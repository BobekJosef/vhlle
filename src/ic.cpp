#include <fstream>
#include <iomanip>
#include <iomanip>
#include <TF1.h>
#include <TF2.h>

#include "fld.h"
#include "eos.h"
#include "ic.h"
#include "inc.h"
#include "s95p.h"
#include "icurqmd.h"

using namespace std;

const double A = 197.0;
const double Ra = 6.37;
const double dlt = 0.54;
const double sigma = 5.1;

const int nphi = 301;

extern char icInputFile[255];
extern int icModel;
extern double s0ScaleFactor;

IC::IC(double e, double impactPar, double a) {
  epsilon = e;
  alpha = a;
  b = impactPar;
}

IC::~IC(void) {}

double IC::eProfile(double x, double y) {
  if (icModel == 1) {  // analytical Glauber for ...
    prms[0] = sqrt((x + b / 2.0) * (x + b / 2.0) + y * y);
    const double tpp = iff->Integral(-3.0 * Ra, 3.0 * Ra, prms, 1.0e-9);
    prms[0] = sqrt((x - b / 2.0) * (x - b / 2.0) + y * y);
    const double tmm = iff->Integral(-3.0 * Ra, 3.0 * Ra, prms, 1.0e-9);
    return epsilon *
           pow(1. / rho0 * (tpp * (1.0 - pow((1.0 - sigma * tmm / A), A)) +
                            tmm * (1.0 - pow((1.0 - sigma * tpp / A), A))),
               1.0);
  } else if (icModel == 2) {
    return s95p::getSongEps(x, y);
  }
}

void IC::findRPhi(void) {
  _rphi = new double[nphi];
  for (int iphi = 0; iphi < nphi; iphi++) {
    double phi = iphi * C_PI * 2. / (nphi - 1);
    double r, r1 = 0., r2 = 2. * Ra;
    while (fabs((r2 - r1) / r2) > 0.001 && r2 > 0.001) {
      r = 0.5 * (r1 + r2);
      if (eProfile(r * cos(phi), r * sin(phi)) > 0.5)
        r1 = r;
      else
        r2 = r;
    }
    _rphi[iphi] = r;
  }
}

double IC::rPhi(double phi) {
  const double cpi = C_PI;
  phi = phi - 2. * cpi * floor(phi / 2. / cpi);
  int iphi = (int)(phi / (2. * cpi) * (nphi - 1));
  int iphi1 = iphi + 1;
  if (iphi1 == nphi) iphi = nphi - 2;
  return _rphi[iphi] * (1. - (phi / (2. * cpi) * (nphi - 1) - iphi)) +
         _rphi[iphi1] * (phi / (2. * cpi) * (nphi - 1) - iphi);
}

void IC::setIC(Fluid *f, EoS *eos, double tau) {
  double e, nb, nq, e_corr, p, vx = 0., vy = 0., vz = 0.;
  const double be = 1.;
  const double bn = 1.;
  Cell *c;
  ofstream fvel("velocity_debug.txt");

  TF2 *ff = 0;
  double prms2[2], intgr2;
  cout << "finding normalization constant\n";
  ff = new TF2("ThicknessF", this, &IC::Thickness, -3.0 * Ra, 3.0 * Ra,
               -3.0 * Ra, 3.0 * Ra, 2, "IC", "Thickness");
  prms2[0] = Ra;
  prms2[1] = dlt;
  ff->SetParameters(prms2);
  intgr2 = ff->Integral(-3.0 * Ra, 3.0 * Ra, -3.0 * Ra, 3.0 * Ra, 1.0e-9);
  if (intgr2 == 0.0) {
    cerr << "IC::setICGlauber Error! ff->Integral == 0; Return -1\n";
    delete ff;
    exit(1);
  }
  delete ff;
  cout << "a = " << A / intgr2 << endl;
  prms[1] = A / intgr2;
  prms[2] = Ra;
  prms[3] = dlt;
  iff = new TF1("WoodSaxonDF", this, &IC::WoodSaxon, -3.0 * Ra, 3.0 * Ra, 4,
                "IC", "WoodSaxon");
  prms[0] = 0.0;
  const double tpp = iff->Integral(-3.0 * Ra, 3.0 * Ra, prms, 1.0e-9);
  rho0 = 2.0 * tpp * (1.0 - pow((1.0 - sigma * tpp / A), A));
  if (icModel == 2) {
    s95p::loadSongIC(icInputFile, s0ScaleFactor);
    // s95p::loadSongIC("/home/yura/Downloads/results-CGC-part-0005-id/Initial/KLN00-05_Part.dat2",
    // 17.75) ;
  }
  if (icModel == 3) {
    icurqmd::loadIC(icInputFile);
  }

  findRPhi();  // fill in R(phi) table
  cout << "R(phi) =  ";
  for (int jj = 0; jj < 5; jj++) cout << rPhi(jj * C_PI / 2.) << "  ";  // test
  cout << endl;
  //--------------
  double avv_num = 0., avv_den = 0.;
  double Etotal = 0.0, Stotal = 0.0;

  for (int ix = 0; ix < f->getNX(); ix++)
    for (int iy = 0; iy < f->getNY(); iy++)
      for (int iz = 0; iz < f->getNZ(); iz++) {
        c = f->getCell(ix, iy, iz);
        double x = f->getX(ix);
        double y = f->getY(iy);
        double z = f->getZ(iz);

        if (icModel == 3) {
          icurqmd::getIC(x, y, z, e, nb, nq, vx, vy, vz);
          //    e = e*1.5 ;
          // vz = -tanh(eta) ; // for v_z=0 (in Cartesian coordinates)
          if (e < 1e-10) {
            e = nb = nq = 0.0;
            vx = vy = vz = 0.0;
          }
        } else {
          //double eta1 = fabs(eta) < 1.3 ? 0.0 : fabs(eta) - 1.3;
          //e = eProfile(x, y) * exp(-eta1 * eta1 / 2.1 / 2.1) *
          //    (fabs(eta) < 5.3 ? 1.0 : 0.0);  // see Huovinen,Hirano
          // e = eProfile(x,y)*(fabs(eta)<2.0 ? 1.0 : 0.0) ; // box-like in
          // z-dir
          // nb = nq = 0.0;
          // if (e < 0.5) e = 0;
          // vx = vy = 0.;
          // vz = 0.;
        }
        //------- test
         e=30.*exp( - x*x -y*y - z*z ) ;
         nb = nq = 0.0;
         if(e<0.01) e=0. ;
         vx = vy = vz = 0. ;
        //------------

        avv_num += sqrt(vx * vx + vy * vy) * e;
        avv_den += e;

        c->setPrimVar(eos, e, nb, nq, 0., vx, vy, vz);
        double _p = eos->p(e, nb, nq, 0.);
        double _s = eos->s(e, nb, nq, 0.);
        const double gamma2 = 1.0 / (1.0 - vx * vx - vy * vy - vz * vz);
        Etotal +=
            ((e + _p) * gamma2 - _p );
        c->saveQprev();
        Stotal += _s * sqrt(gamma2);

        if (e > 0.) c->setAllM(1.);
      }
  fvel.close();
  cout << "average initial flow = " << avv_num / avv_den << endl;
  cout << "Etot = " << Etotal *f->getDx() * f->getDy() * f->getDz() * tau
       << "  Stot = " << Stotal * f->getDx() * f->getDy() * f->getDz() * tau
       << endl;
  //	exit(1) ;
}

double IC::Thickness(double *x, double *p) {
  // p[0]: Ra radius; p[1]: delta = 0.54fm
  double intgrl, prms[4];
  TF1 *iff = 0;
  iff = new TF1("WoodSaxonDF", this, &IC::WoodSaxon, -3.0 * p[0], 3.0 * p[0], 4,
                "IC", "WoodSaxon");
  prms[0] = sqrt(x[0] * x[0] + x[1] * x[1]);  //
  prms[1] = 1.0;  // normalization parameter which must be found.
  prms[2] = p[0];
  prms[3] = p[1];
  intgrl = iff->Integral(-3.0 * p[0], 3.0 * p[0], prms, 1.0e-9);
  if (iff) delete iff;
  return intgrl;
}

double IC::WoodSaxon(double *x, double *p) {
  return p[1] / (exp((sqrt(x[0] * x[0] + p[0] * p[0]) - p[2]) / p[3]) + 1.0);
}
