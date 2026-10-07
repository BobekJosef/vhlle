// Reads the KMPST3D file written by KoMPoST/src/Main3D.cpp (KoMPoST3D.exe):
//   char magic[8] = "KMPST3D", int32 version = 1, ncomp = 16, nxy, neta,
//   double tau [fm], xy_max [fm], eta_max, then nxy*nxy*neta*16 doubles ordered
//   ieta (slowest), iy, ix, component (fastest);
//   x, y on nxy points in [-xy_max, xy_max], eta_s on neta points in [-eta_max, eta_max];
//   components e [GeV/fm^3], u^t, u^x, u^y, u^z,
//   pi^{tt}, pi^{tx}, pi^{ty}, pi^{tz}, pi^{xx}, pi^{xy}, pi^{xz}, pi^{yy}, pi^{yz}, pi^{zz},
//   Pi [GeV/fm^3], where z = tau * eta (the convention of Cell::pi and of vz).
// KoMPoST decomposes T^{mu nu} with a conformal EoS, T = e u u - (e/3 + Pi) Delta + pi. To keep the
// same T^{mu nu} with vHLLE's EoS the bulk pressure is set to Pi + e/3 - p(e); with zeta/s = 0 the
// hydro relaxes it to zero in the first step.
// The grid is not interpolated: the fluid must use the file's grid and tau0 = tau.

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "eos.h"
#include "fld.h"
#include "icKompost3d.h"
#include "inc.h"

using namespace std;

IcKompost3d::IcKompost3d(Fluid *f, const char *filename, double tau0) {
  cout << "loading KoMPoST3D IC from " << filename << endl;
  ifstream fin(filename, ios::binary);
  if (!fin.good()) {
    cout << "I/O error with " << filename << endl;
    exit(1);
  }
  char magic[8];
  int32_t head[4];
  double geom[3];
  fin.read(magic, sizeof magic);
  fin.read(reinterpret_cast<char *>(head), sizeof head);
  fin.read(reinterpret_cast<char *>(geom), sizeof geom);
  if (!fin || strncmp(magic, "KMPST3D", 8) != 0 || head[0] != 1 || head[1] != 16) {
    cout << filename << " is not a KMPST3D v1 file with 16 components" << endl;
    exit(1);
  }
  nxy = head[2];
  neta = head[3];
  tau = geom[0];
  const double xy_max = geom[1], eta_max = geom[2];
  cout << "KoMPoST3D grid: " << nxy << " x " << nxy << " x " << neta << ", x,y in [-" << xy_max << ", "
       << xy_max << "], eta_s in [-" << eta_max << ", " << eta_max << "], tau = " << tau << endl;

  // same grid as the fluid, to a small fraction of a cell
  const double tolXY = 1e-4 * f->getDx(), tolZ = 1e-4 * f->getDz();
  const bool same = f->getNX() == nxy && f->getNY() == nxy && f->getNZ() == neta &&
                    fabs(f->getX(0) + xy_max) < tolXY && fabs(f->getX(nxy - 1) - xy_max) < tolXY &&
                    fabs(f->getY(0) + xy_max) < tolXY && fabs(f->getY(nxy - 1) - xy_max) < tolXY &&
                    fabs(f->getZ(0) + eta_max) < tolZ && fabs(f->getZ(neta - 1) - eta_max) < tolZ;
  if (!same) {
    cout << "The hydro grid must equal the KoMPoST3D grid: set nx = ny = " << nxy << ", nz = " << neta
         << ", xmin = ymin = " << -xy_max << ", xmax = ymax = " << xy_max << ", etamin = " << -eta_max
         << ", etamax = " << eta_max << endl;
    exit(1);
  }
  if (fabs(tau - tau0) > 1e-6 * tau) {
    cout << "tau0 = " << tau0 << " must equal the KoMPoST3D output time tOut = " << tau << endl;
    exit(1);
  }

  data.resize(static_cast<size_t>(nxy) * nxy * neta * 16);
  fin.read(reinterpret_cast<char *>(data.data()), data.size() * sizeof(double));
  if (!fin) {
    cout << filename << " is truncated" << endl;
    exit(1);
  }
}

void IcKompost3d::setIC(Fluid *f, EoS *eos) {
  // pi component order in the file, upper triangle row by row
  static const int pi_mu[10] = {0, 0, 0, 0, 1, 1, 1, 2, 2, 3};
  static const int pi_nu[10] = {0, 1, 2, 3, 1, 2, 3, 2, 3, 3};
  const double dx = f->getDx(), dy = f->getDy(), dz = f->getDz();
  double E = 0, Eideal = 0, Pz = 0, Px = 0, Py = 0, S = 0;
  int ncapped = 0;

  for (int iz = 0; iz < neta; iz++)
    for (int iy = 0; iy < nxy; iy++)
      for (int ix = 0; ix < nxy; ix++) {
        const double *r = &data[((static_cast<size_t>(iz) * nxy + iy) * nxy + ix) * 16];
        Cell *c = f->getCell(ix, iy, iz);
        const double e = r[0] > 0 ? r[0] : 0.0;
        double vx = 0, vy = 0, vz = 0;
        if (e > 0) {
          vx = r[2] / r[1];
          vy = r[3] / r[1];
          vz = r[4] / r[1];
          const double v2 = vx * vx + vy * vy + vz * vz;
          if (v2 > 0.99 * 0.99) {  // as ickw.cpp does
            const double scale = 0.99 / sqrt(v2);
            vx *= scale;
            vy *= scale;
            vz *= scale;
            ncapped++;
          }
        }
        c->setPrimVar(eos, tau, e, 0., 0., 0., vx, vy, vz);
        c->saveQprev();

        const double p = eos->p(e, 0., 0., 0.);
        const double Pi = e > 0 ? r[15] + e / 3. - p : 0.0;
        double pi[4][4] = {};
        if (e > 0)
          for (int k = 0; k < 10; k++) pi[pi_mu[k]][pi_nu[k]] = pi[pi_nu[k]][pi_mu[k]] = r[5 + k];
        for (int i = 0; i < 4; i++)
          for (int j = 0; j <= i; j++) {
            c->setpi(i, j, pi[i][j]);
            c->setpiH(i, j, pi[i][j]);
            c->setpi0(i, j, 0.);
            c->setpiH0(i, j, 0.);
          }
        c->setPi(Pi);
        c->setPiH(Pi);
        c->setPi0(0.);
        c->setPiH0(0.);
        c->setAllM(e > 0 ? 1. : 0.);
        if (e <= 0) continue;

        // conserved quantities with and without the viscous part
        const double eta = f->getZ(iz), ch = cosh(eta), sh = sinh(eta);
        const double gamma = 1. / sqrt(1. - vx * vx - vy * vy - vz * vz);
        const double u[4] = {gamma, gamma * vx, gamma * vy, gamma * vz};
        const double P = p + Pi;
        const double Ttt = (e + P) * u[0] * u[0] - P + pi[0][0];
        const double Ttz = (e + P) * u[0] * u[3] + pi[0][3];
        const double dV = tau * dx * dy * dz;
        E += (Ttt * ch + Ttz * sh) * dV;
        Pz += (Ttt * sh + Ttz * ch) * dV;
        Px += ((e + P) * u[0] * u[1] + pi[0][1]) * dV;
        Py += ((e + P) * u[0] * u[2] + pi[0][2]) * dV;
        Eideal += ((e + p) * u[0] * (u[0] * ch + u[3] * sh) - p * ch) * dV;
        S += eos->s(e, 0., 0., 0.) * u[0] * dV;
      }
  data.clear();
  data.shrink_to_fit();
  if (ncapped) cout << "flow velocity capped at |v| = 0.99 in " << ncapped << " cells" << endl;
  cout << "hydrodynamic E = " << E << "  Pz = " << Pz << "  (ideal part E = " << Eideal << ")" << endl
       << "  Px = " << Px << "  Py = " << Py << endl;
  cout << "initial_entropy S_ini = " << S << endl;
}
