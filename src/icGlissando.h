
class Fluid;
class EoS;

class IcGlissando {
private:
 double eta0; // midrapidity plateau
 double sigEta; // diffuseness of rapidity profile
 double ybeam; // beam rapidity
 double alphaMix; // WN/binary mixing
 double Rg; // Gaussian smearing in transverse dir
 double sNorm; // normalization of initial entropy profile
 double A; // initial shear flow constant
 int nx, ny, nz, nevents;
 double xmin, xmax, ymin, ymax, zmin, zmax;
 double dx, dy, dz;
 double ***rho;
 static const int NP = 10000;  // dimension for particle arrays
 // auxiliary particle arrays
 double X[NP], Y[NP], W[NP];
 int C[NP];

 double tau0;
 int nsmoothx;  // smoothly distribute to +- this many cells
 int nsmoothy;
 int nsmoothz;
 void makeSmoothTable(int npart);

public:
 IcGlissando(Fluid *f, const char *filename, double tau0, const char* setup);
 ~IcGlissando();
 void setIC(Fluid *f, EoS *eos);
};
