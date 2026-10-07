#pragma once
#include <vector>

class Fluid;
class EoS;

// Initial state from KoMPoST3D (KoMPoST run on every eta_s slice of an MC-EKRT T^{mu nu}): a binary
// KMPST3D file with e, u^mu, pi^{mu nu} and Pi on the hydro grid at tau0. See icKompost3d.cpp.
class IcKompost3d {
 public:
  // reads the file and checks that its grid and time equal the fluid's grid and tau0
  IcKompost3d(Fluid *f, const char *filename, double tau0);
  void setIC(Fluid *f, EoS *eos);

 private:
  int nxy, neta;
  double tau;
  std::vector<double> data;  // ieta, iy, ix, 16 components
};
