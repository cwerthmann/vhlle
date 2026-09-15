//#include <TError.h>
//#include <TApplication.h>
//#include <TGraph.h>
//#include <TCanvas.h>
//#include <TMath.h>
//#include <TGraph.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
//#include <TF1.h>

#include "inc.h"
#include "eos.h"
#include "eoChiral.h"

using namespace std;

// ---- auxiliary EoS class. Two instances (objects) of this class will be
// created
// ---- to store two EoS tables: chiraleos and chiralsmall
// ---- Not performance-wise, but the code is elegant :)
// One table cell.  The four interpolated quantities are stored interleaved so
// that a cell occupies 32 contiguous bytes and the two neighbours (in, in+1)
// needed by the bilinear interpolation share a cache line.  The original
// layout (five separate double** tables) touched up to ten cache lines per
// call behind two dependent pointer loads each.
struct EoSNode {
 double p, T, mub, mus;
};

class EoSaux {
 double emax, nmax, emin, nmin;
 double de, dn;          // grid spacings, formerly recomputed (two divisions)
 double inv_de, inv_dn;  // on every single get()/p() call
 int ne, nn;
 EoSNode* tab;           // flat, row-major: tab[ie * nn + in]

 // Cell lookup and bilinear weights, shared by get() and p(), which used to
 // carry two verbatim copies of this code.
 inline void locate(double e, double nb, const EoSNode*& r0, const EoSNode*& r1,
                    double& w0e, double& w1e, double& w0n, double& w1n) const {
  int ie = (int)((e - emin) * inv_de);
  int in = (int)((nb - nmin) * inv_dn);
  if (ie < 0) ie = 0;
  if (in < 0) in = 0;
  if (ie > ne - 2) ie = ne - 2;
  if (in > nn - 2) in = nn - 2;
  const double em = (e - emin - ie * de) * inv_de;
  const double nm = (nb - nmin - in * dn) * inv_dn;
  w1e = em;  w0e = 1. - em;
  w1n = nm;  w0n = 1. - nm;
  r0 = tab + (size_t)ie * nn + in;
  r1 = r0 + nn;
 }

public:
 EoSaux(const char* filename, int Ne, int Nn);
 ~EoSaux();
 void get(double e, double nb, double& p, double& T, double& mub, double& mus);
 double p(double e, double nb);
};

EoSaux::EoSaux(const char* filename, int Ne, int Nn) {
 ne = Ne;
 nn = Nn;
 tab = new EoSNode[(size_t)ne * nn];
 double* e = new double[ne];
 double* n = new double[nn];
 ifstream fin(filename);
 double a;
 if (!fin) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 double sdummy;
 for (int in = 0; in < nn; in++)
  for (int ie = 0; ie < ne; ie++) {
   // T [MeV], mu_q [MeV], energy density [e_0], pressure [e_0], baryon
   // density [n_0], entropy density [n_0], mu_S [MeV], placeholder.
   EoSNode& nd = tab[(size_t)ie * nn + in];
   fin >> nd.T >> nd.mub >> e[ie] >> nd.p >> n[in] >> sdummy >> nd.mus >> a;
   nd.T /= 1000.0;          // --> T[GeV]
   nd.mub *= 3.0 / 1000.0;  // --> mub[GeV]
   nd.mus /= 1000.0;        // --> mus[GeV]
   nd.p *= 0.146;           // --> p[GeV/fm3]
   // the entropy-density column is read and discarded, as before: nothing
   // ever used stab, EoS::s() recomputes s from e, p, T and the potentials
  }
 emin = e[0] * 0.146;
 emax = e[ne - 1] * 0.146;
 nmin = n[0] * 0.15;
 nmax = n[nn - 1] * 0.15;
 de = (emax - emin) / (ne - 1);
 dn = (nmax - nmin) / (nn - 1);
 inv_de = 1.0 / de;
 inv_dn = 1.0 / dn;
 cout << "EoSaux: table " << filename
      << " read, [emin,emax,nmin,nmax] = " << emin << "  " << emax << "  "
      << nmin << "  " << nmax << endl;
 delete[] e;
 delete[] n;
}

EoSaux::~EoSaux() { delete[] tab; }

void EoSaux::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {
 if (e <= 0. || e <= mN * nb) {
  T = mub = mus = p = 0.;
  return;
 }
 const EoSNode *r0, *r1;
 double w0e, w1e, w0n, w1n;
 locate(e, nb, r0, r1, w0e, w1e, w0n, w1n);
 const double c00 = w0e * w0n, c01 = w0e * w1n;
 const double c10 = w1e * w0n, c11 = w1e * w1n;
 p = c00 * r0[0].p + c01 * r0[1].p + c10 * r1[0].p + c11 * r1[1].p;
 T = c00 * r0[0].T + c01 * r0[1].T + c10 * r1[0].T + c11 * r1[1].T;
 mub = c00 * r0[0].mub + c01 * r0[1].mub + c10 * r1[0].mub + c11 * r1[1].mub;
 mus = c00 * r0[0].mus + c01 * r0[1].mus + c10 * r1[0].mus + c11 * r1[1].mus;
 if (p < 0.0) p = 0.0;
}

double EoSaux::p(double e, double nb) {
 if (e < 0. || e <= mN * nb) return 0.0;
 const EoSNode *r0, *r1;
 double w0e, w1e, w0n, w1n;
 locate(e, nb, r0, r1, w0e, w1e, w0n, w1n);
 double p = w0e * (w0n * r0[0].p + w1n * r0[1].p) +
            w1e * (w0n * r1[0].p + w1n * r1[1].p);
 if (p < 0.0) p = 0.0;
 return p;
}

EoSChiral::EoSChiral() {
 eosbig = new EoSaux("eos/chiraleos.dat", 2001, 401);
 eossmall = new EoSaux("eos/chiralsmall.dat", 201, 201);
}

EoSChiral::~EoSChiral() {
 delete eosbig;
 delete eossmall;
}

void EoSChiral::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
 if (e < 1.46 && nb < 0.3)
  eossmall->get(e, nb, p, T, mub, mus);
 else if (e < 146. && nb < 6.)
  eosbig->get(e, nb, p, T, mub, mus);
 else if (e>=mN*nb){
  p = 0.2964 * e;
  T = 0.15120476935 * pow(e, 0.25);
  mub = mus = 0.0;
 }
 else {
  p = T = mub = mus = 0.;
 }
 muq = 0.0;  // generally it's not zero, but...but
}

double EoSChiral::p(double e, double nb, double nq, double ns) {
 if (e < 1.46 && nb < 0.3)
  return eossmall->p(e, nb);
 else if (e < 146. && nb < 6.)
  return eosbig->p(e, nb);
 else if (e>=mN*nb)
  return 0.2964 * e;
 else return 0.;
}
