#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include "eos.h"
#include "eoHadron.h"

using namespace std;

EoSHadron::EoSHadron(char *filename) {
 ifstream fin(filename);
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 fin >> ne >> nnb >> nnq >> e0 >> n0 >> logemin >> logemax >> lognmax;

 const size_t ntot = (size_t)ne * nnb * nnq;
 tab = new EoSHadronNode[ntot];
 statustab = new int[ntot];

 double *e = new double[ne];
 double *nb = new double[nnb];
 double *nq = new double[nnq];

 for (int ie = 0; ie < ne; ie++)
  for (int inb = 0; inb < nnb; inb++)
   for (int inq = 0; inq < nnq; inq++) {
    // e  nb  nq  p  T  mub  muq  mus  status
    const int i3 = index3(ie, inb, inq);
    EoSHadronNode &nd = tab[i3];
    fin >> e[ie] >> nb[inb] >> nq[inq] >> nd.p >> nd.T >> nd.mub >> nd.muq >>
        nd.mus >> statustab[i3];
   }
 double emin = e[0];
 double emax = e[ne - 1];
 double nbmin = nb[0];
 double nbmax = nb[nnb - 1];
 double nqmin = nq[0];
 double nqmax = nq[nnq - 1];

 dxe = (logemax - logemin) / (ne - 1);
 dxnb = (2.0 * lognmax) / (nnb - 1);
 dxnq = (2.0 * lognmax) / (nnq - 1);
 inv_dxe = 1.0 / dxe;
 inv_dxnb = 1.0 / dxnb;
 inv_dxnq = 1.0 / dxnq;
 inv_e0 = 1.0 / e0;
 inv_n0 = 1.0 / n0;
 nb_abs_min = n0 * exp(2.0 * dxnb - lognmax);
 nq_abs_min = n0 * exp(2.0 * dxnq - lognmax);
 e_min = e0 * exp(logemin);
 // reference cell used for the linear ramp below e_min
 p_elow = tab[index3(0, nnb / 2, nnq / 2)].p;
 T_elow = tab[index3(0, nnb / 2, nnq / 2)].T;

 if (fabs((e0 * exp(logemin) - emin) / emin) > 1e-5 ||
     fabs((e0 * exp(logemax) - emax) / emax) > 1e-5 ||
     fabs(n0 * exp(lognmax) - nbmax) > 1e-4) {
  cout << "EoSHadron: wrong eps or nb range: " << setw(14) << emin << setw(14)
       << emax << setw(14) << nbmin << setw(14) << nbmax << endl;
  exit(1);
 }
 cout << "EoHadron: table " << filename
      << " read, [emin,emax,nmin,nmax] = " << emin << "  " << emax << "  "
      << nbmin << "  " << nbmax << "  " << nqmin << "  " << nqmax << endl;
 delete[] e;
 delete[] nb;
 delete[] nq;
}

EoSHadron::~EoSHadron() {
 delete[] tab;
 delete[] statustab;
}

// Shared cell lookup + trilinear weights.  Returns false when the cell is
// flagged as unusable (statustab == 1), matching the original early exits.
bool EoSHadron::locate(double e, double nb, double nq,
                       const EoSHadronNode *&base, double w[8]) const {
 const double xe = log(e * inv_e0);
 double xnb, xnq;  // position of the point in the table
 // Inverse log mapping for nb, nq.  The original always evaluated
 // log(|nb|/n0) and only then tested whether the result was above the
 // threshold.  That test, znb >= -lognmax + 2*dxnb, is exactly
 // |nb| >= n0*exp(2*dxnb - lognmax) == nb_abs_min, which we already hold, so
 // the magnitude can be tested first and the logarithm skipped entirely on
 // the linear branch.  Two of the three logs per call disappear whenever the
 // charge densities are small -- the common case at collider energies.
 const double anb = fabs(nb);
 if (anb >= nb_abs_min) {
  const double znb = log(anb * inv_n0);
  xnb = copysign(0.5 * (znb + lognmax), nb);  // was nb/fabs(nb)*... : a divide
 } else if (nb != 0.0) {
  xnb = nb / nb_abs_min * dxnb;
 } else
  xnb = 0.;
 const double anq = fabs(nq);
 if (anq >= nq_abs_min) {
  const double znq = log(anq * inv_n0);
  xnq = copysign(0.5 * (znq + lognmax), nq);
 } else if (nq != 0.0) {
  xnq = nq / nq_abs_min * dxnq;
 } else
  xnq = 0.;

 int ie = (int)((xe - logemin) * inv_dxe);
 int inb = (int)((xnb + lognmax) * inv_dxnb);
 int inq = (int)((xnq + lognmax) * inv_dxnq);
 if (ie < 0) ie = 0;
 if (inb < 0) inb = 0;
 if (inq < 0) inq = 0;
 if (ie > ne - 2) ie = ne - 2;
 if (inb > nnb - 2) inb = nnb - 2;
 if (inq > nnq - 2) inq = nnq - 2;

 const int i3 = index3(ie, inb, inq);
 if (statustab[i3] == 1) return false;

 const double em = (xe - logemin - ie * dxe) * inv_dxe;
 const double nbm = (xnb + lognmax - inb * dxnb) * inv_dxnb;
 const double nqm = (xnq + lognmax - inq * dxnq) * inv_dxnq;

 const double we0 = 1. - em, we1 = em;
 const double wb0 = 1. - nbm, wb1 = nbm;
 const double wq0 = 1. - nqm, wq1 = nqm;

 // eight corner weights, in the (je, jnb, jnq) order of the original loops
 const double a00 = we0 * wb0, a01 = we0 * wb1;
 const double a10 = we1 * wb0, a11 = we1 * wb1;
 w[0] = a00 * wq0;  w[1] = a00 * wq1;
 w[2] = a01 * wq0;  w[3] = a01 * wq1;
 w[4] = a10 * wq0;  w[5] = a10 * wq1;
 w[6] = a11 * wq0;  w[7] = a11 * wq1;

 base = tab + i3;
 return true;
}

void EoSHadron::eos(double e, double nb, double nq, double ns, double &T,
                    double &mub, double &muq, double &mus, double &p) {
 if (e <= 0.) {
  T = mub = muq = mus = p = 0.;
  return;
 }
 if (e < e_min) {
  p = e / e_min * p_elow;
  T = e / e_min * T_elow;
  mub = muq = mus = 0.;
  return;
 }
 const EoSHadronNode *b;
 double w[8];
 if (!locate(e, nb, nq, b, w)) {
  T = mub = muq = mus = p = 0.;
  return;
 }
 // strides: +1 in nq, +nnq in nb, +nnq*nnb in e
 const int sb = nnq, se = nnq * nnb;
 const EoSHadronNode *const c0 = b, *const c2 = b + sb;
 const EoSHadronNode *const c4 = b + se, *const c6 = b + se + sb;
 p = w[0] * c0[0].p + w[1] * c0[1].p + w[2] * c2[0].p + w[3] * c2[1].p +
     w[4] * c4[0].p + w[5] * c4[1].p + w[6] * c6[0].p + w[7] * c6[1].p;
 T = w[0] * c0[0].T + w[1] * c0[1].T + w[2] * c2[0].T + w[3] * c2[1].T +
     w[4] * c4[0].T + w[5] * c4[1].T + w[6] * c6[0].T + w[7] * c6[1].T;
 mub = w[0] * c0[0].mub + w[1] * c0[1].mub + w[2] * c2[0].mub +
       w[3] * c2[1].mub + w[4] * c4[0].mub + w[5] * c4[1].mub +
       w[6] * c6[0].mub + w[7] * c6[1].mub;
 muq = w[0] * c0[0].muq + w[1] * c0[1].muq + w[2] * c2[0].muq +
       w[3] * c2[1].muq + w[4] * c4[0].muq + w[5] * c4[1].muq +
       w[6] * c6[0].muq + w[7] * c6[1].muq;
 mus = w[0] * c0[0].mus + w[1] * c0[1].mus + w[2] * c2[0].mus +
       w[3] * c2[1].mus + w[4] * c4[0].mus + w[5] * c4[1].mus +
       w[6] * c6[0].mus + w[7] * c6[1].mus;
 if (p < 0.0) p = 0.0;
}

double EoSHadron::p(double e, double nb, double nq, double ns) {
 if (e <= 0.) return 0.0;
 if (e < e_min) return e / e_min * p_elow;
 const EoSHadronNode *b;
 double w[8];
 if (!locate(e, nb, nq, b, w)) return 0.0;
 const int sb = nnq, se = nnq * nnb;
 double p = w[0] * b[0].p + w[1] * b[1].p +
            w[2] * b[sb].p + w[3] * b[sb + 1].p +
            w[4] * b[se].p + w[5] * b[se + 1].p +
            w[6] * b[se + sb].p + w[7] * b[se + sb + 1].p;
 if (p < 0.0) p = 0.0;
 return p;
}
