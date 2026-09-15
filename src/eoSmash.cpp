#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include "eos.h"
#include "eoSmash.h"

using namespace std;

EoSSmash::EoSSmash(char* filename, int Ne, int Nnb, int Nq) {
  ifstream fin(filename);
  if (!fin.good()) {
   cout << "I/O error with " << filename << endl;
   exit(1);
  }
 // Set members
 ne = Ne;
 nnb = Nnb;
 nnq = Nq;

 // Allocate the interleaved table of T, p, muB, muQ, muS
 tab = new EoSSmashNode[(size_t)ne * nnb * nnq];

 double *e = new double[ne];
 double *nb = new double[nnb];
 double *nq = new double[nnq];

 // Skip header
 string header_str;
 getline(fin, header_str);

 // Read table from file
 for (int ie = 0; ie < ne; ie++)
  for (int inb = 0; inb < nnb; inb++)
   for (int inq = 0; inq < nnq; inq++) {
    // e nb nq T p mub mus muq
    EoSSmashNode &nd = tab[index3(ie, inb, inq)];
    fin >> e[ie] >> nb[inb] >> nq[inq] >> nd.T >> nd.p >> nd.mub >> nd.mus >>
        nd.muq;
   }

// Find upper and lower bounds of e, nB, nQ
emin = e[0];
emax = e[ne - 1];
nbmin = nb[0];
nbmax = nb[nnb - 1];
nqmin = nq[0];
nqmax = nq[nnq - 1];

// Step widths and their reciprocals: constant for the whole run
de = (emax - emin) / (ne - 1);
dnb = (nbmax - nbmin) / (nnb - 1);
dnq = (nqmax - nqmin) / (nnq - 1);
inv_de = 1.0 / de;
inv_dnb = 1.0 / dnb;
inv_dnq = 1.0 / dnq;
p_elow = tab[index3(0, nnb / 2, nnq / 2)].p;
T_elow = tab[index3(0, nnb / 2, nnq / 2)].T;

cout << "EoSSMASH: table " << filename
      << " read, [emin,emax,nbmin,nbmax,qmin,qmax] = " << emin << "  " << emax << "  "
      << nbmin << "  " << nbmax << " " << nqmin << "  " << nqmax << endl;

// e, nb and nq are not needed anymore
delete[] e;
delete[] nb;
delete[] nq;
}

EoSSmash::~EoSSmash() { delete[] tab; }

// Shared cell lookup and trilinear weights for eos() and p().
void EoSSmash::locate(double e, double nb, double nq,
                      const EoSSmashNode *&base, double w[8]) const {
 // Find appropriate indices
 int ie = (int)((e - emin) * inv_de);
 int inb = (int)((nb - nbmin) * inv_dnb);
 int inq = (int)((nq - nqmin) * inv_dnq);
 if (ie < 0) ie = 0;
 if (inb < 0) inb = 0;
 if (inq < 0) inq = 0;
 if (ie > ne - 2) ie = ne - 2;
 if (inb > nnb - 2) inb = nnb - 2;
 if (inq > nnq - 2) inq = nnq - 2;

 // Fraction of a step in the e, nB and nQ directions
 const double em = (e - emin - ie * de) * inv_de;
 const double nbm = (nb - nbmin - inb * dnb) * inv_dnb;
 const double nqm = (nq - nqmin - inq * dnq) * inv_dnq;

 const double we0 = 1. - em, we1 = em;
 const double wb0 = 1. - nbm, wb1 = nbm;
 const double wq0 = 1. - nqm, wq1 = nqm;

 // Eight corner weights, in the (je, jnb, jnq) order of the original loops
 const double a00 = we0 * wb0, a01 = we0 * wb1;
 const double a10 = we1 * wb0, a11 = we1 * wb1;
 w[0] = a00 * wq0;  w[1] = a00 * wq1;
 w[2] = a01 * wq0;  w[3] = a01 * wq1;
 w[4] = a10 * wq0;  w[5] = a10 * wq1;
 w[6] = a11 * wq0;  w[7] = a11 * wq1;

 base = tab + index3(ie, inb, inq);
}

void EoSSmash::eos(double e, double nb, double nq, double ns, double &T,
                    double &mub, double &muq, double &mus, double &p) {

 // e < 0 is physically impossible
 if (e <= 0.) {
  T = mub = muq = mus = p = 0.;
  return;
 }

 // Linearly scale down lowest tabularised value
 if (e < emin) {
  p = e / emin * p_elow;
  T = e / emin * T_elow;
  mub = muq = mus = 0.;
  return;
 }

 const EoSSmashNode *b;
 double w[8];
 locate(e, nb, nq, b, w);

 // strides: +1 in nq, +nnq in nb, +nnq*nnb in e
 const int sb = nnq, se = nnq * nnb;
 const EoSSmashNode *const c0 = b, *const c2 = b + sb;
 const EoSSmashNode *const c4 = b + se, *const c6 = b + se + sb;

 // Compute T, p, muB, muQ, muS by interpolating between neighbouring
 // tabularised points
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

double EoSSmash::p(double e, double nb, double nq, double ns) {
  if (e <= 0.) return 0.0;
  if (e < emin) return e / emin * p_elow;

  const EoSSmashNode *b;
  double w[8];
  locate(e, nb, nq, b, w);

  const int sb = nnq, se = nnq * nnb;
  double p = w[0] * b[0].p + w[1] * b[1].p +
             w[2] * b[sb].p + w[3] * b[sb + 1].p +
             w[4] * b[se].p + w[5] * b[se + 1].p +
             w[6] * b[se + sb].p + w[7] * b[se + sb + 1].p;

  if (p < 0.0) p = 0.0;

  return p;
}
