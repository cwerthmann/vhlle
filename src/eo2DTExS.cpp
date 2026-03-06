#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <algorithm>

#include "eo2DTExS.h"
#include "eos.h"

using namespace std;

namespace {
constexpr double HBARC = 0.1973269804;               // GeV*fm
constexpr double HBARC3 = HBARC * HBARC * HBARC;
constexpr double GEV4_TO_GEV_FM3 = 1.0 / HBARC3;     // GeV^4 -> GeV/fm^3
constexpr double GEV3_TO_FM3     = GEV4_TO_GEV_FM3;  // GeV^3 -> 1/fm^3
}

EoS2DTExS_aux::EoS2DTExS_aux(const char* filename, int Nt, int Nmb)
  : tildeTMin(1e30), tildeTMax(-1e30), tildeMuBMin(1e30), tildeMuBMax(-1e30)
{
 // Open file and detect grid size; if the first line has Nt Nmb, use them
 ifstream fin(filename);
 // Columns per header:
 // # Ttilde(GeV) muBtilde(GeV) e(GeV^4) nB(GeV^3) T(GeV) muB(GeV) P(GeV^4) s(GeV^3) cs2 chi2(GeV^2) chi
 double tildeT = 0.0, tildeMuB = 0.0, soundSpeedSquared = 0.0, chi2Val = 0.0, chiVal = 0.0;
 double tildeEnergy = 0.0, tildeBaryonDensity = 0.0;
 if (!fin.good()) {
  cout << "I/O error with " << filename << endl;
  exit(1);
 }
 // Try to read the first two integers as grid dimensions
 std::streampos start = fin.tellg();
 int ntFromFile = -1, nmbFromFile = -1;
 if (fin >> ntFromFile >> nmbFromFile) {
  NtGrid = ntFromFile;
  NmbGrid = nmbFromFile;
 } else {
  fin.clear();
  fin.seekg(start);
  NtGrid = Nt;
  NmbGrid = Nmb;
 }
 // Discard the remainder of the line after the dimensions
 fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

 // Allocate once dimensions are known
 ptab = new double*[NtGrid];
 Ttab = new double*[NtGrid];
 mubtab = new double*[NtGrid];
 mustab = new double*[NtGrid];
 stab = new double*[NtGrid];
 for (int i = 0; i < NtGrid; i++) {
  ptab[i] = new double[NmbGrid];
  Ttab[i] = new double[NmbGrid];
  mubtab[i] = new double[NmbGrid];
  mustab[i] = new double[NmbGrid];
  stab[i] = new double[NmbGrid];
 }
 // Table is regular in (Ttilde, muBtilde); iterate Ttilde first, then muBtilde
 std::string line;
 for (int iT = 0; iT < NtGrid; ++iT)
   for (int iMu = 0; iMu < NmbGrid; ++iMu) {
      // Skip comment lines starting with '#'
      while (fin.peek() == '#') std::getline(fin, line);
      // Read in the new column order according to the header
      if (!(fin >> tildeT >> tildeMuB >> tildeEnergy >> tildeBaryonDensity >> Ttab[iT][iMu] >> mubtab[iT][iMu]
                   >> ptab[iT][iMu] >> stab[iT][iMu] >> soundSpeedSquared >> chi2Val >> chiVal)) {
         cout << "Unexpected EoS or format error while reading " << filename << endl;
         exit(1);
      }
      // Convert table units: T, muB stay in GeV; P: GeV^4 -> GeV/fm^3; s: GeV^3 -> 1/fm^3
      ptab[iT][iMu]  *= GEV4_TO_GEV_FM3;
      stab[iT][iMu]  *= GEV3_TO_FM3;
      // No unit conversion: values are already in GeV-based units
      // The table does not provide mu_S; set it to zero
      mustab[iT][iMu] = 0.0;
      // Track tilde grid range
      if (tildeT < tildeTMin) tildeTMin = tildeT;
      if (tildeT > tildeTMax) tildeTMax = tildeT;
      if (tildeMuB < tildeMuBMin) tildeMuBMin = tildeMuB;
      if (tildeMuB > tildeMuBMax) tildeMuBMax = tildeMuB;
    }
 // Grid is regular in tilde coordinates

 cout << "EoS2DTExS: table '" << filename << "'" << endl
         << "  grid size      : NTildeT = " << NtGrid << ", NTildemuB = " << NmbGrid << endl
         << "  Ttilde range   : [" << tildeTMin << ", " << tildeTMax << "]" << endl
         << "  muBtilde range : [" << tildeMuBMin << ", " << tildeMuBMax << "]" << endl;
}

EoS2DTExS_aux::~EoS2DTExS_aux() {
 for (int i = 0; i < NtGrid; i++) {
  delete[] ptab[i];
  delete[] Ttab[i];
  delete[] mubtab[i];
  delete[] mustab[i];
  delete[] stab[i];
 }
 delete[] ptab;
 delete[] Ttab;
 delete[] mubtab;
 delete[] mustab;
 delete[] stab;
}

void EoS2DTExS_aux::get(double e, double nb, double& p, double& T, double& mub,
                 double& mus) {

   //cout << "EoS2DTExS_aux::get called with e=" << e << ", nb=" << nb << endl;
 if (e < 0.) {
  T = mub = mus = p = 0.;
  return;
 }
 // Convert inputs (GeV/fm^3, 1/fm^3) to natural units expected by the table (GeV^4, GeV^3)
 const double energyGeV4  = e  * HBARC3;
 const double baryonDensityGeV3 = nb * HBARC3;

 // Convert (e, nb) to tilde variables based on conformal relations
 const double tildeT = pow((12.0 * energyGeV4) / (19.0 * M_PI * M_PI), 0.25);
 if (tildeT <= 0.0) {
    T = mub = mus = p = 0.0;
    return;
 }
 double tildeMuB = 5.0 * baryonDensityGeV3 / (tildeT * tildeT); // based on https://arxiv.org/abs/2406.11610

 if (NtGrid < 2 || NmbGrid < 2) {
  T = mub = mus = p = 0.0;
  return;
 }


 const double tildeTStep   = (tildeTMax - tildeTMin) / (NtGrid - 1);
 const double tildeMuBStep = (tildeMuBMax - tildeMuBMin) / (NmbGrid - 1);

 double clampedTildeT  = tildeT;
 double clampedTildeMuB = tildeMuB;
 if (clampedTildeT  < tildeTMin) clampedTildeT  = tildeTMin;
 if (clampedTildeT  > tildeTMax) clampedTildeT  = tildeTMax;
 if (clampedTildeMuB < tildeMuBMin) clampedTildeMuB = tildeMuBMin;
 if (clampedTildeMuB > tildeMuBMax) clampedTildeMuB = tildeMuBMax;

 int iT = static_cast<int>((clampedTildeT - tildeTMin) / tildeTStep);
 int iMu = static_cast<int>((clampedTildeMuB - tildeMuBMin) / tildeMuBStep);
 if (iT < 0) iT = 0;
 if (iMu < 0) iMu = 0;
 if (iT > NtGrid - 2) iT = NtGrid - 2;
 if (iMu > NmbGrid - 2) iMu = NmbGrid - 2;

 const double tildeTOffset  = clampedTildeT  - tildeTMin - iT * tildeTStep;
 const double tildeMuBOffset = clampedTildeMuB - tildeMuBMin - iMu * tildeMuBStep;

 double weightsT[2] = {1.0 - tildeTOffset / tildeTStep, tildeTOffset / tildeTStep};
 double weightsMuB[2] = {1.0 - tildeMuBOffset / tildeMuBStep, tildeMuBOffset / tildeMuBStep};

 T = mub = mus = p = 0.0;
 for (int jT = 0; jT < 2; ++jT)
   for (int jMu = 0; jMu < 2; ++jMu) {
    const double weight = weightsT[jT] * weightsMuB[jMu];
    p   += weight * ptab[iT + jT][iMu + jMu];
    T   += weight * Ttab[iT + jT][iMu + jMu];
    mub += weight * mubtab[iT + jT][iMu + jMu];
    mus += weight * mustab[iT + jT][iMu + jMu];
   }
  if (p < 0.0) p = 0.0;
}

double EoS2DTExS_aux::p(double e, double nb) {
 if (e < 0.) return 0.0;
 const double energyGeV4  = e  * HBARC3;
 const double baryonDensityGeV3 = nb * HBARC3;

 const double tildeT = pow((12.0 * energyGeV4) / (19.0 * M_PI * M_PI), 0.25);
 if (tildeT <= 0.0) return 0.0;
 double tildeMuB = 3.0 * baryonDensityGeV3 / (tildeT * tildeT);

 if (NtGrid < 2 || NmbGrid < 2) return 0.0;


 const double tildeTStep   = (tildeTMax - tildeTMin) / (NtGrid - 1);
 const double tildeMuBStep = (tildeMuBMax - tildeMuBMin) / (NmbGrid - 1);



 double clampedTildeT  = tildeT;
 double clampedTildeMuB = tildeMuB;
 if (clampedTildeT  < tildeTMin) clampedTildeT  = tildeTMin;
 if (clampedTildeT  > tildeTMax) clampedTildeT  = tildeTMax;
 if (clampedTildeMuB < tildeMuBMin) clampedTildeMuB = tildeMuBMin;
 if (clampedTildeMuB > tildeMuBMax) clampedTildeMuB = tildeMuBMax;

 int iT = static_cast<int>((clampedTildeT - tildeTMin) / tildeTStep);
 int iMu = static_cast<int>((clampedTildeMuB - tildeMuBMin) / tildeMuBStep);
 if (iT < 0) iT = 0;
 if (iMu < 0) iMu = 0;
 if (iT > NtGrid - 2) iT = NtGrid - 2;
 if (iMu > NmbGrid - 2) iMu = NmbGrid - 2;

 const double tildeTOffset  = clampedTildeT  - tildeTMin - iT * tildeTStep;
 const double tildeMuBOffset = clampedTildeMuB - tildeMuBMin - iMu * tildeMuBStep;

 double weightsT[2] = {1.0 - tildeTOffset / tildeTStep, tildeTOffset / tildeTStep};
 double weightsMuB[2] = {1.0 - tildeMuBOffset / tildeMuBStep, tildeMuBOffset / tildeMuBStep};

 double pres = 0.0;
 for (int jT = 0; jT < 2; ++jT)
    for (int jMu = 0; jMu < 2; ++jMu)
      pres += weightsT[jT] * weightsMuB[jMu] * ptab[iT + jT][iMu + jMu];

 if (pres < 0.0) pres = 0.0;
 return pres;
}

double EoS2DTExS::p(double e, double nb) {
 if (e < 0.) return 0.0;
 return eos_TExS->p(e, nb);
}

EoS2DTExS::EoS2DTExS(const char* filename, int Nt, int Nmb) {
 eos_TExS = new EoS2DTExS_aux(filename, Nt, Nmb);
}

EoS2DTExS::~EoS2DTExS() {
 delete eos_TExS;
}

void EoS2DTExS::eos(double e, double nb, double nq, double ns, double& T,
                    double& mub, double& muq, double& mus, double& p) {
 eos_TExS->get(e, nb, p, T, mub, mus);
 muq = 0.0; 
}

double EoS2DTExS::p(double e, double nb, double nq, double ns) {
 return eos_TExS->p(e, nb);
}
