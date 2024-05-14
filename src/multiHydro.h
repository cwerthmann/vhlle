#include <vector>
#include <cmath>

class Fluid;
class Hydro;
class EoS;
class TransportCoeff;
class CrossSections;
class Cornelius;
class Nucleon;

class MultiHydro {
 Fluid *f_p, *f_t, *f_f; // projectile and target fluids
 Hydro *h_p, *h_t, *h_f; // solutions of hydro eqs. for projectile and target
 EoS *eos;
 Cornelius *cornelius;
 TransportCoeff *trcoeff;
 CrossSections *xsect;
 std::ofstream fmhfreeze_p, fmhfreeze_f, fmhfreeze_t, ffricx, ffricy, ffricz, ffricall, EfIfile;
 std::string EfIfilename;
 double ***MHeps, ***MHepsPrev;
 double **EfITable;
 int NTemp, Nvatilde,xsectparam;
 double Tmax;
 std::vector<std::vector<double>> retardedFriction;
 double ecrit, vEff_p, vEff_t, vEff_f;
 int nx, ny, nz;
 double dx, dy, dz, dtau, tau0, sNN, Etot, Q0min;
 double xi_fa, formationTime, lambda, xi_q, xi_h;
 int frictionModel, decreasingFormTime, unification;
 double dtauf;
 double tau_unification;
 const double gmunu[4][4] = {
     {1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}};
 double EtotSurf[3] = {0., 0., 0.}, EtotSurf_positive[3] = {0., 0., 0.},
     EtotSurf_negative[3] = {0., 0., 0.};
 std::vector<std::vector<Nucleon>> nucleons;

public:
 MultiHydro(Fluid *f_p, Fluid *f_t, Fluid *f_f, Hydro *h_p, Hydro *h_t,
  Hydro* h_f, EoS *eos, TransportCoeff *trcoeff, double dtau, double eCrit, double sNN, double Etot,
  double xi_fa, double lambda, double formationTime, int frictionModel, int decreasingFormTime,
  double xi_q, double xi_h, double alpha, double beta, int unification, double tau_unification, int NTemp, int Nvatilde, double Tmax, int xsectparam, std::vector<std::vector<Nucleon>> nucl);
 ~MultiHydro(void);
 void setDtau(double newdtau);
 void initOutput(const char *dir);
 double EfIeval(double Tf,double vatilde);
 void performStep();
 void frictionSubstep();
 void getEnergyMomentumTensor(double (&T)[4][4], double Q_p[7], double Q_f[7], double Q_t[7]);
 void getDiagonalizedEnergyDensity();
 void getMaxEnergyDensity();
 void getEnergyDensity();
 void getSumEnergyDensity();
 void updateEnergyDensity();
 int  findFreezeout(EoS *eosH);
 void printFreezeout(std::ofstream &fout, double t, double x, double y, double z, double dsigma[4], double uC[4], double TC, double mub, double muq, double mus, double picart[10], double PiC, double dVEff);
 void outputEnergyDensity();
 void resizeMHeps();
 void setFluids(Fluid *f_p, Fluid *f_t, Fluid *f_f, Hydro *h_p, Hydro *h_t,
  Hydro* h_f);
 void addRetardedFriction(double flux, double x, double y, double z, double t, int i);
 double calculateRetardedFriction(double x, double y, double z, double t, int i);
 void clearRetardedFriction();
 double totalScatRate(double px, double T, double mu, double u[4]);
 double calculateScatRates(double px, double T, double mu, double u[4], int particle);
 double pp_total(double mandelstam_s);
 double Fermi(double nb);
 double getLocalEnergyDensity(double x, double y, double eta);
 void evolveSpectators(void);
 void printSpectators(std::ofstream &fout);
 inline double expk2(double z){
    if(z>300.0) return std::sqrt(M_PI/2.0/z)*(1.0+15.0/8.0/z+105.0/128.0/z/z-945.0/3072.0/z/z/z);
    return std::exp(z)*std::cyl_bessel_k(2,z);
 }
};
