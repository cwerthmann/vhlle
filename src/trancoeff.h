class EoS;

// this class contains the information about the transport coefficients
// of the fluid: eta/s, zeta/s and the corresponding relaxation times,
// taupi (\tau_\pi) and tauPi (\tau_\Pi)
class TransportCoeff {

 double etaS0, zetaS0, taupi, tauPi, ah, al, aRho, etaSMin, eEtaSMin, T0;
 double etaH_eta0, etaH_eta2, etaH_eta4;
 EoS *eos;  // EoS instance is needed optionally for zeta/s or eta/s parametrization,
 int zetaSparam, etaSparam;
            // which depends on the speed of sound
 double zeta(double e, double s, double T);
 double eta(double e, double p, double nb, double s, double T, double mub);
public:
 TransportCoeff(double _etaS, double _zetaS, int _zetaSparam, EoS *_eos, int _etaSparam, double _ah, double _al,
   double _aRho, double _T0, double _etaSMin, double _eEtaSMin, double _etaH_eta0, double _etaH_eta2, double _etaH_eta4) :
  etaS0(_etaS), zetaS0(_zetaS), zetaSparam(_zetaSparam), eos(_eos), etaSparam(_etaSparam), ah(_ah), al(_al), aRho(_aRho), T0(_T0),
  etaSMin(_etaSMin), eEtaSMin(_eEtaSMin), etaH_eta0(_etaH_eta0), etaH_eta2(_etaH_eta2), etaH_eta4(_etaH_eta4) {};
  //zetaSparam: 0 - basic ,1 ,2 - arxiv:1910.12930, 3 arxiv:2103.09848
 ~TransportCoeff(){};
 void printZetaT();
 // returns (optionally temperature dependent) eta/s and zeta/s
 void getEta(double e, double p, double nb, double s, double T, double mub, double &_eta, double &_zeta);
 // returns shear and bulk relaxation times
 void getTau(double e, double p, double nb, double s, double T, double mub, double &_taupi, double &_tauPi);
 // deltapipi, taupipi, lambdapiPi * divided by tau_pi * !
 void getOther(double e, double nb, double nq, double ns,
   double &deltapipi, double &taupipi, double &lambdapiPi, double &phi7) {
  deltapipi = 4./3.;  taupipi = 0.;  lambdapiPi = 0.;
  phi7 = 0.;
  if(std::isinf(phi7)) phi7=0.0;
 }
 void getOtherBulk(double e, double nb, double nq, double ns,
   double &delPiPi, double &lamPipi) {
  delPiPi = 2./3.;  lamPipi = 0.;
 }
 // isViscous tells whether the fluid is viscous or inviscid
 inline bool isViscous() {
  if (etaS0 > 0. || zetaS0 > 0. || etaSparam>0)
   return true;
  else
   return false;
 }
};
