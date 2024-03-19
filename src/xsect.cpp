#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <TGraph.h>
#include "xsect.h"
#include "inc.h"

using namespace std;

CrossSections::CrossSections(void)
{
 // E_kin [GeV]
 float Ekin [] =
 {0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,
  2.0,  2.5,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0};
  // sigma_E [fm^2]
 float sigmaE [] =
 {0.0,  0.0,  0.011, 0.04,  0.105, 0.16,  0.22,  0.25,  0.28,  0.4,
  0.49,  0.57,  0.64,  0.73,  0.81,  0.865, 0.91,  0.94,  0.97,  0.99};
  // sigma_T [fm^2]
 float sigmaT [] =
 {1.67, 1.48, 1.32, 1.23, 1.26, 1.37, 1.25, 1.10, 1.01, 0.86,
  0.75,  0.68,  0.62,  0.55,  0.49,  0.45,  0.41,  0.38,  0.35,  0.33};
 gSigmaE = new TGraph(sizeof(Ekin)/sizeof(float), Ekin, sigmaE);
 gSigmaT = new TGraph(sizeof(Ekin)/sizeof(float), Ekin, sigmaT);
 // reading tabular data for pim-p total cross section
 ifstream finpimp("tables/pimprot-pdd.dat") ;//14
 //ifstream finpimp("tables/pimprot-pdd-nodelta.dat") ;//15
 //ifstream finpimp("tables/pimprot-delta.dat") ;//16
 if(!finpimp){ cout << "cannot read tables/pimprot-pdd.dat\n"; exit(1) ; }
 string line ;
 istringstream instream ;
 const int dimTb = 251;
 double plab;
 double sqrtsm [dimTb], sigpimp [dimTb];
 getline(finpimp, line) ;
 getline(finpimp, line) ;
 int i=0;
 while( getline(finpimp, line)){
  if(i>dimTb-1){ cout << "CrossSections: please increase dimTb\n"; exit(1);}
  instream.str(line) ;
  instream.seekg(0) ;
  instream.clear() ; // does not work with gcc 4.1 otherwise
  instream >> sqrtsm[i] >> sigpimp[i] >> plab;
  i++;
 }
 cout << "pim-p cross section: " << i << " lines read.\n";
 gSigmaPimp = new TGraph(i, sqrtsm, sigpimp);

 // reading tabular data for pip-p total cross section
 ifstream finpipp("tables/pipprot-pdd.dat");//14
 //ifstream finpipp("tables/pipprot-pdd-nodelta.dat");//15
 //ifstream finpipp("tables/pipprot-delta.dat");//16
 if(!finpipp){ cout << "cannot read tables/pipprot-pdd.dat\n"; exit(1) ; }
 double sqrtsp [dimTb], sigpipp [dimTb];
 getline(finpipp, line) ;
 getline(finpipp, line) ;
 i=0;
 while(getline(finpipp, line)){
  if(i>dimTb-1){ cout << "CrossSections: please increase dimTb\n"; exit(1); }
  instream.str(line) ;
  instream.seekg(0) ;
  instream.clear() ; // does not work with gcc 4.1 otherwise
  instream >> sqrtsp[i] >> sigpipp[i] >> plab;
  i++;
 }
 cout << "pip-p cross section: " << i << " lines read.\n";
 gSigmaPipp = new TGraph(i, sqrtsp, sigpipp);
}

// cross sections [fm^2] as a function of cm energy sqrt(s)
void CrossSections::NN(double sqrts, double& sigmaNN)
{
 if(sqrts<2*mN){ // catching numerical errors
  cout << "p-t friction: sqrt(s)<2m_N \n";
  exit(1);
 }
 double N=0.0975444498660489, a=2.17159752305036, e=0.770755822129685, b=3.74465647535611,
    c=3.85047465058729, k=40.5487506706665, N2=0.0965904826952028, a2=2.34838850523609,
    e2=1.20275622467399, b2=2.83296079993603, c2=4.22157367676681, k2=40.9003640080888,
    a3=-48.6044134976404, b3=237.401929007043, c3=-248.561728492573, border=2.16020645259424,
    M=470.40583836337, b4=40.8939051741805, a4=0.0362331499096052, a5=8.77989569420941,
    b5=5.13458402261944, a6=3.54812896365638, b6=28.3320853923458, k7=38.4919255573532,
    a8=0.0528869122546821, b8=42.0695534993737, sigmapp, sigmapn;
 sigmaNN=0.0;
 if(sqrts<border){
    sigmaNN=M*pow(sqrts-2*mN,2);
 }else if(sqrts<a){
    sigmapp=M*pow(sqrts-2*mN,2);
    sigmapn=a3*sqrts*sqrts+b3*sqrts+c3;
    sigmaNN=0.5*(sigmapp+sigmapn);
 }else if(sqrts<a2){
    sigmapp=N*pow(sqrts-a,e)*std::exp(-b*(sqrts-c))+k;
    sigmapn=a3*sqrts*sqrts+b3*sqrts+c3;
    sigmaNN=0.5*(sigmapp+sigmapn);
 }else if(sqrts<3.5){
    sigmapp=N*pow(sqrts-a,e)*std::exp(-b*(sqrts-c))+k;
    sigmapn=N2*pow(sqrts-a2,e2)*std::exp(-b2*(sqrts-c2))+k2;
    sigmaNN=0.5*(sigmapp+sigmapn);
 }else if(sqrts<7.0){
    sigmapp=b4*pow(sqrts-2.5,-a4);
    sigmapn=b8*pow(sqrts-2.5,-a8);
    sigmaNN=0.5*(sigmapp+sigmapn);
 }else if(sqrts<18.0){
    sigmaNN=k7;
 }else if(sqrts<85.0){
    sigmaNN=a6*std::log(sqrts)+b6;
 }else if(sqrts<1000.0){
    sigmaNN=a5*std::log(sqrts)+b5;
 }
 sigmaNN*=0.1; // converts [mb] -> [fm^2]
}

void CrossSections::Ivanov(double Ekin, double& sigmaT,
 double& sigmaE, double& sigmaP)
{
 if(Ekin<0.){ // catching numerical errors
  cout << "p-t friction: Ekin<0 \n";
  exit(1);
 }
 if(Ekin<=0.2){
  sigmaT = 1.239 + 0.0448/Ekin + 0.00831/(Ekin*Ekin);
  sigmaE = 0.;
 }
 else if(Ekin>0.2 && Ekin<10){
  sigmaT = gSigmaT->Eval(Ekin);
  sigmaE = gSigmaE->Eval(Ekin);
 }
 else if(Ekin>=10. && Ekin<100.){
  sigmaT = 1.365 * pow(Ekin, -0.623);
  sigmaE = 0.464 + (0.257 - 0.0125*log(Ekin)) * log(Ekin);
 }
 else if(Ekin>=100.){
  sigmaT = 0.865 * pow(Ekin, -0.525);
  sigmaE = 1.403 + (-0.15 + 0.0317 * log(Ekin)) * log(Ekin);
 }
 // ... and same formula for sigmaP for all Ekin
 sigmaP = sigmaT + (1. + 2.*mN/Ekin)*sigmaE;
}

void CrossSections::Ivanovbar(double sqrts, double& sigmaEbar, double& sigmaPbar, double& sigmaRbar)
{
 double Ekin=sqrts*sqrts/2/mN-2*mN, E0=sqrts/2, yalpha=std::acosh(E0/mN);//x=sqrt(E0*E0-mN*mN)/E0;
 double sigmaE, sigmaNN;
 NN(sqrts,sigmaNN);
 if(Ekin<0.){ // catching numerical errors
  cout << "p-t friction: Ekin<0 \n";
  exit(1);
 }

 double Aint=A(yalpha), Bint=B(yalpha), Cint=C(yalpha), Dint=D(yalpha), Fint=F(yalpha);
 sigmaEbar=alpha*sigmaNN/2*(Aint-Bint-Dint)/Aint;
 sigmaPbar=sigmaNN/2/Aint*((Aint-Bint)-(1.0-alpha)*Cint-alpha*Fint);
 sigmaRbar=sigmaNN/2*Bint/Aint;
}

double CrossSections::piN(double sqrts)
{
 if(sqrts<mN+mpi){ // catching numerical errors
  cout << "f-p/t friction: sqrt(s)<m_N+m_pi \n";
  exit(1);
 }
 double xsect=0.0;
 double apip=6.74164804191884, bpip=1.91040774819671,
    cpip=22.8185019412713, apim=7.64349878021023, bpim=1.98522045695842, cpim=24.0705206942523,
    api2=0.574008225964179, bpi2=22.7230259927231, sigmapimp, sigmapipp;
 if (sqrts < 1.74) {
  sigmapimp=gSigmaPimp->Eval(sqrts);
  sigmapipp=gSigmaPipp->Eval(sqrts);
  xsect=0.5*(sigmapimp+sigmapipp);
 }else if(sqrts<1.97){
  sigmapimp=apim*pow(std::log(sqrts)-bpim,2)+cpim;
  sigmapipp=gSigmaPipp->Eval(sqrts);
  xsect=0.5*(sigmapimp+sigmapipp);
 }else if(sqrts<10.0){
  sigmapimp=apim*pow(std::log(sqrts)-bpim,2)+cpim;
  sigmapipp=apip*pow(std::log(sqrts)-bpip,2)+cpip;
  xsect=0.5*(sigmapimp+sigmapipp);
 }else if(sqrts<1000.0){
  xsect=api2*std::log(sqrts)+bpi2;
 }
 return 0.1*xsect;  // converts [mb] -> [fm^2]
}
