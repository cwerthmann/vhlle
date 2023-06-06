#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <TGraph.h>
#include "xsect.h"

using namespace std;

CrossSections::CrossSections(void)
{
 // reading tabular data for pim-p total cross section
 ifstream finpimp("tables/pimprot-pdd.dat") ;
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
  sigpimp[i] = sigpimp[i] * 0.1; // converts [mb] -> [fm^2]
  i++;
 }
 cout << "pim-p cross section: " << i << " lines read.\n";
 gSigmaPimp = new TGraph(i, sqrtsm, sigpimp);

 // reading tabular data for pip-p total cross section
 ifstream finpipp("tables/pipprot-pdd.dat");
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
  sigpipp[i] = sigpipp[i] * 0.1; // converts [mb] -> [fm^2]
  i++;
 }
 cout << "pip-p cross section: " << i << " lines read.\n";
 gSigmaPipp = new TGraph(i, sqrtsp, sigpipp);
}

// cross sections [fm^2] as a function of cm energy sqrt(s)
void CrossSections::NN(double sqrts, double& sigmaNN)
{
 const double mN = 0.939; // nucleon mass [GeV]
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

double CrossSections::piN(double sqrts)
{
 double mN = 0.939;
 double mPi = 0.1396;
 if(sqrts<mN+mPi){ // catching numerical errors
  cout << "f-p/t friction: sqrt(s)<m_N+m_pi \n";
  exit(1);
 }
 double apip=6.74164804191884, bpip=1.91040774819671,
    cpip=22.8185019412713, apim=7.64349878021023, bpim=1.98522045695842, cpim=24.0705206942523,
    api2=0.574008225964179, bpi2=22.7230259927231, sigmapimp, sigmapipp;
 if (sqrts < 1.74) {
  sigmapimp=gSigmaPimp->Eval(sqrts);
  sigmapipp=gSigmaPipp->Eval(sqrts);
  return 0.5*(sigmapimp+sigmapipp);
 }else if(sqrts<1.97){
  sigmapimp=apim*pow(std::log(sqrts)-bpim,2)+cpim;
  sigmapipp=gSigmaPipp->Eval(sqrts);
  return 0.5*(sigmapimp+sigmapipp);
 }else if(sqrts<10.0){
  sigmapimp=apim*pow(std::log(sqrts)-bpim,2)+cpim;
  sigmapipp=apip*pow(std::log(sqrts)-bpip,2)+cpip;
  return 0.5*(sigmapimp+sigmapipp);
 }else if(sqrts<1000.0){
  return api2*std::log(sqrts)+bpi2;
 }
 return 0.0;
}
