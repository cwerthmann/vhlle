#include <cmath>
#include <fstream>
#include <iostream>
#include <TMatrixDEigen.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <iomanip>
#include <Math/Functor.h>
#include <Math/GaussIntegrator.h>
#include <omp.h>

#include "multiHydro.h"
#include "hdo.h"
#include "fld.h"
#include "eos.h"
#include "rmn.h"
#include "trancoeff.h"
#include "cll.h"
#include "xsect.h"
#include "cornelius.h"
#include "nucleon.h"
#include "EfI.h"

#define OUTPI

using namespace std;

MultiHydro::MultiHydro(Fluid *_f_p, Fluid *_f_t, Fluid *_f_f, Hydro *_h_p,
 Hydro *_h_t, Hydro *_h_f, EoS *_eos, TransportCoeff *_trcoeff, double _dtau,
 double eCrit, double _sNN, double _Etot, double _xi_fa, double _lambda, double _formationTime,
 int _frictionModel, int _decreasingFormTime, double _xi_q, double _xi_h, int _unification, double _tau_unification, int _NTemp, int _Nvatilde, double _Tmax, int _xsectparam, std::vector<std::vector<Nucleon>> nucl)
{
 f_p = _f_p;
 f_t = _f_t;
 f_f = _f_f;
 h_p = _h_p;
 h_t = _h_t;
 h_f = _h_f;
 eos = _eos;
 trcoeff = _trcoeff;
 xsect = new CrossSections;
 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 dtau = _dtau;
 tau0 = h_p->getTau();
 sNN = _sNN;
 Etot = _Etot;
 xi_fa = _xi_fa;
 lambda = _lambda;
 formationTime = _formationTime;
 frictionModel = _frictionModel;
 decreasingFormTime = _decreasingFormTime;
 xi_q = _xi_q;
 xi_h = _xi_h;
 unification=_unification;
 tau_unification=_tau_unification;
 dtauf = formationTime / 10.0;
 NTemp=_NTemp;
 Nvatilde=_Nvatilde;
 Tmax=_Tmax;
 xsectparam=_xsectparam;
 EfIfilename="tables/EfItableNTemp";
 EfIfilename.append(to_string(NTemp));
 EfIfilename.append("Nvatilde");
 EfIfilename.append(to_string(Nvatilde));
 EfIfilename.append("Tmax");
 EfIfilename.append(to_string(Tmax));
 EfIfilename.append("xsectparam");
 EfIfilename.append(to_string(xsectparam));
 EfIfilename.append(".dat");

 nucleons = nucl;

 //---- Cornelius init
 double arrayDx[4] = {h_p->getDtau(), f_p->getDx(), f_p->getDy(), f_p->getDz()};
 cornelius = new Cornelius;
 cornelius->init(4, eCrit, arrayDx);
 ecrit = eCrit;
 vEff = 0.;
 vEff_p = 0.;
 vEff_t = 0.;
 vEff_f = 0.;

 // allocate field for oveall energy density
 MHeps = new double**[nx];
 MHepsPrev = new double**[nx];
 for (int ix = 0; ix < nx; ix++) {
  MHeps[ix] = new double*[ny];
  MHepsPrev[ix] = new double*[ny];
  for (int iy = 0; iy < ny; iy++) {
   MHeps[ix][iy] = new double[nz];
   MHepsPrev[ix][iy] = new double[nz];
   for (int iz = 0; iz < nz; iz++) {
    MHeps[ix][iy][iz] = 0.0;
    MHepsPrev[ix][iy][iz] = 0.0;
   }
  }
 }

 Q0min=1e-6*Etot/(dx)/(dy)/(dz)/f_p->getNX()/f_p->getNY()/f_p->getNZ();
 cout<<"min Q0 level: "<<Q0min<<endl;
 ifstream EfIstream(EfIfilename);
 if(!EfIstream.good()){
  //compute table of fireball friction integrals if necessary
  cout<<"Producing EfITable with filename "<< EfIfilename<<endl;
  ofstream EfIoutstream;
  EfIoutstream.open(EfIfilename);
  EfIoutstream<<setprecision(15);

  //double smin=pow(mN+mpi,2), smax;
  double pmin=mpi, pmax;
  ROOT::Math::GaussIntegrator GI;
  GI.SetRelTolerance(1e-5);
  double vatilde=0.0;
  double gammaatilde=1.0;
  double Tf=0.0;
  double zeta3 = 1.20205690315959;

  int progresslast=0;
  int progress=0;

  EfITable=new double*[Nvatilde];
  for(int iv=0;iv<Nvatilde;iv++){
   EfITable[iv]=new double[NTemp];
   vatilde=1.0*iv/Nvatilde;
   gammaatilde=1.0/sqrt(1.0-vatilde*vatilde);
   for(int iT=0;iT<NTemp;iT++){
    Tf=Tmax*iT/(NTemp-1);
    EfIntegrand EfI(xsect,Tf,0.0,vatilde);
    ROOT::Math::Functor1D func(&EfI,&EfIntegrand::EvalNpi);
    GI.SetFunction(func);

    //smax=20*mN*Tf/gammaatilde/(1-vatilde)+smin;
    pmax=10*Tf+pmin;
    EfITable[iv][iT]=GI.Integral(pmin,pmax);
    double saa=mN*mN+mpi*mpi+2.0*mN*mpi*gammaatilde;
    //EfITable[iv][iT]=3.0/12.0*Tf*Tf*mN*mpi*sqrt(gammaatilde-1)*xsect->piN(sqrt(saa))*gevtofm*gevtofm*gevtofm;//13
    //EfITable[iv][iT]=3.0*zeta3*Tf*Tf*Tf/M_PI/M_PI*mN*sqrt(gammaatilde-1)*xsect->piN(sqrt(saa))*gevtofm*gevtofm*gevtofm;//11

    EfIoutstream<<EfITable[iv][iT]<<endl;
   }
   progress=floor(100.0*iv/Nvatilde);
   if(progress>progresslast){
    progresslast=progress;
    cout<<progresslast<<"% done"<<endl;
   }
  }
  EfIoutstream.close();
  cout<<"EfITable written and ready"<<endl;


 }else{
 //read table of fireball friction integrals
 cout<<"Reading EfITable from file "<<EfIfilename<<endl;
 string line;
 EfITable=new double*[Nvatilde];
 for(int iv=0;iv<Nvatilde;iv++){
  EfITable[iv]=new double[NTemp];
  for(int iT=0;iT<NTemp;iT++){
  if(!getline(EfIstream,line)){
    cout<<"Error reading EfITable on line "<<iv*NTemp+iT<<endl;
    exit(1);
  }
  EfITable[iv][iT]=stod(line);
  }
 }
  cout<<"EfITable ready"<<endl;
 }

 /* Debug of scattering rates calculation
 double u[4] = {1, 0, 0, 0};
 for (double p = 0.1; p <= 2; p+=0.1) {
  cout << p << " " << totalScatRate(p, 0.12, 0.3, u) << endl;
 }
 exit(1);*/
}

MultiHydro::~MultiHydro() {
 for (int ix = 0; ix < nx; ix++) {
  for (int iy = 0; iy < ny; iy++) {
   delete[] MHeps[ix][iy];
   delete[] MHepsPrev[ix][iy];
  }
  delete[] MHeps[ix];
  delete[] MHepsPrev[ix];
 }
 delete[] MHeps;
 delete[] MHepsPrev;
}


 void MultiHydro::setDtau(double newdtau){
 h_p->setDtau(newdtau);
 h_t->setDtau(newdtau);
 h_f->setDtau(newdtau);

 dtau=newdtau;

 double arrayDx[4] = {h_p->getDtau(), f_p->getDx(), f_p->getDy(), f_p->getDz()};
 delete cornelius;
 cornelius = new Cornelius;
 cornelius->init(4, ecrit, arrayDx);
}



void MultiHydro::setFluids(Fluid *_f_p, Fluid *_f_t, Fluid *_f_f, Hydro *_h_p,
 Hydro *_h_t, Hydro *_h_f) {
 f_p = _f_p;
 f_t = _f_t;
 f_f = _f_f;
 h_p = _h_p;
 h_t = _h_t;
 h_f = _h_f;
 nx = f_p->getNX();
 ny = f_p->getNY();
 nz = f_p->getNZ();
 dx = f_p->getDx();
 dy = f_p->getDy();
 dz = f_p->getDz();
 dtau = h_p->getDtau();

 //---- Cornelius init
 double arrayDx[4] = {h_p->getDtau(), f_p->getDx(), f_p->getDy(), f_p->getDz()};
 cornelius = new Cornelius;
 cornelius->init(4, ecrit, arrayDx);

 resizeMHeps();
}

void MultiHydro::initOutput(const char *dir) {
 string outfreeze_p = dir, outfreeze_f = dir, outfreeze_t = dir, outfricx=dir, outfricy=dir, outfricz=dir, outfricall=dir;
 outfreeze_p.append("/freezeout_p.dat");
 outfreeze_t.append("/freezeout_t.dat");
 outfreeze_f.append("/freezeout_f.dat");
 outfricx.append("/fricx.dat");
 outfricy.append("/fricy.dat");
 outfricz.append("/fricz.dat");
 outfricall.append("/fricall.dat");
 fmhfreeze_p.open(outfreeze_p.c_str());
 fmhfreeze_t.open(outfreeze_t.c_str());
 fmhfreeze_f.open(outfreeze_f.c_str());
 ffricx.open(outfricx.c_str());
 ffricy.open(outfricy.c_str());
 ffricz.open(outfricz.c_str());
 ffricall.open(outfricall.c_str());
}

double MultiHydro::EfIeval(double Tf, double vatilde){
 int iT=floor((NTemp-1)*Tf/Tmax);
 if((NTemp-1)*Tf>2.1e9){
  iT=NTemp-1;
 }
 int iv=floor(Nvatilde*vatilde);
 double lambdav=vatilde*Nvatilde-1.0*iv;
 double lambdaT=Tf/Tmax*(NTemp-1)-1.0*iT;
 if(iv>=Nvatilde-1&&iT>=NTemp-1){
  return EfITable[Nvatilde-1][NTemp-1];
 }
 if(iv>=Nvatilde-1){
  double low=EfITable[Nvatilde-1][iT];
  double high=EfITable[Nvatilde-1][iT+1];
  return lambdaT*high+(1.0-lambdaT)*low;
 }
 if(iT>=NTemp-1){
  double low=EfITable[iv][NTemp-1];
  double high=EfITable[iv+1][NTemp-1];
  return lambdav*high+(1.0-lambdav)*low;
 }
 double lowlow=EfITable[iv][iT];
 double lowhigh=EfITable[iv][iT+1];
 double highlow=EfITable[iv+1][iT];
 double highhigh=EfITable[iv+1][iT+1];
 return lambdav*lambdaT*highhigh+lambdav*(1.0-lambdaT)*highlow+(1.0-lambdav)*lambdaT*lowhigh+(1.0-lambdav)*(1.0-lambdaT)*lowlow;
}

void MultiHydro::performStep()
{
 #pragma omp parallel
 {
 #pragma omp sections
 {
 #pragma omp section
 {
 h_p->performStep();
 }
 #pragma omp section
 {
 h_t->performStep();
 }
 #pragma omp section
 {
 h_f->performStep();
 }
 }
 }
 frictionSubstep();

}


void MultiHydro::frictionSubstep()
{
 double mindtaufric=dtau;
 int NLimitedFriction=0;
 int Nunphys=0;
 double ELimitedFriction=0.0;
 double EtotFriction=0.0;
 double EtotLimited=0.0;
 int NSkip=0;
 int Nloop=0;
 // here it is assumed that projectile and target grids
 // have same dimensions and physical sizes
 //#pragma omp parallel for num_threads(3) collapse(3)
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    double dtaufric;
    double dtaufrictot=0.0;

    while(dtaufrictot<dtau){
    dtaufric=dtau-dtaufrictot;
    double flux_p [4] = {0.}, flux_t [4] = {0.},
           flux_pf [4] = {0.}, flux_tf [4] = {0.},
           nbflux_p=0.0, nbflux_t=0.0, nbflux_pf=0.0, nbflux_tf=0.0;
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    double _Q_p[7], _Q_t[7], _Q_f[7], Q_p_new[7]={0}, Q_t_new[7]={0}, Q_f_new[7]={0};
    c_p->getQ(_Q_p);
    c_t->getQ(_Q_t);
    c_f->getQ(_Q_f);


    double taup = h_p->getTau()+dtaufrictot;
    double taut = h_t->getTau()+dtaufrictot;
    double tauf = h_f->getTau()+dtaufrictot;
    double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
    double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
    double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
    c_p->getPrimVar(eos, taup, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
    c_t->getPrimVar(eos, taut, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
    c_f->getPrimVar(eos, tauf, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);
    if(_Q_p[0]>Q0min||_Q_t[0]>Q0min){
    double TCp, mubCp, muqCp, musCp, pCp;
    double TCt, mubCt, muqCt, musCt, pCt;
    double TCf, mubCf, muqCf, musCf, pCf;
    eos->eos(ep, nbp, nqp, nsp, TCp, mubCp, muqCp, musCp, pCp);
    eos->eos(et, nbt, nqt, nst, TCt, mubCt, muqCt, musCt, pCt);
    eos->eos(ef, nbf, nqf, nsf, TCf, mubCf, muqCf, musCf, pCf);
    // 4-velocities, u_p and u_t
    double gammap = 1.0/sqrt(1.0-vxp*vxp-vyp*vyp-vzp*vzp);
    double up [4] = {gammap,gammap*vxp,gammap*vyp,gammap*vzp};
    double gammat = 1.0/sqrt(1.0-vxt*vxt-vyt*vyt-vzt*vzt);
    double ut [4] = {gammat,gammat*vxt,gammat*vyt,gammat*vzt};
    double gammaf = 1.0/sqrt(1.0-vxf*vxf-vyf*vyf-vzf*vzf);
    double uf [4] = {gammaf,gammaf*vxf,gammaf*vyf,gammaf*vzf};
    double upt_abs = sqrt(pow(up[0]+ut[0],2)-pow(up[1]+ut[1],2)-pow(up[2]+ut[2],2)-pow(up[3]+ut[3],2));
    double U_F[4] = {(up[0]+ut[0])/upt_abs, (up[1]+ut[1])/upt_abs, (up[2]+ut[2])/upt_abs, (up[3]+ut[3])/upt_abs};
    double vfsq=vxf*vxf+vyf*vyf+vzf*vzf;
    double vpsq=vxp*vxp+vyp*vyp+vzp*vzp;
    double vtsq=vxt*vxt+vyt*vyt+vzt*vzt;
    double vfvp=vxf*vxp+vyf*vyp+vzf*vzp;
    double vfvt=vxf*vxt+vyf*vyt+vzf*vzt;
    double vTf=2.0*TCf/mpi*(1.0+TCf/mpi)/expk2(mpi/TCf);
    double vptilde=sqrt(max(0.0,1.0-(1.0-vfsq)*(1.0-vpsq)/(1.0-vfvp)/(1.0-vfvp)));//sqrt(abs(vfsq+vpsq-2.0*vfvp+vfvp*vfvp-vfsq*vpsq))/abs(1.0-vfvp);
    double vttilde=sqrt(max(0.0,1.0-(1.0-vfsq)*(1.0-vtsq)/(1.0-vfvt)/(1.0-vfvt)));//sqrt(abs(vfsq+vtsq-2.0*vfvt+vfvt*vfvt-vfsq*vtsq))/abs(1.0-vfvt);
    double unification_factor_vp=0.0;
    if(vptilde/vTf<1e-5){
        unification_factor_vp=1-vptilde*vptilde/vTf/vTf-0.5*pow(vptilde/vTf,4);
    }else if(vptilde<vTf){
        unification_factor_vp=exp(1.0/(1.0-vTf*vTf/vptilde/vptilde));
    }
    double unification_factor_vt=0.0;
    if(vttilde/vTf<1e-5){
        unification_factor_vt=1-vttilde*vttilde/vTf/vTf-0.5*pow(vttilde/vTf,4);
    }else if(vttilde<vTf){
        unification_factor_vt=exp(1.0/(1.0-vTf*vTf/vttilde/vttilde));
    }
    double uput = gammap*gammat*(1.0 - vxp*vxt - vyp*vyt - vzp*vzt);
    double savg = 2.0*mN*mN*(1.0 + uput);


     // 1. projectile-target friction
    if (ep>0. && et>0.) {
    // qgb densities
    double zeta3 = 1.20205690315959;
    double nquark_p = 18 * zeta3 * pow(TCp, 3) / pow(M_PI, 2) + 2*pow(mubCp/3, 3);
    double ngluon_p = 16 * zeta3 * pow(TCp, 3) / pow(M_PI, 2);
    double nquark_t = 18 * zeta3 * pow(TCt, 3) / pow(M_PI, 2) + 2*pow(mubCt/3, 3);
    double ngluon_t = 16 * zeta3 * pow(TCt, 3) / pow(M_PI, 2);
    double dens_p, dens_t;
    //if (ep < 0.7) {
     //dens_p = xi_h*pow(2*mN/sqrt(savg), 0.5)*nbp;
     dens_p = xi_h*nbp;
    /*} else {
     dens_p = xi_q*pow(2*mN/sqrt(s), 0.5)*(nquark_p+ngluon_p)/3;
    }*/
    //if (et < 0.7) {
     //dens_t = xi_h*pow(2*mN/sqrt(savg), 0.5)*nbt;
     dens_t = xi_h*nbt;
    /*} else {
     dens_t = xi_q*pow(2*mN/sqrt(s), 0.5)*(nquark_t+ngluon_t)/3;
    }*/
    // Moeller factor
    double Vrel = sqrt(uput*uput - 1.0);



    if (frictionModel == 1||frictionModel==2) {
    double sigmaNN;
    xsect->NN(std::sqrt(savg),sigmaNN);
    // friction coefficient
    double D_N = mN*Vrel*sigmaNN;
     for(int i=0; i<4; i++){
      // Csernai Tmunu friction terms
      flux_p[i] += -dens_p*dens_t*up[i]*D_N;
      flux_t[i] += -dens_p*dens_t*ut[i]*D_N;

      if (formationTime > 0) {
       addRetardedFriction((-flux_p[i]-flux_t[i])*f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau(),
         f_p->getX(ix)+U_F[1]*formationTime, f_p->getY(iy)+U_F[2]*formationTime,
         f_p->getZ(iz)+U_F[3]*formationTime, h_p->getTau()+U_F[0]*formationTime, i);
       //flux_f[i] += calculateRetardedFriction(f_p->getX(ix), f_p->getY(iy), f_p->getZ(iz),
       //            h_p->getTau(), i)/(f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau());
      } else {
       //flux_f[i] += -flux_p[i]-flux_t[i];
      }
     }
     //Csernai nb friction terms
      nbflux_p += -dens_p*dens_t/mN*D_N;
      nbflux_t += -dens_p*dens_t/mN*D_N;
      if (formationTime > 0) {

      } else {
       //nbflux_f += -nbflux_p-nbflux_t;
      }
    } else if(frictionModel==3){

    double Ekin = savg/(2.0*mN) - 2.0*mN;
    double sigmaT, sigmaE, sigmaP;
    if (Ekin <= 0.) {
     Ekin = 0;
     sigmaT = 0;
     sigmaE = 0;
     sigmaP = 0;
    } else {
     xsect->Ivanov(Ekin, sigmaT, sigmaE, sigmaP);
    }
    double D_P = mN*Vrel*sigmaP;
    double D_E = mN*Vrel*sigmaE;



    for(int i=0; i<4; i++){
        // Ivanov's friction terms
      flux_p[i] += -dens_p*dens_t*(D_P*(up[i] - ut[i]) + D_E*(up[i] + ut[i]));
      flux_t[i] += -dens_p*dens_t*(D_P*(ut[i] - up[i]) + D_E*(up[i] + ut[i]));

      if (formationTime > 0) {
       addRetardedFriction((-flux_p[i]-flux_t[i])*f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau(),
         f_p->getX(ix)+U_F[1]*formationTime, f_p->getY(iy)+U_F[2]*formationTime,
         f_p->getZ(iz)+U_F[3]*formationTime, h_p->getTau()+U_F[0]*formationTime, i);
       //flux_f[i] += calculateRetardedFriction(f_p->getX(ix), f_p->getY(iy), f_p->getZ(iz),
       //            h_p->getTau(), i)/(f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau());
      } else {
       //flux_f[i] += -flux_p[i]-flux_t[i];
      }
    }
    nbflux_p=0.0;nbflux_t=0.0;

    } else if(frictionModel==4){
    double sigmaRbar, sigmaEbar, sigmaPbar;
    if (sqrt(savg) <= 2*mN) {
     sigmaRbar = 0;
     sigmaEbar = 0;
     sigmaPbar = 0;
    } else {
     xsect->Ivanovbar(sqrt(savg), sigmaEbar, sigmaPbar, sigmaRbar);
    }
    double D_Pbar = mN*Vrel*sigmaPbar;
    double D_Ebar = mN*Vrel*sigmaEbar;
    double D_Rbar = mN*Vrel*sigmaRbar;

     for(int i=0; i<4; i++){
        // Ivanov's friction terms
      flux_p[i] += -dens_p*dens_t*(D_Pbar*(up[i] - ut[i]) + D_Ebar*(up[i] + ut[i])+D_Rbar*up[i]);
      flux_t[i] += -dens_p*dens_t*(D_Pbar*(ut[i] - up[i]) + D_Ebar*(up[i] + ut[i])+D_Rbar*ut[i]);

      if (formationTime > 0) {
       addRetardedFriction((-flux_p[i]-flux_t[i])*f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau(),
         f_p->getX(ix)+U_F[1]*formationTime, f_p->getY(iy)+U_F[2]*formationTime,
         f_p->getZ(iz)+U_F[3]*formationTime, h_p->getTau()+U_F[0]*formationTime, i);
       //flux_f[i] += calculateRetardedFriction(f_p->getX(ix), f_p->getY(iy), f_p->getZ(iz),
       //            h_p->getTau(), i)/(f_p->getDx()*f_p->getDy()*f_p->getDz()*h_p->getTau());
      } else {
       //flux_f[i] += -flux_p[i]-flux_t[i];
      }
    }
      nbflux_p += -dens_p*dens_t/mN*D_Rbar;
      nbflux_t += -dens_p*dens_t/mN*D_Rbar;
      if (formationTime > 0) {

      } else {
       //nbflux_f += -nbflux_p-nbflux_t;
      }

    } else {

    }
   }
   // 2. projectile-fireball friction
   if(ep>0. && ef>0.) {
    double dens_p = xi_fa*nbp;
    //double gammaptilde=1.0/std::sqrt(1.0-vptilde*vptilde);
    double EfNpi=MultiHydro::EfIeval(TCf,vptilde);
    double EfNN=0.0;
    //cerr << EfNpi<<" "<< dens_p <<" "<< TCf<<" " << mubCf<<" "<<vptilde<<" "<< xsect->piN(2*smin) <<endl;
    /*if(frictionModel==2){
        ROOT::Math::Functor1D func2(&EfI,&EfIntegrand::EvalNN);
        GI.SetFunction(func2);
        smin=pow(mN+mN,2);
        EfNN=GI.Integral(smin,smax);
    }*/

    if(xsectparam<17){
    for(int i=0; i<4; i++){//14
     flux_pf[i] += -dens_p*up[i]*(EfNpi+EfNN);
    }
    nbflux_pf += -dens_p/mN*(EfNpi+EfNN);

    }else if(xsectparam==17){
    for(int i=0; i<4; i++){//17
     flux_pf[i] += dens_p/mN*uf[i]*(EfNpi+EfNN);
    }
    nbflux_pf += 0.0;
   }

   }

    /*double upuf = gammap*gammaf*(1.0 - vxp*vxf - vyp*vyf - vzp*vzf);
    double savgpf = 2.0*mN*mpi*uput+mN*mN+mpi*mpi;
    double sigmaNpi=xsect->piN(std::sqrt(savg));
    double Vrel = sqrt(uput*uput - 1.0);
    // friction coefficient
    double D_N = mN*Vrel*sigmaNpi;
    for(int i=0; i<4; i++){
      // Csernai Tmunu friction terms
      flux_p[i] += -xi_fa*nb*up[i]*D_N*h_p->getDtau();
     }
     //Csernai nb friction terms
      nbflux_p += -dens_p*dens_t/mN*D_N*h_p->getDtau();*/



   // 3. target-fireball friction
   if(et>0. && ef>0.) {
    double dens_t = xi_fa*nbt;
    //double gammattilde=1.0/std::sqrt(1.0-vttilde*vttilde);
    double EfNpi=MultiHydro::EfIeval(TCf,vttilde);
    double EfNN=0.0;
    /*if(frictionModel==2){
        ROOT::Math::Functor1D func2(&EfI,&EfIntegrand::EvalNN);
        GI.SetFunction(func2);
        smin=pow(mN+mN,2);
        EfNN=GI.Integral(smin,smax);
    }*/
    if(xsectparam<17){
    for(int i=0; i<4; i++){//14
     flux_tf[i] += -dens_t*ut[i]*(EfNpi+EfNN);
    }
    nbflux_tf += -dens_t/mN*(EfNpi+EfNN);
   }else if(xsectparam==17){

    for(int i=0; i<4; i++){//17
     flux_tf[i] += dens_t/mN*uf[i]*(EfNpi+EfNN);
    }
    nbflux_tf += 0.0;
    }


   }





    //timestep scaling according to friction size
   if(_Q_p[0]>Q0min){
    dtaufric=min(dtaufric,0.1*_Q_p[0]/taup/abs(flux_p[0]+flux_pf[0]));
   }
   if(_Q_t[0]>Q0min){
    dtaufric=min(dtaufric,0.1*_Q_t[0]/taut/abs(flux_t[0]+flux_tf[0]));
   }
   if(dtaufric<mindtaufric&&dtaufric<dtau-dtaufrictot){
    mindtaufric=dtaufric;
   }


   for(int i=0;i<4;i++){
    flux_p[i]*=dtaufric;
    flux_pf[i]*=dtaufric;
    flux_t[i]*=dtaufric;
    flux_tf[i]*=dtaufric;
   }
   nbflux_p*=dtaufric;
   nbflux_pf*=dtaufric;
   nbflux_t*=dtaufric;
   nbflux_tf*=dtaufric;


    double unification_factor_t=dtaufric/tau_unification;
    if(unification_factor_t>1e-5) unification_factor_t=1.0-exp(-dtaufric/tau_unification);



    //friction limiter

    /*if(ix==f_p->getNX()/2&&iy==f_p->getNY()/2&&abs(iz-f_p->getNZ()/2)<1){
        cout << iz << setw(14) << ep << setw(14) << e_p_new << setw(14) <<nbp <<setw(14)<<nb_p_new<< setw(14) << ep/mN/nbp << setw(14) << e_p_new/mN/nb_p_new << setw(14) << pp <<setw(14) <<p_p_new<< endl;
        cout << flux_p[0] << setw(14)<< flux_p[1] << setw(14)<< flux_p[2] << setw(14)<< flux_p[3] << setw(14) << nbflux_p << endl;
        cout << up[0] << setw(14) << up[1] << setw(14) << up[2] <<setw(14)<<up[3]<<setw(14)<< 1/mN <<endl;
        cout << _Q_p[0] << setw(14) << _Q_p[1] << setw(14) << _Q_p[2] << setw(14) << _Q_p[3] << setw(14) << _Q_p[NB_] << endl;
        cout << taup*((ep+pp)*up[0]*up[0]-pp) << setw(14) << taup*((ep+pp)*up[0]*up[1]) << setw(14)<< taup*((ep+pp)*up[0]*up[2]) << setw(14)<< taup*((ep+pp)*up[0]*up[3]) << setw(14) << taup*gammap*nbp << endl;
    }*/

          EtotFriction+=abs(flux_p[0]+flux_pf[0])+abs(flux_t[0]+flux_tf[0]);

   //double energy_balance=min(e_p_new-1.2*mN*nb_p_new,min(e_t_new-1.2*mN*nb_t_new,e_f_new-1.2*mN*nb_f_new));
   //if (energy_balance >= 0 &&
       //if(_Q_p[T_] + (flux_p[0]+flux_pf[0])*taup >= 0 &&
      // _Q_t[T_] + (flux_t[0]+flux_tf[0])*taut >= 0 &&
      // _Q_f[T_] + (-flux_pf[0]-flux_tf[0]+flux_f[0])*tauf >= 0) {
        if(_Q_p[T_] + (flux_p[0]+flux_pf[0])*taup < 0||_Q_t[T_] + (flux_t[0]+flux_tf[0])*taut < 0||_Q_f[T_] + (-flux_pf[0]-flux_tf[0]-flux_p[0]-flux_t[0])*tauf < 0){
          EtotLimited+=abs(flux_p[0]+flux_pf[0])+abs(flux_t[0]+flux_tf[0]);
          NLimitedFriction++;
         }
        if(_Q_p[T_] + (flux_p[0]+flux_pf[0])*taup < 0){
         if(_Q_p[T_] + (flux_p[0])*taup < 0){
          ELimitedFriction+=abs(flux_p[0]+flux_pf[0])-0.5*_Q_p[T_]/taup;
          for(int i=0;i<4;i++){
          flux_p[i]=-0.5*_Q_p[i]/taup;
          flux_pf[i]=0;
         }
          nbflux_p=-0.5*_Q_p[NB_]/taup;
          nbflux_pf=0;
         }else{
          ELimitedFriction+=abs(flux_pf[0]);
          for(int i=0;i<4;i++){
          flux_pf[i]=0;
          }
          nbflux_pf=0;
         }
        }
        if(_Q_t[T_] + (flux_t[0]+flux_tf[0])*taut < 0){
         if(_Q_t[T_] + (flux_t[0])*taup < 0){
          ELimitedFriction+=abs(flux_t[0]+flux_tf[0])-0.5*_Q_t[T_]/taut;
          for(int i=0;i<4;i++){
          flux_t[i]=-0.5*_Q_t[i]/taut;
          flux_tf[i]=0;
         }
          nbflux_t=-0.5*_Q_t[NB_]/taut;
          nbflux_tf=0;
         }else{
          ELimitedFriction+=abs(flux_tf[0]);
          for(int i=0;i<4;i++){
          flux_tf[i]=0;
          }
          nbflux_tf=0;
         }
        }
        if(_Q_f[T_] + (-flux_pf[0]-flux_tf[0]-flux_p[0]-flux_t[0])*tauf < 0){
         ELimitedFriction+=abs(flux_t[0]+flux_tf[0]);
         ELimitedFriction+=abs(flux_p[0]+flux_pf[0]);
         for(int i=0;i<4;i++){
          flux_t[i]=0;
          flux_p[i]=0;
          flux_tf[i]=0;
          flux_pf[i]=0;
         }
         nbflux_p=0;
         nbflux_t=0;
         nbflux_tf=0;
         nbflux_pf=0;
         }


    for(int i=0;i<7;i++){
    Q_p_new[i]=_Q_p[i]/taup;
    Q_t_new[i]=_Q_t[i]/taup;
    Q_f_new[i]=_Q_f[i]/taup;
   }
   for(int i=0;i<4;i++){
    Q_p_new[i]+=(flux_p[i]+flux_pf[i]);
    Q_t_new[i]+=(flux_t[i]+flux_tf[i]);
    Q_f_new[i]+=(-flux_p[i]-flux_t[i]-flux_pf[i]-flux_tf[i]);
   }
    Q_p_new[4]+=(nbflux_p+nbflux_pf);
    Q_t_new[4]+=(nbflux_t+nbflux_tf);
    Q_f_new[4]+=(-nbflux_p-nbflux_t-nbflux_pf-nbflux_tf);
    double e_p_new, e_t_new, e_f_new, nb_p_new, nb_t_new, nb_f_new, p_p_new, p_t_new, p_f_new, nq_new, ns_new, vx_p_new, vy_p_new, vz_p_new, vx_t_new, vy_t_new, vz_t_new, vx_f_new, vy_f_new, vz_f_new, eratio;
    transformPV(eos,Q_p_new,e_p_new,p_p_new,nb_p_new,nq_new,ns_new,vx_p_new,vy_p_new,vz_p_new,false);
    transformPV(eos,Q_t_new,e_t_new,p_t_new,nb_t_new,nq_new,ns_new,vx_t_new,vy_t_new,vz_t_new,false);
    transformPV(eos,Q_f_new,e_f_new,p_f_new,nb_f_new,nq_new,ns_new,vx_f_new,vy_f_new,vz_f_new,false);



       if(e_p_new<mN*nb_p_new||e_t_new<mN*nb_t_new){
         EtotLimited+=abs(flux_p[0]+flux_pf[0])+abs(flux_t[0]+flux_tf[0]);
         NLimitedFriction++;
       }


        if(e_p_new<mN*nb_p_new){
         Nunphys++;
         eratio=(flux_p[0]+flux_pf[0])/_Q_p[0];

         for(int i=1;i<4;i++){
          flux_p[i]=eratio*_Q_p[i];
          flux_pf[i]=0;
         }
          nbflux_p=eratio*_Q_p[NB_];
          nbflux_pf=0;
        }
        if(e_t_new<mN*nb_t_new){
         Nunphys++;
         eratio=(flux_t[0]+flux_tf[0])/_Q_t[0];

         for(int i=0;i<4;i++){
          flux_t[i]=eratio*_Q_t[i];
          flux_tf[i]=0;
         }
          nbflux_t=eratio*_Q_p[NB_];
          nbflux_tf=0;
        }

    //unification
    if(unification==1){
    for(int i=0;i<4;i++){
     flux_p[i]*=(1.0-unification_factor_vp);
     flux_t[i]*=(1.0-unification_factor_vt);
     flux_pf[i]*=(1.0-unification_factor_vp);
     flux_tf[i]*=(1.0-unification_factor_vt);
    }
    nbflux_p*=(1.0-unification_factor_vp);
    nbflux_t*=(1.0-unification_factor_vt);
    nbflux_pf*=(1.0-unification_factor_vp);
    nbflux_tf*=(1.0-unification_factor_vt);
    for(int i=0;i<4;i++){
        flux_p[i]-=_Q_p[i]/taup*unification_factor_vp*unification_factor_t;
        flux_t[i]-=_Q_t[i]/tauf*unification_factor_vt*unification_factor_t;
    }
    nbflux_p-=_Q_p[4]/taup*unification_factor_vp*unification_factor_t;
    nbflux_t-=_Q_t[4]/tauf*unification_factor_vt*unification_factor_t;
    }

    c_p->addFlux((flux_p[0]+flux_pf[0])*taup, (flux_p[1]+flux_pf[1])*taup,
     (flux_p[2]+flux_pf[2])*taup, (flux_p[3]+flux_pf[3])*taup,(nbflux_p+nbflux_pf)*taup, 0., 0.);
    c_t->addFlux((flux_t[0]+flux_tf[0])*taut, (flux_t[1]+flux_tf[1])*taut,
     (flux_t[2]+flux_tf[2])*taut, (flux_t[3]+flux_tf[3])*taut, (nbflux_t+nbflux_tf)*taut, 0., 0.);
    c_f->addFlux((-flux_pf[0]-flux_tf[0]-flux_p[0]-flux_t[0])*tauf, (-flux_pf[1]-flux_tf[1]-flux_p[1]-flux_t[1])*tauf,
     (-flux_pf[2]-flux_tf[2]-flux_p[2]-flux_t[2])*tauf, (-flux_pf[3]-flux_tf[3]-flux_p[3]-flux_t[3])*tauf, (-nbflux_pf-nbflux_tf-nbflux_p-nbflux_t)*tauf, 0., 0.);
    c_p->updateByFlux();
    c_t->updateByFlux();
    c_f->updateByFlux();
    c_p->clearFlux();
    c_t->clearFlux();
    c_f->clearFlux();


   /*} else {
        NLimitedFriction++;
        ELimitedFriction+=abs(flux_t[0]+flux_tf[0])+abs(flux_p[0]+flux_pf[0]);
        for(int i=0;i<4;i++){
            flux_f[i]=0;
            flux_p[i]=0;
            flux_t[i]=0;
            flux_pf[i]=0;
            flux_tf[i]=0;
        }
        nbflux_f=0;
        nbflux_p=0;
        nbflux_t=0;
        nbflux_pf=0;
        nbflux_tf=0;
   }*/
   }else{
    if(dtaufric==dtau){
    NSkip++;
    }
   }
   //friction output
   //X direction
   /*if(iy==ny/2&&iz==nz/2){
    double x=f_f->getX(ix);
    ffricx << setw(14) << tauf << setw(14) << dtaufric << setw(14) << x;
    for(int i=0;i<4;i++){
        ffricx << setw(14) << flux_p[i] << setw(14) << flux_t[i] << setw(14) << flux_pf[i] << setw(14) << flux_tf[i];
    }
    ffricx << setw(14) << nbflux_p << setw(14) << nbflux_t << setw(14) << nbflux_pf << setw(14) << nbflux_tf;
    ffricx << setw(14) << ep << setw(14) << et << setw(14) << nbp << setw(14) << nbt << setw(14) << vzp << setw(14) << vzt << endl;
   }
    //Y direction
   if(ix==nx/2&&iz==nz/2){
    double y=f_f->getY(iy);
    ffricy << setw(14) << tauf << setw(14) << dtaufric << setw(14) << y;
    for(int i=0;i<4;i++){
        ffricy << setw(14) << flux_p[i] << setw(14) << flux_t[i] << setw(14) << flux_pf[i] << setw(14) << flux_tf[i];
    }
    ffricy << setw(14) << nbflux_p << setw(14) << nbflux_t << setw(14) << nbflux_pf << setw(14) << nbflux_tf;
    ffricy << setw(14) << ep << setw(14) << et << setw(14) << nbp << setw(14) << nbt << setw(14) << vzp << setw(14) << vzt << endl;
   }
    //Z direction
   if(iy==ny/2&&ix==nx/2){
    double z=f_f->getZ(iz);
    ffricz << setw(14) << tauf << setw(14) << dtaufric << setw(14) << z;
    for(int i=0;i<4;i++){
        ffricz << setw(14) << flux_p[i] << setw(14) << flux_t[i] << setw(14) << flux_pf[i] << setw(14) << flux_tf[i];
    }
    ffricz << setw(14) << nbflux_p << setw(14) << nbflux_t << setw(14) << nbflux_pf << setw(14) << nbflux_tf;
    ffricz << setw(14) << ep << setw(14) << et << setw(14) << nbp << setw(14) << nbt << setw(14) << vzp << setw(14) << vzt << endl;
   }
   //central region
   if(iy>ny/2-5&&iy<6*ny/10&&ix>3*nx/10&&ix<7*nx/10&&iz>9*nz/20&&iz<11*nz/20){
   double x=f_f->getX(ix);
   double y=f_f->getY(iy);
   double z=f_f->getZ(iz);
    ffricall << setw(14) << tauf << setw(14) << dtaufric << setw(14) << x << setw(14) << y << setw(14) << z;
    for(int i=0;i<4;i++){
        ffricall << setw(14) << flux_p[i] << setw(14) << flux_t[i] << setw(14) << flux_pf[i] << setw(14) << flux_tf[i];
    }
    ffricall << setw(14) << nbflux_p << setw(14) << nbflux_t << setw(14) << nbflux_pf << setw(14) << nbflux_tf;
    ffricall << setw(14) << ep << setw(14) << et << setw(14) << nbp << setw(14) << nbt << setw(14) << vzp << setw(14) << vzt << endl;
   }*/


   if(-flux_p[0]-flux_t[0] > 0. && c_f->getMaxM()<0.01)
    c_f->setAllM(1.0);

    Nloop++;
    dtaufrictot+=dtaufric;
   } // end dtau while loop
   } // end cell loop
   ffricx << endl;
   ffricy << endl;
   ffricz << endl;
 clearRetardedFriction();
 cout << "Friction update done at tau="<<h_p->getTau()<<", average loop number: "<<1.0*(Nloop-NSkip)/(f_p->getNX()*f_p->getNY()*f_p->getNZ()-NSkip)<<", smallest dtaufric: "<<mindtaufric<<"."<<endl;
 cout << "skipped friction due to small energy density of target and projectile in "<< NSkip << " cells ("<< 100.0*NSkip/f_p->getNX()/f_p->getNY()/f_p->getNZ() << "%)"<<endl;
 cout << "friction limited "<<NLimitedFriction<<" times (" << 100.0*NLimitedFriction/(Nloop-NSkip) << "%) with " << Nunphys << "unphysical fluid cells, total energy transfer of "<< ELimitedFriction*h_p->getTau()*dx*dy*dz << " (" <<100.0*ELimitedFriction/EtotFriction<<"%)"<< endl;
// <<"), of which only partially dropped: "<< 100.0*NPartiallyLimitedFriction/Nloop <<"% ("<<NPartiallyLimitedFriction<<" times, total energy transfer of "<< EPartiallyLimitedFriction*h_p->getTau()*dx*dy*dz <<")"<<endl;
 if(EtotLimited>0){
 cout << "average local fraction of dropped energy transfer: "<< 100.0*ELimitedFriction/EtotLimited << "%"<<endl;
 }
 if (decreasingFormTime == 1) {
  formationTime -= dtau * dtauf;
  if (formationTime < 0) formationTime = 0;
 }
}

void MultiHydro::addRetardedFriction(double flux, double x, double y, double z, double t, int i)
{
 vector<double> v = {flux, x, y, z, t, (double)i};
 if (abs(flux) > 1e-3) retardedFriction.push_back(v);
}


double MultiHydro::calculateRetardedFriction(double x, double y, double z, double t, int i)
{
 double total_flux = 0.;
 for (int it = 0; it < retardedFriction.size(); it++)
 {
  if ((abs(retardedFriction[it][1]-x) < dx) && (abs(retardedFriction[it][2]-y) < dy) &&
   (abs(retardedFriction[it][3]-z) < dz) && (abs(retardedFriction[it][4]-t) < dtau) &&
   (retardedFriction[it][5] == (double)i))
   {
    double wCenX = 1-abs(retardedFriction[it][1]-x)/dx;
    double wCenY = 1-abs(retardedFriction[it][2]-y)/dy;
    double wCenZ = 1-abs(retardedFriction[it][3]-z)/dz;
    double wCenT = 1-abs(retardedFriction[it][4]-t)/dtau;
    total_flux += retardedFriction[it][0] * wCenX * wCenY * wCenZ * wCenT;
   }
 }
 return total_flux;
}

void MultiHydro::clearRetardedFriction()
{
 retardedFriction.erase(
  remove_if(
   retardedFriction.begin(),
   retardedFriction.end(),
   [&](vector<double> element){
    if (element[4] < h_p->getTau()) return true;
    else return false;
   }
  ),
  retardedFriction.end()
 );
 // next part is only for printing energy stored in the retardedFriction vector
 /*double storedEnergy = 0.;
 for (int it = 0; it < retardedFriction.size(); it++)
 {
  if (abs(retardedFriction[it][4]-h_p->getTau()) < dtau)
  {
   storedEnergy += retardedFriction[it][0]*(abs(retardedFriction[it][4]-h_p->getTau())/dtau);
  }
  else
  {
   storedEnergy += retardedFriction[it][0];
  }
 }
 cout << "Stored energy: " << storedEnergy << endl;*/
}

void MultiHydro::getEnergyMomentumTensor(double (&T)[4][4], double Q_p[7], double Q_f[7], double Q_t[7])
{
 const double delta[4][4] = {
     {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
 double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
 double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
 double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
 transformPV(eos, Q_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
 transformPV(eos, Q_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
 transformPV(eos, Q_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);
 // 4-velocities, u_p and u_t
 double gammap = 1.0/sqrt(1.0-vxp*vxp-vyp*vyp-vzp*vzp);
 double up [4] = {gammap,gammap*vxp,gammap*vyp,gammap*vzp};
 double gammat = 1.0/sqrt(1.0-vxt*vxt-vyt*vyt-vzt*vzt);
 double ut [4] = {gammat,gammat*vxt,gammat*vyt,gammat*vzt};
 double gammaf = 1.0/sqrt(1.0-vxf*vxf-vyf*vyf-vzf*vzf);
 double uf [4] = {gammaf,gammaf*vxf,gammaf*vyf,gammaf*vzf};

 // calculation of the energy-momentum tensor
 for (int i=0; i<4; i++){
  for (int j=0; j<4; j++){
   T[i][j] = (ep + pp) * up[i] * up[j] - pp * gmunu[i][j]
    + (et + pt) * ut[i] * ut[j] - pt * gmunu[i][j]
    + (ef + pf) * uf[i] * uf[j] - pf * gmunu[i][j];
  }
 }
}

void MultiHydro::getDiagonalizedEnergyDensity()
{
 double Q_p[7], Q_f[7], Q_t[7];
 double Ttemp[4][4];
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    c_p->getQ(Q_p);
    c_f->getQ(Q_f);
    c_t->getQ(Q_t);
    for (int i = 0; i < 7; i++) {
     Q_p[i] = Q_p[i]/h_p->getTau();
     Q_t[i] = Q_t[i]/h_t->getTau();
     Q_f[i] = Q_f[i]/h_f->getTau();
    }
    getEnergyMomentumTensor(Ttemp, Q_p, Q_f, Q_t);

    // calculation of the energy-momentum tensor
    TMatrixD T(4,4);
    for (int i=0; i<4; i++)
     for (int j=0; j<4; j++){
      T[i][j] = Ttemp[i][j]*gmunu[j][j];
    }
    if (T[0][0] == 0 && T[1][1] == 0 && T[2][2] == 0 && T[3][3] == 0)
     MHeps[ix][iy][iz] = 0;
    else
    {
     // diagonalization of the energy-momentum tensor
     TMatrixDEigen Te(T);
     TMatrixD eigenValues = Te.GetEigenValues();
     TMatrixD eigenVectors = Te.GetEigenVectors();

     double energyDensity;
     TVectorD v(4);
     for (int i=0; i<4; i++) {
      double vmuvmu = 0;
      energyDensity = eigenValues[i][i];
      v = TMatrixDColumn(eigenVectors,i);
      for (int j=0; j<4; j++) {
       vmuvmu += v[j]*v[j]*gmunu[j][j];
      }
      if (vmuvmu > 0 && energyDensity >= 0) {
       break;
      }
      else if (i == 3) {
       cout << "Multihydro: None of the eigenvectors is time-like, ";
       cout << "using largest eigenvalue for energy density." << endl;
       energyDensity = eigenValues[0][0];
       v = TMatrixDColumn(eigenVectors,0);
       break;
      }
     }

     // save computed energy density into private field
     MHeps[ix][iy][iz] = energyDensity;
    }
   }
}

void MultiHydro::getMaxEnergyDensity()
{
 double Q_p[7], Q_f[7], Q_t[7];
 double Ttemp[4][4];
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    c_p->getQ(Q_p);
    c_f->getQ(Q_f);
    c_t->getQ(Q_t);
    for (int i = 0; i < 7; i++) {
     Q_p[i] = Q_p[i]/h_p->getTau();
     Q_t[i] = Q_t[i]/h_t->getTau();
     Q_f[i] = Q_f[i]/h_f->getTau();
    }
    double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
    double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
    double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
    transformPV(eos, Q_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
    transformPV(eos, Q_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
    transformPV(eos, Q_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);

    MHeps[ix][iy][iz] = max(ep, max(et, ef));
   }
}

void MultiHydro::getSumEnergyDensity()
{
 double Q_p[7], Q_f[7], Q_t[7];
 double Ttemp[4][4];
 for (int iy = 0; iy < f_p->getNY(); iy++)
  for (int iz = 0; iz < f_p->getNZ(); iz++)
   for (int ix = 0; ix < f_p->getNX(); ix++) {
    Cell *c_p = f_p->getCell(ix, iy, iz);
    Cell *c_t = f_t->getCell(ix, iy, iz);
    Cell *c_f = f_f->getCell(ix, iy, iz);
    c_p->getQ(Q_p);
    c_f->getQ(Q_f);
    c_t->getQ(Q_t);
    for (int i = 0; i < 7; i++) {
     Q_p[i] = Q_p[i]/h_p->getTau();
     Q_t[i] = Q_t[i]/h_t->getTau();
     Q_f[i] = Q_f[i]/h_f->getTau();
    }
    double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
    double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
    double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
    transformPV(eos, Q_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
    transformPV(eos, Q_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
    transformPV(eos, Q_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);

    MHeps[ix][iy][iz] = ep + et + ef;
   }
}

void MultiHydro::getEnergyDensity()
{
 // In this function we decide which energy density will be used
 // to find the freezeout hypersurface.
 // - use getMaxEnergyDensity() to use maximum energy of the three fluids
 // - use getDiagonalizedEnergyDensity() to use energy density obtained
 //   by diagonalization of sum of the three energy-momentum tensors
 // getMaxEnergyDensity();
 getDiagonalizedEnergyDensity();
 // getSumEnergyDensity();
}

void MultiHydro::updateEnergyDensity()
{
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    MHepsPrev[ix][iy][iz] = MHeps[ix][iy][iz];
    MHeps[ix][iy][iz] = 0.0;
   }
}

void MultiHydro::outputEnergyDensity()
{
 ofstream fenergy;
 fenergy.open("energy.dat");
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    fenergy << h_p->getTau() << "\t" << f_p->getX(ix) << "\t" << f_p->getY(iy) << "\t" << f_p->getZ(iz) << "\t" << MHeps[ix][iy][iz] << endl;
   }
 fenergy.close();
}

void MultiHydro::resizeMHeps()
{
 double temp[nx][ny][nz];
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    if (2*(ix-(nx-1)/2)+(nx-1)/2 >= 0 && 2*(ix-(nx-1)/2)+(nx-1)/2 < nx && 2*(iy-(ny-1)/2)+(ny-1)/2 >=0 && 2*(iy-(ny-1)/2)+(ny-1)/2 < ny) {
     temp[ix][iy][iz] = MHeps[2*(ix-(nx-1)/2)+(nx-1)/2][2*(iy-(ny-1)/2)+(ny-1)/2][iz];
    } else {
     temp[ix][iy][iz] = 0.;
    }
   }
 for (int ix = 0; ix < nx; ix++)
  for (int iy = 0; iy < ny; iy++)
   for (int iz = 0; iz < nz; iz++) {
    MHeps[ix][iy][iz] = temp[ix][iy][iz];
   }
}

int MultiHydro::findFreezeout(EoS* eosH)
{
 updateEnergyDensity();
 getEnergyDensity();
 swap(eos, eosH);
 int nelements = 0;
 int ne_pos = 0;
 double E=0., Efull = 0.;

 // allocating corner points for Cornelius
 double ****ccube = new double ***[2];
 for (int i1 = 0; i1 < 2; i1++) {
  ccube[i1] = new double **[2];
  for (int i2 = 0; i2 < 2; i2++) {
   ccube[i1][i2] = new double *[2];
   for (int i3 = 0; i3 < 2; i3++) {
    ccube[i1][i2][i3] = new double[2];
   }
  }
 }

 for (int ix = 2; ix < nx - 2; ix++)
  for (int iy = 2; iy < ny - 2; iy++)
   for (int iz = 2; iz < nz - 2; iz++) {
    double QCube_p[2][2][2][2][7], QCube_f[2][2][2][2][7], QCube_t[2][2][2][2][7];
    // array for storing full energy-momentum tensor of all three fluids at corners
    double piSquare_p[2][2][2][10], PiSquare_p[2][2][2];
    double piSquare_t[2][2][2][10], PiSquare_t[2][2][2];
    double piSquare_f[2][2][2][10], PiSquare_f[2][2][2];

    // fill all corner cell with energy-momentum tensor
    for (int jx = 0; jx < 2; jx++)
     for (int jy = 0; jy < 2; jy++)
      for (int jz = 0; jz < 2; jz++) {
       ccube[0][jx][jy][jz] = MHepsPrev[ix + jx][iy + jy][iz + jz];
       ccube[1][jx][jy][jz] = MHeps[ix + jx][iy + jy][iz + jz];
       Cell *cc_p = f_p->getCell(ix + jx, iy + jy, iz + jz);
       Cell *cc_t = f_t->getCell(ix + jx, iy + jy, iz + jz);
       Cell *cc_f = f_f->getCell(ix + jx, iy + jy, iz + jz);
       double Qc_p[7], Qc_f[7], Qc_t[7];
       cc_p->getQ(QCube_p[1][jx][jy][jz]);
       cc_f->getQ(QCube_f[1][jx][jy][jz]);
       cc_t->getQ(QCube_t[1][jx][jy][jz]);
       for (int i = 0; i < 7; i++) {
        QCube_p[1][jx][jy][jz][i] = QCube_p[1][jx][jy][jz][i]/h_p->getTau();
        QCube_t[1][jx][jy][jz][i] = QCube_t[1][jx][jy][jz][i]/h_t->getTau();
        QCube_f[1][jx][jy][jz][i] = QCube_f[1][jx][jy][jz][i]/h_f->getTau();
       }
       cc_p->getQprev(QCube_p[0][jx][jy][jz]);
       cc_f->getQprev(QCube_f[0][jx][jy][jz]);
       cc_t->getQprev(QCube_t[0][jx][jy][jz]);
       for (int i = 0; i < 7; i++) {
        QCube_p[0][jx][jy][jz][i] = QCube_p[0][jx][jy][jz][i]/h_p->getTau();
        QCube_t[0][jx][jy][jz][i] = QCube_t[0][jx][jy][jz][i]/h_t->getTau();
        QCube_f[0][jx][jy][jz][i] = QCube_f[0][jx][jy][jz][i]/h_f->getTau();
       }
       for (int ii = 0; ii < 4; ii++)
        for (int jj = 0; jj <= ii; jj++) {
         piSquare_p[jx][jy][jz][index44(ii, jj)] = cc_p->getpi(ii, jj);
         piSquare_t[jx][jy][jz][index44(ii, jj)] = cc_t->getpi(ii, jj);
         piSquare_f[jx][jy][jz][index44(ii, jj)] = cc_f->getpi(ii, jj);
       }
       PiSquare_p[jx][jy][jz] = cc_p->getPi();
       PiSquare_t[jx][jy][jz] = cc_t->getPi();
       PiSquare_f[jx][jy][jz] = cc_f->getPi();
    }

    // cornelius
    cornelius->find_surface_4d(ccube);
    const int Nsegm = cornelius->get_Nelements();
    for (int isegm = 0; isegm < Nsegm; isegm++) {
     nelements++;

     // interpolation procedure
     double vxC = 0., vyC = 0., vzC = 0., TC = 0., mubC = 0., muqC = 0.,
            musC = 0., PiC_p = 0., PiC_t = 0., PiC_f = 0., nbC = 0., nqC = 0.;
     double piC_p[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     double piC_t[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     double piC_f[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     double QC_p[7] = {0., 0., 0., 0., 0., 0., 0.};
     double QC_t[7] = {0., 0., 0., 0., 0., 0., 0.};
     double QC_f[7] = {0., 0., 0., 0., 0., 0., 0.};
     for (int ii = 0; ii < 10; ii++) {
      piC_p[ii] = 0.0;
      piC_t[ii] = 0.0;
      piC_f[ii] = 0.0;
     }
     double wCenT[2] = {1. - cornelius->get_centroid_elem(isegm, 0) / h_p->getDtau(),
                        cornelius->get_centroid_elem(isegm, 0) / h_p->getDtau()};
     double wCenX[2] = {1. - cornelius->get_centroid_elem(isegm, 1) / dx,
                        cornelius->get_centroid_elem(isegm, 1) / dx};
     double wCenY[2] = {1. - cornelius->get_centroid_elem(isegm, 2) / dy,
                        cornelius->get_centroid_elem(isegm, 2) / dy};
     double wCenZ[2] = {1. - cornelius->get_centroid_elem(isegm, 3) / dz,
                        cornelius->get_centroid_elem(isegm, 3) / dz};

     for (int jt = 0; jt < 2; jt++)
      for (int jx = 0; jx < 2; jx++)
       for (int jy = 0; jy < 2; jy++)
        for (int jz = 0; jz < 2; jz++) {
         for (int i = 0; i < 4; i++)
         for (int i = 0; i < 7; i++) {
          QC_p[i] += QCube_p[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
          QC_t[i] += QCube_t[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
          QC_f[i] += QCube_f[jt][jx][jy][jz][i] * wCenT[jt] * wCenX[jx] *
                   wCenY[jy] * wCenZ[jz];
         }
     }

     /*for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
       TmunuC[i][j] = TmunuC[i][j] / (h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0));
     for (int i = 0; i < 7; i++) {
      QC_p[i] = QC_p[i] / (h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0));
      QC_t[i] = QC_t[i] / (h_t->getTau() - h_t->getDtau() + cornelius->get_centroid_elem(isegm, 0));
      QC_f[i] = QC_f[i] / (h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0));
     }*/
     double _ns = 0.0;
     double ep, pp, nbp, nqp, nsp, vxp, vyp, vzp;
     double et, pt, nbt, nqt, nst, vxt, vyt, vzt;
     double ef, pf, nbf, nqf, nsf, vxf, vyf, vzf;
     transformPV(eos, QC_p, ep, pp, nbp, nqp, nsp, vxp, vyp, vzp);
     transformPV(eos, QC_t, et, pt, nbt, nqt, nst, vxt, vyt, vzt);
     transformPV(eos, QC_f, ef, pf, nbf, nqf, nsf, vxf, vyf, vzf);

     //transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC);
     nbC = nbp + nbt + nbf;
     nqC = nqp + nqt + nqf;
     _ns = nsp + nst + nsf;
     double TCp, mubCp, muqCp, musCp, pCp;
     double TCt, mubCt, muqCt, musCt, pCt;
     double TCf, mubCf, muqCf, musCf, pCf;
     eos->eos(ep, nbp, nqp, nsp, TCp, mubCp, muqCp, musCp, pCp);
     eos->eos(et, nbt, nqt, nst, TCt, mubCt, muqCt, musCt, pCt);
     eos->eos(ef, nbf, nqf, nsf, TCf, mubCf, muqCf, musCf, pCf);
     /*if (TC > 0.4 || abs(mubC) > 0.85) {
      cout << "#### Error (multifluid surface): high T/mu_b (T=" << TC << "/mu_b=" << mubC << ") ####\n";
     }*/
     for (int jx = 0; jx < 2; jx++)
      for (int jy = 0; jy < 2; jy++)
       for (int jz = 0; jz < 2; jz++) {
        for (int ii = 0; ii < 10; ii++) {
         piC_p[ii] += piSquare_p[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
         piC_t[ii] += piSquare_t[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
         piC_f[ii] += piSquare_f[jx][jy][jz][ii] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
        }
        PiC_p += PiSquare_p[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
        PiC_t += PiSquare_t[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
        PiC_f += PiSquare_f[jx][jy][jz] * wCenX[jx] * wCenY[jy] * wCenZ[jz];
       }
     double v2C = vxC * vxC + vyC * vyC + vzC * vzC;
     if (v2C > 1.) {
      vxC *= sqrt(0.99 / v2C);
      vyC *= sqrt(0.99 / v2C);
      vzC *= sqrt(0.99 / v2C);
      v2C = 0.99;
     }
     double etaC = f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3);
     transformToLab(etaC, vxp, vyp, vzp);
     transformToLab(etaC, vxt, vyt, vzt);
     transformToLab(etaC, vxf, vyf, vzf);
     double gammaC_p = 1. / sqrt(1. - vxp * vxp - vyp * vyp - vzp * vzp);
     double gammaC_t = 1. / sqrt(1. - vxt * vxt - vyt * vyt - vzt * vzt);
     double gammaC_f = 1. / sqrt(1. - vxf * vxf - vyf * vyf - vzf * vzf);
     double uC_p[4] = {gammaC_p, gammaC_p * vxp, gammaC_p * vyp, gammaC_p * vzp};
     double uC_t[4] = {gammaC_t, gammaC_t * vxt, gammaC_t * vyt, gammaC_t * vzt};
     double uC_f[4] = {gammaC_f, gammaC_f * vxf, gammaC_f * vyf, gammaC_f * vzf};
     const double tauC = h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0);
     double dsigma[4], dsds;
     // ---- transform dsigma to lab.frame :
     const double ch = cosh(etaC);
     const double sh = sinh(etaC);
     dsigma[0] = tauC * (ch * cornelius->get_normal_elem(0, 0) -
                         sh / tauC * cornelius->get_normal_elem(0, 3));
     dsigma[3] = tauC * (-sh * cornelius->get_normal_elem(0, 0) +
                         ch / tauC * cornelius->get_normal_elem(0, 3));
     dsigma[1] = tauC * cornelius->get_normal_elem(0, 1);
     dsigma[2] = tauC * cornelius->get_normal_elem(0, 2);
     dsds = dsigma[0]*dsigma[0] - dsigma[1]*dsigma[1] - dsigma[2]*dsigma[2] - dsigma[3]*dsigma[3];
     double dVEff_p = 0.0, dVEff_t = 0.0, dVEff_f = 0.0, dVEff = 0.0;
     for (int ii = 0; ii < 4; ii++) {
      dVEff_p += dsigma[ii] * uC_p[ii];  // normalize for Delta eta=1
      dVEff_t += dsigma[ii] * uC_t[ii];
      dVEff_f += dsigma[ii] * uC_f[ii];
     }
     if (dVEff_p > 0) ne_pos++;
     vEff_p += dVEff_p;
     vEff_t += dVEff_t;
     vEff_f += dVEff_f;

     double picart_p[10];
     double picart_t[10];
     double picart_f[10];
#ifdef OUTPI
     // projectile fluid
     /*pi00*/ picart_p[index44(0, 0)] = ch * ch * piC_p[index44(0, 0)] +
                                      2. * ch * sh * piC_p[index44(0, 3)] +
                                      sh * sh * piC_p[index44(3, 3)];
     /*pi01*/ picart_p[index44(0, 1)] =
         ch * piC_p[index44(0, 1)] + sh * piC_p[index44(3, 1)];
     /*pi02*/ picart_p[index44(0, 2)] =
         ch * piC_p[index44(0, 2)] + sh * piC_p[index44(3, 2)];
     /*pi03*/ picart_p[index44(0, 3)] =
         ch * sh * (piC_p[index44(0, 0)] + piC_p[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC_p[index44(0, 3)];
     /*pi11*/ picart_p[index44(1, 1)] = piC_p[index44(1, 1)];
     /*pi12*/ picart_p[index44(1, 2)] = piC_p[index44(1, 2)];
     /*pi13*/ picart_p[index44(1, 3)] =
         sh * piC_p[index44(0, 1)] + ch * piC_p[index44(3, 1)];
     /*pi22*/ picart_p[index44(2, 2)] = piC_p[index44(2, 2)];
     /*pi23*/ picart_p[index44(2, 3)] =
         sh * piC_p[index44(0, 2)] + ch * piC_p[index44(3, 2)];
     /*pi33*/ picart_p[index44(3, 3)] = sh * sh * piC_p[index44(0, 0)] +
                                      ch * ch * piC_p[index44(3, 3)] +
                                      2. * sh * ch * piC_p[index44(0, 3)];
     // target fluid
     /*pi00*/ picart_t[index44(0, 0)] = ch * ch * piC_t[index44(0, 0)] +
                                      2. * ch * sh * piC_t[index44(0, 3)] +
                                      sh * sh * piC_t[index44(3, 3)];
     /*pi01*/ picart_t[index44(0, 1)] =
         ch * piC_t[index44(0, 1)] + sh * piC_t[index44(3, 1)];
     /*pi02*/ picart_t[index44(0, 2)] =
         ch * piC_t[index44(0, 2)] + sh * piC_t[index44(3, 2)];
     /*pi03*/ picart_t[index44(0, 3)] =
         ch * sh * (piC_t[index44(0, 0)] + piC_t[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC_t[index44(0, 3)];
     /*pi11*/ picart_t[index44(1, 1)] = piC_t[index44(1, 1)];
     /*pi12*/ picart_t[index44(1, 2)] = piC_t[index44(1, 2)];
     /*pi13*/ picart_t[index44(1, 3)] =
         sh * piC_t[index44(0, 1)] + ch * piC_t[index44(3, 1)];
     /*pi22*/ picart_t[index44(2, 2)] = piC_t[index44(2, 2)];
     /*pi23*/ picart_t[index44(2, 3)] =
         sh * piC_t[index44(0, 2)] + ch * piC_t[index44(3, 2)];
     /*pi33*/ picart_t[index44(3, 3)] = sh * sh * piC_t[index44(0, 0)] +
                                      ch * ch * piC_t[index44(3, 3)] +
                                      2. * sh * ch * piC_t[index44(0, 3)];
     // fireball fluid
     /*pi00*/ picart_f[index44(0, 0)] = ch * ch * piC_f[index44(0, 0)] +
                                      2. * ch * sh * piC_f[index44(0, 3)] +
                                      sh * sh * piC_f[index44(3, 3)];
     /*pi01*/ picart_f[index44(0, 1)] =
         ch * piC_f[index44(0, 1)] + sh * piC_f[index44(3, 1)];
     /*pi02*/ picart_f[index44(0, 2)] =
         ch * piC_f[index44(0, 2)] + sh * piC_f[index44(3, 2)];
     /*pi03*/ picart_f[index44(0, 3)] =
         ch * sh * (piC_f[index44(0, 0)] + piC_f[index44(3, 3)]) +
         (ch * ch + sh * sh) * piC_f[index44(0, 3)];
     /*pi11*/ picart_f[index44(1, 1)] = piC_f[index44(1, 1)];
     /*pi12*/ picart_f[index44(1, 2)] = piC_f[index44(1, 2)];
     /*pi13*/ picart_f[index44(1, 3)] =
         sh * piC_f[index44(0, 1)] + ch * piC_f[index44(3, 1)];
     /*pi22*/ picart_f[index44(2, 2)] = piC_f[index44(2, 2)];
     /*pi23*/ picart_f[index44(2, 3)] =
         sh * piC_f[index44(0, 2)] + ch * piC_f[index44(3, 2)];
     /*pi33*/ picart_f[index44(3, 3)] = sh * sh * piC_f[index44(0, 0)] +
                                      ch * ch * piC_f[index44(3, 3)] +
                                      2. * sh * ch * piC_f[index44(0, 3)];
#endif

    double PiC=0.0;
    double picart[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

     double dEtotSurf[3] = {0., 0., 0.};
     dEtotSurf[0] = (ep + pCp) * uC_p[0] * dVEff_p - pCp * dsigma[0]; // projectile
     dEtotSurf[1] = (et + pCt) * uC_t[0] * dVEff_t - pCt * dsigma[0]; // target
     dEtotSurf[2] = (ef + pCf) * uC_f[0] * dVEff_f - pCf * dsigma[0]; // fireball
     EtotSurf[0] += dEtotSurf[0];
     EtotSurf[1] += dEtotSurf[1];
     EtotSurf[2] += dEtotSurf[2];
     if (dEtotSurf[0] > 0) EtotSurf_positive[0] += dEtotSurf[0];
     else EtotSurf_negative[0] += dEtotSurf[0];
     if (dEtotSurf[1] > 0) EtotSurf_positive[1] += dEtotSurf[1];
     else EtotSurf_negative[1] += dEtotSurf[1];
     if (dEtotSurf[2] > 0) EtotSurf_positive[2] += dEtotSurf[2];
     else EtotSurf_negative[2] += dEtotSurf[2];

     // exclude segments which fulfills dSigma_0 < 0 & dSigma^2 > 0 - those cells have energy flow into >
     if (dsds > 0) {
      if (dsigma[0] > 0) {
       printFreezeout(
        fmhfreeze_p,
        h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC_p, TCp, mubCp, muqCp, musCp, picart_p, PiC_p, dVEff_p
       );
       printFreezeout(
        fmhfreeze_t,
        h_t->getTau() - h_t->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_t->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_t->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_t->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC_t, TCt, mubCt, muqCt, musCt, picart_t, PiC_t, dVEff_t
       );
       printFreezeout(
        fmhfreeze_f,
        h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0),
        f_f->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
        f_f->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
        f_f->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
        dsigma, uC_f, TCf, mubCf, muqCf, musCf, picart_f, PiC_f, dVEff_f
       );
      }
     } else {
      if (dEtotSurf[0] > 0 && dVEff_p > 0) printFreezeout(
       fmhfreeze_p,
       h_p->getTau() - h_p->getDtau() + cornelius->get_centroid_elem(isegm, 0),
       f_p->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
       f_p->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
       f_p->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
       dsigma, uC_p, TCp, mubCp, muqCp, musCp, picart_p, PiC_p, dVEff_p
      );
      if (dEtotSurf[1] > 0 && dVEff_t > 0) printFreezeout(
       fmhfreeze_t,
       h_t->getTau() - h_t->getDtau() + cornelius->get_centroid_elem(isegm, 0),
       f_t->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
       f_t->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
       f_t->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
       dsigma, uC_t, TCt, mubCt, muqCt, musCt, picart_t, PiC_t, dVEff_t
      );
      if (dEtotSurf[2] > 0 && dVEff_f > 0) printFreezeout(
       fmhfreeze_f,
       h_f->getTau() - h_f->getDtau() + cornelius->get_centroid_elem(isegm, 0),
       f_f->getX(ix) + cornelius->get_centroid_elem(isegm, 1),
       f_f->getY(iy) + cornelius->get_centroid_elem(isegm, 2),
       f_f->getZ(iz) + cornelius->get_centroid_elem(isegm, 3),
       dsigma, uC_f, TCf, mubCf, muqCf, musCf, picart_f, PiC_f, dVEff_f
      );
     }

    }
 }

 cout << setw(10) << h_p->getTau() << setw(10) << nelements << "\t" << ne_pos << "\t"
      << EtotSurf[0] << "\t" << EtotSurf_positive[0] << "\t" << EtotSurf_negative[0] << "\t"
      << EtotSurf[1] << "\t" << EtotSurf_positive[1] << "\t" << EtotSurf_negative[1] << "\t"
      << EtotSurf[2] << "\t" << EtotSurf_positive[2] << "\t" << EtotSurf_negative[2] << endl;
 swap(eos, eosH); // get back to the hydrodynamic EoS
 for (int i1 = 0; i1 < 2; i1++) {
  for (int i2 = 0; i2 < 2; i2++) {
   for (int i3 = 0; i3 < 2; i3++) {
    delete[] ccube[i1][i2][i3];
   }
   delete[] ccube[i1][i2];
  }
  delete[] ccube[i1];
 }
 delete[] ccube;
 if (nelements == 0 && h_p->getTau() > 5 + tau0)
  return 1;   // particlization surface ended - return 1 for the evolution to stop
  //return 0;
 else
  return 0;   // return 0 for the evolution to continue
}

void MultiHydro::printFreezeout(std::ofstream &fout, double t, double x, double y, double z, double dsigma[4], double uC[4], double TC, double mub, double muq, double mus, double picart[10], double PiC, double dVEff)
{
 fout.precision(15);
 fout << setw(24) << t << setw(24) << x << setw(24) << y << setw(24) << z;
 for (int ii = 0; ii < 4; ii++) {
  fout << setw(24) << dsigma[ii];
 }
 for (int ii = 0; ii < 4; ii++) {
  fout << setw(24) << uC[ii];
 }
 fout << setw(24) << TC << setw(24) << mub << setw(24) << muq << setw(24) << mus;
#ifdef OUTPI
 for (int ii = 0; ii < 10; ii++) {
  fout << setw(24) << picart[ii];
 }
 fout << setw(24) << PiC;
#else
 fout << setw(24) << dVEff;
#endif
 fout << endl;
}


double MultiHydro::totalScatRate(double px, double T, double mu, double u[4])
{
 return 2*calculateScatRates(px, T, mu, u, 1) + 3*calculateScatRates(px, T, mu, u, 0);
}


double MultiHydro::calculateScatRates(double px, double T, double mu, double u[4], int particle)
{
 /*
   This function calculates scattering rates of nucleon on nucleons in fluids
 */
 double m1 = 0.938; // scattered particle
 double m2 = particle == 0 ? 0.1396 : 0.938; // scatterers in cloud
 TLorentzVector pLV;
 pLV.SetPxPyPzE(px, 0, 0, sqrt(px*px+m1*m1));
 pLV.Boost(-u[1]/u[0], -u[2]/u[0], -u[3]/u[0]);
 double p = sqrt(pow(pLV.Px(),2) + pow(pLV.Py(),2) + pow(pLV.Pz(),2));
 double gamma = u[0];
 double R;
 double Epi = sqrt(m1*m1+p*p);
 double scrit = m1*m1 + m2*m2 + 2*Epi*m2;
 double ds = 0.1;
 double smin = pow(m1 + m2,2);
 double smax = 10;
 double g = 2;
 for (double s = smin + ds/2; s < smax; s += ds)
 {
  //double Ekin = s/(2.0*mN) - 2.0*mN;
  //double pLab = sqrt(s*(s-4*mN*mN)) / (2.0*mN);
  double sa = pow(m1+m2,2);
  double sb = pow(m1-m2,2);
  double sprime = s - m1*m1 - m2*m2;
  double ebar, deltaE;
  if (s <= scrit) {
   ebar = 0.5*((Epi*sprime+p*sqrt(sprime*sprime-4*m1*m1*m2*m2))/(2*m1*m1)+2*Epi*m2*m2/sprime);
   deltaE = 0.5*((Epi*sprime+p*sqrt(sprime*sprime-4*m1*m1*m2*m2))/(2*m1*m1)-2*Epi*m2*m2/sprime);
  } else {
   ebar = Epi*sprime/(2*m1*m1);
   deltaE = p*sqrt(sprime*sprime-4*m1*m1*m2*m2)/(2*m1*m1);
  }
  double dR = sqrt(s-sa)*sqrt(s-sb)*exp(-ebar/T)*sinh(deltaE/T)/pow(0.197,3);
  // case of proton-proton scattering
  if (particle == 1) dR *= 0.1 * pp_total(s);
  // case of pion-proton scattering
  if (particle == 0) dR *= xsect->piN(s);
  //cout << s << " " << exp(ebar/T) << " " << exp(-sqrt(m1*m1+p*p)*(s - m1*m1 - m2*m2)/(2*T*m1*m1)) << " " << dR << endl;
  if (dR > 1e-20 && dR < 1e20) R += dR * ds;
 }
 R *= g*T*exp(mu/T)/(8*pow(M_PI,2)*p*Epi);
 R *= gamma;
 //exit(1);
 return R;
}

double MultiHydro::pp_total(double mandelstam_s) {
  double mN = 0.938;
  const double p_lab = sqrt(mandelstam_s*(mandelstam_s-4*mN*mN)) / (2.0*mN);
  if (p_lab < 0.4) {
    return 34 * pow(p_lab / 0.4, -2.104);
  } else if (p_lab < 0.8) {
    return 23.5 + 1000 * pow(p_lab - 0.7, 4);
  } else if (p_lab < 1.5) {
    return 23.5 + 24.6 / (1 + exp(-(p_lab - 1.2) / 0.1));
  } else if (p_lab < 5.0) {
    return 41 + 60 * (p_lab - 0.9) * exp(-1.2 * p_lab);
  } else {
    const auto logp = log(p_lab);
    return 48.0 + 0.522 * logp * logp - 4.51 * logp;
  }
}

double MultiHydro::Fermi(double nb) {
 return pow(0.197,2)*pow(6*nb/M_PI,2/3)*M_PI*M_PI/(2*0.938);
}

double MultiHydro::getLocalEnergyDensity(double x, double y, double eta) {
// this function returns effective energy density, pre-computed in the MHeps[][][] array,
// at an arbitrary position using tri-linear interpolation over the table(array)
 int ix = (int)((x - f_p->getX(0))/dx);
 int iy = (int)((y - f_p->getY(0))/dy);
 int iz = (int)((eta - f_p->getZ(0))/dz);
 if(ix<0 || ix>nx-1 || iy<0 || iy>ny-1 || iz<0  || iz>nz-1) {
  cout << "MultiHydro::getLocalEnergyDensity: out of bounds\n";
  exit(77);
 }
 double xm = (x - f_p->getX(0) - ix*dx) / dx;
 double ym = (y - f_p->getY(0) - iy*dy) / dy;
 double zm = (eta - f_p->getZ(0) - iz*dz) / dz;
 double wx [2] = {1. - xm, xm};
 double wy [2] = {1. - ym, ym};
 double wz [2] = {1. - zm, zm};
 double interpolatedEps = 0.;
 for(int i=0; i<2; i++)
 for(int j=0; j<2; j++)
 for(int k=0; k<2; k++)
  interpolatedEps += wx[i]*wy[j]*wz[k]*MHeps[ix+i][iy+j][iz+k];
 return interpolatedEps;
}

void MultiHydro::evolveSpectators(void) {
 for(int iev=0; iev<nucleons.size(); iev++) {
 for(int i=0; i<nucleons[iev].size(); i++) {
  nucleons[iev][i].eta = nucleons[iev][i].rap +
   asinh(h_f->getTau()/(h_f->getTau()+dtau)*sinh(nucleons[iev][i].eta - nucleons[iev][i].rap));
 }
 }
 for(int iev=0; iev<nucleons.size(); iev++) {
  for (auto it = nucleons[iev].begin(); it != nucleons[iev].end();) {
   if(getLocalEnergyDensity((*it).x, (*it).y, (*it).eta)>ecrit) {
    it = nucleons[iev].erase(it);
   } else
   it++;
  }
 }
}

void MultiHydro::printSpectators(std::ofstream &fout) {
 for(int iev=0; iev<nucleons.size(); iev++) {
  fout << "#event " << iev << endl;
  for(int i=0; i<nucleons[iev].size(); i++) {
  // id offset for the output is big enough so that there is no overlap
  // with ids of hadrons, generated by smash-hadron-sampler
  const int id = i + 50000;
  Nucleon &n = nucleons[iev][i];
  const int pdg = n.charge==1 ? 2212 : 2112;
  const double mass = 0.938;
  fout << tau0*cosh(n.getEta0()) << " " << n.x
   << " " << n.y << " " << tau0*sinh(n.getEta0())
   << " " << mass << " " << mass*cosh(n.rap)
   << " " << 0 << " " << 0 << " " << mass*sinh(n.rap)
   << " " << pdg << " " << id << " " << n.charge << endl;
  }
 }
 fout << "#event end" << endl;
}
