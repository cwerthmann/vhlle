#include "xsect.h"
#include "EfI.h"
#include "inc.h"
#include <iostream>

using namespace std;

EfIntegrand::EfIntegrand(CrossSections *_xsect, double _Tf, double _mubf, double _valphatilde){
    xsect=_xsect;
    Tf=_Tf;
    mubf=_mubf;
    valphatilde=_valphatilde;
    gammaalphatilde=1.0/sqrt(1.0-_valphatilde*_valphatilde);
}

/*double EfIntegrand::EvalNpi(double s){//1
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    if(valphatilde<0.001){
        return 9.0*M_PI/mN/mN*sigmaNpi*(seff*seff-4.0*mN*mN*mpi*mpi)/(exp(seff*gammaalphatilde/2.0/mN/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
    }
    return -9.0*M_PI/mN/mN*sigmaNpi*moller*(2.0*moller-Tf*mN/gammaalphatilde/valphatilde*log((exp(p0max/Tf)-1.0)/(exp(p0min/Tf)-1.0)))*gevtofm*gevtofm*gevtofm;
}*/

/*double EfIntegrand::EvalNpi(double s){//4&10
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    if(valphatilde<0.001){
        return 9.0*M_PI/mN*sigmaNpi*moller*seff/mN*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff)/(exp((p0max+p0min)/2/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
    }
    return 9.0*M_PI/mN*Tf/gammaalphatilde/valphatilde*sigmaNpi*moller*(log((exp(-p0max/Tf)-1.0)/(exp(-p0min/Tf)-1.0)))*gevtofm*gevtofm*gevtofm;
}*/

/*double EfIntegrand::EvalNpi(double p){//12
    double seff=2.0*p*mN*gammaalphatilde;
    double sigmaNpi=xsect-> piN(sqrt(seff+mN*mN+mpi*mpi));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    return 36.0*M_PI*sqrt(p*p-mpi*mpi)*sigmaNpi*moller/(exp(p/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
}*/

/*double EfIntegrand::EvalNpi(double p){//14
    double seff=2.0*p*mN*gammaalphatilde;
    double sigmaNpi=xsect-> piN(sqrt(seff+mN*mN+mpi*mpi));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    return 1.5/M_PI/M_PI*sqrt(p*p-mpi*mpi)*sigmaNpi*moller/(exp(p/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
}*/

double EfIntegrand::EvalNpi(double p){//17
    double seff=2.0*p*mN*gammaalphatilde;
    double sigmaNpi=xsect-> piN(sqrt(seff+mN*mN+mpi*mpi));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    return 1.5/M_PI/M_PI*p*sqrt(p*p-mpi*mpi)*sigmaNpi*moller/(exp(p/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
}



/*double EfIntegrand::EvalNpi(double s){//8
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    if(valphatilde<0.001){
        return 9.0*M_PI/mN/mN*sigmaNpi*(seff*seff-4.0*mN*mN*mpi*mpi)/(exp(seff*gammaalphatilde/2.0/mN/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
    }
    return 9.0*M_PI/mN*Tf/gammaalphatilde/valphatilde*sigmaNpi*moller*(log(1.0-exp(-p0max/Tf))-log(1.0-exp(-p0min/Tf)))*gevtofm*gevtofm*gevtofm;
}*/

/*double EfIntegrand::EvalNpi(double s){//5
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    if(valphatilde<0.001){
        return 9.0*M_PI/mN/mN*sigmaNpi*(seff*seff-4.0*mN*mN*mpi*mpi)/(exp(seff*gammaalphatilde/2.0/mN/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
    }
    return -4.5*M_PI/mN/gammaalphatilde/valphatilde*sigmaNpi*moller*(p0max-p0min-Tf*log((cosh(p0max/Tf)-1.0)/(cosh(p0min/Tf)-1.0)))*gevtofm*gevtofm*gevtofm;
}*/


/*double EfIntegrand::EvalNpi(double s){//7
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    if(valphatilde<0.001){
        return 9.0*M_PI/mN/mN*sigmaNpi*(seff*seff-4.0*mN*mN*mpi*mpi)/(exp(seff*gammaalphatilde/2.0/mN/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
    }
    return 9.0*M_PI/mN/gammaalphatilde/valphatilde*sigmaNpi*moller*(-p0max+p0min+Tf*log(1.0-exp(p0max/Tf)/(1.0-exp(p0min/Tf))))*gevtofm*gevtofm*gevtofm;
}*/

/*double EfIntegrand::EvalNpi(double s){//6&9
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    if(valphatilde<0.001){
        return 9.0*M_PI/mN/mN*sigmaNpi*(seff*seff-4.0*mN*mN*mpi*mpi)/(exp(seff*gammaalphatilde/2.0/mN/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
    }
    return 9.0*M_PI/mN/gammaalphatilde/valphatilde*sigmaNpi*moller*(p0max-p0min)/(exp((p0max+p0min)/2/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
}*/

/*double EfIntegrand::EvalNpi(double s){//3
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    if(valphatilde<0.01){
        return 9.0*M_PI/mN/mN*sigmaNpi*(seff*seff-4.0*mN*mN*mpi*mpi)/(exp(seff*gammaalphatilde/2.0/mN/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
    }
    return -9.0*M_PI*Tf/gammaalphatilde/valphatilde/mN*sigmaNpi*moller*log((exp(seff*gammaalphatilde/mN/Tf)-exp(p0max/Tf))/(exp(seff*gammaalphatilde/mN/Tf)-exp(p0min/Tf)))*gevtofm*gevtofm*gevtofm;
}*/

/*double EfIntegrand::EvalNpi(double s){//2
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    return 9.0*M_PI/mN/mN*sigmaNpi*(seff*seff-4.0*mN*mN*mpi*mpi)/(exp(seff*gammaalphatilde/2.0/mN/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
}*/

double EfIntegrand::EvalNN(double s){
    double sigmaNN;
    xsect-> NN(sqrt(s),sigmaNN);
    double seff=s-mN*mN-mN*mN;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mN*mN/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mN*mN/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mN*mN);
    return -4.0*M_PI/valphatilde/gammaalphatilde/mN*sigmaNN*moller*log((exp(-(p0max-mubf)/Tf)+1.0)/(exp(-(p0min-mubf)/Tf)+1.0))*gevtofm*gevtofm*gevtofm;
}

