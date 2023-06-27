#include "xsect.h"
#include "EfI.h"
#include "inc.h"

using namespace std;

EfIntegrand::EfIntegrand(CrossSections *_xsect, double _Tf, double _mubf, double _valphatilde){
    xsect=_xsect;
    Tf=_Tf;
    mubf=_mubf;
    valphatilde=_valphatilde;
    gammaalphatilde=1.0/sqrt(1.0-_valphatilde*_valphatilde);
}

double EfIntegrand::EvalNpi(double s){
    double sigmaNpi=xsect-> piN(sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    return 9*sigmaNpi*moller*Tf*log((exp(-p0min/Tf)-1.0)/(exp(-p0max/Tf)-1.0))*gevtofm*gevtofm*gevtofm;
}

double EfIntegrand::EvalNN(double s){
    double sigmaNN;
    xsect-> NN(sqrt(s),sigmaNN);
    double seff=s-mN*mN-mN*mN;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*sqrt(1.0-4.0*mN*mN*mN*mN/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*sqrt(1.0-4.0*mN*mN*mN*mN/seff/seff));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mN*mN);
    return 4*sigmaNN*moller*Tf*log((exp(-(p0min-mubf)/Tf)+1.0)/(exp(-(p0max-mubf)/Tf)+1.0))*gevtofm*gevtofm*gevtofm;
}

