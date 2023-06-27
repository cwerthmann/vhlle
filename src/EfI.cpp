#include "xsect.h"
#include "EfI.h"
#include "inc.h"

using namespace std;

EfIntegrand::EfIntegrand(CrossSections *_xsect, double _Tf, double _valphatilde){
    xsect=_xsect;
    Tf=_Tf;
    valphatilde=_valphatilde;
    gammaalphatilde=1.0/std::sqrt(1.0-_valphatilde*_valphatilde);
}

double EfIntegrand::Eval(double s){
    double sigmaNpi=xsect-> piN(std::sqrt(s));
    double seff=s-mN*mN-mpi*mpi;
    double p0min=0.5*seff*gammaalphatilde/mN*(1.0-valphatilde*std::sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double p0max=0.5*seff*gammaalphatilde/mN*(1.0+valphatilde*std::sqrt(1.0-4.0*mN*mN*mpi*mpi/seff/seff));
    double moller=0.5*std::sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    return sigmaNpi*moller*Tf*std::log((std::exp(-p0min/Tf)-1.0)/(std::exp(-p0max/Tf)-1.0))*gevtofm*gevtofm*gevtofm;
}

