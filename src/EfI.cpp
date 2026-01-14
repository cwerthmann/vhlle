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


double EfIntegrand::EvalNpi_to_f(double p){
    double seff=2.0*p*mN*gammaalphatilde;
    double sigmaNpi=xsect-> piN(sqrt(seff+mN*mN+mpi*mpi));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    return 1.5/M_PI/M_PI*sqrt(p*p-mpi*mpi)*sigmaNpi*moller/(exp(p/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
}

double EfIntegrand::EvalNpi_to_pt(double p){
    double seff=2.0*p*mN*gammaalphatilde;
    double sigmaNpi=xsect-> piN(sqrt(seff+mN*mN+mpi*mpi));
    double moller=0.5*sqrt(seff*seff-4.0*mN*mN*mpi*mpi);
    return 1.5/M_PI/M_PI*p/mN*sqrt(p*p-mpi*mpi)*sigmaNpi*moller/(exp(p/Tf)-1.0)*gevtofm*gevtofm*gevtofm;
}


