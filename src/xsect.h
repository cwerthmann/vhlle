#include <cmath>

class TGraph;

class CrossSections {
 TGraph *gSigmaT, *gSigmaE, *gSigmaPimp, *gSigmaPipp;
 double alpha=0.5, beta=0.5;
public:
 CrossSections(void);
 ~CrossSections(void);
 void Ivanov(double Ekin, double& sigmaT, double& sigmaE, double& sigmaP);
 void Ivanovbar(double sqrts, double& sigmaEbar, double& sigmaPbar, double& sigmaRbar);
 void NN(double Ekin, double& sigmaNN);
 double piN(double sqrts);
 inline double A(double yalpha){return yalpha+0.6*std::tanh(10*yalpha);}
 inline double B(double yalpha){return beta*yalpha+0.6*std::sinh(10*beta*yalpha)/std::cosh(10*yalpha);}
 inline double C(double yalpha){return (std::log(std::cosh(yalpha)/std::cosh(beta*yalpha)))/std::tanh(yalpha)*(1.0-6.0/std::cosh(10*yalpha))
        +(6*std::cosh(2*yalpha)-6*std::cosh(2*beta*yalpha)-3*std::cosh(4*yalpha)+3*std::cosh(4*beta*yalpha)+2*std::cosh(6*yalpha)-2*std::cosh(6*beta*yalpha)
          -1.5*std::cosh(8*yalpha)+1.5*std::cosh(8*beta*yalpha)+0.6*std::cosh(10*yalpha)-0.6*std::cosh(10*beta*yalpha))/std::cosh(10*yalpha)/std::tanh(yalpha);}
 inline double D(double yalpha){return (std::sinh(yalpha)-std::sinh(beta*yalpha)+(1.0/3.0*(std::sinh(9*yalpha)-std::sinh(9*beta*yalpha))+3.0/11.0*(std::sinh(11*yalpha)-std::sinh(11*beta*yalpha)))/std::cosh(10*yalpha))/std::cosh(yalpha);}
 inline double F(double yalpha){return (std::cosh(yalpha)-std::cosh(beta*yalpha)+(-1.0/3.0*(std::cosh(9*yalpha)-std::cosh(9*beta*yalpha))+3.0/11.0*(std::cosh(11*yalpha)-std::cosh(11*beta*yalpha)))/std::cosh(10*yalpha))/std::sinh(yalpha);}
};
