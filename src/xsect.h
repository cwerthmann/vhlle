#include <cmath>

class TGraph;

class CrossSections {
 TGraph *gSigmaT, *gSigmaE, *gSigmaPimp, *gSigmaPipp;
public:
 CrossSections(void);
 ~CrossSections(void);
 void Ivanov(double Ekin, double& sigmaT, double& sigmaE, double& sigmaP);
 void Ivanovbar(double sqrts, double& sigmaEbar, double& sigmaPbar, double& sigmaRbar);
 void NN(double Ekin, double& sigmaNN);
 double piN(double sqrts);
};
