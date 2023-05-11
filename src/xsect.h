class TGraph;

class CrossSections {
 TGraph *gSigmaPimp, *gSigmaPipp;
public:
 CrossSections(void);
 ~CrossSections(void);
 void NN(double Ekin, double& sigmaNN);
 double piN(double sqrts);
};
