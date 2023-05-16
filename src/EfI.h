#include <cmath>
class CrossSections;

class EfIntegrand{
private:
    double Tf, valphatilde, gammaalphatilde;
    CrossSections *xsect;
public:
    EfIntegrand(CrossSections *_xsect, double _Tf, double _valphatilde);
    double Eval(double s);
};
