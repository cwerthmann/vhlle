#include <cmath>
class CrossSections;

class EfIntegrand{
private:
    double Tf, mubf, valphatilde, gammaalphatilde;
    CrossSections *xsect;
public:
    EfIntegrand(CrossSections *_xsect, double _Tf, double _mubf, double _valphatilde);
    double EvalNpi(double s);
    double EvalNN(double s);
};
