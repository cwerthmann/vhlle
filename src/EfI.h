#include <cmath>
class CrossSections;

class EfIntegrand{
private:
    double Tf, mubf, valphatilde, gammaalphatilde;
    CrossSections *xsect;
public:
    EfIntegrand(CrossSections *_xsect, double _Tf, double _mubf, double _valphatilde);
    double EvalNpi_to_f(double p);
    double EvalNpi_to_pt(double p);
};
