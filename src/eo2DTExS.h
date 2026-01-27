// EoS2DTExS.h
#pragma once

#include "eos.h"

// Forward declaration of the auxiliary tilde-grid reader/interpolator
class EoS2DTExS_aux;

// Wrapper class implementing EoS using a 2D T-μB table
class EoS2DTExS : public EoS {
public:
    // Nt = number of Ttilde points, Nmb = number of muBtilde points
    EoS2DTExS(const char* filename = "eos/EoS2DTExS.dat", int Nt = 500, int Nmb = 500);
    virtual ~EoS2DTExS();

    // EoS interface
    virtual void eos(double e, double nb, double nq, double ns,
                     double& _T, double& _mub, double& _muq, double& _mus,
                     double& _p) override;

    virtual double p(double e, double nb, double nq, double ns) override;

    // Pressure helper without extra charges
    double p(double e, double nb);

private:
    EoS2DTExS_aux* eos_TExS; // pointer to the auxiliary interpolator
};

// Auxiliary class that stores the tilde-grid table and performs interpolation
class EoS2DTExS_aux {
private:
    // Physical bounds
    double emin, emax;
    double nmin, nmax;

    // Tilde-grid bounds
    double tildeTMin, tildeTMax;
    double tildeMuBMin, tildeMuBMax;

    // Grid sizes
    int NtGrid, NmbGrid;

    // Tables on the tilde grid
    double **ptab;    // pressure
    double **Ttab;    // temperature
    double **mubtab;  // mu_B
    double **mustab;  // mu_S
    double **stab;    // entropy density

public:
    // Constructor: load table from file
    EoS2DTExS_aux(const char* filename, int Nt, int Nmb);

    // Destructor: free allocated memory
    ~EoS2DTExS_aux();

    // Interpolated access to thermodynamic quantities
    void get(double e, double nb, double& p, double& T, double& mub, double& mus);

    // Convenience wrapper: pressure only
    double p(double e, double nb);
};
