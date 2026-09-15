class EoS;

// One table cell: the five interpolated quantities interleaved into a single
// 40-byte record, so a cell and its (inq, inq+1) neighbour sit together.  The
// original kept five separate double arrays, so one call touched up to 40
// scattered cache lines instead of ~8.
struct EoSSmashNode {
 double p, T, mub, muq, mus;
};

class EoSSmash : public EoS {
  // Bounds (upper and lower) of energy density, baryon density and charge density
  double emax, emin, nbmax, nbmin, nqmax, nqmin;
  // Number of points in energy density, baryon density, charge density
  int ne, nnb, nnq;
  // Grid spacings and their reciprocals, precomputed once.  The original
  // recomputed de, dnb, dnq -- three divisions -- on every eos()/p() call and
  // then used six more divisions to build the interpolation weights.
  double de, dnb, dnq, inv_de, inv_dnb, inv_dnq;
  // Tabularised values for pressure, Temperature, muB, muQ, muS
  // Dimension: ne * nnb * nnq
  EoSSmashNode *tab;
  double p_elow, T_elow;  // p and T at the (0, nnb/2, nnq/2) reference cell

  // Get absolute index from given indices in energy density, baryon density
  // and charge density
  inline int index3(int ie, int inb, int inq) const {
   return inq + nnq * inb + nnq * nnb * ie;
  }
  // Locate the cell and form the eight trilinear weights.  Shared by eos()
  // and p(), which previously held two verbatim copies of this code.
  void locate(double e, double nb, double nq, const EoSSmashNode *&base,
              double w[8]) const;

public:
 EoSSmash(char* filename, int Ne, int Nnb, int Nq);
 ~EoSSmash(void);

 // Set (T, muB, muQ, muS, p) from (e, nb, nq, ns) according to EoS
 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 // Calculate pressure from (e, nb, nq, ns)
 virtual double p(double e, double nb, double nq, double ns);
};
