class EoS;

// One table cell: the five interpolated quantities are stored interleaved so
// that a cell occupies 40 contiguous bytes and the (inq, inq+1) pair needed by
// the trilinear interpolation sits together.  The original layout kept five
// separate double arrays, so a single call touched up to 40 scattered cache
// lines instead of ~8.
struct EoSHadronNode {
 double p, T, mub, muq, mus;
};

class EoSHadron : public EoS {
 double e0, n0, logemin, logemax, lognmax;
 int ne, nnb, nnq;
 double nb_abs_min, nq_abs_min, e_min;
 // grid spacings and their reciprocals, precomputed once.  The original
 // recomputed dxe, dxnb, dxnq -- three divisions -- on every eos()/p() call,
 // then used six more divisions to form the interpolation weights.
 double dxe, dxnb, dxnq, inv_dxe, inv_dxnb, inv_dxnq;
 double inv_e0, inv_n0;
 EoSHadronNode *tab;
 int *statustab;
 double p_elow, T_elow;   // ptab/Ttab at the (0, nnb/2, nnq/2) reference cell
 inline int index3(int ie, int inb, int inq) const {
  return inq + nnq * inb + nnq * nnb * ie;
 }
 // Locate the cell and form the eight trilinear weights.  Shared by eos() and
 // p(), which previously held two verbatim copies of this code.
 // Returns false if the cell is flagged unusable in statustab.
 bool locate(double e, double nb, double nq, const EoSHadronNode *&base,
             double w[8]) const;

public:
 EoSHadron(char *filename);
 ~EoSHadron(void);

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);
};
