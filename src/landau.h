#pragma once
#include <cmath>

// ---------------------------------------------------------------------------
//  Landau-frame energy density of a sum of three ideal fluids.
//
//  MultiHydro used to build the mixed energy-momentum tensor
//
//      T^mu_nu = sum_I (e_I + p_I) u_I^mu u_{I,nu} - P delta^mu_nu ,
//                                                    P = sum_I p_I
//
//  as a ROOT TMatrixD and hand it to TMatrixDEigen -- a general, complex-
//  capable 4x4 eigensolver with heap allocation -- once per cell per timestep.
//  On a 121x121x161 grid that is 2.4 million eigendecompositions per call to
//  findFreezeout().
//
//  The decomposition does not need a general solver.  The first term is a sum
//  of only three rank-1 pieces, so write
//
//      A^mu_nu = sum_I w_I u_I^mu u_{I,nu} ,   w_I = e_I + p_I ,
//      T^mu_nu = A^mu_nu - P delta^mu_nu .
//
//  A and T share eigenvectors and their eigenvalues differ by the constant P.
//  Acting with A on a vector of the form v = sum_J c_J u_J gives
//
//      A v = sum_I w_I u_I (u_I . v) = sum_I u_I sum_J [ w_I (u_I . u_J) ] c_J ,
//
//  so every eigenvalue of A outside its 4-3 dimensional kernel is an
//  eigenvalue of the 3x3 matrix
//
//      B_IJ = w_I G_IJ ,   G_IJ = u_I . u_J   (Minkowski, +---, G_II = 1).
//
//  B = W G with W = diag(w_I) is not symmetric, but it is similar to
//  S = W^1/2 G W^1/2, which is: S = W^-1/2 B W^1/2.  So the eigenvalues are
//  real, a symmetric 3x3 solver suffices, and the eigenvectors map back as
//  c = W^1/2 y.  This is exact, not an approximation of the 4x4 problem.
//
//  Fluids with w_I = 0 need no special casing: row/column I of S vanishes, its
//  eigenvector is the unit vector e_I, and c = W^1/2 e_I = 0 gives v = 0, which
//  fails the time-likeness test below exactly as it should.
//
//  The kernel of A (eigenvalue 0 of A, i.e. -P for T) is the Minkowski
//  orthogonal complement of span{u_I}.  That span contains a time-like vector,
//  so its complement is space-like and is rejected by the same test -- which is
//  what makes "first eigenvector that is time-like" a well-posed selection rule.
// ---------------------------------------------------------------------------

namespace landau {

// Jacobi eigenvalue iteration for a real symmetric 3x3 matrix.
// On return d[] holds the eigenvalues and column k of V[][] the corresponding
// eigenvector.  Chosen over a closed-form cubic because it stays accurate for
// nearly degenerate spectra, which is precisely the regime of this problem:
// the three fluid velocities become equal wherever unification kicks in.
inline void jacobi3(double a[3][3], double d[3], double V[3][3]) {
 double A[3][3];
 for (int i = 0; i < 3; i++)
  for (int j = 0; j < 3; j++) {
   A[i][j] = a[i][j];
   V[i][j] = (i == j) ? 1.0 : 0.0;
  }
 for (int sweep = 0; sweep < 20; sweep++) {
  double off = fabs(A[0][1]) + fabs(A[0][2]) + fabs(A[1][2]);
  if (off < 1e-300) break;
  double scale = fabs(A[0][0]) + fabs(A[1][1]) + fabs(A[2][2]);
  if (off <= 1e-16 * scale) break;
  for (int p = 0; p < 2; p++)
   for (int q = p + 1; q < 3; q++) {
    const double apq = A[p][q];
    if (fabs(apq) < 1e-300) continue;
    const double theta = 0.5 * (A[q][q] - A[p][p]) / apq;
    double t = (theta >= 0. ? 1.0 : -1.0) /
               (fabs(theta) + sqrt(theta * theta + 1.0));
    const double c = 1.0 / sqrt(t * t + 1.0);
    const double s = t * c;
    const double tauR = s / (1.0 + c);
    const double h = t * apq;
    A[p][p] -= h;
    A[q][q] += h;
    A[p][q] = A[q][p] = 0.0;
    for (int r = 0; r < 3; r++) {
     if (r != p && r != q) {
      const double arp = A[r][p], arq = A[r][q];
      A[r][p] = A[p][r] = arp - s * (arq + arp * tauR);
      A[r][q] = A[q][r] = arq + s * (arp - arq * tauR);
     }
     const double vrp = V[r][p], vrq = V[r][q];
     V[r][p] = vrp - s * (vrq + vrp * tauR);
     V[r][q] = vrq + s * (vrp - vrq * tauR);
    }
   }
 }
 for (int i = 0; i < 3; i++) d[i] = A[i][i];
}

// Landau energy density of the three-fluid mixture.
//   w[I]  = e_I + p_I  (enthalpy density of fluid I, clamped at 0)
//   P     = p_p + p_t + p_f
//   u[I]  = contravariant 4-velocity of fluid I, normalised u.u = 1
// notTimelike is set when no eigenvector came out time-like, so the caller can
// reproduce the diagnostic the ROOT version printed.
inline double energyDensity(const double w[3], double P,
                            const double u[3][4], bool &notTimelike) {
 notTimelike = false;
 const double gmunu[4] = {1., -1., -1., -1.};

 // How many fluids actually carry enthalpy here?  Over most of the grid the
 // answer is 0, 1 or 2, and those cases have closed-form answers -- worth
 // separating out, because the iterative branch below is ~10x their cost and
 // the three-fluid overlap region is a small part of the volume.
 int act[3], nAct = 0;
 for (int I = 0; I < 3; I++)
  if (w[I] > 0.) act[nAct++] = I;

 if (nAct == 0) return -P;  // no matter: T^mu_nu = -P delta^mu_nu

 if (nAct == 1) {
  // A is rank 1 with eigenvalue w_I along u_I, which is time-like by
  // construction.  w_J = e_J + p_J = 0 with both non-negative forces p_J = 0,
  // so P is that fluid's own pressure and the result is just its e.
  const double eDens = w[act[0]] - P;
  if (eDens >= 0.) return eDens;
 }

 // G_IJ = u_I . u_J, then S_IJ = sqrt(w_I w_J) G_IJ
 double sw[3];
 for (int I = 0; I < 3; I++) sw[I] = sqrt(w[I] > 0. ? w[I] : 0.);

 if (nAct == 2) {
  const int I = act[0], J = act[1];
  double g = 0.;
  for (int m = 0; m < 4; m++) g += gmunu[m] * u[I][m] * u[J][m];
  const double b = sw[I] * sw[J] * g;
  const double half = 0.5 * (w[I] + w[J]);
  const double dif = 0.5 * (w[I] - w[J]);
  const double disc = sqrt(dif * dif + b * b);
  // Two time-like unit vectors have u_I.u_J >= 1, so det S = w_I w_J (1-g^2)
  // is non-positive: exactly one eigenvalue is positive and it carries the
  // time-like eigenvector.  The test below is still applied, not assumed.
  const double lamPair[2] = {half + disc, half - disc};
  for (int k = 0; k < 2; k++) {
   const double lam = lamPair[k];
   // eigenvector of the 2x2 block, taken in whichever form is better scaled
   double yI, yJ;
   if (fabs(lam - w[I]) > fabs(lam - w[J])) {
    yI = lam - w[J];
    yJ = b;
   } else {
    yI = b;
    yJ = lam - w[I];
   }
   double v[4];
   const double cI = sw[I] * yI, cJ = sw[J] * yJ;
   for (int m = 0; m < 4; m++) v[m] = cI * u[I][m] + cJ * u[J][m];
   double vv = 0.;
   for (int m = 0; m < 4; m++) vv += gmunu[m] * v[m] * v[m];
   const double eDens = lam - P;
   if (vv > 0. && eDens >= 0.) return eDens;
  }
  // the kernel directions give -P; fall through to the general path only if
  // neither in-plane eigenvector qualified
  notTimelike = true;
  return (fabs(lamPair[0]) > fabs(P) ? lamPair[0] : 0.) - P;
 }

 double S[3][3];
 for (int I = 0; I < 3; I++)
  for (int J = I; J < 3; J++) {
   double g = 0.;
   for (int m = 0; m < 4; m++) g += gmunu[m] * u[I][m] * u[J][m];
   S[I][J] = S[J][I] = sw[I] * sw[J] * g;
  }

 double lam[3], Y[3][3];
 jacobi3(S, lam, Y);

 // Walk the eigenvalues in order of decreasing magnitude, exactly as the ROOT
 // version did (TMatrixDEigen sorts that way), and take the first one whose
 // eigenvector is time-like and whose energy density is non-negative.
 int order[3] = {0, 1, 2};
 for (int i = 0; i < 2; i++)
  for (int j = i + 1; j < 3; j++)
   if (fabs(lam[order[j]]) > fabs(lam[order[i]])) {
    const int tmp = order[i];
    order[i] = order[j];
    order[j] = tmp;
   }

 double fallback = lam[order[0]] - P;
 for (int k = 0; k < 3; k++) {
  const int c = order[k];
  const double eDens = lam[c] - P;
  // v^mu = sum_J c_J u_J^mu with c = W^1/2 y
  double v[4] = {0., 0., 0., 0.};
  for (int J = 0; J < 3; J++) {
   const double cJ = sw[J] * Y[J][c];
   for (int m = 0; m < 4; m++) v[m] += cJ * u[J][m];
  }
  double vv = 0.;
  for (int m = 0; m < 4; m++) vv += gmunu[m] * v[m] * v[m];
  if (vv > 0. && eDens >= 0.) return eDens;
 }
 notTimelike = true;
 return fallback;
}

}  // namespace landau
