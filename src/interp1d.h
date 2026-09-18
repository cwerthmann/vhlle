#pragma once
#include <cstddef>
#include <vector>

// ---------------------------------------------------------------------------
//  Drop-in replacement for the ROOT TGraph objects that eos.cpp and xsect.cpp
//  used purely as one-dimensional lookup tables.
//
//  Two reasons for it, in order of importance:
//
//  1. Thread safety.  TGraph::Eval() is a non-const virtual member of a ROOT
//     object; calling it concurrently from an OpenMP region is not supported.
//     EoSs::p() sits at the bottom of transformPV()'s Newton iteration and is
//     therefore called from every parallel loop in the code, and
//     CrossSections::Ivanov()/piN() are called from the parallel friction
//     substep.  This class holds no mutable state at all, so any number of
//     threads may evaluate it at once.
//
//  2. Speed.  Eval() costs a virtual dispatch plus a binary search even when
//     the abscissae are equally spaced, which they are for every table shipped
//     with the code.  Uniform spacing is detected once at construction and the
//     lookup then becomes a multiply and a truncation.
//
//  The interpolation reproduces TGraph::Eval(x, 0, "") exactly: locate the
//  largest i with x[i] <= x, clamp the bracket to the ends of the table (so
//  out-of-range arguments are linearly extrapolated, as TGraph does), and
//  evaluate with the identical expression so the result is bit-for-bit the
//  same as before.  Abscissae must be sorted ascending, which was already a
//  precondition of TGraph::Eval.
// ---------------------------------------------------------------------------
class Interp1D {
 std::vector<double> xs, ys;
 double x0 = 0., invdx = 0.;
 bool uniform = false;

public:
 Interp1D() {}

 template <typename TX, typename TY>
 Interp1D(int n, const TX *x, const TY *y) {
  init(n, x, y);
 }

 template <typename TX, typename TY>
 void init(int n, const TX *x, const TY *y) {
  xs.resize(n);
  ys.resize(n);
  for (int i = 0; i < n; i++) {
   xs[i] = (double)x[i];
   ys[i] = (double)y[i];
  }
  uniform = false;
  if (n > 2) {
   const double dx = (xs[n - 1] - xs[0]) / (n - 1);
   if (dx > 0.) {
    uniform = true;
    // tolerance is relative to the span, so a table written with a fixed
    // number of decimals still counts as uniform
    const double tol = 1e-9 * (xs[n - 1] - xs[0]);
    for (int i = 1; i < n - 1; i++) {
     const double d = xs[i] - (xs[0] + i * dx);
     if (d > tol || d < -tol) {
      uniform = false;
      break;
     }
    }
    if (uniform) {
     x0 = xs[0];
     invdx = 1.0 / dx;
    }
   }
  }
 }

 int size() const { return (int)xs.size(); }

 double eval(double x) const {
  const int n = (int)xs.size();
  if (n == 0) return 0.;
  if (n == 1) return ys[0];
  int low;
  if (uniform) {
   low = (int)((x - x0) * invdx);
  } else {
   // largest index with xs[low] <= x, or -1; same contract as
   // TMath::BinarySearch, which is what TGraph::Eval used
   int lo = 0, hi = n - 1;
   if (x < xs[0])
    low = -1;
   else {
    while (hi - lo > 1) {
     const int mid = (lo + hi) / 2;
     if (xs[mid] <= x)
      lo = mid;
     else
      hi = mid;
    }
    low = lo;
   }
  }
  if (low < 0) low = 0;
  if (low > n - 2) low = n - 2;
  const int up = low + 1;
  if (xs[low] == xs[up]) return ys[low];
  return ys[up] + (x - xs[up]) * (ys[low] - ys[up]) / (xs[low] - xs[up]);
 }
};
