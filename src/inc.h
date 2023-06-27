#pragma once

// Q[] component indexes
#define T_ 0
#define X_ 1
#define Y_ 2
#define Z_ 3
#define NB_ 4
#define NQ_ 5
#define NS_ 6

#define PREDICT 0
#define CORRECT 1

#define C_1D 0
#define C_2D 1

const double C_PI = 3.14159265358979312;

#define BAG_ASYMPT  // necessary if you use HIRANO EoS !!!

const double gevtofm = 5.067728853; //funny way of writing (hbar*c)^(-1)
const double mN = 0.939; //nucleon mass [GeV]
const double mpi = 0.1396; //pion mass [GeV]


#ifndef _DEBUG
//#define UI
#endif
