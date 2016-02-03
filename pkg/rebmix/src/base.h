#ifndef BASE_H_INCLUDED
#define BASE_H_INCLUDED

#include <float.h>

#ifndef _MEMORY_LEAK_SWITCH
#define _MEMORY_LEAK_SWITCH 0 
#endif

#ifndef FLOAT
#define FLOAT double 
#endif

#ifndef FLOAT_MIN
#define FLOAT_MIN DBL_MIN 
#endif

#ifndef FLOAT_MAX
#define FLOAT_MAX DBL_MAX 
#endif

#ifndef FLOAT_EPSILON
#define FLOAT_EPSILON DBL_EPSILON  
#endif

#ifndef Sqrt2
#define Sqrt2 (FLOAT)1.4142135623730950488016887242097
#endif               

#ifndef Pi            
#define Pi (FLOAT)3.1415926535897932384626433832795
#endif

#ifndef SqrtPi            
#define SqrtPi (FLOAT)1.7724538509055160272981674833411
#endif

#ifndef Sqrt2Pi
#define Sqrt2Pi (FLOAT)2.506628274631000502415765284811
#endif

#ifndef LogPi
#define LogPi (FLOAT)1.1447298858494001741434273513531
#endif

#ifndef LogSqrt2Pi
#define LogSqrt2Pi (FLOAT)0.91893853320467274178032973640562
#endif 

#ifndef Euler       
#define Euler (FLOAT)0.5772156649015328606065120900824
#endif

#ifndef Phi       
#define Phi (FLOAT)1.6180339887498948482045868343656
#endif  

#ifndef Eps
#define Eps (FLOAT)0.000001
#endif

#ifndef ItMax
#define ItMax 1000
#endif

#ifndef BufInc
#define BufInc 1000
#endif

#ifndef IA
#define IA 16807
#endif

#ifndef IM
#define IM 2147483647
#endif

#ifndef AM
#define AM (FLOAT)(1.0 / IM)
#endif

#ifndef IQ
#define IQ 127773
#endif

#ifndef IR
#define IR 2836
#endif

#ifndef NTAB
#define NTAB 32
#endif

#ifndef NDIV
#define NDIV (1 + (IM - 1) / NTAB)
#endif

#ifndef RNMX
#define RNMX (FLOAT)(1.0 - 1.2E-7)
#endif

#define Min(x, y) ((x < y) ? x : y) 

#define Max(x, y) ((x > y) ? x : y) 

#define IsNan(x) ((x) != (x)) 

#define IsInf(x) (!IsNan(x) && IsNan((x) - (x)))

typedef enum {
    pfNormal,    // Normal distribution.
    pfLognormal, // Lognormal distribution.
    pfWeibull,   // Weibull distribution.
    pfGamma,     // Gamma distribution.
    pfBinomial,  // Binomial distribution.
    pfPoisson,   // Poisson distribution.
    pfDirac      // Dirac distribution.
} ParametricFamilyType_e;

typedef struct componetdistributiontype {
    ParametricFamilyType_e *pdf;    // Parametric family types.
    FLOAT                  *Theta1; // Component parameters.
    FLOAT                  *Theta2; // Component parameters.
} ComponentDistributionType;

// length_pdf_ = d_ for class Rebmix.
// length_Theta_ = 2 for class Rebmix.
// length_theta_[0] = d_ and length_theta_[1] = d_ for class Rebmix.

class Base {
public:
    // Members.
    int Trigger_;       // Trigger.
    int length_pdf_;    // Length of pdf_. 
    int length_Theta_;  // Length of Theta_.
    int *length_theta_; // Length of Theta_[i].
    // Constructor.
    Base();
    // Destructor.
    ~Base();
}; // Base

FLOAT Ran1(int *IDum);

// Inserts y into ascending list Y of length n. Set n = 0 initially.

void Insert(FLOAT y,   // Inserted value.
            int   *n,  // Length of Y.
            FLOAT *Y); // Pointer to Y = [y0,...,yn-1].

// Returns the value log(Gamma(y)) for y > 0. See http://www.nr.com/.

FLOAT Gammaln(FLOAT y);

// Returns the digamma for y > 0. See http://www.nr.com/.

int Digamma(FLOAT y, FLOAT *Psi);

// Returns the inverse of the binomial c.d.f. for the specified n and p.

FLOAT BinomialInv(FLOAT Fy, FLOAT n, FLOAT p);

// Returns the inverse of the Poisson c.d.f. for the specified Theta.

FLOAT PoissonInv(FLOAT Fy, FLOAT Theta);

// Returns the incomplete gamma function P(a, y) evaluated by its series
// representation as GamSer. Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

int GammaSer(FLOAT a,       // Constant a > 0.
             FLOAT y,       // Variable y > 0.
             FLOAT *GamSer, // Incomplete gamma function.
             FLOAT *Gamln); // Log(Gamma(a)).

// Returns the incomplete gamma function Q(a, y) evaluated by its continued
// fraction representation as GamCfg. Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

int GammaCfg(FLOAT a,       // Constant a > 0.
             FLOAT y,       // Variable y > 0.
             FLOAT *GamCfg, // Incomplete gamma function.
             FLOAT *Gamln); // Log(Gamma(a)).

// Returns the incomplete gamma function P(a, y). Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

int GammaP(FLOAT a,       // Constant a > 0.
           FLOAT y,       // Variable y > 0.
           FLOAT *GamP,   // Incomplete gamma function.
           FLOAT *Gamln); // Log(Gamma(a)).

// Returns the inverse of the gamma c.d.f. for the specified Theta and Beta. See http://www.nr.com/.

int GammaInv(FLOAT Fy, FLOAT Theta, FLOAT Beta, FLOAT *y);

// Returns the error function erf(y). See http://www.nr.com/.

int ErrorF(FLOAT y,     // Variable y.
           FLOAT *ErF); // Error function.

// Returns the determinant and the inverse matrix of A. See http://www.nr.com/

int LUinvdet(int   n,      // Size of square matrix.
             FLOAT *A,     // Pointer to the square matrix A.
             FLOAT *Ainv,  // Pointer to the inverse matrix of A.
             FLOAT *Adet); // Pointer to the determinant of A.

// Returns the Cholesky decomposition of matrix A. See http://www.nr.com/ 

int Choldc(int   n,   // Size of square matrix.
           FLOAT *A,  // Pointer to the square matrix A.
           FLOAT *L); // Lower triangular factors.

#endif