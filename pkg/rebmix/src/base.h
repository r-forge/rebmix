#ifndef BASE_H_INCLUDED
#define BASE_H_INCLUDED

#ifdef _MSC_VER
#pragma warning(disable: 4514)
#pragma warning(disable: 4710)
#pragma warning(disable: 4820)
#pragma warning(disable: 5045)
#endif

#include <float.h>

#ifndef _MEMORY_LEAK_SWITCH
#define _MEMORY_LEAK_SWITCH 0
#endif

#ifndef _MAINTAIN_SWITCH
#define _MAINTAIN_SWITCH 0
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

#ifndef Pi2
#define Pi2 (FLOAT)6.2831853071795864769252867665590
#endif

#ifndef SqrtPi
#define SqrtPi (FLOAT)1.7724538509055160272981674833411
#endif

#ifndef SqrtPi2
#define SqrtPi2 (FLOAT)2.506628274631000502415765284811
#endif

#ifndef LogPi
#define LogPi (FLOAT)1.1447298858494001741434273513531
#endif

#ifndef LogPi2
#define LogPi2 (FLOAT)1.8378770664093454835606594728112
#endif

#ifndef LogSqrtPi2
#define LogSqrtPi2 (FLOAT)0.91893853320467274178032973640562
#endif

#ifndef Euler
#define Euler (FLOAT)0.5772156649015328606065120900824
#endif

#ifndef SqrPi6
#define SqrPi6 (FLOAT)1.644934066848226436472415166646
#endif

#ifndef Phi
#define Phi (FLOAT)1.6180339887498948482045868343656
#endif

#ifndef Eps
#define Eps (FLOAT)0.00001
#endif

#ifndef TLA_MAX_LN
#define TLA_MAX_LN 46.0
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
    pfGumbel,    // Gumbel distribution.
    pfvonMises,  // Von Mises distribution.
    pfBinomial,  // Binomial distribution.
    pfPoisson,   // Poisson distribution.
    pfDirac,     // Dirac distribution.
    pfUniform    // Uniform distribution.
} ParametricFamilyType_e;

typedef enum {
    varEM,     // Classic Expectation-Maximization algorithm (Soft (continious posterior)).
    varECM,    // Expectation - Conditional - Maximization algorithm (Hard (posterior 0 or 1) EM).
} EmVariantType_e;

typedef enum {
    acc_fixed,  // Fixed constant multiplier acceleration of the EM algorithm.
    acc_line,   // Line search for multiplier acceleration of the EM algorithm.
    acc_golden  // Golden search for multiplier accleration of the EM algorithm.
} EmAccelerationType_e;

typedef enum {
    strategy_none,       // EM algorithm is not employed for estimation of mixture model parameters.
    strategy_single,     // Single REBMIX + EM strategy.
    strategy_best,       // Best REBMIX + EM strategy.
    strategy_exhaustive  // Exhaustive REBMIX + EM strategy.
} EmStrategyType_e;

typedef struct summaryparametertype {
    int   c;     // Optimal number of components.
    int   k;     // Optimal v or optimal k.
    FLOAT *y0;   // Optimal origins of length d.
    FLOAT *ymin; // Minimum observations.
    FLOAT *ymax; // Maximum observations.
    FLOAT *h;    // Optimal class widths of length d.
    FLOAT IC;    // Optimal information criterion.
    FLOAT logL;  // Log-likelihood.
    int   M;     // Degrees of freedom.
} SummaryParameterType;

typedef struct additinalparametertype {
    int Bracket; // 1 for bracketing and 0 for golden section.
    int a;       // Golden section constant.
    int b;       // Golden section constant.
    int c;       // Golden section constant.
    int d;       // Golden section constant.
} AdditionalParameterType;

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

class CompnentDistribution : public Base {
public:
    // Members.
    Base                   *owner_;  // Owner object.
    ParametricFamilyType_e *pdf_;    // Parametric family types.
    FLOAT                  **Theta_; // Component parameters.
    // Constructor.
    CompnentDistribution(Base *owner);
    // Destructor.
    ~CompnentDistribution();
    // Methods.
    int Realloc(int length_pdf, int length_Theta, int *length_theta);
    int Memmove(CompnentDistribution *CmpTheta);
}; // CompnentDistribution

typedef struct mixtureparametertype {
    FLOAT                *W;          // Pointer to weight.
    CompnentDistribution **MixTheta;  // Pointer to mixture parameters.
    int                  c;           // Number of components in mixture.
    FLOAT                logL;        // Estimated value of log likelihood.
    FLOAT                logV;        // Logaritmic value of V.
    int                  k;           // Bin number for histogram preprocessing or smothing parameter for kernel density estimation and k-nearest neighbour preprocessing.
    FLOAT                *h;          // Bin widths for histogram preprocessing or smothing parameter for kernel density estimation and k-nearest neighbour preprocessing.
    FLOAT                *y0;         // Bin origins for histogram preprocessing or NULL for kernel density estimation and k-nearest neighbour preprocessing.
    FLOAT                *ymin;       // Minimum values from data.
    FLOAT                *ymax;       // Maximum values from data.
    int                  n_iter_em;   // Number of performed iterations of EM algorithm.
    int                  initialized; // Boolean indicator if struct contains mixture model parameters. 
} MixtureParameterType;

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

// Returns the inverse of the Weibull c.d.f. for the specified Theta and Beta.

FLOAT WeibullInv(FLOAT Fy, FLOAT Theta, FLOAT Beta);

// Returns the inverse of the Gumbel c.d.f. for the specified Mean and Beta.

FLOAT GumbelInv(FLOAT Fy, FLOAT Mean, FLOAT Beta);

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

// Returns the determinant and the inverse matrix of A. See http ://www.nr.com/

int Cholinvdet(int   n,         // Size of square matrix.
               FLOAT *A,        // Pointer to the symmetric square matrix A.
               FLOAT *Ainv,     // Pointer to the inverse matrix of A.
               FLOAT *logAdet); // Pointer to the logarithm of determinant of A.

// Returns modified Bessel function of order 0. See http://people.math.sfu.ca/~cbm/aands/page_378.htm

FLOAT BesselI0(FLOAT y);

// Returns modified Bessel function of order 1. See http://people.math.sfu.ca/~cbm/aands/page_378.htm

FLOAT BesselI1(FLOAT y);

// Returns the inverse of the von Mises c.d.f. for the specified Mean and Kappa.

FLOAT vonMisesInv(FLOAT Fy, FLOAT Mean, FLOAT Kappa);

// Returns x * log(x).

FLOAT xlogx(FLOAT x);

#endif
