#ifdef __cplusplus /* If this is a C++ compiler, use C linkage. */
extern "C" {
#endif

#ifndef REBMIXF_H_INCLUDED
#define REBMIXF_H_INCLUDED

#include <stdlib.h>
#include <string.h>

#ifndef _MEMORY_LEAK_SWITCH
#define _MEMORY_LEAK_SWITCH 0 
#endif

#ifndef _TIME_LEFT_SWITCH
#define _TIME_LEFT_SWITCH 0 
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

#define Min(x, y) ((x < y) ? x : y) 
#define Max(x, y) ((x > y) ? x : y) 
#define IsNan(x) ((x) != (x)) 
#define IsInf(x) (!IsNan(x) && IsNan((x) - (x))) 

typedef enum {
    poHistogram,         /* Histogram approach. */
    poParzenWindow,      /* Parzen window. */
    poKNearestNeighbour  /* K-nearest neighbour. */
} PreprocessingType_e;

typedef enum {
    vtContinuous,  /* Continuous variable. */
    vtDiscrete     /* Ordered or non-ordered binary or discrete variable. */
} VariablesType_e;

typedef enum {
    pfNormal,      /* Normal distribution. */
    pfLognormal,   /* Lognormal distribution. */
    pfWeibull,     /* Weibull distribution. */
    pfGamma,       /* Gamma distribution. */
    pfBinomial,    /* Binomial distribution. */
    pfPoisson,     /* Poisson distribution. */
    pfDirac        /* Dirac distribution. */
} ParametricFamilyType_e;

typedef enum {
    rtRigid,  /* Rigid restraints. */
    rtLoose   /* Loose restraints. */
} PestraintsType_e;

typedef struct marginaldistributiontype {
    ParametricFamilyType_e ParFamType; /* Parametric family. */
    FLOAT                  Par0;       /* Parameter 0. */ 
    FLOAT                  Par1;       /* Parameter 0. */ 
} MarginalDistributionType;

typedef enum {
    icAIC,    /* AIC - Akaike information criterion Akaike (1973). */ 
    icAIC3,   /* AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
    icAIC4,   /* AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
    icAICc,   /* AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989). */
    icBIC,    /* BIC - Bayesian information criterion Schwarz (1978). */
    icCAIC,   /* CAIC - Consistent Akaike information criterion Bozdogan (1987). */
    icHQC,    /* HQC - Hannan-Quinn information criterion Hannan & Quinn (1979). */
    icMDL2,   /* MDL2 - Minimum description length Liang et al. (1992). */
    icMDL5,   /* MDL5 - Minimum description length Liang et al. (1992). */
    icAWE,    /* AWE - Approximate weight of evidence criterion Banfield & Raftery (1993). */ 
    icCLC,    /* CLC - Classification likelihood criterion Biernacki & Govaert (1997). */
    icICL,    /* ICL - Integrated classification likelihood Biernacki et al. (1998). */
    icPC,     /* PC - Partition coefficient Bezdek (1981). */
    icICLBIC, /* ICL-BIC - Integrated classification likelihood criterion Biernacki et al. (1998). */ 
    icD,      /* D - Total of positive relative deviations Nagode & Fajdiga (2011). */
    icSSE,    /* SSE - Sum of squares error Bishop (1998). */
} InformationCriterionType_e;

typedef struct roughparametertype {
    FLOAT ym;    /* Mode position. */
    FLOAT flm;   /* Component conditional empirical density. */
    FLOAT klm;   /* Component conditional total number of observations. */
} RoughParameterType;

typedef struct inputrebmixparametertype {
    char                       *curr;       /* Path to the currently open data file. */
    int                        o;           /* Number of paths. */ 
    char                       **open;      /* Paths to open data files. */
    PreprocessingType_e        PreType;     /* Preprocessing of observations. */
    int                        cmax;        /* Maximum number of components. */  
    InformationCriterionType_e ICType;      /* Information criterion types. */
    int                        d;           /* Number of independent random variables. */ 
    VariablesType_e            *VarType;    /* Types of variables. */
    ParametricFamilyType_e     *ParFamType; /* Parametric family types. */
    FLOAT                      *Par0;       /* Initial component parameters. */
    FLOAT                      *Par1;       /* Initial component parameters. */
    int                        kmax;        /* Number of classes or k-nearest neighbours to be processed. */
    int                        *K;          /* Numbers of classes or k-nearest neighbours. */
    FLOAT                      *y0;         /* Origins. */
    FLOAT                      *ymin;       /* Minimum observations. */
    FLOAT                      *ymax;       /* Maximum observations. */
    FLOAT                      ar;          /* Acceleration rate. */
    PestraintsType_e           ResType;     /* Restraints type. */
    char                       *save;       /* Path to the save data file. */
} InputREBMIXParameterType;

typedef struct outputrebmixparametertype {
    int                      k;         /* Optimal number of classes or k-nearest neighbours. */
    FLOAT                    *h;        /* Optimal class widths. */ 
    FLOAT                    *y0;       /* Optimal origins. */
    FLOAT                    IC;        /* Optimal information criterion. */
    FLOAT                    logL;      /* Log-likelihood. */
    int                      M;         /* Degrees of freedom. */
    int                      c;         /* Optimal number of components. */ 
    FLOAT                    *W;        /* Optimal component weights. */
    MarginalDistributionType **Theta;   /* Optimal parameters. */
    int                      n;         /* Total number of independent observations. */
    FLOAT                    **X;       /* Pointer to the input observations [x0,...,xd-1]. */
} OutputREBMIXParameterType;

typedef struct optrebmixparametertype {
    int   Imax;   /* Number of iterations for optimal number of classes or k-nearest neighbours. */
    int   *c;     /* Numbers of components for optimal number of classes or k-nearest neighbours. */ 
    FLOAT *IC;    /* Information criteria for optimal number of classes or k-nearest neighbours. */  
    FLOAT *logL;  /* Log-likelihoods for optimal number of classes or k-nearest neighbours. */
    FLOAT *D;     /* Totals of positive relative deviations for optimal number of classes or k-nearest neighbours. */ 
} OptREBMIXParameterType;

typedef struct allrebmixparametertype {
    int   Bracket; /* 1 for bracketing and 0 for golden section.*/
    int   a;       /* Golden section constant. */
    int   b;       /* Golden section constant. */
    int   c;       /* Golden section constant. */
    int   d;       /* Golden section constant. */
    int   kmax;    /* Number of classes or k-nearest neighbours. */
    int   *K;      /* Numbers of classes or k-nearest neighbours. */
    FLOAT *IC;     /* Information criteria for numbers of classes or k-nearest neighbours. */  
} AllREBMIXParameterType;

/* Returns the value log(Gamma(y)) for y > 0. See http://www.nrbook.com/a/bookcpdf/c6-1.pdf */

FLOAT Gammaln(FLOAT y);

/* Returns the digamma for y > 0. */

int Digamma(FLOAT y, FLOAT *Psi);

/* Returns the inverse of the gamma c.d.f. for the specified Theta and Beta. */

int GammaInv(FLOAT Fy, FLOAT Theta, FLOAT Beta, FLOAT *y);

/* Returns component p.d.f or c.d.f. */ 

int ComponentDist(int                      d,            /* Number of independent random variables. */
                  FLOAT                    *Y,           /* Pointer to the input point [y0,...,yd-1]. */
                  MarginalDistributionType *MrgDistType, /* Marginal distribution type. */
                  FLOAT                    *CmpDist,     /* Component distribution. */
                  int                      Cumulative);  /* Set 1 if c.d.f. or 0 if p.d.f. */

/* Returns mixture p.d.f or c.d.f. */ 

int MixtureDist(int                      d,             /* Number of independent random variables. */
                FLOAT                    *Y,            /* Pointer to the input point [y0,...,yd-1]. */
                int                      c,             /* Number of components. */ 
                FLOAT                    *W,            /* Component weights. */
                MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                FLOAT                    *MixDist,      /* Mixture distribution. */
                int                      Cumulative);   /* Set 1 if c.d.f. or 0 if p.d.f. */

/* REBMIX algorithm. */

int REBMIX(InputREBMIXParameterType   *InpParType,  /* Input parameters. */ 
           OutputREBMIXParameterType  *OutParType,  /* Output parameters. */
           OptREBMIXParameterType     *OptParType,  /* Optimal parameters. */ 
           AllREBMIXParameterType     *AllParType); /* All parameters. */ 

/* Reads input data from the file stream. */

int ReadREBMIXDataFile(InputREBMIXParameterType  *InpParType,  /* Input parameters. */ 
                       OutputREBMIXParameterType *OutParType); /* Output parameters. */

/* Writes input and output parameters into the file stream. */

int WriteREBMIXParameterFile(InputREBMIXParameterType  *InpParType,   /* Input parameters. */ 
                             OutputREBMIXParameterType *OutParType);  /* Output parameters. */

/* Runs REBMIX template file stream. */

int RunREBMIXTemplateFile(char *file); /* File stream. */

/* Preprocessing of observations for k-nearest neighbour. */

int PreprocessingKNN(int   k,    /* k-nearest neighbours. */
                     FLOAT *h,   /* Normalizing vector. */
                     int   n,    /* Total number of independent observations. */
                     int   d,    /* Number of independent random variables. */ 
                     FLOAT **Y); /* Pointer to the input array [y0,...,yd-1,kl,V,R]. */

/* Preprocessing of observations for Parzen window. */

int PreprocessingPW(FLOAT *h,   /* Sides of the hypersquare. */
                    int   n,    /* Total number of independent observations. */
                    int   d,    /* Number of independent random variables. */ 
                    FLOAT **Y); /* Pointer to the input array [y0,...,yd-1,kl,k]. */

/* Preprocessing of observations for histogram. */

int PreprocessingH(FLOAT                  *h,          /* Sides of the hypersquare. */
                   FLOAT                  *y0,         /* Origins. */
                   ParametricFamilyType_e *ParFamType, /* Parametric family types. */
                   int                    *k,          /* Total number of bins. */
                   int                    n,           /* Total number of independent observations. */
                   int                    d,           /* Number of independent random variables. */ 
                   FLOAT                  **X,         /* Pointer to the input points [x0,...,xd-1]. */
                   FLOAT                  **Y);        /* Pointer to the input array [y0,...,yd-1,kl]. */

/* Returns information criterion for k-nearest neighbour. */ 

int InformationCriterionKNN(InformationCriterionType_e ICType,        /* Information criterion type. */
                            int                        k,             /* k-nearest neighbours. */
                            int                        n,             /* Total number of independent observations. */
                            int                        d,             /* Number of independent random variables. */ 
                            FLOAT                      **Y,           /* Pointer to the input points [y0,...,yd-1,kl,V,R]. */
                            int                        c,             /* Number of components. */ 
                            FLOAT                      *W,            /* Component weights. */
                            MarginalDistributionType   **MrgDistType, /* Marginal distribution type. */
                            FLOAT                      *IC,           /* Information criterion. */
                            FLOAT                      *logL,         /* log-likelihood. */
                            int                        *M,            /* Degrees of freedom. */
                            FLOAT                      *D);           /* Total of positive relative deviations. */

/* Returns information criterion for Parzen window. */ 

int InformationCriterionPW(InformationCriterionType_e ICType,        /* Information criterion type. */
                           FLOAT                      V,             /* Volume of the hypersquare. */
                           int                        n,             /* Total number of independent observations. */
                           int                        d,             /* Number of independent random variables. */ 
                           FLOAT                      **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                           int                        c,             /* Number of components. */ 
                           FLOAT                      *W,            /* Component weights. */
                           MarginalDistributionType   **MrgDistType, /* Marginal distribution type. */
                           FLOAT                      *IC,           /* Information criterion. */
                           FLOAT                      *logL,         /* log-likelihood. */
                           int                        *M,            /* Degrees of freedom. */
                           FLOAT                      *D);           /* Total of positive relative deviations. */

/* Returns information criterion for histogram. */ 

int InformationCriterionH(InformationCriterionType_e ICType,        /* Information criterion type. */
                          FLOAT                      V,             /* Volume of the hypersquare. */
                          int                        k,             /* Total number of bins. */
                          int                        n,             /* Total number of independent observations. */
                          int                        d,             /* Number of independent random variables. */ 
                          FLOAT                      **Y,           /* Pointer to the input points [y0,...,yd-1,kl]. */
                          int                        c,             /* Number of components. */ 
                          FLOAT                      *W,            /* Component weights. */
                          MarginalDistributionType   **MrgDistType, /* Marginal distribution type. */
                          FLOAT                      *IC,           /* Information criterion. */
                          FLOAT                      *logL,         /* log-likelihood. */
                          int                        *M,            /* Degrees of freedom. */
                          FLOAT                      *D);           /* Total of positive relative deviations. */

#endif

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage. */
}
#endif
