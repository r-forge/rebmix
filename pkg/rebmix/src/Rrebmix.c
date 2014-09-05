#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "rngmixf.h"
#include "rebmixf.h"

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

/* Runs RNGMIX in R. */

void RRNGMIX(int    *IDum,        /* Random seed. */
             int    *d,           /* Number of independent random variables. */
             int    *c,           /* Number of components. */
             int    *N,           /* Total number of observations in class l. */
             char   **ParFamType, /* Parametric family types. */
             double *Par0,        /* Component parameters. */
             double *Par1,        /* Component parameters. */
             int    *n,           /* Total number of independent observations. */
             double *X,           /* Pointer to the output array X. */
             int    *Error)       /* Error code. */
{
    InputRNGMIXParameterType  InpParType;
    OutputRNGMIXParameterType OutParType;
    int                       i, j, k;

    memset(&InpParType, 0, sizeof(InputRNGMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputRNGMIXParameterType));

    InpParType.IDum = *IDum;
    InpParType.d = *d;
    InpParType.c = *c;

    InpParType.N = (int*)malloc(InpParType.c * sizeof(int));

    *Error = NULL == InpParType.N; if (*Error) goto E0;

    InpParType.Theta = (MarginalDistributionType**)malloc(InpParType.c * sizeof(MarginalDistributionType*));

    *Error = NULL == InpParType.Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < InpParType.c; j++) {
        InpParType.N[j] = N[j];

        InpParType.Theta[j] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        *Error = NULL == InpParType.Theta[j]; if (*Error) goto E0;

        for (k = 0; k < InpParType.d; k++) {
            if (!strcmp(ParFamType[i], "normal")) {
                InpParType.Theta[j][k].ParFamType = pfNormal;            
                
                InpParType.Theta[j][k].Par0 = Par0[i];
                InpParType.Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "lognormal")) {
                InpParType.Theta[j][k].ParFamType = pfLognormal;

                InpParType.Theta[j][k].Par0 = Par0[i];
                InpParType.Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Weibull")) {
                InpParType.Theta[j][k].ParFamType = pfWeibull;

                InpParType.Theta[j][k].Par0 = Par0[i];
                InpParType.Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "gamma")) {
                InpParType.Theta[j][k].ParFamType = pfGamma;

                InpParType.Theta[j][k].Par0 = Par0[i];
                InpParType.Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "binomial")) {
                InpParType.Theta[j][k].ParFamType = pfBinomial;

                InpParType.Theta[j][k].Par0 = Par0[i];
                InpParType.Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Poisson")) {
                InpParType.Theta[j][k].ParFamType = pfPoisson;

                InpParType.Theta[j][k].Par0 = Par0[i];
            }
            else
            if (!strcmp(ParFamType[i], "Dirac")) {
                InpParType.Theta[j][k].ParFamType = pfDirac;

                InpParType.Theta[j][k].Par0 = Par0[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    *Error = RNGMIX(&InpParType, &OutParType);

    if (*Error) goto E0;

    *n = OutParType.n; i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (k = 0; k < OutParType.n; k++) {
            X[i] = OutParType.X[k][j]; i++;
        }
    }

E0:    if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }

        free(OutParType.X);
    }
    
    if (InpParType.N) free(InpParType.N);
    
    if (InpParType.Theta) {
        for (i = 0; i < InpParType.c; i++) {
            if (InpParType.Theta[i]) free(InpParType.Theta[i]);
        }

        free(InpParType.Theta);
    }
} /* RRNGMIX */

/* Runs REBMIX in R. */

void RREBMIX(char   **PreType,       /* Preprocessing of observations. */
             double *D,              /* Total of positive relative deviations. */
             int    *cmax,           /* Maximum number of components. */
             char   **ICType,        /* Information criterion types. */
             int    *d,              /* Number of independent random variables. */ 
             char   **VarType,       /* Types of variables. */
             char   **IniParFamType, /* Initial parametric family types. */
             int    *length_IniPar0, /* Length of IniPar0. */
             double *IniPar0,        /* Initial component parameters. */
             int    *length_IniPar1, /* Length of IniPar1. */
             double *IniPar1,        /* Initial component parameters. */
             int    *kmax,           /* Number of classes or k-nearest neighbours to be processed. */
             int    *K,              /* Number of classes or k-nearest neighbours. */
             int    *length_Iniy0,   /* Length of Iniy0. */
             double *Iniy0,          /* Initial origins. */
             int    *length_ymin,    /* Length of ymin. */
             double *ymin,           /* Minimum observations. */
             int    *length_ymax,    /* Length of ymax. */
             double *ymax,           /* Maximum observations. */
             double *b,              /* Minimum weight multiplier. */
             double *ar,             /* Acceleration rate. */
             char   **ResType,       /* Restraints type. */
             int    *n,              /* Total number of independent observations. */
             double *X,              /* Pointer to the input array X. */
             int    *k,              /* Optimal number of classes or k-nearest neighbours. */
             double *h,              /* Optimal class widths. */ 
             double *y0,             /* Optimal origins. */ 
             double *IC,             /* Optimal information criterion. */
             double *logL,           /* Log-likelihood. */
             int    *M,              /* Degrees of freedom. */ 
             int    *c,              /* Optimal number of components. */ 
             double *W,              /* Optimal component weights. */
             char   **ParFamType,    /* Parametric family types. */
             double *Par0,           /* Component parameters. */
             double *Par1,           /* Component parameters. */
             int    *all_Imax,       /* Number of iterations. */
             int    *all_c,          /* All numbers of components. */ 
             double *all_IC,         /* All information criteria. */  
             double *all_logL,       /* All log-likelihoods. */
             double *all_D,          /* All totals of positive relative deviations. */ 
             int    *Error)          /* Error code. */
{
    InputREBMIXParameterType   InpParType;
    OutputREBMIXParameterType  OutParType;
    HistoryREBMIXParameterType HisParType;
    int                        i, j, l;

    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));
    memset(&HisParType, 0, sizeof(HistoryREBMIXParameterType));

    if (!strcmp(PreType[0], "histogram")) {
        InpParType.PreType = poHistogram;
    }
    else
    if (!strcmp(PreType[0], "Parzen window")) {
        InpParType.PreType = poParzenWindow;
    }
    else
    if (!strcmp(PreType[0], "k-nearest neighbour")) {
        InpParType.PreType = poKNearestNeighbour;
    }
    else {
        *Error = 1; goto E0;
    }

    InpParType.D = *D;

    InpParType.cmax = *cmax;

    if (!strcmp(ICType[0], "AIC"))
        InpParType.ICType = icAIC; 
    else
    if (!strcmp(ICType[0], "AIC3"))
        InpParType.ICType = icAIC3;
    else
    if (!strcmp(ICType[0], "AIC4"))
        InpParType.ICType = icAIC4;
    else
    if (!strcmp(ICType[0], "AICc")) {
        InpParType.ICType = icAICc;
    }
    else
    if (!strcmp(ICType[0], "BIC"))
        InpParType.ICType = icBIC;
    else
    if (!strcmp(ICType[0], "CAIC"))
        InpParType.ICType = icCAIC;
    else
    if (!strcmp(ICType[0], "HQC"))
        InpParType.ICType = icHQC;
    else
    if (!strcmp(ICType[0], "MDL2"))
        InpParType.ICType = icMDL2;
    else
    if (!strcmp(ICType[0], "MDL5"))
        InpParType.ICType = icMDL5;
    else
    if (!strcmp(ICType[0], "AWE"))
        InpParType.ICType = icAWE;
    else
    if (!strcmp(ICType[0], "CLC"))
        InpParType.ICType = icCLC;
    else
    if (!strcmp(ICType[0], "ICL"))
        InpParType.ICType = icICL;
    else
    if (!strcmp(ICType[0], "PC"))
        InpParType.ICType = icPC;
    else
    if (!strcmp(ICType[0], "ICL-BIC"))
        InpParType.ICType = icICLBIC;
    else
    if (!strcmp(ICType[0], "D"))
        InpParType.ICType = icD;
    else
    if (!strcmp(ICType[0], "SSE"))
        InpParType.ICType = icSSE;
    else {
        *Error = 1; goto E0;
    }

    InpParType.d = *d;

    InpParType.VarType = (VariablesType_e*)malloc(InpParType.d * sizeof(VariablesType_e));

    *Error = NULL == InpParType.VarType; if (*Error) goto E0;

    for (i = 0; i < InpParType.d; i++) {
        if (!strcmp(VarType[i], "continuous")) {
            InpParType.VarType[i] = vtContinuous;
        }
        else
        if (!strcmp(VarType[i], "discrete")) {
            InpParType.VarType[i] = vtDiscrete;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    InpParType.ParFamType = (ParametricFamilyType_e*)malloc(InpParType.d * sizeof(ParametricFamilyType_e));

    *Error = NULL == InpParType.ParFamType; if (*Error) goto E0;

    for (i = 0; i < InpParType.d; i++) {
        if (!strcmp(IniParFamType[i], "normal")) {
            InpParType.ParFamType[i] = pfNormal;
        }
        else
        if (!strcmp(IniParFamType[i], "lognormal")) {
            InpParType.ParFamType[i] = pfLognormal;
        }
        else
        if (!strcmp(IniParFamType[i], "Weibull")) {
            InpParType.ParFamType[i] = pfWeibull;
        }
        else
        if (!strcmp(IniParFamType[i], "gamma")) {
            InpParType.ParFamType[i] = pfGamma;
        }
        else
        if (!strcmp(IniParFamType[i], "binomial")) {
            InpParType.ParFamType[i] = pfBinomial;
        }
        else
        if (!strcmp(IniParFamType[i], "Poisson")) {
            InpParType.ParFamType[i] = pfPoisson;
        }
        else
        if (!strcmp(IniParFamType[i], "Dirac")) {
            InpParType.ParFamType[i] = pfDirac;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    if (*length_IniPar0 > 0) {
        InpParType.Par0 = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == InpParType.Par0; if (*Error) goto E0;

        for (i = 0; i < InpParType.d; i++) {
            InpParType.Par0[i] = IniPar0[i];
        }
    }
    else {
        InpParType.Par0 = NULL;
    }

    if (*length_IniPar1 > 0) {
        InpParType.Par1 = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == InpParType.Par1; if (*Error) goto E0;

        for (i = 0; i < InpParType.d; i++) {
            InpParType.Par1[i] = IniPar1[i];
        }
    }
    else {
        InpParType.Par1 = NULL;
    }

    InpParType.kmax = *kmax;

    InpParType.K = (int*)malloc(InpParType.kmax * sizeof(int));

    *Error = NULL == InpParType.K; if (*Error) goto E0;

    for (i = 0; i < InpParType.kmax; i++) {
        InpParType.K[i] = K[i];
    }

    if (*length_Iniy0 > 0) {
        InpParType.y0 = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == InpParType.y0; if (*Error) goto E0;

        for (i = 0; i < InpParType.d; i++) {
            InpParType.y0[i] = Iniy0[i];
        }
    }
    else {
        InpParType.y0 = NULL;
    }

    if (*length_ymin > 0) {
        InpParType.ymin = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == InpParType.ymin; if (*Error) goto E0;

        for (i = 0; i < InpParType.d; i++) {
            InpParType.ymin[i] = ymin[i];
        }
    }
    else {
        InpParType.ymin = NULL;
    }

    if (*length_ymax > 0) {
        InpParType.ymax = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == InpParType.ymax; if (*Error) goto E0;

        for (i = 0; i < InpParType.d; i++) {
            InpParType.ymax[i] = ymax[i];
        }
    }
    else {
        InpParType.ymax = NULL;
    }

    InpParType.b = *b;
    InpParType.ar = *ar;

    if (!strcmp(ResType[0], "rigid")) {
        InpParType.ResType = rtRigid;
    }
    else
    if (!strcmp(ResType[0], "loose")) {
        InpParType.ResType = rtLoose;
    }
    else {
        *Error = 1; goto E0;
    }

    OutParType.n = *n;

    OutParType.X = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == OutParType.X; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        OutParType.X[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == OutParType.X[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (l = 0; l < OutParType.n; l++) {
            OutParType.X[l][j] = X[i]; i++;
        }
    }

    *Error = REBMIX(&InpParType, &OutParType, &HisParType);

     if (*Error) goto E0;

    *k = OutParType.k;

    if (OutParType.h) for (i = 0; i < InpParType.d; i++) {
        h[i] = OutParType.h[i];
    }

    if (OutParType.y0) for (i = 0; i < InpParType.d; i++) {
        y0[i] = OutParType.y0[i];
    }

    *IC = OutParType.IC;
    *logL = OutParType.logL;
    *M = OutParType.M;
    *c = OutParType.c;

    i = 0;

    for (j = 0; j < OutParType.c; j++) {
        W[j] = OutParType.W[j];

        for (l = 0; l < InpParType.d; l++) {
            switch (OutParType.Theta[j][l].ParFamType) {
            case pfNormal:
                strcpy(ParFamType[i], "normal");

                Par0[i] = OutParType.Theta[j][l].Par0;
                Par1[i] = OutParType.Theta[j][l].Par1;

                break;
            case pfLognormal:
                strcpy(ParFamType[i], "lognormal");

                Par0[i] = OutParType.Theta[j][l].Par0;
                Par1[i] = OutParType.Theta[j][l].Par1;

                break;
            case pfWeibull:
                strcpy(ParFamType[i], "Weibull");

                Par0[i] = OutParType.Theta[j][l].Par0;
                Par1[i] = OutParType.Theta[j][l].Par1;

                break;
            case pfGamma:
                strcpy(ParFamType[i], "gamma");

                Par0[i] = OutParType.Theta[j][l].Par0;
                Par1[i] = OutParType.Theta[j][l].Par1;

                break;
            case pfBinomial:
                strcpy(ParFamType[i], "binomial");

                Par0[i] = OutParType.Theta[j][l].Par0;
                Par1[i] = OutParType.Theta[j][l].Par1;

                break;
            case pfPoisson:
                strcpy(ParFamType[i], "Poisson");

                Par0[i] = OutParType.Theta[j][l].Par0;

                break;
            case pfDirac:
                strcpy(ParFamType[i], "Dirac");

                Par0[i] = OutParType.Theta[j][l].Par0;
            } 

            i++;
        }
    }

    *all_Imax = HisParType.Imax;

    for (i = 0; i < HisParType.Imax; i++) {
        all_c[i] = HisParType.c[i];
        all_IC[i] = HisParType.IC[i];
        all_logL[i] = HisParType.logL[i];
        all_D[i] = HisParType.D[i];
    }

E0: if (InpParType.save) free(InpParType.save);

    if (InpParType.ymax) free(InpParType.ymax);

    if (InpParType.ymin) free(InpParType.ymin);

    if (InpParType.y0) free(InpParType.y0);

    if (InpParType.K) free(InpParType.K);

    if (InpParType.Par1) free(InpParType.Par1);

    if (InpParType.Par0) free(InpParType.Par0);

    if (InpParType.ParFamType) free(InpParType.ParFamType);

    if (InpParType.VarType) free(InpParType.VarType);
    
    if (InpParType.open) {
        for (i = 0; i < InpParType.o; i++) {
            if (InpParType.open[i]) free(InpParType.open[i]);
        }
         
        free(InpParType.open);
    }

    if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }

    if (OutParType.Theta) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (OutParType.Theta[i]) free(OutParType.Theta[i]);
        }
         
        free(OutParType.Theta);
    }

    if (OutParType.W) free(OutParType.W);

    if (OutParType.y0) free(OutParType.y0);

    if (OutParType.h) free(OutParType.h);

    if (HisParType.D) free(HisParType.D);

    if (HisParType.logL) free(HisParType.logL);

    if (HisParType.IC) free(HisParType.IC);

    if (HisParType.c) free(HisParType.c);
} /* RREBMIX */

/* Returns k-nearest neighbour empirical densities in R. */

void RdensKNearestNeighbourXY(int    *n,     /* Total number of independent observations. */
                              double *x,     /* Pointer to the input array x. */
                              double *y,     /* Pointer to the input array y. */
                              double *p,     /* Pointer to the output array p. */
                              int    *k,     /* k-nearest neighbours. */
                              double *hx,    /* Normalizing vector. */
                              double *hy,    /* Normalizing vector. */
                              int    *Error) /* Error code. */
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, C;
    int   i, j, K, l, m, q;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1; 
    
    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * Pi * (*hx) * (*hy));

    for (i = 0; i < *n; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (x[i] - x[j]) / (*hx); Dc = R * R;
            R = (y[i] - y[j]) / (*hy); Dc += R * R;

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < K; l++) {
                if (Dc < Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != K - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = (FLOAT)sqrt(Dk[K - 1]);

        if (q >= K) R *= (FLOAT)sqrt((K + (FLOAT)1.0) / (q + (FLOAT)2.0));

        p[i] = C / (R * R);
    }

E0: if (Dk) free(Dk);
} /* RdensKNearestNeighbourXY */

/* Returns Parzen window empirical densities in R. */

void RdensParzenWindowXY(int    *n,     /* Total number of independent observations. */
                         double *x,     /* Pointer to the input array x. */
                         double *y,     /* Pointer to the input array y. */
                         double *p,     /* Pointer to the output array p. */
                         double *hx,    /* Side of the hypersquare. */
                         double *hy,    /* Side of the hypersquare. */
                         int    *Error) /* Error code. */
{
    int   i, j;
    FLOAT C, rx, ry;

    *Error = *n < 1; if (*Error) return;

    C = (FLOAT)1.0 / (*hx) / (*hy) / (*n); rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    for (i = 0; i < *n; i++) {
        p[i] = (FLOAT)0.0;
    }

    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if (((FLOAT)fabs(x[j] - x[i]) <= rx) && ((FLOAT)fabs(y[j] - y[i]) <= ry)) {
                p[i] += C; if (i != j) p[j] += C;
            }
        }
    }
} /* RdensParzenWindowXY */

/* Returns histogram empirical densities in R. */

void RdensHistogramXY(int    *k,     /* Total number of bins. */ 
                      int    *n,     /* Total number of independent observations. */
                      double *x,     /* Pointer to the input array x. */
                      double *y,     /* Pointer to the input array y. */
                      double *p,     /* Pointer to the output array p. */
                      double *x0,    /* Origin. */ 
                      double *y0,    /* Origin. */ 
                      double *hx,    /* Side of the hypersquare. */
                      double *hy,    /* Side of the hypersquare. */
                      int    *cx,    /* If x discrete then cx = 1 else cx = 0. */
                      int    *cy,    /* If y discrete then cy = 1 else cy = 0. */
                      char   **px,   /* Parametric family type. */
                      char   **py,   /* Parametric family type. */
                      int    *Error) /* Error code. */
{
    int                    i, j;
    ParametricFamilyType_e pdfx, pdfy;
    FLOAT                  C, rx, ry;

    *Error = *n < 1; if (*Error) return;

    if (!strcmp(px[0], "normal")) {
        pdfx = pfNormal;
    }
    else
    if (!strcmp(px[0], "lognormal")) {
        pdfx = pfLognormal;
    }
    else
    if (!strcmp(px[0], "Weibull")) {
        pdfx = pfWeibull;
    }
    else
    if (!strcmp(px[0], "gamma")) {
        pdfx = pfGamma;
    }
    else
    if (!strcmp(px[0], "binomial")) {
        pdfx = pfBinomial;
    }
    else
    if (!strcmp(px[0], "Poisson")) {
        pdfx = pfPoisson;
    }
    else
    if (!strcmp(px[0], "Dirac")) {
        pdfx = pfDirac;
    }
    else {
      *Error = 1; return;
    }

    if (!strcmp(py[0], "normal")) {
        pdfy = pfNormal;
    }
    else
    if (!strcmp(py[0], "lognormal")) {
        pdfy = pfLognormal;
    }
    else
    if (!strcmp(py[0], "Weibull")) {
        pdfy = pfWeibull;
    }
    else
    if (!strcmp(py[0], "gamma")) {
        pdfy = pfGamma;
    }
    else
    if (!strcmp(py[0], "binomial")) {
        pdfy = pfBinomial;
    }
    else
    if (!strcmp(py[0], "Poisson")) {
        pdfy = pfPoisson;
    }
    else
    if (!strcmp(py[0], "Dirac")) {
        pdfy = pfDirac;
    }
    else {
      *Error = 1; return;
    }

    C = (FLOAT)1.0 / (*hx) / (*hy) / (*n); rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5); 
        
        x[*k] = (*x0) + j * (*hx);

        switch (pdfx) {
        default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (x[*k] <= FLOAT_MIN) x[*k] += (*hx);
        }

        j = (int)floor((y[i] - (*y0)) / (*hy) + (FLOAT)0.5); 

        y[*k] = (*y0) + j * (*hy);

        switch (pdfy) {
        default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (y[*k] <= FLOAT_MIN) y[*k] += (*hy);
        }

        for (j = 0; j < *k; j++) {
            if (((FLOAT)fabs(x[j] - x[*k]) > rx) || ((FLOAT)fabs(y[j] - y[*k]) > ry)) goto S0;
                
            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++; 
S1:;}
} /* RdensHistogramXY */

/* Returns k-nearest neighbour empirical densities in R. */

void RdensKNearestNeighbourX(int    *n,     /* Total number of independent observations. */
                             double *x,     /* Pointer to the input array x. */
                             double *p,     /* Pointer to the output array p. */
                             int    *k,     /* k-nearest neighbours. */
                             double *hx,    /* Normalizing vector. */
                             int    *Error) /* Error code. */
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, C;
    int   i, j, K, l, m, q;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1; 
    
    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * (FLOAT)2.0 * (*hx));

    for (i = 0; i < *n; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (FLOAT)fabs((x[i] - x[j]) / (*hx)); Dc = R;

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < K; l++) {
                if (Dc < Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != K - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = Dk[K - 1];

        if (q >= K) R *= (FLOAT)((K + (FLOAT)1.0) / (q + (FLOAT)2.0));

        p[i] = C / R;
    }

E0: if (Dk) free(Dk);
} /* RdensKNearestNeighbourX */

/* Returns Parzen window empirical densities in R. */

void RdensParzenWindowX(int    *n,     /* Total number of independent observations. */
                        double *x,     /* Pointer to the input array x. */
                        double *p,     /* Pointer to the output array p. */
                        double *hx,    /* Side of the hypersquare. */
                        int    *Error) /* Error code. */
{
    int   i, j;
    FLOAT C, rx;

    *Error = *n < 1; if (*Error) return;

    C = (FLOAT)1.0 / (*hx) / (*n); rx = (FLOAT)0.5 * (*hx);

    for (i = 0; i < *n; i++) {
        p[i] = (FLOAT)0.0;
    }
    
    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if ((FLOAT)fabs(x[j] - x[i]) <= rx) {
                p[i] += C; if (i != j) p[j] += C;
            }
        }
    }
} /* RdensParzenWindowX */

/* Returns histogram empirical densities in R. */

void RdensHistogramX(int    *k,     /* Total number of bins. */ 
                     int    *n,     /* Total number of independent observations. */
                     double *x,     /* Pointer to the input array x. */
                     double *p,     /* Pointer to the output array p. */
                     double *x0,    /* Origin. */ 
                     double *hx,    /* Side of the hypersquare. */
                     int    *cx,    /* If x discrete then cx = 1 else cx = 0. */
                     char   **px,   /* Parametric family type. */
                     int    *Error) /* Error code. */
{
    int                    i, j;
    ParametricFamilyType_e pdfx;
    FLOAT                  C, rx;

    *Error = *n < 1; if (*Error) return;

    if (!strcmp(px[0], "normal")) {
        pdfx = pfNormal;
    }
    else
    if (!strcmp(px[0], "lognormal")) {
        pdfx = pfLognormal;
    }
    else
    if (!strcmp(px[0], "Weibull")) {
        pdfx = pfWeibull;
    }
    else
    if (!strcmp(px[0], "gamma")) {
        pdfx = pfGamma;
    }
    else
    if (!strcmp(px[0], "binomial")) {
        pdfx = pfBinomial;
    }
    else
    if (!strcmp(px[0], "Poisson")) {
        pdfx = pfPoisson;
    }
    else
    if (!strcmp(px[0], "Dirac")) {
        pdfx = pfDirac;
    }
    else {
      *Error = 1; return;
    }
    
    C = (FLOAT)1.0 / (*hx) / (*n); rx = (FLOAT)0.5 * (*hx);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5); 
        
        x[*k] = (*x0) + j * (*hx);

        switch (pdfx) {
        default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (x[*k] <= FLOAT_MIN) x[*k] += (*hx);
        }

        for (j = 0; j < *k; j++) {
            if ((FLOAT)fabs(x[j] - x[*k]) > rx) goto S0;
                
            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++; 
S1:;}
} /* RdensHistogramX */

/* Returns classified observations in R. */

void RCLSMIX(int    *n,           /* Total number of independent observations. */
             double *X,           /* Pointer to the input array X. */
             int    *s,           /* Number of classes. */ 
             int    *o,           /* Number of input REBMIX objects. */    
             int    *d,           /* Number of independent random variables in REBMIX objects. */
             int    *c,           /* Number of components. */
             double *W,           /* Component weights. */
             char   **ParFamType, /* Parametric family types. */
             double *Par0,        /* Component parameters. */
             double *Par1,        /* Component parameters. */
             double *P,           /* Prior probabilities. */  
             int    *Z,           /* Pointer to the output array Z. */
             int    *Error)       /* Error code. */
{
    int                      **C; 
    FLOAT                    ***Q = NULL;
    FLOAT                    *Y = NULL;
    MarginalDistributionType ****Theta = NULL; 
    FLOAT                    CmpDist, MixDist, MaxMixDist;
    int                      i, j, k, l, m;

    C = (int**)malloc(*s * sizeof(int*));

    *Error = NULL == C; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        C[j] = (int*)malloc(*o * sizeof(int));

        *Error = NULL == C[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            C[j][k] = c[i]; i++;
        }
    }

    Q = (FLOAT***)malloc(*s * sizeof(FLOAT**));

    *Error = NULL == Q; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        Q[j] = (FLOAT**)malloc(*o * sizeof(FLOAT*));

        *Error = NULL == Q[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            Q[j][k] = (FLOAT*)malloc(C[j][k] * sizeof(FLOAT));

            *Error = NULL == Q[j][k]; if (*Error) goto E0;

            for (l = 0; l < C[j][k]; l++) {
                Q[j][k][l] = W[i]; i++;
            }
        }
    }

    Theta = (MarginalDistributionType****)malloc(*s * sizeof(MarginalDistributionType***));

    *Error = NULL == Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        Theta[j] = (MarginalDistributionType***)malloc(*o * sizeof(MarginalDistributionType**));

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            Theta[j][k] = (MarginalDistributionType**)malloc(C[j][k] * sizeof(MarginalDistributionType*));

            *Error = NULL == Theta[j][k]; if (*Error) goto E0;

            for (l = 0; l < C[j][k]; l++) {
                Theta[j][k][l] = (MarginalDistributionType*)malloc(d[k] * sizeof(MarginalDistributionType));

                *Error = NULL == Theta[j][k][l]; if (*Error) goto E0;

                for (m = 0; m < d[k]; m++) {
                    if (!strcmp(ParFamType[i], "normal")) {
                        Theta[j][k][l][m].ParFamType = pfNormal;

                        Theta[j][k][l][m].Par0 = Par0[i];
                        Theta[j][k][l][m].Par1 = Par1[i];
                    }
                    else
                    if (!strcmp(ParFamType[i], "lognormal")) {
                        Theta[j][k][l][m].ParFamType = pfLognormal;

                        Theta[j][k][l][m].Par0 = Par0[i];
                        Theta[j][k][l][m].Par1 = Par1[i];
                    }
                    else
                    if (!strcmp(ParFamType[i], "Weibull")) {
                        Theta[j][k][l][m].ParFamType = pfWeibull;

                        Theta[j][k][l][m].Par0 = Par0[i];
                        Theta[j][k][l][m].Par1 = Par1[i];
                    }
                    else
                    if (!strcmp(ParFamType[i], "gamma")) {
                        Theta[j][k][l][m].ParFamType = pfGamma;

                        Theta[j][k][l][m].Par0 = Par0[i];
                        Theta[j][k][l][m].Par1 = Par1[i];
                    }
                    else
                    if (!strcmp(ParFamType[i], "binomial")) {
                        Theta[j][k][l][m].ParFamType = pfBinomial;

                        Theta[j][k][l][m].Par0 = Par0[i];
                        Theta[j][k][l][m].Par1 = Par1[i];
                    }
                    else
                    if (!strcmp(ParFamType[i], "Poisson")) {
                        Theta[j][k][l][m].ParFamType = pfPoisson;

                        Theta[j][k][l][m].Par0 = Par0[i];
                    }
                    else
                    if (!strcmp(ParFamType[i], "Dirac")) {
                        Theta[j][k][l][m].ParFamType = pfDirac;

                        Theta[j][k][l][m].Par0 = Par0[i];
                    }
                    else {
                        *Error = 1; goto E0;
                    }

                    i++;
                }
            }
        }
    }

    i = d[0]; for (j = 1; j < *o; j++) if (d[j] > i) i = d[j]; 

    Y = (FLOAT*)malloc(i * sizeof(FLOAT));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < *n; i++) {
        Z[i] = 0; MaxMixDist = (FLOAT)0.0; 
         
        for (j = 0; j < *s; j++) {
            k = 0; MixDist = (FLOAT)1.0;
            
            for (l = 0; l < *o; l++) {
                for (m = 0; m < d[l]; m++) {
                    Y[m] = X[i + (*n) * (m + k)];
                }

                *Error = MixtureDist(d[l], Y, C[j][l], Q[j][l], Theta[j][l], &CmpDist, 0);

                if (*Error) goto E0;

                k += d[l]; MixDist *= CmpDist; 
            }

            MixDist *= P[j];

            if (MixDist > MaxMixDist) {
                Z[i] = j; MaxMixDist = MixDist;
            }
        }
    }

E0: if (Y) free(Y);

    if (Theta) {
        for (i = 0; i < *s; i++) {
            if (Theta[i]) {
                for (j = 0; j < *o; j++) {
                    if (Theta[i][j]) {
                        for (k = 0; k < C[i][j]; k++) {
                            if (Theta[i][j][k]) free(Theta[i][j][k]);
                        }
                        
                        free(Theta[i][j]);
                    }
                }

                free(Theta[i]);
            }
        }

        free(Theta);
    }

    if (Q) {
        for (i = 0; i < *s; i++) {
            if (Q[i]) {
                for (j = 0; j < *o; j++) if (Q[i][j]) free(Q[i][j]);
                
                free(Q[i]);
            }
        }

        free(Q);
    } 

    if (C) {
        for (i = 0; i < *s; i++) {
            if (C[i]) free(C[i]);
        }

        free(C);
    }
} /* RCLSMIX */

/* Returns clustered observations in R. */

void RCLRMIX(int    *n,           /* Total number of independent observations. */
             double *X,           /* Pointer to the input array X. */
             int    *d,           /* Number of independent random variables. */
             int    *c,           /* Number of components. */ 
             double *W,           /* Component weights. */
             char   **ParFamType, /* Parametric family types. */
             double *Par0,        /* Component parameters. */
             double *Par1,        /* Component parameters. */
             double *Z,           /* Pointer to the output array Z. */
             int    *Error)       /* Error code. */
{
    FLOAT                    *Y = NULL;
    MarginalDistributionType **Theta = NULL; 
    FLOAT                    CmpDist, MaxCmpDist;
    int                      i, j, k;

    Theta = (MarginalDistributionType**)malloc(*c * sizeof(MarginalDistributionType*));

    *Error = NULL == Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *c; j++) {
        Theta[j] = (MarginalDistributionType*)malloc(*d * sizeof(MarginalDistributionType));

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        for (k = 0; k < *d; k++) {
            if (!strcmp(ParFamType[i], "normal")) {
                Theta[j][k].ParFamType = pfNormal;

                Theta[j][k].Par0 = Par0[i];
                Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "lognormal")) {
                Theta[j][k].ParFamType = pfLognormal;

                Theta[j][k].Par0 = Par0[i];
                Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Weibull")) {
                Theta[j][k].ParFamType = pfWeibull;

                Theta[j][k].Par0 = Par0[i];
                Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "gamma")) {
                Theta[j][k].ParFamType = pfGamma;

                Theta[j][k].Par0 = Par0[i];
                Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "binomial")) {
                Theta[j][k].ParFamType = pfBinomial;

                Theta[j][k].Par0 = Par0[i];
                Theta[j][k].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Poisson")) {
                Theta[j][k].ParFamType = pfPoisson;

                Theta[j][k].Par0 = Par0[i];
            }
            else
            if (!strcmp(ParFamType[i], "Dirac")) {
                Theta[j][k].ParFamType = pfDirac;

                Theta[j][k].Par0 = Par0[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    Y = (FLOAT*)malloc(*d * sizeof(FLOAT));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *d; j++) {
            Y[j] = X[i + (*n) * j];
        }

        Z[i] = (FLOAT)0.0; MaxCmpDist = (FLOAT)0.0;
         
        for (j = 0; j < *c; j++) {
            *Error = ComponentDist(*d, Y, Theta[j], &CmpDist, 0);

            if (*Error) goto E0;

            CmpDist *= W[j];

            if (CmpDist > MaxCmpDist) {
                Z[i] = (FLOAT)j; MaxCmpDist = CmpDist;
            }
        }
    }

E0: if (Y) free(Y);

    if (Theta) {
        for (i = 0; i < *c; i++) {
            if (Theta[i]) free(Theta[i]);
        }

        free(Theta);
    }
} /* RCLRMIX */

void RPreprocessingKNN(int    *k,     /* k-nearest neighbours. */
                       double *h,     /* Normalizing vector. */
                       int    *n,     /* Total number of independent observations. */
                       int    *d,     /* Number of independent random variables. */ 
                       double *x,     /* Pointer to the input array x. */
                       double *y,     /* Pointer to the output array y. */
                       int    *Error) /* Error code. */
{
    InputREBMIXParameterType  InpParType;
    OutputREBMIXParameterType OutParType;
    int                       i, j, l;

    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));

    OutParType.n = *n;
    InpParType.d = *d;

    OutParType.X = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == OutParType.X; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        OutParType.X[i] = (FLOAT*)malloc((InpParType.d + 3) * sizeof(FLOAT));

        *Error = NULL == OutParType.X[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (l = 0; l < OutParType.n; l++) {
            OutParType.X[l][j] = x[i]; i++;
        }
    }

    OutParType.k = *k;

    *Error = PreprocessingKNN(OutParType.k, h, OutParType.n, InpParType.d, OutParType.X);

    if (*Error) goto E0;

    i = 0;
    
    for (j = 0; j < InpParType.d + 3; j++) {
        for (l = 0; l < OutParType.n; l++) {
            y[i] = OutParType.X[l][j]; i++;
        }
    }

E0: if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }
} /* RPreprocessingKNN */

void RPreprocessingPW(double *h,     /* Sides of the hypersquare. */
                      int    *n,     /* Total number of independent observations. */
                      int    *d,     /* Number of independent random variables. */ 
                      double *x,     /* Pointer to the input array x. */
                      double *y,     /* Pointer to the output array y. */
                      int    *Error) /* Error code. */
{
    InputREBMIXParameterType  InpParType;
    OutputREBMIXParameterType OutParType;
    int                       i, j, l;

    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));

    OutParType.n = *n;
    InpParType.d = *d;

    OutParType.X = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == OutParType.X; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        OutParType.X[i] = (FLOAT*)malloc((InpParType.d + 2) * sizeof(FLOAT));

        *Error = NULL == OutParType.X[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (l = 0; l < OutParType.n; l++) {
            OutParType.X[l][j] = x[i]; i++;
        }
    }

    *Error = PreprocessingPW(h, OutParType.n, InpParType.d, OutParType.X); 

    if (*Error) goto E0;
    
    i = 0;
    
    for (j = 0; j < InpParType.d + 2; j++) {
        for (l = 0; l < OutParType.n; l++) {
            y[i] = OutParType.X[l][j]; i++;
        }
    }

E0: if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }
} /* RPreprocessingPW */

void RPreprocessingH(double *h,           /* Sides of the hypersquare. */
                     double *y0,          /* Origins. */
                     char   **ParFamType, /* Parametric family types. */
                     int    *k,           /* Total number of bins. */
                     int    *n,           /* Total number of independent observations. */
                     int    *d,           /* Number of independent random variables. */ 
                     double *x,           /* Pointer to the input array x. */
                     double *y,           /* Pointer to the output array y. */
                     int    *Error)       /* Error code. */
{
    InputREBMIXParameterType  InpParType;
    OutputREBMIXParameterType OutParType;
    FLOAT                     **Y = NULL;
    int                       i, j, l;

    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));

    InpParType.d = *d;

    InpParType.ParFamType = (ParametricFamilyType_e*)malloc(InpParType.d * sizeof(ParametricFamilyType_e));

    *Error = NULL == InpParType.ParFamType; if (*Error) goto E0;

    for (i = 0; i < InpParType.d; i++) {
        if (!strcmp(ParFamType[i], "normal")) {
            InpParType.ParFamType[i] = pfNormal;
        }
        else
        if (!strcmp(ParFamType[i], "lognormal")) {
            InpParType.ParFamType[i] = pfLognormal;
        }
        else
        if (!strcmp(ParFamType[i], "Weibull")) {
            InpParType.ParFamType[i] = pfWeibull;
        }
        else
        if (!strcmp(ParFamType[i], "gamma")) {
            InpParType.ParFamType[i] = pfGamma;
        }
        else
        if (!strcmp(ParFamType[i], "binomial")) {
            InpParType.ParFamType[i] = pfBinomial;
        }
        else
        if (!strcmp(ParFamType[i], "Poisson")) {
            InpParType.ParFamType[i] = pfPoisson;
        }
        else
        if (!strcmp(ParFamType[i], "Dirac")) {
            InpParType.ParFamType[i] = pfDirac;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    OutParType.n = *n;

    OutParType.X = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == OutParType.X; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        OutParType.X[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == OutParType.X[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (l = 0; l < OutParType.n; l++) {
            OutParType.X[l][j] = x[i]; i++;
        }
    }
    
    Y = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType.d + 1) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    OutParType.k = *k;

    *Error = PreprocessingH(h, y0, InpParType.ParFamType, &OutParType.k, OutParType.n, InpParType.d, OutParType.X, Y);

    if (*Error) goto E0;

    i = 0;
    
    for (j = 0; j < InpParType.d + 1; j++) {
        for (l = 0; l < OutParType.k; l++) {
            y[i] = Y[l][j]; i++;
        }
    }

    *k = OutParType.k;

E0: if (Y) {
        for (i = 0; i < OutParType.n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }

    if (InpParType.ParFamType) free(InpParType.ParFamType);
} /* RPreprocessingH */

void RInformationCriterionKNN(int    *k,           /* Total number of bins. */
                              double *h,           /* Sides of the hypersquare. */
                              int    *n,           /* Total number of independent observations. */
                              int    *d,           /* Number of independent random variables. */ 
                              double *x,           /* Pointer to the input array x. */
                              char   **ICType,     /* Information criterion type. */ 
                              int    *c,           /* Number of components. */
                              double *W,           /* Component weights. */
                              char   **ParFamType, /* Parametric family types. */
                              double *Par0,        /* Component parameters. */
                              double *Par1,        /* Component parameters. */
                              double *IC,          /* Information criterion. */
                              double *logL,        /* log-likelihood. */
                              int    *M,           /* Degrees of freedom. */
                              double *D,           /* Total of positive relative deviations. */
                              int    *Error)       /* Error code. */
{
    InputREBMIXParameterType  InpParType;
    OutputREBMIXParameterType OutParType;
    int                       i, j, l;
    
    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));

    OutParType.n = *n;
    InpParType.d = *d;

    OutParType.X = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == OutParType.X; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        OutParType.X[i] = (FLOAT*)malloc((InpParType.d + 3) * sizeof(FLOAT));

        *Error = NULL == OutParType.X[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (l = 0; l < OutParType.n; l++) {
            OutParType.X[l][j] = x[i]; i++;
        }
    }

    OutParType.k = *k;

    *Error = PreprocessingKNN(OutParType.k, h, OutParType.n, InpParType.d, OutParType.X);

    if (*Error) goto E0;

    if (!strcmp(ICType[0], "AIC"))
        InpParType.ICType = icAIC;
    else
    if (!strcmp(ICType[0], "AIC3"))
        InpParType.ICType = icAIC3;
    else
    if (!strcmp(ICType[0], "AIC4"))
        InpParType.ICType = icAIC4;
    else
    if (!strcmp(ICType[0], "AICc")) {
        InpParType.ICType = icAICc;
    }
    else
    if (!strcmp(ICType[0], "BIC"))
        InpParType.ICType = icBIC;
    else
    if (!strcmp(ICType[0], "CAIC"))
        InpParType.ICType = icCAIC;
    else
    if (!strcmp(ICType[0], "HQC"))
        InpParType.ICType = icHQC;
    else
    if (!strcmp(ICType[0], "MDL2"))
        InpParType.ICType = icMDL2;
    else
    if (!strcmp(ICType[0], "MDL5"))
        InpParType.ICType = icMDL5;
    else
    if (!strcmp(ICType[0], "AWE"))
        InpParType.ICType = icAWE;
    else
    if (!strcmp(ICType[0], "CLC"))
        InpParType.ICType = icCLC;
    else
    if (!strcmp(ICType[0], "ICL"))
        InpParType.ICType = icICL;
    else
    if (!strcmp(ICType[0], "PC"))
        InpParType.ICType = icPC;
    else
    if (!strcmp(ICType[0], "ICL-BIC"))
        InpParType.ICType = icICLBIC;
    else
    if (!strcmp(ICType[0], "D"))
        InpParType.ICType = icD;
    else
    if (!strcmp(ICType[0], "SSE"))
        InpParType.ICType = icSSE;
    else {
        *Error = 1; goto E0;
    }

    OutParType.c = *c;

    OutParType.W = W;

    OutParType.Theta = (MarginalDistributionType**)malloc(OutParType.c * sizeof(MarginalDistributionType*));

    *Error = NULL == OutParType.Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < OutParType.c; j++) {
        OutParType.Theta[j] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        *Error = NULL == OutParType.Theta[j]; if (*Error) goto E0;

        for (l = 0; l < InpParType.d; l++) {
            if (!strcmp(ParFamType[i], "normal")) {
                OutParType.Theta[j][l].ParFamType = pfNormal;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "lognormal")) {
                OutParType.Theta[j][l].ParFamType = pfLognormal;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Weibull")) {
                OutParType.Theta[j][l].ParFamType = pfWeibull;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "gamma")) {
                OutParType.Theta[j][l].ParFamType = pfGamma;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "binomial")) {
                OutParType.Theta[j][l].ParFamType = pfBinomial;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Poisson")) {
                OutParType.Theta[j][l].ParFamType = pfPoisson;

                OutParType.Theta[j][l].Par0 = Par0[i];
            }
            else
            if (!strcmp(ParFamType[i], "Dirac")) {
                OutParType.Theta[j][l].ParFamType = pfDirac;

                OutParType.Theta[j][l].Par0 = Par0[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    *Error = InformationCriterionKNN(InpParType.ICType, 
                                     OutParType.k, 
                                     OutParType.n, 
                                     InpParType.d, 
                                     OutParType.X, 
                                     OutParType.c, 
                                     OutParType.W, 
                                     OutParType.Theta, 
                                     IC, 
                                     logL, 
                                     M,
                                     D);

    if (*Error) goto E0;

E0:    if (OutParType.Theta) {
        for (i = 0; i < OutParType.c; i++) {
            if (OutParType.Theta[i]) free(OutParType.Theta[i]);
        }
         
        free(OutParType.Theta);
    }
    
    if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }
} /* RInformationCriterionKNN */

void RInformationCriterionPW(double *h,           /* Sides of the hypersquare. */
                             int    *n,           /* Total number of independent observations. */
                             int    *d,           /* Number of independent random variables. */ 
                             double *x,           /* Pointer to the input array x. */
                             char   **ICType,     /* Information criterion type. */ 
                             int    *c,           /* Number of components. */
                             double *W,           /* Component weights. */
                             char   **ParFamType, /* Parametric family types. */
                             double *Par0,        /* Component parameters. */
                             double *Par1,        /* Component parameters. */
                             double *IC,          /* Information criterion. */
                             double *logL,        /* log-likelihood. */
                             int    *M,           /* Degrees of freedom. */
                             double *D,           /* Total of positive relative deviations. */
                             int    *Error)       /* Error code. */
{
    InputREBMIXParameterType  InpParType;
    OutputREBMIXParameterType OutParType;
    FLOAT                     V;
    int                       i, j, l;
    
    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));

    OutParType.n = *n;
    InpParType.d = *d;

    OutParType.X = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == OutParType.X; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        OutParType.X[i] = (FLOAT*)malloc((InpParType.d + 2) * sizeof(FLOAT));

        *Error = NULL == OutParType.X[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (l = 0; l < OutParType.n; l++) {
            OutParType.X[l][j] = x[i]; i++;
        }
    }

    *Error = PreprocessingPW(h, OutParType.n, InpParType.d, OutParType.X); 

    if (*Error) goto E0;

    if (!strcmp(ICType[0], "AIC"))
        InpParType.ICType = icAIC;
    else
    if (!strcmp(ICType[0], "AIC3"))
        InpParType.ICType = icAIC3;
    else
    if (!strcmp(ICType[0], "AIC4"))
        InpParType.ICType = icAIC4;
    else
    if (!strcmp(ICType[0], "AICc")) {
        InpParType.ICType = icAICc;
    }
    else
    if (!strcmp(ICType[0], "BIC"))
        InpParType.ICType = icBIC;
    else
    if (!strcmp(ICType[0], "CAIC"))
        InpParType.ICType = icCAIC;
    else
    if (!strcmp(ICType[0], "HQC"))
        InpParType.ICType = icHQC;
    else
    if (!strcmp(ICType[0], "MDL2"))
        InpParType.ICType = icMDL2;
    else
    if (!strcmp(ICType[0], "MDL5"))
        InpParType.ICType = icMDL5;
    else
    if (!strcmp(ICType[0], "AWE"))
        InpParType.ICType = icAWE;
    else
    if (!strcmp(ICType[0], "CLC"))
        InpParType.ICType = icCLC;
    else
    if (!strcmp(ICType[0], "ICL"))
        InpParType.ICType = icICL;
    else
    if (!strcmp(ICType[0], "PC"))
        InpParType.ICType = icPC;
    else
    if (!strcmp(ICType[0], "ICL-BIC"))
        InpParType.ICType = icICLBIC;
    else
    if (!strcmp(ICType[0], "D"))
        InpParType.ICType = icD;
    else
    if (!strcmp(ICType[0], "SSE"))
        InpParType.ICType = icSSE;
    else {
        *Error = 1; goto E0;
    }

    OutParType.c = *c;

    OutParType.W = W;

    OutParType.Theta = (MarginalDistributionType**)malloc(OutParType.c * sizeof(MarginalDistributionType*));

    *Error = NULL == OutParType.Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < OutParType.c; j++) {
        OutParType.Theta[j] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        *Error = NULL == OutParType.Theta[j]; if (*Error) goto E0;

        for (l = 0; l < InpParType.d; l++) {
            if (!strcmp(ParFamType[i], "normal")) {
                OutParType.Theta[j][l].ParFamType = pfNormal;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "lognormal")) {
                OutParType.Theta[j][l].ParFamType = pfLognormal;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Weibull")) {
                OutParType.Theta[j][l].ParFamType = pfWeibull;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "gamma")) {
                OutParType.Theta[j][l].ParFamType = pfGamma;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "binomial")) {
                OutParType.Theta[j][l].ParFamType = pfBinomial;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Poisson")) {
                OutParType.Theta[j][l].ParFamType = pfPoisson;

                OutParType.Theta[j][l].Par0 = Par0[i];
            }
            else
            if (!strcmp(ParFamType[i], "Dirac")) {
                OutParType.Theta[j][l].ParFamType = pfDirac;

                OutParType.Theta[j][l].Par0 = Par0[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    V = (FLOAT)1.0;

    for (i = 0; i < InpParType.d; i++) {
        V *= h[i];
    }

    *Error = InformationCriterionPW(InpParType.ICType, 
                                    V, 
                                    OutParType.n, 
                                    InpParType.d, 
                                    OutParType.X, 
                                    OutParType.c, 
                                    OutParType.W, 
                                    OutParType.Theta, 
                                    IC, 
                                    logL, 
                                    M,
                                    D);

    if (*Error) goto E0;

E0:    if (OutParType.Theta) {
        for (i = 0; i < OutParType.c; i++) {
            if (OutParType.Theta[i]) free(OutParType.Theta[i]);
        }
         
        free(OutParType.Theta);
    }
    
    if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }
} /* RInformationCriterionPW */

void RInformationCriterionH(double *h,              /* Sides of the hypersquare. */
                            double *y0,             /* Origins. */
                            char   **IniParFamType, /* Initial parametric family types. */
                            int    *k,              /* Total number of bins. */
                            int    *n,              /* Total number of independent observations. */
                            int    *d,              /* Number of independent random variables. */ 
                            double *x,              /* Pointer to the input array x. */
                            char   **ICType,        /* Information criterion type. */ 
                            int    *c,              /* Number of components. */
                            double *W,              /* Component weights. */
                            char   **ParFamType,    /* Parametric family types. */
                            double *Par0,           /* Component parameters. */
                            double *Par1,           /* Component parameters. */
                            double *IC,             /* Information criterion. */
                            double *logL,           /* log-likelihood. */
                            int    *M,              /* Degrees of freedom. */
                            double *D,              /* Total of positive relative deviations. */
                            int    *Error)          /* Error code. */
{
    InputREBMIXParameterType  InpParType;
    OutputREBMIXParameterType OutParType;
    FLOAT                     **Y = NULL;
    FLOAT                     V;
    int                       i, j, l;
    
    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));
    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));

    InpParType.d = *d;

    InpParType.ParFamType = (ParametricFamilyType_e*)malloc(InpParType.d * sizeof(ParametricFamilyType_e));

    *Error = NULL == InpParType.ParFamType; if (*Error) goto E0;

    for (i = 0; i < InpParType.d; i++) {
        if (!strcmp(IniParFamType[i], "normal")) {
            InpParType.ParFamType[i] = pfNormal;
        }
        else
        if (!strcmp(IniParFamType[i], "lognormal")) {
            InpParType.ParFamType[i] = pfLognormal;
        }
        else
        if (!strcmp(IniParFamType[i], "Weibull")) {
            InpParType.ParFamType[i] = pfWeibull;
        }
        else
        if (!strcmp(IniParFamType[i], "gamma")) {
            InpParType.ParFamType[i] = pfGamma;
        }
        else
        if (!strcmp(IniParFamType[i], "binomial")) {
            InpParType.ParFamType[i] = pfBinomial;
        }
        else
        if (!strcmp(IniParFamType[i], "Poisson")) {
            InpParType.ParFamType[i] = pfPoisson;
        }
        else
        if (!strcmp(IniParFamType[i], "Dirac")) {
            InpParType.ParFamType[i] = pfDirac;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    OutParType.n = *n;

    OutParType.X = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == OutParType.X; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        OutParType.X[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        *Error = NULL == OutParType.X[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < InpParType.d; j++) {
        for (l = 0; l < OutParType.n; l++) {
            OutParType.X[l][j] = x[i]; i++;
        }
    }
    
    Y = (FLOAT**)malloc(OutParType.n * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < OutParType.n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType.d + 1) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    OutParType.k = *k;

    *Error = PreprocessingH(h, y0, InpParType.ParFamType, &OutParType.k, OutParType.n, InpParType.d, OutParType.X, Y);

    if (*Error) goto E0;

    *k = OutParType.k;

    if (!strcmp(ICType[0], "AIC"))
        InpParType.ICType = icAIC;
    else
    if (!strcmp(ICType[0], "AIC3"))
        InpParType.ICType = icAIC3;
    else
    if (!strcmp(ICType[0], "AIC4"))
        InpParType.ICType = icAIC4;
    else
    if (!strcmp(ICType[0], "AICc")) {
        InpParType.ICType = icAICc;
    }
    else
    if (!strcmp(ICType[0], "BIC"))
        InpParType.ICType = icBIC;
    else
    if (!strcmp(ICType[0], "CAIC"))
        InpParType.ICType = icCAIC;
    else
    if (!strcmp(ICType[0], "HQC"))
        InpParType.ICType = icHQC;
    else
    if (!strcmp(ICType[0], "MDL2"))
        InpParType.ICType = icMDL2;
    else
    if (!strcmp(ICType[0], "MDL5"))
        InpParType.ICType = icMDL5;
    else
    if (!strcmp(ICType[0], "AWE"))
        InpParType.ICType = icAWE;
    else
    if (!strcmp(ICType[0], "CLC"))
        InpParType.ICType = icCLC;
    else
    if (!strcmp(ICType[0], "ICL"))
        InpParType.ICType = icICL;
    else
    if (!strcmp(ICType[0], "PC"))
        InpParType.ICType = icPC;
    else
    if (!strcmp(ICType[0], "ICL-BIC"))
        InpParType.ICType = icICLBIC;
    else
    if (!strcmp(ICType[0], "D"))
        InpParType.ICType = icD;
    else
    if (!strcmp(ICType[0], "SSE"))
        InpParType.ICType = icSSE;
    else {
        *Error = 1; goto E0;
    }

    OutParType.c = *c;

    OutParType.W = W;

    OutParType.Theta = (MarginalDistributionType**)malloc(OutParType.c * sizeof(MarginalDistributionType*));

    *Error = NULL == OutParType.Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < OutParType.c; j++) {
        OutParType.Theta[j] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        *Error = NULL == OutParType.Theta[j]; if (*Error) goto E0;

        for (l = 0; l < InpParType.d; l++) {
            if (!strcmp(ParFamType[i], "normal")) {
                OutParType.Theta[j][l].ParFamType = pfNormal;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "lognormal")) {
                OutParType.Theta[j][l].ParFamType = pfLognormal;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Weibull")) {
                OutParType.Theta[j][l].ParFamType = pfWeibull;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "gamma")) {
                OutParType.Theta[j][l].ParFamType = pfGamma;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "binomial")) {
                OutParType.Theta[j][l].ParFamType = pfBinomial;

                OutParType.Theta[j][l].Par0 = Par0[i];
                OutParType.Theta[j][l].Par1 = Par1[i];
            }
            else
            if (!strcmp(ParFamType[i], "Poisson")) {
                OutParType.Theta[j][l].ParFamType = pfPoisson;

                OutParType.Theta[j][l].Par0 = Par0[i];
            }
            else
            if (!strcmp(ParFamType[i], "Dirac")) {
                OutParType.Theta[j][l].ParFamType = pfDirac;

                OutParType.Theta[j][l].Par0 = Par0[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    V = (FLOAT)1.0;

    for (i = 0; i < InpParType.d; i++) {
        V *= h[i];
    }

    *Error = InformationCriterionH(InpParType.ICType, 
                                   V, 
                                   OutParType.k, 
                                   OutParType.n, 
                                   InpParType.d, 
                                   Y, 
                                   OutParType.c, 
                                   OutParType.W, 
                                   OutParType.Theta, 
                                   IC, 
                                   logL, 
                                   M,
                                   D);

    if (*Error) goto E0;

E0: if (OutParType.Theta) {
        for (i = 0; i < OutParType.c; i++) {
            if (OutParType.Theta[i]) free(OutParType.Theta[i]);
        }
         
        free(OutParType.Theta);
    }
    
    if (Y) {
        for (i = 0; i < OutParType.n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }

    if (InpParType.ParFamType) free(InpParType.ParFamType);
} /* RInformationCriterionH */