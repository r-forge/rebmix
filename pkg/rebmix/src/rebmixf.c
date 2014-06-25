#if (_MEMORY_LEAK_SWITCH)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include <math.h>
#include <float.h> 

#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "rebmixf.h"

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

#if (_TIME_LEFT_SWITCH)
char Progress[2048];
int  ProgressLength = 0;
#endif

/* Inserts y into ascending list Y of length n. Set n = 0 initially. */

void Insert(FLOAT y,  /* Inserted value. */
            int   *n, /* Length of Y. */
            FLOAT *Y) /* Pointer to Y = [y0,...,yn-1]. */
{
    int i, j;

    Y[*n] = y;

    for (i = 0; i < *n; i++) {
        if (y < Y[i]) {
            for (j = *n; j > i; j--) Y[j] = Y[j - 1]; Y[i] = y;
                    
            break;
        }
    }

    *n += 1;
} /* Insert */

/* Returns the value log(Gamma(y)) for y > 0. See http://www.nrbook.com/a/bookcpdf/c6-1.pdf */

FLOAT Gammaln(FLOAT y)
{
    FLOAT x, z, Tmp, Ser;
    int   j;

    static FLOAT Cof[6] = {(FLOAT)76.18009172947146, -(FLOAT)86.50532032941677,
                           (FLOAT)24.01409824083091, -(FLOAT)1.231739572450155,
                           (FLOAT)0.1208650973866179E-2, -(FLOAT)0.5395239384953E-5};

    static FLOAT Stp = (FLOAT)2.5066282746310005;

    z = x = y; Tmp = x + (FLOAT)5.5; Tmp -= (x + (FLOAT)0.5) * (FLOAT)log(Tmp);

    Ser = (FLOAT)1.000000000190015;

    for (j = 0; j < 6; j++) Ser += Cof[j] / ++z;

    return (-Tmp + (FLOAT)log(Stp * Ser / x));
} /* Gammaln */

/* Returns the digamma for y > 0. */

int Digamma(FLOAT y, FLOAT *Psi)
{
    static FLOAT piov4 = (FLOAT)0.785398163397448;
    static FLOAT dy0 = (FLOAT)1.461632144968362341262659542325721325;
    static FLOAT p1[7] = {(FLOAT)0.0089538502298197, (FLOAT)4.77762828042627, (FLOAT)142.441585084029,
                          (FLOAT)1186.45200713425, (FLOAT)3633.51846806499, (FLOAT)4138.10161269013, 
                          (FLOAT)1305.60269827897};
    static FLOAT q1[6] = {(FLOAT)44.8452573429826, (FLOAT)520.752771467162, (FLOAT)2210.0079924783,
                          (FLOAT)3641.27349079381, (FLOAT)1908.310765963, (FLOAT)6.91091682714533e-6};
    static FLOAT p2[4] = {-(FLOAT)2.12940445131011, -(FLOAT)7.01677227766759,
                          -(FLOAT)4.48616543918019, -(FLOAT)0.648157123766197};
    static FLOAT q2[4] = {(FLOAT)32.2703493791143, (FLOAT)89.2920700481861,
                          (FLOAT)54.6117738103215, (FLOAT)7.77788548522962};
    int   i, m, n, nq;
    FLOAT d2;
    FLOAT w, z;
    FLOAT den, aug, sgn, ymy0, ymax, upper, ymin;
    int Error = 0;

    ymax = (FLOAT)INT_MAX; d2 = (FLOAT)1.0 / DBL_EPSILON; if (ymax > d2) ymax = d2; ymin = (FLOAT)1E-9; aug = (FLOAT)0.0;

    if (y < (FLOAT)0.5) {
        if ((FLOAT)fabs(y) <= ymin) {
              if (y == (FLOAT)0.0) {
                Error = 1; goto E0;
            }
        
            aug = -(FLOAT)1.0 / y;
        }
        else { 
            w = -y; sgn = piov4;

            if (w <= (FLOAT)0.0) {
                w = -w; sgn = -sgn;
            }

            if (w >= ymax) {
                Error = 1; goto E0;
            }

            nq = (int)w; w -= (FLOAT)nq; nq = (int)(w * (FLOAT)4.0); w = (w - (FLOAT) nq * (FLOAT)0.25) * (FLOAT)4.0; n = nq / 2;

            if (n + n != nq) {
                w = (FLOAT)1.0 - w;
            }

            z = piov4 * w; m = n / 2;

            if (m + m != n) {
                sgn = -sgn;
            }

            n = (nq + 1) / 2; m = n / 2; m += m;
            
            if (m == n) {
                if (z == 0.0) {
                    Error = 1; goto E0;
                }

                aug = sgn * (cos(z) / sin(z) * (FLOAT)4.0);
            } 
            else { 
                aug = sgn * ((FLOAT)sin(z) / (FLOAT)cos(z) * (FLOAT)4.0);
            }
        }

        y = (FLOAT)1.0 - y;
    }

    if (y <= (FLOAT)3.0) {
        den = y; upper = p1[0] * y;

        for (i = 1; i <= 5; ++i) {
            den = (den + q1[i - 1]) * y; upper = (upper + p1[i]) * y;
        }

        den = (upper + p1[6]) / (den + q1[5]); ymy0 = y - dy0;

        *Psi = den * ymy0 + aug;
    }
    else
    if (y < ymax) {
        w = (FLOAT)1.0 / (y * y); den = w; upper = p2[0] * w;

        for (i = 1; i <= 3; ++i) {
            den = (den + q2[i - 1]) * w; upper = (upper + p2[i]) * w;
        }

        aug = upper / (den + q2[3]) - (FLOAT)0.5 / y + aug;

        *Psi = aug + (FLOAT)log(y);
    }

E0: return (Error);
} /* Digamma */

/* Returns the inverse of the binomial c.d.f. for the specified n and p. */

FLOAT BinomialInv(FLOAT Fy, FLOAT n, FLOAT p)
{
    FLOAT Sum, y, ypb;

    Sum = ypb = (FLOAT)pow((FLOAT)1.0 - p, n); y = (FLOAT)0.0;

    while ((Sum < Fy) && (ypb > FLOAT_MIN)) {
        y++; ypb *= (n - y + (FLOAT)1.0) * p / y / ((FLOAT)1.0 - p); Sum += ypb;
    }

    if ((Fy < (FLOAT)0.5) && (y > (FLOAT)0.0)) y--;

    return (y);
} /* BinomialInv */

/* Returns the inverse of the Poisson c.d.f. for the specified Theta. */

FLOAT PoissonInv(FLOAT Fy, FLOAT Theta)
{
    FLOAT Sum, y, ypb;

    Sum = ypb = (FLOAT)exp(-Theta); y = (FLOAT)0.0;

    while ((Sum < Fy) && (ypb > FLOAT_MIN)) {
        y++; ypb *= Theta / y; Sum += ypb;
    }

    if ((Fy < (FLOAT)0.5) && (y > (FLOAT)0.0)) y--;

    return (y);
} /* PoissonInv */

/* Returns the incomplete gamma function P(a, y) evaluated by its series
   representation as GamSer. Also returns log(Gamma(a)) as Gamln. */

int GammaSer(FLOAT a,       /* Constant a > 0. */
             FLOAT y,       /* Variable y > 0. */
             FLOAT *GamSer, /* Incomplete gamma function. */
             FLOAT *Gamln)  /* Log(Gamma(a)). */
{
    int   i;
    FLOAT Sum, Del, ap;
    int   Error = 0;

    *Gamln = Gammaln(a);

    if (y <= FLOAT_MIN) {
        if (y < (FLOAT)0.0) Error = 1; if (Error) goto E0;

        *GamSer = (FLOAT)0.0;
    }
    else {
        ap = a; Sum = (FLOAT)1.0 / a; Del = Sum; Error = 1; i = 1;

        while ((i <= ItMax) && Error) {
            ap += (FLOAT)1.0; Del *= y / ap; Sum += Del;

            if ((FLOAT)fabs(Del / Sum) < Eps) Error = 0;

            i++;
        }

        if (Error) goto E0;

        *GamSer = Sum * (FLOAT)exp(-y + a * log(y) - *Gamln);
    }

E0: return (Error);
} /* GammaSer */

/* Returns the incomplete gamma function Q(a, y) evaluated by its continued
   fraction representation as GamCfg. Also returns log(Gamma(a)) as Gamln. */

int GammaCfg(FLOAT a,       /* Constant a > 0. */
             FLOAT y,       /* Variable y > 0. */
             FLOAT *GamCfg, /* Incomplete gamma function. */
             FLOAT *Gamln)  /* Log(Gamma(a)). */
{
    int   i;
    FLOAT Gold, G, Fac, a0, a1, b0, b1, aif, aia, ai;
    int   Error = 0;

    *Gamln = Gammaln(a); 

    if (y <= FLOAT_MIN) {
        if (y < (FLOAT)0.0) Error = 1; if (Error) goto E0;

        *GamCfg = (FLOAT)0.0;
    }
    else {
        G = (FLOAT)0.0; Gold = (FLOAT)0.0; Fac = (FLOAT)1.0;

        a0 = (FLOAT)1.0; a1 = y; b0 = (FLOAT)0.0; b1 = (FLOAT)1.0; Error = 1; i = 1;

        while ((i <= ItMax) && Error) {
            ai = (FLOAT)1.0 * i; aia = ai - a; aif = ai * Fac;

            a0 = (a1 + a0 * aia) * Fac;
            b0 = (b1 + b0 * aia) * Fac;

            a1 = y * a0 + aif * a1;
            b1 = y * b0 + aif * b1;

            if (a1 != (FLOAT)0.0) {
                Fac = (FLOAT)1.0 / a1; G = b1 * Fac;

                if ((FLOAT)fabs((G - Gold) / G) < Eps) Error = 0; else Gold = G;
            }

            i++;
        }

        if (Error) goto E0;

        *GamCfg = (FLOAT)exp(-y + a * log(y) - *Gamln) * G;
    }

E0: return (Error);
} /* GammaCfg */

/* Returns the incomplete gamma function P(a, y). Also returns log(Gamma(a)) as Gamln. See http://www.nrbook.com/a/bookcpdf/c6-2.pdf */

int GammaP(FLOAT a,       /* Constant a > 0. */
           FLOAT y,       /* Variable y > 0. */
           FLOAT *GamP,   /* Incomplete gamma function. */
           FLOAT *Gamln)  /* Log(Gamma(a)). */
{
    FLOAT GamSer, GamCfg;
    int   Error = 0;

    if ((y < (FLOAT)0.0) || (a <= FLOAT_MIN)) {
        *GamP = (FLOAT)0.0; Error = 1; if (Error) goto E0;
    }
    else
    if (y < a + (FLOAT)1.0) {
        Error = GammaSer(a, y, &GamSer, Gamln); 

        if (Error) goto E0;
        
        *GamP = GamSer;
    }
    else {
        Error = GammaCfg(a, y, &GamCfg, Gamln); 

        if (Error) goto E0;
        
        *GamP = (FLOAT)1.0 - GamCfg;
    }

E0: return (Error);
} /* GammaP */

/* Returns the inverse of the gamma c.d.f. for the specified Theta and Beta. */

int GammaInv(FLOAT Fy, FLOAT Theta, FLOAT Beta, FLOAT *y)
{
    FLOAT dy, GamP, Gamln, Tmp;
    int   i;
    int   Error = 0;

    if (Beta > (FLOAT)1.0)
        *y = (Beta - (FLOAT)1.0) * Theta + Eps;
    else
        *y = Eps;

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        if (GammaP(Beta, *y / Theta, &GamP, &Gamln)) goto E0;

        Tmp = *y / Theta;

        dy = (GamP - Fy) / ((FLOAT)exp(Beta * (FLOAT)log(Tmp) - Tmp - Gamln) / (*y));

        *y -= dy;

        #if (_REBMIXEXE || _REBMIXR)
        if (IsNan(dy) || IsInf(dy)) {
            Error = 1; goto E0;
        }
        else
        if (*y < Eps) {
            *y = Eps; Error = 0;
        }
        #endif

        if ((FLOAT)fabs(dy / (*y)) < Eps) Error = 0; 

        i++;
    }

E0: return (Error);
} /* GammaInv */

/* Returns the error function erf(y). */

int ErrorF(FLOAT y,     /* Variable y. */
           FLOAT *ErF)  /* Error function. */
{
    FLOAT GamP, Gamln;
    int   Error = 0;

    Error = GammaP((FLOAT)0.5, y * y, &GamP, &Gamln);

    if (Error) goto E0;

    if (y < (FLOAT)0.0)
        *ErF = -GamP;
    else
        *ErF = +GamP;

E0: return (Error);
} /* ErrorF */

/* Returns component p.d.f or c.d.f. */ 

int ComponentDist(int                      d,            /* Number of independent random variables. */
                  FLOAT                    *Y,           /* Pointer to the input point [y0,...,yd-1]. */
                  MarginalDistributionType *MrgDistType, /* Marginal distribution type. */
                  FLOAT                    *CmpDist,     /* Component distribution. */
                  int                      Cumulative)   /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT y, ypb, ErF, Sum, p, Theta, GamP, Gamln;
    int   i, j, k, n;
    int   Error = 0;

    *CmpDist = (FLOAT)1.0;

    if (Cumulative) {
        for (i = 0; i < d; i++) {
            switch (MrgDistType[i].ParFamType) {
            case pfNormal:
                y = (Y[i] - MrgDistType[i].Par0) / (Sqrt2 * MrgDistType[i].Par1);

                Error = ErrorF(y, &ErF);

                if (Error) goto E0;

                *CmpDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);

                break;
            case pfLognormal:
                if (Y[i] > FLOAT_MIN) {
                    y = ((FLOAT)log(Y[i]) - MrgDistType[i].Par0) / (Sqrt2 * MrgDistType[i].Par1);

                    Error = ErrorF(y, &ErF);

                    if (Error) goto E0;

                    *CmpDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);
                }
                else {
                    *CmpDist = (FLOAT)0.0;
                }
                
                break;
            case pfWeibull:
                if (Y[i] > FLOAT_MIN) {
                    ypb = (FLOAT)exp(MrgDistType[i].Par1 * log(Y[i] / MrgDistType[i].Par0));

                    *CmpDist *= ((FLOAT)1.0 - (FLOAT)exp(-ypb));
                }
                else {
                    *CmpDist = (FLOAT)0.0;
                }

                break;
            case pfGamma:
                if (Y[i] > FLOAT_MIN) {
                    Error = GammaP(MrgDistType[i].Par1, Y[i] / MrgDistType[i].Par0, &GamP, &Gamln);

                    if (Error) goto E0;

                    *CmpDist *= GamP;
                }
                else {
                    *CmpDist = (FLOAT)0.0;
                }

                break;
            case pfBinomial:
                k = (int)Y[i]; n = (int)MrgDistType[i].Par0; p = MrgDistType[i].Par1;

                if (k < 0)
                    *CmpDist *= (FLOAT)0.0;
                else
                if (k == 0)
                    *CmpDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
                else
                if (k >= n)
                    *CmpDist *= (FLOAT)1.0;
                else {
                    Sum = ypb = (FLOAT)pow((FLOAT)1.0 - p, n);
    
                    for (j = 1; j <= k; j++) {
                        ypb *= (n - j + (FLOAT)1.0) * p / j / ((FLOAT)1.0 - p); Sum += ypb; 
                    }

                    *CmpDist *= Sum;
                }

                break;
            case pfPoisson:
                k = (int)Y[i]; Theta = MrgDistType[i].Par0;

                Sum = ypb = (FLOAT)exp(-Theta);
    
                for (j = 1; j <= k; j++) {
                   ypb *= Theta / j; Sum += ypb;
                }

                *CmpDist *= Sum;

                break;
            case pfDirac:
                if (Y[i] < MrgDistType[i].Par0) {
                    *CmpDist *= (FLOAT)0.0;
                }
            }
        }
    }
    else {
        for (i = 0; i < d; i++) {
            switch (MrgDistType[i].ParFamType) {
            case pfNormal:
                y = (Y[i] - MrgDistType[i].Par0) / (Sqrt2 * MrgDistType[i].Par1);

                *CmpDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Par1);

                break;
            case pfLognormal:
                if (Y[i] > FLOAT_MIN) {
                    y = ((FLOAT)log(Y[i]) - MrgDistType[i].Par0)/(Sqrt2 * MrgDistType[i].Par1);

                    *CmpDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Par1) / Y[i];
                }
                else {
                    *CmpDist = (FLOAT)0.0;
                }

                break;
            case pfWeibull:
                if (Y[i] > FLOAT_MIN) {
                    ypb = (FLOAT)exp(MrgDistType[i].Par1 * log(Y[i] /  MrgDistType[i].Par0));

                    *CmpDist *= MrgDistType[i].Par1 * ypb * (FLOAT)exp(-ypb) / Y[i];
                }
                else {
                    *CmpDist = (FLOAT)0.0;
                }

                break;
            case pfGamma:
                if (Y[i] > FLOAT_MIN) {
                    ypb = Y[i] /  MrgDistType[i].Par0;

                    *CmpDist *= (FLOAT)exp(MrgDistType[i].Par1 * log(ypb) - ypb - Gammaln(MrgDistType[i].Par1)) / Y[i];
                }
                else {
                    *CmpDist = (FLOAT)0.0;
                }

                break;
            case pfBinomial:
                k = (int)Y[i]; n = (int)MrgDistType[i].Par0; p = MrgDistType[i].Par1;

                if (k < 0)
                    *CmpDist *= (FLOAT)0.0;
                else
                if (k == 0)
                    *CmpDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
                else
                if (k == n)
                    *CmpDist *= (FLOAT)pow(p, n);
                else
                if (k > n)
                    *CmpDist *= (FLOAT)0.0;
                else
                    *CmpDist *= (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(k + (FLOAT)1.0) - Gammaln(n - k + (FLOAT)1.0)) *
                                (FLOAT)pow(p, k) * (FLOAT)pow((FLOAT)1.0 - p, n - k);

                break;
            case pfPoisson:
                k = (int)Y[i]; Theta = MrgDistType[i].Par0;

                   *CmpDist *= (FLOAT)exp(k * log(Theta) - Theta - Gammaln(k + (FLOAT)1.0));

                break;
            case pfDirac:
                if ((FLOAT)fabs(Y[i] - MrgDistType[i].Par0) > FLOAT_MIN) {
                    *CmpDist *= (FLOAT)0.0;
                }
            }
        }
    }

E0: return (Error);
} /* ComponentDist */

/* Returns component marginal p.d.f or c.d.f. */ 

int ComponentMarginalDist(int                      i,            /* Index of variable y. */
                          FLOAT                    *Y,           /* Pointer to the input point [y0,...,yd-1]. */
                          MarginalDistributionType *MrgDistType, /* Marginal distribution type. */
                          FLOAT                    *CmpMrgDist,  /* Component marginal distribution. */
                          int                      Cumulative)   /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT y, ypb, Sum, p, Theta, GamP, Gamln;
    int   j, k, n;
    FLOAT ErF;
    int   Error = 0;
    
    *CmpMrgDist = (FLOAT)1.0;

    if (Cumulative) {
        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            y = (Y[i] - MrgDistType[i].Par0) / (Sqrt2 * MrgDistType[i].Par1);

            Error = ErrorF(y, &ErF);

            if (Error) goto E0;

            *CmpMrgDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);

            break;
        case pfLognormal:
            if (Y[i] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i]) - MrgDistType[i].Par0) / (Sqrt2 * MrgDistType[i].Par1);

                Error = ErrorF(y, &ErF);

                if (Error) goto E0;

                *CmpMrgDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);
            }
            else {
                *CmpMrgDist = (FLOAT)0.0;
            }
                
            break;
        case pfWeibull:
            if (Y[i] > FLOAT_MIN) {
                ypb = (FLOAT)exp(MrgDistType[i].Par1 * log(Y[i] / MrgDistType[i].Par0));

                *CmpMrgDist *= ((FLOAT)1.0 - (FLOAT)exp(-ypb));
            }
            else {
                *CmpMrgDist = (FLOAT)0.0;
            }

            break;
        case pfGamma:
            if (Y[i] > FLOAT_MIN) {
                Error = GammaP(MrgDistType[i].Par1, Y[i] / MrgDistType[i].Par0, &GamP, &Gamln);

                if (Error) goto E0;

                *CmpMrgDist *= GamP;
            }
            else {
                *CmpMrgDist = (FLOAT)0.0;
            }

            break;
        case pfBinomial:
            k = (int)Y[i]; n = (int)MrgDistType[i].Par0; p = MrgDistType[i].Par1;

            if (k < 0)
                *CmpMrgDist *= (FLOAT)0.0;
            else
            if (k == 0)
                *CmpMrgDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
            else
            if (k >= n)
                *CmpMrgDist *= (FLOAT)1.0;
            else {
                Sum = ypb = (FLOAT)pow((FLOAT)1.0 - p, n);
    
                for (j = 1; j <= k; j++) {
                    ypb *= (n - j + (FLOAT)1.0) * p / j / ((FLOAT)1.0 - p); Sum += ypb; 
                }

                *CmpMrgDist *= Sum;
            }

            break;
        case pfPoisson:
            k = (int)Y[i]; Theta = MrgDistType[i].Par0;

            Sum = ypb = (FLOAT)exp(-Theta);
    
            for (j = 1; j <= k; j++) {
                ypb *= Theta / j; Sum += ypb; 
            }

            *CmpMrgDist *= Sum;

            break;
        case pfDirac:
            if (Y[i] < MrgDistType[i].Par0) {
                *CmpMrgDist *= (FLOAT)0.0;
            }
        }
    }
    else {
        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            y = (Y[i] - MrgDistType[i].Par0) / (Sqrt2 * MrgDistType[i].Par1);

            *CmpMrgDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Par1);

            break;
        case pfLognormal:
            if (Y[i] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i]) - MrgDistType[i].Par0) / (Sqrt2 * MrgDistType[i].Par1);

                *CmpMrgDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Par1) / Y[i];
            }
            else {
                *CmpMrgDist = (FLOAT)0.0;
            }

            break;
        case pfWeibull:
            if (Y[i] > FLOAT_MIN) {
                ypb = (FLOAT)exp(MrgDistType[i].Par1 * log(Y[i] /  MrgDistType[i].Par0));

                *CmpMrgDist *= MrgDistType[i].Par1 * ypb * (FLOAT)exp(-ypb) / Y[i];
            }
            else {
                *CmpMrgDist = (FLOAT)0.0;
            }

            break;
        case pfGamma:
            if (Y[i] > FLOAT_MIN) {
                ypb = Y[i] /  MrgDistType[i].Par0;

                *CmpMrgDist *= (FLOAT)exp(MrgDistType[i].Par1 * log(ypb) - ypb - Gammaln(MrgDistType[i].Par1)) / Y[i];
            }
            else {
                *CmpMrgDist = (FLOAT)0.0;
            }

            break;
        case pfBinomial:
            k = (int)Y[i]; n = (int)MrgDistType[i].Par0; p = MrgDistType[i].Par1;

            if (k < 0)
                *CmpMrgDist *= (FLOAT)0.0;
            else
            if (k == 0)
                *CmpMrgDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
            else
            if (k == n)
                *CmpMrgDist *= (FLOAT)pow(p, n);
            else
            if (k > n)
                *CmpMrgDist *= (FLOAT)0.0;
            else
                *CmpMrgDist *= (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(k + (FLOAT)1.0) - Gammaln(n - k + (FLOAT)1.0)) *
                               (FLOAT)pow(p, k) * (FLOAT)pow((FLOAT)1.0 - p, n - k);

            break;
        case pfPoisson:
            k = (int)Y[i]; Theta = MrgDistType[i].Par0;

               *CmpMrgDist *= (FLOAT)exp(k * log(Theta) - Theta - Gammaln(k + (FLOAT)1.0));

            break;
        case pfDirac:
            if ((FLOAT)fabs(Y[i] - MrgDistType[i].Par0) > FLOAT_MIN) {
                *CmpMrgDist *= (FLOAT)0.0;
            }
        }
    }

E0: return (Error);
} /* ComponentMarginalDist */

/* Returns mixture p.d.f or c.d.f. */ 

int MixtureDist(int                      d,             /* Number of independent random variables. */
                FLOAT                    *Y,            /* Pointer to the input point [y0,...,yd-1]. */
                int                      c,             /* Number of components. */ 
                FLOAT                    *W,            /* Component weights. */
                MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                FLOAT                    *MixDist,      /* Mixture distribution. */
                int                      Cumulative)    /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT CmpDist;
    int   i;
    int   Error = 0;

    *MixDist = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = ComponentDist(d, Y, MrgDistType[i], &CmpDist, Cumulative);

        if (Error) goto E0;

        *MixDist += W[i] * CmpDist;
    }

E0: return (Error);
} /* MixtureDist */

/* Returns mixture marginal p.d.f or c.d.f. */ 

int MixtureMarginalDist(int                      i,             /* Index of variable y. */  
                        FLOAT                    *Y,            /* Pointer to the input point [y0,...,yd-1]. */
                        int                      c,             /* Number of components. */ 
                        FLOAT                    *W,            /* Component weights. */
                        MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                        FLOAT                    *MixMrgDist,   /* Mixture marginal distribution. */
                        int                      Cumulative)    /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT CmpMrgDist;
    int   j;
    int   Error = 0;

    *MixMrgDist = (FLOAT)0.0;

    for (j = 0; j < c; j++) {
        Error = ComponentMarginalDist(i, Y, MrgDistType[j], &CmpMrgDist, Cumulative);

        if (Error) goto E0;

        *MixMrgDist += W[j] * CmpMrgDist;
    }

E0: return (Error);
} /* MixtureMarginalDist */

int DegreesOffreedom(int                      d,             /* Number of independent random variables. */
                     int                      c,             /* Number of components. */ 
                     MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                     int                      *M)            /* Degrees of freedom. */
{
    int i, j;
    int Error = 0;

    *M = c - 1;

    for (i = 0; i < c; i++) {
        for (j = 0; j < d; j++) {
            switch (MrgDistType[i][j].ParFamType) {
            case pfNormal:
                *M += 2;

                break;
            case pfLognormal:
                *M += 2;
                
                break;
            case pfWeibull:
                *M += 2;

                break;
            case pfGamma:
                *M += 2;

                break;
            case pfBinomial:
                *M += 1;

                break;
            case pfPoisson:
                *M += 1;

                break;
            case pfDirac:
                *M += 1;
            }
        }
    }

    return (Error);
} /* DegreesOffreedom */

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
                            FLOAT                      *D)            /* Total of positive relative deviations. */
{
    int   i, j;
    FLOAT E, SSE, EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    Error = DegreesOffreedom(d, c, MrgDistType, M);

    if (Error) goto E0;

    *IC = *logL = EN = *D = SSE = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < n; i++) {
        Error = MixtureDist(d, Y[i], c, W, MrgDistType, &MixDist, 0);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *logL += (FLOAT)log(MixDist);
        }

        E = Y[i][d] / n - MixDist * Y[i][d + 1] / k;
    
        if (E > (FLOAT)0.0) {
            *D += E;
        }

        switch (ICType) {
        case icAWE: case icCLC: case icICL: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                tau = W[j] * CmpDist / MixDist; 
                
                if (tau > FLOAT_MIN) {
                    EN -= tau * (FLOAT)log(tau); PC += tau * tau; 
                }
            }

            break;
        case icSSE:
            E = Y[i][d] - n * MixDist * Y[i][d + 1] / k; E *= E;
    
             SSE += E;

            break;
        default:
            EN = SSE = PC = (FLOAT)0.0;
        }
    }

    switch (ICType) {
    case icAIC: /* AIC - Akaike information criterion Akaike (1973). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M);

        break;
    case icAIC3: /* AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * (*M);

        break;
    case icAIC4: /* AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * (*M);

        break;
    case icAICc: /* AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * ((FLOAT)1.0 + ((*M) + 1) / (n - (*M) - (FLOAT)1.0));

        break;
    case icBIC: /* BIC - Bayesian information criterion Schwarz (1978). */
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icCAIC: /* CAIC - Consistent Akaike information criterion Bozdogan (1987). */
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * ((FLOAT)log((FLOAT)n) + (FLOAT)1.0);

        break;
    case icHQC: /* HQC - Hannan-Quinn information criterion Hannan & Quinn (1979). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log(log((FLOAT)n));

        break;
    case icMDL2: /* MDL2 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icMDL5: /* MDL5 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icAWE: /* AWE - Approximate weight of evidence criterion Banfield & Raftery (1993). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * (*M) * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n));

        break;
    case icCLC: /* CLC - Classification likelihood criterion Biernacki & Govaert (1997). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: /* ICL - Integrated classification likelihood Biernacki et al. (1998). */
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n * PW - (FLOAT)2.0 * K + ((*M) - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n);

        break;
    case icPC: /* PC - Partition coeficient Bezdek (1981). */
        *IC = PC; 

        break;
    case icICLBIC: /* ICL-BIC - Integrated classification likelihood criterion Biernacki et al. (1998). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icD: /* D - Total of positive relative deviations Nagode & Fajdiga (2011). */
        *IC = *D;

        break;
    case icSSE: /* SSE - Sum of squares error Bishop (1998). */
        *IC = (FLOAT)0.5 * SSE;
    }

E0: return (Error);
} /* InformationCriterionKNN */

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
                           FLOAT                      *D)            /* Total of positive relative deviations. */
{
    int   i, j;
    FLOAT E, SSE, EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    Error = DegreesOffreedom(d, c, MrgDistType, M);

    if (Error) goto E0;

    *IC = *logL = EN = *D = SSE = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < n; i++) {
        Error = MixtureDist(d, Y[i], c, W, MrgDistType, &MixDist, 0);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *logL += (FLOAT)log(MixDist);
        }

        E =  Y[i][d] / n - MixDist * V / Y[i][d + 1];
    
        if (E > (FLOAT)0.0) {
            *D += E;
        }

        switch (ICType) {
        case icAWE: case icCLC: case icICL: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                tau = W[j] * CmpDist / MixDist; 
                
                if (tau > FLOAT_MIN) {
                    EN -= tau * (FLOAT)log(tau); PC += tau * tau; 
                }
            }

            break;
        case icSSE:
            E = Y[i][d] - n * MixDist * V / Y[i][d + 1]; E *= E;
    
             SSE += E;

            break;
        default:
            EN = SSE = PC = (FLOAT)0.0;
        }
    }

    switch (ICType) {
    case icAIC: /* AIC - Akaike information criterion Akaike (1973). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M);

        break;
    case icAIC3: /* AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * (*M);

        break;
    case icAIC4: /* AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * (*M);

        break;
    case icAICc: /* AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * ((FLOAT)1.0 + ((*M) + 1) / (n - (*M) - (FLOAT)1.0));

        break;
    case icBIC: /* BIC - Bayesian information criterion Schwarz (1978). */
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icCAIC: /* CAIC - Consistent Akaike information criterion Bozdogan (1987). */
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * ((FLOAT)log((FLOAT)n) + (FLOAT)1.0);

        break;
    case icHQC: /* HQC - Hannan-Quinn information criterion Hannan & Quinn (1979). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log(log((FLOAT)n));

        break;
    case icMDL2: /* MDL2 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icMDL5: /* MDL5 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icAWE: /* AWE - Approximate weight of evidence criterion Banfield & Raftery (1993). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * (*M) * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n));

        break;
    case icCLC: /* CLC - Classification likelihood criterion Biernacki & Govaert (1997). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: /* ICL - Integrated classification likelihood Biernacki et al. (1998). */
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n * PW - (FLOAT)2.0 * K + ((*M) - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n);

        break;
    case icPC: /* PC - Partition coeficient Bezdek (1981). */
        *IC = PC; 

        break;
    case icICLBIC: /* ICL-BIC - Integrated classification likelihood criterion Biernacki et al. (1998). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icD: /* D - Total of positive relative deviations Nagode & Fajdiga (2011). */
        *IC = *D;

        break;
    case icSSE: /* SSE - Sum of squares error Bishop (1998). */
        *IC = (FLOAT)0.5 * SSE;
    }

E0: return (Error);
} /* InformationCriterionPW */

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
                          FLOAT                      *D)            /* Total of positive relative deviations. */
{
    int   i, j;
    FLOAT E, SSE, EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    Error = DegreesOffreedom(d, c, MrgDistType, M);

    if (Error) goto E0;

    *IC = *logL = EN = *D = SSE = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < k; i++) {
        Error = MixtureDist(d, Y[i], c, W, MrgDistType, &MixDist, 0);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *logL += Y[i][d] * (FLOAT)log(MixDist);
        }

        E = Y[i][d] / n - MixDist * V;
    
        if (E > (FLOAT)0.0) {
            *D += E;
        }

        switch (ICType) {
        case icAWE: case icCLC: case icICL: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                tau = W[j] * CmpDist / MixDist; 

                if (tau > FLOAT_MIN) {
                    EN -= Y[i][d] * tau * (FLOAT)log(tau); PC += Y[i][d] * tau * tau; 
                }
            }

            break;
        case icSSE:
            E = Y[i][d] - n * MixDist * V; E *= E;
    
            SSE += E;

            break;
        default:
            EN = SSE = PC = (FLOAT)0.0;
        }
    }

    switch (ICType) {
    case icAIC: /* AIC - Akaike information criterion Akaike (1973). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M);

        break;
    case icAIC3: /* AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * (*M);

        break;
    case icAIC4: /* AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * (*M);

        break;
    case icAICc: /* AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * ((FLOAT)1.0 + ((*M) + 1) / (n - (*M) - (FLOAT)1.0));

        break;
    case icBIC: /* BIC - Bayesian information criterion Schwarz (1978). */
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icCAIC: /* CAIC - Consistent Akaike information criterion Bozdogan (1987). */
        *IC = -(FLOAT)2.0 * (*logL) + (*M) * ((FLOAT)log((FLOAT)n) + (FLOAT)1.0);

        break;
    case icHQC: /* HQC - Hannan-Quinn information criterion Hannan & Quinn (1979). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log(log((FLOAT)n));

        break;
    case icMDL2: /* MDL2 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icMDL5: /* MDL5 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icAWE: /* AWE - Approximate weight of evidence criterion Banfield & Raftery (1993). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * (*M) * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n));

        break;
    case icCLC: /* CLC - Classification likelihood criterion Biernacki & Govaert (1997). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: /* ICL - Integrated classification likelihood Biernacki et al. (1998). */
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n * PW - (FLOAT)2.0 * K + ((*M) - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n);

        break;
    case icPC: /* PC - Partition coeficient Bezdek (1981). */
        *IC = PC; 

        break;
    case icICLBIC: /* ICL-BIC - Integrated classification likelihood criterion Biernacki et al. (1998). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (*M) * (FLOAT)log((FLOAT)n);

        break;
    case icD: /* D - Total of positive relative deviations Nagode & Fajdiga (2011). */
        *IC = *D;

        break;
    case icSSE: /* SSE - Sum of squares error Bishop (1998). */
        *IC = (FLOAT)0.5 * SSE;
    }

E0: return (Error);
} /* InformationCriterionH */

/* Preprocessing of observations for k-nearest neighbour. */

int PreprocessingKNN(int   k,    /* k-nearest neighbours. */
                     FLOAT *h,   /* Normalizing vector. */
                     int   n,    /* Total number of independent observations. */
                     int   d,    /* Number of independent random variables. */ 
                     FLOAT **Y)  /* Pointer to the input array [y0,...,yd-1,kl,V,R]. */
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, V, Vn;
    int   i, j, l, m, q;
    int   Error = n < 1;

    if (Error) goto E0;

    if (k > 1) k -= 1; else k = 1; 
    
    Dk = (FLOAT*)malloc(k * sizeof(FLOAT));

    Error = NULL == Dk; if (Error) goto E0;

    Vn = (FLOAT)exp(d * LogPi / (FLOAT)2.0 - Gammaln((FLOAT)1.0 + d / (FLOAT)2.0));

    for (i = 0; i < n; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < n; j++) if (i != j) {
            Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                R = (Y[i][l] - Y[j][l]) / h[l]; Dc += R * R;
            }

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < k; l++) {
                if (Dc < Dk[l]) {
                    for (m = k - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != k - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = (FLOAT)sqrt(Dk[k - 1]); 

        if (q >= k) R *= (FLOAT)exp(log((k + (FLOAT)1.0) / (q + (FLOAT)2.0)) / d);

        V = Vn * (FLOAT)exp(d * log(R));

        for (j = 0; j < d; j++) V *= h[j];

        Y[i][d] = (FLOAT)1.0; Y[i][d + 1] = V; Y[i][d + 2] = R;
    }

E0: if (Dk) free(Dk);

    return (Error);
} /* PreprocessingKNN */

/* Preprocessing of observations for Parzen window. */

int PreprocessingPW(FLOAT *h,   /* Sides of the hypersquare. */
                    int   n,    /* Total number of independent observations. */
                    int   d,    /* Number of independent random variables. */ 
                    FLOAT **Y)  /* Pointer to the input array [y0,...,yd-1,kl,k]. */
{
    int i, j, k;
    int Error = n < 1;

    if (Error) goto E0;

    for (i = 0; i < n; i++) {
        Y[i][d] = (FLOAT)1.0; Y[i][d + 1] = (FLOAT)0.0;
    }

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            for (k = 0; k < d; k++) if ((FLOAT)fabs(Y[i][k] - Y[j][k]) > (FLOAT)0.5 * h[k]) goto S0;
            
            Y[i][d + 1] += (FLOAT)1.0; if (i != j) Y[j][d + 1] += (FLOAT)1.0;
S0:;    }
    }

E0: return (Error);
} /* PreprocessingPW */

/* Preprocessing of observations for histogram. */

int PreprocessingH(FLOAT           *h,          /* Sides of the hypersquare. */
                   FLOAT           *y0,         /* Origin. */
                   VariablesType_e *VarType,    /* Types of variables. */
                   int             *k,          /* Total number of bins. */
                   int             n,           /* Total number of independent observations. */
                   int             d,           /* Number of independent random variables. */ 
                   FLOAT           **X,         /* Pointer to the input points [x0,...,xd-1]. */
                   FLOAT           **Y)         /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int i, j, l, m = *k - 1;
    int Error = n < 1;

    if (Error) goto E0;

    *k = 0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            l = (int)floor((X[i][j] - y0[j]) / h[j] + (FLOAT)0.5);

            switch (VarType[j]) {
            case vtContinuous:
                if (l < 0) l = 0; else if (l > m) l = m;

                break;
            case vtDiscrete:
                break;
            }

            Y[*k][j] = y0[j] + l * h[j];
        }

        for (j = 0; j < *k; j++) {
            for (l = 0; l < d; l++) if ((FLOAT)fabs(Y[j][l] - Y[*k][l]) > (FLOAT)0.5 * h[l]) goto S0;

            Y[j][d] += (FLOAT)1.0; goto S1;
S0:;    }

        Y[*k][d] = (FLOAT)1.0; (*k)++;
S1:;}

E0: return (Error);
} /* PreprocessingH */

/* Global mode detection for k-nearest neighbour. */

int GlobalModeKNN(int   *m,   /* Global mode. */
                  int   n,    /* Total number of independent observations. */
                  int   d,    /* Number of independent random variables. */ 
                  FLOAT **Y)  /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int i, j;
    int Error = 0;

    j = 0; 

    for (i = 1; i < n; i++) if (Y[i][d] > FLOAT_MIN) {
        if (Y[i][d] / Y[i][d + 1] > Y[j][d] / Y[j][d + 1]) {
            j = i;
        }
    }

    *m = j;

    return (Error);
} /* GlobalModeKNN */

/* Global mode detection for Parzen window. */

int GlobalModePW(int   *m,    /* Global mode. */
                 int   n,     /* Total number of independent observations. */
                 int   d,     /* Number of independent random variables. */ 
                 FLOAT **Y)   /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int i, j;
    int Error = 0;

    j = 0; 

    for (i = 1; i < n; i++) if (Y[i][d] > FLOAT_MIN) {
        if (Y[i][d] * Y[i][d + 1] > Y[j][d] * Y[j][d + 1]) {
            j = i;
        }
    }

    *m = j;

    return (Error);
} /* GlobalModePW */

/* Global mode detection for histogram. */

int GlobalModeH(int   *m,    /* Global mode. */
                int   k,     /* Total number of bins. */
                int   d,     /* Number of independent random variables. */ 
                FLOAT **Y)   /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int i, j;
    int Error = 0;

    j = 0; 

    for (i = 1; i < k; i++) if (Y[i][d] > FLOAT_MIN) {
        if (Y[i][d] > Y[j][d]) {
            j = i;
        }
    }

    *m = j;

    return (Error);
} /* GlobalModeH */

/* Returns rough normal parameters. */

int RoughNormalParameters(FLOAT ym,   
                          FLOAT fm,
                          FLOAT *Mean,
                          FLOAT *Stdev)
{
    int Error = 0;

    *Mean = ym; *Stdev = (FLOAT)1.0 / (Sqrt2Pi * fm);

    return (Error);
} /* RoughNormalParameters */

/* Returns rough lognormal parameters. */

int RoughLognormalParameters(FLOAT ym,   
                             FLOAT fm,
                             FLOAT *Mean,
                             FLOAT *Stdev)
{
    FLOAT Lambda, dLambda;
    FLOAT A[3];
    int   i;
    int   Error = 0;

    Error = ym <= FLOAT_MIN; if (Error) goto E0;

    Lambda = (FLOAT)1.0 + Eps;

    i = 1; A[0] = (FLOAT)2.0 * (FLOAT)log(Sqrt2Pi * ym * fm); Error = 1;
    while ((i <= ItMax) && Error) {
        A[1] = (FLOAT)1.0 / Lambda; A[2] = Lambda - (FLOAT)1.0;

        dLambda = ((FLOAT)1.0 - A[1] + (FLOAT)log(Lambda * A[2]) + A[0]) / (A[1] * ((FLOAT)1.0 + A[1]) + (FLOAT)1.0 / A[2]);

        Lambda -= dLambda;

        #if (_REBMIXEXE || _REBMIXR)
        if (IsNan(dLambda) || IsInf(dLambda)) {
            Error = 1; goto E0;
        }
        else
        if (Lambda < (FLOAT)1.0 + Eps) {
            Lambda = (FLOAT)1.0 + Eps; Error = 0;
        }
        #endif

        if ((FLOAT)fabs(dLambda / Lambda) < Eps) Error = 0;

        i++;
    }

    if (Error) goto E0;

    *Mean = Lambda - (FLOAT)1.0 + (FLOAT)log(ym);

    *Stdev = (FLOAT)pow(Lambda * (Lambda - (FLOAT)1.0), (FLOAT)0.5);

E0: return (Error);
} /* RoughLognormalParameters */

/* Returns rough Weibull parameters. */

int RoughWeibullParameters(FLOAT ym,   
                           FLOAT fm,
                           FLOAT *Theta,
                           FLOAT *Beta)
{
    FLOAT Alpha, dAlpha;
    FLOAT A[4];
    int   i;
    int   Error = 0;

    Error = ym <= FLOAT_MIN; if (Error) goto E0;
                     
    Alpha = (FLOAT)1.3349695;

    i = 1; A[0] = (FLOAT)exp((FLOAT)1.0) * ym * fm; Error = 1;
    while ((i <= ItMax) && Error) {
        A[1] = Alpha - (FLOAT)1.0; 

        A[2] = (FLOAT)1.0 + (Euler + (FLOAT)log(A[1] / Alpha)) / Alpha;

        A[3] = (FLOAT)exp((FLOAT)1.0 / Alpha); 
        
        dAlpha = (A[2] * A[1] * A[3] - A[0]) / (A[3] * ((FLOAT)1.0 - (A[1] - A[2]) / Alpha / Alpha));

        Alpha -= dAlpha;

        #if (_REBMIXEXE || _REBMIXR)                                
        if (IsNan(dAlpha) || IsInf(dAlpha)) {
            Error = 1; goto E0;
        }
        else
        if (Alpha < (FLOAT)1.234332) {
            Alpha = (FLOAT)1.234332; Error = 0;
        }
        #endif

        if ((FLOAT)fabs(dAlpha / Alpha) < Eps) Error = 0; 

        i++;
    }

    if (Error) goto E0;

    *Beta = Alpha + Euler + (FLOAT)log((FLOAT)1.0 - (FLOAT)1.0 / Alpha);

    *Theta = ym * (FLOAT)pow(Alpha / (Alpha - (FLOAT)1.0), (FLOAT)1.0 / (*Beta));

E0: return (Error);
} /* RoughWeibullParameters */

/* Returns rough Gamma parameters. */

int RoughGammaParameters(FLOAT ym,   
                         FLOAT fm,
                         FLOAT *Theta,
                         FLOAT *Beta)
{

    FLOAT Alpha, dAlpha;
    FLOAT A[5];
    int   i;
    int   Error = 0;

    Error = ym <= FLOAT_MIN; if (Error) goto E0;

    Alpha = (FLOAT)1.00032;

    i = 1; A[0] = (FLOAT)log(ym * fm * Sqrt2Pi); Error = 1;
    while ((i <= ItMax) && Error) {
        A[1] = (FLOAT)log((FLOAT)1.0 - (FLOAT)1.0 / Alpha);

        A[2] = A[1] + (FLOAT)1.0 / Alpha; 

        A[3] = Euler * ((FLOAT)1.0 + Alpha) / (Euler - (FLOAT)1.0 - Alpha * A[1]);

        A[4] = A[3] * ((FLOAT)1.0 + A[3] * (A[1] + (FLOAT)1.0 / (Alpha - (FLOAT)1.0)) / Euler) / ((FLOAT)1.0 + Alpha);

        dAlpha = (A[3] * A[2] + (FLOAT)0.5 * (FLOAT)log(A[3]) - A[0]) / (A[4] * (A[2] + (FLOAT)0.5 / A[3]) + A[3] / (Alpha - (FLOAT)1.0) / Alpha / Alpha);

        Alpha -= dAlpha;

        #if (_REBMIXEXE || _REBMIXR)
        if (IsNan(dAlpha) || IsInf(dAlpha)) {
            Error = 1; goto E0;
        }
        else
        if (Alpha < 1.00032) {
            Alpha = 1.00032; Error = 0;
        }
        #endif

        if ((FLOAT)fabs(dAlpha / (Alpha)) < Eps) Error = 0; 

        i++;
    }

    if (Error) goto E0;

    *Beta = Euler * ((FLOAT)1.0 + Alpha) / (Euler - (FLOAT)1.0 - Alpha * (FLOAT)log((FLOAT)1.0 - (FLOAT)1.0 / Alpha));

    *Theta = ym * Alpha / (Alpha - (FLOAT)1.0) / (*Beta);

E0: return (Error);
} /* RoughGammaParameters */

/* Returns rough binomial parameters. */

int RoughBinomialParameters(FLOAT ym,
                            FLOAT fm,
                            FLOAT n,
                            FLOAT *p)
{
    int Error = 0;

    if ((int)ym == 0) {
        *p = (fm < (FLOAT)1.0) ? (FLOAT)1.0 - (FLOAT)pow(fm, (FLOAT)1.0 / n) : (FLOAT)0.0;
    }
    else
    if ((int)ym == (int)n) {
        *p = (fm < (FLOAT)1.0) ? (FLOAT)pow(fm, (FLOAT)1.0 / n) : (FLOAT)1.0;
    }
    else {
        *p = ym / n;
    }

    return (Error);
} /* RoughBinomialParameters */

/* Returns rough Poisson parameters. */

int RoughPoissonParameters(FLOAT ym,
                           FLOAT fm,
                           FLOAT *Theta)
{
    int Error = 0;

    if ((int)ym == 0) {
        *Theta = (fm < (FLOAT)1.0) ? -(FLOAT)log(fm) : (FLOAT)0.0;
    }
    else {
        *Theta = ym;
    }

    return (Error);
} /* RoughPoissonParameters */

/* Rough component parameter estimation for k-nearest neighbours. */

int RoughEstimationKNN(int                      n,             /* Total number of independent observations. */
                       int                      d,             /* Number of independent random variables. */ 
                       FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,V,R]. */
                       int                      k,             /* k-nearest neighbours. */
                       FLOAT                    *h,            /* Normalizing vector. */
                       FLOAT                    nl,            /* Total number of observations in class l. */
                       int                      m,             /* Mode index. */
                       MarginalDistributionType *MrgDistType,  /* Marginal distribution type. */
                       PestraintsType_e         ResType)       /* Restraints type. */  
{
    int                i, j, l, o, p;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, flm, flmin, flmax, Dlm, Dlmin, Dc, R;
    int                Error = 0, Stop;

    Mode = (RoughParameterType*)malloc(d * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    /* Rigid restraints. */

    flm = (FLOAT)1.0;

    for (i = 0; i < d; i++) {
        if (d > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < n; j++) {
                Dc = (FLOAT)0.0;

                for (l = 0; l < d; l++) if (i != l) {
                    R = (Y[j][l] - Y[m][l]) / h[l]; Dc += R * R;
                }

                R = (FLOAT)sqrt(Dc);

                if (R > Y[m][d + 2]) goto S0;

                Mode[i].klm += Y[j][d];
S0:;        }
        }
        else
            Mode[i].klm = nl;

        Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][d] * k / (Mode[i].klm * (FLOAT)2.0 * Y[m][d + 2] * h[i]); flm *= Mode[i].flm;
    }

    epsilon = (FLOAT)exp(log(Y[m][d] * k / (nl * Y[m][d + 1] * flm)) / d);

    for (i = 0; i < d; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].flm *= epsilon;

        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfGamma:
            Error = RoughGammaParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].flm, MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfPoisson:
            Error = RoughPoissonParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0);

            if (Error) goto E0;

            break;
        case pfDirac:
            MrgDistType[i].Par0 = Mode[i].ym;
        }
    }

    if (ResType == rtRigid) goto E0;

    /* Loose restraints. */

    for (i = 0; i < d; i++) {
        if (MrgDistType[i].ParFamType == pfDirac) goto E1;

        /* Bracketing. */

        Dlm = (FLOAT)0.0; 

        for (o = 0; o < n; o++) if (Y[o][d] > FLOAT_MIN) {
            Dc = (FLOAT)0.0;

            for (p = 0; p < d; p++) if (i != p) {
                R = (Y[o][p] - Y[m][p]) / h[p]; Dc += R * R;
            }

            R = (FLOAT)sqrt(Dc);

            if (R > Y[m][d + 2]) goto S1;

            Error = ComponentMarginalDist(i, Y[o], MrgDistType, &CmpMrgDist, 0);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * (FLOAT)2.0 * Y[o][d + 2] * h[i] / k;
S1:;    } 

        Dlm += (FLOAT)0.998;

        if (Dlm > (FLOAT)0.0) goto E1; 

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)0.998;
        flmax = (FLOAT)Mode[i].flm;

        /* Bisection. */

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            switch (MrgDistType[i].ParFamType) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;
                
                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfGamma:
                Error = RoughGammaParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].ym, flm, MrgDistType[i].Par0, &MrgDistType[i].Par1);
                                                       
                if (Error) goto E0;

                break;
            case pfPoisson:
                Error = RoughPoissonParameters(Mode[i].ym, flm, &MrgDistType[i].Par0);
                                                       
                if (Error) goto E0;

                break;
            case pfDirac:
                break;
            }

            Dlm = (FLOAT)0.0; 

            for (o = 0; o < n; o++) if (Y[o][d] > FLOAT_MIN) {
                Dc = (FLOAT)0.0;

                for (p = 0; p < d; p++) if (i != p) {
                    R = (Y[o][p] - Y[m][p]) / h[p]; Dc += R * R;
                }

                R = (FLOAT)sqrt(Dc);

                if (R > Y[m][d + 2]) goto S2;

                Error = ComponentMarginalDist(i, Y[o], MrgDistType, &CmpMrgDist, 0);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * (FLOAT)2.0 * Y[o][d + 2] * h[i] / k;
S2:;        } 

            Dlm += (FLOAT)0.998;

            if (((FLOAT)fabs(Dlm) < Eps) || ((flmax - flmin) / (FLOAT)2.0 < Eps)) {
                Stop = 1;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm; 
                }
                else {
                    flmax = flm;
                }
            }
        }
E1:;}

E0: if (Mode) free(Mode);

    return (Error);
} /* RoughEstimationKNN */

/* Rough component parameter estimation for Parzen window. */

int RoughEstimationPW(int                      n,             /* Total number of independent observations. */
                      int                      d,             /* Number of independent random variables. */ 
                      FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                      FLOAT                    *h,            /* Sides of the hypersquare. */
                      FLOAT                    nl,            /* Total number of observations in class l. */
                      int                      m,             /* Mode index. */
                      MarginalDistributionType *MrgDistType,  /* Marginal distribution type. */
                      PestraintsType_e         ResType)       /* Restraints type. */  
{
    int                i, j, l, o, p;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, flm, flmin, flmax, V, Dlm, Dlmin;
    int                Error = 0, Stop;

    Mode = (RoughParameterType*)malloc(d * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    /* Rigid restraints. */

    flm = (FLOAT)1.0; V = (FLOAT)1.0;

    for (i = 0; i < d; i++) {
        V *= h[i];

        if (d > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < n; j++) {
                for (l = 0; l < d; l++) if ((i != l) && ((FLOAT)fabs(Y[j][l] - Y[m][l]) > (FLOAT)0.5 * h[l])) goto S0;

                Mode[i].klm += Y[j][d];
S0:;        }
        }
        else
            Mode[i].klm = nl;

        Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][d] * Y[m][d + 1] / (Mode[i].klm * h[i]); flm *= Mode[i].flm;
    }

    epsilon = (FLOAT)exp(log(Y[m][d] * Y[m][d + 1] / (nl * V * flm)) / d);

    for (i = 0; i < d; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].flm *= epsilon;

        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfGamma:
            Error = RoughGammaParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].flm, MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfPoisson:
            Error = RoughPoissonParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0);

            if (Error) goto E0;

            break;
        case pfDirac:
            MrgDistType[i].Par0 = Mode[i].ym;
        }
    }

    if (ResType == rtRigid) goto E0;

    /* Loose restraints. */

    for (i = 0; i < d; i++) {
        if (MrgDistType[i].ParFamType == pfDirac) goto E1;

        /* Bracketing. */

        Dlm = (FLOAT)0.0;

        for (o = 0; o < n; o++) if (Y[o][d] > FLOAT_MIN) {
            for (p = 0; p < d; p++) if ((i != p) && ((FLOAT)fabs(Y[o][p] - Y[m][p]) > (FLOAT)0.5 * h[p])) goto S1;

            Error = ComponentMarginalDist(i, Y[o], MrgDistType, &CmpMrgDist, 0);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * h[i] / Y[o][d + 1];
S1:;    }

        Dlm += (FLOAT)0.998;

        if (Dlm > (FLOAT)0.0) goto E1; 

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)0.998;
        flmax = (FLOAT)Mode[i].flm;

        /* Bisection. */

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            switch (MrgDistType[i].ParFamType) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;
                
                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfGamma:
                Error = RoughGammaParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].ym, flm, MrgDistType[i].Par0, &MrgDistType[i].Par1);
                                                       
                if (Error) goto E0;

                break;
            case pfPoisson:
                Error = RoughPoissonParameters(Mode[i].ym, flm, &MrgDistType[i].Par0);
                                                       
                if (Error) goto E0;

                break;
            case pfDirac:
                break;
            }

            Dlm = (FLOAT)0.0; 

            for (o = 0; o < n; o++) if (Y[o][d] > FLOAT_MIN) {
                for (p = 0; p < d; p++) if ((i != p) && ((FLOAT)fabs(Y[o][p] - Y[m][p]) > (FLOAT)0.5 * h[p])) goto S2;

                Error = ComponentMarginalDist(i, Y[o], MrgDistType, &CmpMrgDist, 0);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * h[i] / Y[o][d + 1];
S2:;        }

            Dlm += (FLOAT)0.998;

            if (((FLOAT)fabs(Dlm) < Eps) || ((flmax - flmin) / (FLOAT)2.0 < Eps)) {
                Stop = 1;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm; 
                }
                else {
                    flmax = flm;
                }
            }
        }
E1:;}

E0: if (Mode) free(Mode);

    return (Error);
} /* RoughEstimationPW */

/* Rough component parameter estimation for histogram. */

int RoughEstimationH(int                      k,             /* Total number of bins. */
                     int                      d,             /* Number of independent random variables. */ 
                     FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl]. */
                     FLOAT                    *h,            /* Sides of the hypersquare. */
                     FLOAT                    nl,            /* Total number of observations in class l. */
                     int                      m,             /* Mode index. */
                     MarginalDistributionType *MrgDistType,  /* Marginal distribution type. */
                     PestraintsType_e         ResType)       /* Restraints type. */  
{
    int                i, j, l, o, p;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, flm, flmin, flmax, V, Dlm, Dlmin;
    int                Error = 0, Stop;

    Mode = (RoughParameterType*)malloc(d * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    /* Rigid restraints. */

    flm = (FLOAT)1.0; V = (FLOAT)1.0;

    for (i = 0; i < d; i++) {
        V *= h[i];

        if (d > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < k; j++) {
                for (l = 0; l < d; l++) if ((i != l) && (Y[j][l] != Y[m][l])) goto S0;

                Mode[i].klm += Y[j][d];
S0:;        }
        }
        else
            Mode[i].klm = nl;

        Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][d] / (Mode[i].klm * h[i]); flm *= Mode[i].flm;
    }

    epsilon = (FLOAT)exp(log(Y[m][d] / (nl * V * flm)) / d);

    for (i = 0; i < d; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].flm *= epsilon;

        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfGamma:
            Error = RoughGammaParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].flm, MrgDistType[i].Par0, &MrgDistType[i].Par1);

            if (Error) goto E0;

            break;
        case pfPoisson:
            Error = RoughPoissonParameters(Mode[i].ym, Mode[i].flm, &MrgDistType[i].Par0);

            if (Error) goto E0;

            break;
        case pfDirac:
            MrgDistType[i].Par0 = Mode[i].ym;
        }
    }

    if (ResType == rtRigid) goto E0;

    /* Loose restraints. */

    for (i = 0; i < d; i++) {
        if (MrgDistType[i].ParFamType == pfDirac) goto E1;

        /* Bracketing. */

        Dlm = (FLOAT)0.0; 

        for (o = 0; o < k; o++) if (Y[o][d] > FLOAT_MIN) {
            for (p = 0; p < d; p++) if ((i != p) && (Y[o][p] != Y[m][p])) goto S1;

            Error = ComponentMarginalDist(i, Y[o], MrgDistType, &CmpMrgDist, 0);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * h[i];
S1:;    }

        Dlm += (FLOAT)0.998;

        if (Dlm > (FLOAT)0.0) goto E1; 

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)0.998;
        flmax = (FLOAT)Mode[i].flm;

        /* Bisection. */

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            switch (MrgDistType[i].ParFamType) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;
                
                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfGamma:
                Error = RoughGammaParameters(Mode[i].ym, flm, &MrgDistType[i].Par0, &MrgDistType[i].Par1);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].ym, flm, MrgDistType[i].Par0, &MrgDistType[i].Par1);
                                                       
                if (Error) goto E0;

                break;
            case pfPoisson:
                Error = RoughPoissonParameters(Mode[i].ym, flm, &MrgDistType[i].Par0);
                                                       
                if (Error) goto E0;

                break;
            case pfDirac:
                break;
            }

            Dlm = (FLOAT)0.0; 

            for (o = 0; o < k; o++) if (Y[o][d] > FLOAT_MIN) {
                for (p = 0; p < d; p++) if ((i != p) && (Y[o][p] != Y[m][p])) goto S2;

                Error = ComponentMarginalDist(i, Y[o], MrgDistType, &CmpMrgDist, 0);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * h[i];
S2:;        }

            Dlm += (FLOAT)0.998;

            if (((FLOAT)fabs(Dlm) < Eps) || ((flmax - flmin) / (FLOAT)2.0 < Eps)) {
                Stop = 1;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm; 
                }
                else {
                    flmax = flm;
                }
            }
        }
E1:;}

E0: if (Mode) free(Mode);

    return (Error);
} /* RoughEstimationH */

/* Enhanced component parameter estimation for k-nearest neighbours. */

int EnhancedEstimationKNN(int                      n,             /* Total number of independent observations. */
                          int                      d,             /* Number of independent random variables. */ 
                          FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,V,R]. */
                          FLOAT                    nl,            /* Total number of observations in class l. */
                          MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    A[4], T[2];
    int                      i, j, l;
    FLOAT                    dP, MrgVar, TmpVar;
    int                      Error = 0;

    TmpDistType = (MarginalDistributionType*)calloc(d, sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            TmpDistType[i].ParFamType = pfNormal;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                TmpDistType[i].Par0 += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Par0 /= nl;

            for(j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] = Y[j][i] - TmpDistType[i].Par0; 

                TmpDistType[i].Par1 += Y[j][d] * T[0] * T[0];
            }

            TmpDistType[i].Par1 /= nl;

            if (TmpDistType[i].Par1 <= FLOAT_MIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Par1 = (FLOAT)sqrt(TmpDistType[i].Par1);

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par1;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par1;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfLognormal:
            TmpDistType[i].ParFamType = pfLognormal;

            for (j = 0; j < n; j++) {
                if ((Y[j][d] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][d] * (FLOAT)log(Y[j][i]);

                    TmpDistType[i].Par0 += T[0]; 
                    TmpDistType[i].Par1 += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            TmpDistType[i].Par0 /= nl; 
            TmpDistType[i].Par1 = TmpDistType[i].Par1 / nl - TmpDistType[i].Par0 * TmpDistType[i].Par0;

            if (TmpDistType[i].Par1 <= FLOAT_MIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Par1 = (FLOAT)sqrt(TmpDistType[i].Par1);

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par1;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par1;

            TmpVar = ((FLOAT)exp(TmpVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * TmpDistType[i].Par0 + TmpVar);
            MrgVar = ((FLOAT)exp(MrgVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * MrgDistType[i].Par0 + MrgVar);

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfWeibull:
            TmpDistType[i].ParFamType = pfWeibull;

            TmpDistType[i].Par1 = (FLOAT)1.0;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < n; l++) {
                    if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * TmpDistType[i].Par1);

                        A[0] += Y[l][d] * T[0];
                        A[1] += Y[l][d] * T[1] * T[0];
                        A[2] += Y[l][d] * T[1];
                        A[3] += Y[l][d] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = TmpDistType[i].Par1 * TmpDistType[i].Par1;

                dP = ((FLOAT)1.0 / TmpDistType[i].Par1 + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                TmpDistType[i].Par1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Par1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Par0 = (FLOAT)exp(log(A[2]) / TmpDistType[i].Par1);

            if ((TmpDistType[i].Par0 <= FLOAT_MIN) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0 * TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par0 * MrgDistType[i].Par0;

            TmpVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / TmpDistType[i].Par1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / TmpDistType[i].Par1)); 
            MrgVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / MrgDistType[i].Par1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / MrgDistType[i].Par1)); 

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfGamma:
            TmpDistType[i].ParFamType = pfGamma;

            TmpDistType[i].Par1 = (FLOAT)1.0 + Eps;

            memset(&A, 0, 2 * sizeof(FLOAT));

            for (l = 0; l < n; l++) {
                if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                    A[0] += Y[l][d] * Y[l][i];
                    A[1] += Y[l][d] * (FLOAT)log(Y[l][i]);
                }
            }

            A[0] /= nl; A[1] /= nl; 

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                if (Digamma(TmpDistType[i].Par1, &T[0]) || Digamma(TmpDistType[i].Par1 + Eps, &T[1])) goto E0;

                dP = ((FLOAT)log(TmpDistType[i].Par1) - T[0] - (FLOAT)log(A[0]) + A[1]) / ((FLOAT)1.0 / TmpDistType[i].Par1 - (T[1] - T[0]) / Eps);

                TmpDistType[i].Par1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Par1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Par0 = A[0] / TmpDistType[i].Par1;

            if ((TmpDistType[i].Par0 <= FLOAT_MIN) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par0 * TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par0 * MrgDistType[i].Par0;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfBinomial:
            TmpDistType[i].ParFamType = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }
            
            TmpDistType[i].Par0 = MrgDistType[i].Par0;

            TmpDistType[i].Par1 = T[0] / TmpDistType[i].Par0 / nl;

            if ((TmpDistType[i].Par0 < (FLOAT)0.0) || (TmpDistType[i].Par1 < (FLOAT)0.0) || (TmpDistType[i].Par1 > (FLOAT)1.0)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0 * TmpDistType[i].Par1 * ((FLOAT)1.0 - TmpDistType[i].Par1);
            MrgVar = MrgDistType[i].Par0 * MrgDistType[i].Par1 * ((FLOAT)1.0 - MrgDistType[i].Par1);

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfPoisson:
            TmpDistType[i].ParFamType = pfPoisson;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Par0 = T[0] / nl;

            TmpDistType[i].Par1 = (FLOAT)0.0;

            if (TmpDistType[i].Par0 < (FLOAT)0.0) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par0;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfDirac:
            TmpDistType[i].ParFamType = pfDirac;

            TmpDistType[i].Par0 = MrgDistType[i].Par0;
        }
    }

    memmove(MrgDistType, TmpDistType, d * sizeof(MarginalDistributionType));

E0: if (TmpDistType) free(TmpDistType);

    return (Error);
} /* EnhancedEstimationKNN */

/* Enhanced component parameter estimation for Parzen window. */

int EnhancedEstimationPW(int                      n,             /* Total number of independent observations. */
                         int                      d,             /* Number of independent random variables. */ 
                         FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                         FLOAT                    nl,            /* Total number of observations in class l. */
                         MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    A[4], T[2];
    int                      i, j, l;
    FLOAT                    dP, MrgVar, TmpVar;
    int                      Error = 0;

    TmpDistType = (MarginalDistributionType*)calloc(d, sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            TmpDistType[i].ParFamType = pfNormal;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                TmpDistType[i].Par0 += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Par0 /= nl;

            for(j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] = Y[j][i] - TmpDistType[i].Par0; 

                TmpDistType[i].Par1 += Y[j][d] * T[0] * T[0];
            }

            TmpDistType[i].Par1 /= nl;

            if (TmpDistType[i].Par1 <= FLOAT_MIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Par1 = (FLOAT)sqrt(TmpDistType[i].Par1);

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par1;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par1;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfLognormal:
            TmpDistType[i].ParFamType = pfLognormal;

            for (j = 0; j < n; j++) {
                if ((Y[j][d] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][d] * (FLOAT)log(Y[j][i]);

                    TmpDistType[i].Par0 += T[0]; 
                    TmpDistType[i].Par1 += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            TmpDistType[i].Par0 /= nl; 
            TmpDistType[i].Par1 = TmpDistType[i].Par1 / nl - TmpDistType[i].Par0 * TmpDistType[i].Par0;

            if (TmpDistType[i].Par1 <= FLOAT_MIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Par1 = (FLOAT)sqrt(TmpDistType[i].Par1);

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par1;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par1;

            TmpVar = ((FLOAT)exp(TmpVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * TmpDistType[i].Par0 + TmpVar);
            MrgVar = ((FLOAT)exp(MrgVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * MrgDistType[i].Par0 + MrgVar);

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfWeibull:
            TmpDistType[i].ParFamType = pfWeibull;

            TmpDistType[i].Par1 = (FLOAT)1.0;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < n; l++) {
                    if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * TmpDistType[i].Par1);

                        A[0] += Y[l][d] * T[0];
                        A[1] += Y[l][d] * T[1] * T[0];
                        A[2] += Y[l][d] * T[1];
                        A[3] += Y[l][d] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = TmpDistType[i].Par1 * TmpDistType[i].Par1;

                dP = ((FLOAT)1.0 / TmpDistType[i].Par1 + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                TmpDistType[i].Par1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Par1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Par0 = (FLOAT)exp(log(A[2]) / TmpDistType[i].Par1);

            if ((TmpDistType[i].Par0 <= FLOAT_MIN) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0 * TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par0 * MrgDistType[i].Par0;

            TmpVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / TmpDistType[i].Par1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / TmpDistType[i].Par1)); 
            MrgVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / MrgDistType[i].Par1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / MrgDistType[i].Par1)); 

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfGamma:
            TmpDistType[i].ParFamType = pfGamma;

            TmpDistType[i].Par1 = (FLOAT)1.0 + Eps;

            memset(&A, 0, 2 * sizeof(FLOAT));

            for (l = 0; l < n; l++) {
                if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                    A[0] += Y[l][d] * Y[l][i];
                    A[1] += Y[l][d] * (FLOAT)log(Y[l][i]);
                }
            }

            A[0] /= nl; A[1] /= nl; 

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                if (Digamma(TmpDistType[i].Par1, &T[0]) || Digamma(TmpDistType[i].Par1 + Eps, &T[1])) goto E0;

                dP = ((FLOAT)log(TmpDistType[i].Par1) - T[0] - (FLOAT)log(A[0]) + A[1]) / ((FLOAT)1.0 / TmpDistType[i].Par1 - (T[1] - T[0]) / Eps);

                TmpDistType[i].Par1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Par1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Par0 = A[0] / TmpDistType[i].Par1;

            if ((TmpDistType[i].Par0 <= FLOAT_MIN) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par0 * TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par0 * MrgDistType[i].Par0;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfBinomial:
            TmpDistType[i].ParFamType = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }
            
            TmpDistType[i].Par0 = MrgDistType[i].Par0;

            TmpDistType[i].Par1 = T[0] / TmpDistType[i].Par0 / nl;

            if ((TmpDistType[i].Par0 < (FLOAT)0.0) || (TmpDistType[i].Par1 < (FLOAT)0.0) || (TmpDistType[i].Par1 > (FLOAT)1.0)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0 * TmpDistType[i].Par1 * ((FLOAT)1.0 - TmpDistType[i].Par1);
            MrgVar = MrgDistType[i].Par0 * MrgDistType[i].Par1 * ((FLOAT)1.0 - MrgDistType[i].Par1);

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfPoisson:
            TmpDistType[i].ParFamType = pfPoisson;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Par0 = T[0] / nl;

            TmpDistType[i].Par1 = (FLOAT)0.0;

            if (TmpDistType[i].Par0 < (FLOAT)0.0) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par0;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfDirac:
            TmpDistType[i].ParFamType = pfDirac;

            TmpDistType[i].Par0 = MrgDistType[i].Par0;
        }
    }

    memmove(MrgDistType, TmpDistType, d * sizeof(MarginalDistributionType));

E0: if (TmpDistType) free(TmpDistType);

    return (Error);
} /* EnhancedEstimationPW */

/* Enhanced component parameter estimation for histogram. */

int EnhancedEstimationH(int                      k,             /* Total number of bins. */
                        int                      d,             /* Number of independent random variables. */ 
                        FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                        FLOAT                    nl,            /* Total number of observations in class l. */
                        MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    A[4], T[2];
    int                      i, j, l;
    FLOAT                    dP, MrgVar, TmpVar;
    int                      Error = 0;

    TmpDistType = (MarginalDistributionType*)calloc(d, sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            TmpDistType[i].ParFamType = pfNormal;

            for (j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                TmpDistType[i].Par0 += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Par0 /= nl;

            for(j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] = Y[j][i] - TmpDistType[i].Par0; 

                TmpDistType[i].Par1 += Y[j][d] * T[0] * T[0];
            }

            TmpDistType[i].Par1 /= nl;

            if (TmpDistType[i].Par1 <= FLOAT_MIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Par1 = (FLOAT)sqrt(TmpDistType[i].Par1);

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par1;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par1;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }
            
            break;
        case pfLognormal:
            TmpDistType[i].ParFamType = pfLognormal;

            for (j = 0; j < k; j++) {
                if ((Y[j][d] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][d] * (FLOAT)log(Y[j][i]);

                    TmpDistType[i].Par0 += T[0]; 
                    TmpDistType[i].Par1 += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            TmpDistType[i].Par0 /= nl; 
            TmpDistType[i].Par1 = TmpDistType[i].Par1 / nl - TmpDistType[i].Par0 * TmpDistType[i].Par0;

            if (TmpDistType[i].Par1 <= FLOAT_MIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Par1 = (FLOAT)sqrt(TmpDistType[i].Par1);

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par1;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par1;

            TmpVar = ((FLOAT)exp(TmpVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * TmpDistType[i].Par0 + TmpVar);
            MrgVar = ((FLOAT)exp(MrgVar) - (FLOAT)1.0) * (FLOAT)exp((FLOAT)2.0 * MrgDistType[i].Par0 + MrgVar);

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfWeibull:
            TmpDistType[i].ParFamType = pfWeibull;

            TmpDistType[i].Par1 = (FLOAT)1.0;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < k; l++) {
                    if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * TmpDistType[i].Par1);

                        A[0] += Y[l][d] * T[0];
                        A[1] += Y[l][d] * T[1] * T[0];
                        A[2] += Y[l][d] * T[1];
                        A[3] += Y[l][d] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = TmpDistType[i].Par1 * TmpDistType[i].Par1;

                dP = ((FLOAT)1.0 / TmpDistType[i].Par1 + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                TmpDistType[i].Par1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Par1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Par0 = (FLOAT)exp(log(A[2]) / TmpDistType[i].Par1);

            if ((TmpDistType[i].Par0 <= FLOAT_MIN) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                Error = 1; if (Error) goto E0;
            }
            
            TmpVar = TmpDistType[i].Par0 * TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par0 * MrgDistType[i].Par0;

            TmpVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / TmpDistType[i].Par1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / TmpDistType[i].Par1)); 
            MrgVar *= (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / MrgDistType[i].Par1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / MrgDistType[i].Par1)); 

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfGamma:
            TmpDistType[i].ParFamType = pfGamma;

            TmpDistType[i].Par1 = (FLOAT)1.0 + Eps;

            memset(&A, 0, 2 * sizeof(FLOAT));

            for (l = 0; l < k; l++) {
                if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                    A[0] += Y[l][d] * Y[l][i];
                    A[1] += Y[l][d] * (FLOAT)log(Y[l][i]);
                }
            }

            A[0] /= nl; A[1] /= nl; 

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                if (Digamma(TmpDistType[i].Par1, &T[0]) || Digamma(TmpDistType[i].Par1 + Eps, &T[1])) goto E0;

                dP = ((FLOAT)log(TmpDistType[i].Par1) - T[0] - (FLOAT)log(A[0]) + A[1]) / ((FLOAT)1.0 / TmpDistType[i].Par1 - (T[1] - T[0]) / Eps);

                TmpDistType[i].Par1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Par1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Par0 = A[0] / TmpDistType[i].Par1;

            if ((TmpDistType[i].Par0 <= FLOAT_MIN) || (TmpDistType[i].Par1 <= FLOAT_MIN)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par1 * TmpDistType[i].Par0 * TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par1 * MrgDistType[i].Par0 * MrgDistType[i].Par0;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfBinomial:
            TmpDistType[i].ParFamType = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Par0 = MrgDistType[i].Par0;

            TmpDistType[i].Par1 = T[0] / TmpDistType[i].Par0 / nl;

            if ((TmpDistType[i].Par0 < (FLOAT)0.0) || (TmpDistType[i].Par1 < (FLOAT)0.0) || (TmpDistType[i].Par1 > (FLOAT)1.0)) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0 * TmpDistType[i].Par1 * ((FLOAT)1.0 - TmpDistType[i].Par1);
            MrgVar = MrgDistType[i].Par0 * MrgDistType[i].Par1 * ((FLOAT)1.0 - MrgDistType[i].Par1);

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfPoisson:
            TmpDistType[i].ParFamType = pfPoisson;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Par0 = T[0] / nl;

            TmpDistType[i].Par1 = (FLOAT)0.0;

            if (TmpDistType[i].Par0 < (FLOAT)0.0) {
                Error = 1; if (Error) goto E0;
            }

            TmpVar = TmpDistType[i].Par0;
            MrgVar = MrgDistType[i].Par0;

            if (TmpVar < MrgVar) {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfDirac:
            TmpDistType[i].ParFamType = pfDirac;

            TmpDistType[i].Par0 = MrgDistType[i].Par0;
        }
    }

    memmove(MrgDistType, TmpDistType, d * sizeof(MarginalDistributionType));

E0: if (TmpDistType) free(TmpDistType);

    return (Error);
} /* EnhancedEstimationH */

/* Moments calculation. */

int MomentsCalculation(int                      d,            /* Number of independent random variables. */ 
                       MarginalDistributionType *MrgDistType, /* Marginal distribution type. */
                       FLOAT                    *FirstM,      /* First moment. */
                       FLOAT                    *SecondM)     /* Second moment. */
{
    int i;
    int Error = 0;
    
    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParFamType) {
        case pfNormal:
            FirstM[i] = MrgDistType[i].Par0;

            SecondM[i] = MrgDistType[i].Par1 * MrgDistType[i].Par1 + MrgDistType[i].Par0 * MrgDistType[i].Par0;

            break;
        case pfLognormal:
            FirstM[i] = (FLOAT)exp(MrgDistType[i].Par0 + (FLOAT)0.5 * MrgDistType[i].Par1 * MrgDistType[i].Par1); 
            
            SecondM[i] = (FLOAT)exp((FLOAT)2.0 * (MrgDistType[i].Par0 + MrgDistType[i].Par1 * MrgDistType[i].Par1));

            break;
        case pfWeibull:
            FirstM[i] = MrgDistType[i].Par0 * (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)1.0 / MrgDistType[i].Par1));

            SecondM[i] = MrgDistType[i].Par0 * MrgDistType[i].Par0 * (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / MrgDistType[i].Par1));

            break;
        case pfGamma:
            FirstM[i] = MrgDistType[i].Par1 * MrgDistType[i].Par0;

            SecondM[i] = MrgDistType[i].Par1 * MrgDistType[i].Par0 * MrgDistType[i].Par0 * ((FLOAT)1.0 + MrgDistType[i].Par1);

            break;
        case pfBinomial:
            FirstM[i] = MrgDistType[i].Par0 * MrgDistType[i].Par1;

            SecondM[i] = (FLOAT)0.0;

            break;
        case pfPoisson:
            FirstM[i] = MrgDistType[i].Par0;

            SecondM[i] = (FLOAT)0.0;

            break;
        case pfDirac:
            FirstM[i] = MrgDistType[i].Par0;

            SecondM[i] = (FLOAT)0.0;
        }
    }

    return(Error);
} /* MomentsCalculation */

/* Returns Bayes Weibull parameters. */

int BayesWeibullParameters(FLOAT                    FirstM,        /* First moment. */
                           FLOAT                    SecondM,       /* Second moment. */
                           MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    FLOAT A;
    FLOAT xl, xm, xh, fl, fm, fh, dx;
    FLOAT i;
    int   Error = 0;
    
    A = (FLOAT)log(SecondM / FirstM / FirstM); xl = 0.001; xh = (FLOAT)10.0;

    fl = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xl) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xl);
    fh = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xh) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xh);

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        if (fl * fh < (FLOAT)0.0)
            Error = 0;
        else
        if ((FLOAT)fabs(fl) < (FLOAT)fabs(fh)) {
            xl += (FLOAT)1.6 * (xl - xh);
            fl = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xl) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xl);
        }
        else {
            xh += (FLOAT)1.6 * (xh - xl);
            fh = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xh) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xh);
        }

        i++;
    }

    if (Error) goto E0;

    /* Root must be bracketed for bisection. */

    if (fl < (FLOAT)0.0) {
        MrgDistType->Par1 = xl; dx = xh - xl;
    }
    else {
        MrgDistType->Par1 = xh; dx = xl - xh;
    }

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        dx = (FLOAT)0.5 * dx; xm = MrgDistType->Par1 + dx;

        fm = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xm) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xm);

        if (fm <= FLOAT_MIN) MrgDistType->Par1 = xm;

        if ((FLOAT)fabs(dx / MrgDistType->Par1) < Eps) Error = 0;

        i++;
    }

    if (Error) goto E0;

    MrgDistType->Par0 = FirstM / (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)1.0 / MrgDistType->Par1));

E0: return (Error);
} /* BayesWeibullParameters */

/* Bayes classification of the remaining observations for k-nearest neighbour. */

int BayesClassificationKNN(int                      n,             /* Total number of independent observations. */
                           int                      d,             /* Number of independent random variables. */
                           FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1]. */
                           int                      c,             /* Number of components. */ 
                           FLOAT                    *W,            /* Component weights. */
                           MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                           FLOAT                    **FirstM,      /* First moments. */
                           FLOAT                    **SecondM)     /* Second moments. */
{
    int   i, j, l;
    FLOAT CmpDist, Max, Tmp, dW;
    int   Error = 0;

    for (i = 0; i < n; i++) {
        if (Y[i][d] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(d, Y[i], MrgDistType[l], &CmpDist, 0);

            if (Error) goto E0;

            Max = W[l] * CmpDist;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) { 
                    l = j; Max = Tmp; 
                }
            }

            dW = Y[i][d] / n; W[l] += dW;

            for (j = 0; j < d; j++) {
                FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                SecondM[l][j] += dW * (Y[i][j] * Y[i][j] - SecondM[l][j]) / W[l];
            }
        }
    }

    for (i = 0; i < c; i++) for (j = 0; j < d; j++) { 
        switch (MrgDistType[i][j].ParFamType) {
        case pfNormal: 
            MrgDistType[i][j].Par0 = FirstM[i][j]; 
            
            MrgDistType[i][j].Par1 = (FLOAT)sqrt(SecondM[i][j] - MrgDistType[i][j].Par0 * MrgDistType[i][j].Par0);

            break;
        case pfLognormal: 
            MrgDistType[i][j].Par0 = (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]) - (FLOAT)0.5 * (FLOAT)log(SecondM[i][j]);
            
           
            MrgDistType[i][j].Par1 = (FLOAT)sqrt(log(SecondM[i][j]) - (FLOAT)2.0 * log(FirstM[i][j]));

            break;
        case pfWeibull:
            BayesWeibullParameters(FirstM[i][j], SecondM[i][j], &MrgDistType[i][j]);

            break;
        case pfGamma:
            MrgDistType[i][j].Par1 = (FLOAT)1.0 / (SecondM[i][j] / FirstM[i][j] / FirstM[i][j] - (FLOAT)1.0);

            MrgDistType[i][j].Par0 = FirstM[i][j] / MrgDistType[i][j].Par1;

            break;
        case pfBinomial:
            MrgDistType[i][j].Par1 = FirstM[i][j] / MrgDistType[i][j].Par0;

            break;
        case pfPoisson:
            MrgDistType[i][j].Par0 = FirstM[i][j];

            break;
        case pfDirac:
            break;
        }
    }

E0: return (Error);
} /* BayesClassificationKNN */

/* Bayes classification of the remaining observations for Parzen window. */

int BayesClassificationPW(int                      n,             /* Total number of independent observations. */
                          int                      d,             /* Number of independent random variables. */
                          FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1]. */
                          int                      c,             /* Number of components. */ 
                          FLOAT                    *W,            /* Component weights. */
                          MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                          FLOAT                    **FirstM,      /* First moments. */
                          FLOAT                    **SecondM)     /* Second moments. */
{
    int   i, j, l;
    FLOAT CmpDist, Max, Tmp, dW;
    int   Error = 0;

    for (i = 0; i < n; i++) {
        if (Y[i][d] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(d, Y[i], MrgDistType[l], &CmpDist, 0);

            if (Error) goto E0;

            Max = W[l] * CmpDist;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) { 
                    l = j; Max = Tmp; 
                }
            }

            dW = Y[i][d] / n; W[l] += dW;

            for (j = 0; j < d; j++) {
                FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                SecondM[l][j] += dW * (Y[i][j] * Y[i][j] - SecondM[l][j]) / W[l];
            }
        }
    }

    for (i = 0; i < c; i++) for (j = 0; j < d; j++) { 
        switch (MrgDistType[i][j].ParFamType) {
        case pfNormal: 
            MrgDistType[i][j].Par0 = FirstM[i][j]; 
            
            MrgDistType[i][j].Par1 = (FLOAT)sqrt(SecondM[i][j] - MrgDistType[i][j].Par0 * MrgDistType[i][j].Par0);

            break;
        case pfLognormal: 
            MrgDistType[i][j].Par0 = (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]) - (FLOAT)0.5 * (FLOAT)log(SecondM[i][j]);
            
           
            MrgDistType[i][j].Par1 = (FLOAT)sqrt(log(SecondM[i][j]) - (FLOAT)2.0 * log(FirstM[i][j]));

            break;
        case pfWeibull:
            BayesWeibullParameters(FirstM[i][j], SecondM[i][j], &MrgDistType[i][j]);

            break;
        case pfGamma:
            MrgDistType[i][j].Par1 = (FLOAT)1.0 / (SecondM[i][j] / FirstM[i][j] / FirstM[i][j] - (FLOAT)1.0);

            MrgDistType[i][j].Par0 = FirstM[i][j] / MrgDistType[i][j].Par1;

            break;
        case pfBinomial:
            MrgDistType[i][j].Par1 = FirstM[i][j] / MrgDistType[i][j].Par0;

            break;
        case pfPoisson:
            MrgDistType[i][j].Par0 = FirstM[i][j];

            break;
        case pfDirac:
            break;
        }
    }

E0: return (Error);
} /* BayesClassificationPW */

/* Bayes classification of the remaining observations for histogram. */

int BayesClassificationH(int                      k,             /* Total number of bins. */
                         int                      n,             /* Total number of independent observations. */
                         int                      d,             /* Number of independent random variables. */
                         FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1]. */
                         int                      c,             /* Number of components. */ 
                         FLOAT                    *W,            /* Component weights. */
                         MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                         FLOAT                    **FirstM,      /* First moments. */
                         FLOAT                    **SecondM)     /* Second moments. */
{
    int   i, j, l;
    FLOAT CmpDist, Max, Tmp, dW;
    int   Error = 0;

    for (i = 0; i < k; i++) {
        if (Y[i][d] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(d, Y[i], MrgDistType[l], &CmpDist, 0);

            if (Error) goto E0;

            Max = W[l] * CmpDist;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) { 
                    l = j; Max = Tmp; 
                }
            }

            dW = Y[i][d] / n; W[l] += dW;

            for (j = 0; j < d; j++) {
                FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                SecondM[l][j] += dW * (Y[i][j] * Y[i][j] - SecondM[l][j]) / W[l];
            }
        }
    }

    for (i = 0; i < c; i++) for (j = 0; j < d; j++) { 
        switch (MrgDistType[i][j].ParFamType) {
        case pfNormal: 
            MrgDistType[i][j].Par0 = FirstM[i][j]; 
            
            MrgDistType[i][j].Par1 = (FLOAT)sqrt(SecondM[i][j] - MrgDistType[i][j].Par0 * MrgDistType[i][j].Par0);

            break;
        case pfLognormal: 
            MrgDistType[i][j].Par0 = (FLOAT)2.0 * (FLOAT)log(FirstM[i][j]) - (FLOAT)0.5 * (FLOAT)log(SecondM[i][j]);
            
           
            MrgDistType[i][j].Par1 = (FLOAT)sqrt(log(SecondM[i][j]) - (FLOAT)2.0 * log(FirstM[i][j]));

            break;
        case pfWeibull:
            BayesWeibullParameters(FirstM[i][j], SecondM[i][j], &MrgDistType[i][j]);

            break;
        case pfGamma:
            MrgDistType[i][j].Par1 = (FLOAT)1.0 / (SecondM[i][j] / FirstM[i][j] / FirstM[i][j] - (FLOAT)1.0);

            MrgDistType[i][j].Par0 = FirstM[i][j] / MrgDistType[i][j].Par1;

            break;
        case pfBinomial:
            MrgDistType[i][j].Par1 = FirstM[i][j] / MrgDistType[i][j].Par0;

            break;
        case pfPoisson:
            MrgDistType[i][j].Par0 = FirstM[i][j];

            break;
        case pfDirac:
            break;
        }
    }

E0: return (Error);
} /* BayesClassificationH */

/* REBMIX algorithm for k-nearest neighbours. */

int REBMIXKNN(InputREBMIXParameterType   *InpParType,  /* Input parameters. */ 
              OutputREBMIXParameterType  *OutParType,  /* Output parameters. */
              HistoryREBMIXParameterType *HisParType)  /* Output parameters. */ 
{
    FLOAT                      **Y = NULL;
    FLOAT                      *h = NULL;
    FLOAT                      *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                      *W = NULL;
    MarginalDistributionType   **Theta = NULL; 
    FLOAT                      **FirstM = NULL, **SecondM = NULL;
    HistoryREBMIXParameterType TmpParType;
    int                        c = 0, i, I, j, J, l, m, M;
    FLOAT                      Dmin, r, nl, elp, eln, epsilonlmax, fl, Dl, f, IC, logL, D;
    #if (_TIME_LEFT_SWITCH)
    clock_t                    Start;
    FLOAT                      TimeLeft;
    #endif
    int                        Error = 0, Stop = 0, Found = 0;

    /* InpParType allocation and initialisation. */

    if (InpParType->ymin == NULL) {
        InpParType->ymin = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == InpParType->ymin; if (Error) goto E0;

        for (i = 0; i < InpParType->d; i++) {
            InpParType->ymin[i] = OutParType->X[0][i];

            for (j = 1; j < OutParType->n; j++) {
                if (OutParType->X[j][i] < InpParType->ymin[i]) InpParType->ymin[i] = OutParType->X[j][i];
            }
        }
    }

    if (InpParType->ymax == NULL) {
        InpParType->ymax = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == InpParType->ymax; if (Error) goto E0;

        for (i = 0; i < InpParType->d; i++) {
            InpParType->ymax[i] = OutParType->X[0][i];

            for (j = 1; j < OutParType->n; j++) {
                if (OutParType->X[j][i] > InpParType->ymax[i]) InpParType->ymax[i] = OutParType->X[j][i];
            }
        }
    }

    /* OutParType allocation and initialisation. */

    OutParType->IC = FLOAT_MAX;

    OutParType->h = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == OutParType->h; if (Error) goto E0;

    OutParType->y0 = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == OutParType->y0; if (Error) goto E0;

    OutParType->W = (FLOAT*)malloc(InpParType->cmax * sizeof(FLOAT));

    Error = NULL == OutParType->W; if (Error) goto E0;

    OutParType->Theta = (MarginalDistributionType**)malloc(InpParType->cmax * sizeof(MarginalDistributionType*));

    Error = NULL == OutParType->Theta; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        OutParType->Theta[i] = (MarginalDistributionType*)malloc(InpParType->d * sizeof(MarginalDistributionType));

        Error = NULL == OutParType->Theta[i]; if (Error) goto E0;
    }

    /* HisParType allocation and initialisation. */ 

    HisParType->Imax = ItMax;

    HisParType->c = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == HisParType->c; if (Error) goto E0;

    HisParType->IC = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->IC; if (Error) goto E0;

    HisParType->logL = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->logL; if (Error) goto E0;

    HisParType->D = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->D; if (Error) goto E0;

    /* Allocation and initialisation. */

    Y = (FLOAT**)malloc(OutParType->n * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < OutParType->n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType->d + 3) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;

        for (j = 0; j < InpParType->d; j++) Y[i][j] = OutParType->X[i][j];
    }

    h = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    for (i = 0; i < InpParType->d; i++) {
        h[i] = InpParType->ymax[i] - InpParType->ymin[i];
    }

    R = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    W = (FLOAT*)malloc(InpParType->cmax * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    Theta = (MarginalDistributionType**)malloc(InpParType->cmax * sizeof(MarginalDistributionType*));

    Error = NULL == Theta; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        Theta[i] = (MarginalDistributionType*)calloc(InpParType->d, sizeof(MarginalDistributionType));

        Error = NULL == Theta[i]; if (Error) goto E0;

        for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].ParFamType = InpParType->IniFamType[j];
        }

        if (InpParType->Ini0 != NULL) for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].Par0 = InpParType->Ini0[j];
        }

        if (InpParType->Ini1 != NULL) for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].Par1 = InpParType->Ini1[j];
        }
    }

    FirstM = (FLOAT**)malloc(InpParType->cmax * sizeof(FLOAT*));

    Error = NULL == FirstM; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        FirstM[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == FirstM[i]; if (Error) goto E0;
    }

    SecondM = (FLOAT**)malloc(InpParType->cmax * sizeof(FLOAT*));

    Error = NULL == SecondM; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        SecondM[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == SecondM[i]; if (Error) goto E0;
    }

    memset(&TmpParType, 0, sizeof(HistoryREBMIXParameterType));

    TmpParType.Imax = ItMax;

    TmpParType.c = (int*)malloc(TmpParType.Imax * sizeof(int));

    Error = NULL == TmpParType.c; if (Error) goto E0;

    TmpParType.IC = (FLOAT*)malloc(TmpParType.Imax * sizeof(FLOAT));

    Error = NULL == TmpParType.IC; if (Error) goto E0;

    TmpParType.logL = (FLOAT*)malloc(TmpParType.Imax * sizeof(FLOAT));

    Error = NULL == TmpParType.logL; if (Error) goto E0;

    TmpParType.D = (FLOAT*)malloc(TmpParType.Imax * sizeof(FLOAT));

    Error = NULL == TmpParType.D; if (Error) goto E0;

    #if (_TIME_LEFT_SWITCH)
    Start = clock();
    #endif

    for (i = 0; i < InpParType->kmax; i++) {
        /* Preprocessing of observations. */

        Error = PreprocessingKNN(InpParType->K[i], h, OutParType->n, InpParType->d, Y);

        if (Error) goto E0;

        Found = 0; Dmin = (FLOAT)0.25; J = 1;

        /* Outer loop. */

        do {
            l = 0; r = (FLOAT)OutParType->n; nl = (FLOAT)OutParType->n;

            /* Middle loop. */

            while (nl / OutParType->n > (FLOAT)2.0 * Dmin * l) {
                /* Global mode detection. */

                Error = GlobalModeKNN(&m, OutParType->n, InpParType->d, Y);

                if (Error) goto E0;

                I = 1; W[l] = nl / OutParType->n; memset(R, 0, OutParType->n * sizeof(FLOAT));

                /* Inner loop. */

                while (I <= ItMax) {
                    /* Rough component parameter estimation. */

                    Error = RoughEstimationKNN(OutParType->n, InpParType->d, Y, InpParType->K[i], h, nl, m, Theta[l], InpParType->ResType);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0;

                    for (j = 0; j < OutParType->n; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][InpParType->d] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = ComponentDist(InpParType->d, Y[j], Theta[l], &fl, 0);

                            if (Error) goto E0;

                            E[j] = Y[j][InpParType->d] - nl * fl * Y[j][InpParType->d + 1] / InpParType->K[i];

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][InpParType->d]; 
                                
                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j]; 
                                
                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl; epsilonlmax -= InpParType->ar;

                    if (Dl <= Dmin / W[l]) {
                        /* Enhanced component parameter estimation. */

                        EnhancedEstimationKNN(OutParType->n, InpParType->d, Y, nl, Theta[l]);

                        break;
                    }
                    else {
                        for (j = 0; j < OutParType->n; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][InpParType->d] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < OutParType->n; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][InpParType->d] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / OutParType->n;
                    }

                    I++;
                } 

                /* Moments calculation. */

                Error = MomentsCalculation(InpParType->d, Theta[l], FirstM[l], SecondM[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < OutParType->n; j++) Y[j][InpParType->d] = R[j];

                Stop = (c >= OutParType->n) || (c >= InpParType->cmax);

                if (Stop) break;
            }

            /* Bayes classification of the remaining observations. */
            
            Error = BayesClassificationKNN(OutParType->n, InpParType->d, Y, c, W, Theta, FirstM, SecondM);

            if (Error) goto E0;

            for (j = 0; j < OutParType->n; j++) Y[j][InpParType->d] = (FLOAT)1.0;

            Error = InformationCriterionKNN(InpParType->ICType, InpParType->K[i], OutParType->n, InpParType->d, Y, c, W, Theta, &IC, &logL, &M, &D);
            
            if (Error) goto E0;

            if (IC < OutParType->IC) {
                Found = 1;

                OutParType->k = InpParType->K[i];

                memmove(OutParType->h, h, InpParType->d * sizeof(FLOAT));  
                
                OutParType->IC = IC; OutParType->logL = logL; OutParType->M = M; OutParType->c = c; 

                memmove(OutParType->W, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    memmove(OutParType->Theta[j], Theta[j], InpParType->d * sizeof(MarginalDistributionType));  
                }
            }

            j = J - 1; TmpParType.c[j] = c; TmpParType.IC[j] = IC; TmpParType.logL[j] = logL; TmpParType.D[j] = D;

            Stop |= (D <= InpParType->D) || (D <= FLOAT_MIN) || (J >= ItMax); Dmin *= c / (c + (FLOAT)1.0); J++;
        }
        while (!Stop);

        TmpParType.Imax = J - 1;

        if (Found) {
            HisParType->Imax = TmpParType.Imax;

            memmove(HisParType->c, TmpParType.c, HisParType->Imax * sizeof(int));  
            memmove(HisParType->IC, TmpParType.IC, HisParType->Imax * sizeof(FLOAT));  
            memmove(HisParType->logL, TmpParType.logL, HisParType->Imax * sizeof(FLOAT));  
            memmove(HisParType->D, TmpParType.D, HisParType->Imax * sizeof(FLOAT));  
        }

        #if (_TIME_LEFT_SWITCH)
        ProgressLength = (int)strlen(Progress); for (j = 0; j < ProgressLength; j++) Progress[j] = ' '; Progress[ProgressLength] = '\0';

        #if (_REBMIXEXE)
        printf("\r%s\r", Progress); 
        #elif (_REBMIXR)
        Rprintf("\r%s\r", Progress);
        R_FlushConsole();
        R_ProcessEvents();
        #endif
        
        TimeLeft = (FLOAT)(InpParType->kmax - i) * (clock() - Start) / CLOCKS_PER_SEC / (i + 1);

        sprintf(Progress, "Time left %2.1f sec", TimeLeft);

        #if (_REBMIXEXE)
        printf("%s", Progress); 
        #elif (_REBMIXR)
        Rprintf("%s", Progress);
        R_FlushConsole();
        R_ProcessEvents();
        #endif
        #endif
    }

E0:;

    #if (_TIME_LEFT_SWITCH)
    ProgressLength = (int)strlen(Progress); for (j = 0; j < ProgressLength; j++) Progress[j] = ' '; Progress[ProgressLength] = '\0';

    #if (_REBMIXEXE)
    printf("\r%s\r", Progress); 
    #elif (_REBMIXR)
    Rprintf("\r%s\r", Progress);
    R_FlushConsole();
    R_ProcessEvents();
    #endif
    #endif

    if (TmpParType.D) free(TmpParType.D);

    if (TmpParType.logL) free(TmpParType.logL);

    if (TmpParType.IC) free(TmpParType.IC);

    if (TmpParType.c) free(TmpParType.c);
    
    if (SecondM) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (SecondM[i]) free(SecondM[i]);
        }
         
        free(SecondM);
    }

    if (FirstM) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (FirstM[i]) free(FirstM[i]);
        }
         
        free(FirstM);
    }

    if (Theta) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (Theta[i]) free(Theta[i]);
        }
         
        free(Theta);
    }

    if (W) free(W);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if (h) free(h);

    if (Y) {
        for (i = 0; i < OutParType->n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    return (Error);
} /* REBMIXKNN */

/* REBMIX algorithm for Parzen window. */

int REBMIXPW(InputREBMIXParameterType   *InpParType,  /* Input parameters. */ 
             OutputREBMIXParameterType  *OutParType,  /* Output parameters. */
             HistoryREBMIXParameterType *HisParType)  /* History parameters. */ 
{
    FLOAT                      **Y = NULL;
    FLOAT                      *h = NULL;
    FLOAT                      *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                      *W = NULL;
    MarginalDistributionType   **Theta = NULL; 
    FLOAT                      **FirstM = NULL, **SecondM = NULL;
    HistoryREBMIXParameterType TmpParType;
    int                        c = 0, i, I, j, J, l, m, M;
    FLOAT                      V, Dmin, r, nl, elp, eln, epsilonlmax, fl, Dl, f, IC, logL, D;
    #if (_TIME_LEFT_SWITCH)
    clock_t                    Start;
    FLOAT                      TimeLeft;
    #endif
    int                        Error = 0, Stop = 0, Found = 0;

    /* InpParType allocation and initialisation. */

    if (InpParType->ymin == NULL) {
        InpParType->ymin = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == InpParType->ymin; if (Error) goto E0;

        for (i = 0; i < InpParType->d; i++) {
            InpParType->ymin[i] = OutParType->X[0][i];

            for (j = 1; j < OutParType->n; j++) {
                if (OutParType->X[j][i] < InpParType->ymin[i]) InpParType->ymin[i] = OutParType->X[j][i];
            }
        }
    }

    if (InpParType->ymax == NULL) {
        InpParType->ymax = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == InpParType->ymax; if (Error) goto E0;

        for (i = 0; i < InpParType->d; i++) {
            InpParType->ymax[i] = OutParType->X[0][i];

            for (j = 1; j < OutParType->n; j++) {
                if (OutParType->X[j][i] > InpParType->ymax[i]) InpParType->ymax[i] = OutParType->X[j][i];
            }
        }
    }

    /* OutParType allocation and initialisation. */

    OutParType->IC = FLOAT_MAX;

    OutParType->h = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == OutParType->h; if (Error) goto E0;

    OutParType->y0 = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == OutParType->y0; if (Error) goto E0;

    OutParType->W = (FLOAT*)malloc(InpParType->cmax * sizeof(FLOAT));

    Error = NULL == OutParType->W; if (Error) goto E0;

    OutParType->Theta = (MarginalDistributionType**)malloc(InpParType->cmax * sizeof(MarginalDistributionType*));

    Error = NULL == OutParType->Theta; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        OutParType->Theta[i] = (MarginalDistributionType*)malloc(InpParType->d * sizeof(MarginalDistributionType));

        Error = NULL == OutParType->Theta[i]; if (Error) goto E0;
    }

    /* HisParType allocation and initialisation. */ 

    HisParType->Imax = ItMax;

    HisParType->c = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == HisParType->c; if (Error) goto E0;

    HisParType->IC = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->IC; if (Error) goto E0;

    HisParType->logL = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->logL; if (Error) goto E0;

    HisParType->D = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->D; if (Error) goto E0;

    /* Allocation and initialisation. */

    Y = (FLOAT**)malloc(OutParType->n * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < OutParType->n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType->d + 2) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;

        for (j = 0; j < InpParType->d; j++) Y[i][j] = OutParType->X[i][j];
    }

    h = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    R = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    W = (FLOAT*)malloc(InpParType->cmax * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    Theta = (MarginalDistributionType**)malloc(InpParType->cmax * sizeof(MarginalDistributionType*));

    Error = NULL == Theta; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        Theta[i] = (MarginalDistributionType*)calloc(InpParType->d, sizeof(MarginalDistributionType));

        Error = NULL == Theta[i]; if (Error) goto E0;

        for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].ParFamType = InpParType->IniFamType[j];
        }

        if (InpParType->Ini0 != NULL) for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].Par0 = InpParType->Ini0[j];
        }

        if (InpParType->Ini1 != NULL) for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].Par1 = InpParType->Ini1[j];
        }
    }

    FirstM = (FLOAT**)malloc(InpParType->cmax * sizeof(FLOAT*));

    Error = NULL == FirstM; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        FirstM[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == FirstM[i]; if (Error) goto E0;
    }

    SecondM = (FLOAT**)malloc(InpParType->cmax * sizeof(FLOAT*));

    Error = NULL == SecondM; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        SecondM[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == SecondM[i]; if (Error) goto E0;
    }

    memset(&TmpParType, 0, sizeof(HistoryREBMIXParameterType));

    TmpParType.Imax = ItMax;

    TmpParType.c = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == TmpParType.c; if (Error) goto E0;

    TmpParType.IC = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == TmpParType.IC; if (Error) goto E0;

    TmpParType.logL = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == TmpParType.logL; if (Error) goto E0;

    TmpParType.D = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == TmpParType.D; if (Error) goto E0;

    #if (_TIME_LEFT_SWITCH)
    Start = clock();
    #endif

    for (i = 0; i < InpParType->kmax; i++) {
        /* Preprocessing of observations. */

        V = (FLOAT)1.0;
        
        for (j = 0; j < InpParType->d; j++) {
            switch (InpParType->VarType[j]) {
            case vtContinuous:
                h[j] = (InpParType->ymax[j] - InpParType->ymin[j]) / InpParType->K[i]; V *= h[j]; 

                break;
            case vtDiscrete:
                h[j] = (FLOAT)1.0;
            }
        }

        Error = PreprocessingPW(h, OutParType->n, InpParType->d, Y);

        if (Error) goto E0;

        Found = 0; Dmin = (FLOAT)0.25; J = 1;

        /* Outer loop. */

        do {
            l = 0; r = (FLOAT)OutParType->n; nl = (FLOAT)OutParType->n;

            /* Middle loop. */

            while (nl / OutParType->n > (FLOAT)2.0 * Dmin * l) {
                /* Global mode detection. */

                Error = GlobalModePW(&m, OutParType->n, InpParType->d, Y);

                if (Error) goto E0;

                I = 1; W[l] = nl / OutParType->n; memset(R, 0, OutParType->n * sizeof(FLOAT));

                /* Inner loop. */

                while (I <= ItMax) {
                    /* Rough component parameter estimation. */

                    Error = RoughEstimationPW(OutParType->n, InpParType->d, Y, h, nl, m, Theta[l], InpParType->ResType);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0;

                    for (j = 0; j < OutParType->n; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][InpParType->d] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = ComponentDist(InpParType->d, Y[j], Theta[l], &fl, 0);

                            if (Error) goto E0;

                            E[j] = Y[j][InpParType->d] - nl * fl * V / Y[j][InpParType->d + 1];

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][InpParType->d];

                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j]; 
                                
                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl; epsilonlmax -= InpParType->ar;

                    if (Dl <= Dmin / W[l]) {
                        /* Enhanced component parameter estimation. */

                        EnhancedEstimationPW(OutParType->n, InpParType->d, Y, nl, Theta[l]);

                        break;
                    }
                    else {
                        for (j = 0; j < OutParType->n; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][InpParType->d] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < OutParType->n; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][InpParType->d] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / OutParType->n;
                    }

                    I++;
                }

                /* Moments calculation. */

                Error = MomentsCalculation(InpParType->d, Theta[l], FirstM[l], SecondM[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < OutParType->n; j++) Y[j][InpParType->d] = R[j];

                Stop = (c >= OutParType->n) || (c >= InpParType->cmax);

                if (Stop) break;
            }

            /* Bayes classification of the remaining observations. */

            Error = BayesClassificationPW(OutParType->n, InpParType->d, Y, c, W, Theta, FirstM, SecondM);

            if (Error) goto E0;

            for (j = 0; j < OutParType->n; j++) Y[j][InpParType->d] = (FLOAT)1.0;

            Error = InformationCriterionPW(InpParType->ICType, V, OutParType->n, InpParType->d, Y, c, W, Theta, &IC, &logL, &M, &D);
            
            if (Error) goto E0;

            if (IC < OutParType->IC) {
                Found = 1;

                OutParType->k = InpParType->K[i];

                memmove(OutParType->h, h, InpParType->d * sizeof(FLOAT));  
                
                OutParType->IC = IC; OutParType->logL = logL; OutParType->M = M; OutParType->c = c; 

                memmove(OutParType->W, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    memmove(OutParType->Theta[j], Theta[j], InpParType->d * sizeof(MarginalDistributionType));  
                }
            }

            j = J - 1; TmpParType.c[j] = c; TmpParType.IC[j] = IC; TmpParType.logL[j] = logL; TmpParType.D[j] = D;

            Stop |= (D <= InpParType->D) || (D <= FLOAT_MIN) || (J >= ItMax); Dmin *= c / (c + (FLOAT)1.0); J++; 
        }
        while (!Stop);

        TmpParType.Imax = J - 1;

        if (Found) {
            HisParType->Imax = TmpParType.Imax;

            memmove(HisParType->c, TmpParType.c, HisParType->Imax * sizeof(int));  
            memmove(HisParType->IC, TmpParType.IC, HisParType->Imax * sizeof(FLOAT));  
            memmove(HisParType->logL, TmpParType.logL, HisParType->Imax * sizeof(FLOAT));  
            memmove(HisParType->D, TmpParType.D, HisParType->Imax * sizeof(FLOAT));  
        }

        #if (_TIME_LEFT_SWITCH)
        ProgressLength = (int)strlen(Progress); for (j = 0; j < ProgressLength; j++) Progress[j] = ' '; Progress[ProgressLength] = '\0';

        #if (_REBMIXEXE)
        printf("\r%s\r", Progress); 
        #elif (_REBMIXR)
        Rprintf("\r%s\r", Progress);
        R_FlushConsole();
        R_ProcessEvents();
        #endif
        
        TimeLeft = (FLOAT)(InpParType->kmax - i) * (clock() - Start) / CLOCKS_PER_SEC / (i + 1);

        sprintf(Progress, "Time left %2.1f sec", TimeLeft);

        #if (_REBMIXEXE)
        printf("%s", Progress); 
        #elif (_REBMIXR)
        Rprintf("%s", Progress);
        R_FlushConsole();
        R_ProcessEvents();
        #endif
        #endif
    }

E0:; 
    
    #if (_TIME_LEFT_SWITCH)
    ProgressLength = (int)strlen(Progress); for (j = 0; j < ProgressLength; j++) Progress[j] = ' '; Progress[ProgressLength] = '\0';

    #if (_REBMIXEXE)
    printf("\r%s\r", Progress); 
    #elif (_REBMIXR)
    Rprintf("\r%s\r", Progress);
    R_FlushConsole();
    R_ProcessEvents();
    #endif
    #endif

    if (TmpParType.D) free(TmpParType.D);

    if (TmpParType.logL) free(TmpParType.logL);

    if (TmpParType.IC) free(TmpParType.IC);

    if (TmpParType.c) free(TmpParType.c);
    
    if (SecondM) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (SecondM[i]) free(SecondM[i]);
        }
         
        free(SecondM);
    }

    if (FirstM) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (FirstM[i]) free(FirstM[i]);
        }
         
        free(FirstM);
    }

    if (Theta) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (Theta[i]) free(Theta[i]);
        }
         
        free(Theta);
    }

    if (W) free(W);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if (h) free(h);

    if (Y) {
        for (i = 0; i < OutParType->n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    return (Error);
} /* REBMIXPW */

/* REBMIX algorithm for histogram. */

int REBMIXH(InputREBMIXParameterType   *InpParType,  /* Input parameters. */ 
            OutputREBMIXParameterType  *OutParType,  /* Output parameters. */
            HistoryREBMIXParameterType *HisParType)  /* History parameters. */ 
{
    FLOAT                      **Y = NULL;
    FLOAT                      *h = NULL, *y0 = NULL;
    FLOAT                      *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                      *K = NULL;
    FLOAT                      *W = NULL;
    MarginalDistributionType   **Theta = NULL; 
    FLOAT                      **FirstM = NULL, **SecondM = NULL;
    HistoryREBMIXParameterType TmpParType;
    int                        c = 0, i, I, j, J, k, l, m, M;
    FLOAT                      V, Dmin, r, nl, elp, eln, epsilonlmax, fl, Dl, f, IC, logL, D;
    #if (_TIME_LEFT_SWITCH)
    clock_t                    Start;
    FLOAT                      TimeLeft;
    #endif
    int                        Error = 0, Stop = 0, Found = 0;

    /* InpParType allocation and initialisation. */

    if (InpParType->ymin == NULL) {
        InpParType->ymin = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == InpParType->ymin; if (Error) goto E0;

        for (i = 0; i < InpParType->d; i++) {
            InpParType->ymin[i] = OutParType->X[0][i];

            for (j = 1; j < OutParType->n; j++) {
                if (OutParType->X[j][i] < InpParType->ymin[i]) InpParType->ymin[i] = OutParType->X[j][i];
            }
        }
    }

    if (InpParType->ymax == NULL) {
        InpParType->ymax = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == InpParType->ymax; if (Error) goto E0;

        for (i = 0; i < InpParType->d; i++) {
            InpParType->ymax[i] = OutParType->X[0][i];

            for (j = 1; j < OutParType->n; j++) {
                if (OutParType->X[j][i] > InpParType->ymax[i]) InpParType->ymax[i] = OutParType->X[j][i];
            }
        }
    }

    /* OutParType allocation and initialisation. */

    OutParType->IC = FLOAT_MAX;

    OutParType->h = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == OutParType->h; if (Error) goto E0;

    OutParType->y0 = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == OutParType->y0; if (Error) goto E0;

    OutParType->W = (FLOAT*)malloc(InpParType->cmax * sizeof(FLOAT));

    Error = NULL == OutParType->W; if (Error) goto E0;

    OutParType->Theta = (MarginalDistributionType**)malloc(InpParType->cmax * sizeof(MarginalDistributionType*));

    Error = NULL == OutParType->Theta; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        OutParType->Theta[i] = (MarginalDistributionType*)malloc(InpParType->d * sizeof(MarginalDistributionType));

        Error = NULL == OutParType->Theta[i]; if (Error) goto E0;
    }

    /* HisParType allocation and initialisation. */ 

    HisParType->Imax = ItMax;

    HisParType->c = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == HisParType->c; if (Error) goto E0;

    HisParType->IC = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->IC; if (Error) goto E0;

    HisParType->logL = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->logL; if (Error) goto E0;

    HisParType->D = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == HisParType->D; if (Error) goto E0;

    /* Allocation and initialisation. */

    Y = (FLOAT**)malloc(OutParType->n * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < OutParType->n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType->d + 1) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;
    }

    /* Allocation and normalizing vector h calculation. */

    h = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    y0 = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

    Error = NULL == y0; if (Error) goto E0;

    R = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    K = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == K; if (Error) goto E0;

    W = (FLOAT*)malloc(InpParType->cmax * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    Theta = (MarginalDistributionType**)malloc(InpParType->cmax * sizeof(MarginalDistributionType*));

    Error = NULL == Theta; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        Theta[i] = (MarginalDistributionType*)calloc(InpParType->d, sizeof(MarginalDistributionType));

        Error = NULL == Theta[i]; if (Error) goto E0;

        for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].ParFamType = InpParType->IniFamType[j];
        }

        if (InpParType->Ini0 != NULL) for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].Par0 = InpParType->Ini0[j];
        }

        if (InpParType->Ini1 != NULL) for (j = 0; j < InpParType->d; j++) {
            Theta[i][j].Par1 = InpParType->Ini1[j];
        }
    }

    FirstM = (FLOAT**)malloc(InpParType->cmax * sizeof(FLOAT*));

    Error = NULL == FirstM; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        FirstM[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == FirstM[i]; if (Error) goto E0;
    }

    SecondM = (FLOAT**)malloc(InpParType->cmax * sizeof(FLOAT*));

    Error = NULL == SecondM; if (Error) goto E0;

    for (i = 0; i < InpParType->cmax; i++) {
        SecondM[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));

        Error = NULL == SecondM[i]; if (Error) goto E0;
    }

    memset(&TmpParType, 0, sizeof(HistoryREBMIXParameterType));

    TmpParType.Imax = ItMax;

    TmpParType.c = (int*)malloc(ItMax * sizeof(int));

    Error = NULL == TmpParType.c; if (Error) goto E0;

    TmpParType.IC = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == TmpParType.IC; if (Error) goto E0;

    TmpParType.logL = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == TmpParType.logL; if (Error) goto E0;

    TmpParType.D = (FLOAT*)malloc(ItMax * sizeof(FLOAT));

    Error = NULL == TmpParType.D; if (Error) goto E0;

    #if (_TIME_LEFT_SWITCH)
    Start = clock();
    #endif

    for (i = 0; i < InpParType->kmax; i++) {
        /* Preprocessing of observations. */

        k = InpParType->K[i]; V = (FLOAT)1.0; 
        
        for (j = 0; j < InpParType->d; j++) {
            switch (InpParType->VarType[j]) {
            case vtContinuous:
                h[j] = (InpParType->ymax[j] - InpParType->ymin[j]) / InpParType->K[i]; y0[j] = InpParType->ymin[j] + (FLOAT)0.5 * h[j]; V *= h[j]; 

                break;
            case vtDiscrete:
                h[j] = (FLOAT)1.0; y0[j] = InpParType->ymin[j];
            }
        }

        Error = PreprocessingH(h, y0, InpParType->VarType, &InpParType->K[i], OutParType->n, InpParType->d, OutParType->X, Y);

        if (Error) goto E0;

        for (j = 0; j < InpParType->K[i]; j++) K[j] = Y[j][InpParType->d];

        Found = 0; Dmin = (FLOAT)0.25; J = 1;

        /* Outer loop. */

        do {
            l = 0; r = (FLOAT)OutParType->n; nl = (FLOAT)OutParType->n;

            /* Middle loop. */

            while (nl / OutParType->n > (FLOAT)2.0 * Dmin * l) {
                /* Global mode detection. */

                Error = GlobalModeH(&m, InpParType->K[i], InpParType->d, Y);

                if (Error) goto E0;

                I = 1; W[l] = nl / OutParType->n; memset(R, 0, InpParType->K[i] * sizeof(FLOAT));

                /* Inner loop. */

                while (I <= ItMax) { 
                    /* Rough component parameter estimation. */

                    Error = RoughEstimationH(InpParType->K[i], InpParType->d, Y, h, nl, m, Theta[l], InpParType->ResType);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0;

                    for (j = 0; j < InpParType->K[i]; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][InpParType->d] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = ComponentDist(InpParType->d, Y[j], Theta[l], &fl, 0);

                            if (Error) goto E0;

                            E[j] = Y[j][InpParType->d] - nl * fl * V;

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][InpParType->d]; 
                                
                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j]; 
                                
                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl; epsilonlmax -= InpParType->ar;
                    
                    if (Dl <= Dmin / W[l]) {
                        /* Enhanced component parameter estimation. */

                        EnhancedEstimationH(InpParType->K[i], InpParType->d, Y, nl, Theta[l]);

                        break;
                    }
                    else {
                        for (j = 0; j < InpParType->K[i]; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][InpParType->d] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < InpParType->K[i]; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][InpParType->d] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / OutParType->n;
                    }

                    I++;
                }

                /* Moments calculation. */

                Error = MomentsCalculation(InpParType->d, Theta[l], FirstM[l], SecondM[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < InpParType->K[i]; j++) Y[j][InpParType->d] = R[j];

                Stop = (c >= InpParType->K[i]) || (c >= InpParType->cmax);

                if (Stop) break;
            }

            /* Bayes classification of the remaining observations. */

            Error = BayesClassificationH(InpParType->K[i], OutParType->n, InpParType->d, Y, c, W, Theta, FirstM, SecondM);

            if (Error) goto E0;

            for (j = 0; j < InpParType->K[i]; j++) Y[j][InpParType->d] = K[j];

            Error = InformationCriterionH(InpParType->ICType, V, InpParType->K[i], OutParType->n, InpParType->d, Y, c, W, Theta, &IC, &logL, &M, &D);
            
            if (Error) goto E0;

            if (IC < OutParType->IC) {
                Found = 1; 

                OutParType->k = k;

                memmove(OutParType->h, h, InpParType->d * sizeof(FLOAT));

                memmove(OutParType->y0, y0, InpParType->d * sizeof(FLOAT)); 
                
                OutParType->IC = IC; OutParType->logL = logL; OutParType->M = M; OutParType->c = c; 

                memmove(OutParType->W, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    memmove(OutParType->Theta[j], Theta[j], InpParType->d * sizeof(MarginalDistributionType));  
                }
            }

            j = J - 1; TmpParType.c[j] = c; TmpParType.IC[j] = IC; TmpParType.logL[j] = logL; TmpParType.D[j] = D;

            Stop |= (D <= InpParType->D) || (D <= FLOAT_MIN) || (J >= ItMax); Dmin *= c / (c + (FLOAT)1.0); J++;
        }
        while (!Stop);

        TmpParType.Imax = J - 1;

        if (Found) {
            HisParType->Imax = TmpParType.Imax;

            memmove(HisParType->c, TmpParType.c, HisParType->Imax * sizeof(int));  
            memmove(HisParType->IC, TmpParType.IC, HisParType->Imax * sizeof(FLOAT));  
            memmove(HisParType->logL, TmpParType.logL, HisParType->Imax * sizeof(FLOAT));  
            memmove(HisParType->D, TmpParType.D, HisParType->Imax * sizeof(FLOAT));  
        }

        InpParType->K[i] = k;

        #if (_TIME_LEFT_SWITCH)
        ProgressLength = (int)strlen(Progress); for (j = 0; j < ProgressLength; j++) Progress[j] = ' '; Progress[ProgressLength] = '\0';

        #if (_REBMIXEXE)
        printf("\r%s\r", Progress); 
        #elif (_REBMIXR)
        Rprintf("\r%s\r", Progress);
        R_FlushConsole();
        R_ProcessEvents();
        #endif
        
        TimeLeft = (FLOAT)(InpParType->kmax - i) * (clock() - Start) / CLOCKS_PER_SEC / (i + 1);

        sprintf(Progress, "Time left %2.1f sec", TimeLeft);

        #if (_REBMIXEXE)
        printf("%s", Progress); 
        #elif (_REBMIXR)
        Rprintf("%s", Progress);
        R_FlushConsole();
        R_ProcessEvents();
        #endif
        #endif
    }

E0:;

    #if (_TIME_LEFT_SWITCH)
    ProgressLength = (int)strlen(Progress); for (j = 0; j < ProgressLength; j++) Progress[j] = ' '; Progress[ProgressLength] = '\0';

    #if (_REBMIXEXE)
    printf("\r%s\r", Progress); 
    #elif (_REBMIXR)
    Rprintf("\r%s\r", Progress);
    R_FlushConsole();
    R_ProcessEvents();
    #endif
    #endif

    if (TmpParType.D) free(TmpParType.D);

    if (TmpParType.logL) free(TmpParType.logL);

    if (TmpParType.IC) free(TmpParType.IC);

    if (TmpParType.c) free(TmpParType.c);
    
    if (SecondM) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (SecondM[i]) free(SecondM[i]);
        }
         
        free(SecondM);
    }

    if (FirstM) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (FirstM[i]) free(FirstM[i]);
        }
         
        free(FirstM);
    }

    if (Theta) {
        for (i = 0; i < InpParType->cmax; i++) {
            if (Theta[i]) free(Theta[i]);
        }
         
        free(Theta);
    }

    if (W) free(W);

    if (K) free(K);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if (y0) free(y0);

    if (h) free(h);

    if (Y) {
        for (i = 0; i < OutParType->n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    return (Error);
} /* REBMIXH */

/* Reads input data from the file stream. */

int ReadREBMIXDataFile(InputREBMIXParameterType  *InpParType,  /* Input parameters. */ 
                       OutputREBMIXParameterType *OutParType)  /* Output parameters. */
{
    char line[65536];
    char *pchar = NULL;
    int  i, j, BufSize = 0;
    FILE *fp = NULL;
    int  Error = 0;

    memset(OutParType, 0, sizeof(OutputREBMIXParameterType));

    if ((fp = fopen(InpParType->curr, "r")) == NULL) {
        Error = 1; goto E0;
    }

    OutParType->X = (FLOAT**)malloc((BufSize + BufInc) * sizeof(FLOAT*));

    for (i = BufSize; i < BufSize + BufInc; i++) {
        OutParType->X[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));
    }

    BufSize += BufInc;

    OutParType->n = 0;

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n");

        if (!pchar) goto S0;

        j = 0;
        for (i = 0; i < (int)strlen(pchar); i++) {
            if (pchar[i] == ',') {
                line[j] = '.'; j++;
            }
            else
            if (pchar[i] != ' ') {
                line[j] = pchar[i]; j++;
            }
        }

        line[j] = '\0';

        if (!j) goto S0;

        pchar = strtok(pchar, "\t");

        if (OutParType->n == BufSize) {
            OutParType->X = (FLOAT**)realloc(OutParType->X, (BufSize + BufInc) * sizeof(FLOAT*));

            for (i = BufSize; i < BufSize + BufInc; i++) {
                OutParType->X[i] = (FLOAT*)malloc(InpParType->d * sizeof(FLOAT));
            }

            BufSize += BufInc;
        }

        i = 0;
        while (pchar) {
            OutParType->X[OutParType->n][i] = (FLOAT)atof(pchar); pchar = strtok(NULL, "\t"); i++;
        }

        OutParType->n++;
    }

    for (i = OutParType->n; i < BufSize; i++) {
        if (OutParType->X[i]) free(OutParType->X[i]); 
    }

    OutParType->X = (FLOAT**)realloc(OutParType->X, OutParType->n * sizeof(FLOAT*));

E0: if (fp) fclose(fp);

    return (Error);
} /* ReadREBMIXDataFile */

/* Writes input and output parameters into the file stream. */

int WriteREBMIXParameterFile(InputREBMIXParameterType  *InpParType,  /* Input parameters. */ 
                             OutputREBMIXParameterType *OutParType)  /* Output parameters. */
{
    int  i, j;
    char line[65536];
    char mode[2];
    char path[FILENAME_MAX];
    char ext[FILENAME_MAX];
    char *pchar = NULL;
    FILE *fp0 = NULL, *fp1 = NULL;
    int  Error = 0;

    if (InpParType->curr == InpParType->open[0])
        strcpy(mode, "w");
    else
        strcpy(mode, "a");

    strcpy(path, InpParType->save); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(path, "%s%s%s", path, "_1", ext);

    if ((fp0 = fopen(path, mode)) == NULL) {
        Error = 1; goto E0;
    }

    strcpy(path, InpParType->save); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(path, "%s%s%s", path, "_2", ext);

    if ((fp1 = fopen(path, mode)) == NULL) {
        Error = 1; goto E0;
    }

    if (!strcmp(mode, "w")) {
        fprintf(fp0, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "Dataset",
                                                           "Preprocessing",
                                                           "D",
                                                           "cmax",
                                                           "Criterion",
                                                           "ar",
                                                           "Restraints", 
                                                           "c",
                                                           "b");

        switch (InpParType->PreType) {
        case poHistogram:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < InpParType->d; i++) {
                if (InpParType->d == 1)
                    sprintf(line, "%s", "y0");
                else
                    sprintf(line, "%s%d", "y0", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            for (i = 0; i < InpParType->d; i++) {
                if (InpParType->d == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            break;
        case poParzenWindow:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < InpParType->d; i++) {
                if (InpParType->d == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            break;
        case poKNearestNeighbour:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < InpParType->d; i++) {
                if (InpParType->d == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }
        }

        fprintf(fp0, "\t%s\t%s\n", "IC",
                                   "logL");

        fprintf(fp1, "%s\t%s", "Dataset",
                               "w");

        for (i = 0; i < InpParType->d; i++) {
            switch (InpParType->IniFamType[i]) {
            case pfNormal: case pfLognormal: case pfWeibull: case pfGamma: case pfBinomial:
                if (InpParType->d == 1)
                    fprintf(fp1, "\t%s\t%s\t%s", "pdf", "theta1", "theta2");
                else
                    fprintf(fp1, "\t%s%d\t%s%d\t%s%d", "pdf", i + 1, "theta1.", i + 1, "theta2.", i + 1);

                break;
            case pfPoisson: case pfDirac:
                if (InpParType->d == 1)
                    fprintf(fp1, "\t%s\t%s", "pdf", "theta1");
                else
                    fprintf(fp1, "\t%s%d\t%s%d", "pdf", i + 1, "theta1.", i + 1);
            }
        }

        fprintf(fp1, "\n");
    }

    strcpy(path, InpParType->curr); 

    pchar = strrchr(path, '\\');

    if (!pchar) {
        pchar = strrchr(path, '/');
    }

    if (pchar) {
        strcpy(path, pchar + 1);
    }

    pchar = strrchr(path, '.'); 
        
    if (pchar) pchar[0] = '\0';

    fprintf(fp0, "%s", path);

    switch (InpParType->PreType) {
    case poHistogram: 
        strcpy(line, "histogram");

        break;
    case poParzenWindow: 
        strcpy(line, "Parzen window");

        break;
    case poKNearestNeighbour:
        strcpy(line, "k-nearest neighbour");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%E\t%d", InpParType->D,
                             InpParType->cmax);

    switch (InpParType->ICType) {
    case icAIC:
        strcpy(line, "AIC");

        break;
    case icAIC3:
        strcpy(line, "AIC3");

        break;
    case icAIC4:
        strcpy(line, "AIC4");

        break;
    case icAICc:
        strcpy(line, "AICc");

        break;
    case icBIC:
        strcpy(line, "BIC");

        break;
    case icCAIC:
        strcpy(line, "CAIC");

        break;
    case icHQC:
        strcpy(line, "HQC");

        break;
    case icMDL2:
        strcpy(line, "MDL2");

        break;
    case icMDL5:
        strcpy(line, "MDL5");

        break;
    case icAWE:
        strcpy(line, "AWE");

        break;
    case icCLC:
        strcpy(line, "CLC");

        break;
    case icICL:
        strcpy(line, "ICL");

        break;
    case icPC:
        strcpy(line, "PC");

        break;
    case icICLBIC:
        strcpy(line, "ICL-BIC");

        break;
    case icD:
        strcpy(line, "D");
        
        break;
    case icSSE:
        strcpy(line, "SSE");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%E", InpParType->ar);

    switch (InpParType->ResType) {
    case rtRigid:
        strcpy(line, "rigid");

        break;
    case rtLoose:
        strcpy(line, "loose");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%d", OutParType->c);

    fprintf(fp0, "\t%E", InpParType->b);

    switch (InpParType->PreType) {
    case poHistogram:
        fprintf(fp0, "\t%d", OutParType->k);

        for (i = 0; i < InpParType->d; i++) {
            fprintf(fp0, "\t%E", OutParType->y0[i]);
        }

        for (i = 0; i < InpParType->d; i++) {
            fprintf(fp0, "\t%E", OutParType->h[i]);
        }

        break;
    case poParzenWindow:
        fprintf(fp0, "\t%d", OutParType->k);

        for (i = 0; i < InpParType->d; i++) {
            fprintf(fp0, "\t%E", OutParType->h[i]);
        }

        break;
    case poKNearestNeighbour:
        fprintf(fp0, "\t%d", OutParType->k);

        for (i = 0; i < InpParType->d; i++) {
            fprintf(fp0, "\t%E", OutParType->h[i]);
        }
    }

    fprintf(fp0, "\t%E\t%E\n", OutParType->IC,
                               OutParType->logL);

    for (i = 0; i < OutParType->c; i++) {
        fprintf(fp1, "%s\t%E", path,
                               OutParType->W[i]);

        for (j = 0; j < InpParType->d; j++) switch (OutParType->Theta[i][j].ParFamType) {
            case pfNormal:
                fprintf(fp1, "\t%s\t%E\t%E", "normal", 
                                             OutParType->Theta[i][j].Par0,
                                             OutParType->Theta[i][j].Par1);

                break;
            case pfLognormal:
                fprintf(fp1, "\t%s\t%E\t%E", "lognormal", 
                                             OutParType->Theta[i][j].Par0,
                                             OutParType->Theta[i][j].Par1);
                break;
            case pfWeibull:
                fprintf(fp1, "\t%s\t%E\t%E", "Weibull", 
                                             OutParType->Theta[i][j].Par0,
                                             OutParType->Theta[i][j].Par1);
                break;
            case pfGamma:
                fprintf(fp1, "\t%s\t%E\t%E", "gamma", 
                                             OutParType->Theta[i][j].Par0,
                                             OutParType->Theta[i][j].Par1);
                break;
            case pfBinomial:
                fprintf(fp1, "\t%s\t%E\t%E", "binomial", 
                                             OutParType->Theta[i][j].Par0,
                                             OutParType->Theta[i][j].Par1);
                break;
            case pfPoisson:
                fprintf(fp1, "\t%s\t%E", "Poisson", 
                                         OutParType->Theta[i][j].Par0);
                break;
            case pfDirac:
                fprintf(fp1, "\t%s\t%E", "Dirac", 
                                         OutParType->Theta[i][j].Par0);
        }

        fprintf(fp1, "\n");
    }

E0: if (fp0) fclose(fp0);
    if (fp1) fclose(fp1);

    return (Error);
} /* WriteREBMIXParameterFile */

/* REBMIX algorithm. */

int REBMIX(InputREBMIXParameterType   *InpParType,  /* Input parameters. */ 
           OutputREBMIXParameterType  *OutParType,  /* Output parameters. */
           HistoryREBMIXParameterType *HisParType)  /* History parameters. */ 
{
    int  Error = 0;

    switch (InpParType->PreType) {
    case poHistogram:
        Error = REBMIXH(InpParType, OutParType, HisParType);

        if (Error) goto E0;

        break;
    case poParzenWindow:
        Error = REBMIXPW(InpParType, OutParType, HisParType);

        if (Error) goto E0;

        break;
    case poKNearestNeighbour:
        Error = REBMIXKNN(InpParType, OutParType, HisParType);

        if (Error) goto E0;
    }

E0: return (Error);
} /* REBMIX */

/* Runs REBMIX template file stream. */

int RunREBMIXTemplateFile(char *file)
{
    char                       line[65536], ident[65536], list[65536];
    char                       *pchar = NULL, *rchar = NULL;
    int                        i, imin, imax, iinc, j, k, isI;
    FLOAT                      isF;
    InputREBMIXParameterType   InpParType;
    OutputREBMIXParameterType  OutParType;
    HistoryREBMIXParameterType HisParType;
    FILE                       *fp = NULL;
    int                        Error = 0;
    #if (_MEMORY_LEAK_SWITCH)
    _CrtMemState               s1, s2, s3;

    _CrtMemCheckpoint(&s1);
    #endif

    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));

    /* Recommended values. */

    InpParType.D = (FLOAT)0.025;
    InpParType.cmax = 15;
    InpParType.ICType = icAIC;
    InpParType.b = (FLOAT)1.0;
    InpParType.ar = (FLOAT)0.1;
    InpParType.ResType = rtLoose;

    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));
    memset(&HisParType, 0, sizeof(HistoryREBMIXParameterType));

    #if (_REBMIXEXE)
    printf("REBMIX Version 2.6.0\n");
    #endif

    if ((fp = fopen(file, "r")) == NULL) {
        Error = 1; goto E0;
    }

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n"); 
        
        pchar = strtok(pchar, "=");

        if (!pchar) goto S0;

        j = 0;

        for (i = 0; i < (int)strlen(pchar); i++) {
            if (pchar[i] != ' ') {
                ident[j] = (char)toupper(pchar[i]); j++;
            }
        }

        ident[j] = '\0';

        j = 0; list[j] = '\0'; imin = 0; imax = 0;

        while((pchar = strtok(NULL, ",")) != NULL) {
            if (!strcmp(ident, "DATASET") || !strcmp(ident, "SAVE")) {
                for (i = 0; pchar[i] != '\0'; i++) {
                    list[j] = pchar[i]; j++;

                    if (pchar[i] != ' ') {
                        imax = i; if (!imin) imin = i;
                    }
                }

                j = imax + 1 - imin;

                for (i = 0; i < j; i++) {
                    list[i] = list[imin + i];
                }
            }
            else {
                for (i = 0; pchar[i] != '\0'; i++) if (pchar[i] != '[' && pchar[i] != ']') {
                    if (pchar[i] == ' ') {
                    }
                    else {
                        list[j] = (char)toupper(pchar[i]); j++;
                    }
                }
            }

            list[j] = '\t'; j++;
        }

        if (!j) goto S0; else list[j - 1] = '\0';

        pchar = strtok(list, "\t");

        if (!strcmp(ident, "RUN")) {
            for (k = 0; k < InpParType.o; k++) {
                InpParType.curr = InpParType.open[k];

                Error = ReadREBMIXDataFile(&InpParType, &OutParType);

                if (Error) goto E0;

                #if (_REBMIXEXE)
                printf("Dataset = %s\n", InpParType.curr);
                #endif

                Error = REBMIX(&InpParType, &OutParType, &HisParType);

                if (Error) goto E0;

                Error = WriteREBMIXParameterFile(&InpParType, &OutParType);

                if (Error) goto E0;
            }
        }
        else
        if (!strcmp(ident, "DATASET")) {
            InpParType.open = (char**)realloc(InpParType.open, (InpParType.o + 1) * sizeof(char*));

            Error = NULL == InpParType.open; if (Error) goto E0;

            InpParType.open[InpParType.o] = (char*)malloc((strlen(pchar) + 1) * sizeof(char));

            Error = NULL == InpParType.open[InpParType.o]; if (Error) goto E0;

            strcpy(InpParType.open[InpParType.o], pchar); InpParType.o++;
        }
        else
        if (!strcmp(ident, "PREPROCESSING")) {
            if (!strcmp(pchar, "HISTOGRAM"))
                InpParType.PreType = poHistogram;
            else
            if (!strcmp(pchar, "PARZENWINDOW"))
                InpParType.PreType = poParzenWindow;
            else
            if (!strcmp(pchar, "K-NEARESTNEIGHBOUR"))
                InpParType.PreType = poKNearestNeighbour;
            else {
                Error = 1; goto E0;
            }
        } else
        if (!strcmp(ident, "D")) {
            InpParType.D = isF = (FLOAT)atof(pchar);

            Error = (isF < (FLOAT)0.0) || (isF > (FLOAT)1.0); if (Error) goto E0;
        } else
        if (!strcmp(ident, "CMAX")) {
            InpParType.cmax = isI = (int)atol(pchar);

            Error = isI <= 0; if (Error) goto E0;
        } else
        if (!strcmp(ident, "CRITERION")) {
            if (!strcmp(pchar, "AIC"))
                InpParType.ICType = icAIC;
            else
            if (!strcmp(pchar, "AIC3"))
                InpParType.ICType = icAIC3;
            else
            if (!strcmp(pchar, "AIC4"))
                InpParType.ICType = icAIC4;
            else
            if (!strcmp(pchar, "AICC"))
                InpParType.ICType = icAICc;
            else
            if (!strcmp(pchar, "BIC"))
                InpParType.ICType = icBIC;
            else
            if (!strcmp(pchar, "CAIC"))
                InpParType.ICType = icCAIC;
            else
            if (!strcmp(pchar, "HQC"))
                InpParType.ICType = icHQC;
            else
            if (!strcmp(pchar, "MDL2"))
                InpParType.ICType = icMDL2;
            else
            if (!strcmp(pchar, "MDL5"))
                InpParType.ICType = icMDL5;
            else
            if (!strcmp(pchar, "AWE"))
                InpParType.ICType = icAWE;
            else
            if (!strcmp(pchar, "CLC"))
                InpParType.ICType = icCLC;
            else
            if (!strcmp(pchar, "ICL"))
                InpParType.ICType = icICL;
            else
            if (!strcmp(pchar, "PC"))
                InpParType.ICType = icPC;
            else
            if (!strcmp(pchar, "ICL-BIC"))
                InpParType.ICType = icICLBIC;
            else
            if (!strcmp(pchar, "D"))
                InpParType.ICType = icD;
            else
            if (!strcmp(pchar, "SSE"))
                InpParType.ICType = icSSE;
            else {
                Error = 1; goto E0;
            }
        } else
        if (!strcmp(ident, "VARIABLES")) {
            i = 0;

            while (pchar) {
                InpParType.VarType = (VariablesType_e*)realloc(InpParType.VarType, (i + 1) * sizeof(VariablesType_e));

                Error = NULL == InpParType.VarType; if (Error) goto E0;

                if (!strcmp(pchar, "CONTINUOUS"))
                    InpParType.VarType[i] = vtContinuous; 
                else
                if (!strcmp(pchar, "DISCRETE"))
                    InpParType.VarType[i] = vtDiscrete;
                else {
                    Error = 1; goto E0;
                }
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "PDF")) {
            i = 0;

            while (pchar) {
                InpParType.IniFamType = (ParametricFamilyType_e*)realloc(InpParType.IniFamType, (i + 1) * sizeof(ParametricFamilyType_e));

                Error = NULL == InpParType.IniFamType; if (Error) goto E0;

                if (!strcmp(pchar, "NORMAL"))
                    InpParType.IniFamType[i] = pfNormal; 
                else
                if (!strcmp(pchar, "LOGNORMAL"))
                    InpParType.IniFamType[i] = pfLognormal;
                else
                if (!strcmp(pchar, "WEIBULL"))
                    InpParType.IniFamType[i] = pfWeibull;
                else
                if (!strcmp(pchar, "GAMMA"))
                    InpParType.IniFamType[i] = pfGamma;
                else
                if (!strcmp(pchar, "BINOMIAL"))
                    InpParType.IniFamType[i] = pfBinomial;
                else
                if (!strcmp(pchar, "POISSON"))
                    InpParType.IniFamType[i] = pfPoisson;
                else
                if (!strcmp(pchar, "DIRAC"))
                    InpParType.IniFamType[i] = pfDirac;
                else {
                    Error = 1; goto E0;
                }
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "THETA1")) {
            i = 0;

            while (pchar) {
                InpParType.Ini0 = (FLOAT*)realloc(InpParType.Ini0, (i + 1) * sizeof(FLOAT));

                Error = NULL == InpParType.Ini0; if (Error) goto E0;

                InpParType.Ini0[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "THETA2")) {
            i = 0;

            while (pchar) {
                InpParType.Ini1 = (FLOAT*)realloc(InpParType.Ini1, (i + 1) * sizeof(FLOAT));

                Error = NULL == InpParType.Ini1; if (Error) goto E0;

                InpParType.Ini1[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "K")) {
            i = 0;

            while (pchar != NULL) {
                if ((rchar = strrchr(pchar, '-')) != NULL) {
                    imin = (int)atol(pchar); imax = (int)atol(rchar + 1);

                    if (imin > imax) {
                        j = imin; imin = imax; imax = j;
                    }

                    if ((rchar = strrchr(pchar, '+')) != NULL)
                        iinc = (int)atol(rchar + 1);
                    else
                        iinc = 1;

                    InpParType.K = (int*)realloc(InpParType.K, (i + (imax - imin) / iinc + 1) * sizeof(int));

                    Error = NULL == InpParType.K; if (Error) goto E0;

                    for (j = imin; j <= imax; j += iinc) {
                        InpParType.K[i] = isI = j; 

                        Error = isI <= 0; if (Error) goto E0;
                        
                        i++;
                    }
                
                    InpParType.kmax = i;
                }
                else {
                    InpParType.K = (int*)realloc(InpParType.K, (i + 1) * sizeof(int));

                    Error = NULL == InpParType.K; if (Error) goto E0;

                    InpParType.K[i] = isI = (int)atol(pchar);

                    Error = isI <= 0; if (Error) goto E0;
                
                    InpParType.kmax = ++i;
                }

                pchar = strtok(NULL, "\t"); 
            }
        } else
        if (!strcmp(ident, "YMIN")) {
            i = 0;

            while (pchar) {
                InpParType.ymin = (FLOAT*)realloc(InpParType.ymin, (i + 1) * sizeof(FLOAT));

                Error = NULL == InpParType.ymin; if (Error) goto E0;

                InpParType.ymin[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "YMAX")) {
            i = 0;

            while (pchar) {
                InpParType.ymax = (FLOAT*)realloc(InpParType.ymax, (i + 1) * sizeof(FLOAT));

                Error = NULL == InpParType.ymax; if (Error) goto E0;

                InpParType.ymax[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "B")) {
            InpParType.b = isF = (FLOAT)atof(pchar);

            Error = (isF < (FLOAT)0.0) || (isF > (FLOAT)1.0); if (Error) goto E0;
        } else
        if (!strcmp(ident, "AR")) {
            InpParType.ar = isF = (FLOAT)atof(pchar);

            Error = (isF <= (FLOAT)0.0) || (isF > (FLOAT)1.0); if (Error) goto E0;

        } else
        if (!strcmp(ident, "RESTRAINTS")) {
            if (!strcmp(pchar, "RIGID"))
                InpParType.ResType = rtRigid;
            else
            if (!strcmp(pchar, "LOOSE"))
                InpParType.ResType = rtLoose;
            else {
                Error = 1; goto E0;
            }
        } else
        if (!strcmp(ident, "SAVE")) {
            InpParType.save = (char*)realloc(InpParType.save, (strlen(pchar) + 1) * sizeof(char));

            Error = NULL == InpParType.save; if (Error) goto E0;

            strcpy(InpParType.save, pchar);
        }
    }

E0: if (fp) fclose(fp);

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

    if (InpParType.save) free(InpParType.save);

    if (InpParType.ymax) free(InpParType.ymax);

    if (InpParType.ymin) free(InpParType.ymin);

    if (InpParType.K) free(InpParType.K);

    if (InpParType.Ini1) free(InpParType.Ini1);

    if (InpParType.Ini0) free(InpParType.Ini0);

    if (InpParType.IniFamType) free(InpParType.IniFamType);

    if (InpParType.VarType) free(InpParType.VarType);
    
    if (InpParType.open) {
        for (i = 0; i < InpParType.o; i++) {
            if (InpParType.open[i]) free(InpParType.open[i]);
        }
         
        free(InpParType.open);
    }

    if (HisParType.D) free(HisParType.D);

    if (HisParType.logL) free(HisParType.logL);

    if (HisParType.IC) free(HisParType.IC);

    if (HisParType.c) free(HisParType.c);

    #if (_MEMORY_LEAK_SWITCH)
    _CrtMemCheckpoint(&s2);

    if (_CrtMemDifference(&s3, &s1, &s2)) _CrtMemDumpStatistics(&s3);
    #endif

    return (Error);
} /* RunREBMIXTemplateFile */