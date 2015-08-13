#include <math.h>
#include <limits.h>

#include "base.h"

#if (_MEMORY_LEAK_SWITCH)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

// Inserts y into ascending list Y of length n. Set n = 0 initially.

void Insert(FLOAT y,   // Inserted value.
            int   *n,  // Length of Y.
            FLOAT *Y)  // Pointer to Y = [y0,...,yn-1].
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
} // Insert

// Returns the value log(Gamma(y)) for y > 0. See http://www.nr.com/.

FLOAT Gammaln(FLOAT y)
{
    FLOAT x, z, Tmp, Ser;
    int   j;

    static FLOAT Cof[6] = { (FLOAT)76.18009172947146, -(FLOAT)86.50532032941677,
        (FLOAT)24.01409824083091, -(FLOAT)1.231739572450155,
        (FLOAT)0.1208650973866179E-2, -(FLOAT)0.5395239384953E-5 };

    static FLOAT Stp = (FLOAT)2.5066282746310005;

    z = x = y; Tmp = x + (FLOAT)5.5; Tmp -= (x + (FLOAT)0.5) * (FLOAT)log(Tmp);

    Ser = (FLOAT)1.000000000190015;

    for (j = 0; j < 6; j++) Ser += Cof[j] / ++z;

    return (-Tmp + (FLOAT)log(Stp * Ser / x));
} // Gammaln

// Returns the digamma for y > 0. See http://www.nr.com/.

int Digamma(FLOAT y, FLOAT *Psi)
{
    static FLOAT piov4 = (FLOAT)0.785398163397448;
    static FLOAT dy0 = (FLOAT)1.461632144968362341262659542325721325;
    static FLOAT p1[7] = { (FLOAT)0.0089538502298197, (FLOAT)4.77762828042627, (FLOAT)142.441585084029,
        (FLOAT)1186.45200713425, (FLOAT)3633.51846806499, (FLOAT)4138.10161269013,
        (FLOAT)1305.60269827897 };
    static FLOAT q1[6] = { (FLOAT)44.8452573429826, (FLOAT)520.752771467162, (FLOAT)2210.0079924783,
        (FLOAT)3641.27349079381, (FLOAT)1908.310765963, (FLOAT)6.91091682714533e-6 };
    static FLOAT p2[4] = { -(FLOAT)2.12940445131011, -(FLOAT)7.01677227766759,
        -(FLOAT)4.48616543918019, -(FLOAT)0.648157123766197 };
    static FLOAT q2[4] = { (FLOAT)32.2703493791143, (FLOAT)89.2920700481861,
        (FLOAT)54.6117738103215, (FLOAT)7.77788548522962 };
    int   i, m, n, nq;
    FLOAT d2;
    FLOAT w, z;
    FLOAT den, aug, sgn, ymy0, ymax, upper, ymin;
    int Error = 0;

    ymax = (FLOAT)INT_MAX; d2 = (FLOAT)1.0 / FLOAT_EPSILON; if (ymax > d2) ymax = d2; ymin = (FLOAT)1E-9; aug = (FLOAT)0.0;

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

            nq = (int)w; w -= (FLOAT)nq; nq = (int)(w * (FLOAT)4.0); w = (w - (FLOAT)nq * (FLOAT)0.25) * (FLOAT)4.0; n = nq / 2;

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
} // Digamma

// Returns the inverse of the binomial c.d.f. for the specified n and p.

FLOAT BinomialInv(FLOAT Fy, FLOAT n, FLOAT p)
{
    FLOAT Sum, y, ypb;

    Sum = ypb = (FLOAT)pow((FLOAT)1.0 - p, n); y = (FLOAT)0.0;

    while ((Sum < Fy) && (ypb > FLOAT_MIN)) {
        y++; ypb *= (n - y + (FLOAT)1.0) * p / y / ((FLOAT)1.0 - p); Sum += ypb;
    }

    if ((Fy < (FLOAT)0.5) && (y >(FLOAT)0.0)) y--;

    return (y);
} // BinomialInv

// Returns the inverse of the Poisson c.d.f. for the specified Theta.

FLOAT PoissonInv(FLOAT Fy, FLOAT Theta)
{
    FLOAT Sum, y, ypb;

    Sum = ypb = (FLOAT)exp(-Theta); y = (FLOAT)0.0;

    while ((Sum < Fy) && (ypb > FLOAT_MIN)) {
        y++; ypb *= Theta / y; Sum += ypb;
    }

    if ((Fy < (FLOAT)0.5) && (y >(FLOAT)0.0)) y--;

    return (y);
} // PoissonInv

// Returns the incomplete gamma function P(a, y) evaluated by its series
// representation as GamSer. Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

int GammaSer(FLOAT a,       // Constant a > 0.
             FLOAT y,       // Variable y > 0.
             FLOAT *GamSer, // Incomplete gamma function.
             FLOAT *Gamln)  // Log(Gamma(a)).
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

            if ((FLOAT)fabs(Del) < FLOAT_EPSILON) Error = 0;

            i++;
        }

        if (Error) goto E0;

        *GamSer = Sum * (FLOAT)exp(-y + a * log(y) - *Gamln);
    }

E0: return (Error);
} // GammaSer

// Returns the incomplete gamma function Q(a, y) evaluated by its continued
// fraction representation as GamCfg. Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

int GammaCfg(FLOAT a,       // Constant a > 0.
             FLOAT y,       // Variable y > 0.
             FLOAT *GamCfg, // Incomplete gamma function.
             FLOAT *Gamln)  // Log(Gamma(a)).
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

                if ((FLOAT)fabs(G - Gold) < FLOAT_EPSILON) Error = 0; else Gold = G;
            }

            i++;
        }

        if (Error) goto E0;

        *GamCfg = (FLOAT)exp(-y + a * log(y) - *Gamln) * G;
    }

E0: return (Error);
} // GammaCfg

// Returns the incomplete gamma function P(a, y). Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

int GammaP(FLOAT a,       // Constant a > 0.
           FLOAT y,       // Variable y > 0.
           FLOAT *GamP,   // Incomplete gamma function.
           FLOAT *Gamln)  // Log(Gamma(a)).
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
} // GammaP

// Returns the inverse of the gamma c.d.f. for the specified Theta and Beta. See http://www.nr.com/.

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

        if ((FLOAT)fabs(dy) < Eps) Error = 0;

        i++;
    }

E0: return (Error);
} // GammaInv

// Returns the error function erf(y). See http://www.nr.com/.

int ErrorF(FLOAT y,     // Variable y.
           FLOAT *ErF)  // Error function.
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
} // ErrorF

