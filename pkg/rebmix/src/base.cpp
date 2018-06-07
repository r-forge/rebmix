#ifdef _MSC_VER
#pragma warning(disable: 4514)
#endif

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "base.h"

#if (_MEMORY_LEAK_SWITCH)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

// Base constructor.

Base::Base()
{
    Trigger_ = 0;
    length_pdf_ = 0;
    length_Theta_ = 0;
    length_theta_ = NULL;
} // Base

// Base destructor.

Base::~Base()
{
    if (length_theta_) free(length_theta_);
} // ~Base

static long IY = 0;
static long IV[NTAB];

// Minimal random number generator of Park and Miller with Bays-Durham shuffle and added
// safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
// values). Call with IDum a negative integer to initialize; thereafter do not alter IDum
// between successive deviates in a sequence.RNMX should approximate the largest floating
// value that is less than 1. See http://www.nrbook.com/a/bookcpdf/c7-1.pdf

FLOAT Ran1(int *IDum)
{
    int   j, k;
    FLOAT Tmp;

    if (*IDum <= 0 || !IY) {
        *IDum = (-(*IDum) < 1) ? 1 : -(*IDum);

        for (j = NTAB + 7; j >= 0; j--) {
            k = *IDum / IQ;

            *IDum = IA * (*IDum - k * IQ) - IR * k;

            if (*IDum < 0) *IDum += IM;

            if (j < NTAB) IV[j] = *IDum;
        }

        IY = IV[0];
    }

    k = *IDum / IQ; *IDum = IA * (*IDum - k * IQ) - IR * k;

    if (*IDum < 0) *IDum += IM;

    j = IY / NDIV; IY = IV[j]; IV[j] = *IDum;

    if ((Tmp = AM * IY) > RNMX) return (RNMX); else return (Tmp);
} // Ran1

// Inserts y into ascending list Y of length n.Set n = 0 initially.

void Insert(FLOAT y,   // Inserted value.
            int   *n,  // Length of Y.
            FLOAT *Y)  // Pointer to Y = [y0,...,yn-1].
{
    int i, j;

    Y[*n] = y;

    for (i = 0; i < *n; i++) {
        if (y < Y[i]) {
            for (j = *n; j > i; j--) Y[j] = Y[j - 1];
            
            Y[i] = y;

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

    if ((Fy < (FLOAT)0.5) && (y > (FLOAT)0.0)) y--;

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
		*GamSer = (FLOAT)0.0;
    }
    else {
        ap = a; Sum = (FLOAT)1.0 / a; Del = Sum; Error = 1; i = 1;

        while ((i <= ItMax) && Error) {
            ap += (FLOAT)1.0; Del *= y / ap; Sum += Del;

            if ((FLOAT)fabs(Del) < Eps) Error = 0;

            i++;
        }

        if (Error) Error = 0; // ItMax too small.

        *GamSer = Sum * (FLOAT)exp(-y + a * (FLOAT)log(y) - *Gamln);
    }

    return (Error);
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

                if ((FLOAT)fabs(G - Gold) < Eps) Error = 0; else Gold = G;
            }

            i++;
        }

        if (Error) Error = 0; // ItMax too small.

        *GamCfg = (FLOAT)exp(-y + a * (FLOAT)log(y) - *Gamln) * G;
    }

    return (Error);
} // GammaCfg

// Returns the incomplete gamma function P(a, y). Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

int GammaP(FLOAT a,       // Constant a > 0.
           FLOAT y,       // Variable y > 0.
           FLOAT *GamP,   // Incomplete gamma function.
           FLOAT *Gamln)  // Log(Gamma(a)).
{
    FLOAT GamSer, GamCfg;
    int   Error = 0;

	if ((y <= FLOAT_MIN) || (a <= FLOAT_MIN)) {
		*GamP = (FLOAT)0.0;
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

	if (Beta > (FLOAT)1.0) {
		*y = (Beta - (FLOAT)1.0) * Theta + Eps;
	}
    else
        *y = Eps;

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        if (GammaP(Beta, *y / Theta, &GamP, &Gamln)) goto E0;

        Tmp = *y / Theta;

        dy = (GamP - Fy) / ((FLOAT)exp(Beta * (FLOAT)log(Tmp) - Tmp - Gamln) / (*y));

        *y -= dy;

        if (IsNan(dy) || IsInf(dy)) {
			Error = 1; goto E0;
        }
        else
        if (*y < Eps) {
            *y = Eps; Error = 0;
        }

        if ((FLOAT)fabs(dy) < Eps) Error = 0;

        i++;
    }

	if (Error) Error = 0; // ItMax too small.

E0: return (Error);
} // GammaInv

// Returns the inverse of the Weibull c.d.f. for the specified Theta and Beta.

FLOAT WeibullInv(FLOAT Fy, FLOAT Theta, FLOAT Beta)
{
    FLOAT y;

    y = Theta * (FLOAT)pow(-(FLOAT)log((FLOAT)1.0 - Fy), (FLOAT)1.0 / Beta);

    return (y);
} // WeibullInv

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

// Returns the LU decomposition of matrix A. See http://www.nr.com/ 

int LUdcmp(int   n,     // Size of square matrix.
           FLOAT *A,    // Pointer to the square matrix A.
           int   *indx, // Pointer to the permutation vector.
           FLOAT *det)  // Determinant.
{
    int   i, imax, j, k;
    FLOAT Big, Tmp;
    FLOAT *V;
    int   Error = 0;

    V = (FLOAT*)malloc(n * sizeof(FLOAT));

    Error = NULL == V; if (Error) goto E0;

    for (i = 0; i < n; i++) {
        Big = (FLOAT)0.0;

        for (j = 0; j < n; j++) {
            if ((Tmp = (FLOAT)fabs(A[i * n + j])) > Big) Big = Tmp;
        }

        if ((FLOAT)fabs(Big) <= FLOAT_MIN) {
            Error = 1; goto E0;
        }

        V[i] = (FLOAT)1.0 / Big;
    }

    *det = (FLOAT)1.0;

    for (k = 0; k < n; k++) {
        Big = (FLOAT)0.0; imax = k;

        for (i = k; i < n; i++) {
            Tmp = V[i] * (FLOAT)fabs(A[i * n + k]);

            if (Tmp > Big) {
                Big = Tmp; imax = i;
            }
        }

        if (k != imax) {
            for (j = 0; j < n; j++) {
                Tmp = A[imax * n + j]; A[imax * n + j] = A[k * n + j]; A[k * n + j] = Tmp;
            }

            *det = -(*det); V[imax] = V[k];
        }

        indx[k] = imax;

        if ((FLOAT)fabs(A[k * n + k]) <= FLOAT_MIN) A[k * n + k] = FLOAT_MIN;

        for (i = k + 1; i < n; i++) {
            Tmp = A[i * n + k] /= A[k * n + k];

            for (j = k + 1; j < n; j++) A[i * n + j] -= Tmp * A[k * n + j];
        }
    }

    for (i = 0; i < n; i++) *det *= A[i * n + i];

    if (IsNan(*det) || ((FLOAT)fabs(*det) <= FLOAT_MIN)) {
        Error = 1; goto E0;
    }

E0: if (V) free(V);

    return (Error);
} // LUdcmp 

// Solves the set of n linear equations A x = b. See See http://www.nr.com/ 

int LUbksb(int   n,     // Size of square matrix.
           FLOAT *A,    // Pointer to the square matrix A.
           int   *indx, // Pointer to the permutation vector.
           FLOAT *b)    // Pointer to the solution vector.
{
    int   i, ii = 0, ip, j;
    FLOAT Sum;
    int   Error = 0;

    for (i = 0; i < n; i++) {
        ip = indx[i]; Sum = b[ip]; b[ip] = b[i];

        if (ii) {
            for (j = ii - 1; j < i; j++) Sum -= A[i * n + j] * b[j];
        }
        else
        if (Sum) {
            ii = i + 1;
        }

        b[i] = Sum;
    }

    for (i = n - 1; i >= 0; i--) {
        Sum = b[i];

        for (j = i + 1; j < n; j++) Sum -= A[i * n + j] * b[j];

        b[i] = Sum / A[i * n + i];
    }

    return (Error);
} // LUbksb 

// Returns the determinant and the inverse matrix of A. See http://www.nr.com/ 

int LUinvdet(int   n,     // Size of square matrix.
             FLOAT *A,    // Pointer to the square matrix A.
             FLOAT *Ainv, // Pointer to the inverse matrix of A.
             FLOAT *Adet) // Pointer to the determinant of A.
{
    int   i, *indx = NULL, j;
    FLOAT *b = NULL, *B = NULL;
    int   Error = 0;

    indx = (int*)calloc((size_t)n, sizeof(int));

    Error = NULL == indx; if (Error) goto E0;

    b = (FLOAT*)malloc(n * sizeof(FLOAT));

    Error = NULL == b; if (Error) goto E0;

    B = (FLOAT*)malloc(n * n * sizeof(FLOAT));

    Error = NULL == B; if (Error) goto E0;

    memmove(B, A, n * n * sizeof(FLOAT));

    Error = LUdcmp(n, B, indx, Adet);

    if (Error) goto E0;

    for (j = 0; j < n; j++) {
        memset(b, 0, n * sizeof(FLOAT));

        b[j] = (FLOAT)1.0;

        Error = LUbksb(n, B, indx, b);

        if (Error) goto E0;

        for (i = 0; i < n; i++) Ainv[i * n + j] = b[i];
    }

E0: if (B) free(B);	
 
    if (b) free(b);

    if (indx) free(indx);

    return (Error);
} // LUinvdet

// Returns the Cholesky decomposition of matrix A. See http://www.nr.com/ 

int Choldc(int   n,   // Size of square matrix.
           FLOAT *A,  // Pointer to the square matrix A.
           FLOAT *L)  // Lower triangular factors.
{
    int   i, j, k;
    FLOAT Sum;
    FLOAT *p = NULL;
    int   Error = 0;

    memmove(L, A, n * n * sizeof(FLOAT));

    p = (FLOAT*)malloc(n * sizeof(FLOAT));

    Error = NULL == p; if (Error) goto E0;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            Sum = L[i * n + j];

            for (k = 0; k < i; k++) Sum -= L[i * n + k] * L[j * n + k];

            if (i == j) {
                if (Sum < Eps) {
                    A[i * n + j] = Eps - Sum; Sum = Eps;
                }

                p[i] = (FLOAT)sqrt(Sum);
            }
            else {
                L[j * n + i] = Sum / p[i];
            }
        }
    }

    for (i = 0; i < n; i++) {
        L[i * n + i] = p[i]; for (j = 0; j < i; j++) L[j * n + i] = (FLOAT)0.0;
    }

    if (p) free(p);

E0: return (Error);
} // Choldc

// Returns the determinant and the inverse matrix of A. See http ://www.nr.com/ 

int Cholinvdet(int   n,         // Size of square matrix.
               FLOAT *A,        // Pointer to the symmetric square matrix A.
               FLOAT *Ainv,     // Pointer to the inverse matrix of A.
               FLOAT *logAdet)  // Pointer to the logarithm of determinant of A.
{
    int   i, j, k;
    FLOAT *L = NULL, Sum;
    FLOAT *p = NULL;
    int   Error = 0;

    L = (FLOAT*)malloc(n * n * sizeof(FLOAT));

    Error = NULL == L; if (Error) goto E0;

    memmove(L, A, n * n * sizeof(FLOAT));

    p = (FLOAT*)malloc(n * sizeof(FLOAT));

    Error = NULL == p; if (Error) goto E0;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            Sum = L[i * n + j];

            for (k = 0; k < i; k++) Sum -= L[i * n + k] * L[j * n + k];

            if (i == j) {
                if (Sum < Eps) {
                    A[i * n + j] = Eps - Sum; Sum = Eps;
                }

                p[i] = (FLOAT)sqrt(Sum);
            }
            else {
                L[j * n + i] = Sum / p[i];
            }
        }
    }

    *logAdet = (FLOAT)0.0;

    for (i = 0; i < n; i++) {
        L[i * n + i] = (FLOAT)1.0 / p[i]; *logAdet += (FLOAT)log(p[i]);

        for (j = i - 1; j >= 0; j--) {
            Sum = (FLOAT)0.0;

            for (k = j; k < i; k++) Sum -= L[i * n + k] * L[j * n + k];

            L[j * n + i] = Sum / p[i];
        }
    }

    *logAdet *= (FLOAT)2.0;

    for (i = 0; i < n; i++) {
       for (j = i; j < n; j++) {
           Sum = (FLOAT)0.0;

           for (k = j; k < n; k++) Sum += L[i * n + k] * L[j * n + k];

           Ainv[i * n + j] = Ainv[j * n + i] = Sum;
       }
    }

E0: if (p) free(p);

    if (L) free(L);

    return (Error);
} // Cholinvdet

// Returns modified Bessel function of order 0. See http://people.math.sfu.ca/~cbm/aands/page_378.htm 

FLOAT BesselI0(FLOAT y)
{
	FLOAT t, I0;

	y = (FLOAT)fabs(y); t =  y / (FLOAT)3.75;

	if (y <= FLOAT_MIN) {
		I0 = (FLOAT)1.0;
	}
	else
	if (y <= (FLOAT)3.75) {
		I0 = (FLOAT)1.0 +
		 	 (FLOAT)3.5156229 * (FLOAT)pow(t, 2) +
			 (FLOAT)3.0899424 * (FLOAT)pow(t, 4) +
			 (FLOAT)1.2067492 * (FLOAT)pow(t, 6) +
			 (FLOAT)0.2659732 * (FLOAT)pow(t, 8) +
			 (FLOAT)0.0360768 * (FLOAT)pow(t, 10) +
			 (FLOAT)0.0045813 * (FLOAT)pow(t, 12);
	}
	else {
		I0 = (FLOAT)0.39894228 +
			 (FLOAT)0.01328592 * (FLOAT)pow(t, -1) +
			 (FLOAT)0.00225319 * (FLOAT)pow(t, -2) -
			 (FLOAT)0.00157565 * (FLOAT)pow(t, -3) +
			 (FLOAT)0.00916281 * (FLOAT)pow(t, -4) -
			 (FLOAT)0.02057706 * (FLOAT)pow(t, -5) +
			 (FLOAT)0.02635537 * (FLOAT)pow(t, -6) -
			 (FLOAT)0.01647633 * (FLOAT)pow(t, -7) +
			 (FLOAT)0.00392377 * (FLOAT)pow(t, -8);

		I0 *= (FLOAT)exp(y) / (FLOAT)sqrt(y);
	}

    return (I0);
} // BesselI0

// Returns modified Bessel function of order 1. See http://people.math.sfu.ca/~cbm/aands/page_378.htm 

FLOAT BesselI1(FLOAT y)
{
	FLOAT t, I1, sgn;

	if (y < (FLOAT)0.0) {
		sgn = -(FLOAT)1.0; y = -y;
	}
	else {
		sgn = (FLOAT)1.0;
	}

	t = y / (FLOAT)3.75;

	if (y <= FLOAT_MIN) {
		I1 = (FLOAT)0.0;
	}
	else
	if (y <= (FLOAT)3.75) {
		I1 = (FLOAT)0.5 +
			 (FLOAT)0.87890594 * (FLOAT)pow(t, 2) +
			 (FLOAT)0.51498869 * (FLOAT)pow(t, 4) +
			 (FLOAT)0.15084934 * (FLOAT)pow(t, 6) +
			 (FLOAT)0.02658733 * (FLOAT)pow(t, 8) +
			 (FLOAT)0.00301532 * (FLOAT)pow(t, 10) +
			 (FLOAT)0.00032411 * (FLOAT)pow(t, 12);
		
		I1 *= sgn * y;
	}
	else {
		I1 = (FLOAT)0.39894228 -
			 (FLOAT)0.03988024 * (FLOAT)pow(t, -1) -
			 (FLOAT)0.00362018 * (FLOAT)pow(t, -2) +
			 (FLOAT)0.00163801 * (FLOAT)pow(t, -3) -
			 (FLOAT)0.01031555 * (FLOAT)pow(t, -4) +
			 (FLOAT)0.02282967 * (FLOAT)pow(t, -5) -
			 (FLOAT)0.02895312 * (FLOAT)pow(t, -6) +
			 (FLOAT)0.01787654 * (FLOAT)pow(t, -7) -
			 (FLOAT)0.00420059 * (FLOAT)pow(t, -8);

		I1 *= sgn * (FLOAT)exp(y) / (FLOAT)sqrt(y);
	}

	return (I1);
} // BesselI1

// Returns the von Mises c.d.f. for the specified Mean and Kappa.

FLOAT vonMisesCdf(FLOAT y, FLOAT Mean, FLOAT Kappa)
{
	FLOAT A[3], Io, In, Fy;
	int   i, Error;

	if (y > Pi2) {
		Fy = (FLOAT)1.0;
	}
	else
	if (y < (FLOAT)0.0) {
		Fy = (FLOAT)0.0;
	}
	else {
		Io = BesselI0(Kappa); In = BesselI1(Kappa);

		A[0] = (FLOAT)1.0 / Pi2; A[1] = (FLOAT)2.0 * A[0] / Io;

		i = 1; Fy = A[0] * y; Error = 1;
		while ((i <= ItMax) && Error) {
			Fy += A[1] * In * ((FLOAT)sin(i * (y - Mean)) + (FLOAT)sin(i * Mean)) / i;

			if (In < Eps) Error = 0;

			A[2] = Io - (FLOAT)2.0 * i * In / Kappa; Io = In; In = A[2];

			i++;
		}

		if (Fy >(FLOAT)1.0) {
			Fy = (FLOAT)1.0;
		}
		else
		if (Fy < (FLOAT)0.0) {
			Fy = (FLOAT)0.0;
		}
	}

	return (Fy);
} // vonMisesCdf

// Returns the inverse of the von Mises c.d.f. for the specified Mean and Kappa.

FLOAT vonMisesInv(FLOAT Fy, FLOAT Mean, FLOAT Kappa)
{
	FLOAT yl, ym, yh, fl, fm;
	int   Stop;

	if (Fy >= (FLOAT)1.0) {
		ym = Pi2;
	}
	else
	if (Fy <= (FLOAT)0.0) {
		ym = (FLOAT)0.0;
	}
	else {
		yl = (FLOAT)0.0; fl = Fy - vonMisesCdf(yl, Mean, Kappa);
		yh = Pi2;

		// Bisection.

		Stop = 0;

		while (!Stop) {
			ym = (yh + yl) / (FLOAT)2.0;

			fm = Fy - vonMisesCdf(ym, Mean, Kappa);

			if (((FLOAT)fabs(fm) < Eps) || (yh - yl < Eps)) {
				Stop = 1;
			}
			else {
				if (fm * fl >(FLOAT)0.0) {
					yl = ym; fl = fm;
				}
				else {
					yh = ym;
				}
			}
		}
	}

	return (ym);
} // vonMisesInv

FLOAT xlogx(FLOAT x)
{
    if (x > FLOAT_MIN) {
        x *= (FLOAT)log(x);
    }
    else {
        x = (FLOAT)0.0;
    }

    return(x);
} // xlogx