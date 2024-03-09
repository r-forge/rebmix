#include "base.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

char _e_line_[65536] = "";
char _w_line_[2][65536] = { "", "" };

void E_begin()
{
    strcpy(_e_line_, ""); strcpy(_w_line_[0], ""); strcpy(_w_line_[1], "");
} // E_begin

void Print_e_line_(const char *file, INT line, INT error)
{
    sprintf(_e_line_, "File = %s; line = %d; error = %d.", file, line, error);
} // Print_e_line_

void Print_w_line_(INT idx)
{
    if (!strcmp(_w_line_[idx], "")) strcpy(_w_line_[idx], _e_line_);
    
    strcpy(_e_line_, "");
} // Print_w_line_

void Print_e_list_(char *elist)
{
    sprintf(elist, "%s\n%s\n%s\n", _e_line_, _w_line_[0], _w_line_[1]);
} // Print_e_list_

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

static INT IY = 0;
static INT IV[NTAB];

// Minimal random number generator of Park and Miller with Bays-Durham shuffle and added
// safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
// values). Call with IDum a negative integer to initialize; thereafter do not alter IDum
// between successive deviates in a sequence.RNMX should approximate the largest floating
// value that is less than 1. See http://www.nrbook.com/a/bookcpdf/c7-1.pdf

FLOAT Ran1(INT *IDum)
{
    FLOAT Tmp;
    INT   j, k;

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

    if ((Tmp = AM * IY) > RNMX) return RNMX; else return Tmp;
} // Ran1

// Inserts y into ascending list Y of length n.Set n = 0 initially.

void Insert(FLOAT y,   // Inserted value.
            INT   *n,  // Length of Y.
            FLOAT *Y)  // Pointer to Y = [y0,...,yn-1].
{
    INT i, j;

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
    static FLOAT Cof[6] = { (FLOAT)76.18009172947146, -(FLOAT)86.50532032941677,
        (FLOAT)24.01409824083091, -(FLOAT)1.231739572450155,
        (FLOAT)0.1208650973866179E-2, -(FLOAT)0.5395239384953E-5 };
    static FLOAT Stp = (FLOAT)2.5066282746310005;

    FLOAT Ser, Tmp, x, z;
    INT   j;

    z = x = y; Tmp = x + (FLOAT)5.5; Tmp -= (x + (FLOAT)0.5) * (FLOAT)log(Tmp);

    Ser = (FLOAT)1.000000000190015;

    for (j = 0; j < 6; j++) Ser += Cof[j] / ++z;

    return (-Tmp + (FLOAT)log(Stp * Ser / x));
} // Gammaln

// Returns the digamma for y > 0. See http://www.nr.com/.

INT Digamma(FLOAT y, FLOAT *Psi)
{
    static FLOAT piov4 = (FLOAT)0.785398163397448;
    static FLOAT dy0 = (FLOAT)1.461632144968362341262659542325721325;
    static FLOAT p1[7] = {(FLOAT)0.0089538502298197, (FLOAT)4.77762828042627, (FLOAT)142.441585084029,
        (FLOAT)1186.45200713425, (FLOAT)3633.51846806499, (FLOAT)4138.10161269013,
        (FLOAT)1305.60269827897};
    static FLOAT q1[6] = {(FLOAT)44.8452573429826, (FLOAT)520.752771467162, (FLOAT)2210.0079924783,
        (FLOAT)3641.27349079381, (FLOAT)1908.310765963, (FLOAT)6.91091682714533e-6 };
    static FLOAT p2[4] = {-(FLOAT)2.12940445131011, -(FLOAT)7.01677227766759,
        -(FLOAT)4.48616543918019, -(FLOAT)0.648157123766197};
    static FLOAT q2[4] = {(FLOAT)32.2703493791143, (FLOAT)89.2920700481861,
        (FLOAT)54.6117738103215, (FLOAT)7.77788548522962};

    FLOAT aug, d2, den, sgn, upper, w, ymax, ymin, ymy0, z;
    INT   i, m, n, nq, Error = E_OK;

    ymax = (FLOAT)INT_MAX; d2 = (FLOAT)1.0 / FLOAT_EPSILON; if (ymax > d2) ymax = d2; ymin = (FLOAT)1E-9; aug = (FLOAT)0.0;

    if (y < (FLOAT)0.5) {
        if ((FLOAT)fabs(y) <= ymin) {
            E_CHECK(y == (FLOAT)0.0, E_ARG);

            aug = -(FLOAT)1.0 / y;
        }
        else {
            w = -y; sgn = piov4;

            if (w <= (FLOAT)0.0) {
                w = -w; sgn = -sgn;
            }

            E_CHECK(w >= ymax, E_ARG);

            nq = (INT)w; w -= (FLOAT)nq; nq = (INT)(w * (FLOAT)4.0); w = (w - (FLOAT)nq * (FLOAT)0.25) * (FLOAT)4.0; n = nq / 2;

            if (n + n != nq) {
                w = (FLOAT)1.0 - w;
            }

            z = piov4 * w; m = n / 2;

            if (m + m != n) {
                sgn = -sgn;
            }

            n = (nq + 1) / 2; m = n / 2; m += m;

            if (m == n) {
                E_CHECK(z == (FLOAT)0.0, E_ARG);

                aug = sgn * ((FLOAT)cos(z) / (FLOAT)sin(z) * (FLOAT)4.0);
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

EEXIT: 
    
    E_RETURN(Error);
} // Digamma

// Returns binomial c.d.f. for the specified n and p. See http://www.nr.com/.

FLOAT BinomialCdf(INT k, INT n, FLOAT p)
{
    FLOAT Fy, ypb;
    INT   y;

    if (k < 0)
        Fy = (FLOAT)0.0;
    else
    if (k == 0)
        Fy = (FLOAT)pow((FLOAT)1.0 - p, n);
    else
    if (k == n)
        Fy = (FLOAT)1.0;
    else
    if (k > n)
        Fy = (FLOAT)0.0;
    else {
        Fy = ypb = (FLOAT)pow((FLOAT)1.0 - p, n); y = 0;

        while ((y < k) && (ypb > FLOAT_MIN)) {
            y++; ypb *= (n - y + (FLOAT)1.0) * p / y / ((FLOAT)1.0 - p); Fy += ypb;
        }
    }

    return Fy;
} // BinomialCdf

// Returns the inverse of the binomial c.d.f. for the specified n and p. See http://www.nr.com/.

INT BinomialInv(FLOAT Fy, INT n, FLOAT p)
{
    FLOAT Sum, ypb;
    INT   y;

    Sum = ypb = (FLOAT)pow((FLOAT)1.0 - p, n); y = 0;

    while ((Sum < Fy) && (ypb > FLOAT_MIN)) {
        y++; ypb *= (n - y + (FLOAT)1.0) * p / y / ((FLOAT)1.0 - p); Sum += ypb;
    }

    if ((Fy < (FLOAT)0.5) && (y > 0)) y--;

    return y;
} // BinomialInv

// Returns the Poisson c.d.f. for the specified Theta.

FLOAT PoissonCdf(INT k, FLOAT Theta)
{
    FLOAT Fy, ypb;
    INT   y;

    if (k < 0)
        Fy = (FLOAT)0.0;
    else
    if (k == 0)
        Fy = (FLOAT)exp(-Theta);
    else {
        Fy = ypb = (FLOAT)exp(-Theta); y = 0;

        while ((y < k) && (ypb > FLOAT_MIN)) {
            y++; ypb *= Theta / y; Fy += ypb;
        }
    }

    return Fy;
} // PoissonCdf

// Returns the inverse of the Poisson c.d.f. for the specified Theta.

INT PoissonInv(FLOAT Fy, FLOAT Theta)
{
    FLOAT Sum, ypb;
    INT   y;

    Sum = ypb = (FLOAT)exp(-Theta); y = 0;

    while ((Sum < Fy) && (ypb > FLOAT_MIN)) {
        y++; ypb *= Theta / y; Sum += ypb;
    }

    if ((Fy < (FLOAT)0.5) && (y > 0)) y--;

    return y;
} // PoissonInv

// Returns the incomplete gamma function P(a, y) evaluated by its series
// representation as GamSer. Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

INT GammaSer(FLOAT a,       // Constant a > 0.
             FLOAT y,       // Variable y > 0.
             FLOAT *GamSer, // Incomplete gamma function.
             FLOAT *Gamln)  // Log(Gamma(a)).
{
    FLOAT ap, Del, Sum;
    INT   i, Error = E_OK;

    *Gamln = Gammaln(a);

    if (y <= FLOAT_MIN) {
        *GamSer = (FLOAT)0.0;
    }
    else {
        i = 1; Error = E_CON; ap = a; Sum = (FLOAT)1.0 / a; Del = Sum;

        while ((i <= ItMax) && (Error != E_OK)) {
            ap += (FLOAT)1.0; Del *= y / ap; Sum += Del;

            if ((FLOAT)fabs(Del) < Eps) Error = E_OK;

            i++;
        }

        if (Error != E_OK) Error = E_OK; // ItMax too small.

        *GamSer = Sum * (FLOAT)exp(-y + a * (FLOAT)log(y) - *Gamln);
    }

    E_RETURN(Error);
} // GammaSer

// Returns the incomplete gamma function Q(a, y) evaluated by its continued
// fraction representation as GamCfg. Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

INT GammaCfg(FLOAT a,       // Constant a > 0.
             FLOAT y,       // Variable y > 0.
             FLOAT *GamCfg, // Incomplete gamma function.
             FLOAT *Gamln)  // Log(Gamma(a)).
{
    FLOAT a0, a1, ai, aia, aif, b0, b1, Fac, G, Gold;
    INT   i, Error = E_OK;

    *Gamln = Gammaln(a);

    if (y <= FLOAT_MIN) {
        *GamCfg = (FLOAT)0.0;
    }
    else {
        G = (FLOAT)0.0; Gold = (FLOAT)0.0; Fac = (FLOAT)1.0;

        i = 1; Error = E_CON; a0 = (FLOAT)1.0; a1 = y; b0 = (FLOAT)0.0; b1 = (FLOAT)1.0;

        while ((i <= ItMax) && (Error != E_OK)) {
            ai = (FLOAT)1.0 * i; aia = ai - a; aif = ai * Fac;

            a0 = (a1 + a0 * aia) * Fac;
            b0 = (b1 + b0 * aia) * Fac;

            a1 = y * a0 + aif * a1;
            b1 = y * b0 + aif * b1;

            if (a1 != (FLOAT)0.0) {
                Fac = (FLOAT)1.0 / a1; G = b1 * Fac;

                if ((FLOAT)fabs(G - Gold) < Eps) Error = E_OK; else Gold = G;
            }

            i++;
        }

        if (Error != E_OK) Error = E_OK; // ItMax too small.

        *GamCfg = (FLOAT)exp(-y + a * (FLOAT)log(y) - *Gamln) * G;
    }

    E_RETURN(Error);
} // GammaCfg

// Returns the incomplete gamma function P(a, y). Also returns log(Gamma(a)) as Gamln. See http://www.nr.com/.

INT GammaP(FLOAT a,       // Constant a > 0.
           FLOAT y,       // Variable y > 0.
           FLOAT *GamP,   // Incomplete gamma function.
           FLOAT *Gamln)  // Log(Gamma(a)).
{
    FLOAT GamCfg, GamSer;
    INT   Error = E_OK;

    if ((y <= FLOAT_MIN) || (a <= FLOAT_MIN)) {
        *GamP = (FLOAT)0.0;
    }
    else
    if (y < a + (FLOAT)1.0) {
        Error = GammaSer(a, y, &GamSer, Gamln);

        E_CHECK(Error != E_OK, Error);

        *GamP = GamSer;
    }
    else {
        Error = GammaCfg(a, y, &GamCfg, Gamln);

        E_CHECK(Error != E_OK, Error);

        *GamP = (FLOAT)1.0 - GamCfg;
    }

EEXIT:

    E_RETURN(Error);
} // GammaP

// Returns the inverse of the gamma c.d.f. for the specified Theta and Beta. See http://www.nr.com/.

INT GammaInv(FLOAT Fy, FLOAT Theta, FLOAT Beta, FLOAT *y)
{
    FLOAT dx, dy, Gamln, GamP, Tmp;
    INT   i, error, Error = E_OK;

    if (Beta > (FLOAT)1.0) {
        *y = (Beta - (FLOAT)1.0) * Theta + Eps;
    }
    else
        *y = Eps;

    i = 1; Error = E_CON; dx = (FLOAT)0.0; 

    while ((i <= ItMax) && (Error != E_OK)) {
        error = GammaP(Beta, *y / Theta, &GamP, &Gamln);

        E_CHECK(error != E_OK, error);

        Tmp = *y / Theta;

        dy = (GamP - Fy) / ((FLOAT)exp(Beta * (FLOAT)log(Tmp) - Tmp - Gamln) / (*y));

        *y -= dy;

        E_CHECK(IsNan(dy) || IsInf(dy), E_CON);

        if (*y < Eps) {
            *y = Eps; Error = E_OK;
        }

        if (((FLOAT)fabs(dy) < Eps) || ((FLOAT)fabs(dx + dy) < Eps)) Error = E_OK; 
        
        dx = dy;

        i++;
    }

EEXIT:

    E_RETURN(Error);
} // GammaInv

// Returns the inverse of the Weibull c.d.f. for the specified Theta and Beta.

FLOAT WeibullInv(FLOAT Fy, FLOAT Theta, FLOAT Beta)
{
    FLOAT y;

    y = Theta * (FLOAT)pow(-(FLOAT)log((FLOAT)1.0 - Fy), (FLOAT)1.0 / Beta);

    return y;
} // WeibullInv

// Returns the inverse of the Gumbel c.d.f. for the specified Mean, Sigma and Xi.

FLOAT GumbelInv(FLOAT Fy, FLOAT Mean, FLOAT Sigma, FLOAT Xi)
{
    FLOAT y;

    if (Xi > Eps) {
        y = Mean + Sigma * (FLOAT)log((FLOAT)log((FLOAT)1.0 / ((FLOAT)1.0 - Fy)));
    }
    else {
        y = Mean - Sigma * (FLOAT)log((FLOAT)log((FLOAT)1.0 / Fy));
    }

    return y;
} // GumbelInv

// Returns the error function erf(y). See http://www.nr.com/.

INT ErrorF(FLOAT y,     // Variable y.
           FLOAT *ErF)  // Error function.
{
    FLOAT Gamln, GamP;
    INT   Error = E_OK;

    Error = GammaP((FLOAT)0.5, y * y, &GamP, &Gamln);

    E_CHECK(Error != E_OK, Error);

    if (y < (FLOAT)0.0)
        *ErF = -GamP;
    else
        *ErF = GamP;

EEXIT:

    E_RETURN(Error);
} // ErrorF

// Returns the LU decomposition of matrix A. See http://www.nr.com/

INT LUdcmp(INT   n,     // Size of square matrix.
           FLOAT *A,    // Pointer to the square matrix A.
           INT   *indx, // Pointer to the permutation vector.
           FLOAT *det)  // Determinant.
{
    FLOAT Big, Tmp, *V = NULL;
    INT   i, imax, j, k, Error = E_OK;

    V = (FLOAT*)malloc(n * sizeof(FLOAT));

    E_CHECK(NULL == V, E_MEM);

    for (i = 0; i < n; i++) {
        Big = (FLOAT)0.0;

        for (j = 0; j < n; j++) {
            if ((Tmp = (FLOAT)fabs(A[i * n + j])) > Big) Big = Tmp;
        }

        E_CHECK((FLOAT)fabs(Big) <= FLOAT_MIN, E_CON);

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

    E_CHECK(IsNan(*det) || ((FLOAT)fabs(*det) <= FLOAT_MIN), E_CON);

EEXIT:

    if (V) free(V);

    E_RETURN(Error);
} // LUdcmp

// Solves the set of n linear equations A x = b. See See http://www.nr.com/

INT LUbksb(INT   n,     // Size of square matrix.
           FLOAT *A,    // Pointer to the square matrix A.
           INT   *indx, // Pointer to the permutation vector.
           FLOAT *b)    // Pointer to the solution vector.
{
    FLOAT Sum;
    INT   i, ii = 0, ip, j, Error = E_OK;

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

    E_RETURN(Error);
} // LUbksb

// Returns the determinant and the inverse matrix of A. See http://www.nr.com/

INT LUinvdet(INT   n,     // Size of square matrix.
             FLOAT *A,    // Pointer to the square matrix A.
             FLOAT *Ainv, // Pointer to the inverse matrix of A.
             FLOAT *Adet) // Pointer to the determinant of A.
{
    FLOAT *b = NULL, *B = NULL;
    INT   i, *indx = NULL, j, Error = E_OK;

    indx = (INT*)calloc((size_t)n, sizeof(INT));

    E_CHECK(NULL == indx, E_MEM);

    b = (FLOAT*)malloc(n * sizeof(FLOAT));

    E_CHECK(NULL == b, E_MEM);

    B = (FLOAT*)malloc(n * n * sizeof(FLOAT));

    E_CHECK(NULL == B, E_MEM);

    memmove(B, A, n * n * sizeof(FLOAT));

    Error = LUdcmp(n, B, indx, Adet);

    E_CHECK(Error != E_OK, Error);

    for (j = 0; j < n; j++) {
        memset(b, 0, n * sizeof(FLOAT));

        b[j] = (FLOAT)1.0;

        Error = LUbksb(n, B, indx, b);

        E_CHECK(Error != E_OK, Error);

        for (i = 0; i < n; i++) Ainv[i * n + j] = b[i];
    }

EEXIT:

    if (B) free(B);

    if (b) free(b);

    if (indx) free(indx);

    E_RETURN(Error);
} // LUinvdet

// Returns the Cholesky decomposition of matrix A. See http://www.nr.com/

INT Choldc(INT   n,   // Size of square matrix.
           FLOAT *A,  // Pointer to the square matrix A.
           FLOAT *L)  // Lower triangular factors.
{
    FLOAT *p = NULL, Sum;
    INT   i, j, k, Error = E_OK;

    memmove(L, A, n * n * sizeof(FLOAT));

    p = (FLOAT*)malloc(n * sizeof(FLOAT));

    E_CHECK(NULL == p, E_MEM);

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

EEXIT:

    if (p) free(p);

    E_RETURN(Error);
} // Choldc

// Returns the determinant and the inverse matrix of A. See http://www.nr.com/.

INT Cholinvdet(INT   n,         // Size of square matrix.
               FLOAT *A,        // Pointer to the symmetric square matrix A.
               FLOAT *Ainv,     // Pointer to the inverse matrix of A.
               FLOAT *logAdet)  // Pointer to the logarithm of determinant of A.
{
    FLOAT *L = NULL, *p = NULL, Sum;
    INT   i, j, k, Error = E_OK;

    L = (FLOAT*)malloc(n * n * sizeof(FLOAT));

    E_CHECK(NULL == L, E_MEM);

    memmove(L, A, n * n * sizeof(FLOAT));

    p = (FLOAT*)malloc(n * sizeof(FLOAT));

    E_CHECK(NULL == p, E_MEM);

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

EEXIT:

    if (p) free(p);

    if (L) free(L);

    E_RETURN(Error);
} // Cholinvdet

// Returns modified Bessel function of order 0. See http://www.nr.com/.

FLOAT BesselI0(FLOAT x)
{
    FLOAT I0, y;

    if (x < (FLOAT)3.75) {
        y = x / (FLOAT)3.75; y *= y;
        
        I0 = (FLOAT)1.0 + y * ((FLOAT)3.5156229 + y * ((FLOAT)3.0899424 + y * ((FLOAT)1.2067492
            + y * ((FLOAT)0.2659732 + y * ((FLOAT)0.360768e-1 + y * (FLOAT)0.45813e-2)))));
    }
    else {
        y = 3.75 / x;

        I0 = ((FLOAT)exp(x) / (FLOAT)sqrt(x)) * ((FLOAT)0.39894228 + y * ((FLOAT)0.1328592e-1
            + y * ((FLOAT)0.225319e-2 + y * (-(FLOAT)0.157565e-2 + y * ((FLOAT)0.916281e-2
            + y * (-(FLOAT)0.2057706e-1 + y * ((FLOAT)0.2635537e-1 + y * (-(FLOAT)0.1647633e-1
            + y * (FLOAT)0.392377e-2))))))));
    }

    return I0;
} // BesselI0

// Returns modified Bessel function of order 1. See http://www.nr.com/.

FLOAT BesselI1(FLOAT x)
{
    FLOAT I1, y;
    INT   sgn = 1;

    if (x < (FLOAT)0.0) {
        sgn = -1; x = -x;
    }
    else {
        sgn = 1;
    }

    if (x < (FLOAT)3.75) {
        y = x / (FLOAT)3.75; y *= y;

        I1 = x * ((FLOAT)0.5 + y * ((FLOAT)0.87890594 + y * ((FLOAT)0.51498869 + y * ((FLOAT)0.15084934
            + y * ((FLOAT)0.2658733e-1 + y * ((FLOAT)0.301532e-2 + y * (FLOAT)0.32411e-3))))));
    }
    else {
        y = (FLOAT)3.75 / x;

        I1 = (FLOAT)0.2282967e-1 + y * (-(FLOAT)0.2895312e-1 + y * ((FLOAT)0.1787654e-1
            - y * (FLOAT)0.420059e-2));

        I1 = (FLOAT)0.39894228 + y * (-(FLOAT)0.3988024e-1 + y * (-(FLOAT)0.362018e-2
            + y * ((FLOAT)0.163801e-2 + y * (-(FLOAT)0.1031555e-1 + y * I1))));

        I1 *= ((FLOAT)exp(I1) / (FLOAT)sqrt(x));
    }

    return sgn < 0 ? -I1 : I1;
} // BesselI1

INT vonMisesCdf(FLOAT y, FLOAT Mean, FLOAT Kappa, FLOAT *Fy)
{
    FLOAT A[3], Io, In;
    INT   i, Error = E_OK;

    if (y > Pi2) {
        *Fy = (FLOAT)1.0;
    }
    else
    if (y < (FLOAT)0.0) {
        *Fy = (FLOAT)0.0;
    }
    else {
        Io = BesselI0(Kappa); In = BesselI1(Kappa);

        A[0] = (FLOAT)1.0 / Pi2; A[1] = (FLOAT)2.0 * A[0] / Io;

        i = 1; Error = E_CON; *Fy = A[0] * y;

        while ((i <= ItMax) && (Error != E_OK)) {
            *Fy += A[1] * In * ((FLOAT)sin(i * (y - Mean)) + (FLOAT)sin(i * Mean)) / i;

            A[2] = Io - (FLOAT)2.0 * i * In / Kappa; Io = In; In = A[2];

            if (In < Eps) Error = E_OK;

            i++;
        }

        if (*Fy > (FLOAT)1.0) {
             *Fy = (FLOAT)1.0;
        }
        else
        if (*Fy < (FLOAT)0.0) {
            *Fy = (FLOAT)0.0;
        }
    }

    E_RETURN(Error);
} // vonMisesCdf

// Returns the inverse of the von Mises c.d.f. for the specified Mean and Kappa.

INT vonMisesInv(FLOAT Fy, FLOAT Mean, FLOAT Kappa, FLOAT *y)
{
    FLOAT fl, fm, Fyt, yl, yh;
    INT   i, Error = E_OK;

    if (Fy >= (FLOAT)1.0) {
        *y = Pi2;
    }
    else
    if (Fy <= (FLOAT)0.0) {
        *y = (FLOAT)0.0;
    }
    else {
        yl = (FLOAT)0.0;

        Error = vonMisesCdf(yl, Mean, Kappa, &Fyt);

        E_CHECK(Error != E_OK, Error);

        fl = Fy - Fyt; yh = Pi2;

        // Bisection.

        i = 1; Error = E_CON;

        while ((i <= ItMax) && (Error != E_OK)) {
            *y = (yh + yl) / (FLOAT)2.0;

            Error = vonMisesCdf(*y, Mean, Kappa, &Fyt);

            E_CHECK(Error != E_OK, Error);

            fm = Fy - Fyt;

            if (((FLOAT)fabs(fm) < Eps) || (yh - yl < Eps)) {
                Error = E_OK;
            }
            else {
                if (fm * fl > (FLOAT)0.0) {
                    yl = *y; fl = fm;
                }
                else {
                    yh = *y;
                }
            }

            i++;
        }
    }

EEXIT:

    E_RETURN(Error);
} // vonMisesInv

FLOAT xlogx(FLOAT x)
{
    if (x > FLOAT_MIN) {
        x *= (FLOAT)log(x);
    }
    else {
        x = (FLOAT)0.0;
    }

    return x;
} // xlogx

// Returns merged intervals.

void MergeIntervals(FLOAT    ym,  // Mode position. 
                    INT      *n,  // Total number of intervals.
                    Interval *X)  // Pointer to the intervals.
{
    Interval Tmp;
    INT      i, j, k;

    if (*n <= 1) return;

    // Bubble sort.

    for (i = 0; i < *n - 1; i++) {
        for (j = 0; j < *n - i - 1; j++) {
            if (X[j].a > X[j + 1].a) {
                Tmp = X[j]; X[j] = X[j + 1]; X[j + 1] = Tmp;
            }
        }
    }

    // Merge intervals.

    k = 0;

    for (i = 1; i < *n; i++) {
        if (X[i].a <= X[k].b) {
            if (X[i].b > X[k].b) X[k].b = X[i].b;
        }
        else {
            k++; X[k] = X[i];
        }
    }
    
    *n = ++k;

    // Side of interval.

    for (i = 0; i < *n; i++) {
        if (X[i].b <= ym) {
            X[i].s = 0;
        }
        else
        if (X[i].a >= ym) {
            X[i].s = 1;
        }
        else {
            X[k].b = X[i].b; X[i].b = X[k].a = ym; X[i].s = 0; X[k].s = 1;
            
            k++;
        }
    }

    *n = k;
} // MergeIntervals
