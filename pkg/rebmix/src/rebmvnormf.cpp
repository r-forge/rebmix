#include "rebmvnormf.h"

#include <math.h>

// Perform necessary initializations.

INT Rebmvnorm::Initialize()
{
    INT Error = E_OK;

    p_value_ = (FLOAT)0.0001;

    min_dist_mul_ = (FLOAT)2.5;

    var_mul_ = (FLOAT)0.0625;

    kmax_ = (INT)floor(((FLOAT)1.0 + (FLOAT)1.0 / length_pdf_) * (FLOAT)pow((FLOAT)n_, (FLOAT)1.0 / ((FLOAT)1.0 + (FLOAT)1.0 / length_pdf_)));

    Error = GammaInv((FLOAT)1.0 - (FLOAT)2.0 * p_value_, (FLOAT)2.0, length_pdf_ / (FLOAT)2.0, &ChiSqr_);

    E_CHECK(Error != E_OK, Error);

EEXIT:

    E_RETURN(Error);
} // Initialize

INT Rebmvnorm::ComponentConditionalDist(INT                  i,           // Index of variable y.
                                        FLOAT                Y,           // Variable Y[i].
                                        FLOAT                *Cinv,       // Inverse correlation matrix.
                                        CompnentDistribution *CmpTheta,   // Component distribution type.
                                        FLOAT                *CmpConDist) // Component conditional distribution.
{
    FLOAT Mean, Stdev, y;
    INT   o, Error = E_OK;

    Mean = CmpTheta->Theta_[0][i];

    o = i * length_pdf_ + i;

    Stdev = (FLOAT)sqrt(CmpTheta->Theta_[1][o] / Cinv[o]);

    y = (Y - Mean) / (Sqrt2 * Stdev); y *= y;

    *CmpConDist = (FLOAT)exp(-y) / (SqrtPi2 * Stdev);

    E_RETURN(Error);
} // ComponentConditionalDist

INT Rebmvnorm::ComponentConditionalCdf(INT                  i,           // Index of variable y.
                                       FLOAT                Y,           // Variable Y[i].
                                       FLOAT                *Cinv,       // Inverse correlation matrix.
                                       CompnentDistribution *CmpTheta,   // Component distribution type.
                                       FLOAT                *CmpConCdf)  // Component conditional cumulative distribution.
{
    FLOAT Mean, Stdev, y;
    INT   o, Error = E_OK;

    Mean = CmpTheta->Theta_[0][i];

    o = i * length_pdf_ + i;

    Stdev = (FLOAT)sqrt(CmpTheta->Theta_[1][o] / Cinv[o]);

    y = (Y - Mean) / (Sqrt2 * Stdev);

    Error = ErrorF(y, CmpConCdf);

    E_CHECK(Error != E_OK, Error);

    *CmpConCdf = (FLOAT)0.5 * ((FLOAT)1.0 + *CmpConCdf);

EEXIT:

    E_RETURN(Error);
} // ComponentConditionalCdf

// Rough component parameter estimation for k-nearest neighbours.

INT Rebmvnorm::RoughEstimationKNN(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,logV,R].
                                  INT                  k,           // k-nearest neighbours.
                                  FLOAT                *h,          // Normalizing vector.
                                  FLOAT                nl,          // Total number of observations in class l.
                                  INT                  m,           // Mode index.
                                  CompnentDistribution *RigidTheta, // Rigid parameters.
                                  CompnentDistribution *LooseTheta) // Loose parameters.
{
    RoughParameterType *Mode = NULL;
    Interval           *X = NULL;
    FLOAT              *C = NULL, *Cinv = NULL, CmpConCdf[2], *D = NULL, Dll, Dlm, Dlmin, Dlr, epsilon, flm, flmax, flmin, logCdet, logflm, Mean, Stdev, Sum, Tmp, ym;
    INT                i, I, ii, j, jj, l, *N = NULL, o, p, q, r, s, error, Error = E_OK;

    // Global mode.

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    E_CHECK(NULL == Mode, E_MEM);

    N = (INT*)malloc(length_pdf_ * sizeof(INT));

    E_CHECK(NULL == N, E_MEM);

    D = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == D, E_MEM);

    for (i = 0; i < length_pdf_; i++) {
        N[i] = 0; D[i] = (FLOAT)2.0 * Y[length_pdf_ + 2][m] * h[i];

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                for (l = 0; l < length_pdf_; l++) if ((i != l) && ((FLOAT)fabs(Y[l][j] - Y[l][m]) > (FLOAT)0.5 * D[l])) goto S0;

                Mode[i].klm += Y[length_pdf_][j]; X_[i][N[i]] = (FLOAT)j; N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                X_[i][N[i]] = (FLOAT)j; N[i] += 1;
            }
        }

        Mode[i].ym = Y[i][m]; Mode[i].flm = Y[length_pdf_][m] * k / (Mode[i].klm * D[i]);
    }

    logflm = (FLOAT)log(Y[length_pdf_][m] * k / nl) - Y[length_pdf_ + 1][m];

    // Variance-covariance matrix.

    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == C, E_MEM);

    Cinv = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == Cinv, E_MEM);

    if (nl > length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            Sum = (FLOAT)0.0;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                Sum += Y[length_pdf_][j] * (Y[i][j] - Mode[i].ym) * (Y[i][j] - Mode[i].ym);
            }

            C[i * length_pdf_ + i] = Sum / nl;

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;

                for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                    Sum += Y[length_pdf_][j] * (Y[i][j] - Mode[i].ym) * (Y[ii][j] - Mode[ii].ym);
                }

                C[i * length_pdf_ + ii] = C[ii * length_pdf_ + i] = Sum / nl;
            }
        }

        // Correlation matrix.

        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i;

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

                C[p] = C[q] /= (FLOAT)sqrt(C[o] * C[r]);

                if (IsNan(C[q]) || IsInf(C[q])) {
                    C[p] = C[q] = (FLOAT)0.0;
                }
            }
        }

        for (i = 0; i < length_pdf_; i++) {
            C[i * length_pdf_ + i] = (FLOAT)1.0;
        }

        Error = Cholinvdet(length_pdf_, C, Cinv, &logCdet);

        E_CHECK(Error != E_OK, Error);
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0;

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
            }
        }

        logCdet = (FLOAT)0.0;
    }

    // Rigid restraints.

    RigidTheta->Theta_[3][0] = logCdet;

    for (i = 0; i < length_pdf_; i++) {
        RigidTheta->Theta_[0][i] = Mode[i].ym;

        Stdev = (FLOAT)1.0 / (SqrtPi2 * Mode[i].flm); Stdev *= Stdev;

        o = i * length_pdf_ + i;

        RigidTheta->Theta_[1][o] = Cinv[o] * Stdev; RigidTheta->Theta_[2][o] = (FLOAT)1.0 / Stdev;

        RigidTheta->Theta_[3][0] += (FLOAT)log(RigidTheta->Theta_[1][o]);

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(RigidTheta->Theta_[1][o] * RigidTheta->Theta_[1][r]);

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] = C[q] * Stdev;

            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

    epsilon = (FLOAT)exp(-(FLOAT)2.0 * (LogSqrtPi2 + logflm / length_pdf_) - RigidTheta->Theta_[3][0] / length_pdf_);

    RigidTheta->Theta_[3][0] += length_pdf_ * (FLOAT)log(epsilon);

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        RigidTheta->Theta_[1][o] *= epsilon;
        RigidTheta->Theta_[2][o] /= epsilon;

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] *= epsilon;
            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] /= epsilon;
        }
    }

    Error = LooseTheta->Memmove(RigidTheta);

    E_CHECK(Error != E_OK, Error);

    if ((Restraints_ == rtRigid) || (nl <= length_pdf_)) goto EEXIT;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        X = (Interval*)malloc((N[i] + 1) * sizeof(Interval));

        E_CHECK(NULL == X, E_MEM);

        // Bracketing.

        Tmp = ((FLOAT)1.0 + Eps) * h[i];

        for (j = 0; j < N[i]; j++) {
            l = (INT)X_[i][j];

            X[j].a = Y[i][l] - Y[length_pdf_ + 2][l] * Tmp;
            X[j].b = Y[i][l] + Y[length_pdf_ + 2][l] * Tmp;
        }

        I = N[i]; MergeIntervals(Mode[i].ym, &I, X);

        Error = ComponentConditionalCdf(i, Mode[i].ym, Cinv, LooseTheta, &Tmp);

        E_CHECK(Error != E_OK, Error);

        Dll = Tmp; Dlr = (FLOAT)1.0 - Tmp;

        for (j = 0; j < I; j++) {
            Error = ComponentConditionalCdf(i, X[j].a, Cinv, LooseTheta, &CmpConCdf[0]);

            E_CHECK(Error != E_OK, Error);

            Error = ComponentConditionalCdf(i, X[j].b, Cinv, LooseTheta, &CmpConCdf[1]);

            E_CHECK(Error != E_OK, Error);

            Tmp = CmpConCdf[1] - CmpConCdf[0];

            if (X[j].s == 1) Dlr -= Tmp; else Dll -= Tmp;
        }

        Dlm = Dll + Dlr;

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0; flmax = Mode[i].flm;

        // Bisection.

        ii = 1; Error = E_CON;

        while ((ii <= ItMax) && (Error != E_OK)) {
            ym = Mode[i].ym; flm = (flmax + flmin) / (FLOAT)2.0;

            jj = 1; s = 0;

            while (jj <= ItMax) {
                Mean = ym;

                LooseTheta->Theta_[0][i] = Mean;

                Stdev = (FLOAT)1.0 / (SqrtPi2 * flm); Stdev *= Stdev;

                o = i * length_pdf_ + i;

                LooseTheta->Theta_[1][o] = Cinv[o] * Stdev;

                error = ComponentConditionalCdf(i, Mode[i].ym, Cinv, LooseTheta, &Tmp);

                E_CHECK(error != E_OK, error);

                Dll = Tmp; Dlr = (FLOAT)1.0 - Tmp;

                for (j = 0; j < I; j++) {
                    error = ComponentConditionalCdf(i, X[j].a, Cinv, LooseTheta, &CmpConCdf[0]);

                    E_CHECK(error != E_OK, error);

                    error = ComponentConditionalCdf(i, X[j].b, Cinv, LooseTheta, &CmpConCdf[1]);

                    E_CHECK(error != E_OK, error);

                    Tmp = CmpConCdf[1] - CmpConCdf[0];

                    if (X[j].s == 1) Dlr -= Tmp; else Dll -= Tmp;
                }

                if (Dll < Dlr) {
                    if (s < 0) break; else s = 1; ym -= Mode[i].h;
                }
                else {
                    if (s > 0) break; else s = -1; ym += Mode[i].h;
                }

                jj++;
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
                Error = E_OK;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm;
                }
                else {
                    flmax = flm;
                }
            }

            ii++;
        }

E1:     if (X) free(X);
    }

    LooseTheta->Theta_[3][0] = logCdet;

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        LooseTheta->Theta_[3][0] += (FLOAT)log(LooseTheta->Theta_[1][o]); LooseTheta->Theta_[2][o] = Cinv[o] / LooseTheta->Theta_[1][o];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(LooseTheta->Theta_[1][o] * LooseTheta->Theta_[1][r]);

            LooseTheta->Theta_[1][p] = LooseTheta->Theta_[1][q] = C[q] * Stdev;

            LooseTheta->Theta_[2][p] = LooseTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

EEXIT:

    if (Cinv) free(Cinv);

    if (C) free(C);

    if (D) free(D);

    if (N) free(N);

    if (Mode) free(Mode);

    E_RETURN(Error);
} // RoughEstimationKNN

// Rough component parameter estimation for kernel density estimation.

INT Rebmvnorm::RoughEstimationKDE(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                  FLOAT                *h,          // Sides of the hypersquare.
                                  FLOAT                nl,          // Total number of observations in class l.
                                  INT                  m,           // Mode index.
                                  CompnentDistribution *RigidTheta, // Rigid parameters.
                                  CompnentDistribution *LooseTheta) // Loose parameters.
{
    RoughParameterType *Mode = NULL;
    Interval           *X = NULL;
    FLOAT              *C = NULL, *Cinv = NULL, CmpConCdf[2], Dll, Dlm, Dlmin, Dlr, epsilon, flm, flmax, flmin, logCdet, logflm, logV, Mean, Stdev, Sum, Tmp, ym;
    INT                i, I, ii, j, jj, l, *N = NULL, o, p, q, r, s, error, Error = E_OK;

    // Global mode.

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    E_CHECK(NULL == Mode, E_MEM);

    N = (INT*)malloc(length_pdf_ * sizeof(INT));

    E_CHECK(NULL == N, E_MEM);

    logV = (FLOAT)0.0;

    for (i = 0; i < length_pdf_; i++) {
        logV += (FLOAT)log(h[i]); N[i] = 0;

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                for (l = 0; l < length_pdf_; l++) if ((i != l) && ((FLOAT)fabs(Y[l][j] - Y[l][m]) > (FLOAT)0.5 * h[l])) goto S0;

                Mode[i].klm += Y[length_pdf_][j]; X_[i][N[i]] = (FLOAT)j; N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                X_[i][N[i]] = (FLOAT)j; N[i] += 1;
            }
        }

        Mode[i].ym = Y[i][m]; Mode[i].flm = Y[length_pdf_][m] * Y[length_pdf_ + 1][m] / (Mode[i].klm * h[i]);
    }

    logflm = (FLOAT)log(Y[length_pdf_][m] * Y[length_pdf_ + 1][m] / nl) - logV;

    // Variance-covariance matrix.

    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == C, E_MEM);

    Cinv = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == Cinv, E_MEM);

    if (nl > length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            Sum = (FLOAT)0.0;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                Sum += Y[length_pdf_][j] * (Y[i][j] - Mode[i].ym) * (Y[i][j] - Mode[i].ym);
            }

            C[i * length_pdf_ + i] = Sum / nl;

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;

                for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                    Sum += Y[length_pdf_][j] * (Y[i][j] - Mode[i].ym) * (Y[ii][j] - Mode[ii].ym);
                }

                C[i * length_pdf_ + ii] = C[ii * length_pdf_ + i] = Sum / nl;
            }
        }

        // Correlation matrix.

        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i;

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

                C[p] = C[q] /= (FLOAT)sqrt(C[o] * C[r]);

                if (IsNan(C[q]) || IsInf(C[q])) {
                    C[p] = C[q] = (FLOAT)0.0;
                }
            }
        }

        for (i = 0; i < length_pdf_; i++) {
            C[i * length_pdf_ + i] = (FLOAT)1.0;
        }

        Error = Cholinvdet(length_pdf_, C, Cinv, &logCdet);

        E_CHECK(Error != E_OK, Error);
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0;

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
            }
        }

        logCdet = (FLOAT)0.0;
    }

    // Rigid restraints.

    RigidTheta->Theta_[3][0] = logCdet;

    for (i = 0; i < length_pdf_; i++) {
        RigidTheta->Theta_[0][i] = Mode[i].ym;

        Stdev = (FLOAT)1.0 / (SqrtPi2 * Mode[i].flm); Stdev *= Stdev;

        o = i * length_pdf_ + i;

        RigidTheta->Theta_[1][o] = Cinv[o] * Stdev; RigidTheta->Theta_[2][o] = (FLOAT)1.0 / Stdev;

        RigidTheta->Theta_[3][0] += (FLOAT)log(RigidTheta->Theta_[1][o]);

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(RigidTheta->Theta_[1][o] * RigidTheta->Theta_[1][r]);

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] = C[q] * Stdev;

            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

    epsilon = (FLOAT)exp(-(FLOAT)2.0 * (LogSqrtPi2 + logflm / length_pdf_) - RigidTheta->Theta_[3][0] / length_pdf_);

    RigidTheta->Theta_[3][0] += length_pdf_ * (FLOAT)log(epsilon);

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        RigidTheta->Theta_[1][o] *= epsilon;
        RigidTheta->Theta_[2][o] /= epsilon;

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] *= epsilon;
            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] /= epsilon;
        }
    }

    Error = LooseTheta->Memmove(RigidTheta);

    E_CHECK(Error != E_OK, Error);

    if ((Restraints_ == rtRigid) || (nl <= length_pdf_)) goto EEXIT;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        X = (Interval*)malloc((N[i] + 1) * sizeof(Interval));

        E_CHECK(NULL == X, E_MEM);

        // Bracketing.

        Tmp = (FLOAT)0.5 * h[i] * ((FLOAT)1.0 + Eps);

        for (j = 0; j < N[i]; j++) {
            l = (INT)X_[i][j];

            X[j].a = Y[i][l] - Tmp;
            X[j].b = Y[i][l] + Tmp;
        }

        I = N[i]; MergeIntervals(Mode[i].ym, &I, X);

        Error = ComponentConditionalCdf(i, Mode[i].ym, Cinv, LooseTheta, &Tmp);

        E_CHECK(Error != E_OK, Error);

        Dll = Tmp; Dlr = (FLOAT)1.0 - Tmp;

        for (j = 0; j < I; j++) {
            Error = ComponentConditionalCdf(i, X[j].a, Cinv, LooseTheta, &CmpConCdf[0]);

            E_CHECK(Error != E_OK, Error);

            Error = ComponentConditionalCdf(i, X[j].b, Cinv, LooseTheta, &CmpConCdf[1]);

            E_CHECK(Error != E_OK, Error);

            Tmp = CmpConCdf[1] - CmpConCdf[0];

            if (X[j].s == 1) Dlr -= Tmp; else Dll -= Tmp;
        }

        Dlm = Dll + Dlr;

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0; flmax = Mode[i].flm;

        // Bisection.

        ii = 1; Error = E_CON;

        while ((ii <= ItMax) && (Error != E_OK)) {
            ym = Mode[i].ym; flm = (flmax + flmin) / (FLOAT)2.0;

            jj = 1; s = 0;

            while (jj <= ItMax) {
                Mean = ym;

                LooseTheta->Theta_[0][i] = Mean;

                Stdev = (FLOAT)1.0 / (SqrtPi2 * flm); Stdev *= Stdev;

                o = i * length_pdf_ + i;

                LooseTheta->Theta_[1][o] = Cinv[o] * Stdev;

                error = ComponentConditionalCdf(i, Mode[i].ym, Cinv, LooseTheta, &Tmp);

                E_CHECK(error != E_OK, error);

                Dll = Tmp; Dlr = (FLOAT)1.0 - Tmp;

                for (j = 0; j < I; j++) {
                    error = ComponentConditionalCdf(i, X[j].a, Cinv, LooseTheta, &CmpConCdf[0]);

                    E_CHECK(error != E_OK, error);

                    error = ComponentConditionalCdf(i, X[j].b, Cinv, LooseTheta, &CmpConCdf[1]);

                    E_CHECK(error != E_OK, error);

                    Tmp = CmpConCdf[1] - CmpConCdf[0];

                    if (X[j].s == 1) Dlr -= Tmp; else Dll -= Tmp;
                }

                if (Dll < Dlr) {
                    if (s < 0) break; else s = 1; ym -= Mode[i].h;
                }
                else {
                    if (s > 0) break; else s = -1; ym += Mode[i].h;
                }

                jj++;
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
                Error = E_OK;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm;
                }
                else {
                    flmax = flm;
                }
            }

            ii++;
        }

E1:     if (X) free(X);
    }

    LooseTheta->Theta_[3][0] = logCdet;

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        LooseTheta->Theta_[3][0] += (FLOAT)log(LooseTheta->Theta_[1][o]); LooseTheta->Theta_[2][o] = Cinv[o] / LooseTheta->Theta_[1][o];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(LooseTheta->Theta_[1][o] * LooseTheta->Theta_[1][r]);

            LooseTheta->Theta_[1][p] = LooseTheta->Theta_[1][q] = C[q] * Stdev;

            LooseTheta->Theta_[2][p] = LooseTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

EEXIT:

    if (Cinv) free(Cinv);

    if (C) free(C);

    if (N) free(N);

    if (Mode) free(Mode);

    E_RETURN(Error);
} // RoughEstimationKDE

// Rough component parameter estimation for histogram.

INT Rebmvnorm::RoughEstimationH(INT                  k,           // Total number of bins.
                                FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl].
                                FLOAT                *h,          // Sides of the hypersquare.
                                FLOAT                nl,          // Total number of observations in class l.
                                INT                  m,           // Mode index.
                                CompnentDistribution *RigidTheta, // Rigid parameters.
                                CompnentDistribution *LooseTheta) // Loose parameters.
{
    RoughParameterType *Mode = NULL;
    Interval           *X = NULL;
    FLOAT              *C = NULL, *Cinv = NULL, CmpConCdf[2], Dll, Dlm, Dlmin, Dlr, epsilon, flm, flmax, flmin, logCdet, logflm, Mean, logV, Stdev, Sum, Tmp, ym;
    INT                i, I, ii, j, jj, l, *N = NULL, o, p, q, r, s, error, Error = E_OK;

    // Global mode.

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    E_CHECK(NULL == Mode, E_MEM);

    N = (INT*)malloc(length_pdf_ * sizeof(INT));

    E_CHECK(NULL == N, E_MEM);

    logV = (FLOAT)0.0;

    for (i = 0; i < length_pdf_; i++) {
        logV += (FLOAT)log(h[i]); N[i] = 0;

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                for (l = 0; l < length_pdf_; l++) if ((i != l) && (Y[l][j] != Y[l][m])) goto S0;

                Mode[i].klm += Y[length_pdf_][j]; X_[i][N[i]] = (FLOAT)j; N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < k; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                X_[i][N[i]] = (FLOAT)j; N[i] += 1;
            }
        }

        Mode[i].ym = Y[i][m]; Mode[i].flm = Y[length_pdf_][m] / (Mode[i].klm * h[i]);
    }

    logflm = (FLOAT)log(Y[length_pdf_][m] / nl) - logV;

    // Variance-covariance matrix.

    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == C, E_MEM);

    Cinv = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == Cinv, E_MEM);

    if (nl > length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            Sum = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                Sum += Y[length_pdf_][j] * (Y[i][j] - Mode[i].ym) * (Y[i][j] - Mode[i].ym);
            }

            C[i * length_pdf_ + i] = Sum / nl;

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;

                for (j = 0; j < k; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                    Sum += Y[length_pdf_][j] * (Y[i][j] - Mode[i].ym) * (Y[ii][j] - Mode[ii].ym);
                }

                C[i * length_pdf_ + ii] = C[ii * length_pdf_ + i] = Sum / nl;
            }
        }

        // Correlation matrix.

        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i;

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

                C[p] = C[q] /= (FLOAT)sqrt(C[o] * C[r]);

                if (IsNan(C[q]) || IsInf(C[q])) {
                    C[p] = C[q] = (FLOAT)0.0;
                }
            }
        }

        for (i = 0; i < length_pdf_; i++) {
            C[i * length_pdf_ + i] = (FLOAT)1.0;
        }

        Error = Cholinvdet(length_pdf_, C, Cinv, &logCdet);

        E_CHECK(Error != E_OK, Error);
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0;

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
            }
        }

        logCdet = (FLOAT)0.0;
    }

    // Rigid restraints.

    RigidTheta->Theta_[3][0] = logCdet;

    for (i = 0; i < length_pdf_; i++) {
        RigidTheta->Theta_[0][i] = Mode[i].ym;

        Stdev = (FLOAT)1.0 / (SqrtPi2 * Mode[i].flm); Stdev *= Stdev;

        o = i * length_pdf_ + i;

        RigidTheta->Theta_[1][o] = Cinv[o] * Stdev; RigidTheta->Theta_[2][o] = (FLOAT)1.0 / Stdev;

        RigidTheta->Theta_[3][0] += (FLOAT)log(RigidTheta->Theta_[1][o]);

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(RigidTheta->Theta_[1][o] * RigidTheta->Theta_[1][r]);

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] = C[q] * Stdev;

            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

    epsilon = (FLOAT)exp(-(FLOAT)2.0 * (LogSqrtPi2 + logflm / length_pdf_) - RigidTheta->Theta_[3][0] / length_pdf_);

    RigidTheta->Theta_[3][0] += length_pdf_ * (FLOAT)log(epsilon);

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        RigidTheta->Theta_[1][o] *= epsilon;
        RigidTheta->Theta_[2][o] /= epsilon;

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] *= epsilon;
            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] /= epsilon;
        }
    }
    
    Error = LooseTheta->Memmove(RigidTheta);

    E_CHECK(Error != E_OK, Error);

    if ((Restraints_ == rtRigid) || (nl <= length_pdf_)) goto EEXIT;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        X = (Interval*)malloc((N[i] + 1) * sizeof(Interval));

        E_CHECK(NULL == X, E_MEM);

        // Bracketing.

        Tmp = (FLOAT)0.5 * h[i] * ((FLOAT)1.0 + Eps);

        for (j = 0; j < N[i]; j++) {
            l = (INT)X_[i][j];

            X[j].a = Y[i][l] - Tmp;
            X[j].b = Y[i][l] + Tmp;
        }

        I = N[i]; MergeIntervals(Mode[i].ym, &I, X);

        Error = ComponentConditionalCdf(i, Mode[i].ym, Cinv, LooseTheta, &Tmp);

        E_CHECK(Error != E_OK, Error);

        Dll = Tmp; Dlr = (FLOAT)1.0 - Tmp;

        for (j = 0; j < I; j++) {
            Error = ComponentConditionalCdf(i, X[j].a, Cinv, LooseTheta, &CmpConCdf[0]);

            E_CHECK(Error != E_OK, Error);

            Error = ComponentConditionalCdf(i, X[j].b, Cinv, LooseTheta, &CmpConCdf[1]);

            E_CHECK(Error != E_OK, Error);

            Tmp = CmpConCdf[1] - CmpConCdf[0];

            if (X[j].s == 1) Dlr -= Tmp; else Dll -= Tmp;
        }

        Dlm = Dll + Dlr;

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0; flmax = Mode[i].flm;

        // Bisection.

        ii = 1; Error = E_CON;

        while ((ii <= ItMax) && (Error != E_OK)) {
            ym = Mode[i].ym; flm = (flmax + flmin) / (FLOAT)2.0;

            jj = 1; s = 0;

            while (jj <= ItMax) {
                Mean = ym;

                LooseTheta->Theta_[0][i] = Mean;

                Stdev = (FLOAT)1.0 / (SqrtPi2 * flm); Stdev *= Stdev;

                o = i * length_pdf_ + i;

                LooseTheta->Theta_[1][o] = Cinv[o] * Stdev;

                error = ComponentConditionalCdf(i, Mode[i].ym, Cinv, LooseTheta, &Tmp);

                E_CHECK(error != E_OK, error);

                Dll = Tmp; Dlr = (FLOAT)1.0 - Tmp;

                for (j = 0; j < I; j++) {
                    error = ComponentConditionalCdf(i, X[j].a, Cinv, LooseTheta, &CmpConCdf[0]);

                    E_CHECK(error != E_OK, error);

                    error = ComponentConditionalCdf(i, X[j].b, Cinv, LooseTheta, &CmpConCdf[1]);

                    E_CHECK(error != E_OK, error);

                    Tmp = CmpConCdf[1] - CmpConCdf[0];

                    if (X[j].s == 1) Dlr -= Tmp; else Dll -= Tmp;
                }

                if (Dll < Dlr) {
                    if (s < 0) break; else s = 1; ym -= Mode[i].h;
                }
                else {
                    if (s > 0) break; else s = -1; ym += Mode[i].h;
                }

                jj++;
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
                Error = E_OK;
            }
            else {
                if (Dlm * Dlmin > (FLOAT)0.0) {
                    flmin = flm; Dlmin = Dlm;
                }
                else {
                    flmax = flm;
                }
            }

            ii++;
        }

E1:     if (X) free(X);
    }

    LooseTheta->Theta_[3][0] = logCdet;

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        LooseTheta->Theta_[3][0] += (FLOAT)log(LooseTheta->Theta_[1][o]); LooseTheta->Theta_[2][o] = Cinv[o] / LooseTheta->Theta_[1][o];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(LooseTheta->Theta_[1][o] * LooseTheta->Theta_[1][r]);

            LooseTheta->Theta_[1][p] = LooseTheta->Theta_[1][q] = C[q] * Stdev;

            LooseTheta->Theta_[2][p] = LooseTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

EEXIT:

    if (Cinv) free(Cinv);

    if (C) free(C);

    if (N) free(N);

    if (Mode) free(Mode);

    E_RETURN(Error);
} // RoughEstimationH

// Enhanced component parameter estimation for k-nearest neighbours.

INT Rebmvnorm::EnhancedEstimationKNN(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,logV,R].
                                     FLOAT                nl,          // Total number of observations in class l.
                                     CompnentDistribution *RigidTheta, // Rigid parameters.
                                     CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                Sum;
    INT                  i, ii, j, o, Error = E_OK;

    EnhanTheta = new CompnentDistribution(this);

    E_CHECK(NULL == EnhanTheta, E_MEM);

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    E_CHECK(Error != E_OK, Error);

    E_CHECK(nl <= (FLOAT)1.0, E_ARG);

    for (i = 0; i < length_pdf_; i++) {
        EnhanTheta->pdf_[i] = pfNormal;

        Sum = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
            Sum += Y[length_pdf_][j] * Y[i][j];
        }

        EnhanTheta->Theta_[0][i] = Sum / nl;

        o = i * length_pdf_ + i;

        Sum = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
            Sum += Y[length_pdf_][j] * (Y[i][j] - EnhanTheta->Theta_[0][i]) * (Y[i][j] - EnhanTheta->Theta_[0][i]);
        }

        EnhanTheta->Theta_[1][o] = Sum / nl;

        for (ii = 0; ii < i; ii++) {
            Sum = (FLOAT)0.0;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                Sum += Y[length_pdf_][j] * (Y[i][j] - EnhanTheta->Theta_[0][i]) * (Y[ii][j] - EnhanTheta->Theta_[0][ii]);
            }

            EnhanTheta->Theta_[1][i * length_pdf_ + ii] = EnhanTheta->Theta_[1][ii * length_pdf_ + i] = Sum / nl;
        }
    }

    Error = Cholinvdet(length_pdf_, EnhanTheta->Theta_[1], EnhanTheta->Theta_[2], EnhanTheta->Theta_[3]);

    E_CHECK(Error != E_OK, Error);

    E_CHECK(*EnhanTheta->Theta_[3] < *RigidTheta->Theta_[3] + (FLOAT)log(var_mul_), E_CON);

    Error = LooseTheta->Memmove(EnhanTheta);

    E_CHECK(Error != E_OK, Error);

EEXIT:

    if (EnhanTheta) delete EnhanTheta;

    E_RETURN(Error);
} // EnhancedEstimationKNN

// Enhanced component parameter estimation for kernel density estimation.

INT Rebmvnorm::EnhancedEstimationKDE(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                     FLOAT                nl,          // Total number of observations in class l.
                                     CompnentDistribution *RigidTheta, // Rigid parameters.
                                     CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                Sum;
    INT                  i, ii, j, o, Error = E_OK;

    EnhanTheta = new CompnentDistribution(this);

    E_CHECK(NULL == EnhanTheta, E_MEM);

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    E_CHECK(Error != E_OK, Error);

    E_CHECK(nl <= (FLOAT)1.0, E_ARG);

    for (i = 0; i < length_pdf_; i++) {
        EnhanTheta->pdf_[i] = pfNormal;

        Sum = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
            Sum += Y[length_pdf_][j] * Y[i][j];
        }

        EnhanTheta->Theta_[0][i] = Sum / nl;

        o = i * length_pdf_ + i;

        Sum = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
            Sum += Y[length_pdf_][j] * (Y[i][j] - EnhanTheta->Theta_[0][i]) * (Y[i][j] - EnhanTheta->Theta_[0][i]);
        }

        EnhanTheta->Theta_[1][o] = Sum / nl;

        for (ii = 0; ii < i; ii++) {
            Sum = (FLOAT)0.0;

            for (j = 0; j < nr_; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                Sum += Y[length_pdf_][j] * (Y[i][j] - EnhanTheta->Theta_[0][i]) * (Y[ii][j] - EnhanTheta->Theta_[0][ii]);
            }

            EnhanTheta->Theta_[1][i * length_pdf_ + ii] = EnhanTheta->Theta_[1][ii * length_pdf_ + i] = Sum / nl;
        }
    }

    Error = Cholinvdet(length_pdf_, EnhanTheta->Theta_[1], EnhanTheta->Theta_[2], EnhanTheta->Theta_[3]);

    E_CHECK(Error != E_OK, Error);

    E_CHECK(*EnhanTheta->Theta_[3] < *RigidTheta->Theta_[3] + (FLOAT)log(var_mul_), E_CON);

    Error = LooseTheta->Memmove(EnhanTheta);

    E_CHECK(Error != E_OK, Error);

EEXIT:

    if (EnhanTheta) delete EnhanTheta;

    E_RETURN(Error);
} // EnhancedEstimationKDE

// Enhanced component parameter estimation for histogram.

INT Rebmvnorm::EnhancedEstimationH(INT                  k,           // Total number of bins.
                                   FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                   FLOAT                nl,          // Total number of observations in class l.
                                   FLOAT                *h,          // Sides of the hypersquare.
                                   CompnentDistribution *RigidTheta, // Rigid parameters.
                                   CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                Sum;
    INT                  i, ii, j, o, Error = E_OK;

    (void)h;

    EnhanTheta = new CompnentDistribution(this);

    E_CHECK(NULL == EnhanTheta, E_MEM);

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    E_CHECK(Error != E_OK, Error);

    E_CHECK(nl <= (FLOAT)1.0, E_ARG);

    for (i = 0; i < length_pdf_; i++) {
        EnhanTheta->pdf_[i] = pfNormal;

        Sum = (FLOAT)0.0;

        for (j = 0; j < k; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
            Sum += Y[length_pdf_][j] * Y[i][j];
        }

        EnhanTheta->Theta_[0][i] = Sum / nl;

        o = i * length_pdf_ + i;

        Sum = (FLOAT)0.0;

        for (j = 0; j < k; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
            Sum += Y[length_pdf_][j] * (Y[i][j] - EnhanTheta->Theta_[0][i]) * (Y[i][j] - EnhanTheta->Theta_[0][i]);
        }

        EnhanTheta->Theta_[1][o] = Sum / nl;

        for (ii = 0; ii < i; ii++) {
            Sum = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[length_pdf_][j] > FLOAT_MIN) {
                Sum += Y[length_pdf_][j] * (Y[i][j] - EnhanTheta->Theta_[0][i]) * (Y[ii][j] - EnhanTheta->Theta_[0][ii]);
            }

            EnhanTheta->Theta_[1][i * length_pdf_ + ii] = EnhanTheta->Theta_[1][ii * length_pdf_ + i] = Sum / nl;
        }
    }

    Error = Cholinvdet(length_pdf_, EnhanTheta->Theta_[1], EnhanTheta->Theta_[2], EnhanTheta->Theta_[3]);

    E_CHECK(Error != E_OK, Error);

    E_CHECK(*EnhanTheta->Theta_[3] < *RigidTheta->Theta_[3] + (FLOAT)log(var_mul_), E_CON);

    Error = LooseTheta->Memmove(EnhanTheta);

    E_CHECK(Error != E_OK, Error);

EEXIT:

    if (EnhanTheta) delete EnhanTheta;

    E_RETURN(Error);
} // EnhancedEstimationH

// Moments calculation.

INT Rebmvnorm::MomentsCalculation(CompnentDistribution *CmpTheta, // Component parameters.
                                  FLOAT                *FirstM,   // First moment.
                                  FLOAT                *SecondM)  // Second moment.
{
    INT i, ii, o, p, q, Error = E_OK;

    for (i = 0; i < length_pdf_; i++) {
        FirstM[i] = CmpTheta->Theta_[0][i];

        o = i * length_pdf_ + i;

        SecondM[o] = CmpTheta->Theta_[1][o] + CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][i];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

            SecondM[p] = SecondM[q] = CmpTheta->Theta_[1][p] + CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][ii];
        }
    }

    E_RETURN(Error);
} // MomentsCalculation

// Bayes classification of the remaining observations for k-nearest neighbour.

INT Rebmvnorm::BayesClassificationKNN(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                      INT                  c,          // Number of components.
                                      FLOAT                *W,         // Component weights.
                                      CompnentDistribution **MixTheta, // Mixture parameters.
                                      FLOAT                **FirstM,   // First moments.
                                      FLOAT                **SecondM)  // Second moments.
{
    FLOAT CmpDist, dW, Max, N = (FLOAT)0.0, Tmp;
    INT   i, j, jj, l, o, outlier, Outlier = 0, p, q, Error = E_OK;

    for (i = 0; i < nr_; i++) {
        if (Y[length_pdf_][i] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(i, Y, MixTheta[l], &CmpDist, &outlier);

            E_CHECK(Error != E_OK, Error);

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(i, Y, MixTheta[j], &CmpDist, &outlier);

                E_CHECK(Error != E_OK, Error);

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[length_pdf_][i];
            }
            else {
                dW = Y[length_pdf_][i] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    FirstM[l][j] += dW * (Y[j][i] - FirstM[l][j]) / W[l];

                    o = j * length_pdf_ + j;

                    SecondM[l][o] += dW * (Y[j][i] * Y[j][i] - SecondM[l][o]) / W[l];

                    for (jj = 0; jj < j; jj++) {
                        p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                        SecondM[l][p] = SecondM[l][q] += dW * (Y[j][i] * Y[jj][i] - SecondM[l][q]) / W[l];
                    }
                }
            }
        }
    }

    for (i = 0; i < c; i++) {
        W[i] *= n_ / (n_ - N);

        for (j = 0; j < length_pdf_; j++) {
            MixTheta[i]->Theta_[0][j] = FirstM[i][j];

            o = j * length_pdf_ + j;

            MixTheta[i]->Theta_[1][o] = SecondM[i][o] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][j];

            for (jj = 0; jj < j; jj++) {
                p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                MixTheta[i]->Theta_[1][p] = MixTheta[i]->Theta_[1][q] = SecondM[i][p] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][jj];
            }
        }

        Error = Cholinvdet(length_pdf_, MixTheta[i]->Theta_[1], MixTheta[i]->Theta_[2], MixTheta[i]->Theta_[3]);

        E_CHECK(Error != E_OK, Error);
    }

EEXIT:

    E_RETURN(Error);
} // BayesClassificationKNN

// Bayes classification of the remaining observations for kernel density estimation.

INT Rebmvnorm::BayesClassificationKDE(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                      INT                  c,          // Number of components.
                                      FLOAT                *W,         // Component weights.
                                      CompnentDistribution **MixTheta, // Mixture parameters.
                                      FLOAT                **FirstM,   // First moments.
                                      FLOAT                **SecondM)  // Second moments.
{
    FLOAT CmpDist, dW, Max, N = (FLOAT)0.0, Tmp;
    INT   i, j, jj, l, o, outlier, Outlier = 0, p, q, Error = E_OK;

    for (i = 0; i < nr_; i++) {
        if (Y[length_pdf_][i] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(i, Y, MixTheta[l], &CmpDist, &outlier);

            E_CHECK(Error != E_OK, Error);

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(i, Y, MixTheta[j], &CmpDist, &outlier);

                E_CHECK(Error != E_OK, Error);

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[length_pdf_][i];
            }
            else {
                dW = Y[length_pdf_][i] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    FirstM[l][j] += dW * (Y[j][i] - FirstM[l][j]) / W[l];

                    o = j * length_pdf_ + j;

                    SecondM[l][o] += dW * (Y[j][i] * Y[j][i] - SecondM[l][o]) / W[l];

                    for (jj = 0; jj < j; jj++) {
                        p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                        SecondM[l][p] = SecondM[l][q] += dW * (Y[j][i] * Y[jj][i] - SecondM[l][q]) / W[l];
                    }
                }
            }
        }
    }

    for (i = 0; i < c; i++) {
        W[i] *= n_ / (n_ - N);

        for (j = 0; j < length_pdf_; j++) {
            MixTheta[i]->Theta_[0][j] = FirstM[i][j];

            o = j * length_pdf_ + j;

            MixTheta[i]->Theta_[1][o] = SecondM[i][o] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][j];

            for (jj = 0; jj < j; jj++) {
                p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                MixTheta[i]->Theta_[1][p] = MixTheta[i]->Theta_[1][q] = SecondM[i][p] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][jj];
            }
        }

        Error = Cholinvdet(length_pdf_, MixTheta[i]->Theta_[1], MixTheta[i]->Theta_[2], MixTheta[i]->Theta_[3]);

        E_CHECK(Error != E_OK, Error);
    }

EEXIT:

    E_RETURN(Error);
} // BayesClassificationKDE

// Bayes classification of the remaining observations for histogram.

INT Rebmvnorm::BayesClassificationH(INT                  k,          // Total number of bins.
                                    FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                    INT                  c,          // Number of components.
                                    FLOAT                *W,         // Component weights.
                                    CompnentDistribution **MixTheta, // Mixture parameters.
                                    FLOAT                **FirstM,   // First moments.
                                    FLOAT                **SecondM)  // Second moments.
{
    FLOAT CmpDist, dW, Max, N = (FLOAT)0.0, Tmp;
    INT   i, j, jj, l, o, outlier, Outlier = 0, p, q, Error = E_OK;

    for (i = 0; i < k; i++) {
        if (Y[length_pdf_][i] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(i, Y, MixTheta[l], &CmpDist, &outlier);

            E_CHECK(Error != E_OK, Error);

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(i, Y, MixTheta[j], &CmpDist, &outlier);

                E_CHECK(Error != E_OK, Error);

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[length_pdf_][i];
            }
            else {
                dW = Y[length_pdf_][i] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    FirstM[l][j] += dW * (Y[j][i] - FirstM[l][j]) / W[l];

                    o = j * length_pdf_ + j;

                    SecondM[l][o] += dW * (Y[j][i] * Y[j][i] - SecondM[l][o]) / W[l];

                    for (jj = 0; jj < j; jj++) {
                        p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                        SecondM[l][p] = SecondM[l][q] += dW * (Y[j][i] * Y[jj][i] - SecondM[l][q]) / W[l];
                    }
                }
            }
        }
    }

    for (i = 0; i < c; i++) {
        W[i] *= n_ / (n_ - N);

        for (j = 0; j < length_pdf_; j++) {
            MixTheta[i]->Theta_[0][j] = FirstM[i][j];

            o = j * length_pdf_ + j;

            MixTheta[i]->Theta_[1][o] = SecondM[i][o] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][j];

            for (jj = 0; jj < j; jj++) {
                p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                MixTheta[i]->Theta_[1][p] = MixTheta[i]->Theta_[1][q] = SecondM[i][p] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][jj];
            }
        }

        Error = Cholinvdet(length_pdf_, MixTheta[i]->Theta_[1], MixTheta[i]->Theta_[2], MixTheta[i]->Theta_[3]);

        E_CHECK(Error != E_OK, Error);
    }

EEXIT:

    E_RETURN(Error);
} // BayesClassificationH

// Returns component p.d.f..

INT Rebmvnorm::ComponentDist(INT                  j,         // Indey of observation.  
                             FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                             CompnentDistribution *CmpTheta, // Component parameters.
                             FLOAT                *CmpDist,  // Component distribution.
                             INT                  *Outlier)  // 1 if outlier otherwise 0.
{
    FLOAT y, yi, yk;
    INT   i, k, Error = E_OK;

    y = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        yi = Y[i][j] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + i] * yi * yi;

        for (k = i + 1; k < CmpTheta->length_pdf_; k++) {
            yk = Y[k][j] - CmpTheta->Theta_[0][k]; y += CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + k] * yi * yk;
        }
    }

    if (Outlier) {
        *Outlier = (FLOAT)2.0 * y > ChiSqr_;
    }

    *CmpDist = (FLOAT)exp(-y - CmpTheta->length_pdf_ * LogSqrtPi2 - (FLOAT)0.5 * CmpTheta->Theta_[3][0]);

    E_RETURN(Error);
} // ComponentDist

// Returns logarithm of component p.d.f..

INT Rebmvnorm::LogComponentDist(INT                  j,         // Indey of observation.  
                                FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                                CompnentDistribution *CmpTheta, // Component parameters.
                                FLOAT                *CmpDist,  // Component distribution.
                                INT                  *Outlier)  // 1 if outlier otherwise 0.
{
    FLOAT y, yi, yk;
    INT   i, k, Error = E_OK;

    y = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        yi = Y[i][j] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + i] * yi * yi;

        for (k = i + 1; k < CmpTheta->length_pdf_; k++) {
            yk = Y[k][j] - CmpTheta->Theta_[0][k]; y += CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + k] * yi * yk;
        }
    }

    if (Outlier) {
        *Outlier = (FLOAT)2.0 * y > ChiSqr_;
    }

    *CmpDist = -y - CmpTheta->length_pdf_ * LogSqrtPi2 - (FLOAT)0.5 * CmpTheta->Theta_[3][0];

    E_RETURN(Error);
} // LogComponentDist

INT Rebmvnorm::DegreesOffreedom(INT c,                  // Number of components.
                                CompnentDistribution**, // Mixture parameters.
                                INT *M)                 // Degrees of freedom.
{
    INT i, Error = E_OK;

    *M = c - 1;

    for (i = 0; i < c; i++) {
        *M += length_pdf_ + (length_pdf_ + 1) * length_pdf_ / 2;
    }

    E_RETURN(Error);
} // DegreesOffreedom

/// Panic Branislav
// Runs the EM algorithm or its variant for the initial parameters assesed by the REBMIX algorithm.

INT Rebmvnorm::EMInitialize()
{
    INT Error = E_OK;

    EM_ = new Emmvnorm();

    E_CHECK(NULL == EM_, E_MEM);

    Error = EM_->Initialize(n_,
        nr_,
        nc_,
        Y_,
        cmax_,
        length_pdf_,
        length_Theta_,
        length_theta_,
        EM_TOL_,
        EM_am_,
        EM_max_iter_,
        EM_K_,
        EM_strategy_,
        EM_variant_,
        EM_accel_);

    E_CHECK(Error != E_OK, Error);

EEXIT:

    E_RETURN(Error);
} // EMInitialize
/// End

