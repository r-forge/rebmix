#include <math.h>
#include <float.h> 

#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "base.h"
#include "rebmvnormf.h"

// Perform necessary initializations.

int Rebmvnorm::Initialize()
{
    int Error = 0;

    p_value_ = (FLOAT)0.0001;

    min_dist_mul_ = (FLOAT)2.5;

    var_mul_ = (FLOAT)0.25;

    Error = GammaInv((FLOAT)1.0 - (FLOAT)2.0 * p_value_, (FLOAT)2.0, length_pdf_ / (FLOAT)2.0, &ChiSqr_);

    return (Error);
} // Initialize

int Rebmvnorm::ComponentConditionalDist(int                  i,           // Index of variable y.
                                        FLOAT                *Y,          // Pointer to the input point [y0,...,yd-1].
                                        FLOAT                *Cinv,       // Inverse correlation matrix.
                                        CompnentDistribution *CmpTheta,   // Component distribution type.
                                        FLOAT                *CmpMrgDist) // Component marginal distribution.
{
    int o;
    FLOAT y, Mean, Stdev;
    int   Error = 0;

    Mean = CmpTheta->Theta_[0][i]; 

    o = i * length_pdf_ + i;

    Stdev = (FLOAT)sqrt(CmpTheta->Theta_[1][o] / Cinv[o]);

    y = (Y[i] - Mean) / (Sqrt2 * Stdev); y *= y;

    *CmpMrgDist = (FLOAT)exp(-y) / (SqrtPi2 * Stdev);

    return Error;
} // ComponentConditionalDist

// Rough component parameter estimation for k-nearest neighbours.

int Rebmvnorm::RoughEstimationKNN(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,V,R].
                                  int                  k,           // k-nearest neighbours.
                                  FLOAT                *h,          // Normalizing vector.
                                  FLOAT                nl,          // Total number of observations in class l.
                                  int                  m,           // Mode index.
                                  CompnentDistribution *RigidTheta, // Rigid parameters.
                                  CompnentDistribution *LooseTheta) // Loose parameters.
{
    int                i, ii, j, l, o, p, q, r, *N = NULL;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, flm, flmin, flmax, Dlm, Dlmin, Sum, Stdev, Dc, R, *D = NULL;
    FLOAT              *C = NULL, *Cinv = NULL, Cdet;
    int                Error = 0, Stop;

    // Global mode.

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    N = (int*)malloc(length_pdf_ * sizeof(int));

    Error = NULL == N; if (Error) goto E0;

    D = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == D; if (Error) goto E0;

    for (i = 0; i < length_pdf_; i++) {
        N[i] = 0; D[i] = (FLOAT)2.0 * Y[m][length_pdf_ + 2] * h[i];

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Dc = (FLOAT)0.0;

                for (l = 0; l < length_pdf_; l++) if (i != l) {
                    R = (Y[j][l] - Y[m][l]) / h[l]; Dc += R * R;
                }

                R = (FLOAT)sqrt(Dc);

                if (R > Y[m][length_pdf_ + 2]) goto S0;

                Mode[i].klm += Y[j][length_pdf_];

                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / D[i] + (FLOAT)0.5) * D[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * D[i]) goto S0;
                }

                N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / D[i] + (FLOAT)0.5) * D[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * D[i]) goto S1;
                }

                N[i] += 1;
S1:;
            }
        }

        Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][length_pdf_] * k / (Mode[i].klm * D[i]);
    }

    flm = Y[m][length_pdf_] * k / (nl * Y[m][length_pdf_ + 1]);

    // Variance-covariance matrix.

    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    Error = NULL == C; if (Error) goto E0;

    Cinv = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    Error = NULL == Cinv; if (Error) goto E0;

    if (nl >= length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            Sum = (FLOAT)0.0;
        
            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * (Y[j][i] - Mode[i].ym) * (Y[j][i] - Mode[i].ym);
            }

            C[i * length_pdf_ + i] = Sum / nl;

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;
        
                for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                    Sum += Y[j][length_pdf_] * (Y[j][i] - Mode[i].ym) * (Y[j][ii] - Mode[ii].ym);
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

        Stop = LUinvdet(length_pdf_, C, Cinv, &Cdet);

        if (Stop || (Cdet < (FLOAT)0.05)) {
            for (i = 0; i < length_pdf_; i++) {
                o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0; 

                for (ii = 0; ii < i; ii++) {
                    p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                    C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
                }
            }

            Cdet = (FLOAT)1.0;
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0; 

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
            }
        }

        Cdet = (FLOAT)1.0;
    }

    // Rigid restraints.

    RigidTheta->Theta_[3][0] = Cdet;

    for (i = 0; i < length_pdf_; i++) {
        RigidTheta->Theta_[0][i] = Mode[i].ym;

        Stdev = (FLOAT)1.0 / (SqrtPi2 * Mode[i].flm); Stdev *= Stdev;

        o = i * length_pdf_ + i;

        RigidTheta->Theta_[3][0] *= RigidTheta->Theta_[1][o] = Cinv[o] * Stdev; RigidTheta->Theta_[2][o] = (FLOAT)1.0 / Stdev;

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(RigidTheta->Theta_[1][o] * RigidTheta->Theta_[1][r]);

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] = C[q] * Stdev;

            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

    epsilon = (FLOAT)exp(-(FLOAT)2.0 * (LogSqrtPi2 + (FLOAT)log(flm) / length_pdf_) - (FLOAT)log(RigidTheta->Theta_[3][0]) / length_pdf_);

    if (epsilon > (FLOAT)1.0) {
        RigidTheta->Theta_[3][0] *= (FLOAT)exp(length_pdf_ * (FLOAT)log(epsilon));

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
    }

    Error = LooseTheta->Memmove(RigidTheta);

    if (Error) goto E0;

    if (Restraints_ == rtRigid) goto E0;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        // Bracketing.

        Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

        for (j = 0; j < N[i]; j++) {
            Error = ComponentConditionalDist(i, X_[j], Cinv, LooseTheta, &CmpMrgDist);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * D[i];
        }

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0 - (FLOAT)2.0 * p_value_; flmax = Mode[i].flm;

        // Bisection.

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            Stdev = (FLOAT)1.0 / (SqrtPi2 * flm); Stdev *= Stdev;

            o = i * length_pdf_ + i;

            LooseTheta->Theta_[1][o] = Cinv[o] * Stdev;

            Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

            for (j = 0; j < N[i]; j++) {
                Error = ComponentConditionalDist(i, X_[j], Cinv, LooseTheta, &CmpMrgDist);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * D[i];
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
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
E1:;
    }

    LooseTheta->Theta_[3][0] = Cdet;

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        LooseTheta->Theta_[3][0] *= LooseTheta->Theta_[1][o]; LooseTheta->Theta_[2][o] = Cinv[o] / LooseTheta->Theta_[1][o];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(LooseTheta->Theta_[1][o] * LooseTheta->Theta_[1][r]);

            LooseTheta->Theta_[1][p] = LooseTheta->Theta_[1][q] = C[q] * Stdev;

            LooseTheta->Theta_[2][p] = LooseTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

E0: if (Cinv) free(Cinv);

    if (C) free(C);

    if (D) free(D);

    if (N) free(N);

    if (Mode) free(Mode);

    return Error;
} // RoughEstimationKNN

// Rough component parameter estimation for Parzen window.

int Rebmvnorm::RoughEstimationPW(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                 FLOAT                *h,          // Sides of the hypersquare.
                                 FLOAT                nl,          // Total number of observations in class l.
                                 int                  m,           // Mode index.
                                 CompnentDistribution *RigidTheta, // Rigid parameters.
                                 CompnentDistribution *LooseTheta) // Loose parameters.
{
    int                i, ii, j, l, o, p, q, r, *N = NULL;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, flm, flmin, flmax, V, Dlm, Dlmin, Sum, Stdev;
    FLOAT              *C = NULL, *Cinv = NULL, Cdet;
    int                Error = 0, Stop;

    // Global mode.

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    N = (int*)malloc(length_pdf_ * sizeof(int));

    Error = NULL == N; if (Error) goto E0;

    V = (FLOAT)1.0;

    for (i = 0; i < length_pdf_; i++) {
        V *= h[i]; N[i] = 0;

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                for (l = 0; l < length_pdf_; l++) if ((i != l) && ((FLOAT)fabs(Y[j][l] - Y[m][l]) > (FLOAT)0.5 * h[l])) goto S0;

                Mode[i].klm += Y[j][length_pdf_]; 

                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / h[i] + (FLOAT)0.5) * h[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * h[i]) goto S0;
                }

                N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                X_[N[i]][i] = Y[m][i] + (int)floor((Y[j][i] - Y[m][i]) / h[i] + (FLOAT)0.5) * h[i];

                for (ii = 0; ii < N[i]; ii++) {
                    if ((FLOAT)fabs(X_[N[i]][i] - X_[ii][i]) < (FLOAT)0.5 * h[i]) goto S1;
                }

                N[i] += 1;
S1:;
            }
        }

        Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][length_pdf_] * Y[m][length_pdf_ + 1] / (Mode[i].klm * h[i]);
    }

    flm = Y[m][length_pdf_] * Y[m][length_pdf_ + 1] / (nl * V);

    // Variance-covariance matrix.

    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    Error = NULL == C; if (Error) goto E0;

    Cinv = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    Error = NULL == Cinv; if (Error) goto E0;

    if (nl >= length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            Sum = (FLOAT)0.0;
        
            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * (Y[j][i] - Mode[i].ym) * (Y[j][i] - Mode[i].ym);
            }

            C[i * length_pdf_ + i] = Sum / nl;

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;
        
                for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                    Sum += Y[j][length_pdf_] * (Y[j][i] - Mode[i].ym) * (Y[j][ii] - Mode[ii].ym);
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

        Stop = LUinvdet(length_pdf_, C, Cinv, &Cdet);

        if (Stop || (Cdet < (FLOAT)0.05)) {
            for (i = 0; i < length_pdf_; i++) {
                o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0; 

                for (ii = 0; ii < i; ii++) {
                    p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                    C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
                }
            }

            Cdet = (FLOAT)1.0;
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0; 

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
            }
        }

        Cdet = (FLOAT)1.0;
    }

    // Rigid restraints.

    RigidTheta->Theta_[3][0] = Cdet;

    for (i = 0; i < length_pdf_; i++) {
        RigidTheta->Theta_[0][i] = Mode[i].ym;

        Stdev = (FLOAT)1.0 / (SqrtPi2 * Mode[i].flm); Stdev *= Stdev;

        o = i * length_pdf_ + i;

        RigidTheta->Theta_[3][0] *= RigidTheta->Theta_[1][o] = Cinv[o] * Stdev; RigidTheta->Theta_[2][o] = (FLOAT)1.0 / Stdev;

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(RigidTheta->Theta_[1][o] * RigidTheta->Theta_[1][r]);

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] = C[q] * Stdev;

            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

    epsilon = (FLOAT)exp(-(FLOAT)2.0 * (LogSqrtPi2 + (FLOAT)log(flm) / length_pdf_) - (FLOAT)log(RigidTheta->Theta_[3][0]) / length_pdf_);

    if (epsilon > (FLOAT)1.0) {
        RigidTheta->Theta_[3][0] *= (FLOAT)exp(length_pdf_ * (FLOAT)log(epsilon));

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
    }

    Error = LooseTheta->Memmove(RigidTheta);

    if (Error) goto E0;

    if (Restraints_ == rtRigid) goto E0;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        // Bracketing.

        Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

        for (j = 0; j < N[i]; j++) {
            Error = ComponentConditionalDist(i, X_[j], Cinv, LooseTheta, &CmpMrgDist);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * h[i];
        }

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0 - (FLOAT)2.0 * p_value_; flmax = Mode[i].flm;

        // Bisection.

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            Stdev = (FLOAT)1.0 / (SqrtPi2 * flm); Stdev *= Stdev;

            o = i * length_pdf_ + i;

            LooseTheta->Theta_[1][o] = Cinv[o] * Stdev;

            Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

            for (j = 0; j < N[i]; j++) {
                Error = ComponentConditionalDist(i, X_[j], Cinv, LooseTheta, &CmpMrgDist);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * h[i];
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
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
E1:;
    }

    LooseTheta->Theta_[3][0] = Cdet;

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        LooseTheta->Theta_[3][0] *= LooseTheta->Theta_[1][o]; LooseTheta->Theta_[2][o] = Cinv[o] / LooseTheta->Theta_[1][o];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(LooseTheta->Theta_[1][o] * LooseTheta->Theta_[1][r]);

            LooseTheta->Theta_[1][p] = LooseTheta->Theta_[1][q] = C[q] * Stdev;

            LooseTheta->Theta_[2][p] = LooseTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

E0: if (Cinv) free(Cinv);

    if (C) free(C);

    if (N) free(N);

    if (Mode) free(Mode);

    return Error;
} // RoughEstimationPW

// Rough component parameter estimation for histogram.

int Rebmvnorm::RoughEstimationH(int                  k,           // Total number of bins.
                                FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl].
                                FLOAT                *h,          // Sides of the hypersquare.
                                FLOAT                nl,          // Total number of observations in class l.
                                int                  m,           // Mode index.
                                CompnentDistribution *RigidTheta, // Rigid parameters.
                                CompnentDistribution *LooseTheta) // Loose parameters.
{
    int                i, ii, j, l, o, p, q, r, *N = NULL;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, flm, flmin, flmax, V, Dlm, Dlmin, Sum, Stdev;
    FLOAT              *C = NULL, *Cinv = NULL, Cdet;
    int                Error = 0, Stop;

    // Global mode.

    Mode = (RoughParameterType*)malloc(length_pdf_ * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    N = (int*)malloc(length_pdf_ * sizeof(int));

    Error = NULL == N; if (Error) goto E0;

    V = (FLOAT)1.0;

    for (i = 0; i < length_pdf_; i++) {
        V *= h[i]; N[i] = 0;

        if (length_pdf_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                for (l = 0; l < length_pdf_; l++) if ((i != l) && (Y[j][l] != Y[m][l])) goto S0;

                Mode[i].klm += Y[j][length_pdf_]; X_[N[i]][i] = Y[j][i]; N[i] += 1;
S0:;
            }
        }
        else {
            Mode[i].klm = nl;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                X_[N[i]][i] = Y[j][i]; N[i] += 1;
            }
        }

        Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][length_pdf_] / (Mode[i].klm * h[i]);
    }

    flm = Y[m][length_pdf_] / (nl * V);

    // Variance-covariance matrix.

    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    Error = NULL == C; if (Error) goto E0;

    Cinv = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));

    Error = NULL == Cinv; if (Error) goto E0;

    if (nl >= length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            Sum = (FLOAT)0.0;
        
            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * (Y[j][i] - Mode[i].ym) * (Y[j][i] - Mode[i].ym);
            }

            C[i * length_pdf_ + i] = Sum / nl;

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;
        
                for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                    Sum += Y[j][length_pdf_] * (Y[j][i] - Mode[i].ym) * (Y[j][ii] - Mode[ii].ym);
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

        Stop = LUinvdet(length_pdf_, C, Cinv, &Cdet);

        if (Stop || (Cdet < (FLOAT)0.05)) {
            for (i = 0; i < length_pdf_; i++) {
                o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0; 

                for (ii = 0; ii < i; ii++) {
                    p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                    C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
                }
            }

            Cdet = (FLOAT)1.0;
        }
    }
    else {
        for (i = 0; i < length_pdf_; i++) {
            o = i * length_pdf_ + i; C[o] = Cinv[o] = (FLOAT)1.0; 

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                C[p] = C[q] = Cinv[p] = Cinv[q] = (FLOAT)0.0;
            }
        }

        Cdet = (FLOAT)1.0;
    }

    // Rigid restraints.

    RigidTheta->Theta_[3][0] = Cdet;

    for (i = 0; i < length_pdf_; i++) {
        RigidTheta->Theta_[0][i] = Mode[i].ym;

        Stdev = (FLOAT)1.0 / (SqrtPi2 * Mode[i].flm); Stdev *= Stdev;

        o = i * length_pdf_ + i;

        RigidTheta->Theta_[3][0] *= RigidTheta->Theta_[1][o] = Cinv[o] * Stdev; RigidTheta->Theta_[2][o] = (FLOAT)1.0 / Stdev;

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(RigidTheta->Theta_[1][o] * RigidTheta->Theta_[1][r]);

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] = C[q] * Stdev;

            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

    epsilon = (FLOAT)exp(-(FLOAT)2.0 * (LogSqrtPi2 + (FLOAT)log(flm) / length_pdf_) - (FLOAT)log(RigidTheta->Theta_[3][0]) / length_pdf_);

    if (epsilon > (FLOAT)1.0) {
        RigidTheta->Theta_[3][0] *= (FLOAT)exp(length_pdf_ * (FLOAT)log(epsilon));

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
    }

    Error = LooseTheta->Memmove(RigidTheta);

    if (Error) goto E0;

    if (Restraints_ == rtRigid) goto E0;

    // Loose restraints.

    for (i = 0; i < length_pdf_; i++) if (N[i] > 1) {
        // Bracketing.

        Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

        for (j = 0; j < N[i]; j++)  {
            Error = ComponentConditionalDist(i, X_[j], Cinv, LooseTheta, &CmpMrgDist);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * h[i];
        }

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)1.0 - (FLOAT)2.0 * p_value_; flmax = Mode[i].flm;

        // Bisection.

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            Stdev = (FLOAT)1.0 / (SqrtPi2 * flm); Stdev *= Stdev;

            o = i * length_pdf_ + i;

            LooseTheta->Theta_[1][o] = Cinv[o] * Stdev;

            Dlm = (FLOAT)1.0 - (FLOAT)2.0 * p_value_;

            for (j = 0; j < N[i]; j++) {
                Error = ComponentConditionalDist(i, X_[j], Cinv, LooseTheta, &CmpMrgDist);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * h[i];
            }

            if (((FLOAT)fabs(Dlm) < Eps) || (flmax - flmin < Eps)) {
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
E1:;
    }

    LooseTheta->Theta_[3][0] = Cdet;

    for (i = 0; i < length_pdf_; i++) {
        o = i * length_pdf_ + i;

        LooseTheta->Theta_[3][0] *= LooseTheta->Theta_[1][o]; LooseTheta->Theta_[2][o] = Cinv[o] / LooseTheta->Theta_[1][o];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i; r = ii * length_pdf_ + ii;

            Stdev = (FLOAT)sqrt(LooseTheta->Theta_[1][o] * LooseTheta->Theta_[1][r]);

            LooseTheta->Theta_[1][p] = LooseTheta->Theta_[1][q] = C[q] * Stdev;

            LooseTheta->Theta_[2][p] = LooseTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

E0: if (Cinv) free(Cinv);

    if (C) free(C);

    if (N) free(N);

    if (Mode) free(Mode);

    return Error;
} // RoughEstimationH

// Enhanced component parameter estimation for k-nearest neighbours.

int Rebmvnorm::EnhancedEstimationKNN(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,V,R].
                                     FLOAT                nl,          // Total number of observations in class l.
                                     CompnentDistribution *RigidTheta, // Rigid parameters.
                                     CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                Sum;
    int                  i, ii, j, o;
    int                  Error = 0;

    EnhanTheta = new CompnentDistribution(this);

    Error = NULL == EnhanTheta; if (Error) goto E0;

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    if (Error) goto E0;

    if (nl >= length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            EnhanTheta->pdf_[i] = pfNormal;

            Sum = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = Sum / nl;

            o = i * length_pdf_ + i;

            Sum = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][i] - EnhanTheta->Theta_[0][i]);
            }

            EnhanTheta->Theta_[1][o] = Sum / nl;

            if (EnhanTheta->Theta_[1][o] < RigidTheta->Theta_[1][o] * var_mul_) {
                Error = 1; goto E0;
            }

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;
        
                for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                    Sum += Y[j][length_pdf_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][ii] - EnhanTheta->Theta_[0][ii]);
                }

                EnhanTheta->Theta_[1][i * length_pdf_ + ii] = EnhanTheta->Theta_[1][ii * length_pdf_ + i] = Sum / nl;
            }
        }

        Error = LUinvdet(length_pdf_, EnhanTheta->Theta_[1], EnhanTheta->Theta_[2], EnhanTheta->Theta_[3]);

        if (Error) goto E0;

        Error = LooseTheta->Memmove(EnhanTheta);

        if (Error) goto E0;
    }
    else {
        Error = 1; goto E0;
    }

E0: if (EnhanTheta) delete EnhanTheta;

    return Error;
} // EnhancedEstimationKNN

// Enhanced component parameter estimation for Parzen window.

int Rebmvnorm::EnhancedEstimationPW(FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                    FLOAT                nl,          // Total number of observations in class l.
                                    CompnentDistribution *RigidTheta, // Rigid parameters.
                                    CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                Sum;
    int                  i, ii, j, o;
    int                  Error = 0;

    EnhanTheta = new CompnentDistribution(this);

    Error = NULL == EnhanTheta; if (Error) goto E0;

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    if (Error) goto E0;

    if (nl >= length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            EnhanTheta->pdf_[i] = pfNormal;

            Sum = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = Sum / nl;

            o = i * length_pdf_ + i;

            Sum = (FLOAT)0.0;

            for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][i] - EnhanTheta->Theta_[0][i]);
            }

            EnhanTheta->Theta_[1][o] = Sum / nl;

            if (EnhanTheta->Theta_[1][o] < RigidTheta->Theta_[1][o] * var_mul_) {
                Error = 1; goto E0;
            }

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;
        
                for (j = 0; j < n_; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                    Sum += Y[j][length_pdf_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][ii] - EnhanTheta->Theta_[0][ii]);
                }

                EnhanTheta->Theta_[1][i * length_pdf_ + ii] = EnhanTheta->Theta_[1][ii * length_pdf_ + i] = Sum / nl;
            }
        }

        Error = LUinvdet(length_pdf_, EnhanTheta->Theta_[1], EnhanTheta->Theta_[2], EnhanTheta->Theta_[3]);

        if (Error) goto E0;

        Error = LooseTheta->Memmove(EnhanTheta);

        if (Error) goto E0;
    }
    else {
        Error = 1; goto E0;
    }

E0: if (EnhanTheta) delete EnhanTheta;

    return Error;
} // EnhancedEstimationPW

// Enhanced component parameter estimation for histogram.

int Rebmvnorm::EnhancedEstimationH(int                  k,           // Total number of bins.
                                   FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                   FLOAT                nl,          // Total number of observations in class l.
                                   CompnentDistribution *RigidTheta, // Rigid parameters.
                                   CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta = NULL;
    FLOAT                Sum;
    int                  i, ii, j, o;
    int                  Error = 0;

    EnhanTheta = new CompnentDistribution(this);

    Error = NULL == EnhanTheta; if (Error) goto E0;

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    if (Error) goto E0;

    if (nl >= length_pdf_) {
        for (i = 0; i < length_pdf_; i++) {
            EnhanTheta->pdf_[i] = pfNormal;

            Sum = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * Y[j][i];
            }

            EnhanTheta->Theta_[0][i] = Sum / nl;

            o = i * length_pdf_ + i;

            Sum = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                Sum += Y[j][length_pdf_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][i] - EnhanTheta->Theta_[0][i]);
            }

            EnhanTheta->Theta_[1][o] = Sum / nl;

            if (EnhanTheta->Theta_[1][o] < RigidTheta->Theta_[1][o] * var_mul_) {
                Error = 1; goto E0;
            }

            for (ii = 0; ii < i; ii++) {
                Sum = (FLOAT)0.0;
        
                for (j = 0; j < k; j++) if (Y[j][length_pdf_] > FLOAT_MIN) {
                    Sum += Y[j][length_pdf_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][ii] - EnhanTheta->Theta_[0][ii]);
                }

                EnhanTheta->Theta_[1][i * length_pdf_ + ii] = EnhanTheta->Theta_[1][ii * length_pdf_ + i] = Sum / nl;
            }
        }

        Error = LUinvdet(length_pdf_, EnhanTheta->Theta_[1], EnhanTheta->Theta_[2], EnhanTheta->Theta_[3]);

        if (Error) goto E0;

        Error = LooseTheta->Memmove(EnhanTheta);

        if (Error) goto E0;
    }
    else {
        Error = 1; goto E0;
    }

E0: if (EnhanTheta) delete EnhanTheta;

    return Error;
} // EnhancedEstimationH

// Moments calculation.

int Rebmvnorm::MomentsCalculation(CompnentDistribution *CmpTheta, // Component parameters.
                                  FLOAT                *FirstM,   // First moment.
                                  FLOAT                *SecondM)  // Second moment.
{
    int i, ii, o, p, q;
    int Error = 0;

    for (i = 0; i < length_pdf_; i++) {
        FirstM[i] = CmpTheta->Theta_[0][i];

        o = i * length_pdf_ + i;

        SecondM[o] = CmpTheta->Theta_[1][o] + CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][i];

        for (ii = 0; ii < i; ii++) {
            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

            SecondM[p] = SecondM[q] = CmpTheta->Theta_[1][p] + CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][ii];
        }
    }

    return Error;
} // MomentsCalculation

// Bayes classification of the remaining observations for k-nearest neighbour.

int Rebmvnorm::BayesClassificationKNN(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                      int                  c,          // Number of components.
                                      FLOAT                *W,         // Component weights.
                                      CompnentDistribution **MixTheta, // Mixture parameters.
                                      FLOAT                **FirstM,   // First moments.
                                      FLOAT                **SecondM)  // Second moments.
{
    int   i, j, jj, l, o, p, q, outlier, Outlier = 0;
    FLOAT CmpDist, Max, Tmp, dW, N = (FLOAT)0.0;
    int   Error = 0;

    for (i = 0; i < n_; i++) {
        if (Y[i][length_pdf_] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(Y[i], MixTheta[l], &CmpDist, &outlier);

            if (Error) goto E0;

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, &outlier);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[i][length_pdf_];
            }
            else { 
                dW = Y[i][length_pdf_] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                    o = j * length_pdf_ + j;

                    SecondM[l][o] += dW * (Y[i][j] * Y[i][j] - SecondM[l][o]) / W[l];

                    for (jj = 0; jj < j; jj++) {
                        p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                        SecondM[l][p] = SecondM[l][q] += dW * (Y[i][j] * Y[i][jj] - SecondM[l][q]) / W[l];
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

        Error = LUinvdet(length_pdf_, MixTheta[i]->Theta_[1], MixTheta[i]->Theta_[2], MixTheta[i]->Theta_[3]);

        if (Error) goto E0;
    }

E0: return Error;
} // BayesClassificationKNN

// Bayes classification of the remaining observations for Parzen window.

int Rebmvnorm::BayesClassificationPW(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                     int                  c,          // Number of components.
                                     FLOAT                *W,         // Component weights.
                                     CompnentDistribution **MixTheta, // Mixture parameters.
                                     FLOAT                **FirstM,   // First moments.
                                     FLOAT                **SecondM)  // Second moments.
{
    int   i, j, jj, l, o, p, q, outlier, Outlier = 0;
    FLOAT CmpDist, Max, Tmp, dW, N = (FLOAT)0.0;
    int   Error = 0;

    for (i = 0; i < n_; i++) {
        if (Y[i][length_pdf_] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(Y[i], MixTheta[l], &CmpDist, &outlier);

            if (Error) goto E0;

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, &outlier);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier;
                }
            }

            if (Outlier) {
                N += Y[i][length_pdf_];
            }
            else { 
                dW = Y[i][length_pdf_] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                    o = j * length_pdf_ + j;

                    SecondM[l][o] += dW * (Y[i][j] * Y[i][j] - SecondM[l][o]) / W[l];

                    for (jj = 0; jj < j; jj++) {
                        p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                        SecondM[l][p] = SecondM[l][q] += dW * (Y[i][j] * Y[i][jj] - SecondM[l][q]) / W[l];
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

        Error = LUinvdet(length_pdf_, MixTheta[i]->Theta_[1], MixTheta[i]->Theta_[2], MixTheta[i]->Theta_[3]);

        if (Error) goto E0;
    }

E0: return Error;
} // BayesClassificationPW

// Bayes classification of the remaining observations for histogram.

int Rebmvnorm::BayesClassificationH(int                  k,          // Total number of bins.
                                    FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                    int                  c,          // Number of components.
                                    FLOAT                *W,         // Component weights.
                                    CompnentDistribution **MixTheta, // Mixture parameters.
                                    FLOAT                **FirstM,   // First moments.
                                    FLOAT                **SecondM)  // Second moments.
{
    int   i, j, jj, l, o, p, q, outlier, Outlier = 0;
    FLOAT CmpDist, Max, Tmp, dW, N = (FLOAT)0.0;
    int   Error = 0;

    for (i = 0; i < k; i++) {
        if (Y[i][length_pdf_] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(Y[i], MixTheta[l], &CmpDist, &outlier);

            if (Error) goto E0;

            Max = W[l] * CmpDist; Outlier = outlier;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist, &outlier);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp; Outlier = outlier; 
                }
            }

            if (Outlier) {
                N += Y[i][length_pdf_];
            }
            else { 
                dW = Y[i][length_pdf_] / n_; W[l] += dW;

                for (j = 0; j < length_pdf_; j++) {
                    FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                    o = j * length_pdf_ + j;

                    SecondM[l][o] += dW * (Y[i][j] * Y[i][j] - SecondM[l][o]) / W[l];

                    for (jj = 0; jj < j; jj++) {
                        p = j * length_pdf_ + jj; q = jj * length_pdf_ + j;

                        SecondM[l][p] = SecondM[l][q] += dW * (Y[i][j] * Y[i][jj] - SecondM[l][q]) / W[l];
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

        Error = LUinvdet(length_pdf_, MixTheta[i]->Theta_[1], MixTheta[i]->Theta_[2], MixTheta[i]->Theta_[3]);

        if (Error) goto E0;
    }

E0: return Error;
} // BayesClassificationH

// Returns component p.d.f..

int Rebmvnorm::ComponentDist(FLOAT                *Y,        // Pointer to the input point [y0,...,yd-1].
                             CompnentDistribution *CmpTheta, // Component parameters.
                             FLOAT                *CmpDist,  // Component distribution.
                             int                  *Outlier)  // 1 if outlier otherwise 0.
{
    FLOAT y, yi, yj;
    int   i, j;
    int   Error = 0;

    y = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        yi = Y[i] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + i] * yi * yi;

        for (j = i + 1; j < CmpTheta->length_pdf_; j++) {
            yj = Y[j] - CmpTheta->Theta_[0][j]; y += CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + j] * yi * yj;
        }
    }

    if (Outlier) {
        *Outlier = (FLOAT)2.0 * y > ChiSqr_;
    }

    *CmpDist = (FLOAT)exp(-y) / (FLOAT)sqrt((FLOAT)pow(Pi2, CmpTheta->length_pdf_) * CmpTheta->Theta_[3][0]);

    return Error;
} // ComponentDist

int Rebmvnorm::DegreesOffreedom(int c,                  // Number of components.
                                CompnentDistribution**, // Mixture parameters.
                                int *M)                 // Degrees of freedom.
{
    int i;
    int Error = 0;

    *M = c - 1;

    for (i = 0; i < c; i++) {
        *M += length_pdf_ + (length_pdf_ + 1) * length_pdf_ / 2;
    }

    return Error;
} // DegreesOffreedom
