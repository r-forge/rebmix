#include <math.h>
#include <float.h> 

#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "base.h"
#include "rebmvnormf.h"

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

    o = i * d_ + i;

    Stdev = (FLOAT)sqrt(CmpTheta->Theta_[1][o] / Cinv[o]);

    y = (Y[i] - Mean) / (Sqrt2 * Stdev);

    *CmpMrgDist = (FLOAT)exp(-(y * y)) / (Sqrt2Pi * Stdev);

    return Error;
} // ComponentConditionalDist 

// Rough component parameter estimation for histogram.

int Rebmvnorm::RoughEstimationH(int                  k,           // Total number of bins.
                                FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl].
                                FLOAT                *h,          // Sides of the hypersquare.
                                FLOAT                nl,          // Total number of observations in class l.
                                int                  m,           // Mode index.
                                CompnentDistribution *RigidTheta, // Rigid parameters.
                                CompnentDistribution *LooseTheta) // Loose parameters.
{
    int                i, ii, j, l, o, p, q, r;
    RoughParameterType *Mode = NULL;
    FLOAT              CmpMrgDist, epsilon, flm, flmin, flmax, V, Dlm, Dlmin, Sum, Stdev;
    FLOAT              *C = NULL, *Cinv = NULL, Cdet;
    int                Error = 0, Stop;

    // Global mode.

    Mode = (RoughParameterType*)malloc(d_ * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    flm = (FLOAT)1.0; V = (FLOAT)1.0;

    for (i = 0; i < d_; i++) {
        V *= h[i];

        if (d_ > 1) {
            Mode[i].klm = (FLOAT)0.0;

            for (j = 0; j < k; j++) {
                for (l = 0; l < d_; l++) if ((i != l) && (Y[j][l] != Y[m][l])) goto S0;

                Mode[i].klm += Y[j][d_];
S0:;
            }
        }
        else
            Mode[i].klm = nl;

        Mode[i].ym = Y[m][i]; Mode[i].flm = Y[m][d_] / (Mode[i].klm * h[i]);
    }

    flm = Y[m][d_] / (nl * V);

    // Variance-covariance matrix.

    C = (FLOAT*)malloc(d_ * d_ * sizeof(FLOAT));

    Error = NULL == C; if (Error) goto E0;

    for (i = 0; i < d_; i++) {
        Sum = (FLOAT)0.0;
        
        for (j = 0; j < k; j++) if (Y[j][d_] > FLOAT_MIN) {
            Sum += Y[j][d_] * (Y[j][i] - Mode[i].ym) * (Y[j][i] - Mode[i].ym);
        }

        C[i * d_ + i] = Sum / nl; 

        for (ii = 0; ii < i; ii++) {
            Sum = (FLOAT)0.0;
        
            for (j = 0; j < k; j++) if (Y[j][d_] > FLOAT_MIN) {
                Sum += Y[j][d_] * (Y[j][i] - Mode[i].ym) * (Y[j][ii] - Mode[ii].ym);
            }

            C[i * d_ + ii] = C[ii * d_ + i] = Sum / nl;
        }
    }

    // Correlation matrix.

    for (i = 0; i < d_; i++) {
        o = i * d_ + i; 

        for (ii = 0; ii < i; ii++) {
            p = i * d_ + ii; q = ii * d_ + i; r = ii * d_ + ii;

            if (((FLOAT)fabs(C[q]) < FLOAT_MIN) || (C[o] < FLOAT_MIN) || (C[r] < FLOAT_MIN)) {
                C[p] = C[q] = (FLOAT)0.0;
            }
            else {
                C[p] = C[q] /= (FLOAT)sqrt(C[o] * C[r]);
            }
        }
    }

    for (i = 0; i < d_; i++) {
        C[i * d_ + i] = (FLOAT)1.0; 
    }

    Cinv = (FLOAT*)malloc(d_ * d_ * sizeof(FLOAT));

    Error = NULL == Cinv; if (Error) goto E0;

    Error = LUinvdet(d_, C, Cinv, &Cdet);

    if (Error) goto E0;

    // Rigid restraints.

    RigidTheta->Theta_[3][0] = Cdet;

    for (i = 0; i < d_; i++) {
        RigidTheta->Theta_[0][i] = Mode[i].ym;

        Stdev = (FLOAT)1.0 / (Sqrt2Pi * Mode[i].flm); Stdev *= Stdev;

        o = i * d_ + i;

        RigidTheta->Theta_[3][0] *= RigidTheta->Theta_[1][o] = Cinv[o] * Stdev; RigidTheta->Theta_[2][o] = Cinv[o] / Stdev;

        for (ii = 0; ii < i; ii++) {
            p = i * d_ + ii; q = ii * d_ + i; r = ii * d_ + ii;

            Stdev = (FLOAT)sqrt(RigidTheta->Theta_[1][o] * RigidTheta->Theta_[1][r]);

            RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] = C[q] * Stdev;

            RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

    epsilon = (FLOAT)exp(-(FLOAT)2.0 * (LogSqrt2Pi + (FLOAT)log(flm) / d_) - (FLOAT)log(RigidTheta->Theta_[3][0]) / d_);

    if (epsilon > (FLOAT)1.0) {
        RigidTheta->Theta_[3][0] *= (FLOAT)exp(d_ * (FLOAT)log(epsilon));

        for (i = 0; i < d_; i++) {
            o = i * d_ + i;

            RigidTheta->Theta_[1][o] *= epsilon;
            RigidTheta->Theta_[2][o] /= epsilon;

            for (ii = 0; ii < i; ii++) {
                p = i * d_ + ii; q = ii * d_ + i;

                RigidTheta->Theta_[1][p] = RigidTheta->Theta_[1][q] *= epsilon;
                RigidTheta->Theta_[2][p] = RigidTheta->Theta_[2][q] /= epsilon;
            }
        }
    }

    Error = LooseTheta->Memmove(RigidTheta);

    if (Error) goto E0;

    if (Restraints_ == rtRigid) goto E0;

    // Loose restraints.

    for (i = 0; i < d_; i++) {
        // Bracketing.

        Dlm = (FLOAT)0.0;

        for (o = 0; o < k; o++) if (Y[o][d_] > FLOAT_MIN) {
            for (p = 0; p < d_; p++) if ((i != p) && (Y[o][p] != Y[m][p])) goto S1;

            Error = ComponentConditionalDist(i, Y[o], Cinv, LooseTheta, &CmpMrgDist);

            if (Error) goto E0;

            Dlm -= CmpMrgDist * h[i];
S1:;
        }

        Dlm += (FLOAT)0.998;

        if (Dlm > (FLOAT)0.0) goto E1;

        flmin = (FLOAT)0.0; Dlmin = (FLOAT)0.998; flmax = Mode[i].flm;

        // Bisection.

        Stop = 0;

        while (!Stop) {
            flm = (flmax + flmin) / (FLOAT)2.0;

            Stdev = (FLOAT)1.0 / (Sqrt2Pi * flm); Stdev *= Stdev;

            o = i * d_ + i;

            LooseTheta->Theta_[1][o] = Cinv[o] * Stdev;

            Dlm = (FLOAT)0.0;

            for (o = 0; o < k; o++) if (Y[o][d_] > FLOAT_MIN) {
                for (p = 0; p < d_; p++) if ((i != p) && (Y[o][p] != Y[m][p])) goto S2;

                Error = ComponentConditionalDist(i, Y[o], Cinv, LooseTheta, &CmpMrgDist);

                if (Error) goto E0;

                Dlm -= CmpMrgDist * h[i];
S2:;
            }

            Dlm += (FLOAT)0.998;

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

    for (i = 0; i < d_; i++) {
        o = i * d_ + i;

        LooseTheta->Theta_[3][0] *= LooseTheta->Theta_[1][o]; LooseTheta->Theta_[2][o] = Cinv[o] / LooseTheta->Theta_[1][o];

        for (ii = 0; ii < i; ii++) {
            p = i * d_ + ii; q = ii * d_ + i; r = ii * d_ + ii;

            Stdev = (FLOAT)sqrt(LooseTheta->Theta_[1][o] * LooseTheta->Theta_[1][r]);

            LooseTheta->Theta_[1][p] = LooseTheta->Theta_[1][q] = C[q] * Stdev;

            LooseTheta->Theta_[2][p] = LooseTheta->Theta_[2][q] = Cinv[q] / Stdev;
        }
    }

E0: if (Cinv) free(Cinv);

    if (C) free(C);

    if (Mode) free(Mode);

    return Error;
} // RoughEstimationH

// Enhanced component parameter estimation for histogram.

int Rebmvnorm::EnhancedEstimationH(int                  k,           // Total number of bins.
                                   FLOAT                **Y,         // Pointer to the input points [y0,...,yd-1,kl,k].
                                   FLOAT                nl,          // Total number of observations in class l.
                                   CompnentDistribution *RigidTheta, // Rigid parameters.
                                   CompnentDistribution *LooseTheta) // Loose parameters.
{
    CompnentDistribution *EnhanTheta;
    FLOAT                Sum;
    int                  i, ii, j, o;
    int                  Error = 0;

    EnhanTheta = new CompnentDistribution(this);

    Error = NULL == EnhanTheta; if (Error) goto E0;

    Error = EnhanTheta->Realloc(length_pdf_, length_Theta_, length_theta_);

    if (Error) goto E0;

    for (i = 0; i < d_; i++) {
        EnhanTheta->pdf_[i] = pfNormal;

        Sum = (FLOAT)0.0;

        for (j = 0; j < k; j++) if (Y[j][d_] > FLOAT_MIN) {
            Sum += Y[j][d_] * Y[j][i];
        }

        EnhanTheta->Theta_[0][i] = Sum / nl;

        o = i * d_ + i;

        Sum = (FLOAT)0.0;

        for (j = 0; j < k; j++) if (Y[j][d_] > FLOAT_MIN) {
            Sum += Y[j][d_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][i] - EnhanTheta->Theta_[0][i]);
        }

        EnhanTheta->Theta_[1][o] = Sum / nl;

        if (EnhanTheta->Theta_[1][o] <= FLOAT_MIN) {
            Error = 1; if (Error) goto E0;
        }

        if (EnhanTheta->Theta_[1][o] < RigidTheta->Theta_[1][o]) {
            Error = 1; if (Error) goto E0;
        }

        for (ii = 0; ii < i; ii++) {
            Sum = (FLOAT)0.0;
        
            for (j = 0; j < k; j++) if (Y[j][d_] > FLOAT_MIN) {
                Sum += Y[j][d_] * (Y[j][i] - EnhanTheta->Theta_[0][i]) * (Y[j][ii] - EnhanTheta->Theta_[0][ii]);
            }

            EnhanTheta->Theta_[1][i * d_ + ii] = EnhanTheta->Theta_[1][ii * d_ + i] = Sum / nl;
        }
    }

    Error = LUinvdet(d_, EnhanTheta->Theta_[1], EnhanTheta->Theta_[2], EnhanTheta->Theta_[3]);

    if (Error) goto E0;

    Error = LooseTheta->Memmove(EnhanTheta);

    if (Error) goto E0;

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

    for (i = 0; i < d_; i++) {
        FirstM[i] = CmpTheta->Theta_[0][i];

        o = i * d_ + i;

        SecondM[o] = CmpTheta->Theta_[1][o] + CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][i];

        for (ii = 0; ii < i; ii++) {
            p = i * d_ + ii; q = ii * d_ + i;

            SecondM[p] = SecondM[q] = CmpTheta->Theta_[1][p] + CmpTheta->Theta_[0][i] * CmpTheta->Theta_[0][ii];
        }
    }

    return Error;
} // MomentsCalculation

// Bayes classification of the remaining observations for histogram.

int Rebmvnorm::BayesClassificationH(int                  k,          // Total number of bins.
                                    FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                                    int                  c,          // Number of components.
                                    FLOAT                *W,         // Component weights.
                                    CompnentDistribution **MixTheta, // Mixture parameters.
                                    FLOAT                **FirstM,   // First moments.
                                    FLOAT                **SecondM)  // Second moments.
{
    int   i, j, jj, l, o, p, q;
    FLOAT CmpDist, Max, Tmp, dW;
    int   Error = 0;

    for (i = 0; i < k; i++) {
        if (Y[i][d_] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(Y[i], MixTheta[l], &CmpDist);

            if (Error) goto E0;

            Max = W[l] * CmpDist;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(Y[i], MixTheta[j], &CmpDist);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) {
                    l = j; Max = Tmp;
                }
            }

            dW = Y[i][d_] / n_; W[l] += dW;

            for (j = 0; j < d_; j++) {
                FirstM[l][j] += dW * (Y[i][j] - FirstM[l][j]) / W[l];

                o = j * d_ + j;

                SecondM[l][o] += dW * (Y[i][j] * Y[i][j] - SecondM[l][o]) / W[l];

                for (jj = 0; jj < j; jj++) {
                    p = j * d_ + jj; q = jj * d_ + j;

                    SecondM[l][p] = SecondM[l][q] += dW * (Y[i][j] * Y[i][jj] - SecondM[l][p]) / W[l];
                }
            }
        }
    }

    for (i = 0; i < c; i++) for (j = 0; j < d_; j++) {
        MixTheta[i]->Theta_[0][j] = FirstM[i][j];

        o = j * d_ + j;

        MixTheta[i]->Theta_[1][o] = SecondM[i][o] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][j];

        for (jj = 0; jj < j; jj++) {
            p = j * d_ + jj; q = jj * d_ + j;

            MixTheta[i]->Theta_[1][p] = MixTheta[i]->Theta_[1][q] = SecondM[i][p] - MixTheta[i]->Theta_[0][j] * MixTheta[i]->Theta_[0][jj];
        }
    }

E0: return Error;
} // BayesClassificationH

// Returns component p.d.f..

int Rebmvnorm::ComponentDist(FLOAT                *Y,        // Pointer to the input point [y0,...,yd-1].
                             CompnentDistribution *CmpTheta, // Component parameters.
                             FLOAT                *CmpDist)  // Component distribution.
{
    FLOAT y, yi, yj;
    int   i, j;
    int   Error = 0;

    y = (FLOAT)0.0;

    for (i = 0; i < d_; i++) {
        yi = Y[i] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * d_ + i] * yi * yi;

        for (j = i + 1; j < d_; j++) {
            yj = Y[j] - CmpTheta->Theta_[0][j]; y += CmpTheta->Theta_[2][i * d_ + j] * yi * yj;
        }
    }

    *CmpDist = (FLOAT)exp(-y) / (FLOAT)sqrt((FLOAT)pow((FLOAT)2.0 * Pi, d_) * CmpTheta->Theta_[3][0]);

    return Error;
} // ComponentDist

int Rebmvnorm::DegreesOffreedom(int                  c,          // Number of components.
                                CompnentDistribution **MixTheta, // Mixture parameters.
                                int                  *M)         // Degrees of freedom.
{
    int Error = 0;

    *M = c * (d_ + (d_ + 1) * d_ / 2) + c - 1;

    return Error;
} // DegreesOffreedom
