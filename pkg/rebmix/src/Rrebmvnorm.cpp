#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "rngmvnormf.h"
#include "rebmvnormf.h"

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

extern "C" {

// Runs RRNGMVNORM in R.

void RRNGMVNORM(int    *IDum,         // Random seed.
                int    *d,            // Number of independent random variables.
                int    *c,            // Number of components.
                int    *N,            // Numbers of observations.
                int    *length_pdf,   // Length of pdf.
                int    *length_Theta, // Length of Theta.
                int    *length_theta, // Length of Theta[i].
                char   **pdf,         // Parametric family types.
                double *Theta,        // Component parameters.
                int    *n,            // Number of observations.
                double *Y,            // Dataset.
                int    *Error)        // Error code.
{
    Rngmvnorm *rngmvnorm;
    int       i, j, k, l;

    rngmvnorm = new Rngmvnorm;

    *Error = NULL == rngmvnorm; if (*Error) goto E0;

    rngmvnorm->IDum_ = *IDum;
    rngmvnorm->d_ = *d;
    rngmvnorm->c_ = *c;

    rngmvnorm->N_ = (int*)malloc(rngmvnorm->c_ * sizeof(int));

    *Error = NULL == rngmvnorm->N_; if (*Error) goto E0;

    for (i = 0; i < rngmvnorm->c_; i++) rngmvnorm->N_[i] = N[i];

    rngmvnorm->length_pdf_ = *length_pdf;

    rngmvnorm->length_Theta_ = *length_Theta;

    rngmvnorm->length_theta_ = (int*)malloc(rngmvnorm->length_Theta_ * sizeof(int));

    *Error = NULL == rngmvnorm->length_theta_; if (*Error) goto E0;

    for (i = 0; i < rngmvnorm->length_Theta_; i++) {
        rngmvnorm->length_theta_[i] = length_theta[i];
    }

    rngmvnorm->IniTheta_ = new CompnentDistribution(rngmvnorm);

    *Error = NULL == rngmvnorm->IniTheta_; if (*Error) goto E0;

    *Error = rngmvnorm->IniTheta_->Realloc(rngmvnorm->length_pdf_, rngmvnorm->length_Theta_, rngmvnorm->length_theta_);

    if (*Error) goto E0;

    for (i = 0; i < rngmvnorm->length_pdf_; i++) {
        rngmvnorm->IniTheta_->pdf_[i] = pfNormal;
    }

    rngmvnorm->MixTheta_ = new CompnentDistribution*[rngmvnorm->c_];

    *Error = NULL == rngmvnorm->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < rngmvnorm->c_; i++) {
        rngmvnorm->MixTheta_[i] = new CompnentDistribution(rngmvnorm);

        *Error = NULL == rngmvnorm->MixTheta_[i]; if (*Error) goto E0;

        *Error = rngmvnorm->MixTheta_[i]->Realloc(rngmvnorm->length_pdf_, rngmvnorm->length_Theta_, rngmvnorm->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < rngmvnorm->c_; i++) {
        for (j = 0; j < rngmvnorm->length_pdf_; j++) {
            rngmvnorm->MixTheta_[i]->pdf_[j] = rngmvnorm->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rngmvnorm->length_Theta_; j++) if (rngmvnorm->IniTheta_->Theta_[j]) {
        for (k = 0; k < rngmvnorm->c_; k++) {
            for (l = 0; l < rngmvnorm->length_theta_[j]; l++) {
                rngmvnorm->MixTheta_[k]->Theta_[j][l] = Theta[i];

                i++;
            }
        }
    }

    *Error = rngmvnorm->RNGMIX();

    if (*Error) goto E0;

    *n = rngmvnorm->n_; i = 0;

    for (j = 0; j < rngmvnorm->d_; j++) {
        for (k = 0; k < rngmvnorm->n_; k++) {
            Y[i] = rngmvnorm->Y_[k][j]; i++;
        }
    }

E0: if (rngmvnorm) delete rngmvnorm;
} // RRNGMVNORM

// Runs RREBMVNORM in R.

void RREBMVNORM(char   **Preprocessing, // Preprocessing type.
                int    *cmax,           // Maximum number of components.
                char   **Criterion,     // Infromation criterion type.
                int    *d,              // Number of independent random variables.
                char   **Variables,     // Types of variables.
                int    *length_pdf,     // Length of pdf.
                char   **pdf,           // Parametric family types.
                int    *length_Theta,   // Length of Theta.
                int    *length_theta,   // Length of Theta[i].
                double *Theta,          // Component parameters.
                int    *length_K,       // Length of K.
                int    *K,              // Numbers of bins v or numbers of nearest neighbours k.
                int    *length_y0,      // Length of y0.
                double *y0,             // Origins.
                int    *length_ymin,    // Length of ymin.
                double *ymin,           // Minimum observations.
                int    *length_ymax,    // Length of ymax.
                double *ymax,           // Maximum observations.
                double *ar,             // Acceleration rate.
                char   **Restraints,    // Restraints type.
                int    *n,              // Number of observations.
                double *Y,              // Dataset.
                int    *summary_k,      // Optimal v or optimal k.
                double *summary_h,      // Optimal class widths of length d.
                double *summary_y0,     // Optimal origins of length d.
                double *summary_IC,     // Optimal information criterion.
                double *summary_logL,   // Log-likelihood.
                int    *summary_M,      // Degrees of freedom.
                int    *summary_c,      // Optimal number of components.
                double *W,              // Component weights.
                double *theta1,         // Component parameters.
                double *theta2,         // Component parameters.
                int    *opt_length,     // Length of opt_c, opt_IC, opt_logL and opt_D.
                int    *opt_c,          // Numbers of components for optimal v or for optimal k.
                double *opt_IC,         // Information criteria for optimal v or for optimal k.
                double *opt_logL,       // Log-likelihoods for optimal v or for optimal k.
                double *opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
                int    *all_length,     // Length of all_K and all_IC.
                int    *all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                double *all_IC,         // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
                int    *Error)          // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    int       i, j, k, l;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    if (!strcmp(Preprocessing[0], "histogram")) {
        rebmvnorm->Preprocessing_ = poHistogram;
    }
    else
    if (!strcmp(Preprocessing[0], "Parzen window")) {
        rebmvnorm->Preprocessing_ = poParzenWindow;
    }
    else
    if (!strcmp(Preprocessing[0], "k-nearest neighbour")) {
        rebmvnorm->Preprocessing_ = poKNearestNeighbour;
    }
    else {
        *Error = 1; goto E0;
    }

    rebmvnorm->cmax_ = *cmax;

    if (!strcmp(Criterion[0], "AIC"))
        rebmvnorm->Criterion_ = icAIC; 
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmvnorm->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmvnorm->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmvnorm->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmvnorm->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmvnorm->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmvnorm->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmvnorm->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmvnorm->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmvnorm->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmvnorm->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmvnorm->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmvnorm->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmvnorm->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmvnorm->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmvnorm->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmvnorm->d_ = *d;

    rebmvnorm->Variables_ = (VariablesType_e*)malloc(rebmvnorm->d_ * sizeof(VariablesType_e));

    *Error = NULL == rebmvnorm->Variables_; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->d_; i++) {
        if (!strcmp(Variables[i], "continuous")) {
            rebmvnorm->Variables_[i] = vtContinuous;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    *Error = NULL == rebmvnorm->IniTheta_; if (*Error) goto E0;

    rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (int*)malloc(rebmvnorm->length_Theta_ * sizeof(int));

    *Error = NULL == rebmvnorm->length_theta_; if (*Error) goto E0;

    *Error = rebmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmvnorm->IniTheta_->pdf_[i] = pfNormal;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_Theta_; j++) if (rebmvnorm->IniTheta_->Theta_[j]) {
        for (k = 0; k < rebmvnorm->length_theta_[j]; k++) {
            rebmvnorm->IniTheta_->Theta_[j][k] = Theta[i];

            i++;
        }
    }

    rebmvnorm->length_K_ = *length_K;

    rebmvnorm->K_ = (int*)malloc(rebmvnorm->length_K_ * sizeof(int));

    *Error = NULL == rebmvnorm->K_; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_K_; i++) {
        rebmvnorm->K_[i] = K[i];
    }

    if (*length_y0 > 0) {
        rebmvnorm->y0_ = (FLOAT*)malloc(rebmvnorm->d_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->y0_; if (*Error) goto E0;

        for (i = 0; i < rebmvnorm->d_; i++) {
            rebmvnorm->y0_[i] = y0[i];
        }
    }
    else {
        rebmvnorm->y0_ = NULL;
    }

    if (*length_ymin > 0) {
        rebmvnorm->ymin_ = (FLOAT*)malloc(rebmvnorm->d_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->ymin_; if (*Error) goto E0;

        for (i = 0; i < rebmvnorm->d_; i++) {
            rebmvnorm->ymin_[i] = ymin[i];
        }
    }
    else {
        rebmvnorm->ymin_ = NULL;
    }

    if (*length_ymax > 0) {
        rebmvnorm->ymax_ = (FLOAT*)malloc(rebmvnorm->d_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->ymax_; if (*Error) goto E0;

        for (i = 0; i < rebmvnorm->d_; i++) {
            rebmvnorm->ymax_[i] = ymax[i];
        }
    }
    else {
        rebmvnorm->ymax_ = NULL;
    }

    rebmvnorm->ar_ = *ar;

    if (!strcmp(Restraints[0], "rigid")) {
        rebmvnorm->Restraints_ = rtRigid;
    }
    else
    if (!strcmp(Restraints[0], "loose")) {
        rebmvnorm->Restraints_ = rtLoose;
    }
    else {
        *Error = 1; goto E0;
    }

    rebmvnorm->n_ = *n;

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->n_ * sizeof(FLOAT*));

    *Error = NULL == rebmvnorm->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->n_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->d_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmvnorm->d_; j++) {
        for (l = 0; l < rebmvnorm->n_; l++) {
            rebmvnorm->Y_[l][j] = Y[i]; i++;
        }
    }

    *Error = rebmvnorm->REBMIX();

    if (*Error) goto E0;

    *summary_k = rebmvnorm->summary_.k;

    if (rebmvnorm->summary_.h) for (i = 0; i < rebmvnorm->d_; i++) {
        summary_h[i] = rebmvnorm->summary_.h[i];
    }

    if (rebmvnorm->summary_.y0) for (i = 0; i < rebmvnorm->d_; i++) {
        summary_y0[i] = rebmvnorm->summary_.y0[i];
    }

    *summary_IC = rebmvnorm->summary_.IC;
    *summary_logL = rebmvnorm->summary_.logL;
    *summary_M = rebmvnorm->summary_.M;
    *summary_c = rebmvnorm->summary_.c;

    for (j = 0; j < rebmvnorm->summary_.c; j++) {
        W[j] = rebmvnorm->W_[j];
    }

    i = 0;

    for (j = 0; j < rebmvnorm->summary_.c; j++) {
        for (l = 0; l < rebmvnorm->length_theta_[0]; l++) {
            theta1[i] = rebmvnorm->MixTheta_[j]->Theta_[0][l];

            i++;
        }
    }

    i = 0;

    for (j = 0; j < rebmvnorm->summary_.c; j++) {
        for (l = 0; l < rebmvnorm->length_theta_[1]; l++) {
            theta2[i] = rebmvnorm->MixTheta_[j]->Theta_[1][l];

            i++;
        }
    }

    *opt_length = rebmvnorm->opt_length_;

    for (i = 0; i < rebmvnorm->opt_length_; i++) {
        opt_c[i] = rebmvnorm->opt_c_[i];
        opt_IC[i] = rebmvnorm->opt_IC_[i];
        opt_logL[i] = rebmvnorm->opt_logL_[i];
        opt_D[i] = rebmvnorm->opt_D_[i];
    }

    i = 0;

    for (j = 0; j < rebmvnorm->all_length_; j++) if (rebmvnorm->all_K_[j]) {
        all_K[i] = rebmvnorm->all_K_[j]; all_IC[i] = rebmvnorm->all_IC_[j]; 

        i++;
    }

    *all_length = i;

E0: if (rebmvnorm) delete rebmvnorm;
} // RREBMVNORM

// Returns classified observations in R.

void RCLSMVNORM(int    *n,      // Total number of independent observations.
                double *X,      // Pointer to the input array X.
                int    *s,      // Number of classes.
                int    *o,      // Number of input REBMVNORM objects.
                int    *d,      // Number of independent random variables in REBMVNORM objects.
                int    *c,      // Number of components.
                double *W,      // Component weights.
                char   **pdf,   // Component parameters.
                double *theta1, // Component parameters.
                double *theta2, // Component parameters.
                double *P,      // Prior probabilities.
                int    *Z,      // Pointer to the output array Z.
                int    *Error)  // Error code.
{
    Rebmvnorm            *rebmvnorm = NULL;
    int                  **C = NULL; 
    int                  A[2];
    FLOAT                ***Q = NULL;
    FLOAT                *Y = NULL;
    CompnentDistribution ****Theta = NULL; 
    FLOAT                CmpDist, MixDist, MaxMixDist;
    int                  i, j, k, l, m;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

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

    Theta = new CompnentDistribution*** [*s];

    *Error = NULL == Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        Theta[j] = new CompnentDistribution** [*o];

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            Theta[j][k] = new CompnentDistribution* [C[j][k]];

            *Error = NULL == Theta[j][k]; if (*Error) goto E0;

            for (l = 0; l < C[j][k]; l++) {
                Theta[j][k][l] = new CompnentDistribution(rebmvnorm);

                *Error = NULL == Theta[j][k][l]; if (*Error) goto E0;

                A[0] = A[1] = d[k];

                *Error = Theta[j][k][l]->Realloc(d[k], 2, A);

                if (*Error) goto E0;

                for (m = 0; m < d[k]; m++) {
                    if (!strcmp(pdf[i], "normal")) {
                        Theta[j][k][l]->pdf_[m] = pfNormal;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
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

                *Error = rebmvnorm->MixtureDist(Y, C[j][l], Q[j][l], Theta[j][l], &CmpDist);

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
                            if (Theta[i][j][k]) delete Theta[i][j][k];
                        }
                        
                        delete Theta[i][j];
                    }
                }

                delete Theta[i];
            }
        }

        delete Theta;
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

    if (rebmvnorm) delete rebmvnorm;
} // RCLSMVNORM 

} // extern "C"