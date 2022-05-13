#include "rngmvnormf.h"
#include "rebmvnormf.h"

#include <math.h>

extern "C" {

// Runs RRNGMVNORM in R.

void RRNGMVNORM(INT    *IDum,         // Random seed.
                INT    *d,            // Number of independent random variables.
                INT    *c,            // Number of components.
                INT    *N,            // Numbers of observations.
                INT    *length_pdf,   // Length of pdf.
                INT    *length_Theta, // Length of Theta.
                INT    *length_theta, // Length of Theta[i].
                double *Theta,        // Component parameters.
                INT    *n,            // Number of observations.
                double *Y,            // Dataset.
                INT    *Z,            // Component membership.
                INT    *Error)        // Error code.
{
    Rngmvnorm *rngmvnorm = NULL;
    INT       i, j, k, l;

    rngmvnorm = new Rngmvnorm;

    *Error = NULL == rngmvnorm; if (*Error) goto E0;

    rngmvnorm->IDum_ = *IDum;
    rngmvnorm->length_pdf_ = *d;
    rngmvnorm->c_ = *c;

    rngmvnorm->N_ = (INT*)malloc(rngmvnorm->c_ * sizeof(INT));

    *Error = NULL == rngmvnorm->N_; if (*Error) goto E0;

    for (i = 0; i < rngmvnorm->c_; i++) rngmvnorm->N_[i] = N[i];

    rngmvnorm->IniTheta_ = new CompnentDistribution(rngmvnorm);

    *Error = NULL == rngmvnorm->IniTheta_; if (*Error) goto E0;

    rngmvnorm->length_pdf_ = *length_pdf;

    rngmvnorm->length_Theta_ = *length_Theta;

    rngmvnorm->length_theta_ = (INT*)malloc(rngmvnorm->length_Theta_ * sizeof(INT));

    *Error = NULL == rngmvnorm->length_theta_; if (*Error) goto E0;

    *Error = rngmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rngmvnorm->length_pdf_; i++) {
        rngmvnorm->IniTheta_->pdf_[i] = pfNormal;
    }

    rngmvnorm->MixTheta_ = new CompnentDistribution*[(unsigned INT)rngmvnorm->c_];

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

    for (j = 0; j < rngmvnorm->length_pdf_; j++) {
        for (k = 0; k < rngmvnorm->n_; k++) {
            Y[i] = rngmvnorm->Y_[j][k]; i++;
        }
    }

    for (i = 0; i < rngmvnorm->n_; i++) {
        Z[i] = rngmvnorm->Z_[i];
    }

E0: if (rngmvnorm) delete rngmvnorm;
} // RRNGMVNORM

// Runs RREBMVNORM in R.

void RREBMVNORM(char   **Preprocessing, // Preprocessing type.
                INT    *cmax,           // Maximum number of components.
                INT    *cmin,           // Minimum number of components.
                char   **Criterion,     // Information criterion type.
                INT    *d,              // Number of independent random variables.
                char   **Variables,     // Types of variables.
                INT    *length_pdf,     // Length of pdf.
                char   **pdf,           // Parametric family types.
                INT    *length_Theta,   // Length of Theta.
                INT    *length_theta,   // Length of Theta[i].
                double *Theta,          // Component parameters.
                INT    *length_K,       // Length of K.
                INT    *K,              // Numbers of bins v or numbers of nearest neighbours k.
                INT    *length_ymin,    // Length of ymin.
                double *ymin,           // Minimum observations.
                INT    *length_ymax,    // Length of ymax.
                double *ymax,           // Maximum observations.
                INT    *length_h,       // Length of h.
                double *h,              // Sides of the hypersquare.
                double *ar,             // Acceleration rate.
                char   **Restraints,    // Restraints type.
                INT    *n,              // Number of observations.
                double *Y,              // Dataset.
                INT    *Y_type,         // Dataset type.
/// Panic Branislav
                char   **EMStrategy,       // Strategy for EM algorithm.
                char   **EMVariant,        // EM algorithm variant.
                char   **EMAcceleration,   // Acceleration for the standard EM algorithm.
                double *EMTolerance,       // Tolerance for EM algortihm.
                double *EMAccelerationMul, // Acceleration rate for Em algorithm.
                INT    *EMMaxIter,         // Maximum number of iterations in EM algorithm.
                INT    *EMK,               // Number of bins for histogram EM algorithm.
                INT    *n_iter,            // Number of iterations for optimal case.
                INT    *n_iter_sum,        // Number of iterations in whole run.
/// End
                INT    *summary_k,      // Optimal v or optimal k.
                double *summary_h,      // Optimal class widths of length d.
                double *summary_y0,     // Optimal origins of length d.
                double *summary_ymin,   // Optimal minimum observations of length d.
                double *summary_ymax,   // Optimal maximum observations of length d.
                double *summary_IC,     // Optimal information criterion.
                double *summary_logL,   // Log-likelihood.
                INT    *summary_M,      // Degrees of freedom.
                INT    *summary_c,      // Optimal number of components.
                double *W,              // Component weights.
                double *theta1,         // Component parameters.
                double *theta2,         // Component parameters.
                INT    *opt_length,     // Length of opt_c, opt_IC, opt_logL and opt_D.
                INT    *opt_c,          // Numbers of components for optimal v or for optimal k.
                double *opt_IC,         // Information criteria for optimal v or for optimal k.
                double *opt_logL,       // Log-likelihoods for optimal v or for optimal k.
                double *opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
                INT    *all_length,     // Length of all_K and all_IC.
                INT    *all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                double *all_IC,         // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
                INT    *Error)          // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    *Error = rebmvnorm->Set(Preprocessing,     // Preprocessing type.
                            cmax,              // Maximum number of components.
                            cmin,              // Minimum number of components.
                            Criterion,         // Information criterion type.
                            d,                 // Number of independent random variables.
                            Variables,         // Types of variables.
                            length_pdf,        // Length of pdf.
                            pdf,               // Parametric family types.
                            length_Theta,      // Length of Theta.
                            length_theta,      // Length of Theta[i].
                            Theta,             // Component parameters.
                            length_K,          // Length of K.
                            K,                 // Numbers of bins v or numbers of nearest neighbours k.
                            length_ymin,       // Length of ymin.
                            ymin,              // Minimum observations.
                            length_ymax,       // Length of ymax.
                            ymax,              // Maximum observations.
                            length_h,          // Length of h.
                            h,                 // Sides of the hypersquare.
                            ar,                // Acceleration rate.
                            Restraints,        // Restraints type.
                            n,                 // Number of observations.
                            Y,                 // Dataset.
                            Y_type,            // Dataset type. 
                            EMStrategy,        // Strategy for EM algorithm.
                            EMVariant,         // EM algorithm variant.
                            EMAcceleration,    // Acceleration for the standard EM algorithm.
                            EMTolerance,       // Tolerance for EM algortihm.
                            EMAccelerationMul, // Acceleration rate for Em algorithm.
                            EMMaxIter,         // Maximum number of iterations in EM algorithm.
                            EMK,               // Number of bins for histogram EM algorithm.
                            NULL,              // Component weights.
                            NULL);             // Mixture parameters.
    
    if (*Error) goto E0;

    *Error = rebmvnorm->REBMIX();

    if (*Error) goto E0;

/// Panic Branislav
    *Error = rebmvnorm->Get(n_iter,       // Number of iterations for optimal case.
                            n_iter_sum,   // Number of iterations in whole run.
/// End
                            summary_k,    // Optimal v or optimal k.
                            summary_h,    // Optimal class widths of length d.
                            summary_y0,   // Optimal origins of length d.
                            summary_ymin, // Optimal minimum observations of length d.
                            summary_ymax, // Optimal maximum observations of length d.
                            summary_IC,   // Optimal information criterion.
                            summary_logL, // Log-likelihood.
                            summary_M,    // Degrees of freedom.
                            summary_c,    // Optimal number of components.
                            W,            // Component weights.
                            theta1,       // Component parameters.
                            theta2,       // Component parameters.
                            NULL,         // Component parameters.
                            opt_length,   // Length of opt_c, opt_IC, opt_logL and opt_D.
                            opt_c,        // Numbers of components for optimal v or for optimal k.
                            opt_IC,       // Information criteria for optimal v or for optimal k.
                            opt_logL,     // Log-likelihoods for optimal v or for optimal k.
                            opt_D,        // Totals of positive relative deviations for optimal v or for optimal k.
                            all_length,   // Length of all_K and all_IC.
                            all_K,        // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                            all_IC);      // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.

    if (*Error) goto E0;

E0: if (rebmvnorm) delete rebmvnorm;
} // RREBMVNORM

// Returns classified observations in R.

void RCLSMVNORM(INT    *n,      // Total number of independent observations.
                double *X,      // Pointer to the input array X.
                INT    *s,      // Number of classes.
                INT    *o,      // Number of input REBMVNORM objects.
                INT    *d,      // Number of independent random variables in REBMVNORM objects.
                INT    *c,      // Number of components.
                double *W,      // Component weights.
                char   **pdf,   // Component parameters.
                double *theta1, // Component parameters.
                double *theta2, // Component parameters.
                double *P,      // Prior probabilities.
                INT    *Z,      // Pointer to the output array Z.
                INT    *Error)  // Error code.
{
    Rebmvnorm            *rebmvnorm = NULL;
    INT                  **C = NULL;
    INT                  A[4];
    FLOAT                ***Q = NULL;
    FLOAT                **Y = NULL;
    CompnentDistribution ****Theta = NULL;
    FLOAT                CmpDist, MixDist, MaxMixDist;
    INT                  i, j, k, l, m, dmax = 0;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    C = (INT**)malloc(*s * sizeof(INT*));

    *Error = NULL == C; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *s; j++) {
        C[j] = (INT*)malloc(*o * sizeof(INT));

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

    Theta = new CompnentDistribution*** [(unsigned INT)(*s)];

    *Error = NULL == Theta; if (*Error) goto E0;

    for (j = 0; j < *s; j++) {
        Theta[j] = new CompnentDistribution** [(unsigned INT)(*o)];

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        for (k = 0; k < *o; k++) {
            Theta[j][k] = new CompnentDistribution* [(unsigned INT)C[j][k]];

            *Error = NULL == Theta[j][k]; if (*Error) goto E0;

            for (l = 0; l < C[j][k]; l++) {
                Theta[j][k][l] = new CompnentDistribution(rebmvnorm);

                *Error = NULL == Theta[j][k][l]; if (*Error) goto E0;

                A[0] = d[k];
                A[1] = A[2] = d[k] * d[k];
                A[3] = 1;

                *Error = Theta[j][k][l]->Realloc(d[k], 4, A);

                if (*Error) goto E0;
            }
        }
    }

    i = 0;

    for (j = 0; j < *s; j++) {
        for (k = 0; k < *o; k++) {
            for (l = 0; l < C[j][k]; l++) {
                for (m = 0; m < d[k]; m++) {
                    if (!strcmp(pdf[i], "normal")) {
                        Theta[j][k][l]->pdf_[m] = pfNormal;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                    }
                    else {
                        *Error = 1; goto E0;
                    }

                    i++;
                }
            }
        }
    }

    i = 0;

    for (j = 0; j < *s; j++) {
        for (k = 0; k < *o; k++) {
            for (l = 0; l < C[j][k]; l++) {
                for (m = 0; m < d[k] * d[k]; m++) {
                    Theta[j][k][l]->Theta_[1][m] = theta2[i];

                    i++;
                }
            }
        }
    }

    for (j = 0; j < *s; j++) {
        for (k = 0; k < *o; k++) {
            for (l = 0; l < C[j][k]; l++) {
                *Error = Cholinvdet(d[k], Theta[j][k][l]->Theta_[1], Theta[j][k][l]->Theta_[2], Theta[j][k][l]->Theta_[3]);

                if (*Error) goto E0;
            }
        }
    }

    dmax = d[0]; for (i = 1; i < *o; i++) if (d[i] > dmax) dmax = d[i];

    Y = (FLOAT**)malloc(dmax * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < dmax; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    for (i = 0; i < *n; i++) {
        Z[i] = 1; MaxMixDist = (FLOAT)0.0;

        for (j = 0; j < *s; j++) {
            k = 0; MixDist = (FLOAT)1.0;

            for (l = 0; l < *o; l++) {
                for (m = 0; m < d[l]; m++) {
                    Y[m][0] = X[i + (*n) * (m + k)];
                }

                *Error = rebmvnorm->MixtureDist(0, Y, C[j][l], Q[j][l], Theta[j][l], &CmpDist);

                if (*Error) goto E0;

                k += d[l]; MixDist *= CmpDist;
            }

            MixDist *= P[j];

            if (MixDist > MaxMixDist) {
                Z[i] = j + 1; MaxMixDist = MixDist;
            }
        }
    }

E0: if (Y) {
        for (i = 0; i < dmax; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (Theta) {
        for (i = 0; i < *s; i++) {
            if (Theta[i]) {
                for (j = 0; j < *o; j++) {
                    if (Theta[i][j]) {
                        for (k = 0; k < C[i][j]; k++) {
                            if (Theta[i][j][k]) delete Theta[i][j][k];
                        }

                        delete[] Theta[i][j];
                    }
                }

                delete[] Theta[i];
            }
        }

        delete[] Theta;
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

// Returns clustered observations in R.

void RCLRMVNORM(INT    *n,      // Total number of independent observations.
                double *X,      // Pointer to the input array X.
                INT    *d,      // Number of independent random variables.
                INT    *c,      // Number of components.
                double *W,      // Component weights.
                char   **pdf,   // Component parameters.
                double *theta1, // Component parameters.
                double *theta2, // Component parameters.
                INT    *Z,      // Pointer to the output array Z.
                INT    *Error)  // Error code.
{
    Rebmvnorm            *rebmvnorm = NULL;
    FLOAT                **Y = NULL;
    INT                  A[4];
    CompnentDistribution **Theta = NULL;
    FLOAT                CmpDist, MaxCmpDist;
    INT                  i, j, k;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    rebmvnorm->length_pdf_ = *d;

    Theta = new CompnentDistribution* [(unsigned INT)(*c)];

    *Error = NULL == Theta; if (*Error) goto E0;

    A[0] = *d;
    A[1] = A[2] = (*d) * (*d);
    A[3] = 1;

    for (j = 0; j < *c; j++) {
        Theta[j] = new CompnentDistribution(rebmvnorm);

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        *Error = Theta[j]->Realloc(*d, 4, A);

        if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < *c; j++) {
        for (k = 0; k < *d; k++) {
            if (!strcmp(pdf[i], "normal")) {
                Theta[j]->pdf_[k] = pfNormal;

                Theta[j]->Theta_[0][k] = theta1[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    i = 0;

    for (j = 0; j < *c; j++) {
        for (k = 0; k < (*d) * (*d); k++) {
            Theta[j]->Theta_[1][k] = theta2[i];

            i++;
        }
    }

    for (j = 0; j < *c; j++) {
        *Error = Cholinvdet(*d, Theta[j]->Theta_[1], Theta[j]->Theta_[2], Theta[j]->Theta_[3]);

        if (*Error) goto E0;
    }

    Y = (FLOAT**)malloc(*d * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < *d; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *d; j++) {
            Y[j][0] = X[i + (*n) * j];
        }

        Z[i] = 1; MaxCmpDist = (FLOAT)0.0;

        for (j = 0; j < *c; j++) {
            *Error = rebmvnorm->ComponentDist(0, Y, Theta[j], &CmpDist, NULL);

            if (*Error) goto E0;

            CmpDist *= W[j];

            if (CmpDist > MaxCmpDist) {
                Z[i] = j + 1; MaxCmpDist = CmpDist;
            }
        }
    }

E0: if (Y) {
        for (i = 0; i < *d; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (Theta) {
        for (i = 0; i < *c; i++) {
            if (Theta[i]) delete Theta[i];
        }

        delete[] Theta;
    }

    if (rebmvnorm) delete rebmvnorm;
} // RCLRMVNORM

void RPreprocessingKNNMVNORM(INT    *k,     // k-nearest neighbours.
                             double *h,     // Normalizing vector.
                             INT    *n,     // Total number of independent observations.
                             INT    *d,     // Number of independent random variables.
                             double *x,     // Pointer to the input array x.
                             double *y,     // Pointer to the output array y.
                             INT    *Error) // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    FLOAT     **Y = NULL;
    INT       i, j, l;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;
    rebmvnorm->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 3) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmvnorm->PreprocessingKNN(*k, h, Y);

    if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 3; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;
} // RPreprocessingKNNMVNORM

void RPreprocessingKDEMVNORM(double *h,     // Sides of the hypersquare.
                             INT    *n,     // Total number of independent observations.
                             INT    *d,     // Number of independent random variables.
                             double *x,     // Pointer to the input array x.
                             double *y,     // Pointer to the output array y.
                             INT    *Error) // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    FLOAT     **Y = NULL;
    INT       i, j, l;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;
    rebmvnorm->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 2) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmvnorm->PreprocessingKDE(h, Y);

    if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 2; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;
} // RPreprocessingKDEMVNORM

void RPreprocessingHMVNORM(double *h,          // Sides of the hypersquare.
                           double *y0,         // Origins.
                           double *ymin,       // Minimum observations.
                           double *ymax,       // Maximum observations.
                           INT    *k,          // Total number of bins.
                           INT    *n,          // Total number of independent observations.
                           INT    *d,          // Number of independent random variables.
                           double *x,          // Pointer to the input array x.
                           double *y,          // Pointer to the output array y.
                           INT    *Error)      // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    FLOAT     **Y = NULL;
    INT       i, j, l;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;
    rebmvnorm->length_pdf_ = *d;

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmvnorm->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 1) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmvnorm->PreprocessingH(h, y0, ymin, ymax, k, Y);

    if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 1; j++) {
        for (l = 0; l < *k; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;
} // RPreprocessingHMVNORM

void RInformationCriterionKNNMVNORM(double *h,            // Sides of the hypersquare.
                                    INT    *k,            // k-nearest neighbours.
                                    char   **Criterion,   // Information criterion type.
                                    INT    *c,            // Number of components.
                                    double *W,            // Component weights.
                                    INT    *length_pdf,   // Length of pdf.
                                    INT    *length_Theta, // Length of Theta.
                                    INT    *length_theta, // Length of Theta[i].
                                    char   **pdf,         // Parametric family types.
                                    double *Theta,        // Component parameters.
                                    INT    *n,            // Number of observations.
                                    double *x,            // Dataset.
                                    double *IC,           // Information criterion.
                                    double *logL,         // log-likelihood.
                                    INT    *M,            // Degrees of freedom.
                                    double *D,            // Total of positive relative deviations.
                                    INT    *Error)        // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    FLOAT     **Y = NULL;
    INT       i, j, l, m;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

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

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmvnorm->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    *Error = NULL == rebmvnorm->IniTheta_; if (*Error) goto E0;

    rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

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

    rebmvnorm->MixTheta_ = new CompnentDistribution*[(unsigned INT)(*c)];

    *Error = NULL == rebmvnorm->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        *Error = NULL == rebmvnorm->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmvnorm->length_pdf_; j++) {
            rebmvnorm->MixTheta_[i]->pdf_[j] = rebmvnorm->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_Theta_; j++) if (rebmvnorm->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmvnorm->length_theta_[j]; m++) {
                rebmvnorm->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmvnorm->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 3) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;

        if (i < rebmvnorm->length_pdf_) {
            for (j = 0; j < rebmvnorm->nr_; j++) Y[i][j] = rebmvnorm->Y_[i][j];
        }
    }

    *Error = rebmvnorm->PreprocessingKNN(*k, h, Y);

    if (*Error) goto E0;

    rebmvnorm->cmax_ = *c;

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        if (*Error) goto E0;
    }

    *Error = rebmvnorm->InformationCriterionKNN(*k,
                                                Y,
                                                *c,
                                                rebmvnorm->W_,
                                                rebmvnorm->MixTheta_,
                                                IC,
                                                logL,
                                                M,
                                                D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;
} // RInformationCriterionKNNMVNORM

void RInformationCriterionKDEMVNORM(double *h,            // Sides of the hypersquare.
                                    char   **Criterion,   // Information criterion type.
                                    INT    *c,            // Number of components.
                                    double *W,            // Component weights.
                                    INT    *length_pdf,   // Length of pdf.
                                    INT    *length_Theta, // Length of Theta.
                                    INT    *length_theta, // Length of Theta[i].
                                    char   **pdf,         // Parametric family types.
                                    double *Theta,        // Component parameters.
                                    INT    *n,            // Number of observations.
                                    double *x,            // Dataset.
                                    double *IC,           // Information criterion.
                                    double *logL,         // log-likelihood.
                                    INT    *M,            // Degrees of freedom.
                                    double *D,            // Total of positive relative deviations.
                                    INT    *Error)        // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    FLOAT     **Y = NULL;
    FLOAT     logV;
    INT       i, j, l, m;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

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

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmvnorm->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    *Error = NULL == rebmvnorm->IniTheta_; if (*Error) goto E0;

    rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

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

    rebmvnorm->MixTheta_ = new CompnentDistribution*[(unsigned INT)(*c)];

    *Error = NULL == rebmvnorm->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        *Error = NULL == rebmvnorm->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmvnorm->length_pdf_; j++) {
            rebmvnorm->MixTheta_[i]->pdf_[j] = rebmvnorm->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_Theta_; j++) if (rebmvnorm->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmvnorm->length_theta_[j]; m++) {
                rebmvnorm->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmvnorm->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 2) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;

        if (i < rebmvnorm->length_pdf_) {
            for (j = 0; j < rebmvnorm->nr_; j++) Y[i][j] = rebmvnorm->Y_[i][j];
        }
    }

    *Error = rebmvnorm->PreprocessingKDE(h, Y);

    if (*Error) goto E0;

    rebmvnorm->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        if (*Error) goto E0;
    }

    *Error = rebmvnorm->InformationCriterionKDE(logV,
                                                Y,
                                                *c,
                                                rebmvnorm->W_,
                                                rebmvnorm->MixTheta_,
                                                IC,
                                                logL,
                                                M,
                                                D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;
} // RInformationCriterionKDEMVNORM

void RInformationCriterionHMVNORM(double *h,            // Sides of the hypersquare.
                                  double *y0,           // Origins.
                                  double *ymin,         // Minimum observations.
                                  double *ymax,         // Maximum observations.
                                  INT    *k,            // Total number of bins.
                                  char   **Criterion,   // Information criterion type.
                                  INT    *c,            // Number of components.
                                  double *W,            // Component weights.
                                  INT    *length_pdf,   // Length of pdf.
                                  INT    *length_Theta, // Length of Theta.
                                  INT    *length_theta, // Length of Theta[i].
                                  char   **pdf,         // Parametric family types.
                                  double *Theta,        // Component parameters.
                                  INT    *n,            // Number of observations.
                                  double *x,            // Dataset.
                                  double *IC,           // Information criterion.
                                  double *logL,         // log-likelihood.
                                  INT    *M,            // Degrees of freedom.
                                  double *D,            // Total of positive relative deviations.
                                  INT    *Error)        // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    FLOAT     **Y = NULL;
    FLOAT     logV;
    INT       i, j, l, m;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

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

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmvnorm->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    *Error = NULL == rebmvnorm->IniTheta_; if (*Error) goto E0;

    rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

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

    rebmvnorm->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    *Error = NULL == rebmvnorm->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        *Error = NULL == rebmvnorm->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmvnorm->length_pdf_; j++) {
            rebmvnorm->MixTheta_[i]->pdf_[j] = rebmvnorm->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_Theta_; j++) if (rebmvnorm->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmvnorm->length_theta_[j]; m++) {
                rebmvnorm->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->length_pdf_ * sizeof(FLOAT*));

    *Error = NULL == rebmvnorm->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == rebmvnorm->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 1) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmvnorm->PreprocessingH(h, y0, ymin, ymax, k, Y);

    if (*Error) goto E0;

    rebmvnorm->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        if (*Error) goto E0;
    }

    *Error = rebmvnorm->InformationCriterionH(logV,
                                              *k,
                                              Y,
                                              *c,
                                              rebmvnorm->W_,
                                              rebmvnorm->MixTheta_,
                                              IC,
                                              logL,
                                              M,
                                              D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;
} // RInformationCriterionHMVNORM

void RInformationCriterionKMVNORM(double *h,            // Sides of the hypersquare.
                                  char   **Criterion,   // Information criterion type.
                                  INT    *c,            // Number of components.
                                  double *W,            // Component weights.
                                  INT    *length_pdf,   // Length of pdf.
                                  INT    *length_Theta, // Length of Theta.
                                  INT    *length_theta, // Length of Theta[i].
                                  char   **pdf,         // Parametric family types.
                                  double *Theta,        // Component parameters.
                                  INT    *n,            // Number of observations.
                                  double *x,            // Dataset.
                                  double *IC,           // Information criterion.
                                  double *logL,         // log-likelihood.
                                  INT    *M,            // Degrees of freedom.
                                  double *D,            // Total of positive relative deviations.
                                  INT    *Error)        // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    FLOAT     **Y = NULL;
    FLOAT     logV;
    INT       i, j, l, m;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

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

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmvnorm->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    *Error = NULL == rebmvnorm->IniTheta_; if (*Error) goto E0;

    rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

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

    rebmvnorm->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    *Error = NULL == rebmvnorm->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        *Error = NULL == rebmvnorm->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmvnorm->length_pdf_; j++) {
            rebmvnorm->MixTheta_[i]->pdf_[j] = rebmvnorm->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_Theta_; j++) if (rebmvnorm->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmvnorm->length_theta_[j]; m++) {
                rebmvnorm->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 1) * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 1; j++) {
        Y[j] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        *Error = NULL == Y[j]; if (*Error) goto E0;

        for (l = 0; l < rebmvnorm->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    rebmvnorm->n_ = 0;

    for (l = 0; l < rebmvnorm->nr_; l++) {
        rebmvnorm->n_ += (INT)Y[rebmvnorm->length_pdf_][l];
    }

    rebmvnorm->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        if (*Error) goto E0;
    }

    rebmvnorm->Y_type_ = 1;

    *Error = rebmvnorm->InformationCriterionH(logV,
                                              *n,
                                              Y,
                                              *c,
                                              rebmvnorm->W_,
                                              rebmvnorm->MixTheta_,
                                              IC,
                                              logL,
                                              M,
                                              D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;
} // RInformationCriterionKMVNORM

void RCombineComponentsMVNORM(INT    *c,            // Number of components.
                              double *W,            // Component weights.
                              INT    *length_pdf,   // Length of pdf.
                              INT    *length_Theta, // Length of Theta.
                              INT    *length_theta, // Length of Theta[i].
                              char   **pdf,         // Parametric family types.
                              double *MixTheta,     // Mixture parameters.
                              INT    *n,            // Number of observations.
                              double *Y,            // Dataset.
                              INT    *Y_type,       // Dataset type.
/// Panic Branislav
                              char   **Rule,        // Mergining rule
/// End
                              double *tau,          // Conditional probabilities.
                              INT    *F,            // From components.
                              INT    *T,            // To components.
                              double *EN,           // Entropy.
                              double *ED,           // Entropy decrease.
                              INT    *Error)        // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    INT       i;

    rebmvnorm = new Rebmvnorm;

    *Error = NULL == rebmvnorm; if (*Error) goto E0;

    rebmvnorm->Set(NULL,         // Preprocessing type.
                   c,            // Maximum number of components.
                   NULL,         // Minimum number of components.
                   NULL,         // Information criterion type.
                   length_pdf,   // Number of independent random variables.
                   NULL,         // Types of variables.
                   length_pdf,   // Length of pdf.
                   pdf,          // Parametric family types.
                   length_Theta, // Length of Theta.
                   length_theta, // Length of Theta[i].
                   NULL,         // Component parameters.
                   NULL,         // Length of K.
                   NULL,         // Numbers of bins v or numbers of nearest neighbours k.
                   NULL,         // Length of ymin.
                   NULL,         // Minimum observations.
                   NULL,         // Length of ymax.
                   NULL,         // Maximum observations.
                   NULL,         // Length of h.
                   NULL,         // Sides of the hypersquare.
                   NULL,         // Acceleration rate.
                   NULL,         // Restraints type.
                   n,            // Number of observations.
                   Y,            // Dataset.
                   Y_type,       // Dataset type. 
                   NULL,         // Strategy for EM algorithm.
                   NULL,         // EM algorithm variant.
                   NULL,         // Acceleration for the standard EM algorithm.
                   NULL,         // Tolerance for EM algortihm.
                   NULL,         // Acceleration rate for Em algorithm.
                   NULL,         // Maximum number of iterations in EM algorithm.
                   NULL,         // Number of bins for histogram EM algorithm.
                   W,            // Component weights.
                   MixTheta);    // Mixture parameters.

    for (i = 0; i < rebmvnorm->cmax_; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        if (*Error) goto E0;
    }

/// Panic Branislav
    if (!strcmp(Rule[0], "Entropy")) {
        *Error = rebmvnorm->CombineComponentsEntropy(*c,
            rebmvnorm->W_,
            rebmvnorm->MixTheta_,
            tau,
            F,
            T,
            EN,
            ED);

        if (*Error) goto E0;
    }
    else
    if (!strcmp(Rule[0], "Demp")) {
        *Error = rebmvnorm->CombineComponentsDemp(*c,
            rebmvnorm->W_,
            rebmvnorm->MixTheta_,
            tau,
            F,
            T,
            EN,
            ED);

        if (*Error) goto E0;
    }
    else {
        *Error = 1; goto E0;
    }
/// End

E0: if (rebmvnorm) delete rebmvnorm;
} // RCombineComponentsMVNORM

void RTvtNormalPdf(INT *n, double *x, double *y, double *Mean, double *Sigma, double *f)
{
    INT   j;
    FLOAT Adet, Ainv[4], xi, yi, z;

    Adet = Sigma[0] * Sigma[3] - Sigma[1] * Sigma[2];

    if (Adet > FLOAT_MIN) {
        Ainv[0] = Sigma[3] / Adet; Ainv[1] = -Sigma[2] / Adet; Ainv[2] = -Sigma[1] / Adet; Ainv[3] = Sigma[0] / Adet;

        for (j = 0; j < *n; j++) {
            xi = x[j] - Mean[0]; yi = y[j] - Mean[1];

            z = -(FLOAT)0.5 * (Ainv[0] * xi * xi + Ainv[3] * yi * yi) - Ainv[1] * xi * yi - (FLOAT)0.5 * (FLOAT)log(Adet) - LogPi2;

            f[j] = (FLOAT)exp(z);
        }
    }
    else {
        for (j = 0; j < *n; j++) {
            f[j] = (FLOAT)0.0;
        }
    }
} // RTvtNormalPdf

void RMvtNormalPdf(INT *n, double *X, INT *d, double *Mean, double *Sigma, double *f)
{
    INT   i, j, k, Error;
    FLOAT logAdet, *Ainv, y, yi, yk;

    Ainv = (FLOAT*)malloc(*d * *d * sizeof(FLOAT));

    Error = NULL == Ainv; if (Error) goto E0;

    Error = Cholinvdet(*d, Sigma, Ainv, &logAdet);

    if (Error) {
        for (j = 0; j < *n; j++) {
            f[j] = (FLOAT)0.0;
        }
    }
    else {
        for (j = 0; j < *n; j++) {
            y = (FLOAT)0.0;

            for (i = 0; i < *d; i++) {
                yi = X[*n * i + j] - Mean[i]; y += (FLOAT)0.5 * Ainv[*d * i + i] * yi * yi;

                for (k = i + 1; k < *d; k++) {
                    yk = X[*n * k + j] - Mean[k]; y += Ainv[*d * k + i] * yi * yk;
                }
            }

            f[j] = (FLOAT)exp(-y - *d * LogSqrtPi2 - (FLOAT)0.5 * logAdet);
        }
    }

E0: if (Ainv) free(Ainv);
} // RMvtNormalPdf

/// Panic Branislav
void REMMVNORM(INT    *d,                 // Number of independent random variables.
               INT    *n,                 // Number of observations.
               double *Y,                 // Dataset.
               INT    *Y_type,            // Dataset type.
               char   **pdf,              // Parametric family types.
               INT    *c,                 // Input number of components.
               double *W,                 // Input weights of components.
               double *Theta1,            // Input parameters theta 1.
               double *Theta2,            // Input parameters theta 2.
               char   **EMVariant,        // EM algorithm variant.
               char   **EMAcceleration,   // Acceleration for the standard EM algorithm.
               double *EMTolerance,       // Tolerance for EM algortihm.
               double *EMAccelerationMul, // Acceleration rate for Em algorithm.
               INT    *EMMaxIter,         // Maximum number of iterations in EM algorithm.
               INT    *EMK,               // Number of bins for histogram EM algorithm.
               INT    *EMMerge,           // Remove zero components.   
               INT    *n_iter,            // Number of iterations for optimal case.
               double *summary_logL,      // Log-likelihood.
               INT    *summary_M,         // Degrees of freedom.
               INT    *Error)             // Error number.)
{

    INT error = 0, i = 0, j = 0, l = 0, length_Theta = 4;

    INT A[4];

    A[0] = *d;
    A[1] = A[2] = (*d) * (*d);
    A[3] = 1;

    Rebmvnorm *rebmvnorm = NULL;

    rebmvnorm = new Rebmvnorm;

    error = NULL == rebmvnorm; if (error) goto E0;

    error = rebmvnorm->Set(NULL,              // Preprocessing type.
                           c,                 // Maximum number of components.
                           c,                 // Minimum number of components.
                           NULL,              // Information criterion type.
                           d,                 // Number of independent random variables.
                           NULL,              // Types of variables.
                           d,                 // Length of pdf.
                           pdf,               // Parametric family types.
                           &length_Theta,     // Length of Theta.
                           A,                 // Length of Theta[i].
                           NULL,              // Component parameters.
                           NULL,              // Length of K.
                           NULL,              // Numbers of bins v or numbers of nearest neighbours k.
                           NULL,              // Length of ymin.
                           NULL,              // Minimum observations.
                           NULL,              // Length of ymax.
                           NULL,              // Maximum observations.
                           NULL,              // Length of h.
                           NULL,              // Sides of the hypersquare.
                           NULL,              // Acceleration rate.
                           NULL,              // Restraints type.
                           n,                 // Number of observations.
                           Y,                 // Dataset.
                           Y_type,            // Dataset type.
                           NULL,              // Strategy for EM algorithm.
                           EMVariant,         // EM algorithm variant.
                           EMAcceleration,    // Acceleration for the standard EM algorithm.
                           EMTolerance,       // Tolerance for EM algortihm.
                           EMAccelerationMul, // Acceleration rate for Em algorithm.
                           EMMaxIter,         // Maximum number of iterations in EM algorithm.
                           EMK,               // Number of bins for histogram EM algorithm.
                           NULL,              // Component weights.
                           NULL);             // Mixture parameters.

    rebmvnorm->EM_strategy_ = strategy_single;

    rebmvnorm->MixTheta_ = new CompnentDistribution*[(unsigned INT)(rebmvnorm->cmax_)];

    error = NULL == rebmvnorm->MixTheta_; if (error) goto E0;

    for (l = 0; l < rebmvnorm->cmax_; l++) {
        rebmvnorm->MixTheta_[l] = new CompnentDistribution(rebmvnorm);

        error = NULL == rebmvnorm->MixTheta_[l]; if (error) goto E0;

        error = rebmvnorm->MixTheta_[l]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        if (error) goto E0;

        for (i = 0; i < rebmvnorm->length_pdf_; i++) {
            rebmvnorm->MixTheta_[l]->pdf_[i] = rebmvnorm->IniTheta_->pdf_[i];
        }
    }

    i = 0;

    for (l = 0; l < *c; l++) {
        for (j = 0; j < rebmvnorm->length_theta_[0]; j++) {
            rebmvnorm->MixTheta_[l]->Theta_[0][j] = Theta1[i]; 
            
            i++;
        }
    }

    i = 0;

    for (l = 0; l < *c; l++) {
        for (j = 0; j < rebmvnorm->length_theta_[1]; j++) {
            rebmvnorm->MixTheta_[l]->Theta_[1][j] = Theta2[i];
            
            i++;
        }
    }

    rebmvnorm->W_ = (double*)malloc(rebmvnorm->cmax_ * sizeof(double));

    error = rebmvnorm->W_ == NULL; if (error) goto E0;

    for (l = 0; l < *c; l++) {
        rebmvnorm->W_[l] = W[l];

        error = Cholinvdet(*d, rebmvnorm->MixTheta_[l]->Theta_[1], rebmvnorm->MixTheta_[l]->Theta_[2], rebmvnorm->MixTheta_[l]->Theta_[3]);

        if (error) goto E0;

    }

    error = rebmvnorm->EMInitialize();

    if (error) goto E0;

    if (*EMMerge) {
        rebmvnorm->EM_->merge_ = merge_naive;
    }

    error = rebmvnorm->EMRun(c, rebmvnorm->W_, rebmvnorm->MixTheta_);

    if (error) goto E0;

    error = rebmvnorm->EM_->LogLikelihood(*c, rebmvnorm->W_, rebmvnorm->MixTheta_, summary_logL);

    if (error) goto E0;

    error = rebmvnorm->DegreesOffreedom(*c, rebmvnorm->MixTheta_, summary_M);

    if (error) goto E0;

    rebmvnorm->summary_.c = *c;

    error = rebmvnorm->Get(n_iter, // Number of iterations for optimal case.
                           NULL,   // Number of iterations in whole run.
                           NULL,   // Optimal v or optimal k.
                           NULL,   // Optimal class widths of length d.
                           NULL,   // Optimal origins of length d.
                           NULL,   // Optimal minimum observations of length d.
                           NULL,   // Optimal maximum observations of length d.
                           NULL,   // Optimal information criterion.
                           NULL,   // Log-likelihood.
                           NULL,   // Degrees of freedom.
                           c,      // Optimal number of components.
                           W,      // Component weights.
                           Theta1, // Component parameters.
                           Theta2, // Component parameters.
                           NULL,   // Component parameters.
                           NULL,   // Length of opt_c, opt_IC, opt_logL and opt_D.
                           NULL,   // Numbers of components for optimal v or for optimal k.
                           NULL,   // Information criteria for optimal v or for optimal k.
                           NULL,   // Log-likelihoods for optimal v or for optimal k.
                           NULL,   // Totals of positive relative deviations for optimal v or for optimal k.
                           NULL,   // Length of all_K and all_IC.
                           NULL,   // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                           NULL);  // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.

E0:
    if (rebmvnorm) delete rebmvnorm;

    *Error = error;
}
/// End
}
