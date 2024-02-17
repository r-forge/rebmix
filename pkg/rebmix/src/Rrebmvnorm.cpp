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

    *Error = EOK;

    rngmvnorm = new Rngmvnorm;

    R_CHECK(NULL == rngmvnorm, ERRNGMVNORM);

    rngmvnorm->IDum_ = *IDum;
    rngmvnorm->length_pdf_ = *d;
    rngmvnorm->c_ = *c;

    rngmvnorm->N_ = (INT*)malloc(rngmvnorm->c_ * sizeof(INT));

    R_CHECK(NULL == rngmvnorm->N_, ERRNGMVNORM);

    for (i = 0; i < rngmvnorm->c_; i++) rngmvnorm->N_[i] = N[i];

    rngmvnorm->IniTheta_ = new CompnentDistribution(rngmvnorm);

    R_CHECK(NULL == rngmvnorm->IniTheta_, ERRNGMVNORM);

    rngmvnorm->length_pdf_ = *length_pdf;

    rngmvnorm->length_Theta_ = *length_Theta;

    rngmvnorm->length_theta_ = (INT*)malloc(rngmvnorm->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rngmvnorm->length_theta_, ERRNGMVNORM);

    *Error = rngmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != EOK, *Error);

    for (i = 0; i < rngmvnorm->length_pdf_; i++) {
        rngmvnorm->IniTheta_->pdf_[i] = pfNormal;
    }

    rngmvnorm->MixTheta_ = new CompnentDistribution*[(unsigned INT)rngmvnorm->c_];

    R_CHECK(NULL == rngmvnorm->MixTheta_, ERRNGMVNORM);

    for (i = 0; i < rngmvnorm->c_; i++) {
        rngmvnorm->MixTheta_[i] = new CompnentDistribution(rngmvnorm);

        R_CHECK(NULL == rngmvnorm->MixTheta_[i], ERRNGMVNORM);

        *Error = rngmvnorm->MixTheta_[i]->Realloc(rngmvnorm->length_pdf_, rngmvnorm->length_Theta_, rngmvnorm->length_theta_);

        R_CHECK(*Error != EOK, *Error);
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

    R_CHECK(*Error != EOK, *Error);

    *n = rngmvnorm->n_; i = 0;

    for (j = 0; j < rngmvnorm->length_pdf_; j++) {
        for (k = 0; k < rngmvnorm->n_; k++) {
            Y[i] = rngmvnorm->Y_[j][k]; i++;
        }
    }

    for (i = 0; i < rngmvnorm->n_; i++) {
        Z[i] = rngmvnorm->Z_[i];
    }

EEXIT:

    if (rngmvnorm) delete rngmvnorm;

    R_RETURN(*Error);
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
                INT    *opt_length,     // Length of opt_c, opt_IC, opt_logL, opt_Dmin and opt_D.
                INT    *opt_c,          // Numbers of components for optimal v or for optimal k.
                double *opt_IC,         // Information criteria for optimal v or for optimal k.
                double *opt_logL,       // Log-likelihoods for optimal v or for optimal k.
                double *opt_Dmin,       // Dmin for optimal v or for optimal k.
                double *opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
                INT    *all_length,     // Length of all_K and all_IC.
                INT    *all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                double *all_IC,         // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
                INT    *Error)          // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERREBMVNORM);

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
    
    R_CHECK(*Error != EOK, *Error);

    *Error = rebmvnorm->REBMIX();

    R_CHECK(*Error != EOK, *Error);

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
                            opt_length,   // Length of opt_c, opt_IC, opt_logL, opt_Dmin and opt_D.
                            opt_c,        // Numbers of components for optimal v or for optimal k.
                            opt_IC,       // Information criteria for optimal v or for optimal k.
                            opt_logL,     // Log-likelihoods for optimal v or for optimal k.
                            opt_Dmin,     // Dmin for optimal v or for optimal k.
                            opt_D,        // Totals of positive relative deviations for optimal v or for optimal k.
                            all_length,   // Length of all_K and all_IC.
                            all_K,        // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                            all_IC);      // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.

    R_CHECK(*Error != EOK, *Error);

EEXIT:

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
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
    CompnentDistribution ****Theta = NULL;
    FLOAT                CmpDist, MixDist, MaxMixDist, ***Q = NULL , **Y = NULL;
    INT                  A[4], **C = NULL, i, j, k, l, m, dmax = 0;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERCLSMVNORM);

    C = (INT**)malloc(*s * sizeof(INT*));

    R_CHECK(NULL == C, ERCLSMVNORM);

    i = 0;

    for (j = 0; j < *s; j++) {
        C[j] = (INT*)malloc(*o * sizeof(INT));

        R_CHECK(NULL == C[j], ERCLSMVNORM);

        for (k = 0; k < *o; k++) {
            C[j][k] = c[i]; i++;
        }
    }

    Q = (FLOAT***)malloc(*s * sizeof(FLOAT**));

    R_CHECK(NULL == Q, ERCLSMVNORM);

    i = 0;

    for (j = 0; j < *s; j++) {
        Q[j] = (FLOAT**)malloc(*o * sizeof(FLOAT*));

        R_CHECK(NULL == Q[j], ERCLSMVNORM);

        for (k = 0; k < *o; k++) {
            Q[j][k] = (FLOAT*)malloc(C[j][k] * sizeof(FLOAT));

            R_CHECK(NULL == Q[j][k], ERCLSMVNORM);

            for (l = 0; l < C[j][k]; l++) {
                Q[j][k][l] = W[i]; i++;
            }
        }
    }

    Theta = new CompnentDistribution*** [(unsigned INT)(*s)];

    R_CHECK(NULL == Theta, ERCLSMVNORM);

    for (j = 0; j < *s; j++) {
        Theta[j] = new CompnentDistribution** [(unsigned INT)(*o)];

        R_CHECK(NULL == Theta[j], ERCLSMVNORM);

        for (k = 0; k < *o; k++) {
            Theta[j][k] = new CompnentDistribution* [(unsigned INT)C[j][k]];

            R_CHECK(NULL == Theta[j][k], ERCLSMVNORM);

            for (l = 0; l < C[j][k]; l++) {
                Theta[j][k][l] = new CompnentDistribution(rebmvnorm);

                R_CHECK(NULL == Theta[j][k][l], ERCLSMVNORM);

                A[0] = d[k];
                A[1] = A[2] = d[k] * d[k];
                A[3] = 1;

                *Error = Theta[j][k][l]->Realloc(d[k], 4, A);

                R_CHECK(*Error != EOK, *Error);
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
                        R_CHECK(1, ERCLSMVNORM);
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

                R_CHECK(*Error != EOK, *Error);
            }
        }
    }

    dmax = d[0]; for (i = 1; i < *o; i++) if (d[i] > dmax) dmax = d[i];

    Y = (FLOAT**)malloc(dmax * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERCLSMVNORM);

    for (i = 0; i < dmax; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERCLSMVNORM);
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

                R_CHECK(*Error != EOK, *Error);

                k += d[l]; MixDist *= CmpDist;
            }

            MixDist *= P[j];

            if (MixDist > MaxMixDist) {
                Z[i] = j + 1; MaxMixDist = MixDist;
            }
        }
    }

EEXIT:

    if (Y) {
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

    R_RETURN(*Error);
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
    CompnentDistribution **Theta = NULL;
    FLOAT                CmpDist, MaxCmpDist, **Y = NULL;
    INT                  A[4], i, j, k;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERCLRMVNORM);

    rebmvnorm->length_pdf_ = *d;

    Theta = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == Theta, ERCLRMVNORM);

    A[0] = *d;
    A[1] = A[2] = (*d) * (*d);
    A[3] = 1;

    for (j = 0; j < *c; j++) {
        Theta[j] = new CompnentDistribution(rebmvnorm);

        R_CHECK(NULL == Theta[j], ERCLRMVNORM);

        *Error = Theta[j]->Realloc(*d, 4, A);

        R_CHECK(*Error != EOK, *Error);
    }

    i = 0;

    for (j = 0; j < *c; j++) {
        for (k = 0; k < *d; k++) {
            if (!strcmp(pdf[i], "normal")) {
                Theta[j]->pdf_[k] = pfNormal;

                Theta[j]->Theta_[0][k] = theta1[i];
            }
            else {
                R_CHECK(1, ERCLRMVNORM);
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

        R_CHECK(*Error != EOK, *Error);
    }

    Y = (FLOAT**)malloc(*d * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERCLRMVNORM);

    for (i = 0; i < *d; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERCLRMVNORM);
    }

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *d; j++) {
            Y[j][0] = X[i + (*n) * j];
        }

        Z[i] = 1; MaxCmpDist = (FLOAT)0.0;

        for (j = 0; j < *c; j++) {
            *Error = rebmvnorm->ComponentDist(0, Y, Theta[j], &CmpDist, NULL);

            R_CHECK(*Error != EOK, *Error);

            CmpDist *= W[j];

            if (CmpDist > MaxCmpDist) {
                Z[i] = j + 1; MaxCmpDist = CmpDist;
            }
        }
    }

EEXIT: 
    
    if (Y) {
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

    R_RETURN(*Error);
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

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERPreprocessingKNNMVNORM);

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;
    rebmvnorm->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 3) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERPreprocessingKNNMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERPreprocessingKNNMVNORM);
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmvnorm->PreprocessingKNN(*k, h, Y);

    R_CHECK(*Error != EOK, *Error);

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 3; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

EEXIT:

    if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
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

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERPreprocessingKDEMVNORM);

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;
    rebmvnorm->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 2) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERPreprocessingKDEMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERPreprocessingKDEMVNORM);
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmvnorm->PreprocessingKDE(h, Y);

    R_CHECK(*Error != EOK, *Error);

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 2; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

EEXIT:

    if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
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

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERPreprocessingHMVNORM);

    rebmvnorm->n_ = rebmvnorm->nr_ = *n;
    rebmvnorm->nc_ = rebmvnorm->length_pdf_ = *d;

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmvnorm->Y_, ERPreprocessingHMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmvnorm->Y_[i], ERPreprocessingHMVNORM);
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 1) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERPreprocessingHMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERPreprocessingHMVNORM);
    }

    *Error = rebmvnorm->PreprocessingH(h, y0, ymin, ymax, k, Y);

    R_CHECK(*Error != EOK, *Error);

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 1; j++) {
        for (l = 0; l < *k; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

EEXIT:

    if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
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

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERInformationCriterionKNNMVNORM);

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
        R_CHECK(1, ERInformationCriterionKNNMVNORM);
    }

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmvnorm->W_, ERInformationCriterionKNNMVNORM);

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    R_CHECK(NULL == rebmvnorm->IniTheta_, ERInformationCriterionKNNMVNORM);

    rebmvnorm->nc_ = rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmvnorm->length_theta_, ERInformationCriterionKNNMVNORM);

    *Error = rebmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != EOK, *Error);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmvnorm->IniTheta_->pdf_[i] = pfNormal;
        }
        else {
            R_CHECK(1, ERInformationCriterionKNNMVNORM);
        }
    }

    rebmvnorm->MixTheta_ = new CompnentDistribution*[(unsigned INT)(*c)];

    R_CHECK(NULL == rebmvnorm->MixTheta_, ERInformationCriterionKNNMVNORM);

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        R_CHECK(NULL == rebmvnorm->MixTheta_[i], ERInformationCriterionKNNMVNORM);

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        R_CHECK(*Error != EOK, *Error);
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

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmvnorm->Y_, ERInformationCriterionKNNMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmvnorm->Y_[i], ERInformationCriterionKNNMVNORM);
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 3) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERInformationCriterionKNNMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERInformationCriterionKNNMVNORM);

        if (i < rebmvnorm->length_pdf_) {
            for (j = 0; j < rebmvnorm->nr_; j++) Y[i][j] = rebmvnorm->Y_[i][j];
        }
    }

    *Error = rebmvnorm->PreprocessingKNN(*k, h, Y);

    R_CHECK(*Error != EOK, *Error);

    rebmvnorm->cmax_ = *c;

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        R_CHECK(*Error != EOK, *Error);
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

    R_CHECK(*Error != EOK, *Error);

EEXIT: 
    
    if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
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
    FLOAT     logV, **Y = NULL;
    INT       i, j, l, m;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERInformationCriterionKDEMVNORM);

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
        R_CHECK(1, ERInformationCriterionKDEMVNORM);
    }

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmvnorm->W_, ERInformationCriterionKDEMVNORM);

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    R_CHECK(NULL == rebmvnorm->IniTheta_, ERInformationCriterionKDEMVNORM);

    rebmvnorm->nc_ = rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmvnorm->length_theta_, ERInformationCriterionKDEMVNORM);

    *Error = rebmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != EOK, *Error);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmvnorm->IniTheta_->pdf_[i] = pfNormal;
        }
        else {
            R_CHECK(1, ERInformationCriterionKDEMVNORM);
        }
    }

    rebmvnorm->MixTheta_ = new CompnentDistribution*[(unsigned INT)(*c)];

    R_CHECK(NULL == rebmvnorm->MixTheta_, ERInformationCriterionKDEMVNORM);

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        R_CHECK(NULL == rebmvnorm->MixTheta_[i], ERInformationCriterionKDEMVNORM);

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        R_CHECK(*Error != EOK, *Error);
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

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmvnorm->Y_, ERInformationCriterionKDEMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmvnorm->Y_[i], ERInformationCriterionKDEMVNORM);
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 2) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERInformationCriterionKDEMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERInformationCriterionKDEMVNORM);

        if (i < rebmvnorm->length_pdf_) {
            for (j = 0; j < rebmvnorm->nr_; j++) Y[i][j] = rebmvnorm->Y_[i][j];
        }
    }

    *Error = rebmvnorm->PreprocessingKDE(h, Y);

    R_CHECK(*Error != EOK, *Error);

    rebmvnorm->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        R_CHECK(*Error != EOK, *Error);
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

    R_CHECK(*Error != EOK, *Error);

EEXIT:

    if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
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
    FLOAT     logV, **Y = NULL;
    INT       i, j, l, m;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERInformationCriterionHMVNORM);

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
        R_CHECK(1, ERInformationCriterionHMVNORM);
    }

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmvnorm->W_, ERInformationCriterionHMVNORM);

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    R_CHECK(NULL == rebmvnorm->IniTheta_, ERInformationCriterionHMVNORM);

    rebmvnorm->nc_ = rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmvnorm->length_theta_, ERInformationCriterionHMVNORM);

    *Error = rebmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != EOK, *Error);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmvnorm->IniTheta_->pdf_[i] = pfNormal;
        }
        else {
            R_CHECK(1, ERInformationCriterionHMVNORM);
        }
    }

    rebmvnorm->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == rebmvnorm->MixTheta_, ERInformationCriterionHMVNORM);

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        R_CHECK(NULL == rebmvnorm->MixTheta_[i], ERInformationCriterionHMVNORM);

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        R_CHECK(*Error != EOK, *Error);
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

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmvnorm->Y_, ERInformationCriterionHMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmvnorm->Y_[i], ERInformationCriterionHMVNORM);
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmvnorm->length_pdf_ + 1) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, ERInformationCriterionHMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], ERInformationCriterionHMVNORM);
    }

    *Error = rebmvnorm->PreprocessingH(h, y0, ymin, ymax, k, Y);

    R_CHECK(*Error != EOK, *Error);

    rebmvnorm->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        R_CHECK(*Error != EOK, *Error);
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

    R_CHECK(*Error != EOK, *Error);

EEXIT:

    if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
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
    FLOAT     logV, **Y = NULL;
    INT       i, j, l, m;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERInformationCriterionKMVNORM);

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
        R_CHECK(1, ERInformationCriterionKMVNORM);
    }

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmvnorm->W_, ERInformationCriterionKMVNORM);

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    R_CHECK(NULL == rebmvnorm->IniTheta_, ERInformationCriterionKMVNORM);

    rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmvnorm->length_theta_, ERInformationCriterionKMVNORM);

    *Error = rebmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != EOK, *Error);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmvnorm->IniTheta_->pdf_[i] = pfNormal;
        }
        else {
            R_CHECK(1, ERInformationCriterionKMVNORM);
        }
    }

    rebmvnorm->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == rebmvnorm->MixTheta_, ERInformationCriterionKMVNORM);

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        R_CHECK(NULL == rebmvnorm->MixTheta_[i], ERInformationCriterionKMVNORM);

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        R_CHECK(*Error != EOK, *Error);
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

    R_CHECK(NULL == Y, ERInformationCriterionKMVNORM);

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_ + 1; j++) {
        Y[j] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[j], ERInformationCriterionKMVNORM);

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

        R_CHECK(*Error != EOK, *Error);
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

    R_CHECK(*Error != EOK, *Error);

EEXIT:

    if (Y) {
        for (i = 0; i < rebmvnorm->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
} // RInformationCriterionKMVNORM

void RInformationCriterionMVNORM(char   **Criterion,   // Information criterion type.
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
    INT       i, j, l, m;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERInformationCriterionMVNORM);

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
        R_CHECK(1, ERInformationCriterionMVNORM);
    }

    rebmvnorm->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmvnorm->W_, ERInformationCriterionMVNORM);

    for (i = 0; i < *c; i++) rebmvnorm->W_[i] = W[i];

    rebmvnorm->IniTheta_ = new CompnentDistribution(rebmvnorm);

    R_CHECK(NULL == rebmvnorm->IniTheta_, ERInformationCriterionMVNORM);

    rebmvnorm->nc_ = rebmvnorm->length_pdf_ = *length_pdf;

    rebmvnorm->length_Theta_ = *length_Theta;

    rebmvnorm->length_theta_ = (INT*)malloc(rebmvnorm->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmvnorm->length_theta_, ERInformationCriterionMVNORM);

    *Error = rebmvnorm->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != EOK, *Error);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmvnorm->IniTheta_->pdf_[i] = pfNormal;
        }
        else {
            R_CHECK(1, ERInformationCriterionMVNORM);
        }
    }

    rebmvnorm->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == rebmvnorm->MixTheta_, ERInformationCriterionMVNORM);

    for (i = 0; i < *c; i++) {
        rebmvnorm->MixTheta_[i] = new CompnentDistribution(rebmvnorm);

        R_CHECK(NULL == rebmvnorm->MixTheta_[i], ERInformationCriterionMVNORM);

        *Error = rebmvnorm->MixTheta_[i]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        R_CHECK(*Error != EOK, *Error);
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

    rebmvnorm->Y_ = (FLOAT**)malloc(rebmvnorm->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmvnorm->Y_, ERInformationCriterionMVNORM);

    for (i = 0; i < rebmvnorm->length_pdf_; i++) {
        rebmvnorm->Y_[i] = (FLOAT*)malloc(rebmvnorm->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmvnorm->Y_[i], ERInformationCriterionMVNORM);
    }

    i = 0;

    for (j = 0; j < rebmvnorm->length_pdf_; j++) {
        for (l = 0; l < rebmvnorm->nr_; l++) {
            rebmvnorm->Y_[j][l] = x[i]; i++;
        }
    }

    rebmvnorm->cmax_ = *c;

    for (i = 0; i < *c; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        R_CHECK(*Error != EOK, *Error);
    }

    *Error = rebmvnorm->InformationCriterion(*c,
                                             rebmvnorm->W_,
                                             rebmvnorm->MixTheta_,
                                             IC,
                                             logL,
                                             M,
                                             D);

    R_CHECK(*Error != EOK, *Error);

EEXIT:

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
} // RInformationCriterionMVNORM

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
                              double *PSS,          // Pairwise similarity scores.
                              INT    *Error)        // Error code.
{
    Rebmvnorm *rebmvnorm = NULL;
    INT       i;

    *Error = EOK;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, ERCombineComponentsMVNORM);

    *Error = rebmvnorm->Set(NULL,         // Preprocessing type.
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

    R_CHECK(*Error != EOK, *Error);

    for (i = 0; i < rebmvnorm->cmax_; i++) {
        *Error = Cholinvdet(rebmvnorm->length_pdf_, rebmvnorm->MixTheta_[i]->Theta_[1], rebmvnorm->MixTheta_[i]->Theta_[2], rebmvnorm->MixTheta_[i]->Theta_[3]);

        R_CHECK(*Error != EOK, *Error);
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
            ED,
            PSS);

        R_CHECK(*Error != EOK, *Error);
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
            ED,
            PSS);

        R_CHECK(*Error != EOK, *Error);
    }
    else {
        R_CHECK(1, ERCombineComponentsMVNORM);
    }
/// End

EEXIT: 
    
    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
} // RCombineComponentsMVNORM

void RTvtNormalPdf(INT *n, double *x, double *y, double *Mean, double *Sigma, double *f)
{
    FLOAT Adet, Ainv[4], xi, yi, z;
    INT   j;

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

void RMvtNormalPdf(INT *n, double *X, INT *d, double *Mean, double *Sigma, double *f, INT *Error)
{
    FLOAT logAdet, *Ainv, y, yi, yk;
    INT   i, j, k;

    *Error = EOK;

    Ainv = (FLOAT*)malloc(*d * *d * sizeof(FLOAT));

    R_CHECK(NULL == Ainv, ERMvtNormalPdf);

    *Error = Cholinvdet(*d, Sigma, Ainv, &logAdet);

    if (*Error != EOK) {
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

EEXIT:

    if (Ainv) free(Ainv);

    R_RETURN(*Error);
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
    Rebmvnorm *rebmvnorm = NULL;
    INT       A[4], i = 0, j = 0, l = 0, length_Theta = 4;

    *Error = EOK;

    A[0] = *d; A[1] = A[2] = (*d) * (*d); A[3] = 1;

    rebmvnorm = new Rebmvnorm;

    R_CHECK(NULL == rebmvnorm, EREMMVNORM);

    *Error = rebmvnorm->Set(NULL,              // Preprocessing type.
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

    R_CHECK(NULL == rebmvnorm->MixTheta_, EREMMVNORM);

    for (l = 0; l < rebmvnorm->cmax_; l++) {
        rebmvnorm->MixTheta_[l] = new CompnentDistribution(rebmvnorm);

        R_CHECK(NULL == rebmvnorm->MixTheta_[l], EREMMVNORM);

        *Error = rebmvnorm->MixTheta_[l]->Realloc(rebmvnorm->length_pdf_, rebmvnorm->length_Theta_, rebmvnorm->length_theta_);

        R_CHECK(*Error != EOK, *Error);

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

    R_CHECK(NULL == rebmvnorm->W_, EREMMVNORM);

    for (l = 0; l < *c; l++) {
        rebmvnorm->W_[l] = W[l];

        *Error = Cholinvdet(*d, rebmvnorm->MixTheta_[l]->Theta_[1], rebmvnorm->MixTheta_[l]->Theta_[2], rebmvnorm->MixTheta_[l]->Theta_[3]);

        R_CHECK(*Error != EOK, *Error);

    }

    *Error = rebmvnorm->EMInitialize();

    R_CHECK(*Error != EOK, *Error);

    if (*EMMerge) {
        rebmvnorm->EM_->merge_ = merge_naive;
    }

    *Error = rebmvnorm->EMRun(c, rebmvnorm->W_, rebmvnorm->MixTheta_);

    R_CHECK(*Error != EOK, *Error);

    *Error = rebmvnorm->EM_->LogLikelihood(*c, rebmvnorm->W_, rebmvnorm->MixTheta_, summary_logL);

    R_CHECK(*Error != EOK, *Error);

    *Error = rebmvnorm->DegreesOffreedom(*c, rebmvnorm->MixTheta_, summary_M);

    R_CHECK(*Error != EOK, *Error);

    rebmvnorm->summary_.c = *c;

    *Error = rebmvnorm->Get(n_iter, // Number of iterations for optimal case.
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
                            NULL,   // Length of opt_c, opt_IC, opt_logL, opt_Dmin and opt_D.
                            NULL,   // Numbers of components for optimal v or for optimal k.
                            NULL,   // Information criteria for optimal v or for optimal k.
                            NULL,   // Log-likelihoods for optimal v or for optimal k.
                            NULL,   // Dmin for optimal v or for optimal k.
                            NULL,   // Totals of positive relative deviations for optimal v or for optimal k.
                            NULL,   // Length of all_K and all_IC.
                            NULL,   // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                            NULL);  // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.

    R_CHECK(*Error != EOK, *Error);

EEXIT:

    if (rebmvnorm) delete rebmvnorm;

    R_RETURN(*Error);
}
/// End
}
