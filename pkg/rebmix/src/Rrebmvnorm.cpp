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

} // extern "C"