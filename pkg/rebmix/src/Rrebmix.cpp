#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "rngmixf.h"
#include "rebmixf.h"

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

extern "C" {

// Runs RNGMIX in R.

void RRNGMIX(int    *IDum,         // Random seed.
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
    Rngmix *rngmix;
    int    i, j, k, l;

    rngmix = new Rngmix;

    *Error = NULL == rngmix; if (*Error) goto E0;

    rngmix->IDum_ = *IDum;
    rngmix->length_pdf_ = *d;
    rngmix->c_ = *c;

    rngmix->N_ = (int*)malloc(rngmix->c_ * sizeof(int));

    *Error = NULL == rngmix->N_; if (*Error) goto E0;

    for (i = 0; i < rngmix->c_; i++) rngmix->N_[i] = N[i];

    rngmix->length_pdf_ = *length_pdf;

    rngmix->length_Theta_ = *length_Theta;

    rngmix->length_theta_ = (int*)malloc(rngmix->length_Theta_ * sizeof(int));

    *Error = NULL == rngmix->length_theta_; if (*Error) goto E0;

    for (i = 0; i < rngmix->length_Theta_; i++) {
        rngmix->length_theta_[i] = length_theta[i];
    }

    rngmix->IniTheta_ = new CompnentDistribution(rngmix);

    *Error = NULL == rngmix->IniTheta_; if (*Error) goto E0;

    *Error = rngmix->IniTheta_->Realloc(rngmix->length_pdf_, rngmix->length_Theta_, rngmix->length_theta_);

    if (*Error) goto E0;

    for (i = 0; i < rngmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rngmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rngmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rngmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rngmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rngmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rngmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rngmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rngmix->MixTheta_ = new CompnentDistribution* [rngmix->c_];

    *Error = NULL == rngmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < rngmix->c_; i++) {
        rngmix->MixTheta_[i] = new CompnentDistribution(rngmix);

        *Error = NULL == rngmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rngmix->MixTheta_[i]->Realloc(rngmix->length_pdf_, rngmix->length_Theta_, rngmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < rngmix->c_; i++) {
        for (j = 0; j < rngmix->length_pdf_; j++) {
            rngmix->MixTheta_[i]->pdf_[j] = rngmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rngmix->length_Theta_; j++) if (rngmix->IniTheta_->Theta_[j]) {
        for (k = 0; k < rngmix->c_; k++) {
            for (l = 0; l < rngmix->length_theta_[j]; l++) {
                rngmix->MixTheta_[k]->Theta_[j][l] = Theta[i];

                i++;
            }
        }
    }

    *Error = rngmix->RNGMIX();

    if (*Error) goto E0;

    *n = rngmix->n_; i = 0;

    for (j = 0; j < rngmix->length_pdf_; j++) {
        for (k = 0; k < rngmix->n_; k++) {
            Y[i] = rngmix->Y_[k][j]; i++;
        }
    }

E0: if (rngmix) delete rngmix;
} // RRNGMIX

// Runs REBMIX in R.

void RREBMIX(char   **Preprocessing, // Preprocessing type.
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
    Rebmix *rebmix = NULL;
    int    i, j, k, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    if (!strcmp(Preprocessing[0], "histogram")) {
        rebmix->Preprocessing_ = poHistogram;
    }
    else
    if (!strcmp(Preprocessing[0], "Parzen window")) {
        rebmix->Preprocessing_ = poParzenWindow;
    }
    else
    if (!strcmp(Preprocessing[0], "k-nearest neighbour")) {
        rebmix->Preprocessing_ = poKNearestNeighbour;
    }
    else {
        *Error = 1; goto E0;
    }

    rebmix->cmax_ = *cmax;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC; 
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->length_pdf_ = *d;

    rebmix->Variables_ = (VariablesType_e*)malloc(rebmix->length_pdf_ * sizeof(VariablesType_e));

    *Error = NULL == rebmix->Variables_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(Variables[i], "continuous")) {
            rebmix->Variables_[i] = vtContinuous;
        }
        else
        if (!strcmp(Variables[i], "discrete")) {
            rebmix->Variables_[i] = vtDiscrete;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (k = 0; k < rebmix->length_theta_[j]; k++) {
            rebmix->IniTheta_->Theta_[j][k] = Theta[i];

            i++;
        }
    }

    rebmix->length_K_ = *length_K;

    rebmix->K_ = (int*)malloc(rebmix->length_K_ * sizeof(int));

    *Error = NULL == rebmix->K_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_K_; i++) {
        rebmix->K_[i] = K[i];
    }

    if (*length_y0 > 0) {
        rebmix->y0_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->y0_; if (*Error) goto E0;

        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->y0_[i] = y0[i];
        }
    }
    else {
        rebmix->y0_ = NULL;
    }

    if (*length_ymin > 0) {
        rebmix->ymin_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->ymin_; if (*Error) goto E0;

        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = ymin[i];
        }
    }
    else {
        rebmix->ymin_ = NULL;
    }

    if (*length_ymax > 0) {
        rebmix->ymax_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->ymax_; if (*Error) goto E0;

        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = ymax[i];
        }
    }
    else {
        rebmix->ymax_ = NULL;
    }

    rebmix->ar_ = *ar;

    if (!strcmp(Restraints[0], "rigid")) {
        rebmix->Restraints_ = rtRigid;
    }
    else
    if (!strcmp(Restraints[0], "loose")) {
        rebmix->Restraints_ = rtLoose;
    }
    else {
        *Error = 1; goto E0;
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[l][j] = Y[i]; i++;
        }
    }

    *Error = rebmix->REBMIX();

    if (*Error) goto E0;

    *summary_k = rebmix->summary_.k;

    if (rebmix->summary_.h) for (i = 0; i < rebmix->length_pdf_; i++) {
        summary_h[i] = rebmix->summary_.h[i];
    }

    if (rebmix->summary_.y0) for (i = 0; i < rebmix->length_pdf_; i++) {
        summary_y0[i] = rebmix->summary_.y0[i];
    }

    *summary_IC = rebmix->summary_.IC;
    *summary_logL = rebmix->summary_.logL;
    *summary_M = rebmix->summary_.M;
    *summary_c = rebmix->summary_.c;

    for (j = 0; j < rebmix->summary_.c; j++) {
        W[j] = rebmix->W_[j];
    }

    i = 0;

    for (j = 0; j < rebmix->summary_.c; j++) {
        for (l = 0; l < rebmix->length_theta_[0]; l++) {
            theta1[i] = rebmix->MixTheta_[j]->Theta_[0][l];

            i++;
        }
    }

    i = 0;

    for (j = 0; j < rebmix->summary_.c; j++) {
        for (l = 0; l < rebmix->length_theta_[1]; l++) {
            theta2[i] = rebmix->MixTheta_[j]->Theta_[1][l];

            i++;
        }
    }

    *opt_length = rebmix->opt_length_;

    for (i = 0; i < rebmix->opt_length_; i++) {
        opt_c[i] = rebmix->opt_c_[i];
        opt_IC[i] = rebmix->opt_IC_[i];
        opt_logL[i] = rebmix->opt_logL_[i];
        opt_D[i] = rebmix->opt_D_[i];
    }

    i = 0;

    for (j = 0; j < rebmix->all_length_; j++) if (rebmix->all_K_[j]) {
        all_K[i] = rebmix->all_K_[j]; all_IC[i] = rebmix->all_IC_[j]; 

        i++;
    }

    *all_length = i;

E0: if (rebmix) delete rebmix;
} // RREBMIX

// Returns k-nearest neighbour empirical densities in R.

void RdensKNearestNeighbourXY(int    *n,     // Total number of independent observations.
                              double *x,     // Pointer to the input array x.
                              double *y,     // Pointer to the input array y.
                              double *p,     // Pointer to the output array p.
                              int    *k,     // k-nearest neighbours.
                              double *hx,    // Normalizing vector.
                              double *hy,    // Normalizing vector.
                              int    *Error) // Error code.
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, C;
    int   i, j, K, l, m, q;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1; 
    
    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * Pi * (*hx) * (*hy));

    for (i = 0; i < *n; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (x[i] - x[j]) / (*hx); Dc = R * R;
            R = (y[i] - y[j]) / (*hy); Dc += R * R;

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < K; l++) {
                if (Dc < Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != K - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = (FLOAT)sqrt(Dk[K - 1]);

        if (q >= K) R *= (FLOAT)sqrt((K + (FLOAT)1.0) / (q + (FLOAT)2.0));

        p[i] = C / (R * R);
    }

E0: if (Dk) free(Dk);
} // RdensKNearestNeighbourXY 

// Returns Parzen window empirical densities in R.

void RdensParzenWindowXY(int    *n,     // Total number of independent observations.
                         double *x,     // Pointer to the input array x.
                         double *y,     // Pointer to the input array y.
                         double *p,     // Pointer to the output array p.
                         double *hx,    // Side of the hypersquare.
                         double *hy,    // Side of the hypersquare.
                         int    *Error) // Error code.
{
    int   i, j;
    FLOAT C, rx, ry;

    *Error = *n < 1; if (*Error) return;

    C = (FLOAT)1.0 / (*hx) / (*hy) / (*n); rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    for (i = 0; i < *n; i++) {
        p[i] = (FLOAT)0.0;
    }

    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if (((FLOAT)fabs(x[j] - x[i]) <= rx) && ((FLOAT)fabs(y[j] - y[i]) <= ry)) {
                p[i] += C; if (i != j) p[j] += C;
            }
        }
    }
} // RdensParzenWindowXY 

// Returns histogram empirical densities in R.

void RdensHistogramXY(int    *k,     // Total number of bins.
                      int    *n,     // Total number of independent observations.
                      double *x,     // Pointer to the input array x.
                      double *y,     // Pointer to the input array y.
                      double *p,     // Pointer to the output array p.
                      double *x0,    // Origin.
                      double *y0,    // Origin.
                      double *hx,    // Side of the hypersquare.
                      double *hy,    // Side of the hypersquare.
                      int    *cx,    // If x discrete then cx = 1 else cx = 0.
                      int    *cy,    // If y discrete then cy = 1 else cy = 0.
                      char   **px,   // Parametric family type.
                      char   **py,   // Parametric family type.
                      int    *Error) // Error code.
{
    int                    i, j;
    ParametricFamilyType_e pdfx, pdfy;
    FLOAT                  C, rx, ry;

    *Error = *n < 1; if (*Error) return;

    if (!strcmp(px[0], "normal")) {
        pdfx = pfNormal;
    }
    else
    if (!strcmp(px[0], "lognormal")) {
        pdfx = pfLognormal;
    }
    else
    if (!strcmp(px[0], "Weibull")) {
        pdfx = pfWeibull;
    }
    else
    if (!strcmp(px[0], "gamma")) {
        pdfx = pfGamma;
    }
    else
    if (!strcmp(px[0], "binomial")) {
        pdfx = pfBinomial;
    }
    else
    if (!strcmp(px[0], "Poisson")) {
        pdfx = pfPoisson;
    }
    else
    if (!strcmp(px[0], "Dirac")) {
        pdfx = pfDirac;
    }
    else {
      *Error = 1; return;
    }

    if (!strcmp(py[0], "normal")) {
        pdfy = pfNormal;
    }
    else
    if (!strcmp(py[0], "lognormal")) {
        pdfy = pfLognormal;
    }
    else
    if (!strcmp(py[0], "Weibull")) {
        pdfy = pfWeibull;
    }
    else
    if (!strcmp(py[0], "gamma")) {
        pdfy = pfGamma;
    }
    else
    if (!strcmp(py[0], "binomial")) {
        pdfy = pfBinomial;
    }
    else
    if (!strcmp(py[0], "Poisson")) {
        pdfy = pfPoisson;
    }
    else
    if (!strcmp(py[0], "Dirac")) {
        pdfy = pfDirac;
    }
    else {
      *Error = 1; return;
    }

    C = (FLOAT)1.0 / (*hx) / (*hy) / (*n); rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5); 
        
        x[*k] = (*x0) + j * (*hx);

        switch (pdfx) {
        default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (x[*k] <= FLOAT_MIN) x[*k] += (*hx);
        }

        j = (int)floor((y[i] - (*y0)) / (*hy) + (FLOAT)0.5); 

        y[*k] = (*y0) + j * (*hy);

        switch (pdfy) {
        default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (y[*k] <= FLOAT_MIN) y[*k] += (*hy);
        }

        for (j = 0; j < *k; j++) {
            if (((FLOAT)fabs(x[j] - x[*k]) > rx) || ((FLOAT)fabs(y[j] - y[*k]) > ry)) goto S0;
                
            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++; 
S1:;}
} // RdensHistogramXY 

// Returns k-nearest neighbour empirical densities in R.

void RdensKNearestNeighbourX(int    *n,     // Total number of independent observations.
                             double *x,     // Pointer to the input array x.
                             double *p,     // Pointer to the output array p.
                             int    *k,     // k-nearest neighbours.
                             double *hx,    // Normalizing vector.
                             int    *Error) // Error code.
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, C;
    int   i, j, K, l, m, q;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1; 
    
    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * (FLOAT)2.0 * (*hx));

    for (i = 0; i < *n; i++) {
        Dk[0] = FLOAT_MAX; q = 0;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (FLOAT)fabs((x[i] - x[j]) / (*hx)); Dc = R;

            q += Dc <= FLOAT_MIN;

            for (l = 0; l < K; l++) {
                if (Dc < Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    if ((Dc > FLOAT_MIN) || (l != K - 1)) Dk[l] = Dc;

                    break;
                }
            }
        }

        R = Dk[K - 1];

        if (q >= K) R *= (FLOAT)((K + (FLOAT)1.0) / (q + (FLOAT)2.0));

        p[i] = C / R;
    }

E0: if (Dk) free(Dk);
} // RdensKNearestNeighbourX 

// Returns Parzen window empirical densities in R.

void RdensParzenWindowX(int    *n,     // Total number of independent observations.
                        double *x,     // Pointer to the input array x.
                        double *p,     // Pointer to the output array p.
                        double *hx,    // Side of the hypersquare.
                        int    *Error) // Error code.
{
    int   i, j;
    FLOAT C, rx;

    *Error = *n < 1; if (*Error) return;

    C = (FLOAT)1.0 / (*hx) / (*n); rx = (FLOAT)0.5 * (*hx);

    for (i = 0; i < *n; i++) {
        p[i] = (FLOAT)0.0;
    }
    
    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if ((FLOAT)fabs(x[j] - x[i]) <= rx) {
                p[i] += C; if (i != j) p[j] += C;
            }
        }
    }
} // RdensParzenWindowX 

// Returns histogram empirical densities in R.

void RdensHistogramX(int    *k,     // Total number of bins.
                     int    *n,     // Total number of independent observations.
                     double *x,     // Pointer to the input array x.
                     double *p,     // Pointer to the output array p.
                     double *x0,    // Origin.
                     double *hx,    // Side of the hypersquare.
                     int    *cx,    // If x discrete then cx = 1 else cx = 0.
                     char   **px,   // Parametric family type.
                     int    *Error) // Error code.
{
    int                    i, j;
    ParametricFamilyType_e pdfx;
    FLOAT                  C, rx;

    *Error = *n < 1; if (*Error) return;

    if (!strcmp(px[0], "normal")) {
        pdfx = pfNormal;
    }
    else
    if (!strcmp(px[0], "lognormal")) {
        pdfx = pfLognormal;
    }
    else
    if (!strcmp(px[0], "Weibull")) {
        pdfx = pfWeibull;
    }
    else
    if (!strcmp(px[0], "gamma")) {
        pdfx = pfGamma;
    }
    else
    if (!strcmp(px[0], "binomial")) {
        pdfx = pfBinomial;
    }
    else
    if (!strcmp(px[0], "Poisson")) {
        pdfx = pfPoisson;
    }
    else
    if (!strcmp(px[0], "Dirac")) {
        pdfx = pfDirac;
    }
    else {
      *Error = 1; return;
    }
    
    C = (FLOAT)1.0 / (*hx) / (*n); rx = (FLOAT)0.5 * (*hx);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5); 
        
        x[*k] = (*x0) + j * (*hx);

        switch (pdfx) {
        default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (x[*k] <= FLOAT_MIN) x[*k] += (*hx);
        }

        for (j = 0; j < *k; j++) {
            if ((FLOAT)fabs(x[j] - x[*k]) > rx) goto S0;
                
            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++; 
S1:;}
} // RdensHistogramX 

// Returns classified observations in R.

void RCLSMIX(int    *n,      // Total number of independent observations.
             double *X,      // Pointer to the input array X.
             int    *s,      // Number of classes.
             int    *o,      // Number of input REBMIX objects.
             int    *d,      // Number of independent random variables in REBMIX objects.
             int    *c,      // Number of components.
             double *W,      // Component weights.
             char   **pdf,   // Component parameters.
             double *theta1, // Component parameters.
             double *theta2, // Component parameters.
             double *P,      // Prior probabilities.
             int    *Z,      // Pointer to the output array Z.
             int    *Error)  // Error code.
{
    Rebmix               *rebmix = NULL;
    int                  **C = NULL; 
    int                  A[2];
    FLOAT                ***Q = NULL;
    FLOAT                *Y = NULL;
    CompnentDistribution ****Theta = NULL; 
    FLOAT                CmpDist, MixDist, MaxMixDist;
    int                  i, j, k, l, m;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

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
                Theta[j][k][l] = new CompnentDistribution(rebmix);

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
                    else
                    if (!strcmp(pdf[i], "lognormal")) {
                        Theta[j][k][l]->pdf_[m] = pfLognormal;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "Weibull")) {
                        Theta[j][k][l]->pdf_[m] = pfWeibull;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "gamma")) {
                        Theta[j][k][l]->pdf_[m] = pfGamma;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "binomial")) {
                        Theta[j][k][l]->pdf_[m] = pfBinomial;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "Poisson")) {
                        Theta[j][k][l]->pdf_[m] = pfPoisson;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                    }
                    else
                    if (!strcmp(pdf[i], "Dirac")) {
                        Theta[j][k][l]->pdf_[m] = pfDirac;

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

    i = d[0]; for (j = 1; j < *o; j++) if (d[j] > i) i = d[j]; 

    Y = (FLOAT*)malloc(i * sizeof(FLOAT));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < *n; i++) {
        Z[i] = 1; MaxMixDist = (FLOAT)0.0; 
         
        for (j = 0; j < *s; j++) {
            k = 0; MixDist = (FLOAT)1.0;
            
            for (l = 0; l < *o; l++) {
                for (m = 0; m < d[l]; m++) {
                    Y[m] = X[i + (*n) * (m + k)];
                }

                *Error = rebmix->MixtureDist(Y, C[j][l], Q[j][l], Theta[j][l], &CmpDist);

                if (*Error) goto E0;

                k += d[l]; MixDist *= CmpDist; 
            }

            MixDist *= P[j];

            if (MixDist > MaxMixDist) {
                Z[i] = j + 1; MaxMixDist = MixDist;
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

    if (rebmix) delete rebmix;
} // RCLSMIX 

// Returns clustered observations in R.

void RCLRMIX(int    *n,      // Total number of independent observations.
             double *X,      // Pointer to the input array X.
             int    *d,      // Number of independent random variables.
             int    *c,      // Number of components.
             double *W,      // Component weights.
             char   **pdf,   // Component parameters.
             double *theta1, // Component parameters.
             double *theta2, // Component parameters.
             double *Z,      // Pointer to the output array Z.
             int    *Error)  // Error code.
{
    Rebmix               *rebmix = NULL; 
    FLOAT                *Y = NULL;
    int                  A[2];
    CompnentDistribution **Theta = NULL; 
    FLOAT                CmpDist, MaxCmpDist;
    int                  i, j, k;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    Theta = new CompnentDistribution* [*c];

    *Error = NULL == Theta; if (*Error) goto E0;

    i = 0;

    for (j = 0; j < *c; j++) {
        Theta[j] = new CompnentDistribution(rebmix);

        *Error = NULL == Theta[j]; if (*Error) goto E0;

        A[0] = A[1] = *d;

        *Error = Theta[j]->Realloc(*d, 2, A);

        if (*Error) goto E0;

        for (k = 0; k < *d; k++) {
            if (!strcmp(pdf[i], "normal")) {
                Theta[j]->pdf_[k] = pfNormal;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "lognormal")) {
                Theta[j]->pdf_[k] = pfLognormal;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "Weibull")) {
                Theta[j]->pdf_[k] = pfWeibull;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "gamma")) {
                Theta[j]->pdf_[k] = pfGamma;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "binomial")) {
                Theta[j]->pdf_[k] = pfBinomial;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "Poisson")) {
                Theta[j]->pdf_[k] = pfPoisson;

                Theta[j]->Theta_[0][k] = theta1[i];
            }
            else
            if (!strcmp(pdf[i], "Dirac")) {
                Theta[j]->pdf_[k] = pfDirac;

                Theta[j]->Theta_[0][k] = theta1[i];
            }
            else {
                *Error = 1; goto E0;
            }

            i++;
        }
    }

    Y = (FLOAT*)malloc(*d * sizeof(FLOAT));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *d; j++) {
            Y[j] = X[i + (*n) * j];
        }

        Z[i] = (FLOAT)0.0; MaxCmpDist = (FLOAT)0.0;
         
        for (j = 0; j < *c; j++) {
            rebmix->length_pdf_ = *d;

            *Error = rebmix->ComponentDist(Y, Theta[j], &CmpDist);

            if (*Error) goto E0;

            CmpDist *= W[j];

            if (CmpDist > MaxCmpDist) {
                Z[i] = (FLOAT)j; MaxCmpDist = CmpDist;
            }
        }
    }

E0: if (Y) free(Y);

    if (Theta) {
        for (i = 0; i < *c; i++) {
            if (Theta[i]) delete Theta[i];
        }

        delete Theta;
    }

    if (rebmix) delete rebmix;
} // RCLRMIX

void RPreprocessingKNN(int    *k,     // k-nearest neighbours.
                       double *h,     // Normalizing vector.
                       int    *n,     // Total number of independent observations.
                       int    *d,     // Number of independent random variables.
                       double *x,     // Pointer to the input array x.
                       double *y,     // Pointer to the output array y.
                       int    *Error) // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->n_ = *n;
    rebmix->length_pdf_ = *d;

    Y = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        Y[i] = (FLOAT*)malloc((rebmix->length_pdf_ + 3) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            Y[l][j] = x[i]; i++;
        }
    }

    *Error = rebmix->PreprocessingKNN(*k, h, Y);

    if (*Error) goto E0;

    i = 0;
    
    for (j = 0; j < rebmix->length_pdf_ + 3; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            y[i] = Y[l][j]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmix->n_; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    if (rebmix) delete rebmix;
} // RPreprocessingKNN 

void RPreprocessingPW(double *h,     // Sides of the hypersquare.
                      int    *n,     // Total number of independent observations.
                      int    *d,     // Number of independent random variables.
                      double *x,     // Pointer to the input array x.
                      double *y,     // Pointer to the output array y.
                      int    *Error) // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL; 
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->n_ = *n;
    rebmix->length_pdf_ = *d;

    Y = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        Y[i] = (FLOAT*)malloc((rebmix->length_pdf_ + 2) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            Y[l][j] = x[i]; i++;
        }
    }

    *Error = rebmix->PreprocessingPW(h, Y); 

    if (*Error) goto E0;
    
    i = 0;
    
    for (j = 0; j < rebmix->length_pdf_ + 2; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            y[i] = Y[l][j]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmix->n_; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    if (rebmix) delete rebmix;
} // RPreprocessingPW 

void RPreprocessingH(double *h,          // Sides of the hypersquare.
                     double *y0,         // Origins.
                     int    *length_pdf, // Length of pdf.
                     char   **pdf,       // Parametric family types.
                     int    *k,          // Total number of bins.
                     int    *n,          // Total number of independent observations.
                     int    *d,          // Number of independent random variables.
                     double *x,          // Pointer to the input array x.
                     double *y,          // Pointer to the output array y.
                     int    *Error)      // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    int    i, j, l;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->length_pdf_ = *d;

    rebmix->length_pdf_ = *length_pdf;

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = rebmix->IniTheta_->Realloc(rebmix->length_pdf_, 0, NULL);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->IniTheta_->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal")) {
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        }
        else
        if (!strcmp(pdf[i], "lognormal")) {
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        }
        else
        if (!strcmp(pdf[i], "Weibull")) {
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        }
        else
        if (!strcmp(pdf[i], "gamma")) {
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        }
        else
        if (!strcmp(pdf[i], "binomial")) {
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        }
        else
        if (!strcmp(pdf[i], "Poisson")) {
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        }
        else
        if (!strcmp(pdf[i], "Dirac")) {
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        }
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[l][j] = x[i]; i++;
        }
    }
    
    Y = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        Y[i] = (FLOAT*)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmix->PreprocessingH(h, y0, k, Y);

    if (*Error) goto E0;

    i = 0;
    
    for (j = 0; j < rebmix->length_pdf_ + 1; j++) {
        for (l = 0; l < *k; l++) {
            y[i] = Y[l][j]; i++;
        }
    }

E0: if (Y) {
        for (i = 0; i < rebmix->n_; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    if (rebmix) delete rebmix;
} // RPreprocessingH

void RInformationCriterionKNN(double *h,            // Sides of the hypersquare.
                              int    *k,            // Total number of bins.
                              int    *d,            // Number of independent random variables.
                              char   **Criterion,   // Information criterion type.
                              int    *c,            // Number of components.
                              double *W,            // Component weights.
                              int    *length_pdf,   // Length of pdf.
                              int    *length_Theta, // Length of Theta.
                              int    *length_theta, // Length of Theta[i].
                              char   **pdf,         // Parametric family types.
                              double *Theta,        // Component parameters.
                              int    *n,            // Number of observations.
                              double *x,            // Dataset.
                              double *IC,           // Information criterion.
                              double *logL,         // log-likelihood.
                              int    *M,            // Degrees of freedom.
                              double *D,            // Total of positive relative deviations.
                              int    *Error)        // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    int    i, j, l, m;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->length_pdf_ = *d;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC;
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmix->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_Theta_; i++) {
        rebmix->length_theta_[i] = length_theta[i];
    }

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal"))
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        else
        if (!strcmp(pdf[i], "lognormal"))
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        else
        if (!strcmp(pdf[i], "Weibull"))
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        else
        if (!strcmp(pdf[i], "gamma"))
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        else
        if (!strcmp(pdf[i], "binomial"))
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        else
        if (!strcmp(pdf[i], "Poisson"))
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        else
        if (!strcmp(pdf[i], "Dirac"))
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->MixTheta_ = new CompnentDistribution*[*c];

    *Error = NULL == rebmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        *Error = NULL == rebmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[l][j] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        Y[i] = (FLOAT*)malloc((rebmix->length_pdf_ + 3) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmix->PreprocessingKNN(*k, h, Y);

    if (*Error) goto E0;

    rebmix->cmax_ = *c;

    *Error = rebmix->InformationCriterionKNN(*k, 
                                             Y, 
                                             *c, 
                                             rebmix->W_, 
                                             rebmix->MixTheta_, 
                                             IC, 
                                             logL,
                                             M,
                                             D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmix->n_; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RInformationCriterionKNN 

void RInformationCriterionPW(double *h,            // Sides of the hypersquare.
                             int    *d,            // Number of independent random variables.
                             char   **Criterion,   // Information criterion type.
                             int    *c,            // Number of components.
                             double *W,            // Component weights.
                             int    *length_pdf,   // Length of pdf.
                             int    *length_Theta, // Length of Theta.
                             int    *length_theta, // Length of Theta[i].
                             char   **pdf,         // Parametric family types.
                             double *Theta,        // Component parameters.
                             int    *n,            // Number of observations.
                             double *x,            // Dataset.
                             double *IC,           // Information criterion.
                             double *logL,         // log-likelihood.
                             int    *M,            // Degrees of freedom.
                             double *D,            // Total of positive relative deviations.
                             int    *Error)        // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    FLOAT  V;
    int    i, j, l, m;

    rebmix = new Rebmix;

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->length_pdf_ = *d;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC;
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmix->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_Theta_; i++) {
        rebmix->length_theta_[i] = length_theta[i];
    }

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

    if (*Error) goto E0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal"))
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        else
        if (!strcmp(pdf[i], "lognormal"))
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        else
        if (!strcmp(pdf[i], "Weibull"))
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        else
        if (!strcmp(pdf[i], "gamma"))
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        else
        if (!strcmp(pdf[i], "binomial"))
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        else
        if (!strcmp(pdf[i], "Poisson"))
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        else
        if (!strcmp(pdf[i], "Dirac"))
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->MixTheta_ = new CompnentDistribution*[*c];

    *Error = NULL == rebmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        *Error = NULL == rebmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[l][j] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        Y[i] = (FLOAT*)malloc((rebmix->length_pdf_ + 2) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmix->PreprocessingPW(h, Y); 

    if (*Error) goto E0;
 
    rebmix->cmax_ = *c;

    V = (FLOAT)1.0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        V *= h[i];
    }

    *Error = rebmix->InformationCriterionPW(V, 
                                            Y, 
                                            *c, 
                                            rebmix->W_, 
                                            rebmix->MixTheta_, 
                                            IC,
                                            logL,
                                            M,
                                            D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmix->n_; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;
} // RInformationCriterionPW 

void RInformationCriterionH(double *h,            // Sides of the hypersquare.
                            double *y0,           // Origins.
                            int    *k,            // Total number of bins.
                            int    *d,            // Number of independent random variables.
                            char   **Criterion,   // Information criterion type.
                            int    *c,            // Number of components.
                            double *W,            // Component weights.
                            int    *length_pdf,   // Length of pdf.
                            int    *length_Theta, // Length of Theta.
                            int    *length_theta, // Length of Theta[i].
                            char   **pdf,         // Parametric family types.
                            double *Theta,        // Component parameters.
                            int    *n,            // Number of observations.
                            double *x,            // Dataset.
                            double *IC,           // Information criterion.
                            double *logL,         // log-likelihood.
                            int    *M,            // Degrees of freedom.
                            double *D,            // Total of positive relative deviations.
                            int    *Error)        // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    FLOAT  V;
    int    i, j, l, m;

    rebmix = new Rebmix;    

    *Error = NULL == rebmix; if (*Error) goto E0;

    rebmix->length_pdf_ = *d;

    if (!strcmp(Criterion[0], "AIC"))
        rebmix->Criterion_ = icAIC;
    else
    if (!strcmp(Criterion[0], "AIC3"))
        rebmix->Criterion_ = icAIC3;
    else
    if (!strcmp(Criterion[0], "AIC4"))
        rebmix->Criterion_ = icAIC4;
    else
    if (!strcmp(Criterion[0], "AICc"))
        rebmix->Criterion_ = icAICc;
    else
    if (!strcmp(Criterion[0], "BIC"))
        rebmix->Criterion_ = icBIC;
    else
    if (!strcmp(Criterion[0], "CAIC"))
        rebmix->Criterion_ = icCAIC;
    else
    if (!strcmp(Criterion[0], "HQC"))
        rebmix->Criterion_ = icHQC;
    else
    if (!strcmp(Criterion[0], "MDL2"))
        rebmix->Criterion_ = icMDL2;
    else
    if (!strcmp(Criterion[0], "MDL5"))
        rebmix->Criterion_ = icMDL5;
    else
    if (!strcmp(Criterion[0], "AWE"))
        rebmix->Criterion_ = icAWE;
    else
    if (!strcmp(Criterion[0], "CLC"))
        rebmix->Criterion_ = icCLC;
    else
    if (!strcmp(Criterion[0], "ICL"))
        rebmix->Criterion_ = icICL;
    else
    if (!strcmp(Criterion[0], "PC"))
        rebmix->Criterion_ = icPC;
    else
    if (!strcmp(Criterion[0], "ICL-BIC"))
        rebmix->Criterion_ = icICLBIC;
    else
    if (!strcmp(Criterion[0], "D"))
        rebmix->Criterion_ = icD;
    else
    if (!strcmp(Criterion[0], "SSE"))
        rebmix->Criterion_ = icSSE;
    else {
        *Error = 1; goto E0;
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == rebmix->W_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (int*)malloc(rebmix->length_Theta_ * sizeof(int));

    *Error = NULL == rebmix->length_theta_; if (*Error) goto E0;

    for (i = 0; i < rebmix->length_Theta_; i++) {
        rebmix->length_theta_[i] = length_theta[i];
    }

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    *Error = NULL == rebmix->IniTheta_; if (*Error) goto E0;

    *Error = rebmix->IniTheta_->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

    if (*Error) goto E0;
    
    for (i = 0; i < rebmix->length_pdf_; i++) {
        if (!strcmp(pdf[i], "normal"))
            rebmix->IniTheta_->pdf_[i] = pfNormal;
        else
        if (!strcmp(pdf[i], "lognormal"))
            rebmix->IniTheta_->pdf_[i] = pfLognormal;
        else
        if (!strcmp(pdf[i], "Weibull"))
            rebmix->IniTheta_->pdf_[i] = pfWeibull;
        else
        if (!strcmp(pdf[i], "gamma"))
            rebmix->IniTheta_->pdf_[i] = pfGamma;
        else
        if (!strcmp(pdf[i], "binomial"))
            rebmix->IniTheta_->pdf_[i] = pfBinomial;
        else
        if (!strcmp(pdf[i], "Poisson"))
            rebmix->IniTheta_->pdf_[i] = pfPoisson;
        else
        if (!strcmp(pdf[i], "Dirac"))
            rebmix->IniTheta_->pdf_[i] = pfDirac;
        else {
            *Error = 1; goto E0;
        }
    }

    rebmix->MixTheta_ = new CompnentDistribution* [*c];

    *Error = NULL == rebmix->MixTheta_; if (*Error) goto E0;

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        *Error = NULL == rebmix->MixTheta_[i]; if (*Error) goto E0;

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        if (*Error) goto E0;
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == rebmix->Y_; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

        *Error = NULL == rebmix->Y_[i]; if (*Error) goto E0;
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->n_; l++) {
            rebmix->Y_[l][j] = x[i]; i++;
        }
    }
    
    Y = (FLOAT**)malloc(rebmix->n_ * sizeof(FLOAT*));

    *Error = NULL == Y; if (*Error) goto E0;

    for (i = 0; i < rebmix->n_; i++) {
        Y[i] = (FLOAT*)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT));

        *Error = NULL == Y[i]; if (*Error) goto E0;
    }

    *Error = rebmix->PreprocessingH(h, y0, k, Y);

    if (*Error) goto E0;

    rebmix->cmax_ = *c;

    V = (FLOAT)1.0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        V *= h[i];
    }

    *Error = rebmix->InformationCriterionH(V, 
                                           *k, 
                                           Y, 
                                           *c, 
                                           rebmix->W_, 
                                           rebmix->MixTheta_, 
                                           IC,
                                           logL,
                                           M,
                                           D);

    if (*Error) goto E0;

E0: if (Y) {
        for (i = 0; i < rebmix->n_; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    if (rebmix) delete rebmix;
} // RInformationCriterionH

} // extern "C"