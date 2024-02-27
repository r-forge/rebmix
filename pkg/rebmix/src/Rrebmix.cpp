#include "rngmixf.h"

#include <math.h>

extern "C" {

// Runs RNGMIX in R.

void RRNGMIX(INT    *IDum,         // Random seed.
             INT    *d,            // Number of independent random variables.
             INT    *c,            // Number of components.
             INT    *N,            // Numbers of observations.
             INT    *length_pdf,   // Length of pdf.
             INT    *length_Theta, // Length of Theta.
             INT    *length_theta, // Length of Theta[i].
             char   **pdf,         // Parametric family types.
             double *Theta,        // Component parameters.
             INT    *n,            // Number of observations.
             double *Y,            // Dataset.
             INT    *Z,            // Component membership.
             INT    *Error)        // Error code.
{
    Rngmix *rngmix = NULL;
    INT    i, j, k, l;

    *Error = E_OK;

    rngmix = new Rngmix;

    R_CHECK(NULL == rngmix, E_MEM);

    rngmix->IDum_ = *IDum;
    rngmix->length_pdf_ = *d;
    rngmix->c_ = *c;

    rngmix->N_ = (INT*)malloc(rngmix->c_ * sizeof(INT));

    R_CHECK(NULL == rngmix->N_, E_MEM);

    for (i = 0; i < rngmix->c_; i++) rngmix->N_[i] = N[i];

    rngmix->IniTheta_ = new CompnentDistribution(rngmix);

    R_CHECK(NULL == rngmix->IniTheta_, E_MEM);

    rngmix->length_pdf_ = *length_pdf;

    rngmix->length_Theta_ = *length_Theta;

    rngmix->length_theta_ = (INT*)malloc(rngmix->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rngmix->length_theta_, E_MEM);

    *Error = rngmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != E_OK, *Error);

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
        if (!strcmp(pdf[i], "Gumbel")) {
            rngmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rngmix->IniTheta_->pdf_[i] = pfvonMises;
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
        else
        if (!strcmp(pdf[i], "uniform")) {
            rngmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            R_CHECK(1, E_ARG);
        }
    }

    rngmix->MixTheta_ = new CompnentDistribution* [(unsigned INT)rngmix->c_];

    R_CHECK(NULL == rngmix->MixTheta_, E_MEM);

    for (i = 0; i < rngmix->c_; i++) {
        rngmix->MixTheta_[i] = new CompnentDistribution(rngmix);

        R_CHECK(NULL == rngmix->MixTheta_[i], E_MEM);

        *Error = rngmix->MixTheta_[i]->Realloc(rngmix->length_pdf_, rngmix->length_Theta_, rngmix->length_theta_);

        R_CHECK(*Error != E_OK, *Error);
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

    R_CHECK(*Error != E_OK, *Error);

    *n = rngmix->n_; i = 0;

    for (j = 0; j < rngmix->length_pdf_; j++) {
        for (k = 0; k < rngmix->n_; k++) {
            Y[i] = rngmix->Y_[j][k]; i++;
        }
    }

    for (i = 0; i < rngmix->n_; i++) {
        Z[i] = rngmix->Z_[i];
    }

EEXIT:

    if (rngmix) delete rngmix;

    R_RETURN(*Error);
} // RRNGMIX

// Runs REBMIX in R.

void RREBMIX(char   **Preprocessing, // Preprocessing type.
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
             double *theta3,         // Component parameters.
             INT    *opt_length,     // Length of opt_c, opt_IC, opt_logL, opt_Dmin and opt_D.
             INT    *opt_c,          // Numbers of components for optimal v or for optimal k.
             double *opt_IC,         // Information criteria for optimal v or for optimal k.
             double *opt_logL,       // Log-likelihoods for optimal v or for optimal k.
             double *opt_Dmin,       // Dmin for optimal v or for optimal k.
             double *opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
             INT    *all_length,     // Length of all_K and all_IC.
             INT    *all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
             double *all_IC,         // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
             char   **EList)         // Error list.
{
    Rebmix *rebmix = NULL;
    INT    Error;

    Error = E_OK;

    rebmix = new Rebmix;

    E_CHECK(NULL == rebmix, E_MEM);

    Error = rebmix->Set(Preprocessing,     // Preprocessing type.
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

    E_CHECK(Error != E_OK, Error);

    Error = rebmix->REBMIX();

    E_CHECK(Error != E_OK, Error);

/// Panic Branislav
    Error = rebmix->Get(n_iter,         // Number of iterations for optimal case.
                        n_iter_sum,     // Number of iterations in whole run.
/// End
                        summary_k,      // Optimal v or optimal k.
                        summary_h,      // Optimal class widths of length d.
                        summary_y0,     // Optimal origins of length d.
                        summary_ymin,   // Optimal minimum observations of length d.
                        summary_ymax,   // Optimal maximum observations of length d.
                        summary_IC,     // Optimal information criterion.
                        summary_logL,   // Log-likelihood.
                        summary_M,      // Degrees of freedom.
                        summary_c,      // Optimal number of components.
                        W,              // Component weights.
                        theta1,         // Component parameters.
                        theta2,         // Component parameters.
                        theta3,         // Component parameters.
                        opt_length,     // Length of opt_c, opt_IC, opt_logL, opt_Dmin and opt_D.
                        opt_c,          // Numbers of components for optimal v or for optimal k.
                        opt_IC,         // Information criteria for optimal v or for optimal k.
                        opt_logL,       // Log-likelihoods for optimal v or for optimal k.
                        opt_Dmin,       // Dmin for optimal v or for optimal k.
                        opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
                        all_length,     // Length of all_K and all_IC.
                        all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                        all_IC);        // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
    
    E_CHECK(Error != E_OK, Error);

EEXIT:

    if (rebmix) delete rebmix;

    E_LIST(*EList);
} // RREBMIX

// Returns k-nearest neighbour empirical densities in R.

void RdensKNearestNeighbourXY(INT    *n,     // Total number of independent observations.
                              double *x,     // Pointer to the input array x.
                              double *y,     // Pointer to the input array y.
                              double *p,     // Pointer to the output array p.
                              INT    *k,     // k-nearest neighbours.
                              double *hx,    // Normalizing vector.
                              double *hy,    // Normalizing vector.
                              INT    *Error) // Error code.
{
    FLOAT C, Dc, *Dk = NULL, R;
    INT   i, j, K, l, m, q;

    *Error = E_OK;

    R_CHECK(*n < 1, E_ARG);

    K = *k; if (K > 1) K -= 1; else K = 1;

    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    R_CHECK(NULL == Dk, E_MEM);

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

EEXIT: 
    
    if (Dk) free(Dk);

    R_RETURN(*Error);
} // RdensKNearestNeighbourXY

// Returns kernel density estimation empirical densities in R.

void RdensKDEXY(INT    *n,     // Total number of independent observations.
                double *x,     // Pointer to the input array x.
                double *y,     // Pointer to the input array y.
                double *p,     // Pointer to the output array p.
                double *hx,    // Side of the hypersquare.
                double *hy,    // Side of the hypersquare.
                INT    *Error) // Error code.
{
    FLOAT C, rx, ry;
    INT   i, j;

    *Error = E_OK;

    R_CHECK(*n < 1, E_ARG);

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

EEXIT:

    R_RETURN(*Error);
} // RdensKDEXY

// Returns histogram empirical densities in R.

void RdensHistogramXY(INT    *k,     // Total number of bins.
                      INT    *n,     // Total number of independent observations.
                      double *x,     // Pointer to the input array x.
                      double *y,     // Pointer to the input array y.
                      double *p,     // Pointer to the output array p.
                      double *x0,    // Origin.
                      double *xmin,  // Minimum observation.
                      double *xmax,  // Maximim observation.
                      double *y0,    // Origin.
                      double *ymin,  // Minimum observation.
                      double *ymax,  // Maximim observation.
                      double *hx,    // Side of the hypersquare.
                      double *hy,    // Side of the hypersquare.
                      char   **px,   // Parametric family type.
                      char   **py,   // Parametric family type.
                      INT    *Error) // Error code.
{
    ParametricFamilyType_e pdfx, pdfy;
    FLOAT                  C, rx, ry;
    INT                    i, j;

    *Error = E_OK;

    R_CHECK(*n < 1, E_ARG);

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
    if (!strcmp(px[0], "Gumbel")) {
        pdfx = pfGumbel;
    }
    else
    if (!strcmp(px[0], "vonMises")) {
        pdfx = pfvonMises;
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
    else
    if (!strcmp(px[0], "uniform")) {
        pdfx = pfUniform;
    }
    else {
        R_CHECK(1, E_ARG);
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
    if (!strcmp(py[0], "Gumbel")) {
        pdfy = pfGumbel;
    }
    else
    if (!strcmp(py[0], "vonMises")) {
        pdfy = pfvonMises;
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
    else
    if (!strcmp(py[0], "uniform")) {
        pdfy = pfUniform;
    }
    else {
        R_CHECK(1, E_ARG);
    }

    C = (FLOAT)1.0 / (*hx) / (*hy) / (*n); rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (INT)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5);

        x[*k] = (*x0) + j * (*hx);

        if (x[*k] < *xmin) {
            x[*k] += (*hx);
        }
        else
        if (x[*k] > *xmax) {
            x[*k] -= (*hx);
        }

        switch (pdfx) {
        case pfNormal: case pfTNormal: case pfvonMises: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform: case pfGumbel: default:
            break;
        case pfLognormal: case pfWeibull: case pfGamma:
            if (x[*k] <= FLOAT_MIN) x[*k] += (*hx);
        }

        j = (INT)floor((y[i] - (*y0)) / (*hy) + (FLOAT)0.5);

        y[*k] = (*y0) + j * (*hy);

        if (y[*k] < *ymin) {
            y[*k] += (*hy);
        }
        else
        if (y[*k] > *ymax) {
            y[*k] -= (*hy);
        }

        switch (pdfy) {
        case pfNormal: case pfTNormal: case pfvonMises: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform: case pfGumbel: default:
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

EEXIT:

    R_RETURN(*Error);
} // RdensHistogramXY

// Returns histogram input empirical densities in R.

void RdensKXY(INT    *v,     // Total number of bins.
              double *x,     // Pointer to the input array x.
              double *y,     // Pointer to the input array y.
              double *k,     // Pointer to the input array k.
              double *p,     // Pointer to the output array p.
              double *hx,    // Side of the hypersquare.
              double *hy,    // Side of the hypersquare.
              INT    *Error) // Error code.
{
    FLOAT C, rx, ry;
    INT   i, j, n;

    *Error = E_OK;

    R_CHECK(*v < 1, E_ARG);

    rx = (FLOAT)0.5 * (*hx); ry = (FLOAT)0.5 * (*hy);

    n = 0;

    for (i = 0; i < *v; i++) {
        p[i] = k[i];

        j = i + 1;

        if (j < *v) do {
            if ((x[j] < x[i] + rx) && (x[j] > x[i] - rx) && (y[j] < y[i] + ry) && (y[j] > y[i] - ry)) {
                p[i] += k[j]; (*v)--; 
                
                x[j] = x[*v]; y[j] = y[*v]; k[j] = k[*v];
            }
            else {
                j++;
            }
        } while (j < *v);

        n += (INT)p[i];
    }

    C = (FLOAT)1.0 / (*hx) / (*hy) / (n);

    for (i = 0; i < *v; i++) p[i] *= C;

EEXIT:

    R_RETURN(*Error);
} // RdensKXY

// Returns k-nearest neighbour empirical densities in R.

void RdensKNearestNeighbourX(INT    *n,     // Total number of independent observations.
                             double *x,     // Pointer to the input array x.
                             double *p,     // Pointer to the output array p.
                             INT    *k,     // k-nearest neighbours.
                             double *hx,    // Normalizing vector.
                             INT    *Error) // Error code.
{
    FLOAT Dc, *Dk = NULL, R, C;
    INT   i, j, K, l, m, q;

    *Error = E_OK;

    R_CHECK(*n < 1, E_ARG);

    K = *k; if (K > 1) K -= 1; else K = 1;

    Dk = (FLOAT*)malloc(K * sizeof(FLOAT));

    R_CHECK(NULL == Dk, E_MEM);

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

EEXIT:  
    
    if (Dk) free(Dk);

    R_RETURN(*Error);
} // RdensKNearestNeighbourX

// Returns kernel density estimation empirical densities in R.

void RdensKDEX(INT    *n,     // Total number of independent observations.
               double *x,     // Pointer to the input array x.
               double *p,     // Pointer to the output array p.
               double *hx,    // Side of the hypersquare.
               INT    *Error) // Error code.
{
    FLOAT C, rx;
    INT   i, j;

    *Error = E_OK;

    R_CHECK(*n < 1, E_ARG);

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

EEXIT:

    R_RETURN(*Error);
} // RdensKDEX

// Returns histogram empirical densities in R.

void RdensHistogramX(INT    *k,     // Total number of bins.
                     INT    *n,     // Total number of independent observations.
                     double *x,     // Pointer to the input array x.
                     double *p,     // Pointer to the output array p.
                     double *x0,    // Origin.
                     double *xmin,  // Minimum observation.
                     double *xmax,  // Maximum observation.
                     double *hx,    // Side of the hypersquare.
                     char   **px,   // Parametric family type.
                     INT    *Error) // Error code.
{
    ParametricFamilyType_e pdfx;
    FLOAT                  C, rx;
    INT                    i, j;

    *Error = E_OK;

    R_CHECK(*n < 1, E_ARG);

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
    if (!strcmp(px[0], "Gumbel")) {
        pdfx = pfGumbel;
    }
    else
    if (!strcmp(px[0], "vonMises")) {
        pdfx = pfvonMises;
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
    else
    if (!strcmp(px[0], "uniform")) {
        pdfx = pfUniform;
    }
    else {
        R_CHECK(1, E_ARG);
    }

    C = (FLOAT)1.0 / (*hx) / (*n); rx = (FLOAT)0.5 * (*hx);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (INT)floor((x[i] - (*x0)) / (*hx) + (FLOAT)0.5);

        x[*k] = (*x0) + j * (*hx);

        if (x[*k] < *xmin) {
            x[*k] += (*hx);
        }
        else
        if (x[*k] > *xmax) {
            x[*k] -= (*hx);
        }

        switch (pdfx) {
        case pfNormal: case pfTNormal: case pfvonMises: case pfBinomial: case pfPoisson: case pfDirac: case pfUniform: case pfGumbel: default:
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

EEXIT:

    R_RETURN(*Error);
} // RdensHistogramX

// Returns histogram input empirical densities in R.

void RdensKX(INT    *v,     // Total number of bins.
             double *x,     // Pointer to the input array x.
             double *k,     // Pointer to the input array k.
             double *p,     // Pointer to the output array p.
             double *hx,    // Side of the hypersquare.
             INT    *Error) // Error code.
{
    FLOAT C, rx;
    INT   i, j, n;

    *Error = E_OK;

    R_CHECK(*v < 1, E_ARG);

    rx = (FLOAT)0.5 * (*hx);

    n = 0;

    for (i = 0; i < *v; i++) {
        p[i] = k[i];

        j = i + 1;

        if (j < *v) do {
            if ((x[j] < x[i] + rx) && (x[j] > x[i] - rx)) {
                p[i] += k[j]; (*v)--;

                x[j] = x[*v]; k[j] = k[*v];
            }
            else {
                j++;
            }
        } while (j < *v);

        n += (INT)p[i];
    }

    C = (FLOAT)1.0 / (*hx) / (n);

    for (i = 0; i < *v; i++) p[i] *= C;

EEXIT:

    R_RETURN(*Error);
} // RdensKX

// Returns classified observations in R.

void RCLSMIX(INT    *n,      // Total number of independent observations.
             double *X,      // Pointer to the input array X.
             INT    *s,      // Number of classes.
             INT    *o,      // Number of input REBMIX objects.
             INT    *d,      // Number of independent random variables in REBMIX objects.
             INT    *c,      // Number of components.
             double *W,      // Component weights.
             char   **pdf,   // Component parameters.
             double *theta1, // Component parameters.
             double *theta2, // Component parameters.
             double *theta3, // Component parameters.
             double *P,      // Prior probabilities.
             INT    *Z,      // Pointer to the output array Z.
             INT    *Error)  // Error code.
{
    Rebmix               *rebmix = NULL;
    CompnentDistribution ****Theta = NULL;
    FLOAT                CmpDist, MixDist, MaxMixDist, ***Q = NULL, **Y = NULL;
    INT                  A[3], **C = NULL, i, j, k, l, m, dmax = 0;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    C = (INT**)malloc(*s * sizeof(INT*));

    R_CHECK(NULL == C, E_MEM);

    i = 0;

    for (j = 0; j < *s; j++) {
        C[j] = (INT*)malloc(*o * sizeof(INT));

        R_CHECK(NULL == C[j], E_MEM);

        for (k = 0; k < *o; k++) {
            C[j][k] = c[i]; i++;
        }
    }

    Q = (FLOAT***)malloc(*s * sizeof(FLOAT**));

    R_CHECK(NULL == Q, E_MEM);

    i = 0;

    for (j = 0; j < *s; j++) {
        Q[j] = (FLOAT**)malloc(*o * sizeof(FLOAT*));

        R_CHECK(NULL == Q[j], E_MEM);

        for (k = 0; k < *o; k++) {
            Q[j][k] = (FLOAT*)malloc(C[j][k] * sizeof(FLOAT));

            R_CHECK(NULL == Q[j][k], E_MEM);

            for (l = 0; l < C[j][k]; l++) {
                Q[j][k][l] = W[i]; i++;
            }
        }
    }

    Theta = new CompnentDistribution*** [(unsigned INT)(*s)];

    R_CHECK(NULL == Theta, E_MEM);

    i = 0;

    for (j = 0; j < *s; j++) {
        Theta[j] = new CompnentDistribution** [(unsigned INT)(*o)];

        R_CHECK(NULL == Theta[j], E_MEM);

        for (k = 0; k < *o; k++) {
            Theta[j][k] = new CompnentDistribution* [(unsigned INT)C[j][k]];

            R_CHECK(NULL == Theta[j][k], E_MEM);

            for (l = 0; l < C[j][k]; l++) {
                Theta[j][k][l] = new CompnentDistribution(rebmix);

                R_CHECK(NULL == Theta[j][k][l], E_MEM);

                A[0] = A[1] = A[2] = d[k];

                *Error = Theta[j][k][l]->Realloc(d[k], 3, A);

                R_CHECK(*Error != E_OK, *Error);

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
                    if (!strcmp(pdf[i], "Gumbel")) {
                        Theta[j][k][l]->pdf_[m] = pfGumbel;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                        Theta[j][k][l]->Theta_[2][m] = theta3[i];
                    }
                    else
                    if (!strcmp(pdf[i], "vonMises")) {
                        Theta[j][k][l]->pdf_[m] = pfvonMises;

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
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else
                    if (!strcmp(pdf[i], "uniform")) {
                        Theta[j][k][l]->pdf_[m] = pfUniform;

                        Theta[j][k][l]->Theta_[0][m] = theta1[i];
                        Theta[j][k][l]->Theta_[1][m] = theta2[i];
                    }
                    else {
                        R_CHECK(1, E_ARG);
                    }

                    i++;
                }
            }
        }
    }

    dmax = d[0]; for (i = 1; i < *o; i++) if (d[i] > dmax) dmax = d[i];

    Y = (FLOAT**)malloc(dmax * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < dmax; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }

    for (i = 0; i < *n; i++) {
        Z[i] = 1; MaxMixDist = (FLOAT)0.0;

        for (j = 0; j < *s; j++) {
            k = 0; MixDist = (FLOAT)1.0;

            for (l = 0; l < *o; l++) {
                for (m = 0; m < d[l]; m++) {
                    Y[m][0] = X[i + (*n) * (m + k)];
                }

                *Error = rebmix->MixtureDist(0, Y, C[j][l], Q[j][l], Theta[j][l], &CmpDist);

                R_CHECK(*Error != E_OK, *Error);

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

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RCLSMIX

// Returns clustered observations in R.

void RCLRMIX(INT    *n,      // Total number of independent observations.
             double *X,      // Pointer to the input array X.
             INT    *d,      // Number of independent random variables.
             INT    *c,      // Number of components.
             double *W,      // Component weights.
             char   **pdf,   // Component parameters.
             double *theta1, // Component parameters.
             double *theta2, // Component parameters.
             double *theta3, // Component parameters.
             INT    *Z,      // Pointer to the output array Z.
             INT    *Error)  // Error code.
{
    Rebmix               *rebmix = NULL;
    CompnentDistribution **Theta = NULL;
    FLOAT                CmpDist, MaxCmpDist, **Y = NULL;
    INT                  A[3], i, j, k;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    rebmix->length_pdf_ = *d;

    Theta = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == Theta, E_MEM);

    i = 0;

    for (j = 0; j < *c; j++) {
        Theta[j] = new CompnentDistribution(rebmix);

        R_CHECK(NULL == Theta[j], E_MEM);

        A[0] = A[1] = A[2] = *d;

        *Error = Theta[j]->Realloc(*d, 3, A);

        R_CHECK(*Error != E_OK, *Error);

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
            if (!strcmp(pdf[i], "Gumbel")) {
                Theta[j]->pdf_[k] = pfGumbel;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
                Theta[j]->Theta_[2][k] = theta3[i];
            }
            else
            if (!strcmp(pdf[i], "vonMises")) {
                Theta[j]->pdf_[k] = pfvonMises;

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
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else
            if (!strcmp(pdf[i], "uniform")) {
                Theta[j]->pdf_[k] = pfUniform;

                Theta[j]->Theta_[0][k] = theta1[i];
                Theta[j]->Theta_[1][k] = theta2[i];
            }
            else {
                R_CHECK(1, E_ARG);
            }

            i++;
        }
    }

    Y = (FLOAT**)malloc(*d * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < *d; i++) {
        Y[i] = (FLOAT*)malloc(sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }

    for (i = 0; i < *n; i++) {
        for (j = 0; j < *d; j++) {
            Y[j][0] = X[i + (*n) * j];
        }

        Z[i] = 1; MaxCmpDist = (FLOAT)0.0;

        for (j = 0; j < *c; j++) {
            *Error = rebmix->ComponentDist(0, Y, Theta[j], &CmpDist, NULL);

            R_CHECK(*Error != E_OK, *Error);

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

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RCLRMIX

void RPreprocessingKNNMIX(INT    *k,     // k-nearest neighbours.
                          double *h,     // Normalizing vector.
                          INT    *n,     // Total number of independent observations.
                          INT    *d,     // Number of independent random variables.
                          double *x,     // Pointer to the input array x.
                          double *y,     // Pointer to the output array y.
                          INT    *Error) // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    INT    i, j, l;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    rebmix->n_ = rebmix->nr_ = *n;
    rebmix->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 3) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmix->PreprocessingKNN(*k, h, Y);

    R_CHECK(*Error != E_OK, *Error);

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 3; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RPreprocessingKNNMIX

void RPreprocessingKDEMIX(double *h,     // Sides of the hypersquare.
                          INT    *n,     // Total number of independent observations.
                          INT    *d,     // Number of independent random variables.
                          double *x,     // Pointer to the input array x.
                          double *y,     // Pointer to the output array y.
                          INT    *Error) // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    INT    i, j, l;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    rebmix->n_ = rebmix->nr_ = *n;
    rebmix->length_pdf_ = *d;

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 2) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    *Error = rebmix->PreprocessingKDE(h, Y);

    R_CHECK(*Error != E_OK, *Error);

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 2; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RPreprocessingKDEMIX

void RPreprocessingHMIX(double *h,          // Sides of the hypersquare.
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
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    INT    i, j, l;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    rebmix->n_ = rebmix->nr_ = *n;
    rebmix->nc_ = rebmix->length_pdf_ = *d;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmix->Y_, E_MEM);

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmix->Y_[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }

    *Error = rebmix->PreprocessingH(h, y0, ymin, ymax, k, Y);

    R_CHECK(*Error != E_OK, *Error);

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 1; j++) {
        for (l = 0; l < *k; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RPreprocessingHMIX

void RPreprocessingKMIX(double *h,      // Sides of the hypersquare.
                        INT    *d,      // Number of independent random variables.
                        INT    *n,      // Length of x.
                        double *x,      // Pointer to the array x.
                        INT    *Error)  // Error code.
{
    INT i, j, k, l, dn, jn, kn;

    *Error = E_OK;

    R_CHECK(*n < 1, E_ARG);

    dn = (*d) * (*n);

    l = 0;

    for (i = 0; i < *n; i++) {
        for (j = 0; j <= *d; j++) {
            jn = j * (*n); x[l + jn] = x[i + jn];
        }

        for (j = 0; j < l; j++) {
            for (k = 0; k < *d; k++) {
                kn = k * (*n);

                if ((FLOAT)fabs(x[j + kn] - x[l + kn]) > (FLOAT)0.5 * h[k]) goto S0;
            }

            x[j + dn] += x[l + dn]; goto S1;
S0:;
        }

        (l)++;
S1:;
    }

    *n = l;

EEXIT:

    R_RETURN(*Error);
} // RPreprocessingKMIX

void RInformationCriterionKNNMIX(double *h,            // Sides of the hypersquare.
                                 INT    *k,            // k-nearest neighbours.
                                 char   **Criterion,   // Information criterion type.
                                 INT    *c,            // Number of components.
                                 double *W,            // Component weights.
                                 double *Theta1,       // Input parameters theta 1.
                                 double *Theta2,       // Input parameters theta 2.
                                 double *Theta3,       // Input parameters theta 3. 
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
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    INT    i, j, l, m;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

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
        R_CHECK(1, E_ARG);
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->W_, E_MEM);

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    R_CHECK(NULL == rebmix->IniTheta_, E_MEM);

    rebmix->nc_ = rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (INT*)malloc(rebmix->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmix->length_theta_, E_MEM);

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != E_OK, *Error);

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
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
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
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            R_CHECK(1, E_ARG);
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        rebmix->IniTheta_->Theta_[0][j] = Theta1[i];
        rebmix->IniTheta_->Theta_[1][j] = Theta2[i];
        rebmix->IniTheta_->Theta_[2][j] = Theta3[i];

        i++;
    }

    rebmix->MixTheta_ = new CompnentDistribution*[(unsigned INT)(*c)];

    R_CHECK(NULL == rebmix->MixTheta_, E_MEM);

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        R_CHECK(NULL == rebmix->MixTheta_[i], E_MEM);

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        R_CHECK(*Error != E_OK, *Error);
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = rebmix->nr_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmix->Y_, E_MEM);

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmix->Y_[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 3) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);

        if (i < rebmix->length_pdf_) {
            for (j = 0; j < rebmix->nr_; j++) Y[i][j] = rebmix->Y_[i][j];
        }
    }

    *Error = rebmix->PreprocessingKNN(*k, h, Y);

    R_CHECK(*Error != E_OK, *Error);

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

    R_CHECK(*Error != E_OK, *Error);

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 3; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RInformationCriterionKNNMIX

void RInformationCriterionKDEMIX(double *h,            // Sides of the hypersquare.
                                 char   **Criterion,   // Information criterion type.
                                 INT    *c,            // Number of components.
                                 double *W,            // Component weights.
                                 double *Theta1,       // Input parameters theta 1.
                                 double *Theta2,       // Input parameters theta 2.
                                 double *Theta3,       // Input parameters theta 3. 
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
    Rebmix *rebmix = NULL;
    FLOAT  logV, **Y = NULL;
    INT    i, j, l, m;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

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
        R_CHECK(1, E_ARG);
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->W_, E_MEM);

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    R_CHECK(NULL == rebmix->IniTheta_, E_MEM);

    rebmix->nc_ = rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (INT*)malloc(rebmix->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmix->length_theta_, E_MEM);

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != E_OK, *Error);

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
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
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
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            R_CHECK(1, E_ARG);
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        rebmix->IniTheta_->Theta_[0][j] = Theta1[i];
        rebmix->IniTheta_->Theta_[1][j] = Theta2[i];
        rebmix->IniTheta_->Theta_[2][j] = Theta3[i];

        i++;
    }

    rebmix->MixTheta_ = new CompnentDistribution*[(unsigned INT)(*c)];

    R_CHECK(NULL == rebmix->MixTheta_, E_MEM);

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        R_CHECK(NULL == rebmix->MixTheta_[i], E_MEM);

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        R_CHECK(*Error != E_OK, *Error);
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = rebmix->nr_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmix->Y_, E_MEM);

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmix->Y_[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 2) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);

        if (i < rebmix->length_pdf_) {
            for (j = 0; j < rebmix->nr_; j++) Y[i][j] = rebmix->Y_[i][j];
        }
    }

    *Error = rebmix->PreprocessingKDE(h, Y);

    R_CHECK(*Error != E_OK, *Error);

    rebmix->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    *Error = rebmix->InformationCriterionKDE(logV,
                                             Y,
                                             *c,
                                             rebmix->W_,
                                             rebmix->MixTheta_,
                                             IC,
                                             logL,
                                             M,
                                             D);

    R_CHECK(*Error != E_OK, *Error);

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 2; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RInformationCriterionKDEMIX

void RInformationCriterionHMIX(double *h,            // Sides of the hypersquare.
                               double *y0,           // Origins.
                               double *ymin,         // Minimum observations.
                               double *ymax,         // Maximum observations.
                               INT    *k,            // Total number of bins.
                               char   **Criterion,   // Information criterion type.
                               INT    *c,            // Number of components.
                               double *W,            // Component weights.
                               double *Theta1,       // Input parameters theta 1.
                               double *Theta2,       // Input parameters theta 2.
                               double *Theta3,       // Input parameters theta 3. 
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
    Rebmix *rebmix = NULL;
    FLOAT  logV, **Y = NULL;
    INT    i, j, l, m;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

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
        R_CHECK(1, E_ARG);
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->W_, E_MEM);

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    R_CHECK(NULL == rebmix->IniTheta_, E_MEM);

    rebmix->nc_ = rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (INT*)malloc(rebmix->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmix->length_theta_, E_MEM);

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != E_OK, *Error);

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
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
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
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            R_CHECK(1, E_ARG);
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        rebmix->IniTheta_->Theta_[0][j] = Theta1[i];
        rebmix->IniTheta_->Theta_[1][j] = Theta2[i];
        rebmix->IniTheta_->Theta_[2][j] = Theta3[i];

        i++;
    }

    rebmix->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == rebmix->MixTheta_, E_MEM);

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        R_CHECK(NULL == rebmix->MixTheta_[i], E_MEM);

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        R_CHECK(*Error != E_OK, *Error);
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = rebmix->nr_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmix->Y_, E_MEM);

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmix->Y_[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }

    *Error = rebmix->PreprocessingH(h, y0, ymin, ymax, k, Y);

    R_CHECK(*Error != E_OK, *Error);

    rebmix->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    *Error = rebmix->InformationCriterionH(logV,
                                           *k,
                                           Y,
                                           *c,
                                           rebmix->W_,
                                           rebmix->MixTheta_,
                                           IC,
                                           logL,
                                           M,
                                           D);

    R_CHECK(*Error != E_OK, *Error);

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RInformationCriterionHMIX

void RInformationCriterionKMIX(double *h,            // Sides of the hypersquare.
                               char   **Criterion,   // Information criterion type.
                               INT    *c,            // Number of components.
                               double *W,            // Component weights.
                               double *Theta1,       // Input parameters theta 1.
                               double *Theta2,       // Input parameters theta 2.
                               double *Theta3,       // Input parameters theta 3. 
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
    Rebmix *rebmix = NULL;
    FLOAT  logV, **Y = NULL;
    INT    i, j, l, m;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

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
        R_CHECK(1, E_ARG);
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->W_, E_MEM);

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    R_CHECK(NULL == rebmix->IniTheta_, E_MEM);

    rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (INT*)malloc(rebmix->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmix->length_theta_, E_MEM);

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != E_OK, *Error);

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
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
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
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            R_CHECK(1, E_ARG);
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        rebmix->IniTheta_->Theta_[0][j] = Theta1[i];
        rebmix->IniTheta_->Theta_[1][j] = Theta2[i];
        rebmix->IniTheta_->Theta_[2][j] = Theta3[i];

        i++;
    }

    rebmix->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == rebmix->MixTheta_, E_MEM);

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        R_CHECK(NULL == rebmix->MixTheta_[i], E_MEM);

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        R_CHECK(*Error != E_OK, *Error);
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = rebmix->nr_ = *n;

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 1; j++) {
        Y[j] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[j], E_MEM);

        for (l = 0; l < rebmix->nr_; l++) {
            Y[j][l] = x[i]; i++;
        }
    }

    rebmix->n_ = 0;

    for (l = 0; l < rebmix->nr_; l++) {
        rebmix->n_ += (INT)Y[rebmix->length_pdf_][l];
    }

    rebmix->cmax_ = *c;

    logV = (FLOAT)0.0;

    for (i = 0; i < rebmix->length_pdf_; i++) {
        logV += (FLOAT)log(h[i]);
    }

    rebmix->Y_type_ = 1;

    *Error = rebmix->InformationCriterionH(logV,
                                           *n,
                                           Y,
                                           *c,
                                           rebmix->W_,
                                           rebmix->MixTheta_,
                                           IC,
                                           logL,
                                           M,
                                           D);

    R_CHECK(*Error != E_OK, *Error);

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RInformationCriterionKMIX

void RInformationCriterionMIX(char   **Criterion,   // Information criterion type.
                              INT    *c,            // Number of components.
                              double *W,            // Component weights.
                              double *Theta1,       // Input parameters theta 1.
                              double *Theta2,       // Input parameters theta 2.
                              double *Theta3,       // Input parameters theta 3. 
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
    Rebmix *rebmix = NULL;
    INT    i, j, l, m;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

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
        R_CHECK(1, E_ARG);
    }

    rebmix->W_ = (FLOAT*)malloc(*c * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->W_, E_MEM);

    for (i = 0; i < *c; i++) rebmix->W_[i] = W[i];

    rebmix->IniTheta_ = new CompnentDistribution(rebmix);

    R_CHECK(NULL == rebmix->IniTheta_, E_MEM);

    rebmix->nc_ = rebmix->length_pdf_ = *length_pdf;

    rebmix->length_Theta_ = *length_Theta;

    rebmix->length_theta_ = (INT*)malloc(rebmix->length_Theta_ * sizeof(INT));

    R_CHECK(NULL == rebmix->length_theta_, E_MEM);

    *Error = rebmix->IniTheta_->Realloc(*length_pdf, *length_Theta, length_theta);

    R_CHECK(*Error != E_OK, *Error);

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
        if (!strcmp(pdf[i], "Gumbel")) {
            rebmix->IniTheta_->pdf_[i] = pfGumbel;
        }
        else
        if (!strcmp(pdf[i], "vonMises")) {
            rebmix->IniTheta_->pdf_[i] = pfvonMises;
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
        else
        if (!strcmp(pdf[i], "uniform")) {
            rebmix->IniTheta_->pdf_[i] = pfUniform;
        }
        else {
            R_CHECK(1, E_ARG);
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        rebmix->IniTheta_->Theta_[0][j] = Theta1[i];
        rebmix->IniTheta_->Theta_[1][j] = Theta2[i];
        rebmix->IniTheta_->Theta_[2][j] = Theta3[i];

        i++;
    }

    rebmix->MixTheta_ = new CompnentDistribution* [(unsigned INT)(*c)];

    R_CHECK(NULL == rebmix->MixTheta_, E_MEM);

    for (i = 0; i < *c; i++) {
        rebmix->MixTheta_[i] = new CompnentDistribution(rebmix);

        R_CHECK(NULL == rebmix->MixTheta_[i], E_MEM);

        *Error = rebmix->MixTheta_[i]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        R_CHECK(*Error != E_OK, *Error);
    }

    for (i = 0; i < *c; i++) {
        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[i]->pdf_[j] = rebmix->IniTheta_->pdf_[j];
        }
    }

    i = 0;

    for (j = 0; j < rebmix->length_Theta_; j++) if (rebmix->IniTheta_->Theta_[j]) {
        for (l = 0; l < *c; l++) {
            for (m = 0; m < rebmix->length_theta_[j]; m++) {
                rebmix->MixTheta_[l]->Theta_[j][m] = Theta[i];

                i++;
            }
        }
    }

    rebmix->n_ = rebmix->nr_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmix->Y_, E_MEM);

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmix->Y_[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    rebmix->cmax_ = *c;

    *Error = rebmix->InformationCriterion(*c,
                                          rebmix->W_,
                                          rebmix->MixTheta_,
                                          IC,
                                          logL,
                                          M,
                                          D);

    R_CHECK(*Error != E_OK, *Error);

EEXIT:  
    
    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RInformationCriterionMIX

void RCombineComponentsMIX(INT    *c,            // Number of components.
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
    Rebmix *rebmix = NULL;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    *Error =  rebmix->Set(NULL,         // Preprocessing type.
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

    R_CHECK(*Error != E_OK, *Error);

/// Panic Branislav
    if (!strcmp(Rule[0], "Entropy")) {
        *Error = rebmix->CombineComponentsEntropy(rebmix->cmax_,
            rebmix->W_,
            rebmix->MixTheta_,
            tau,
            F,
            T,
            EN,
            ED,
            PSS);

        R_CHECK(*Error != E_OK, *Error);
    }
    else
    if (!strcmp(Rule[0], "Demp")) {
        *Error = rebmix->CombineComponentsDemp(rebmix->cmax_,
            rebmix->W_,
            rebmix->MixTheta_,
            tau,
            F,
            T,
            EN,
            ED,
            PSS);

        R_CHECK(*Error != E_OK, *Error);
    }
    else {
        R_CHECK(1, E_ARG);
    }
/// End

EEXIT:  
    
    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // RCombineComponentsMIX

void RvonMisesPdf(INT *n, double *y, double *Mean, double *Kappa, double *f)
{
    FLOAT A;
    INT   i;

    A = Pi2 * BesselI0(*Kappa);

    for (i = 0; i < *n; i++) {
        if (y[i] > Pi2) {
            f[i] = (FLOAT)0.0;
        }
        else
        if (y[i] < (FLOAT)0.0) {
            f[i] = (FLOAT)0.0;
        }
        else {
            f[i] = (FLOAT)exp(*Kappa * (FLOAT)cos(y[i] - *Mean)) / A;
        }
    }
} // RvonMisesPdf

void RvonMisesCdf(INT *n, double *y, double *Mean, double *Kappa, double *F, INT *Error)
{
    FLOAT A[3], Io, In, I0, I1;
    INT   i, j;

    *Error = E_OK;

    I0 = BesselI0(*Kappa); I1 = BesselI1(*Kappa);

    A[0] = (FLOAT)1.0 / Pi2; A[1] = (FLOAT)2.0 * A[0] / I0;

    for (i = 0; i < *n; i++) {
        if (y[i] > Pi2) {
            F[i] = (FLOAT)1.0;
        }
        else
        if (y[i] < (FLOAT)0.0) {
            F[i] = (FLOAT)0.0;
        }
        else {
            j = 1; *Error = E_CON; Io = I0; In = I1; F[i] = A[0] * y[i];

            while ((j <= ItMax) && (*Error != E_OK)) {
                F[i] += A[1] * In * ((FLOAT)sin(j * (y[i] - *Mean)) + (FLOAT)sin(j * *Mean)) / j;

                if (In < Eps) *Error = E_OK;

                A[2] = Io - (FLOAT)2.0 * j * In / *Kappa; Io = In; In = A[2];

                j++;
            }

            if (F[i] > (FLOAT)1.0) {
                F[i] = (FLOAT)1.0;
            }
            else
            if (F[i] < (FLOAT)0.0) {
                F[i] = (FLOAT)0.0;
            }
        }
    }

    R_RETURN(*Error);
} // RvonMisesCdf

void RGumbelPdf(INT *n, double *y, double *Mean, double *Sigma, double *Xi, double *f)
{
    FLOAT A;
    INT   i;

    for (i = 0; i < *n; i++) {
        A = *Xi * (y[i] - *Mean) / (*Sigma);

        f[i] = (FLOAT)exp(A - (FLOAT)exp(A)) / (*Sigma);
    }
} // RGumbelPdf

void RGumbelCdf(INT *n, double *y, double *Mean, double *Sigma, double *Xi, double *F)
{
    FLOAT A;
    INT   i;

    for (i = 0; i < *n; i++) {
        if (*Xi > Eps) {
            A = (y[i] - *Mean) / (*Sigma);

            F[i] = (FLOAT)1.0 - (FLOAT)exp(-(FLOAT)exp(A));
        }
        else {
            A = -(y[i] - *Mean) / (*Sigma);

            F[i] = (FLOAT)exp(-(FLOAT)exp(A));
        }
    }
} // RGumbelCdf

/// Panic Branislav & Marko Nagode.
// Optimal number of bins calculation.

void Roptbins(INT    *d,           // Number of independent random variables.
              INT    *n,           // Number of observations.
              double *x,           // Dataset.
              char   **Rule,       // Rule.
              INT    *length_ymin, // Length of ymin.
              double *ymin,        // Minimum observations.
              INT    *length_ymax, // Length of ymax.
              double *ymax,        // Maximum observations.
              INT    *kmin,        // Minimum number of bins.
              INT    *kmax,        // Maximum number of bins.
              INT    *opt_k,       // Optimal number of bins.
              INT    *Error)       // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    INT    i, j, l, *opt_k_tmp = NULL;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    rebmix->nc_ = rebmix->length_pdf_ = *d;

    rebmix->n_ = rebmix->nr_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmix->Y_, E_MEM);

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmix->Y_[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    rebmix->ymin_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->ymin_, E_MEM);

    if (*length_ymin > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = ymin[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->nr_; j++) {
                if (rebmix->Y_[i][j] < rebmix->ymin_[i]) rebmix->ymin_[i] = rebmix->Y_[i][j];
            }
        }
    }

    rebmix->ymax_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->ymax_, E_MEM);

    if (*length_ymax > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = ymax[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->nr_; j++) {
                if (rebmix->Y_[i][j] > rebmix->ymax_[i]) rebmix->ymax_[i] = rebmix->Y_[i][j];
            }
        }
    }

    rebmix->y0_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->y0_, E_MEM);

    rebmix->h_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->h_, E_MEM);

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }

    *Error = rebmix->Initialize();

    R_CHECK(*Error != E_OK, *Error);

    if (!strcmp(Rule[0], "Sturges")) {
        opt_k[0] = (INT)ceil((FLOAT)1.0 + (FLOAT)log2((FLOAT)rebmix->n_));

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "Log10")) {
        opt_k[0] = (INT)ceil((FLOAT)10.0 * (FLOAT)log10((FLOAT)rebmix->n_));

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "RootN")) {
        opt_k[0] = (INT)ceil((FLOAT)2.0 * (FLOAT)sqrt((FLOAT)rebmix->n_));

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "Knuth equal")) {
        FLOAT logp, logpopt, M;
        INT   k;

        logpopt = -FLOAT_MAX;

        for (l = *kmin; l < *kmax + 1; l++) {
            for (j = 0; j < rebmix->length_pdf_; j++) {
                rebmix->h_[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / l;

                rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * rebmix->h_[j];
            }

            *Error = rebmix->PreprocessingH(rebmix->h_, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, &k, Y);

            R_CHECK(*Error != E_OK, *Error);

            if (k > rebmix->kmax_) {
                break;
            }

            M = (FLOAT)1.0;

            for (j = 0; j < rebmix->length_pdf_; j++) {
                M *= (FLOAT)l;
            }

            logp = (FLOAT)rebmix->n_ * (FLOAT)log(M) + (FLOAT)Gammaln((FLOAT)0.5 * M) - M * (FLOAT)Gammaln((FLOAT)0.5) - (FLOAT)Gammaln((FLOAT)rebmix->n_ + (FLOAT)0.5 * M);

            for (j = 0; j < k; j++) {
                logp += Gammaln(Y[rebmix->length_pdf_][j] + (FLOAT)0.5);
            }

            logp += Max(M - (FLOAT)k, (FLOAT)0.0) * Gammaln((FLOAT)0.5);

            if (logp > logpopt) {
                logpopt = logp; opt_k[0] = l;
            }
        }

        for (j = 1; j < rebmix->length_pdf_; j++) {
            opt_k[j] = opt_k[0];
        }
    }
    else
    if (!strcmp(Rule[0], "Knuth unequal")) {
        FLOAT logp, logpopt, M, M_broken;
        INT   k, num_e, lim_broken, div, inc, Converged = 1;
        INT   num_r, MAX_ITER_NUM = 131072;

        logpopt = -FLOAT_MAX;

        opt_k_tmp = (INT*)malloc(rebmix->length_pdf_ * sizeof(INT));

        R_CHECK(NULL == opt_k_tmp, E_MEM);

        for (j = 0; j < rebmix->length_pdf_; j++) {
            opt_k_tmp[j] = 1;
        }

        while (Converged) {
            for (j = 0; j < rebmix->length_pdf_; j++) {
                for (i = 0; i < rebmix->length_pdf_; i++) {
                    if (i != j) {
                        rebmix->h_[i] = (rebmix->ymax_[i] - rebmix->ymin_[i]) / opt_k_tmp[i];

                        rebmix->y0_[i] = rebmix->ymin_[i] + (FLOAT)0.5 * rebmix->h_[i];
                    }
                }

                opt_k[j] = opt_k_tmp[j];

                for (l = *kmin; l < *kmax + 1; l++) {
                    rebmix->h_[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / l;

                    rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * rebmix->h_[j];

                    *Error = rebmix->PreprocessingH(rebmix->h_, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, &k, Y);

                    R_CHECK(*Error != E_OK, *Error);

                    if (k > rebmix->kmax_ && opt_k_tmp[j] != 1) {
                        break;
                    }

                    M = (FLOAT)1.0;

                    for (i = 0; i < rebmix->length_pdf_; i++) {
                        if (i != j) {
                            M *= (FLOAT)opt_k_tmp[i];
                        }
                        else{
                            M *= (FLOAT)l;
                        }
                    }

                    logp = (FLOAT)rebmix->n_ * (FLOAT)log(M) + (FLOAT)Gammaln((FLOAT)0.5 * M) - M * (FLOAT)Gammaln((FLOAT)0.5) - (FLOAT)Gammaln((FLOAT)rebmix->n_ + (FLOAT)0.5 * M);

                    for (i = 0; i < k; i++) {
                        logp += Gammaln(Y[rebmix->length_pdf_][i] + (FLOAT)0.5);
                    }

                    logp += Max(M - (FLOAT)k, (FLOAT)0.0) * Gammaln((FLOAT)0.5);

                    if (logp > logpopt || opt_k_tmp[j] == 1) {
                        logpopt = logp; opt_k_tmp[j] = l;
                    }
                }

                if (j == 0 && opt_k_tmp[j] == opt_k[j]) {
                    for (i = 0; i < rebmix->length_pdf_; i++) opt_k[i] = opt_k_tmp[i];

                    Converged = 0; break;
                }
            }
        }

        *kmin = *kmax = opt_k[0];

        for (j = 1; j < rebmix->length_pdf_; j++) {
            if (opt_k[j] < *kmin) {
                *kmin = opt_k[j];
            }
            else
            if (opt_k[j] > *kmax) {
                *kmax = opt_k[j];
            }
        }

        num_e = *kmax - *kmin + 1;

        num_r = (INT)floor(pow((FLOAT)num_e, rebmix->length_pdf_));

        if (num_r > MAX_ITER_NUM || num_r <= 0) {
            for (j = *kmax - 1; j >= *kmin; j--) {
                *kmax = j;

                num_e = *kmax - *kmin + 1;

                num_r = (INT)floor(pow((FLOAT)num_e, rebmix->length_pdf_));

                if (num_r <= MAX_ITER_NUM && num_r > 0) {
                    break;
                }
            }
        }

        lim_broken = 0; M_broken = (FLOAT)0.0;

        R_CHECK(*kmax == *kmin, E_ARG);

        for (i = 0; i < num_r; i++) {
            for (j = rebmix->length_pdf_ - 1; j >= 0; j--) {
                div = (INT)floor(pow((FLOAT)num_e, j));

                inc = (i / div) % num_e;

                opt_k_tmp[j] = *kmin + inc;
            }

            for (j = 0; j < rebmix->length_pdf_; j++) {
                rebmix->h_[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / opt_k_tmp[j];

                rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * rebmix->h_[j];
            }

            M = (FLOAT)1.0;

            for (j = 0; j < rebmix->length_pdf_; j++) {
                M *= (FLOAT)opt_k_tmp[j];
            }

            if (lim_broken) {
                if (M > M_broken) continue;
            }

            *Error = rebmix->PreprocessingH(rebmix->h_, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, &k, Y);

            R_CHECK(*Error != E_OK, *Error);

            if (k > rebmix->kmax_) {
                lim_broken = 1; M_broken = M; continue;
            }

            logp = (FLOAT)rebmix->n_ * (FLOAT)log(M) + (FLOAT)Gammaln((FLOAT)0.5 * M) - M * (FLOAT)Gammaln((FLOAT)0.5) - (FLOAT)Gammaln((FLOAT)rebmix->n_ + (FLOAT)0.5 * M);

            for (j = 0; j < k; j++) {
                logp += Gammaln(Y[rebmix->length_pdf_][j] + (FLOAT)0.5);
            }

            logp += Max(M - (FLOAT)k, (FLOAT)0.0) * Gammaln((FLOAT)0.5);

            if (logp > logpopt) {
                logpopt = logp;

                for (j = 0; j < rebmix->length_pdf_; j++) opt_k[j] = opt_k_tmp[j];
            }
        }
    }

EEXIT:  
    
    if (opt_k_tmp) free(opt_k_tmp);
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // Roptbins

/// Binning of data.

void Rbins(INT    *d,           // Number of independent random variables.
           INT    *n,           // Number of observations.
           double *x,           // Dataset.
           INT    *length_ymin, // Length of ymin.
           double *ymin,        // Minimum observations.
           INT    *length_ymax, // Length of ymax.
           double *ymax,        // Maximum observations.
           INT    *k,           // Number of bins.
           INT    *length_y,    // Length of y.
           double *y,           // Binned dataset.
           INT    *Error)       // Error code.
{
    Rebmix *rebmix = NULL;
    FLOAT  **Y = NULL;
    INT    i, j, l;

    *Error = E_OK;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    rebmix->nc_ = rebmix->length_pdf_ = *d;

    rebmix->n_ = rebmix->nr_ = *n;

    rebmix->Y_ = (FLOAT**)malloc(rebmix->nc_ * sizeof(FLOAT*));

    R_CHECK(NULL == rebmix->Y_, E_MEM);

    for (i = 0; i < rebmix->length_pdf_; i++) {
        rebmix->Y_[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == rebmix->Y_[i], E_MEM);
    }

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        for (l = 0; l < rebmix->nr_; l++) {
            rebmix->Y_[j][l] = x[i]; i++;
        }
    }

    rebmix->ymin_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->ymin_, E_MEM);

    if (*length_ymin > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = ymin[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymin_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->nr_; j++) {
                if (rebmix->Y_[i][j] < rebmix->ymin_[i]) rebmix->ymin_[i] = rebmix->Y_[i][j];
            }
        }
    }

    rebmix->ymax_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->ymax_, E_MEM);

    if (*length_ymax > 0) {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = ymax[i];
        }
    }
    else {
        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->ymax_[i] = rebmix->Y_[i][0];

            for (j = 1; j < rebmix->nr_; j++) {
                if (rebmix->Y_[i][j] > rebmix->ymax_[i]) rebmix->ymax_[i] = rebmix->Y_[i][j];
            }
        }
    }

    rebmix->y0_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->y0_, E_MEM);

    rebmix->h_ = (FLOAT*)malloc(rebmix->length_pdf_ * sizeof(FLOAT));

    R_CHECK(NULL == rebmix->h_, E_MEM);

    Y = (FLOAT**)malloc((rebmix->length_pdf_ + 1) * sizeof(FLOAT*));

    R_CHECK(NULL == Y, E_MEM);

    for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
        Y[i] = (FLOAT*)malloc(rebmix->nr_ * sizeof(FLOAT));

        R_CHECK(NULL == Y[i], E_MEM);
    }


    for (j = 0; j < rebmix->length_pdf_; j++) {
        rebmix->h_[j] = (rebmix->ymax_[j] - rebmix->ymin_[j]) / k[j];

        rebmix->y0_[j] = rebmix->ymin_[j] + (FLOAT)0.5 * rebmix->h_[j];
    }

    *Error = rebmix->PreprocessingH(rebmix->h_, rebmix->y0_, rebmix->ymin_, rebmix->ymax_, length_y, Y);

    R_CHECK(*Error != E_OK, *Error);

    i = 0;

    for (j = 0; j < rebmix->length_pdf_ + 1; j++) {
        for (l = 0; l < *length_y; l++) {
            y[i] = Y[j][l]; i++;
        }
    }

EEXIT:  
    
    if (Y) {
        for (i = 0; i < rebmix->length_pdf_ + 1; i++) {
            if (Y[i]) free(Y[i]);
        }

        free(Y);
    }

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
} // Rbins

/// Panic Branislav
void REMMIX(INT    *d,                 // Number of independent random variables.
            INT    *n,                 // Number of observations.
            double *Y,                 // Input dataset.
            INT    *Y_type,            // Dataset type.
            char   **pdf,              // Parametric family types.
            INT    *c,                 // Input number of components.
            double *W,                 // Input weights of components.
            double *Theta1,            // Input parameters theta 1.
            double *Theta2,            // Input parameters theta 2.
            double *Theta3,            // Input parameters theta 3. 
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
            INT    *Error)             // Error number.
{
    Rebmix *rebmix = NULL;
    INT    A[3], i = 0, j = 0, l = 0, length_Theta = 3;

    *Error = E_OK;

    A[0] = A[1] = A[2] = *d;

    rebmix = new Rebmix;

    R_CHECK(NULL == rebmix, E_MEM);

    *Error = rebmix->Set(NULL,              // Preprocessing type.
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

    rebmix->EM_strategy_ = strategy_single;

    i = 0;

    for (j = 0; j < rebmix->length_pdf_; j++) {
        rebmix->IniTheta_->Theta_[0][j] = Theta1[i];
        rebmix->IniTheta_->Theta_[1][j] = Theta2[i];
        rebmix->IniTheta_->Theta_[2][j] = Theta3[i];

        i++;
    }

    rebmix->MixTheta_ = new CompnentDistribution*[(unsigned INT)(rebmix->cmax_)];

    R_CHECK(NULL == rebmix->MixTheta_, E_MEM);

    for (l = 0; l < rebmix->cmax_; l++) {
        rebmix->MixTheta_[l] = new CompnentDistribution(rebmix);

        R_CHECK(NULL == rebmix->MixTheta_[l], E_MEM);

        *Error = rebmix->MixTheta_[l]->Realloc(rebmix->length_pdf_, rebmix->length_Theta_, rebmix->length_theta_);

        R_CHECK(*Error != E_OK, *Error);

        for (i = 0; i < rebmix->length_pdf_; i++) {
            rebmix->MixTheta_[l]->pdf_[i] = rebmix->IniTheta_->pdf_[i];
        }
    }

    i = 0;

    rebmix->W_ = (double*)malloc(rebmix->cmax_ * sizeof(double));

    R_CHECK(NULL == rebmix->W_, E_MEM);

    for (l = 0; l < *c; l++) {
        rebmix->W_[l] = W[l];

        for (j = 0; j < rebmix->length_pdf_; j++) {
            rebmix->MixTheta_[l]->Theta_[0][j] = Theta1[i];
            rebmix->MixTheta_[l]->Theta_[1][j] = Theta2[i];
            rebmix->MixTheta_[l]->Theta_[2][j] = Theta3[i]; 
            
            i++;
        }
    }

    *Error = rebmix->EMInitialize();

    R_CHECK(*Error != E_OK, *Error);

    if (*EMMerge) {
        rebmix->EM_->merge_ = merge_naive;
    }

    *Error = rebmix->EMRun(c, rebmix->W_, rebmix->MixTheta_);

    R_CHECK(*Error != E_OK, *Error);

    *Error = rebmix->EM_->LogLikelihood(*c, rebmix->W_, rebmix->MixTheta_, summary_logL);

    R_CHECK(*Error != E_OK, *Error);

    *Error = rebmix->DegreesOffreedom(*c, rebmix->MixTheta_, summary_M);

    R_CHECK(*Error != E_OK, *Error);

    rebmix->summary_.c = *c;

    *Error = rebmix->Get(n_iter,   // Number of iterations for optimal case.
                         NULL,     // Number of iterations in whole run.
                         NULL,     // Optimal v or optimal k.
                         NULL,     // Optimal class widths of length d.
                         NULL,     // Optimal origins of length d.
                         NULL,     // Optimal minimum observations of length d.
                         NULL,     // Optimal maximum observations of length d.
                         NULL,     // Optimal information criterion.
                         NULL,     // Log-likelihood.
                         NULL,     // Degrees of freedom.
                         NULL,     // Optimal number of components.
                         W,        // Component weights.
                         Theta1,   // Component parameters.
                         Theta2,   // Component parameters.
                         Theta3,   // Component parameters.
                         NULL,     // Length of opt_c, opt_IC, opt_logL, opt_Dmin and opt_D.
                         NULL,     // Numbers of components for optimal v or for optimal k.
                         NULL,     // Information criteria for optimal v or for optimal k.
                         NULL,     // Log-likelihoods for optimal v or for optimal k.
                         NULL,     // Dmin for optimal v or for optimal k.
                         NULL,     // Totals of positive relative deviations for optimal v or for optimal k.
                         NULL,     // Length of all_K and all_IC.
                         NULL,     // All processed numbers of bins v or all processed numbers of nearest neighbours k.
                         NULL);    // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.

    R_CHECK(*Error != E_OK, *Error);

EEXIT:

    if (rebmix) delete rebmix;

    R_RETURN(*Error);
}
/// End

// Fast Histogram calculation.

void Rfhistogram(INT    *K,      // Numbers of bins.
                 double *y0,     // Origins.
                 double *h,      // Sides of the hypersquare.
                 INT    *d,      // Number of independent random variables.
                 INT    *nx,     // Length of x.
                 double *x,      // Pointer to the input array x.
                 INT    *ny,     // Length of y.
                 double *y,      // Pointer to the output array y.
                 INT    *shrink, // If 1 the output array is shrinked.
                 INT    *Error)  // Error code.
{
    INT i, j, k, *l = NULL, *m = NULL, dny;

    *Error = E_OK;

    R_CHECK(*nx < 1, E_ARG);

    l = (INT*)malloc(*d * sizeof(INT));

    R_CHECK(NULL == l, E_MEM);

    m = (INT*)malloc(*d * sizeof(INT));

    R_CHECK(NULL == m, E_MEM);

    dny = (*d) * (*ny); m[*d - 1] = 1;

    for (i = *d - 1; i > 0; i--) {
        m[i - 1] = K[i] * m[i];
    }

    for (i = 0; i < *nx; i++) {
        j = 0;

        for (k = 0; k < *d; k++) {
            l[k] = (INT)floor((x[i + k * (*nx)] - y0[k]) / h[k] + (FLOAT)0.5);

            if (l[k] < 0) l[k] = 0; else if (l[k] >= K[k]) l[k] = K[k] - 1;

            j += l[k] * m[k];
        }

        for (k = 0; k < *d; k++) {
            y[j + k * (*ny)] = y0[k] + l[k] * h[k];
        }

        y[j + dny] += (FLOAT)1.0;
    }

    if (*shrink) {
        i = 0;

        for (j = 0; j < *ny; j++) if (y[j + dny] > FLOAT_MIN) {
            if (i != j) for (k = 0; k <= *d; k++) y[i + k * (*ny)] = y[j + k * (*ny)];

            i++;
        }

        *ny = i;
    }

EEXIT:  
    
    if (m) free(m);

    if (l) free(l);

    R_RETURN(*Error);
} // Rfhistogram

// Compact histogram calculation.

void Rchistogram(INT    *K,      // Numbers of bins.
                 INT    *v,      // Number of nonempty bins.
                 double *y0,     // Origins.
                 double *h,      // Sides of the hypersquare.
                 INT    *d,      // Number of independent random variables.
                 INT    *nx,     // Length of x.
                 double *x,      // Pointer to the input array x.
                 INT    *ny,     // Length of y.
                 double *y,      // Pointer to the output array y.
                 INT    *Error)  // Error code.
{
    INT i, j, k, dny, kny;

    *Error = E_OK;

    R_CHECK(*nx < 1, E_ARG);

    dny = (*d) * (*ny);

    for (i = 0; i < *nx; i++) {
        for (j = 0; j < *d; j++) {
            k = (INT)floor((x[i + j * (*nx)] - y0[j]) / h[j] + (FLOAT)0.5);

            if (k < 0) k = 0; else if (k >= K[j]) k = K[j] - 1;

            y[*v + j * (*ny)] = y0[j] + k * h[j];
        }

        for (j = 0; j < *v; j++) {
            for (k = 0; k < *d; k++) {
                kny = k * (*ny);

                if ((FLOAT)fabs(y[j + kny] - y[*v + kny]) > (FLOAT)0.5 * h[k]) goto S0;
            }

            y[j + dny] += (FLOAT)1.0; goto S1;
S0:;
        }

        y[*v + dny] = (FLOAT)1.0; (*v)++;
S1:;
    }

EEXIT:

    R_RETURN(*Error);
} // Rchistogram

}
