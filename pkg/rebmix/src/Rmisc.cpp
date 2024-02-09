#include "base.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

extern "C" {

// Returns numbers of pixels, label moments for clustered image and adjacency matrix in R.

void RLabelMomentsXY(INT    *nx,     // Image width.
                     INT    *ny,     // Image height.
                     double *Zp,     // Pointer to the predictive cluster membership.
                     INT    *c,      // Number of clusters.
                     double *N,      // Numbers of pixels.
                     double *Mx,     // Raw moments Mx.
                     double *My,     // Raw moments My.  
                     double *Mxy,    // Raw moments Mxy.
                     double *A,      // Adjacency matrix.
                     double *Sigma,  // Scaling factor.
                     INT    *Error)  // Error code.
{
    FLOAT **Mij = NULL, *Mean = NULL, Mul, *Stdev = NULL, Tmp;
    INT   i, ii, j, k, kk, l, n = 0;

    *Error = *c < 2; if (*Error) goto E0;

    Mij = (FLOAT**)malloc(4 * sizeof(FLOAT*));

    *Error = NULL == Mij; if (*Error) goto E0;

    for (i = 0; i < 4; i++) {
        Mij[i] = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

        *Error = NULL == Mij[i]; if (*Error) goto E0;
    }

    Mean = (FLOAT*)calloc((size_t)4, sizeof(FLOAT));

    *Error = NULL == Mean; if (*Error) goto E0;

    Stdev = (FLOAT*)calloc((size_t)4, sizeof(FLOAT));

    *Error = NULL == Stdev; if (*Error) goto E0;

    // Numbers of pixels and raw clustered image moments calculation.   

    for (i = 0; i < *nx; i++) {
        j = (*ny) * i;

        for (k = 0; k < *ny; k++) {
            l = (INT)Zp[j + k] - 1;

            if ((l >= 0) && (l < *c)) {
                Mij[0][l] += (FLOAT)1.0; ii = i + 1; kk = k + 1;

                Mij[1][l] += (FLOAT)ii; Mij[2][l] += (FLOAT)kk;

                Mij[3][l] += (FLOAT)(ii * kk);
            }
        }
    }

    // Number of nonempty clusters and raw moments calculation.
    
    for (i = 0; i < *c; i++) {
        N[i] = Mij[0][i]; 
        
        if (N[i] > FLOAT_MIN) {
            Mij[1][i] /= N[i]; Mij[2][i] /= N[i]; n++;
        }
    }

    //  Means calculation.

    for (i = 1; i < 4; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mean[i] += Mij[i][j];
        }

        Mean[i] /= n;
    }

    //  Standard deviations calculation.

    for (i = 1; i < 4; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Tmp = Mij[i][j] - Mean[i]; Stdev[i] += Tmp * Tmp;
        }
        
        Stdev[i] = (FLOAT)sqrt(Stdev[i] / (n - (FLOAT)1.0)); if (Stdev[i] < Eps) Stdev[i] = Eps;
    }

    for (i = 0; i < *c; i++) {
        Mx[i] = Mij[1][i]; My[i] = Mij[2][i]; Mxy[i] = Mij[3][i] / N[i] - Mx[i] * My[i];
    }

    //  Z-scores calculation.

    for (i = 1; i < 4; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mij[i][j] = (Mij[i][j] - Mean[i]) / Stdev[i];
        }
    }

    // Adjacency matrix calculation based on Gaussian kernel function.

    Mul = (FLOAT)0.5 / (*Sigma) / (*Sigma);

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) {
        j = (*c) * i;

        for (k = i + 1; k < *c; k++) if (N[k] > FLOAT_MIN) {
            ii = j + k; kk = i + (*c) * k; A[ii] = (FLOAT)0.0;
            
            for (l = 1; l < 4; l++) {
                Tmp = Mij[l][i] - Mij[l][k]; A[ii] += Tmp * Tmp;
            }

            A[ii] = A[kk] = (FLOAT)exp(-Mul * A[ii]);
        }
    }

E0: if (Stdev) free(Stdev);

    if (Mean) free(Mean);

    if (Mij) {
        for (i = 0; i < 4; i++) {
            if (Mij[i]) free(Mij[i]);
        }

        free(Mij);
    }
} // RLabelMomentsXY

// Returns numbers of voxels, label moments for clustered image and adjacency matrix in R.

void RLabelMomentsXYZ(INT    *nx,     // Image width.
                      INT    *ny,     // Image height.
                      INT    *nz,     // Image depth.
                      double *Zp,     // Pointer to the predictive cluster membership.
                      INT    *c,      // Number of clusters.
                      double *N,      // Numbers of voxels.
                      double *Mx,     // Raw moments Mx.
                      double *My,     // Raw moments My.
                      double *Mz,     // Raw moments Mz.
                      double *Mxyz,   // Raw moments Mxyz.  
                      double *A,      // Adjacency matrix.
                      double *Sigma,  // Scaling factor.
                      INT    *Error)  // Error code.
{
    FLOAT **Mijk = NULL, *Mean = NULL, Mul, *Stdev = NULL, Tmp;
    INT i, ii, j, k, kk, l, m, mm, n = 0;

    *Error = *c < 2; if (*Error) goto E0;

    Mijk = (FLOAT**)malloc(5 * sizeof(FLOAT*));

    *Error = NULL == Mijk; if (*Error) goto E0;

    for (i = 0; i < 5; i++) {
        Mijk[i] = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

        *Error = NULL == Mijk[i]; if (*Error) goto E0;
    }

    Mean = (FLOAT*)calloc((size_t)5, sizeof(FLOAT));

    *Error = NULL == Mean; if (*Error) goto E0;

    Stdev = (FLOAT*)calloc((size_t)5, sizeof(FLOAT));

    *Error = NULL == Stdev; if (*Error) goto E0;

    //  Numbers of voxels and raw clustered image moments calculation.   

    for (i = 0; i < *nz; i++) {
        j = (*nx) * (*ny) * i;

        for (k = 0; k < *nx; k++) {
            l = j + (*ny) * k;

            for (m = 0; m < *ny; m++) {
                n = (INT)Zp[l + m] - 1;

                if ((n >= 0) && (n < *c)) {
                    Mijk[0][n] += (FLOAT)1.0; ii = i + 1; kk = k + 1; mm = m + 1;

                    Mijk[1][n] += (FLOAT)kk; Mijk[2][n] += (FLOAT)mm; Mijk[3][n] += (FLOAT)ii;

                    Mijk[4][n] += (FLOAT)(ii * kk * mm);
                }
            }
        }
    }

    // Number of nonempty clusters and raw moments calculation.

    for (i = 0; i < *c; i++) {
        N[i] = Mijk[0][i]; 
        
        if (N[i] > FLOAT_MIN) {
            Mijk[1][i] /= N[i]; Mijk[2][i] /= N[i]; Mijk[3][i] /= N[i]; n++;
        }
    }

    //  Means calculation.

    for (i = 1; i < 5; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mean[i] += Mijk[i][j];
        }

        Mean[i] /= n;
    }

    //  Standard deviations calculation.

    for (i = 1; i < 5; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Tmp = Mijk[i][j] - Mean[i]; Stdev[i] += Tmp * Tmp;
        }

        Stdev[i] = (FLOAT)sqrt(Stdev[i] / (n - (FLOAT)1.0)); if (Stdev[i] < Eps) Stdev[i] = Eps;
    }

    for (i = 0; i < *c; i++) {
        Mx[i] = Mijk[1][i]; My[i] = Mijk[2][i]; Mz[i] = Mijk[3][i]; Mxyz[i] = Mijk[4][i];
    }

    //  Z-scores calculation.

    for (i = 1; i < 5; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mijk[i][j] = (Mijk[i][j] - Mean[i]) / Stdev[i];
        }
    }

    // Adjacency matrix calculation based on Gaussian kernel function.

    Mul = (FLOAT)0.5 / (*Sigma) / (*Sigma);

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) {
        j = (*c) * i;

        for (k = i + 1; k < *c; k++) if (N[k] > FLOAT_MIN) {
            ii = j + k; kk = i + (*c) * k; A[ii] = (FLOAT)0.0;

            for (l = 1; l < 5; l++) {
                Tmp = Mijk[l][i] - Mijk[l][k]; A[ii] += Tmp * Tmp;
            }

            A[ii] = A[kk] = (FLOAT)exp(-Mul * A[ii]);
        }
    }

E0: if (Stdev) free(Stdev);

    if (Mean) free(Mean);

    if (Mijk) {
        for (i = 0; i < 5; i++) {
            if (Mijk[i]) free(Mijk[i]);
        }

        free(Mijk);
    }
} // RLabelMomentsXYZ

// Returns normalised adjacency muatrix in R.

void RMergeLabels(INT    *n,      // Number of adjacency matrices.
                  double *A,      // Pointer to the adjacency matrices.
                  INT    *c,      // Number of clusters.
                  double *W,      // Weights.
                  double *L,      // Normalised adjacency matrix.
                  INT    *Error)  // Error code.
{
    FLOAT *D = NULL, p;
    INT   i, ii, j, k, kk, l, m;

    *Error = *n < 1; if (*Error) goto E0;
    
    *Error = *c < 2; if (*Error) goto E0;

    D = (FLOAT*)malloc(*c * sizeof(FLOAT));

    *Error = NULL == D; if (*Error) goto E0;

    // Adjacency matrix calculation.

    for (i = 0; i < *n; i++) {
        p = (FLOAT)0.0; j = (*c) * (*c) * i;

        for (k = 0; k < *c; k++) {
            l = j + (*c) * k;

            for (m = k + 1; m < *c; m++) {
                p += A[l + m];
            }
        }

        p = W[i] / p;

        for (k = 0; k < *c; k++) {
            l = j + (*c) * k;

            for (m = k + 1; m < *c; m++) {
                L[k + (*c) * m] = L[m + (*c) * k] += p * A[l + m];
            }
        }
    }

    // Degree matrix calculation.

    for (i = 0; i < *c; i++) {
        j = (*c) * i;

        D[i] = L[j]; for (k = 1; k < *c; k++) D[i] += L[j + k];

        D[i] = (FLOAT)sqrt((FLOAT)1.0 / D[i]);
    }

    // Normalized adjacency matrix calculation.

    for (i = 0; i < *c; i++) {
        j = (*c) * i;

        for (k = i + 1; k < *c; k++) {
            ii = j + k; kk = i + (*c) * k;

            L[ii] = L[kk] = L[ii] * D[i] * D[k];
        }
    }

E0: if (D) free(D);
} // RMergeLabels

}
