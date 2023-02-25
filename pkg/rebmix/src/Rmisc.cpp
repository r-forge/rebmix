#include "base.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

extern "C" {

// Returns numbers of pixels and label moments for clustered image in R.

void RLabelMomentsXY(INT    *nx,     // Image width.
                     INT    *ny,     // Image height.
                     INT    *Zp,     // Pointer to the predictive cluster membership.
                     INT    *c,      // Number of clusters.
                     double *N,      // Numbers of pixels.
                     double *Mx,     // Raw moments Mx.
                     double *My,     // Raw moments My.  
                     double *Mxy,    // Raw moments Mxy.
                     double *L,      // Normalised affinity matrix calculation.
                     INT    *Error)  // Error code.
{
    FLOAT *D = NULL, **Mij = NULL, *Mean = NULL, *Stdev = NULL;
    INT   i, ii, j, k, kk, l, n = 0;

    *Error = *c < 1; if (*Error) goto E0;

    Mij = (FLOAT**)malloc(4 * sizeof(FLOAT*));

    *Error = NULL == Mij; if (*Error) goto E0;

    for (i = 0; i < 4; i++) {
        Mij[i] = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

        *Error = NULL == Mij[i]; if (*Error) goto E0;
    }

    Mean = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

    *Error = NULL == Mean; if (*Error) goto E0;

    Stdev = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

    *Error = NULL == Stdev; if (*Error) goto E0;

    D = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

    *Error = NULL == D; if (*Error) goto E0;
    
    // Numbers of pixels and raw clustered image moments calculation.   

    for (i = 0; i < *nx; i++) {
        j = (*ny) * i;

        for (k = 0; k < *ny; k++) {
            l = Zp[j + k] - 1;

            if (l >= 0) {
                Mij[0][l] += 1; ii = i + 1; kk = k + 1;

                Mij[1][l] += ii; Mij[2][l] += kk;

                Mij[3][l] += ii * kk;
            }
        }
    }

    for (i = 0; i < *c; i++) {
        N[i] = Mij[0][i]; Mx[i] = Mij[1][i]; My[i] = Mij[2][i]; Mxy[i] = Mij[3][i];
    }

    // Number of nonempty clusters calculation.
    
    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) n++;

    //  Means and standard deviations calculation.

    for (i = 1; i < 4; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mean[i] += Mij[i][j]; Stdev[i] += Mij[i][j] * Mij[i][j];
        }

        Mean[i] /= n; Stdev[i] = (FLOAT)sqrt((Stdev[i] - n * Mean[i] * Mean[i]) / (n - (FLOAT)1.0));
    }

    //  Z-scores calculation.

    for (i = 1; i < 4; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mij[i][j] = (Mij[i][j] - Mean[i]) / Stdev[i];
        }
    }

    // Affinity matrix calculation.

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) {
        j = (*c) * i;

        for (k = i + 1; k < *c; k++) if (N[k] > FLOAT_MIN) {
            ii = j + k; kk = i + (*c) * k; L[ii] = (FLOAT)0.0;
            
            for (l = 1; l < 4; l++) L[ii] += (Mij[l][i] - Mij[l][k]) * (Mij[l][i] - Mij[l][k]);

            L[ii] = L[kk] = (FLOAT)exp(-(FLOAT)0.5 * L[ii]);
        }
    }

    // Square root of inverse diagonal degree matrix calculation.

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) {
        j = (*c) * i;

        for (k = 0; k < *c; k++) if (N[k] > FLOAT_MIN) {
            D[i] += L[j + k];
        }

        D[i] = (FLOAT)sqrt((FLOAT)1.0 / D[i]);
    }

    // Normalised affinity matrix calculation.

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) {
        j = (*c) * i;

        for (k = i + 1; k < *c; k++) if (N[k] > FLOAT_MIN) {
            ii = j + k; kk = i + (*c) * k;

            L[ii] = L[kk] = L[ii] * D[i] * D[k];
        }
    }

E0: if (D) free(D);

    if (Stdev) free(Stdev);

    if (Mean) free(Mean);

    if (Mij) {
        for (i = 0; i < 4; i++) {
            if (Mij[i]) free(Mij[i]);
        }

        free(Mij);
    }
} // RLabelMomentsXY

// Returns numbers of voxels and label moments for clustered image in R.

void RLabelMomentsXYZ(INT    *nx,     // Image width.
                      INT    *ny,     // Image height.
                      INT    *nz,     // Image depth.
                      INT    *Zp,     // Pointer to the predictive cluster membership.
                      INT    *c,      // Number of clusters.
                      double *N,      // Numbers of voxels.
                      double *Mx,     // Raw moments Mx.
                      double *My,     // Raw moments My.
                      double *Mz,     // Raw moments Mz.
                      double *Mxyz,   // Raw moments Mxyz.  
                      double *Eud,    // Euclidean distances. 
                      INT    *Error)  // Error code.
{
    FLOAT **Mijk = NULL, *Mean = NULL, *Stdev = NULL, Sum;
    INT i, ii, j, k, kk, l, m, mm, n = 0;

    *Error = *c < 1; if (*Error) goto E0;

    Mijk = (FLOAT**)malloc(5 * sizeof(FLOAT*));

    *Error = NULL == Mijk; if (*Error) goto E0;

    for (i = 0; i < 5; i++) {
        Mijk[i] = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

        *Error = NULL == Mijk[i]; if (*Error) goto E0;
    }

    Mean = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

    *Error = NULL == Mean; if (*Error) goto E0;

    Stdev = (FLOAT*)calloc((size_t)(*c), sizeof(FLOAT));

    *Error = NULL == Stdev; if (*Error) goto E0;

    //  Numbers of voxels and raw clustered image moments calculation.   

    for (i = 0; i < *nz; i++) {
        j = (*nx) * (*ny) * i;

        for (k = 0; k < *nx; k++) {
            l = j + (*ny) * k;

            for (m = 0; m < *ny; m++) {
                n = Zp[l + m] - 1;

                if (n >= 0) {
                    Mijk[0][n] += 1; ii = i + 1; kk = k + 1; mm = m + 1;

                    Mijk[1][n] += kk; Mijk[2][n] += mm; Mijk[3][n] += ii;

                    Mijk[4][n] += ii * kk * mm;
                }
            }
        }
    }

    for (i = 0; i < *c; i++) {
        N[i] = Mijk[0][i]; Mx[i] = Mijk[1][i]; My[i] = Mijk[2][i]; Mz[i] = Mijk[3][i]; Mxyz[i] = Mijk[4][i];
    }

    //  Number of nonempty clusters calculation.

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) n++;

    //  Means and standard deviations calculation.

    for (i = 0; i < 5; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mean[i] += Mijk[i][j]; Stdev[i] += Mijk[i][j] * Mijk[i][j];
        }

        Mean[i] /= n; Stdev[i] = (FLOAT)sqrt((Stdev[i] - n * Mean[i] * Mean[i]) / (n - (FLOAT)1.0));
    }

    //  Z-scores calculation.

    for (i = 0; i < 5; i++) {
        for (j = 0; j < *c; j++) if (N[j] > FLOAT_MIN) {
            Mijk[i][j] = (Mijk[i][j] - Mean[i]) / Stdev[i];
        }
    }

    //  Euclidean distances calculation.

    Sum = (FLOAT)0.0;

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) {
        j = (*c) * i;

        for (k = i + 1; k < *c; k++) if (N[k] > FLOAT_MIN) {
            ii = j + k; kk = i + (*c) * k; Eud[ii] = (FLOAT)0.0;

            for (l = 1; l < 5; l++) Eud[ii] += (Mijk[l][i] - Mijk[l][k]) * (Mijk[l][i] - Mijk[l][k]);

            Eud[ii] = Eud[kk] = (FLOAT)sqrt(Eud[ii]); Sum += Eud[ii];
        }
    }

    //  Normalized Euclidean distances calculation.

    for (i = 0; i < *c; i++) if (N[i] > FLOAT_MIN) {
        j = (*c) * i;

        for (k = i + 1; k < *c; k++) if (N[k] > FLOAT_MIN) {
            ii = j + k; kk = i + (*c) * k;

            Eud[ii] = Eud[kk] = Eud[ii] / Sum;
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

}
