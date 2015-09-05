#include <math.h>
#include <float.h> 

#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "base.h"
#include "mvnorm.h"

static int   NDevISet = 0;
static FLOAT NDevVSet = (FLOAT)0.0;

int Rngmvnorm::InvComponentDist(CompnentDistribution *CmpDist, FLOAT *Y)
{
    FLOAT C[4];
    FLOAT *y, Sum;
    int   i, j;
    int   Error = 0;

    y = (FLOAT*)malloc(d_ * sizeof(FLOAT));

    Error = NULL == y; if (Error) goto E0;

    for (i = 0; i < d_; i++) {
        if (Trigger_) {
            Trigger_ = 0;

            Error = Choldc(d_, CmpDist->Theta_[1], CmpDist->Theta_[2]);

            if (Error) goto E0;
        }

        if (NDevISet == 0) {
            do {
                C[0] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;
                C[1] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;

                C[2] = C[0] * C[0] + C[1] * C[1];
            } while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

            C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * log(C[2]) / C[2]);

            y[i] = C[3] * C[0];

            NDevISet = 1; NDevVSet = C[3] * C[1];
        }
        else {
            y[i] = NDevVSet; NDevISet = 0;
        }
    }

    for (i = 0; i < d_; i++) {
        Sum = (FLOAT)0.0;

        for (j = 0; j <= i; j++) {
            Sum += CmpDist->Theta_[2][i * d_ + j] * y[j];
        }

        Y[i] = Sum + CmpDist->Theta_[0][i];
    }

E0: if (y) free(y);

    return Error;
} // InvComponentDist

// Returns component p.d.f or c.d.f.

int Rebmvnorm::ComponentDist(FLOAT              *Y,        // Pointer to the input point [y0,...,yd-1].
                          CompnentDistribution *CmpTheta, // Component parameters.
                          FLOAT              *CmpDist)  // Component distribution.
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

    *CmpDist = (FLOAT)exp(-y) / (CmpTheta->Theta_[4][0]);

    return Error;
} // ComponentDist

int Rebmvnorm::DegreesOffreedom(int                c,          // Number of components.
                             CompnentDistribution **MixTheta, // Mixture parameters.
                             int                *M)         // Degrees of freedom.
{
    int Error = 0;

    *M = c * (d_ + (d_ + 1) * d_ / 2) + c - 1;

    return Error;
} // DegreesOffreedom
