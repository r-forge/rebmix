#include "rngmvnormf.h"

#include <math.h>

#if (_MAINTAIN_SWITCH)
#include <ctype.h>
#include <stdio.h>
#endif

static int   NDevISet = 0;
static FLOAT NDevVSet = (FLOAT)0.0;

int Rngmvnorm::InvComponentDist(CompnentDistribution *CmpDist, int j, FLOAT **Y)
{
    FLOAT C[4];
    FLOAT *y, Sum;
    int   i, k;
    int   Error = 0;

    y = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = NULL == y; if (Error) goto E0;

    for (i = 0; i < length_pdf_; i++) {
        if (Trigger_) {
            Trigger_ = 0;

            Error = Choldc(length_pdf_, CmpDist->Theta_[1], CmpDist->Theta_[2]);

            if (Error) goto E0;
        }

        if (NDevISet == 0) {
            do {
                C[0] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;
                C[1] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;

                C[2] = C[0] * C[0] + C[1] * C[1];
            } while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

            C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * (FLOAT)log(C[2]) / C[2]);

            y[i] = C[3] * C[0];

            NDevISet = 1; NDevVSet = C[3] * C[1];
        }
        else {
            y[i] = NDevVSet; NDevISet = 0;
        }
    }

    for (i = 0; i < length_pdf_; i++) {
        Sum = (FLOAT)0.0;

        for (k = 0; k <= i; k++) {
            Sum += CmpDist->Theta_[2][i * length_pdf_ + k] * y[k];
        }

        Y[i][j] = Sum + CmpDist->Theta_[0][i];
    }

E0: if (y) free(y);

    return Error;
} // InvComponentDist
