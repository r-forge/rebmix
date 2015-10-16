#include <math.h>
#include <float.h> 

#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "base.h"
#include "rebmvnormf.h"

// Returns component p.d.f or c.d.f.

int Rebmvnorm::ComponentDist(FLOAT                *Y,        // Pointer to the input point [y0,...,yd-1].
                             CompnentDistribution *CmpTheta, // Component parameters.
                             FLOAT                *CmpDist)  // Component distribution.
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

int Rebmvnorm::DegreesOffreedom(int                  c,          // Number of components.
                                CompnentDistribution **MixTheta, // Mixture parameters.
                                int                  *M)         // Degrees of freedom.
{
    int Error = 0;

    *M = c * (d_ + (d_ + 1) * d_ / 2) + c - 1;

    return Error;
} // DegreesOffreedom
