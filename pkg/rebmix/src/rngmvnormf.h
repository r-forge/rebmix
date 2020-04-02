#ifndef RNGMVNORM_H_INCLUDED
#define RNGMVNORM_H_INCLUDED

#include <stdlib.h>
#include <string.h>

#include "base.h"
#include "rngmixf.h"

class Rngmvnorm : public Rngmix {
public:
    int InvComponentDist(CompnentDistribution *CmpDist, int j, FLOAT **Y);
}; // Rngmvnorm

#endif