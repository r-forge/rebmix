#ifndef RNGMVNORM_H_INCLUDED
#define RNGMVNORM_H_INCLUDED

#include "base.h"
#include "rngmixf.h"

class Rngmvnorm : public Rngmix {
public:
    INT InvComponentDist(CompnentDistribution *CmpDist, INT j, FLOAT **Y);
}; // Rngmvnorm

#endif