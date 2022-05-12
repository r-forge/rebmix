#ifndef REBMVNORM_H_INCLUDED
#define REBMVNORM_H_INCLUDED

#include "base.h"
#include "rebmixf.h"

class Rebmvnorm : public Rebmix {
private:
    INT ComponentConditionalDist(INT i, INT j, FLOAT **Y, FLOAT *Cinv, CompnentDistribution *CmpTheta, FLOAT *CmpMrgDist);
public:
    // Methods.
    INT Initialize();
    INT RoughEstimationKNN(FLOAT **Y, INT k, FLOAT *h, FLOAT nl, INT m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    INT RoughEstimationKDE(FLOAT **Y, FLOAT *h, FLOAT nl, INT m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    INT RoughEstimationH(INT k, FLOAT **Y, FLOAT *h, FLOAT nl, INT m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    INT ComponentDist(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, INT *Outlier);
    INT LogComponentDist(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, INT *Outlier);
    INT EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    INT EnhancedEstimationKDE(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    INT EnhancedEstimationH(INT k, FLOAT **Y, FLOAT nl, FLOAT *h, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    INT MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    INT BayesClassificationKNN(FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    INT BayesClassificationKDE(FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    INT BayesClassificationH(INT k, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    INT DegreesOffreedom(INT c, CompnentDistribution **MixTheta, INT *M);
/// Panic Branislav
    INT EMInitialize();
/// End
}; // Rebmvnorm

#endif