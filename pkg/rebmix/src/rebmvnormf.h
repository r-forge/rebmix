#ifndef REBMVNORM_H_INCLUDED
#define REBMVNORM_H_INCLUDED

#include <stdlib.h>
#include <string.h>

#include "base.h"
#include "rebmixf.h"

class Rebmvnorm : public Rebmix {
private:
    int ComponentConditionalDist(int i, FLOAT *Y, FLOAT *Cinv, CompnentDistribution *CmpTheta, FLOAT *CmpMrgDist);
public:
    // Methods.
    int RoughEstimationKNN(FLOAT **Y, int k, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int RoughEstimationPW(FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int RoughEstimationH(int k, FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int ComponentDist(FLOAT *Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    int EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int EnhancedEstimationPW(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int EnhancedEstimationH(int k, FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    int BayesClassificationKNN(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    int BayesClassificationPW(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    int BayesClassificationH(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    int DegreesOffreedom(int c, CompnentDistribution **MixTheta, int *M);
}; // Rebmvnorm

#endif