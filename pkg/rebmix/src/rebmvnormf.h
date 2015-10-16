#ifndef REBMVNORM_H_INCLUDED
#define REBMVNORM_H_INCLUDED

#include <stdlib.h>
#include <string.h>

#include "base.h"
#include "rebmixf.h"

class Rebmvnorm : public Rebmix {
public:
    // Methods.
    //  virtual int RoughEstimationKNN(FLOAT **Y, int k, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    //  virtual int RoughEstimationPW(FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    //    virtual int RoughEstimationH(int k, FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int ComponentDist(FLOAT *Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    //    virtual int EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    //virtual int EnhancedEstimationPW(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    //virtual int EnhancedEstimationH(int k, FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    //virtual int MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    //virtual int BayesClassificationKNN(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    //virtual int BayesClassificationPW(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    //virtual int BayesClassificationH(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    int DegreesOffreedom(int c, CompnentDistribution **MixTheta, int *M);
}; // Rebmvnorm

#endif