#ifndef REBMVNORM_H_INCLUDED
#define REBMVNORM_H_INCLUDED

#include <stdlib.h>
#include <string.h>

#include "base.h"
#include "rebmixf.h"

class Rebmvnorm : public Rebmix {
private:
    int ComponentConditionalDist(int i, int j, FLOAT **Y, FLOAT *Cinv, CompnentDistribution *CmpTheta, FLOAT *CmpMrgDist);
public:
    // Methods.
    int Initialize();
    int RoughEstimationKNN(FLOAT **Y, int k, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int RoughEstimationKDE(FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int RoughEstimationH(int k, FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int ComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, int *Outlier);
    int LogComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, int *Outlier);
    int EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int EnhancedEstimationKDE(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int EnhancedEstimationH(int k, FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    int MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    int BayesClassificationKNN(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    int BayesClassificationKDE(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    int BayesClassificationH(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    int DegreesOffreedom(int c, CompnentDistribution **MixTheta, int *M);
    /// Panic Branislav: Method for invoking EM algorithm. ///
    int ExpectationMaximizationStep(int c, FLOAT *W, CompnentDistribution **MixTheta, int *n_iter);
    /// End                                                ///
}; // Rebmvnorm

#endif