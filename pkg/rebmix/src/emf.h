/*
 *
 * Helper functions(definitions) for the Expectation-Maximization(EM) algorithm for the finite mixture models.
 *
 * Author: Branislav Panic
 *
*/

#ifndef EMF_H_INCLUDED
#define EMF_H_INCLUDED

#ifdef _MSC_VER
#pragma warning(disable: 4514)
#pragma warning(disable: 4820)
#endif

#include <stdlib.h>
#include <string.h>
#include "base.h"

class Em : public Base {
public:
    FLOAT                TOL_;             // Tolerance for stoping criteria.
    FLOAT                ar_;              // Acceleration rate.
    int                  MAX_ITER_;        // Maximum allowed iterations.
    int                  Initialized_;     // Status check for initialization.
    int                  C_;               // Number of components.
    FLOAT                *W_;              // Component weights.
    FLOAT                *dW_;             // Update component weights.
    CompnentDistribution **MixTheta_;      // Mixture parameters.    
    CompnentDistribution **dMixTheta_;     // Update mixture parameters.
    SummaryParameterType summary_;         // Summary.
    EmVariantType_e      variant_t_;       // Type of EM algorithm (variant).
    EmAccelerationType_e acceleration_t_;  // Type of acceleration.
    int                  n_iter_;          // Number of iterations.
    // Constructor.
    Em();
    // Destructor.
    virtual ~Em();
    int Initialize(int C, FLOAT *iniW, CompnentDistribution **IniMixTheta, FLOAT TOL, int MAX_ITER, EmVariantType_e algType, EmAccelerationType_e accelType, FLOAT accel_mul);
    virtual int LogComponentDist(FLOAT *Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    int MixtureDist(FLOAT *Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    int ExpectationStep(FLOAT **Y, int n, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P);
    int ConditionalStep(int n, int c, FLOAT **P);
    virtual int UpdateMixtureParameters(int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT ar);
    int GoldenRatioSearch(FLOAT **Y, int N, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT *ar_opt);
    int LineSearch(FLOAT **Y, int N, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT *ar_opt);
    virtual int MaximizationStep(FLOAT **Y, int n, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P);
    int LogLikelihood(FLOAT **Y, int n, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *LogL);
    int EM(FLOAT **Y, int n);
    int ECM(FLOAT **Y, int n);
    int Run(FLOAT **Y, int n);
}; // EM

class GaussianDiagMixture: public Em{
public:
    int LogComponentDist(FLOAT *Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    virtual int UpdateMixtureParameters(int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT ar);
    int MaximizationStep(FLOAT **Y, int n, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P);
}; // GaussianDiagMixture

#endif


