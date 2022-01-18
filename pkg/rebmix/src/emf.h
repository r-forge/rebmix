/*
 *
 * Helper functions(definitions) for the Expectation-Maximization(EM) algorithm for the finite mixture models.
 *
 * Author: Branislav Panic
 *
*/

#ifndef EMF_H_INCLUDED
#define EMF_H_INCLUDED

#include "base.h"

class Emmix : public Base {
public:
    // Members.
    int                  n_;               // Number of observations.
    int                  k_;               // Total number of bins.
    FLOAT                **Y_;             // Dataset.
    int                  cmax_;            // Maximum number of components.
    FLOAT                TOL_;             // Tolerance for EM algorithm.
    FLOAT                am_;              // Acceleration multiplier for EM algorithm.
    int                  max_iter_;        // Maximum number of iterations of EM algorithm.
    int                  K_;               // Number of bins for histogram EM algorithm.
    EmStrategyType_e     strategy_;        // EM strategy utilization.
    EmVariantType_e      variant_;         // Type of EM variant algorithm.
    EmAccelerationType_e accel_;           // Type of acceleration of standard EM algorithm.
/// Panic Branislav
    EmMergeCompsType_e   merge_;           // Merge components with zero weight.s
/// End
    int                  n_iter_;          // Number of iterations.
    int                  c_;               // Number of components.
    FLOAT                *W_;              // Component weights.
    CompnentDistribution **MixTheta_;      // Mixture parameters.    
    FLOAT                *dW_;             // Update component weights.
    CompnentDistribution **dMixTheta_;     // Update mixture parameters.
    SummaryParameterType summary_;         // Summary.
    FLOAT                **P_;             // Pointer to posterior probabilities.
    // Constructor.
    Emmix();
    // Destructor.
    virtual ~Emmix();
    int Initialize(int n, FLOAT **Y, int cmax, int length_pdf, int length_Theta, int *length_theta, FLOAT TOL, FLOAT am, int max_iter, int EM_K, EmStrategyType_e strategy, EmVariantType_e variant, EmAccelerationType_e accel);
    int Transform(FLOAT **Y);
    int MixtureDist(int j, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    int LogLikelihood(int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *LogL);
    int ExpectationStep();
    int ConditionalStep();
    int GoldenRatioSearch(FLOAT *am_opt);
    int LineSearch(FLOAT *am_opt);
    int EM();
    int ECM();
    int Run(int *c, FLOAT *W, CompnentDistribution **MixTheta);
    virtual int LogComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    virtual int UpdateMixtureParameters(int *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT am);
    virtual int MaximizationStep();
}; // Emmix

class Emmvnorm : public Emmix {
public:
    // Constructor.
    int LogComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    int UpdateMixtureParameters(int *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT am);
    int MaximizationStep();
}; // Emmvnorm

#endif


