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
    INT                  n_;               // Number of observations.
    INT                  nr_;              // Total number of bins.
    INT                  nc_;              // Number of columns.
    FLOAT                **Y_;             // Dataset.
    INT                  cmax_;            // Maximum number of components.
    FLOAT                TOL_;             // Tolerance for EM algorithm.
    FLOAT                am_;              // Acceleration multiplier for EM algorithm.
    INT                  max_iter_;        // Maximum number of iterations of EM algorithm.
    INT                  K_;               // Number of bins for histogram EM algorithm.
    EmStrategyType_e     strategy_;        // EM strategy utilization.
    EmVariantType_e      variant_;         // Type of EM variant algorithm.
    EmAccelerationType_e accel_;           // Type of acceleration of standard EM algorithm.
/// Panic Branislav
    EmMergeCompsType_e   merge_;           // Merge components with zero weights.
/// End
    INT                  n_iter_;          // Number of iterations.
    INT                  c_;               // Number of components.
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
    INT Initialize(INT n, INT nr, INT nc, FLOAT **Y, INT cmax, INT length_pdf, INT length_Theta, INT *length_theta, FLOAT TOL, FLOAT am, INT max_iter, INT EM_K, EmStrategyType_e strategy, EmVariantType_e variant, EmAccelerationType_e accel);
    INT Transform(FLOAT **Y);
    INT MixtureDist(INT j, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    INT LogLikelihood(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *LogL);
    INT ExpectationStep();
    INT ConditionalStep();
    INT GoldenRatioSearch(FLOAT *am_opt);
    INT LineSearch(FLOAT *am_opt);
    INT EM();
    INT ECM();
    INT Run(INT *c, FLOAT *W, CompnentDistribution **MixTheta);
    virtual INT LogComponentDist(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    virtual INT UpdateMixtureParameters(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT am);
    virtual INT MaximizationStep();
}; // Emmix

class Emmvnorm : public Emmix {
public:
    // Constructor.
    INT LogComponentDist(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist);
    INT UpdateMixtureParameters(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT am);
    INT MaximizationStep();
}; // Emmvnorm

#endif


