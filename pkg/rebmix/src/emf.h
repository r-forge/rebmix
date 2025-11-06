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
private:
    static constexpr INT ITER_NBR_LINE_SEARCH = 10; // Hardcoded number of 
    static constexpr INT ITER_NBR_GOLDEN_SEARCH = 20; // Ha
    static constexpr FLOAT MAX_SEARCH_RANGE = (FLOAT)2.0;
    static constexpr FLOAT MIN_SEARCH_RANGE = (FLOAT)1.0;
    static constexpr FLOAT INC_LINE_SEARCH = (FLOAT)0.1;
    static constexpr FLOAT GOLDEN_SEARCH_TOL = (FLOAT)0.001;

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
    EmAccelParamType_e   accel_eq_;        // Calculation for optimal acceleration parameter.
    EmLikelihoodEstimateType_e loglest_;   // Likelihood estimation type.
    EmConvergenceType_e  toltype_;         // Convergence criterion for loglikelihood estimate. 
/// End
    INT                  n_iter_;          // Number of iterations.
    INT                  c_;               // Number of components.
    FLOAT                aam_;             // Average acceleration rate.
    FLOAT                aam_step_;        // Step for linear/exponential growth/decay.
    FLOAT                *W_;              // Component weights.
    CompnentDistribution **MixTheta_;      // Mixture parameters.    
    FLOAT                *dW_;             // Update component weights.
    CompnentDistribution **dMixTheta_;     // Update mixture parameters.
    FLOAT                *TW_;             // Temporary component weights.
    CompnentDistribution **TMixTheta_;     // Temporary mixture parameters.
    SummaryParameterType summary_;         // Summary.
    FLOAT                **P_;             // Pointer to posterior probabilities.
    FLOAT                LogL_[3];         // Last three estimates of LogL for acceleration;

    // Constructor.
    Emmix();
    // Destructor.
    virtual ~Emmix();
    INT Initialize(INT n, INT nr, INT nc, FLOAT **Y, 
                   INT cmax, INT length_pdf, INT length_Theta, INT *length_theta, 
                   FLOAT TOL, FLOAT am, INT max_iter, INT EM_K, 
                   EmStrategyType_e strategy, EmVariantType_e variant, EmAccelerationType_e accel, 
                   EmAccelParamType_e accel_eq, EmLikelihoodEstimateType_e loglest, EmConvergenceType_e  toltype);
    INT Transform(FLOAT **Y);
    INT MixturePdf(INT j, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixPdf);
    INT LogLikelihood(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *LogL);
    INT ExpectationStep(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P);
    INT ConditionalStep(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P);
    INT StochasticStep(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P);
    INT MaximizationStep(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT *TW, CompnentDistribution **TMixTheta, FLOAT **P);
    INT OptSLen(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P, FLOAT *am_opt);
    INT CheckSLenConv(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT *am_opt);
    INT CheckSLenConvFast(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT *am_opt);
    INT GoldenRatioSearch(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT *TW, CompnentDistribution **TMixTheta, FLOAT *am_opt);
    INT LineSearch(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT *TW, CompnentDistribution **TMixTheta, FLOAT *am_opt);
    INT EM();
    INT ECM();
    INT SEM();
    INT ECM_EM();
    INT SEM_EM();
    INT Run(INT *c, FLOAT *W, CompnentDistribution **MixTheta);
    INT EstimateAlphaNaive(INT i, FLOAT am, FLOAT *am_opt);
    INT AitkenAcceleratedLogL(FLOAT *LogL, FLOAT *LogLPrev);
    INT Converged(FLOAT LogLNew, FLOAT LogLOld, INT N);
    virtual INT LogComponentPdf(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpPdf);
    virtual INT UpdateMixParams(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT am);
    virtual INT UpdateMixParamsSt(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT am);
    virtual INT UpdateMixParamsSqr(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT am);
    virtual INT Resids(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P);
    virtual INT OptSLen1(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P, FLOAT *am_opt);
    virtual INT OptSLen2(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P, FLOAT *am_opt);
    virtual INT OptSLen3(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P, FLOAT *am_opt);
    virtual INT F(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT **P);
}; // Emmix

class Emmvnorm : public Emmix {
public:
    // Constructor.
    INT LogComponentPdf(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpPdf);
    INT UpdateMixParams(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT am);
    INT UpdateMixParamsSt(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT am);
    INT UpdateMixParamsSqr(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT am);
    INT Resids(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P);
    INT OptSLen1(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P, FLOAT *am_opt);
    INT OptSLen2(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P, FLOAT *am_opt);
    INT OptSLen3(INT *c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *RW, CompnentDistribution **RMixTheta, FLOAT *VW, CompnentDistribution **VMixTheta, FLOAT **P, FLOAT *am_opt);
    INT F(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT **P);
}; // Emmvnorm

#endif


