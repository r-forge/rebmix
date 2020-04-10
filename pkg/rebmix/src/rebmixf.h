#ifndef REBMIXF_H_INCLUDED
#define REBMIXF_H_INCLUDED

#ifdef _MSC_VER
#pragma warning(disable: 4514)
#pragma warning(disable: 4820)
#endif

#include <stdlib.h>
#include <string.h>
#include "base.h"
#include "emf.h"

typedef enum {
    poHistogram,        // Histogram approach.
    poKDE,              // kernel density estimation.
    poKNearestNeighbour // K-nearest neighbour.
} PreprocessingType_e;

typedef enum {
    vtContinuous, // Continuous variable.
    vtDiscrete    // Ordered or non-ordered binary or discrete variable.
} VariablesType_e;

typedef enum {
    rtRigid, // Rigid restraints.
    rtLoose  // Loose restraints.
} PestraintsType_e;

typedef enum {
    icAIC,    // AIC - Akaike information criterion Akaike (1973).
    icAIC3,   // AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
    icAIC4,   // AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980).
    icAICc,   // AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989).
    icBIC,    // BIC - Bayesian information criterion Schwarz (1978).
    icCAIC,   // CAIC - Consistent Akaike information criterion Bozdogan (1987).
    icHQC,    // HQC - Hannan-Quinn information criterion Hannan & Quinn (1979).
    icMDL2,   // MDL2 - Minimum description length Liang et al.(1992).
    icMDL5,   // MDL5 - Minimum description length Liang et al.(1992).
    icAWE,    // AWE - Approximate weight of evidence criterion Banfield & Raftery (1993).
    icCLC,    // CLC - Classification likelihood criterion Biernacki & Govaert (1997).
    icICL,    // ICL - Integrated classification likelihood Biernacki et al.(1998).
    icPC,     // PC - Partition coefficient Bezdek (1981).
    icICLBIC, // ICL-BIC - Integrated classification likelihood criterion Biernacki et al.(1998).
    icD,      // D - Total of positive relative deviations Nagode & Fajdiga (2011).
    icSSE     // SSE - Sum of squares error Bishop (1998).
} InformationCriterionType_e;

typedef struct roughparametertype {
    FLOAT h;   // Mode class width;
    FLOAT ym;  // Mode position.
    FLOAT flm; // Component conditional empirical density.
    FLOAT klm; // Component conditional total number of observations.
} RoughParameterType;

class Rebmix : public Base {
    // Methods.
    int Golden();
    int GlobalModeKNN(int *m, FLOAT **Y, int *I);
    int GlobalModeKDE(int *m, FLOAT **Y, int *I);
    int GlobalModeH(int *m, int k, FLOAT **Y, int *O);
    int REBMIXKNN();
    int REBMIXKDE();
    int REBMIXH();
    #if (_MAINTAIN_SWITCH)
    int ReadDataFile();
    int WriteDataFile();
    #endif
public:
    // Input members.
    FLOAT                      p_value_;       // Probability of obtaining a result equal to or "more extreme" than what was actually observed.
    FLOAT                      min_dist_mul_;  // Minimum distance multiplier.
    FLOAT                      var_mul_;       // Variance multiplier.
    int                        kmax_;          // Maximum number of nonempty bins.
    FLOAT                      ChiSqr_;        // Critical Chi square value for outlier detection and p = 2.0 * p_value_.
    char                       *curr_;         // Path to the currently open data file.
    int                        o_;             // Number of paths.
    char                       **open_;        // Paths to open data files.
    char                       *save_;         // Path to the save data file.
    PreprocessingType_e        Preprocessing_; // Preprocessing type.
    int                        cmax_;          // Maximum number of components.
    int                        cmin_;          // Maximum number of components.
    InformationCriterionType_e Criterion_;     // Information criterion type.
    VariablesType_e            *Variables_;    // Types of variables.
    CompnentDistribution       *IniTheta_;     // Initial component parameters.
    int                        length_K_;      // Length of K_.
    int                        *K_;            // Numbers of bins v or numbers of nearest neighbours k.
    FLOAT                      *y0_;           // Origins.
    FLOAT                      *ymin_;         // Minimum observations.
    FLOAT                      *ymax_;         // Maximum observations.
    FLOAT                      ar_;            // Acceleration rate.
    PestraintsType_e           Restraints_;    // Restraints type.
/// Panic Branislav: fields for EM algorithm.
    Emmix                      *EM_;           // Object of class Emmix.
    FLOAT                      EM_TOL_;        // Tolerance for EM algorithm.
    FLOAT                      EM_am_;         // Acceleration multiplier for EM algorithm.
    int                        EM_max_iter_;   // Maximum number of iterations of EM algorithm.    
    EmStrategyType_e           EM_strategy_;   // EM strategy utilization.
    EmVariantType_e            EM_variant_;    // Type of EM variant algorithm.
    EmAccelerationType_e       EM_accel_;      // Type of acceleration of standard EM algorithm.
    MixtureParameterType       *OptMixTheta_;  // Best mixture parameters.
/// End
    // Input members.
    int                        n_;             // Number of observations.
    FLOAT                      **Y_;           // Dataset.
    FLOAT                      **X_;           // Temporary dataset.
    // Output members.
    FLOAT                      *W_;            // Component weights.
    CompnentDistribution       **MixTheta_;    // Mixture parameters.
    SummaryParameterType       summary_;       // Summary.
    int                        opt_length_;    // Length of opt_c_, opt_IC_, opt_logL_ and opt_D_.
    int                        *opt_c_;        // Numbers of components for optimal v or for optimal k.
    FLOAT                      *opt_IC_;       // Information criteria for optimal v or for optimal k.
    FLOAT                      *opt_logL_;     // Log-likelihoods for optimal v or for optimal k.
    FLOAT                      *opt_D_;        // Totals of positive relative deviations for optimal v or for optimal k.
    int                        all_length_;    // Length of all_K_ and all_IC_.
    int                        *all_I_;        // Information on processed numbers of bins v or processed numbers of nearest neighbours k. 0 if not processed, 1 if processed, 2 if error.
    int                        *all_K_;        // All processed numbers of bins v or all processed numbers of nearest neighbours k.
    FLOAT                      *all_IC_;       // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
    AdditionalParameterType    additional_;    // Additional parameters.
/// Panic Branislav.
    int                        n_iter_;        // Number of iterations performed by EM algorithm for optimal mixture selected.
    int                        n_iter_sum_;    // Total number of iterations performed by EM algorithm. 
/// End
    // Constructor.
    Rebmix();
    // Destructor.
    virtual ~Rebmix();
    // Methods.
    virtual int Initialize();
    int PreprocessingKNN(int k, FLOAT *h, FLOAT **Y);
    int PreprocessingKDE(FLOAT *h, FLOAT **Y);
    int PreprocessingH(FLOAT *h, FLOAT *y0, int *k, FLOAT **Y);
    int PreprocessingH(FLOAT *h, FLOAT *y0, int *k, FLOAT **Y, int *State);
    virtual int RoughEstimationKNN(FLOAT **Y, int k, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int RoughEstimationKDE(FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int RoughEstimationH(int k, FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int ComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, int *Outlier);
    virtual int LogComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, int *Outlier);
    virtual int EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int EnhancedEstimationKDE(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int EnhancedEstimationH(int k, FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    virtual int BayesClassificationKNN(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int BayesClassificationKDE(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int BayesClassificationH(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int DegreesOffreedom(int c, CompnentDistribution **MixTheta, int *M);
/// Panic Branislav: method for invoking EM algorithm.
    virtual int EMInitialize();
    virtual int EMRun(int c, FLOAT *W, CompnentDistribution **MixTheta);
/// End 
    int MixtureDist(int j, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    int MixtureDist(FLOAT logV, int j, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    int InformationCriterionKNN(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int InformationCriterionKDE(FLOAT logV, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int InformationCriterionH(FLOAT logV, int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int CombineComponents(int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *tau, int *F, int *T, FLOAT *EN, FLOAT *ED);
    int REBMIX();
    #if (_MAINTAIN_SWITCH)
    int RunTemplateFile(char *file);
    #endif
}; // Rebmix

#endif
