#ifndef REBMIXF_H_INCLUDED
#define REBMIXF_H_INCLUDED

#include "base.h"
#include "emf.h"

#include <stdlib.h>
#include <string.h>

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
    FLOAT h;     // Mode class width;
    FLOAT ym;    // Mode position.
    FLOAT ymean; // Mean position.
    FLOAT ymin;  // Minimum position.
    FLOAT ymax;  // Maximum position.
    FLOAT flm;   // Component conditional empirical density.
    FLOAT klm;   // Component conditional total number of observations.
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
    int REBMIXK();
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
    int                        cmin_;          // Minimum number of components.
    InformationCriterionType_e Criterion_;     // Information criterion type.
    VariablesType_e            *Variables_;    // Types of variables.
    CompnentDistribution       *IniTheta_;     // Initial component parameters.
    int                        length_K_;      // Length of K_.
    int                        *K_;            // Numbers of bins v or numbers of nearest neighbours k.
    FLOAT                      *y0_;           // Origins.
    FLOAT                      *ymin_;         // Minimum observations.
    FLOAT                      *ymax_;         // Maximum observations.
    FLOAT                      *h_;            // Sides of the hypersquare.
    FLOAT                      ar_;            // Acceleration rate.
    PestraintsType_e           Restraints_;    // Restraints type.
/// Panic Branislav
    Emmix                      *EM_;           // Object of class Emmix.
    FLOAT                      EM_TOL_;        // Tolerance for EM algorithm.
    FLOAT                      EM_am_;         // Acceleration multiplier for EM algorithm.
    int                        EM_max_iter_;   // Maximum number of iterations of EM algorithm.
    int                        EM_K_;          // Number of bins for histogram EM algorithm.
    EmStrategyType_e           EM_strategy_;   // EM strategy utilization.
    EmVariantType_e            EM_variant_;    // Type of EM variant algorithm.
    EmAccelerationType_e       EM_accel_;      // Type of acceleration of standard EM algorithm.
    MixtureParameterType       *OptMixTheta_;  // Best mixture parameters.
/// End
    // Input members.
    int                        n_;             // Number of observations.
    int                        nr_;            // Number of rows.
    int                        nc_;            // Number of columns.
    FLOAT                      **Y_;           // Dataset.
    int                        Y_type_;        // Dataset type.
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
/// Panic Branislav
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
    int PreprocessingH(FLOAT *h, FLOAT *y0, FLOAT *ymin, FLOAT *ymax, int *k, FLOAT **Y);
    int PreprocessingH(FLOAT *h, FLOAT *y0, FLOAT *ymin, FLOAT *ymax, int *k, FLOAT **Y, int *State);
    virtual int RoughEstimationKNN(FLOAT **Y, int k, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int RoughEstimationKDE(FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int RoughEstimationH(int k, FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int ComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, int *Outlier);
    virtual int LogComponentDist(int j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, int *Outlier);
    virtual int EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int EnhancedEstimationKDE(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int EnhancedEstimationH(int k, FLOAT **Y, FLOAT nl, FLOAT *h, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    virtual int BayesClassificationKNN(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int BayesClassificationKDE(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int BayesClassificationH(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int DegreesOffreedom(int c, CompnentDistribution **MixTheta, int *M);
/// Panic Branislav
    virtual int EMInitialize();
    virtual int EMRun(int *c, FLOAT *W, CompnentDistribution **MixTheta);
/// End 
    int MixtureDist(int j, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    int MixtureDist(FLOAT logV, int j, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    int InformationCriterionKNN(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int InformationCriterionKDE(FLOAT logV, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int InformationCriterionH(FLOAT logV, int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int CombineComponents(int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *tau, int *F, int *T, FLOAT *EN, FLOAT *ED);
/// Panic Branislav
    int CombineComponentsDemp(int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *tau, int *F, int *T, FLOAT *EN, FLOAT *ED);
/// End
    int REBMIX();
    int Set(char  **Preprocessing,    // Preprocessing type.
            int   *cmax,              // Maximum number of components.
            int   *cmin,              // Minimum number of components.
            char  **Criterion,        // Information criterion type.
            int   *d,                 // Number of independent random variables.
            char  **Variables,        // Types of variables.
            int   *length_pdf,        // Length of pdf.
            char  **pdf,              // Parametric family types.
            int   *length_Theta,      // Length of Theta.
            int   *length_theta,      // Length of Theta[i].
            FLOAT *Theta,             // Component parameters.
            int   *length_K,          // Length of K.
            int   *K,                 // Numbers of bins v or numbers of nearest neighbours k.
            int   *length_y0,         // Length of y0.
            FLOAT *y0,                // Origins.
            int   *length_ymin,       // Length of ymin.
            FLOAT *ymin,              // Minimum observations.
            int   *length_ymax,       // Length of ymax.
            FLOAT *ymax,              // Maximum observations.
            int   *length_h,          // Length of h.
            FLOAT *h,                 // Sides of the hypersquare.
            FLOAT *ar,                // Acceleration rate.
            char  **Restraints,       // Restraints type.
            int   *n,                 // Number of observations.
            FLOAT *Y,                 // Dataset.
            int   *Y_type,            // Dataset type. 
            char  **EMStrategy,       // Strategy for EM algorithm.
            char  **EMVariant,        // EM algorithm variant.
            char  **EMAcceleration,   // Acceleration for the standard EM algorithm.
            FLOAT *EMTolerance,       // Tolerance for EM algortihm.
            FLOAT *EMAccelerationMul, // Acceleration rate for Em algorithm.
            int   *EMMaxIter,         // Maximum number of iterations in EM algorithm.
            int   *EMK);              // Number of bins for histogram EM algorithm.
/// Panic Branislav
    int Get(int   *n_iter,         // Number of iterations for optimal case.
            int   *n_iter_sum,     // Number of iterations in whole run.
/// End
            int   *summary_k,      // Optimal v or optimal k.
            FLOAT *summary_h,      // Optimal class widths of length d.
            FLOAT *summary_y0,     // Optimal origins of length d.
            FLOAT *summary_ymin,   // Optimal minimum observations of length d.
            FLOAT *summary_ymax,   // Optimal maximum observations of length d.
            FLOAT *summary_IC,     // Optimal information criterion.
            FLOAT *summary_logL,   // Log-likelihood.
            int   *summary_M,      // Degrees of freedom.
            int   *summary_c,      // Optimal number of components.
            FLOAT *W,              // Component weights.
            FLOAT *theta1,         // Component parameters.
            FLOAT *theta2,         // Component parameters.
            FLOAT *theta3,         // Component parameters.
            int   *opt_length,     // Length of opt_c, opt_IC, opt_logL and opt_D.
            int   *opt_c,          // Numbers of components for optimal v or for optimal k.
            FLOAT *opt_IC,         // Information criteria for optimal v or for optimal k.
            FLOAT *opt_logL,       // Log-likelihoods for optimal v or for optimal k.
            FLOAT *opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
            int   *all_length,     // Length of all_K and all_IC.
            int   *all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
            FLOAT *all_IC);        // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
    #if (_MAINTAIN_SWITCH)
    int RunTemplateFile(char *file);
    #endif
}; // Rebmix

#endif
