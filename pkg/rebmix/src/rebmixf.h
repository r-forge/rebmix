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
    mtAll,         // The modes are determined in decreasing magnitude from all observations.
    mtOutliers,    // The modes are determined in decreasing magnitude from outliers only. Meanwhile, some outliers are reclassified as inliers. Eventually, when all observations are inliers, the procedure is concluded. 
    mtOutliersPlus // The modes are determined in decreasing magnitude from outliers only. Meanwhile, some outliers are reclassified as inliers. Eventually, when all observations are inliers, they are converted to outliers, and the procedure is repeated.
} ModeType_e;

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
    FLOAT fm;    // Component conditional empirical density.
    FLOAT km;    // Component conditional total number of observations.
} RoughParameterType;

class Rebmix : public Base {
    // Methods.
    INT Golden();
    INT GlobalModeKNN(INT *m, FLOAT **Y, INT *I);
    INT GlobalModeKDE(INT *m, FLOAT **Y, INT *I);
    INT GlobalModeH(INT *m, INT k, FLOAT **Y, INT *O);
    INT REBMIXKNN();
    INT REBMIXKDE();
    INT REBMIXH();
    INT REBMIXK();
    #if (_MAINTAIN_SWITCH)
    INT ReadDataFile();
    INT WriteDataFile();
    #endif
public:
    // Input members.
    FLOAT                      p_value_;       // Probability of obtaining a result equal to or "more extreme" than what was actually observed.
    FLOAT                      var_mul_;       // Variance multiplier.
    INT                        kmax_;          // Maximum number of nonempty bins.
    FLOAT                      ChiSqr_;        // Critical Chi square value for outlier detection and p = 2.0 * p_value_.
    char                       *curr_;         // Path to the currently open data file.
    INT                        o_;             // Number of paths.
    char                       **open_;        // Paths to open data files.
    char                       *save_;         // Path to the save data file.
    PreprocessingType_e        Preprocessing_; // Preprocessing type.
    INT                        cmax_;          // Maximum number of components.
    INT                        cmin_;          // Minimum number of components.
    InformationCriterionType_e Criterion_;     // Information criterion type.
    VariablesType_e            *Variables_;    // Types of variables.
    CompnentDistribution       *IniTheta_;     // Initial component parameters.
    INT                        length_K_;      // Length of K_.
    INT                        *K_;            // Numbers of bins v or numbers of nearest neighbours k.
    FLOAT                      *y0_;           // Origins.
    FLOAT                      *ymin_;         // Minimum observations.
    FLOAT                      *ymax_;         // Maximum observations.
    FLOAT                      *h_;            // Sides of the hypersquare.
    FLOAT                      ar_;            // Acceleration rate.
    PestraintsType_e           Restraints_;    // Restraints type.
    ModeType_e                 Mode_;          // Mode type.
/// Panic Branislav
    Emmix                      *EM_;           // Object of class Emmix.
    FLOAT                      EM_TOL_;        // Tolerance for EM algorithm.
    FLOAT                      EM_am_;         // Acceleration multiplier for EM algorithm.
    INT                        EM_max_iter_;   // Maximum number of iterations of EM algorithm.
    INT                        EM_K_;          // Number of bins for histogram EM algorithm.
    EmStrategyType_e           EM_strategy_;   // EM strategy utilization.
    EmVariantType_e            EM_variant_;    // Type of EM variant algorithm.
    EmAccelerationType_e       EM_accel_;      // Type of acceleration of standard EM algorithm.
    MixtureParameterType       *OptMixTheta_;  // Best mixture parameters.
/// End
    // Input members.
    INT                        n_;             // Number of observations.
    INT                        nr_;            // Number of rows.
    INT                        nc_;            // Number of columns.
    FLOAT                      **Y_;           // Dataset.
    INT                        Y_type_;        // Dataset type.
    FLOAT                      **X_;           // Temporary dataset.
    FLOAT                      **Z_;           // Temporary dataset.
    // Output members.
    FLOAT                      *W_;            // Component weights.
    CompnentDistribution       **MixTheta_;    // Mixture parameters.
    SummaryParameterType       summary_;       // Summary.
    INT                        opt_length_;    // Length of opt_c_, opt_IC_, opt_logL_, opt_Dmin_ and opt_D_.
    INT                        *opt_c_;        // Numbers of components for optimal v or for optimal k.
    FLOAT                      *opt_IC_;       // Information criteria for optimal v or for optimal k.
    FLOAT                      *opt_logL_;     // Log-likelihoods for optimal v or for optimal k.
    FLOAT                      *opt_Dmin_;     // Dmin for optimal v or for optimal k.
    FLOAT                      *opt_D_;        // Totals of positive relative deviations for optimal v or for optimal k.
    INT                        all_length_;    // Length of all_K_ and all_IC_.
    INT                        *all_I_;        // Information on processed numbers of bins v or processed numbers of nearest neighbours k. 0 if not processed, 1 if processed, 2 if error.
    INT                        *all_K_;        // All processed numbers of bins v or all processed numbers of nearest neighbours k.
    FLOAT                      *all_IC_;       // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
    AdditionalParameterType    additional_;    // Additional parameters.
/// Panic Branislav
    INT                        n_iter_;        // Number of iterations performed by EM algorithm for optimal mixture selected.
    INT                        n_iter_sum_;    // Total number of iterations performed by EM algorithm. 
/// End
    // Constructor.
    Rebmix();
    // Destructor.
    virtual ~Rebmix();
    // Methods.
    virtual INT Initialize();
    INT PreprocessingKNN(INT k, FLOAT *h, FLOAT *Rm, FLOAT **Y);
    INT PreprocessingKDE(FLOAT *h, FLOAT **Y);
    INT PreprocessingH(FLOAT *h, FLOAT *y0, FLOAT *ymin, FLOAT *ymax, INT *k, FLOAT **Y);
    INT PreprocessingH(FLOAT *h, FLOAT *y0, FLOAT *ymin, FLOAT *ymax, INT *k, FLOAT **Y, INT *State);
    virtual INT RoughEstimationKNN(FLOAT **Y, INT k, FLOAT *h, FLOAT Rm, FLOAT nl, INT m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual INT RoughEstimationKDE(FLOAT **Y, FLOAT *h, FLOAT nl, INT m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual INT RoughEstimationH(INT k, FLOAT **Y, FLOAT *h, FLOAT nl, INT m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual INT ComponentPdf(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpPdf, INT *Outlier);
    virtual INT LogComponentPdf(INT j, FLOAT **Y, CompnentDistribution *CmpTheta, FLOAT *CmpPdf, INT *Outlier);
    virtual INT EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual INT EnhancedEstimationKDE(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual INT EnhancedEstimationH(INT k, FLOAT **Y, FLOAT nl, FLOAT *h, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual INT MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    virtual INT BayesClassificationKNN(FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual INT BayesClassificationKDE(FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual INT BayesClassificationH(INT k, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual INT DegreesOffreedom(INT c, CompnentDistribution **MixTheta, INT *M);
/// Panic Branislav
    virtual INT EMInitialize();
    virtual INT EMRun(INT *c, FLOAT *W, CompnentDistribution **MixTheta);
/// End 
    INT MixturePdf(INT j, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixPdf);
    INT MixturePdf(FLOAT logV, INT j, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixPdf);
    INT InformationCriterionKNN(INT k, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, INT *M, FLOAT *D);
    INT InformationCriterionKDE(FLOAT logV, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, INT *M, FLOAT *D);
    INT InformationCriterionH(FLOAT logV, INT k, FLOAT **Y, INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, INT *M, FLOAT *D);
    INT InformationCriterion(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, INT *M, FLOAT *D);
    INT CombineComponentsEntropy(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *tau, INT *F, INT *T, FLOAT *EN, FLOAT *ED, FLOAT *PSS);
/// Panic Branislav
    INT CombineComponentsDemp(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *tau, INT *F, INT *T, FLOAT *EN, FLOAT *ED, FLOAT *PSS);
/// End
    INT REBMIX();
    INT Set(char  **Preprocessing,    // Preprocessing type.
            INT   *cmax,              // Maximum number of components.
            INT   *cmin,              // Minimum number of components.
            char  **Criterion,        // Information criterion type.
            INT   *d,                 // Number of independent random variables.
            char  **Variables,        // Types of variables.
            INT   *length_pdf,        // Length of pdf.
            char  **pdf,              // Parametric family types.
            INT   *length_Theta,      // Length of Theta.
            INT   *length_theta,      // Length of Theta[i].
            FLOAT *Theta,             // Component parameters.
            INT   *length_K,          // Length of K.
            INT   *K,                 // Numbers of bins v or numbers of nearest neighbours k.
            INT   *length_ymin,       // Length of ymin.
            FLOAT *ymin,              // Minimum observations.
            INT   *length_ymax,       // Length of ymax.
            FLOAT *ymax,              // Maximum observations.
            INT   *length_h,          // Length of h.
            FLOAT *h,                 // Sides of the hypersquare.
            FLOAT *ar,                // Acceleration rate.
            char  **Restraints,       // Restraints type.
            char  **Mode,             // Mode type.
            INT   *n,                 // Number of observations.
            FLOAT *Y,                 // Dataset.
            INT   *Y_type,            // Dataset type. 
            char  **EMStrategy,       // Strategy for EM algorithm.
            char  **EMVariant,        // EM algorithm variant.
            char  **EMAcceleration,   // Acceleration for the standard EM algorithm.
            FLOAT *EMTolerance,       // Tolerance for EM algortihm.
            FLOAT *EMAccelerationMul, // Acceleration rate for Em algorithm.
            INT   *EMMaxIter,         // Maximum number of iterations in EM algorithm.
            INT   *EMK,               // Number of bins for histogram EM algorithm.
            FLOAT *W,                 // Component weights.
            FLOAT *MixTheta);         // Mixture parameters.

/// Panic Branislav
    INT Get(INT   *n_iter,         // Number of iterations for optimal case.
            INT   *n_iter_sum,     // Number of iterations in whole run.
/// End
            INT   *summary_k,      // Optimal v or optimal k.
            FLOAT *summary_h,      // Optimal class widths of length d.
            FLOAT *summary_y0,     // Optimal origins of length d.
            FLOAT *summary_ymin,   // Optimal minimum observations of length d.
            FLOAT *summary_ymax,   // Optimal maximum observations of length d.
            FLOAT *summary_IC,     // Optimal information criterion.
            FLOAT *summary_logL,   // Log-likelihood.
            INT   *summary_M,      // Degrees of freedom.
            INT   *summary_c,      // Optimal number of components.
            FLOAT *W,              // Component weights.
            FLOAT *theta1,         // Component parameters.
            FLOAT *theta2,         // Component parameters.
            FLOAT *theta3,         // Component parameters.
            INT   *opt_length,     // Length of opt_c, opt_IC, opt_logL and opt_D.
            INT   *opt_c,          // Numbers of components for optimal v or for optimal k.
            FLOAT *opt_IC,         // Information criteria for optimal v or for optimal k.
            FLOAT *opt_logL,       // Log-likelihoods for optimal v or for optimal k.
            FLOAT *opt_Dmin,       // Dmin for optimal v or for optimal k.
            FLOAT *opt_D,          // Totals of positive relative deviations for optimal v or for optimal k.
            INT   *all_length,     // Length of all_K and all_IC.
            INT   *all_K,          // All processed numbers of bins v or all processed numbers of nearest neighbours k.
            FLOAT *all_IC);        // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
    #if (_MAINTAIN_SWITCH)
    INT RunTemplateFile(char *file);
    #endif
}; // Rebmix

#endif
