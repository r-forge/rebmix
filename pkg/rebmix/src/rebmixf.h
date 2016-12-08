#ifndef REBMIXF_H_INCLUDED
#define REBMIXF_H_INCLUDED

#ifdef _MSC_VER
#pragma warning(disable: 4514)
#pragma warning(disable: 4820)
#endif

#include <stdlib.h>
#include <string.h>

#include "base.h"

typedef enum {
    poHistogram,        // Histogram approach.
    poParzenWindow,     // Parzen window.
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
    FLOAT ym;  // Mode position.
    FLOAT flm; // Component conditional empirical density.
    FLOAT klm; // Component conditional total number of observations.
} RoughParameterType;

typedef struct summaryparametertype {
    int   c;     // Optimal number of components.
    int   k;     // Optimal v or optimal k.
    FLOAT *y0;   // Optimal origins of length d.
    FLOAT *ymin; // Minimum observations.
    FLOAT *ymax; // Maximum observations.
    FLOAT *h;    // Optimal class widths of length d.
    FLOAT IC;    // Optimal information criterion.
    FLOAT logL;  // Log-likelihood.
    int   M;     // Degrees of freedom.
} SummaryParameterType;

typedef struct additinalparametertype {
    int Bracket; // 1 for bracketing and 0 for golden section.
    int a;       // Golden section constant.
    int b;       // Golden section constant.
    int c;       // Golden section constant.
    int d;       // Golden section constant.
} AdditionalParameterType;

class CompnentDistribution : public Base {
public:
    // Members.
    Base                   *owner_;  // Owner object. 
    ParametricFamilyType_e *pdf_;    // Parametric family types.
    FLOAT                  **Theta_; // Component parameters.
    // Constructor.
    CompnentDistribution(Base *owner);
    // Destructor.
    ~CompnentDistribution();
    // Methods.
    int Realloc(int length_pdf, int length_Theta, int *length_theta);
    int Memmove(CompnentDistribution *CmpTheta);
}; // CompnentDistribution

class Rebmix : public Base {
    // Methods.
    int Golden();
    int GlobalModeKNN(int *m, FLOAT **Y, FLOAT *h, int *I);
    int GlobalModePW(int *m, FLOAT **Y, FLOAT *h, int *I);
    int GlobalModeH(int *m, int k, FLOAT **Y, FLOAT *h, int *I);
    int REBMIXKNN();
    int REBMIXPW();
    int REBMIXH();
    int ReadDataFile();
    int WriteDataFile();
public:
    // Input members.
    FLOAT                      p_value_;       // Probability of obtaining a result equal to or "more extreme" than what was actually observed.
    FLOAT                      min_dist_mul_;  // Minimum distance multiplier.
    FLOAT                      ChiSqr_;        // Critical Chi square value for outlier detection and p = 2.0 * p_value_.
    char                       *curr_;         // Path to the currently open data file.
    int                        o_;             // Number of paths.
    char                       **open_;        // Paths to open data files.
    char                       *save_;         // Path to the save data file.
    PreprocessingType_e        Preprocessing_; // Preprocessing type.
    int                        cmax_;          // Maximum number of components.
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
    int                        *all_K_;        // All processed numbers of bins v or all processed numbers of nearest neighbours k.
    FLOAT                      *all_IC_;       // Information criteria for all processed numbers of bins v or all processed numbers of nearest neighbours k.
    AdditionalParameterType    additional_;    // Additional parameters.
    // Constructor.
    Rebmix();
    // Destructor.
    virtual ~Rebmix();
    // Methods.
    virtual int Initialize();
    int PreprocessingKNN(int k, FLOAT *h, FLOAT **Y);
    int PreprocessingPW(FLOAT *h, FLOAT **Y);
    int PreprocessingH(FLOAT *h, FLOAT *y0, int *k, FLOAT **Y);
    virtual int RoughEstimationKNN(FLOAT **Y, int k, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int RoughEstimationPW(FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int RoughEstimationH(int k, FLOAT **Y, FLOAT *h, FLOAT nl, int m, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int ComponentDist(FLOAT *Y, CompnentDistribution *CmpTheta, FLOAT *CmpDist, int *Outlier);
    virtual int EnhancedEstimationKNN(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int EnhancedEstimationPW(FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int EnhancedEstimationH(int k, FLOAT **Y, FLOAT nl, CompnentDistribution *RigidTheta, CompnentDistribution *LooseTheta);
    virtual int MomentsCalculation(CompnentDistribution *CmpTheta, FLOAT *FirstM, FLOAT *SecondM);
    virtual int BayesClassificationKNN(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int BayesClassificationPW(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int BayesClassificationH(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **FirstM, FLOAT **SecondM);
    virtual int DegreesOffreedom(int c, CompnentDistribution **MixTheta, int *M);
    int MixtureDist(FLOAT *Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *MixDist);
    int InformationCriterionKNN(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int InformationCriterionPW(FLOAT V, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int InformationCriterionH(FLOAT V, int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *IC, FLOAT *logL, int *M, FLOAT *D);
    int CombineComponentsKNN(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, int *F, int *T, FLOAT *EN, FLOAT *ED);
    int CombineComponentsPW(FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, int *F, int *T, FLOAT *EN, FLOAT *ED);
    int CombineComponentsH(int k, FLOAT **Y, int c, FLOAT *W, CompnentDistribution **MixTheta, int *F, int *T, FLOAT *EN, FLOAT *ED);
    int REBMIX();
    int RunTemplateFile(char *file);
}; // Rebmix

#endif