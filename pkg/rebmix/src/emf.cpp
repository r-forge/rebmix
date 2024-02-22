/*
 *
 * Helper functions for the Expectation-Maximization(EM) algorithm for the finite mixture models.
 *
 * Author: Branislav Panic
 *
*/

#include "rebmixf.h"

#include <math.h>

// Emmix constructor.

Emmix::Emmix()
{
    n_ = 0;
    nr_ = 0;
    nc_ = 0;
    Y_ = NULL;
    cmax_ = 0;
    TOL_ = (FLOAT)0.0;
    am_ = (FLOAT)0.0;
    max_iter_ = 0;
    K_ = 0;
    strategy_ = strategy_none;
    variant_ = varEM;
    accel_ = acc_fixed;
 /// Panic Branislav
    merge_ = merge_none;
 /// End
    n_iter_ = 0;
    c_ = 0;
    W_ = NULL;
    MixTheta_ = NULL;
    dW_ = NULL;
    dMixTheta_ = NULL;
    memset(&summary_, 0, sizeof(SummaryParameterType));
    P_ = NULL;
} // Emmix

// Em destructor.

Emmix::~Emmix()
{
    INT i;

    if (P_) {
        for (i = 0; i < cmax_; i++) {
            if (P_[i]) free(P_[i]);
        }

        free(P_);
    }

    if (dMixTheta_) {
        for (i = 0; i < cmax_; i++) {
            if (dMixTheta_[i]) delete dMixTheta_[i];
        }

        delete[] dMixTheta_;
    }

    if (dW_) free(dW_);

    if (MixTheta_) {
        for (i = 0; i < cmax_; i++) {
            if (MixTheta_[i]) delete MixTheta_[i];
        }

        delete[] MixTheta_;
    }
    
    if (W_) free(W_);

    if (Y_) {
        for (i = 0; i < length_pdf_ + 1; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_);
    }
} // ~Emmix

// Emmix initialize. Returns 0 on success, 1 otherwise.

INT Emmix::Initialize(INT                  n,             // Number of observations.
                      INT                  nr,            // Number of rows.
                      INT                  nc,            // Number of columns.
                      FLOAT                **Y,           // Dataset.
                      INT                  cmax,          // Maximum number of components. 
                      INT                  length_pdf,    // Length of pdf.
                      INT                  length_Theta,  // Length of Theta.
                      INT                  *length_theta, // Length of Theta[i].
                      FLOAT                TOL,           // Tolerance for EM algorithm.
                      FLOAT                am,            // Acceleration multiplier for EM algorithm.
                      INT                  max_iter,      // Maximum number of iterations of EM algorithm.
                      INT                  K,             // Number of bins for histogram EM algorithm.
                      EmStrategyType_e     strategy,      // EM strategy utilization.
                      EmVariantType_e      variant,       // Type of EM variant algorithm.
                      EmAccelerationType_e accel)         // Type of acceleration of standard EM algorithm.
{
    INT i, j, Error = EOK;
    
    n_ = n;
    nr_ = nr;
    nc_ = nc;
    cmax_ = cmax;
    length_pdf_ = length_pdf;
    length_Theta_ = length_Theta;

    length_theta_ = (INT*)malloc(length_Theta_ * sizeof(INT));

    E_CHECK(NULL == length_theta_, EEmmixInitialize);

    for (i = 0; i < length_Theta_; i++) {
        length_theta_[i] = (INT)labs(length_theta[i]);
    }

    Y_ = (FLOAT**)malloc((length_pdf_ + 1) * sizeof(FLOAT*));

    E_CHECK(NULL == Y_, EEmmixInitialize);

    for (i = 0; i < length_pdf_ + 1; i++) {
        Y_[i] = (FLOAT*)malloc(nr_ * sizeof(FLOAT));

        E_CHECK(NULL == Y_[i], EEmmixInitialize);
    }
    
    TOL_ = TOL;
    am_ = am;
    max_iter_ = max_iter;
    K_ = K;

    if (nc_ == length_pdf_) {
        if (K_ > 0) {
            Error = Transform(Y);

            E_CHECK(Error != EOK, Error);
        }
        else {
            for (i = 0; i < n_; i++) {
                for (j = 0; j < length_pdf_; j++) {
                    Y_[j][i] = Y[j][i];
                }

                Y_[length_pdf_][i] = (FLOAT)1.0;
            }
        }
    }
    else
    if (nc_ == length_pdf_ + 1) {
        for (i = 0; i < nr_; i++) {
            for (j = 0; j < nc_; j++) {
                Y_[j][i] = Y[j][i];
            }
        }
    }
    else {
        E_CHECK(1, EEmmixInitialize);
    }

    strategy_ = strategy;
    variant_ = variant;
    accel_ = accel;

    W_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    E_CHECK(NULL == W_, EEmmixInitialize);

    MixTheta_ = new CompnentDistribution* [(unsigned INT)cmax_];

    E_CHECK(NULL == MixTheta_, EEmmixInitialize);

    for (i = 0; i < cmax_; i++) {
        MixTheta_[i] = new CompnentDistribution(this);

        E_CHECK(NULL == MixTheta_[i], EEmmixInitialize);

        Error = MixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        E_CHECK(Error != EOK, Error);
    }
    
    dW_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));
    
    E_CHECK(NULL == dW_, EEmmixInitialize);

    dMixTheta_ = new CompnentDistribution* [(unsigned INT)cmax_];

    E_CHECK(NULL == dMixTheta_, EEmmixInitialize);

    for (i = 0; i < cmax_; i++) {
        dMixTheta_[i] = new CompnentDistribution(this);
        
        E_CHECK(NULL == dMixTheta_[i], EEmmixInitialize);
        
        Error = dMixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);
        
        E_CHECK(Error != EOK, Error);
    }

    P_ = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    E_CHECK(NULL == P_, EEmmixInitialize);

    if (nc_ == length_pdf_) {
        for (i = 0; i < cmax_; i++) {
            P_[i] = (FLOAT*)malloc(n_ * sizeof(FLOAT));

            E_CHECK(NULL == P_[i], EEmmixInitialize);
        }
    }
    else
    if (nc_ == length_pdf_ + 1) {
        for (i = 0; i < cmax_; i++) {
            P_[i] = (FLOAT*)malloc(nr_ * sizeof(FLOAT));

            E_CHECK(NULL == P_[i], EEmmixInitialize);
        }
    }

EEXIT:
    
    E_RETURN(Error);
} // Initialize

INT Emmix::Transform(FLOAT **Y)
{
    FLOAT *h = NULL, *y0 = NULL, *ymax = NULL, *ymin = NULL;
    INT   i, j, l, Error = EOK;

    y0 = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == y0, EEmmixTransform);

    ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == ymin, EEmmixTransform);

    for (i = 0; i < length_pdf_; i++) {
        ymin[i] = Y[i][0];

        for (j = 1; j < n_; j++) {
            if (Y[i][j] < ymin[i]) ymin[i] = Y[i][j];
        }
    }

    ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == ymax, EEmmixTransform);

    for (i = 0; i < length_pdf_; i++) {
        ymax[i] = Y[i][0];

        for (j = 1; j < n_; j++) {
            if (Y[i][j] > ymax[i]) ymax[i] = Y[i][j];
        }
    }

    h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == h, EEmmixTransform);

    for (j = 0; j < length_pdf_; j++) {
        h[j] = (ymax[j] - ymin[j]) / K_;

        y0[j] = ymin[j] + (FLOAT)0.5 * h[j];
    }

    nr_ = 0;

    for (i = 0; i < n_; i++) {
        for (j = 0; j < length_pdf_; j++) {
            l = (INT)floor((Y[j][i] - y0[j]) / h[j] + (FLOAT)0.5);

            Y_[j][nr_] = y0[j] + l * h[j];

            if (Y_[j][nr_] < ymin[j]) {
                Y_[j][nr_] += h[j];
            }
            else
            if (Y_[j][nr_] > ymax[j]) {
                Y_[j][nr_] -= h[j];
            }
        }

        for (j = 0; j < nr_; j++) {
            for (l = 0; l < length_pdf_; l++) if ((FLOAT)fabs(Y_[l][j] - Y_[l][nr_]) > (FLOAT)0.5 * h[l]) goto S0;

            Y_[length_pdf_][j] += (FLOAT)1.0; goto S1;
S0:;
        }

        Y_[length_pdf_][nr_] = (FLOAT)1.0; nr_++;
S1:;
    }

EEXIT:

    if (h) free(h);

    if (ymax) free(ymax);

    if (ymin) free(ymin);

    if (y0) free(y0);

    return Error;
} // Transform

// Returns mixture p.d.f..

INT Emmix::MixtureDist(INT                  j,          // Indey of observation.  
                       FLOAT                **Y,        // Pointer to the input array [y0,...,yd-1,...]
                       INT                  c,          // Number of components.
                       FLOAT                *W,         // Component weights.
                       CompnentDistribution **MixTheta, // Mixture parameters.
                       FLOAT                *MixDist)   // Mixture distribution.
{
    FLOAT CmpDist;
    INT   i, Error = EOK;

    *MixDist = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = LogComponentDist(j, Y, MixTheta[i], &CmpDist);

        E_CHECK(Error != EOK, Error);

        *MixDist += W[i] * (FLOAT)exp(CmpDist);
    }

EEXIT:

    E_RETURN(Error);
} // MixtureDist

// Calculates the log likelihood value of current mixture model parameters. 

INT Emmix::LogLikelihood(INT                  c,          // Number of components.
                         FLOAT                *W,         // Component weights.
                         CompnentDistribution **MixTheta, // Mixture parameters.    
                         FLOAT                *LogL)      // Value of log likelihood.
{
    FLOAT MixDist;
    INT   i, Error = EOK;

    *LogL = (FLOAT)0.0;

    for (i = 0; i < nr_; i++) {
        Error = MixtureDist(i, Y_, c, W, MixTheta, &MixDist);

        E_CHECK(Error != EOK, Error);

        if (MixDist > FLOAT_MIN) {
            *LogL += Y_[length_pdf_][i] * (FLOAT)log(MixDist);
        }
        else {
            *LogL += Y_[length_pdf_][i] * (FLOAT)log(FLOAT_MIN);
        }
    }

EEXIT:

    E_RETURN(Error);
} // LogLikelihood

// Calculates posterior probabilities P for all components in mixture model (expectation step of EM algorith). 

INT Emmix::ExpectationStep()
{
    FLOAT CmpDist, *CmpDistArr = NULL, PostProb;
    INT   i, j, Error = EOK;

    CmpDistArr = (FLOAT*)malloc(c_ * sizeof(FLOAT));

    E_CHECK(NULL == CmpDistArr, EEmmixExpectationStep);

    for (i = 0; i < nr_; i++) {
        PostProb = (FLOAT)0.0;

        for (j = 0; j < c_; j++) {
            Error = LogComponentDist(i, Y_, MixTheta_[j], &CmpDist);

            E_CHECK(Error != EOK, Error);

            CmpDist = (FLOAT)exp(CmpDist);

            CmpDistArr[j] = W_[j] * CmpDist;

            PostProb += CmpDistArr[j];
        }

        for (j = 0; j < c_; j++) {
            P_[j][i] = CmpDistArr[j] / (PostProb + FLOAT_MIN);
        }
    }

EEXIT:

    if (CmpDistArr) free(CmpDistArr);

    E_RETURN(Error);
} // ExpectationStep

// Performs hard clustering (k-means like variant of E-step).

INT Emmix::ConditionalStep()
{
    FLOAT TmpVal;
    INT   i, j, MaxPos, Error = EOK;

    for (i = 0; i < nr_; i++) {
        MaxPos = 0; TmpVal = P_[MaxPos][i]; P_[MaxPos][i] = (FLOAT)0.0;

        for (j = 1; j < c_; j++) {
            if (P_[j][i] > TmpVal) {
                MaxPos = j; TmpVal = P_[MaxPos][i];
            }

            P_[j][i] = (FLOAT)0.0;
        }

        P_[MaxPos][i] = (FLOAT)1.0;
    }

    E_RETURN(Error);
} // ConditionalStep

// Performs golden ration search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment.

INT Emmix::GoldenRatioSearch(FLOAT *am_opt) // Optimal acceleration rate.
{
    CompnentDistribution **MixTheta = NULL;
    FLOAT                arLowerBracket = (FLOAT)1.0;
    FLOAT                arUpperBracket = (FLOAT)1.9;
    FLOAT                arUpdateLower = (FLOAT)0.0;
    FLOAT                arUpdateUpper = (FLOAT)0.0;
    FLOAT                LogLLower = (FLOAT)0.0;
    FLOAT                LogLUpper = (FLOAT)0.0;
    FLOAT                *W = NULL;
    INT                  i, j, Error = EOK;

    W = (FLOAT*)malloc(c_ * sizeof(FLOAT));

    E_CHECK(NULL == W, EEmmixGoldenRatioSearch);

    MixTheta = new CompnentDistribution* [(unsigned INT)c_];

    E_CHECK(NULL == MixTheta, EEmmixGoldenRatioSearch);

    for (i = 0; i < c_; i++) {
        W[i] = W_[i];

        MixTheta[i] = new CompnentDistribution(this);

        E_CHECK(NULL == MixTheta[i], EEmmixGoldenRatioSearch);

        Error = MixTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        E_CHECK(Error != EOK, Error);

        for (j = 0; j < length_pdf_; j++) MixTheta[i]->pdf_[j] = MixTheta_[i]->pdf_[j];

        Error = MixTheta[i]->Memmove(MixTheta_[i]);

        E_CHECK(Error != EOK, Error);
    }

    i = 1; Error = EEmmixGoldenRatioSearch;

    while ((i <= max_iter_) && (Error != EOK)) {
        arUpdateLower = (FLOAT)(arUpperBracket - (arUpperBracket - arLowerBracket) * GoldR);

        arUpdateUpper = (FLOAT)(arLowerBracket + (arUpperBracket - arLowerBracket) * GoldR);

        if ((FLOAT)fabs(arUpdateUpper - arUpdateLower) < TOL_) Error = EOK;

        Error = UpdateMixtureParameters(&c_, W, MixTheta, dW_, dMixTheta_, arUpdateLower); 
        
        E_CHECK(Error != EOK, Error);

        Error = LogLikelihood(c_, W, MixTheta, &LogLLower); 
        
        E_CHECK(Error != EOK, Error);

        for (j = 0; j < c_; j++) {
            W[j] = W_[j];

            Error = MixTheta[j]->Memmove(MixTheta_[j]); 
            
            E_CHECK(Error != EOK, Error);
        }

        Error = UpdateMixtureParameters(&c_, W, MixTheta, dW_, dMixTheta_, arUpdateUpper); 
        
        E_CHECK(Error != EOK, Error);

        Error = LogLikelihood(c_, W, MixTheta, &LogLUpper); 
        
        E_CHECK(Error != EOK, Error);

        for (j = 0; j < c_; j++) {
            W[j] = W_[j];

            Error = MixTheta[j]->Memmove(MixTheta_[j]); 
            
            E_CHECK(Error != EOK, Error);
        }

        if (LogLLower < LogLUpper) {
            arUpperBracket = arUpdateUpper;
        }
        else {
            arLowerBracket = arUpdateLower;
        }

        i++;
    }

    *am_opt = (FLOAT)((arUpperBracket + arLowerBracket) / (FLOAT)2.0);

EEXIT:

    if (MixTheta) {
        for (i = 0; i < c_; i++) {
            if (MixTheta[i]) delete MixTheta[i];
        }

        delete[] MixTheta;
    }

    if (W) free(W);

    E_RETURN(Error);
} // GoldenRatioSearch

// Performs line search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment.

INT Emmix::LineSearch(FLOAT *am_opt) // Return value for optimal acceleration rate.
{
    CompnentDistribution **MixTheta = NULL;
    FLOAT                LogL = (FLOAT)0.0;
    FLOAT                LogLUpdate = (FLOAT)0.0;
    FLOAT                am = (FLOAT)1.0;
    FLOAT                *W = NULL;
    INT                  i, j, Error = EOK;

    W = (FLOAT*)malloc(c_ * sizeof(FLOAT));

    E_CHECK(NULL == W, EEmmixLineSearch);

    MixTheta = new CompnentDistribution* [(unsigned INT)c_];

    E_CHECK(NULL == MixTheta, EEmmixLineSearch);

    for (i = 0; i < c_; i++) {
        W[i] = W_[i];

        MixTheta[i] = new CompnentDistribution(this);

        E_CHECK(NULL == MixTheta[i], EEmmixLineSearch);

        Error = MixTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        E_CHECK(Error != EOK, Error);

        for (j = 0; j < length_pdf_; j++) MixTheta[i]->pdf_[j] = MixTheta_[i]->pdf_[j];

        Error = MixTheta[i]->Memmove(MixTheta_[i]); 
        
        E_CHECK(Error != EOK, Error);
    }

    Error = UpdateMixtureParameters(&c_, W, MixTheta, dW_, dMixTheta_, am); 
    
    E_CHECK(Error != EOK, Error);

    Error = LogLikelihood(c_, W, MixTheta, &LogL); 
    
    E_CHECK(Error != EOK, Error);

    *am_opt = am;

    for (j = 0; j < c_; j++) {
        W[j] = W_[j];

        Error = MixTheta[j]->Memmove(MixTheta_[j]); 
        
        E_CHECK(Error != EOK, Error);
    }

    for (i = 0; i < 9; i++) {
        am += (FLOAT)0.1;

        Error = UpdateMixtureParameters(&c_, W, MixTheta, dW_, dMixTheta_, am); 
        
        E_CHECK(Error != EOK, Error);

        Error = LogLikelihood(c_, W, MixTheta, &LogLUpdate); 
        
        E_CHECK(Error != EOK, Error);

        for (j = 0; j < c_; j++) {
            W[j] = W_[j];

            Error = MixTheta[j]->Memmove(MixTheta_[j]); 
            
            E_CHECK(Error != EOK, Error);
        }

        if (LogLUpdate > LogL) {
            LogL = LogLUpdate; *am_opt = am;
        }
    }

EEXIT:

    if (MixTheta) {
        for (i = 0; i < c_; i++) {
            if (MixTheta[i]) delete MixTheta[i];
        }

        delete[] MixTheta;
    }

    if (W) free(W);

    E_RETURN(Error);
} // LineSearch

// Performs the standard EM algorithm (reference).

INT Emmix::EM()
{
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;
    INT   i, Error = EOK;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    E_CHECK(Error != EOK, Error);

    LogLOld = LogLOld / (FLOAT)n_;

    for (i = 0; i < max_iter_; i++) {
        Error = ExpectationStep();

        E_CHECK(Error != EOK, Error);

        Error = MaximizationStep();

        E_CHECK(Error != EOK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != EOK, Error);

        LogLNew = LogLNew / (FLOAT)n_;

        if ((FLOAT)fabs(LogLNew - LogLOld) <= TOL_) break;

        LogLOld = LogLNew;
    }

    n_iter_ = i;

EEXIT:

    E_RETURN(Error);
} // EM

// Performs the Expectation-Conditional-Maximization algorithm (reference).

INT Emmix::ECM()
{
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;
    INT   i, Error = EOK;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    E_CHECK(Error != EOK, Error);

    for (i = 0; i < max_iter_; i++) {
        Error = ExpectationStep();

        E_CHECK(Error != EOK, Error);

        Error = ConditionalStep();

        E_CHECK(Error != EOK, Error);

        Error = MaximizationStep();

        E_CHECK(Error != EOK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != EOK, Error);

        if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        LogLOld = LogLNew;
    }

    n_iter_ = i;

EEXIT:

    E_RETURN(Error);
} // ECM

// Runs the EM algorithm or its variant.

INT Emmix::Run(INT *c, FLOAT *W, CompnentDistribution **MixTheta)
{
    INT i, Error = EOK;

    c_ = *c;

    for (i = 0; i < c_; i++) {
        W_[i] = W[i];

        Error = MixTheta_[i]->Memmove(MixTheta[i]);

        E_CHECK(Error != EOK, Error);
    }

    switch (variant_) {
    case varEM:
        Error = EM();

        E_CHECK(Error != EOK, Error);

        break;
    case varECM:
        Error = ECM();

        E_CHECK(Error != EOK, Error);
    }

    *c = c_;

    for (i = 0; i < *c; i++) {
        W[i] = W_[i];

        Error = MixTheta[i]->Memmove(MixTheta_[i]);

        E_CHECK(Error != EOK, Error);
    }

EEXIT:

    E_RETURN(Error);
} // Run

// Returns logarithm of component p.d.f..

INT Emmix::LogComponentDist(INT                  j,         // Indey of observation.  
                            FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                            CompnentDistribution *CmpTheta, // Component parameters.
                            FLOAT                *CmpDist)  // Component distribution value.
{
    FLOAT p, Theta, y, ypb;
    INT   i, k, n, Error = EOK;

    *CmpDist = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        switch (CmpTheta->pdf_[i]) {
        case pfNormal:
            y = (Y[i][j] - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

            *CmpDist += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]);

            break;
        case pfTNormal:
            break;
        case pfLognormal:
            if (Y[i][j] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i][j]) - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

                *CmpDist += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpDist = -FLOAT_MAX;
            }

            break;
        case pfWeibull:
            if (Y[i][j] > FLOAT_MIN) {
                ypb = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(Y[i][j] / CmpTheta->Theta_[0][i]));

                *CmpDist += (FLOAT)log(CmpTheta->Theta_[1][i]) + (FLOAT)log(ypb) - ypb - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpDist = -FLOAT_MAX;
            }

            break;
        case pfGamma:
            if (Y[i][j] > FLOAT_MIN) {
                ypb = Y[i][j] / CmpTheta->Theta_[0][i];

                *CmpDist += CmpTheta->Theta_[1][i] * (FLOAT)log(ypb) - ypb - Gammaln(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpDist = -FLOAT_MAX;
            }

            break;
        case pfGumbel:
            ypb = CmpTheta->Theta_[2][i] * (Y[i][j] - CmpTheta->Theta_[0][i]) / CmpTheta->Theta_[1][i];

            *CmpDist += ypb - (FLOAT)exp(ypb) - (FLOAT)log(CmpTheta->Theta_[1][i]);

            break;
        case pfvonMises:
            if ((Y[i][j] < (FLOAT)0.0) || (Y[i][j] > Pi2)) {
                *CmpDist = -FLOAT_MAX;
            }
            else {
                *CmpDist += CmpTheta->Theta_[1][i] * (FLOAT)cos(Y[i][j] - CmpTheta->Theta_[0][i]) - LogPi2 - (FLOAT)log(BesselI0(CmpTheta->Theta_[1][i]));
            }

            break;
        case pfBinomial:
            k = (INT)Y[i][j]; n = (INT)CmpTheta->Theta_[0][i]; p = CmpTheta->Theta_[1][i];

            if (k < 0) {
                *CmpDist = -FLOAT_MAX;
            }
            else
            if (k == 0)
                *CmpDist += n * (FLOAT)log((FLOAT)1.0 - p);
            else
            if (k == n)
                *CmpDist += n * (FLOAT)log(p);
            else
            if (k > n) {
                *CmpDist = -FLOAT_MAX;
            }
            else
               *CmpDist += Gammaln(n + (FLOAT)1.0) - Gammaln(k + (FLOAT)1.0) - Gammaln(n - k + (FLOAT)1.0) +
                           k * (FLOAT)log(p) + (n - k) * (FLOAT)log((FLOAT)1.0 - p);

            break;
        case pfPoisson:
            k = (INT)Y[i][j]; Theta = CmpTheta->Theta_[0][i];

            *CmpDist += k * (FLOAT)log(Theta) - Theta - Gammaln(k + (FLOAT)1.0);

            break;
        case pfDirac:
            if ((FLOAT)fabs(Y[i][j] - CmpTheta->Theta_[0][i]) > FLOAT_MIN) {
                *CmpDist = -FLOAT_MAX;
            }
            else {
                *CmpDist += (FLOAT)0.0;
            }

            break;
        case pfUniform:
            if ((Y[i][j] > CmpTheta->Theta_[1][i]) || (Y[i][j] < CmpTheta->Theta_[0][i])) {
                *CmpDist = -FLOAT_MAX;
            }
            else {
                *CmpDist -= (FLOAT)log(CmpTheta->Theta_[1][i] - CmpTheta->Theta_[0][i]);
            }
        }
    }

    E_RETURN(Error);
} // LogComponentDist

// Updates mixture model parameters with appropriate increment.

INT Emmix::UpdateMixtureParameters(INT                  *c,          // Number of components. 
                                   FLOAT                *W,          // Mixture model weight values.
                                   CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                   FLOAT                *dW,         // Update increment of mixture model weights.
                                   CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                                   FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, j, l, Error = EOK;

    for (l = 0; l < *c; l++) {
        W[l] += am * dW[l];

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

/// Panic Branislav
        if (W[l] < FLOAT_MIN) {
            switch (merge_) {
            case merge_naive:
                for (j = l; j < *c - 1; j++) {
                    dW[j] = dW[j + 1];

                    W[j] = W[j + 1];

                    for (i = 0; i < length_pdf_; i++) {
                        MixTheta[j]->Theta_[0][i] = MixTheta[j + 1]->Theta_[0][i];

                        dMixTheta[j]->Theta_[0][i] = dMixTheta[j + 1]->Theta_[0][i];

                        MixTheta[j]->Theta_[1][i] = MixTheta[j + 1]->Theta_[1][i];

                        dMixTheta[j]->Theta_[1][i] = dMixTheta[j + 1]->Theta_[1][i];

                        MixTheta[j]->Theta_[2][i] = MixTheta[j + 1]->Theta_[2][i];

                        dMixTheta[j]->Theta_[2][i] = dMixTheta[j + 1]->Theta_[2][i];
                    }
                }

                (*c)--; l--; goto S1;
                
                break;
            case merge_none:
                break;
            }
        }
/// End

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta[l]->pdf_[i]) {
            case pfNormal:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfWeibull:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[0][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[0][i] = Eps;
                }

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfGamma:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[0][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[0][i] = Eps;
                }

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfGumbel:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfvonMises:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfBinomial:
                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < (FLOAT)0.0) {
                    MixTheta[l]->Theta_[1][i] = (FLOAT)0.0;
                }
                else
                if (MixTheta[l]->Theta_[1][i] > (FLOAT)1.0) {
                    MixTheta[l]->Theta_[1][i] = (FLOAT)1.0;
                }

                break;
            case pfPoisson:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                if (MixTheta[l]->Theta_[0][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[0][i] = Eps;
                }

                break;
            case pfDirac: case pfUniform:
                break;
            }
        }
S1:;
    }

    E_RETURN(Error);
} // UpdateMixtureParameters

// Maximization step of the EM algoritm.

INT Emmix::MaximizationStep()
{
    FLOAT A[5], am_opt = (FLOAT)1.0, *C = NULL, dC, dM, *M = NULL, T[2], W;
    INT   i, j, k, l, Error = EOK;

    M = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == M, EEmmixMaximizationStep);

    C = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == C, EEmmixMaximizationStep);

    for (l = 0; l < c_; l++) {
        W = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) {
            W += Y_[length_pdf_][j] * P_[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta_[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < nr_; j++) {
                    M[i] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
                }

                M[i] = M[i] / (W + FLOAT_MIN);

                dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                for (j = 0; j < nr_; j++)  if (Y_[i][j] > FLOAT_MIN) {
                    M[i] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)log(Y_[i][j]);
                }

                M[i] = M[i] / (W + FLOAT_MIN);

                dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfWeibull:
                M[i] = MixTheta_[l]->Theta_[1][i];

                j = 1; Error = EEmmixMaximizationStep;
                while ((j <= ItMax) && (Error != EOK)) {
                    memset(&A, 0, 5 * sizeof(FLOAT));

                    for (k = 0; k < nr_; k++) if (Y_[i][k] > FLOAT_MIN) {
                        T[0] = (FLOAT)log(Y_[i][k]);
                        T[1] = (FLOAT)exp(T[0] * M[i]);

                        A[1] += Y_[length_pdf_][k] * P_[l][k] * T[0];
                        A[2] += Y_[length_pdf_][k] * P_[l][k] * T[1];
                        A[3] += Y_[length_pdf_][k] * P_[l][k] * T[1] * T[0];
                        A[4] += Y_[length_pdf_][k] * P_[l][k] * T[1] * T[0] * T[0];
                    }

                    A[0] = W + FLOAT_MIN;

                    T[0] = A[0] / A[2];

                    dM = (A[0] / M[i] + A[1] - T[0] * A[3]) / (T[0] * (A[3] * A[3] / A[2] - A[4]) - A[0] / M[i] / M[i]);

                    M[i] -= dM;

                    E_CHECK(IsNan(dM) || IsInf(dM), EEmmixMaximizationStep);

                    if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = EOK;

                    j++;
                }

                E_CHECK(Error != EOK, Error);

                dMixTheta_[l]->Theta_[1][i] = M[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfGamma:
                M[i] = MixTheta_[l]->Theta_[1][i];

                memset(&A, 0, 4 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) if (Y_[i][j] > FLOAT_MIN) {
                    A[1] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
                    A[2] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)log(Y_[i][j]);
                }

                A[0] = W + FLOAT_MIN;

                A[3] = (FLOAT)log(A[1] / A[0]) - A[2] / A[0];

                j = 1; Error = EEmmixMaximizationStep;
                while ((j <= ItMax) && (Error != EOK)) {
                    if ((Digamma(M[i], &T[0]) != EOK) || (Digamma(M[i] + Eps, &T[1]) != EOK)) goto EEXIT;

                    dM = (T[0] + A[3] - (FLOAT)log(M[i])) / ((T[1] - T[0]) / Eps - (FLOAT)1.0 / M[i]);

                    M[i] -= dM;

                    E_CHECK(IsNan(dM) || IsInf(dM), EEmmixMaximizationStep);

                    if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = EOK;

                    j++;
                }

                E_CHECK(Error != EOK, Error);

                dMixTheta_[l]->Theta_[1][i] = M[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfGumbel:
                M[i] = MixTheta_[l]->Theta_[1][i];

                j = 1; Error = EEmmixMaximizationStep;
                while ((j <= ItMax) && (Error != EOK)) {
                    memset(&A, 0, 5 * sizeof(FLOAT));

                    for (k = 0; k < nr_; k++) {
                        T[0] = (FLOAT)exp(MixTheta_[l]->Theta_[2][i] * Y_[i][k] / M[i]);

                        A[1] += Y_[length_pdf_][k] * P_[l][k] * Y_[i][k];
                        A[2] += Y_[length_pdf_][k] * P_[l][k] * T[0];
                        A[3] += Y_[length_pdf_][k] * P_[l][k] * Y_[i][k] * T[0];
                        A[4] += Y_[length_pdf_][k] * P_[l][k] * Y_[i][k] * Y_[i][k] * T[0];
                    }

                    A[0] = W + FLOAT_MIN;

                    T[0] = A[0] / A[2]; T[1] = A[3] / A[2];

                    dM = (M[i] * A[0] + MixTheta_[l]->Theta_[2][i] * (A[1] - T[0] * A[3])) / (A[0] + T[0] * (A[4] - T[1] * A[3]) / M[i] / M[i]);

                    M[i] -= dM;

                    E_CHECK(IsNan(dM) || IsInf(dM), EEmmixMaximizationStep);

                    if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = EOK;

                    j++;
                }

                E_CHECK(Error != EOK, Error);

                dMixTheta_[l]->Theta_[1][i] = M[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfvonMises:
                memset(&A, 0, 3 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    A[0] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)cos(Y_[i][j]);
                    A[1] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)sin(Y_[i][j]);
                }

                A[0] /= n_; A[1] /= n_; A[2] = (FLOAT)sqrt((FLOAT)pow(A[0], (FLOAT)2.0) + (FLOAT)pow(A[1], (FLOAT)2.0));

                if (A[1] > FLOAT_MIN) {
                    M[i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]);
                }
                else
                if (A[1] < -FLOAT_MIN) {
                    M[i] = (FLOAT)2.0 * (FLOAT)atan((A[2] - A[0]) / A[1]) + Pi2;
                }
                else
                if (A[0] > FLOAT_MIN) {
                    M[i] = (FLOAT)0.0;
                }
                else
                if (A[0] < -FLOAT_MIN) {
                    M[i] = Pi;
                }
                else {
                    E_CHECK(1, EEmmixMaximizationStep);
                }

                dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfBinomial:
                dMixTheta_[l]->Theta_[0][i] = (FLOAT)0.0;

                break;
            case pfPoisson:
                for (j = 0; j < nr_; j++) {
                    M[i] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
                }

                M[i] = M[i] / (W + FLOAT_MIN);

                dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfDirac:
                dMixTheta_[l]->Theta_[0][i] = (FLOAT)0.0;

                break;
            case pfUniform:
                dMixTheta_[l]->Theta_[0][i] = (FLOAT)0.0;
            }
        }

        memset(C, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta_[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < nr_; j++) {
                    C[i] += Y_[length_pdf_][j] * P_[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
                }

                C[i] = (FLOAT)sqrt(C[i] / (W + FLOAT_MIN));

                dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                for (j = 0; j < nr_; j++) if (Y_[i][j] > FLOAT_MIN)  {
                    C[i] += Y_[length_pdf_][j] * P_[l][j] * ((FLOAT)log(Y_[i][j]) - M[i]) * ((FLOAT)log(Y_[i][j]) - M[i]);
                }

                C[i] = (FLOAT)sqrt(C[i] / (W + FLOAT_MIN));

                dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfWeibull:
                memset(&A, 0, 2 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) if (Y_[i][j] > FLOAT_MIN) {
                    T[0] = (FLOAT)log(Y_[i][j]);
                    T[1] = (FLOAT)exp(T[0] * M[i]);

                    A[1] += Y_[length_pdf_][j] * P_[l][j] * T[1];
                }

                A[0] = W + FLOAT_MIN;

                T[0] = A[1] / A[0];

                C[i] = (FLOAT)exp((FLOAT)log(T[0]) / M[i]);

                dMixTheta_[l]->Theta_[0][i] = C[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfGamma:
                memset(&A, 0, 2 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) if (Y_[i][j] > FLOAT_MIN) {
                    A[1] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
                }

                A[0] = W + FLOAT_MIN;

                C[i] = A[1] / A[0] / M[i];

                dMixTheta_[l]->Theta_[0][i] = C[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfGumbel:
                memset(&A, 0, 2 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    T[0] = (FLOAT)exp(MixTheta_[l]->Theta_[2][i] * Y_[i][j] / M[i]);

                    A[1] += Y_[length_pdf_][j] * P_[l][j] * T[0];
                }

                A[0] = W + FLOAT_MIN;

                T[0] = A[1] / A[0];

                C[i] = MixTheta_[l]->Theta_[2][i] * M[i] * (FLOAT)log(T[0]);

                dMixTheta_[l]->Theta_[0][i] = C[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfvonMises:
                memset(&A, 0, 3 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    A[1] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)cos(Y_[i][j] - M[i]);
                }

                A[0] = W + FLOAT_MIN; A[2] = A[1] / A[0];

                C[i] = MixTheta_[l]->Theta_[1][i];

                j = 1; Error = EEmmixMaximizationStep;
                while ((j <= ItMax) && (Error != EOK)) {
                    A[0] = BesselI0(C[i]); A[1] = BesselI1(C[i]);

                    dC = (A[1] - A[2] * A[0]) / (A[0] - (A[2] + (FLOAT)1.0 / C[i]) * A[1]);

                    E_CHECK(IsNan(dC) || IsInf(dC), EEmmixMaximizationStep);

                    C[i] -= dC;

                    if ((FLOAT)fabs(dC) < Max(Eps * (FLOAT)fabs(C[i]), Eps)) Error = EOK;

                    j++;
                }

                dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfBinomial:
                for (j = 0; j < nr_; j++) {
                    C[i] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
                }

                C[i] = C[i] / (W + FLOAT_MIN) / MixTheta_[l]->Theta_[0][i];

                dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfPoisson:
                dMixTheta_[l]->Theta_[1][i] = (FLOAT)0.0;

                break;
            case pfDirac:
                dMixTheta_[l]->Theta_[1][i] = (FLOAT)0.0;

                break;
            case pfUniform:
                dMixTheta_[l]->Theta_[1][i] = (FLOAT)0.0;
            }
        }

        dW_[l] = W / n_ - W_[l];
    }

    if (accel_ == acc_golden) {
        Error = GoldenRatioSearch(&am_opt);

        if (Error != EOK) {
            Error = EOK; am_opt = (FLOAT)1.0;
        }
    }
    else
    if (accel_ == acc_line) {
        Error = LineSearch(&am_opt);

        if (Error != EOK) {
            Error = EOK; am_opt = (FLOAT)1.0;
        }
    }
    else
    if (accel_ == acc_fixed) {
        am_opt = am_;
    }
    else {
        am_opt = (FLOAT)1.0;
    }

    Error = UpdateMixtureParameters(&c_, W_, MixTheta_, dW_, dMixTheta_, am_opt);

    E_CHECK(Error != EOK, Error);

EEXIT:

    if (C) free(C); 
    
    if (M) free(M);

    E_RETURN(Error);
} // MaximizationStep

// Returns logarithm of component p.d.f..

INT Emmvnorm::LogComponentDist(INT                  j,         // Indey of observation.  
                               FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                               CompnentDistribution *CmpTheta, // Component parameters.
                               FLOAT                *CmpDist)  // Component distribution value.
{
    FLOAT y, yi, yk;
    INT   i, k, Error = EOK;

    y = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        yi = Y[i][j] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + i] * yi * yi;

        for (k = i + 1; k < CmpTheta->length_pdf_; k++) {
            yk = Y[k][j] - CmpTheta->Theta_[0][k]; y += CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + k] * yi * yk;
        }
    }
    
    *CmpDist = -y - CmpTheta->length_pdf_ * LogSqrtPi2 - (FLOAT)0.5 * CmpTheta->Theta_[3][0];

    E_RETURN(Error);
} // LogComponentDist

// Updates mixture model parameters with appropriate increment.

INT Emmvnorm::UpdateMixtureParameters(INT                  *c,          // Number of components. 
                                      FLOAT                *W,          // Mixture model weight values.
                                      CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                      FLOAT                *dW,         // Update increment of mixture model weights.
                                      CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                                      FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, ii, j, l, p, q, Error = EOK;
    
    for (l = 0; l < *c; l++) {
        W[l] += am * dW[l];

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

/// Panic Branislav
        if (W[l] < FLOAT_MIN) {
            switch (merge_) {
            case merge_naive:
                for (j = l; j < *c - 1; j++) {
                    dW[j] = dW[j + 1];

                    W[j] = W[j + 1];

                    for (i = 0; i < length_pdf_; i++) {
                        MixTheta[j]->Theta_[0][i] = MixTheta[j + 1]->Theta_[0][i];

                        dMixTheta[j]->Theta_[0][i] = dMixTheta[j + 1]->Theta_[0][i];

                        p = i * length_pdf_ + i;

                        MixTheta[j]->Theta_[1][p] = MixTheta[j + 1]->Theta_[1][p];

                        dMixTheta[j]->Theta_[1][p] = dMixTheta[j + 1]->Theta_[1][p];

                        for (ii = 0; ii < i; ii++) {
                            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                            MixTheta[j]->Theta_[1][p] = MixTheta[j + 1]->Theta_[1][p];

                            dMixTheta[j]->Theta_[1][p] = dMixTheta[j + 1]->Theta_[1][p];

                            MixTheta[j]->Theta_[1][q] = MixTheta[j + 1]->Theta_[1][q];

                            dMixTheta[j]->Theta_[1][q] = dMixTheta[j + 1]->Theta_[1][q];
                        }
                    }
                }
                
                (*c)--; l--; goto S1;
                
                break;
            case merge_none:
                break;
            }
        }
/// End

        for (i = 0; i < length_pdf_; i++) {
            MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

            p = i * length_pdf_ + i;

            MixTheta[l]->Theta_[1][p] += am * dMixTheta[l]->Theta_[1][p];

            if (MixTheta[l]->Theta_[1][p] < Eps) {
                W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][p] = Eps;
            }

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                MixTheta[l]->Theta_[1][p] += am * dMixTheta[l]->Theta_[1][p];

                MixTheta[l]->Theta_[1][q] = MixTheta[l]->Theta_[1][p];
            }
        }

        Error = Cholinvdet(length_pdf_, MixTheta[l]->Theta_[1], MixTheta[l]->Theta_[2], MixTheta[l]->Theta_[3]);

        E_CHECK(Error != EOK, Error);
S1:;
    }

EEXIT:

    E_RETURN(Error);
} // UpdateMixtureParameters

// Maximization step of the EM algoritm.

INT Emmvnorm::MaximizationStep()
{
    FLOAT am_opt = (FLOAT)1.0, *C = NULL, *M = NULL, W;
    INT   i, ii, j, l, p, q, Error = EOK;
     
    M = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == M, EEmmvnormMaximizationStep);
    
    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));
    
    E_CHECK(NULL == C, EEmmvnormMaximizationStep);
    
    for (l = 0; l < c_; l++) {
        W = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) {
            W += Y_[length_pdf_][j] * P_[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            for (j = 0; j < nr_; j++) {
                M[i] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
            }

            M[i] = M[i] / (W + FLOAT_MIN);

            dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];
        }

        memset(C, 0, length_pdf_ * length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            p = i * length_pdf_ + i;

            for (j = 0; j < nr_; j++) {
                C[p] += Y_[length_pdf_][j] * P_[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
            }

            dMixTheta_[l]->Theta_[1][p] = C[p] / (W + FLOAT_MIN) - MixTheta_[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii;

                for (j = 0; j < nr_; j++) {
                    C[p] += Y_[length_pdf_][j] * P_[l][j] * (Y_[i][j] - M[i]) * (Y_[ii][j] - M[ii]);
                }

                dMixTheta_[l]->Theta_[1][p] = C[p] / (W + FLOAT_MIN) - MixTheta_[l]->Theta_[1][p];

                q = ii * length_pdf_ + i;

                dMixTheta_[l]->Theta_[1][q] = dMixTheta_[l]->Theta_[1][p];
            }
        }

        dW_[l] = W / n_ - W_[l];
    }

    if (accel_ == acc_golden) {
        Error = GoldenRatioSearch(&am_opt);

        if (Error != EOK) {
            Error = EOK; am_opt = (FLOAT)1.0;
        }
    } 
    else
    if (accel_ == acc_line) {
        Error = LineSearch(&am_opt);

        if (Error != EOK) {
            Error = EOK; am_opt = (FLOAT)1.0;
        }
    } 
    else
    if (accel_ == acc_fixed) {
        am_opt = am_;
    } 
    else {
        am_opt = (FLOAT)1.0;
    }

    Error = UpdateMixtureParameters(&c_, W_, MixTheta_, dW_, dMixTheta_, am_opt);

    E_CHECK(Error != EOK, Error);

EEXIT:
    
    if (C) free(C); 
    
    if (M) free(M);

    E_RETURN(Error);
} // MaximizationStep