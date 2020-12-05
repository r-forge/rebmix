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
    cmax_ = 0;

    TOL_ = (FLOAT)0.0;
    
    am_ = (FLOAT)0.0;

    max_iter_ = 0;

    strategy_ = strategy_none;

    accel_ = acc_fixed;

    variant_ = varEM;

    n_iter_ = 0;

    dW_ = NULL;
    dMixTheta_ = NULL;
    
    memset(&summary_, 0, sizeof(SummaryParameterType));

    P_ = NULL;
} // Emmix

// Em destructor.

Emmix::~Emmix()
{
    int i;

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
} // ~Emmix

// Emmix initialize. Returns 0 on success, 1 otherwise.

int Emmix::Initialize(int                  n,             // Number of observations.
                      FLOAT                **Y,           // Dataset.
                      int                  cmax,          // Maximum number of components. 
                      int                  length_pdf,    // Length of pdf.
                      int                  length_Theta,  // Length of Theta.
                      int                  *length_theta, // Length of Theta[i].
                      FLOAT                TOL,           // Tolerance for EM algorithm.
                      FLOAT                am,            // Acceleration multiplier for EM algorithm.
                      int                  max_iter,      // Maximum number of iterations of EM algorithm.
                      EmStrategyType_e     strategy,      // EM strategy utilization.
                      EmVariantType_e      variant,       // Type of EM variant algorithm.
                      EmAccelerationType_e accel)         // Type of acceleration of standard EM algorithm.
{
    int i, Error = 0;
    
    n_ = n;

    Y_ = Y;

    cmax_ = cmax;

    length_pdf_ = length_pdf;

    length_Theta_ = length_Theta;

    length_theta_ = (int*)malloc(length_Theta_ * sizeof(int));

    Error = NULL == length_theta_; if (Error) goto E0;

    for (i = 0; i < length_Theta_; i++) {
        length_theta_[i] = (int)labs(length_theta[i]);
    }
    
    TOL_ = TOL;

    am_ = am;

    max_iter_ = max_iter;

    strategy_ = strategy;

    variant_ = variant;

    accel_ = accel;

    W_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    Error = NULL == W_; if (Error) goto E0;

    MixTheta_ = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == MixTheta_; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        MixTheta_[i] = new CompnentDistribution(this);

        Error = NULL == MixTheta_[i]; if (Error) goto E0;

        Error = MixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;
    }
    
    dW_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));
    
    Error = NULL == dW_; if (Error) goto E0;

    dMixTheta_ = new CompnentDistribution* [(unsigned int)cmax_];

    Error = NULL == dMixTheta_; if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        dMixTheta_[i] = new CompnentDistribution(this);
        
        Error = NULL == dMixTheta_[i]; if (Error) goto E0;
        
        Error = dMixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);
        
        if (Error) goto E0;
    }

    P_ = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    Error = P_ == NULL;

    if (Error) goto E0;

    for (i = 0; i < cmax_; i++) {
        P_[i] = (FLOAT*)malloc(n_ * sizeof(FLOAT));

        Error = NULL == P_[i];

        if (Error) goto E0;
    }
    
E0: return Error;
} // Initialize

// Returns mixture p.d.f..

int Emmix::MixtureDist(int                  j,          // Indey of observation.  
                       FLOAT                **Y,        // Pointer to the input array [y0,...,yd-1,...]
                       int                  c,          // Number of components.
                       FLOAT                *W,         // Component weights.
                       CompnentDistribution **MixTheta, // Mixture parameters.
                       FLOAT                *MixDist)   // Mixture distribution.
{
    FLOAT CmpDist;
    int   i;
    int   Error = 0;

    *MixDist = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = LogComponentDist(j, Y, MixTheta[i], &CmpDist);

        if (Error) goto E0;

        *MixDist += W[i] * (FLOAT)exp(CmpDist);
    }

E0: return Error;
} // MixtureDist

// Calculates the log likelihood value of current mixture model parameters. 

int Emmix::LogLikelihood(int                  c,          // Number of components.
                         FLOAT                *W,         // Component weights.
                         CompnentDistribution **MixTheta, // Mixture parameters.    
                         FLOAT                *LogL)      // Value of log likelihood.
{
    FLOAT MixDist;
    int   i, Error = 0;

    *LogL = (FLOAT)0.0;

    for (i = 0; i < n_; i++) {
        Error = MixtureDist(i, Y_, c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *LogL += (FLOAT)log(MixDist);
        }
        else {
            *LogL += (FLOAT)log(FLOAT_MIN);
        }
    }

E0: return Error;
} // LogLikelihood

// Calculates posterior probabilities P for all components in mixture model (expectation step of EM algorith). 

int Emmix::ExpectationStep()
{
    FLOAT CmpDist, PostProb;
    FLOAT *CmpDistArr = NULL;
    int   i, j, Error = 0;

    CmpDistArr = (FLOAT*)malloc(c_ * sizeof(FLOAT));

    Error = CmpDistArr == NULL; if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        PostProb = (FLOAT)0.0;

        for (j = 0; j < c_; j++) {
            Error = LogComponentDist(i, Y_, MixTheta_[j], &CmpDist);

            if (Error) goto E0;

            CmpDist = (FLOAT)exp(CmpDist);

            CmpDistArr[j] = W_[j] * CmpDist;

            PostProb += CmpDistArr[j];
        }

        for (j = 0; j < c_; j++) {
            P_[j][i] = CmpDistArr[j] / (PostProb + FLOAT_MIN);
        }
    }

E0: if (CmpDistArr) free(CmpDistArr);

    return Error;
} // ExpectationStep

// Performs hard clustering (k-means like variant of E-step).

int Emmix::ConditionalStep()
{
    int   i, j, MaxPos, Error = 0;
    FLOAT TmpVal;

    for (i = 0; i < n_; i++) {
        MaxPos = 0; TmpVal = P_[MaxPos][i]; P_[MaxPos][i] = (FLOAT)0.0;

        for (j = 1; j < c_; j++) {
            if (P_[j][i] > TmpVal) {
                MaxPos = j; TmpVal = P_[MaxPos][i];
            }

            P_[j][i] = (FLOAT)0.0;
        }

        P_[MaxPos][i] = (FLOAT)1.0;
    }

    return Error;
} // ConditionalStep

// Performs golden ration search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment.

int Emmix::GoldenRatioSearch(FLOAT *am_opt) // Optimal acceleration rate.
{
    FLOAT                arLowerBracket = (FLOAT)1.0;
    FLOAT                arUpperBracket = (FLOAT)1.9;
    FLOAT                arUpdateLower = (FLOAT)0.0;
    FLOAT                arUpdateUpper = (FLOAT)0.0;
    FLOAT                GolderRation = (FLOAT)((sqrt(5.0) + 1.0) / 2.0);
    FLOAT                LogLLower = (FLOAT)0.0;
    FLOAT                LogLUpper = (FLOAT)0.0;
    int                  i, j, max_iter = 100, Error = 0;
    FLOAT                TOL = (FLOAT)0.01;
    FLOAT                *W = NULL;
    CompnentDistribution **MixTheta = NULL;

    W = (FLOAT*)malloc(c_ * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    MixTheta = new CompnentDistribution* [(unsigned int)c_];

    Error = NULL == MixTheta; if (Error) goto E0;

    for (i = 0; i < c_; i++) {
        W[i] = W_[i];

        MixTheta[i] = new CompnentDistribution(this);

        Error = NULL == MixTheta[i]; if (Error) goto E0;

        Error = MixTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        for (j = 0; j < length_pdf_; j++) MixTheta[i]->pdf_[j] = MixTheta_[i]->pdf_[j];

        Error = MixTheta[i]->Memmove(MixTheta_[i]); if (Error) goto E0;
    }

    for (i = 0; i < max_iter; i++) {
        arUpdateLower = (FLOAT)(arUpperBracket - (arUpperBracket - arLowerBracket) / GolderRation);

        arUpdateUpper = (FLOAT)(arLowerBracket + (arUpperBracket - arLowerBracket) / GolderRation);

        if ((FLOAT)fabs(arUpdateUpper - arUpdateLower) < TOL) break;

        Error = UpdateMixtureParameters(c_, W, MixTheta, dW_, dMixTheta_, arUpdateLower); 
        
        if (Error) goto E0;

        Error = LogLikelihood(c_, W, MixTheta, &LogLLower); 
        
        if (Error) goto E0;

        for (j = 0; j < c_; j++) {
            W[j] = W_[j];

            Error = MixTheta[j]->Memmove(MixTheta_[j]); 
            
            if (Error) goto E0;
        }

        Error = UpdateMixtureParameters(c_, W, MixTheta, dW_, dMixTheta_, arUpdateUpper); 
        
        if (Error) goto E0;

        Error = LogLikelihood(c_, W, MixTheta, &LogLUpper); 
        
        if (Error) goto E0;

        for (j = 0; j < c_; j++) {
            W[j] = W_[j];

            Error = MixTheta[j]->Memmove(MixTheta_[j]); 
            
            if (Error) goto E0;
        }

        if (LogLLower < LogLUpper) {
            arUpperBracket = arUpdateUpper;
        }
        else {
            arLowerBracket = arUpdateLower;
        }
    }

    Error = i == max_iter - 1; if (Error) goto E0;

    *am_opt = (FLOAT)((arUpperBracket + arLowerBracket) / (FLOAT)2.0);

E0: if (MixTheta) {
        for (i = 0; i < c_; i++) {
            if (MixTheta[i]) delete MixTheta[i];
        }

        delete[] MixTheta;
    }

    if (W) free(W);

    return Error;
} // GoldenRatioSearch

// Performs line search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment.

int Emmix::LineSearch(FLOAT *am_opt) // Return value for optimal acceleration rate.
{
    int                  i, j, Error = 0;
    FLOAT                LogL = (FLOAT)0.0;
    FLOAT                LogLUpdate = (FLOAT)0.0;
    FLOAT                am = (FLOAT)1.0;
    FLOAT                *W = NULL;
    CompnentDistribution **MixTheta = NULL;

    W = (FLOAT*)malloc(c_ * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    MixTheta = new CompnentDistribution* [(unsigned int)c_];

    Error = NULL == MixTheta; if (Error) goto E0;

    for (i = 0; i < c_; i++) {
        W[i] = W_[i];

        MixTheta[i] = new CompnentDistribution(this);

        Error = NULL == MixTheta[i]; if (Error) goto E0;

        Error = MixTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        if (Error) goto E0;

        for (j = 0; j < length_pdf_; j++) MixTheta[i]->pdf_[j] = MixTheta_[i]->pdf_[j];

        Error = MixTheta[i]->Memmove(MixTheta_[i]); 
        
        if (Error) goto E0;
    }

    Error = UpdateMixtureParameters(c_, W, MixTheta, dW_, dMixTheta_, am); 
    
    if (Error) goto E0;

    Error = LogLikelihood(c_, W, MixTheta, &LogL); 
    
    if (Error) goto E0;

    *am_opt = am;

    for (j = 0; j < c_; j++) {
        W[j] = W_[j];

        Error = MixTheta[j]->Memmove(MixTheta_[j]); 
        
        if (Error) goto E0;
    }

    for (i = 0; i < 9; i++) {
        am += (FLOAT)0.1;

        Error = UpdateMixtureParameters(c_, W, MixTheta, dW_, dMixTheta_, am); 
        
        if (Error) goto E0;

        Error = LogLikelihood(c_, W, MixTheta, &LogLUpdate); 
        
        if (Error) goto E0;

        for (j = 0; j < c_; j++) {
            W[j] = W_[j];

            Error = MixTheta[j]->Memmove(MixTheta_[j]); 
            
            if (Error) goto E0;
        }

        if (LogLUpdate > LogL) {
            LogL = LogLUpdate; *am_opt = am;
        }
    }

E0: if (MixTheta) {
        for (i = 0; i < c_; i++) {
            if (MixTheta[i]) delete MixTheta[i];
        }

        delete[] MixTheta;
    }

    if (W) free(W);

    return Error;
} // LineSearch

// Performs the standard EM algorithm (reference).

int Emmix::EM()
{
    int   i = 0, Error = 0;
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    if (Error) goto E0;

    LogLOld = LogLOld / (FLOAT)n_;

    for (i = 0; i < max_iter_; i++) {
        Error = ExpectationStep();

        if (Error) goto E0;

        Error = MaximizationStep();

        if (Error) goto E0;

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        if (Error) goto E0;

        LogLNew = LogLNew / (FLOAT)n_;

        if ((FLOAT)fabs(LogLNew - LogLOld) <= TOL_) break;

        LogLOld = LogLNew;
    }

    n_iter_ = i;

E0: return Error;
} // EM

// Performs the Expectation-Conditional-Maximization algorithm (reference).

int Emmix::ECM()
{
    int   i, Error = 0;
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    if (Error) goto E0;

    for (i = 0; i < max_iter_; i++) {
        Error = ExpectationStep();

        if (Error) goto E0;

        Error = ConditionalStep();

        if (Error) goto E0;

        Error = MaximizationStep();

        if (Error) goto E0;

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        if (Error) goto E0;

        if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        LogLOld = LogLNew;
    }

    n_iter_ = i;

E0: return Error;
} // ECM

// Runs the EM algorithm or its variant.

int Emmix::Run(int c, FLOAT *W, CompnentDistribution **MixTheta)
{
    int i, Error = 0;

    c_ = c;

    for (i = 0; i < c_; i++) {
        W_[i] = W[i];

        Error = MixTheta_[i]->Memmove(MixTheta[i]);

        if (Error) goto E0;
    }

    switch (variant_) {
    case varEM:
        Error = EM();

        if (Error) goto E0;

        break;
    case varECM:
        Error = ECM();

        if (Error) goto E0;

        break;
    }

    c = c_;

    for (i = 0; i < c; i++) {
        W[i] = W_[i];

        Error = MixTheta[i]->Memmove(MixTheta_[i]);

        if (Error) goto E0;
    }

E0: return Error;
} // Run

// Returns logarithm of component p.d.f..

int Emmix::LogComponentDist(int                  j,         // Indey of observation.  
                            FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                            CompnentDistribution *CmpTheta, // Component parameters.
                            FLOAT                *CmpDist)  // Component distribution value.
{
    FLOAT y, ypb, p, Theta;
    int   i, k, n;
    int   Error = 0;

    *CmpDist = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        switch (CmpTheta->pdf_[i]) {
        case pfNormal:
            y = (Y[i][j] - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

            *CmpDist += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]);

            break;
        case pfLognormal:
            if (Y[i][j] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i][j]) - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

                *CmpDist += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpDist = -FLOAT_MAX; goto E0;
            }

            break;
        case pfWeibull:
            if (Y[i][j] > FLOAT_MIN) {
                ypb = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(Y[i][j] / CmpTheta->Theta_[0][i]));

                *CmpDist += (FLOAT)log(CmpTheta->Theta_[1][i]) + (FLOAT)log(ypb) - ypb - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpDist = -FLOAT_MAX; goto E0;
            }

            break;
        case pfGamma:
            if (Y[i][j] > FLOAT_MIN) {
                ypb = Y[i][j] / CmpTheta->Theta_[0][i];

                *CmpDist += CmpTheta->Theta_[1][i] * (FLOAT)log(ypb) - ypb - Gammaln(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpDist = -FLOAT_MAX; goto E0;
            }

            break;
        case pfGumbel:
            break;
        case pfvonMises:
            if ((Y[i][j] < (FLOAT)0.0) || (Y[i][j] > Pi2)) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else {
                *CmpDist += CmpTheta->Theta_[1][i] * (FLOAT)cos(Y[i][j] - CmpTheta->Theta_[0][i]) - LogPi2 - (FLOAT)log(BesselI0(CmpTheta->Theta_[1][i]));
            }

            break;
        case pfBinomial:
            k = (int)Y[i][j]; n = (int)CmpTheta->Theta_[0][i]; p = CmpTheta->Theta_[1][i];

            if (k < 0) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else
            if (k == 0)
                *CmpDist += n * (FLOAT)log((FLOAT)1.0 - p);
            else
            if (k == n)
                *CmpDist += n * (FLOAT)log(p);
            else
            if (k > n) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else
                *CmpDist += Gammaln(n + (FLOAT)1.0) - Gammaln(k + (FLOAT)1.0) - Gammaln(n - k + (FLOAT)1.0) +
                            k * (FLOAT)log(p) + (n - k) * (FLOAT)log((FLOAT)1.0 - p);

            break;
        case pfPoisson:
            k = (int)Y[i][j]; Theta = CmpTheta->Theta_[0][i];

            *CmpDist += k * (FLOAT)log(Theta) - Theta - Gammaln(k + (FLOAT)1.0);

            break;
        case pfDirac:
            if ((FLOAT)fabs(Y[i][j] - CmpTheta->Theta_[0][i]) > FLOAT_MIN) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else {
                *CmpDist += (FLOAT)0.0;
            }

            break;
        case pfUniform:
            if ((Y[i][j] > CmpTheta->Theta_[1][i]) || (Y[i][j] < CmpTheta->Theta_[0][i])) {
                *CmpDist = -FLOAT_MAX; goto E0;
            }
            else {
                *CmpDist -= (FLOAT)log(CmpTheta->Theta_[1][i] - CmpTheta->Theta_[0][i]);
            }
        }
    }

E0: return Error;
} // LogComponentDist

// Updates mixture model parameters with appropriate increment.

int Emmix::UpdateMixtureParameters(int                  c,           // Number of components. 
                                   FLOAT                *W,          // Mixture model weight values.
                                   CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                   FLOAT                *dW,         // Update increment of mixture model weights.
                                   CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                                   FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    int i, l, Error = 0;

    for (l = 0; l < c; l++) {
        W[l] += am * dW[l];

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta[l]->pdf_[i]) {
            case pfNormal:
                MixTheta[l]->Theta_[0][i] += am * dMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += am * dMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfLognormal:
                break;
            case pfWeibull:
                break;
            case pfGamma:
                break;
            case pfGumbel:
                break;
            case pfvonMises:
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
                break;
            case pfDirac:
                break;
            case pfUniform:
                break;
            }
        }
    }

    return Error;
} // UpdateMixtureParameters

// Maximization step of the EM algoritm.

int Emmix::MaximizationStep()
{
    FLOAT W, am_opt = (FLOAT)1.0;
    FLOAT *M = NULL, *C = NULL;
    int   i, j, l, Error = 0;

    M = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = M == NULL; if (Error) goto E0;

    C = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = C == NULL; if (Error) goto E0;

    for (l = 0; l < c_; l++) {
        W = (FLOAT)0.0;

        for (j = 0; j < n_; j++) {
            W += P_[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta_[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < n_; j++) {
                    M[i] += P_[l][j] * Y_[i][j];
                }

                M[i] = M[i] / (W + FLOAT_MIN);

                break;
            case pfLognormal:
                break;
            case pfWeibull:
                break;
            case pfGamma:
                break;
            case pfGumbel:
                break;
            case pfvonMises:
                break;
            case pfBinomial:
                break;
            case pfPoisson:
                for (j = 0; j < n_; j++) {
                    M[i] += P_[l][j] * Y_[i][j];
                }

                M[i] = M[i] / (W + FLOAT_MIN);

                break;
            case pfDirac:
                break;
            case pfUniform:
                break;
            }

            dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];
        }

        memset(C, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta_[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < n_; j++) {
                    C[i] += P_[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
                }

                C[i] = (FLOAT)sqrt(C[i] / (W + FLOAT_MIN));

                break;
            case pfLognormal:
                break;
            case pfWeibull:
                break;
            case pfGamma:
                break;
            case pfGumbel:
                break;
            case pfvonMises:
                break;
            case pfBinomial:
                for (j = 0; j < n_; j++) {
                    C[i] += P_[l][j] * Y_[i][j];
                }

                C[i] = C[i] / (W + FLOAT_MIN) / MixTheta_[l]->Theta_[0][i];

                break;
            case pfPoisson:
                break;
            case pfDirac:
                break;
            case pfUniform:
                break;
            }

            dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];
        }

        dW_[l] = W / n_ - W_[l];
    }

    if (accel_ == acc_golden) {
        Error = GoldenRatioSearch(&am_opt);

        if (Error) {
            Error = 0; am_opt = (FLOAT)1.0;
        }
    }
    else
    if (accel_ == acc_line) {
        Error = LineSearch(&am_opt);

        if (Error) {
            Error = 0; am_opt = (FLOAT)1.0;
        }
    }
    else
    if (accel_ == acc_fixed) {
        am_opt = am_;
    }
    else {
        am_opt = (FLOAT)1.0;
    }

    Error = UpdateMixtureParameters(c_, W_, MixTheta_, dW_, dMixTheta_, am_opt);

E0: if (M) free(M);

    if (C) free(C);

    return Error;
} // MaximizationStep

// Returns logarithm of component p.d.f..

int Emmvnorm::LogComponentDist(int                  j,         // Indey of observation.  
                               FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                               CompnentDistribution *CmpTheta, // Component parameters.
                               FLOAT                *CmpDist)  // Component distribution value.
{
    FLOAT y, yi, yk;
    int   i, k, Error = 0;

    y = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        yi = Y[i][j] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + i] * yi * yi;

        for (k = i + 1; k < CmpTheta->length_pdf_; k++) {
            yk = Y[k][j] - CmpTheta->Theta_[0][k]; y += CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + k] * yi * yk;
        }
    }
    
    *CmpDist = -y - CmpTheta->length_pdf_ * LogSqrtPi2 - (FLOAT)0.5 * CmpTheta->Theta_[3][0];

    return Error;
} // LogComponentDist

// Updates mixture model parameters with appropriate increment.

int Emmvnorm::UpdateMixtureParameters(int                  c,           // Number of components. 
                                      FLOAT                *W,          // Mixture model weight values.
                                      CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                      FLOAT                *dW,         // Update increment of mixture model weights.
                                      CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                                      FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    int l, i, ii, p, q, Error = 0;
    
    for (l = 0; l < c; l++) {
        W[l] += am * dW[l];

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

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

        if (Error) goto E0;
    }    

E0: return Error;
} // UpdateMixtureParameters

// Maximization step of the EM algoritm.

int Emmvnorm::MaximizationStep()
{
    FLOAT W, am_opt = (FLOAT)1.0;
    FLOAT *M = NULL, *C = NULL;
    int   i, ii, j, l, p, q, Error = 0;
     
    M = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = M == NULL; if (Error) goto E0;
    
    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));
    
    Error = C == NULL; if (Error) goto E0;
    
    for (l = 0; l < c_; l++) {
        W = (FLOAT)0.0;

        for (j = 0; j < n_; j++) {
            W += P_[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            for (j = 0; j < n_; j++) {
                M[i] += P_[l][j] * Y_[i][j];
            }

            M[i] = M[i] / (W + FLOAT_MIN);

            dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];
        }

        memset(C, 0, length_pdf_ * length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            p = i * length_pdf_ + i;

            for (j = 0; j < n_; j++) {
                C[p] += P_[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
            }

            dMixTheta_[l]->Theta_[1][p] = C[p] / (W + FLOAT_MIN) - MixTheta_[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii;

                for (j = 0; j < n_; j++) {
                    C[p] += P_[l][j] * (Y_[i][j] - M[i]) * (Y_[ii][j] - M[ii]);
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

        if (Error) {
            Error = 0; am_opt = (FLOAT)1.0;
        }
    } 
    else
    if (accel_ == acc_line) {
        Error = LineSearch(&am_opt);

        if (Error) {
            Error = 0; am_opt = (FLOAT)1.0;
        }
    } 
    else
    if (accel_ == acc_fixed) {
        am_opt = am_;
    } 
    else {
        am_opt = (FLOAT)1.0;
    }

    Error = UpdateMixtureParameters(c_, W_, MixTheta_, dW_, dMixTheta_, am_opt); 
    
E0: if (M) free(M);

    if (C) free(C);

    return Error;
} // MaximizationStep