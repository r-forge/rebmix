/*
 *
 * Helper functions for the Expectation-Maximization(EM) algorithm for the finite mixture models.
 *
 * Author: Branislav Panic
 *
*/

#include <math.h>
#include <float.h>

#include <stdio.h>
#include <ctype.h>
#include <time.h>

#include "base.h"
#include "emf.h"

// Em constructor.

Em::Em()
{
    TOL_ = (FLOAT)0.0;
    
    ar_ = (FLOAT)0.0;
    
    C_ = 0;
    
    W_ = NULL;
    
    MixTheta_ = NULL;
    
    dW_ = NULL;
    
    dMixTheta_ = NULL;
    
    memset(&summary_, 0, sizeof(SummaryParameterType));
    
    Initialized_ = 0;
    
    acceleration_t_ = acc_fixed;
    
    variant_t_ = varEM;
    
    n_iter_ = 0;
    
} // Em

// Em destructor.

Em::~Em()
{
    if (Initialized_) {
        int i;
        
        if (W_) free(W_);
        
        if (dW_) free(dW_);
        
        if (MixTheta_) {
            for (i = 0; i < C_; i++) if (MixTheta_[i]) delete MixTheta_[i];
                
            delete[] MixTheta_;
        }

        if (dMixTheta_) {
            for (i = 0; i < C_; i++) if (dMixTheta_[i]) delete dMixTheta_[i];
                
            delete[] dMixTheta_;
        }
    }
} // ~Em

// Em initialize. Returns 0 on success, 1 otherwise.

int Em::Initialize(int C,                              // Number of components. 
                   FLOAT *iniW,                        // Initial mixture model component weights values. 
                   CompnentDistribution **IniMixTheta, // Initial mixture model distribution parameter values.
                   FLOAT TOL,                          // Tolerance threshold for convergence of EM algorithm. 
                   int MAX_ITER,                       // Maximum number of iterations of EM algorithm. 
                   EmVariantType_e algType,            // Variant of EM algorithm (EM, ECM ...).
                   EmAccelerationType_e accelType,     // Acceleration multiplier calculation type.
                   FLOAT accel_mul)                    // Acceleration multiplier value if no calculation if performed.
{
    int Error = 0, i = 0, j = 0;
    
    C_ = C;
    
    TOL_ = TOL;
    
    variant_t_ = algType;
    
    acceleration_t_ = accelType;
    
    ar_ = accel_mul;
    
    MAX_ITER_ = MAX_ITER;
    
    Error = C_ < 1; if (Error) goto E0;
    
    if (C == 1) {
        ar_ = (FLOAT)1.0;
        
        acceleration_t_ = acc_fixed;
        
        variant_t_ = varEM;
    }
    
    W_ = (FLOAT*)malloc(C_ * sizeof(FLOAT));
    
    Error = NULL == W_; if (Error) goto E0;
    
    dW_ = (FLOAT*)malloc(C_ * sizeof(FLOAT));
    
    Error = NULL == dW_; if (Error) goto E0;
    
    MixTheta_ = new CompnentDistribution* [(unsigned int)C_];
    
    Error = NULL == MixTheta_; if (Error) goto E0;
    
    for (i = 0; i < C_; i++) {
        W_[i] = iniW[i];
        
        MixTheta_[i] = new CompnentDistribution(this);
        
        Error = NULL == MixTheta_[i]; if (Error) goto E0;
        
        Error = MixTheta_[i]->Realloc(IniMixTheta[i]->length_pdf_, IniMixTheta[i] -> length_Theta_, IniMixTheta[i] -> length_theta_); 

        if (Error) goto E0;
        
        for (j = 0; j < IniMixTheta[i]->length_pdf_; j++) MixTheta_[i] -> pdf_[j] = IniMixTheta[i] -> pdf_[j];

        Error = MixTheta_[i] -> Memmove(IniMixTheta[i]);
        
        if (Error) goto E0;
    }

    dMixTheta_ = new CompnentDistribution* [(unsigned int)C_];

    Error = NULL == dMixTheta_; if (Error) goto E0;

    for (i = 0; i < C_; i++) {
        dW_[i] = (FLOAT)0.0;
        
        dMixTheta_[i] = new CompnentDistribution(this);
        
        Error = NULL == dMixTheta_[i]; if (Error) goto E0;
        
        Error = dMixTheta_[i]->Realloc(IniMixTheta[i]->length_pdf_, IniMixTheta[i] -> length_Theta_, IniMixTheta[i] -> length_theta_);
        
        if (Error) goto E0;

        for (j = 0; j < IniMixTheta[i]->length_pdf_; j++) dMixTheta_[i] -> pdf_[j] = IniMixTheta[i] -> pdf_[j];
    }

    Initialized_ = 1;
    
E0:
    return Error;
} // Initialize

// Calculates the log component distribution for Gaussian mixture model with unrestricted covariance matrix. Return 0 on success, 1 otherwise.

int Em::LogComponentDist(FLOAT                *Y,        // Pointer to the input point [y0,...,yd-1].
                         CompnentDistribution *CmpTheta, // Component parameters.
                         FLOAT                *CmpDist)  // Component distribution value.
{
    FLOAT y, yi, yj;
    int i, j;
    int Error = 0;

    y = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        yi = Y[i] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + i] * yi * yi;

        for (j = i + 1; j < CmpTheta->length_pdf_; j++) {
            yj = Y[j] - CmpTheta->Theta_[0][j]; y += CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + j] * yi * yj;
        }
    }
    
    *CmpDist = -y - CmpTheta->length_pdf_ * LogSqrtPi2 - (FLOAT)0.5 * CmpTheta->Theta_[3][0];

    return Error;
} // LogComponentDist

// Calculates mixture distribution density probability value. Returns 0 on success, 1 otherwise.

int Em::MixtureDist(FLOAT                *Y,         // Pointer to the input point [y0,...,yd-1].
                    int                  c,          // Number of components.
                    FLOAT                *W,         // Component weights.
                    CompnentDistribution **MixTheta, // Mixture parameters.
                    FLOAT                *MixDist)   // Mixture distribution.
{
    FLOAT CmpDist;
    int  i;
    int  Error = 0;
    
    *MixDist = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = LogComponentDist(Y, MixTheta[i], &CmpDist); if (Error) goto E0;

        *MixDist += W[i] * (FLOAT)exp(CmpDist);
    }
E0: 
    return Error;
} // MixtureDist

// Calculates the posterior probabilities t_ij for each component in mixture model (expectation step of EM algorith). 
// Returns 0 on success, 1 otherwise.

int Em::ExpectationStep(FLOAT                 **Y,        // Pointer to the input points [y0,...,yd-1].
                        int                   N,          // Len of data pointer.
                        int                   c,          // Number of components.
                        FLOAT                 *W,         // Component weights.
                        CompnentDistribution  **MixTheta, // Mixture parameters.
                        FLOAT                  **P)       // Pointer to posterior probabilities.
{
    int i, j, Error = 0;
    FLOAT CmpDist = (FLOAT)0.0, PostProb = (FLOAT)0.0;
    FLOAT *CmpDistArr = NULL;

    CmpDistArr = (FLOAT *)malloc(c * sizeof(FLOAT));
    Error = CmpDistArr == NULL; if (Error) goto E0;

    for (i = 0; i < N; i++) {
        PostProb = (FLOAT)0.0;

        for (j = 0; j < c; j++) {
            Error = LogComponentDist(Y[i], MixTheta[j], &CmpDist);

            if (Error) goto E0;

            CmpDist = (FLOAT)exp(CmpDist);

            CmpDistArr[j] = W[j] * CmpDist;

            PostProb += CmpDistArr[j];
        }

        for (j = 0; j < c; j++) {
            P[i][j] = CmpDistArr[j] / (PostProb + FLOAT_MIN);
        }
    } 
E0:
    if (CmpDistArr != NULL) free(CmpDistArr);

    return Error;
} // ExpectationStep

// Performs hard clustering (k-means like variant of E-step). Returns 0 on success, 1 otherwise.

int Em::ConditionalStep(int   n,    // Length of data pointer.
                        int   c,    // Number of components.
                        FLOAT **P)  // Pointer to posterior probabilities.
{
    int i = 0, j = 0, MaxPos = 0, Error = 0;
    FLOAT TmpVal = (FLOAT)0.0;

    for (i = 0; i < n; i++) {
        TmpVal = P[i][0];

        MaxPos = 0;

        P[i][0] = (FLOAT)0.0;

        for (j = 1; j < c; j++) {
            if (P[i][j] > TmpVal) {
                TmpVal = P[i][j];

                MaxPos = j;
            }

            P[i][j] = (FLOAT)0.0;
        }

        P[i][MaxPos] = (FLOAT)1.0;
    }

    return Error;
} // ConditionalStep

// Updates mixture model parameters with appropriate increment. Return 0 on success, 0 otherwise.

int Em::UpdateMixtureParameters(int                  c,           // Number of components. 
                                FLOAT                *W,          // Mixture model weight values.
                                CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                FLOAT                *dW,         // Update increment of mixture model weights.
                                CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                                FLOAT                ar)          // Acceleration value (multiplicator of incremental update of mixture model parameters).
{
    int i, ii, jj;
    int Error = 0;
    
    for (i = 0; i < c; i++) {
        W[i] += ar * dW[i];

        if (W[i] < (FLOAT)0.0) W[i] = (FLOAT)0.0;

        for (ii = 0; ii < MixTheta[i]->length_pdf_; ii++) {
            MixTheta[i]->Theta_[0][ii] += ar * dMixTheta[i]->Theta_[0][ii];

            MixTheta[i]->Theta_[1][ii * MixTheta[i]->length_pdf_ + ii] += ar * dMixTheta[i]->Theta_[1][ii * MixTheta[i]->length_pdf_ + ii];

            if (MixTheta[i]->Theta_[1][ii * MixTheta[i]->length_pdf_ + ii] < Eps) {
                MixTheta[i]->Theta_[1][ii * MixTheta[i]->length_pdf_ + ii] = Eps;
            }

            for (jj = 0; jj < ii; jj++) {
                MixTheta[i]->Theta_[1][ii * MixTheta[i]->length_pdf_ + jj] += ar * dMixTheta[i]->Theta_[1][ii * MixTheta[i]->length_pdf_ + jj];

                MixTheta[i]->Theta_[1][jj * MixTheta[i]->length_pdf_ + ii] = MixTheta[i]->Theta_[1][ii * MixTheta[i]->length_pdf_ + jj];
            }
        }

        Error = Cholinvdet(MixTheta[i]->length_pdf_, MixTheta[i]->Theta_[1], MixTheta[i]->Theta_[2], MixTheta[i]->Theta_[3]);

        if (Error) return Error;
    }    

    return Error;
} // UpdateMixtureParameters

// Performs line search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment. Returns 0 on success, 1 otherwise.

int Em::LineSearch(FLOAT                **Y,          // Pointer to the input points [y0,...,yd-1]. 
                   int                  N,            // Length of data pointer.
                   int                  c,            // Number of components. 
                   FLOAT                *W,           // Mixture model weight values.
                   CompnentDistribution **MixTheta,   // Mixture model distribution parameter values.
                   FLOAT                *dW,          // Update increment of mixture model weights. 
                   CompnentDistribution **dMixTheta,  // Update increment of mixture model distribution parameter values.
                   FLOAT                *ar_opt)      // Return value for optimal acceleration multiplicator.
{
    int i, j, Error = 0;
    FLOAT LogL = (FLOAT)0.0;
    FLOAT LogLUpdate = (FLOAT)0.0;
    FLOAT ar = (FLOAT)1.0;
    FLOAT *TmpW = NULL;
    CompnentDistribution **TmpMixTheta = NULL;

    TmpW = (FLOAT*)malloc(c * sizeof(FLOAT));

    Error = NULL == TmpW; if (Error) goto E0;

    TmpMixTheta = new CompnentDistribution* [(unsigned int)c];

    Error = NULL == TmpMixTheta; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        TmpW[i] = W[i];

        TmpMixTheta[i] = NULL;

        TmpMixTheta[i] = new CompnentDistribution(this);

        Error = NULL == TmpMixTheta[i]; if (Error) goto E0;

        Error = TmpMixTheta[i]->Realloc(MixTheta[i]->length_pdf_, MixTheta[i] -> length_Theta_, MixTheta[i] -> length_theta_);

        if (Error) goto E0;

        for (j = 0; j < MixTheta[i]->length_pdf_; j++) TmpMixTheta[i] -> pdf_[j] = MixTheta[i] -> pdf_[j];

        Error = TmpMixTheta[i] -> Memmove(MixTheta[i]); if (Error) goto E0;
    }

    Error = UpdateMixtureParameters(c, TmpW, TmpMixTheta, dW, dMixTheta, ar); if (Error) goto E0;

    Error = LogLikelihood(Y, N, c, TmpW, TmpMixTheta, &LogL); if (Error) goto E0;

    *ar_opt = ar;

    for (j = 0; j < c; j++) {
        TmpW[j] = W[j];

        Error = TmpMixTheta[j] -> Memmove(MixTheta[j]); if (Error) goto E0;
    }

    for (i = 0; i < 9; i++) {
        ar += (FLOAT)0.1;

        Error = UpdateMixtureParameters(c, TmpW, TmpMixTheta, dW, dMixTheta, ar); if (Error) goto E0;

        Error = LogLikelihood(Y, N, c, TmpW, TmpMixTheta, &LogLUpdate); if (Error) goto E0;

        for (j = 0; j < c; j++) {
            TmpW[j] = W[j];

            Error = TmpMixTheta[j] -> Memmove(MixTheta[j]); if (Error) goto E0;
        }

        if (LogLUpdate > LogL) {
            LogL = LogLUpdate; *ar_opt = ar;
        }
    }

E0:
    if (TmpW != NULL) {
        free(TmpW);
    }

    if (TmpMixTheta != NULL) {
        for (i = 0; i < c; i++) if (TmpMixTheta[i] != NULL) delete TmpMixTheta[i];
        
        delete[] TmpMixTheta;
    }

    return Error;        
} // LineSearch

// Performs golden ration search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment. Returns 0 on success, 1 otherwise.

int Em::GoldenRatioSearch(FLOAT                **Y,          // Pointer to the input points [y0,...,yd-1].
                          int                  N,            // Length of data pointer.
                          int                  c,            // Number of components.
                          FLOAT                *W,           // Mixture model weight values.
                          CompnentDistribution **MixTheta,   // Mixture model distribution parameter values. 
                          FLOAT                *dW,          // Update increment of mixture model weights. 
                          CompnentDistribution **dMixTheta,  // Update increment of mixture model distribution parameter values. 
                          FLOAT                *ar_opt)      // Return value for optimal acceleration multiplicator.
{
    FLOAT arLowerBracket = (FLOAT)1.0;
    FLOAT arUpperBracket = (FLOAT)1.9;
    FLOAT arUpdateLower = (FLOAT)0.0;
    FLOAT arUpdateUpper = (FLOAT)0.0;
    FLOAT GolderRation = (FLOAT)((sqrt(5.0) + 1.0) / 2.0);
    FLOAT LogLLower = (FLOAT)0.0;
    FLOAT LogLUpper = (FLOAT)0.0;
    int i, j, MAX_ITER = 100, Error = 0;
    FLOAT TOL = (FLOAT)0.01;
    FLOAT *TmpW = NULL;
    CompnentDistribution **TmpMixTheta = NULL;

    TmpW = (FLOAT*)malloc(c * sizeof(FLOAT));

    Error = NULL == TmpW; if (Error) goto E0;

    TmpMixTheta = new CompnentDistribution* [(unsigned int)c];

    Error = NULL == TmpMixTheta; if (Error) goto E0;

    for (i = 0; i < c; i++) {
        TmpW[i] = W[i];

        TmpMixTheta[i] = NULL;

        TmpMixTheta[i] = new CompnentDistribution(this);

        Error = NULL == TmpMixTheta[i]; if (Error) goto E0;

        Error = TmpMixTheta[i]->Realloc(MixTheta[i]->length_pdf_, MixTheta[i] -> length_Theta_, MixTheta[i] -> length_theta_);

        if (Error) goto E0;

        for (j = 0; j < MixTheta[i]->length_pdf_; j++) TmpMixTheta[i] -> pdf_[j] = MixTheta[i] -> pdf_[j];

        Error = TmpMixTheta[i] -> Memmove(MixTheta[i]); if (Error) goto E0; 
    }
    
    for (i = 0; i < MAX_ITER; i++) {
        arUpdateLower = (FLOAT)(arUpperBracket - (arUpperBracket - arLowerBracket) / GolderRation);

        arUpdateUpper = (FLOAT)(arLowerBracket + (arUpperBracket - arLowerBracket) / GolderRation);

        if ((FLOAT)fabs(arUpdateUpper - arUpdateLower) < TOL) {
            break;
        }

        Error = UpdateMixtureParameters(c, TmpW, TmpMixTheta, dW, dMixTheta, arUpdateLower); if (Error) goto E0;

        Error = LogLikelihood(Y, N, c, TmpW, TmpMixTheta, &LogLLower); if (Error) goto E0;

        for (j = 0; j < c; j++) {
            TmpW[j] = W[j];

            Error = TmpMixTheta[j] -> Memmove(MixTheta[j]); if (Error) goto E0;
        }

        Error = UpdateMixtureParameters(c, TmpW, TmpMixTheta, dW, dMixTheta, arUpdateUpper); if (Error) goto E0;
        
        Error = LogLikelihood(Y, N, c, TmpW, TmpMixTheta, &LogLUpper); if (Error) goto E0;

        for (j = 0; j < c; j++) {
            TmpW[j] = W[j];

            Error = TmpMixTheta[j] -> Memmove(MixTheta[j]); if (Error) goto E0;
        }
        
        if (LogLLower < LogLUpper) {
            arUpperBracket = arUpdateUpper;

        } else {
            arLowerBracket = arUpdateLower;
        }
    }

    Error = i == MAX_ITER - 1; if (Error) goto E0;

    *ar_opt = (FLOAT)((arUpperBracket + arLowerBracket) / 2.0);
    
E0:
    if (TmpW != NULL) {
        free(TmpW);
    }

    if (TmpMixTheta != NULL) {
        for (i = 0; i < c; i++) {
            if (TmpMixTheta[i] != NULL) {
                delete TmpMixTheta[i];
            }
        }

        delete[] TmpMixTheta;
    }

    return Error;    
} // GoldenRatioSearch

// Maximization step of EM algoritm. Calculates the update increment of the Gaussian mixture model with unrestricted covariance matrix. Return 0 on success, 1 otherwise.

int Em::MaximizationStep(FLOAT                **Y,        // Pointer to the input points [y0,...,yd-1].
                         int                  N,          // Len of data pointer.
                         int                  c,          // Number of components.
                         FLOAT                *W,         // Component weights.
                         CompnentDistribution **MixTheta, // Mixture parameters.
                         FLOAT                **P)        // Pointer to posterior probabilities.
{
    int i, j, ii, jj, Error = 0;
    FLOAT TmpMeanVal, TmpWeightVal, TmpCovVal;
    FLOAT ar_opt = (FLOAT)1.0;
    FLOAT *TmpMeanVec = NULL;
    FLOAT *TmpCovVec = NULL;
     
    TmpMeanVec = (FLOAT*)malloc(MixTheta[0]->length_pdf_ * sizeof(FLOAT));
    
    TmpCovVec = (FLOAT*)malloc(MixTheta[0]->length_pdf_ * MixTheta[0]->length_pdf_ * sizeof(FLOAT));
    
    Error = TmpMeanVec == NULL; if (Error) goto E0;
    
    for (j = 0; j < c; j++) {
        TmpWeightVal = (FLOAT)0.0; TmpMeanVal = (FLOAT)0.0; TmpCovVal = (FLOAT)0.0;
        
        for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) TmpMeanVec[ii] = (FLOAT)0.0;
        
        for (ii = 0; ii < MixTheta[j]->length_pdf_ * MixTheta[j]->length_pdf_; ii++) TmpCovVec[ii] = (FLOAT)0.0;
        
        for (i = 0; i < N; i++) {
            TmpWeightVal += P[i][j];
            
            for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) {
                TmpMeanVal = P[i][j] * Y[i][ii];
                
                TmpMeanVec[ii] += TmpMeanVal; 
            }
        }
        
        for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) { 
            TmpMeanVec[ii] = TmpMeanVec[ii] / (TmpWeightVal + FLOAT_MIN);
            
            dMixTheta_[j]->Theta_[0][ii] = (FLOAT)(-MixTheta[j]->Theta_[0][ii] + TmpMeanVec[ii]);
        }
        
        for (i = 0; i < N; i++) {
            for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) {
                TmpCovVal = P[i][j] * (Y[i][ii] - TmpMeanVec[ii]) * (Y[i][ii] - TmpMeanVec[ii]);
                
                TmpCovVec[ii * MixTheta[j]->length_pdf_ + ii] += TmpCovVal;
                
                for (jj = 0; jj < ii; jj++) {
                    TmpCovVal = P[i][j] * (Y[i][ii] - TmpMeanVec[ii]) * (Y[i][jj] - TmpMeanVec[jj]);
                    
                    TmpCovVec[ii * MixTheta[j]->length_pdf_ + jj] += TmpCovVal;
                    
                    TmpCovVec[jj * MixTheta[j]->length_pdf_ + ii] = TmpCovVec[ii * MixTheta[j]->length_pdf_ + jj]; 
                }
            }
        }

        for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) {
            TmpCovVec[ii * MixTheta[j]->length_pdf_ + ii] = TmpCovVec[ii * MixTheta[j]->length_pdf_ + ii] / (TmpWeightVal + FLOAT_MIN);
            
            dMixTheta_[j]->Theta_[1][ii * MixTheta[j]->length_pdf_ + ii] = (FLOAT)(-MixTheta[j]->Theta_[1][ii * MixTheta[j]->length_pdf_ + ii] + TmpCovVec[ii * MixTheta[j]->length_pdf_ + ii]);
            
            for (jj = 0; jj < ii; jj++) {
                TmpCovVec[ii * MixTheta[j]->length_pdf_ + jj] = TmpCovVec[ii * MixTheta[j]->length_pdf_ + jj] / (TmpWeightVal + FLOAT_MIN);
                
                TmpCovVec[jj * MixTheta[j]->length_pdf_ + ii] = TmpCovVec[ii * MixTheta[j]->length_pdf_ + jj];

                dMixTheta_[j]->Theta_[1][ii * MixTheta[j]->length_pdf_ + jj] = (FLOAT)(-MixTheta[j]->Theta_[1][ii * MixTheta[j]->length_pdf_ + jj] + TmpCovVec[ii * MixTheta[j]->length_pdf_ + jj]);

                dMixTheta_[j]->Theta_[1][jj * MixTheta[j]->length_pdf_ + ii] = dMixTheta_[j]->Theta_[1][ii * MixTheta[j]->length_pdf_ + jj];
            }
        }

        TmpWeightVal = TmpWeightVal / N;

        dW_[j] = (TmpWeightVal - W_[j]);
    }
    if (acceleration_t_ == acc_golden) {
        Error = GoldenRatioSearch(Y, N, c, W, MixTheta, dW_, dMixTheta_, &ar_opt); if (Error) ar_opt = (FLOAT)1.0;
    } else
    if (acceleration_t_ == acc_line) {
        Error = LineSearch(Y, N, c, W, MixTheta, dW_, dMixTheta_, &ar_opt); if (Error) ar_opt = (FLOAT)1.0;
    } else
    if (acceleration_t_ == acc_fixed) {
        ar_opt = ar_;
    } else {
        ar_opt = (FLOAT)1.0;
    }

    Error = UpdateMixtureParameters(c, W, MixTheta, dW_, dMixTheta_, ar_opt); 
    
E0:
    if (TmpMeanVec != NULL) free(TmpMeanVec);
    if (TmpCovVec != NULL) free(TmpCovVec);

    return Error;
} // MaximizationStep

// Calculates the log likelihood value of current mixture model parameters. Return 0 on success, 1 otherwise.

int Em::LogLikelihood(FLOAT                   **Y,     // Pointer to the input points [y0,...,yd-1].
                      int                   N,         // Len of data pointer.
                      int                  c,          // Number of components.
                      FLOAT                *W,         // Component weights.
                      CompnentDistribution **MixTheta, // Mixture parameters.
                      FLOAT                   *LogL)   // Value of log likelihood.
{
    int i, Error = 0;
    FLOAT MixDistVal = (FLOAT)0.0;
    
    *LogL = (FLOAT)0.0;
    
    for (i = 0; i < N; i++) { 
        Error = MixtureDist(Y[i], c, W, MixTheta, &MixDistVal);
    
        if (Error) goto E0;
        
        if (MixDistVal > FLOAT_MIN) {
            *LogL += (FLOAT)log(MixDistVal);
        }
        else {
            *LogL += (FLOAT)log(FLOAT_MIN); 
        }
    }
E0: return Error;
} // LogLikelihood

// Performs the standard EM algorithm (reference). Return 0 on success, 1 otherwise.

int Em::EM(FLOAT **Y, // Pointer to the input points [y0,...,yd-1].
           int   N)   // Length of data pointer.
{    
    int i = 0, Error = 0;
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;
    FLOAT **P = NULL;
        
    P = (FLOAT**)malloc(N * sizeof(FLOAT*));
        
    Error = P == NULL;
        
    if (Error) goto E0;
        
    for (i = 0; i < N; i++) { 
        P[i] = NULL;
            
        P[i] = (FLOAT*)malloc(C_ * sizeof(FLOAT));
            
        Error = NULL == P[i];
            
        if (Error) goto E0;
    }
        
    Error = LogLikelihood(Y, N, C_, W_, MixTheta_, &LogLOld);
        
    LogLOld = LogLOld / (FLOAT)N;
        
    if (Error) goto E0;
        
    for (i = 0; i < MAX_ITER_; i++) {
        Error = ExpectationStep(Y, N, C_, W_, MixTheta_, P);
            
        if (Error) goto E0;
            
        Error = MaximizationStep(Y, N, C_, W_, MixTheta_, P);
            
        if (Error) goto E0;
            
        Error = LogLikelihood(Y, N, C_, W_, MixTheta_, &LogLNew);
            
        LogLNew = LogLNew / (FLOAT)N;
            
        if ((FLOAT)fabs(LogLNew - LogLOld) <= TOL_) break;
            
        LogLOld = LogLNew;
    }
        
    n_iter_ = i;

E0:
    if (P != NULL) {
        for (i = 0; i < N; i++) { 
             if (P[i] != NULL) {
                free(P[i]);
            }
        }
            
        free(P);
    }
    
    return Error;
} // EM

// Performs the Expectation-Conditional-Maximization algorithm (reference ). Return 0 on success, 1 otherwise.

int Em::ECM(FLOAT  **Y, // Pointer to the input points [y0,...,yd-1].
            int    N)   // Len of data pointer.
{
    int i = 0, Error = 0;
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;
    FLOAT **P = NULL;
        
    P = (FLOAT**)malloc(N * sizeof(FLOAT*));
        
    Error = P == NULL;
        
    if (Error) goto E0;
        
    for (i = 0; i < N; i++) { 
        P[i] = NULL;
            
        P[i] = (FLOAT*)malloc(C_ * sizeof(FLOAT));
            
        Error = NULL == P[i];
            
        if (Error) goto E0;
    }
        
    Error = LogLikelihood(Y, N, C_, W_, MixTheta_, &LogLOld);
        
    if (Error) goto E0;
        
    for (i = 0; i < MAX_ITER_; i++) {
        Error = ExpectationStep(Y, N, C_, W_, MixTheta_, P);
            
        if (Error) goto E0;
            
        Error = ConditionalStep(N, C_, P);
            
        if (Error) goto E0;
            
        Error = MaximizationStep(Y, N, C_, W_, MixTheta_, P);
            
        if (Error) goto E0;
            
        Error = LogLikelihood(Y, N, C_, W_, MixTheta_, &LogLNew);
            
        if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;
            
        LogLOld = LogLNew;
    }
    
    n_iter_ = i;

E0:
    if (P != NULL) { 
        for (i = 0; i < N; i++) {
             if (P[i] != NULL) {
                free(P[i]); 
            }
            
            free(P);
        }
    }

    return Error;
} // ECM

// Runs the EM algorithm or its variant. Return 0 on success, 1 otherwise.

int Em::Run(FLOAT  **Y, // Pointer to the input points [[y0,...,yd-1],...].
            int    N)   // Length of data pointer.
{
    int Error = 0;

    switch (variant_t_) {
    case varEM:
        Error = EM(Y, N);

        if (Error) goto E0;

        break;
    case varECM:
        Error = ECM(Y, N);

        if (Error) goto E0;

        break;
    }

E0: 
    return Error;
} // Run

// Calculates the log component distrubition value for Gaussian mixture model with diagonal covariance structure. Return 0 on success, 1 otherwise.

int GaussianDiagMixture::LogComponentDist(FLOAT                *Y,         // Pointer to the input point [y0,...,yd-1].
                                           CompnentDistribution *CmpTheta, // Component parameters.
                                           FLOAT                *CmpDist)  // Component distribution value.
{
    FLOAT y;
    int i;
    int Error = 0;

    *CmpDist = (FLOAT)0.0;
    
    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        y = (Y[i] - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

        *CmpDist += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]);
    }

    return Error;
} // LogComponentDist

// Updates Gaussian mixture model parameters with diagonal covariance structure with appropriate increment. Return 0 on success, 0 otherwise.

int GaussianDiagMixture::UpdateMixtureParameters(int                  c,           // Number of components. 
                                                 FLOAT                *W,          // Mixture model weight values.
                                                 CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                                 FLOAT                *dW,         // Update increment of mixture model weights.
                                                 CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                                                 FLOAT                ar)          // Acceleration value (multiplicator of incremental update of mixture model parameters).
{
    int i, ii;
    int Error = 0;
    for (i = 0; i < c; i++) {
        W[i] += ar * dW[i];

        for (ii = 0; ii < MixTheta[i]->length_pdf_; ii++) {
            MixTheta[i]->Theta_[0][ii] += ar * dMixTheta[i]->Theta_[0][ii];

            MixTheta[i]->Theta_[1][ii] += ar * dMixTheta[i]->Theta_[1][ii];
        }
    }

    return Error;
} // UpdateMixtureParameters

// Maximization step for Gaussian mixture model with diagonal covariance structure. Return 0 on success, 1 otherwise.

int GaussianDiagMixture::MaximizationStep(FLOAT                   **Y,     // Pointer to the input points [y0,...,yd-1].
                                          int                   N,         // Len of data pointer.
                                          int                  c,          // Number of components.
                                          FLOAT                *W,         // Component weights.
                                          CompnentDistribution **MixTheta, // Mixture parameters.
                                          FLOAT                   **P)     // Pointer to posterior probabilities.
{
    int i, j, ii;
    int Error = 0;
    
    FLOAT TmpMeanVal, TmpWeightVal, TmpCovVal, ar_opt = (FLOAT)1.0;
    
    FLOAT *TmpMeanVec = NULL, *TmpCovVec = NULL;
    
    TmpMeanVec = (FLOAT*)malloc(MixTheta[0]->length_pdf_ * sizeof(FLOAT));
    
    Error = TmpMeanVec == NULL; if (Error) goto E0;
    
    TmpCovVec = (FLOAT*)malloc(MixTheta[0]->length_pdf_ * sizeof(FLOAT));
    
    Error = TmpCovVec == NULL; if (Error) goto E0;
    
    for (j = 0; j < c; j++) {
        TmpWeightVal = (FLOAT)0.0; TmpMeanVal = (FLOAT)0.0; TmpCovVal = (FLOAT)0.0;
        
        for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) {
            TmpMeanVec[ii] = (FLOAT)0.0;
            
            TmpCovVec[ii] = (FLOAT)0.0;
        }
        
        for (i = 0; i < N; i++) {
            TmpWeightVal += P[i][j];
            
            for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) {
                TmpMeanVal = P[i][j] * Y[i][ii];
                
                TmpMeanVec[ii] += TmpMeanVal;
                
                TmpCovVal = P[i][j] * (Y[i][ii] - MixTheta[j]->Theta_[0][ii]) * (Y[i][ii] - MixTheta[j]->Theta_[0][ii]);

                TmpCovVec[ii] += TmpCovVal;
            }
        }
        
        for (ii = 0; ii < MixTheta[j]->length_pdf_; ii++) {
            TmpMeanVec[ii] = TmpMeanVec[ii] / (TmpWeightVal + FLOAT_MIN);

            dMixTheta_[j] -> Theta_[0][ii] = (FLOAT)(-MixTheta[j] -> Theta_[0][ii] + TmpMeanVec[ii]); 
            
            TmpCovVec[ii] = (FLOAT)sqrt(TmpCovVec[ii] / (TmpWeightVal + FLOAT_MIN));

            dMixTheta_[j] -> Theta_[1][ii] = (FLOAT)(-MixTheta[j] -> Theta_[1][ii] + TmpCovVec[ii]);
        }

        TmpWeightVal = TmpWeightVal / N;

        dW_[j] = (TmpWeightVal - W_[j]);
    }

    if (acceleration_t_ == acc_golden) {
        Error = GoldenRatioSearch(Y, N, c, W, MixTheta, dW_, dMixTheta_, &ar_opt); if (Error) ar_opt = (FLOAT)1.0;
    } else
    if (acceleration_t_ == acc_line) {
        Error = LineSearch(Y, N, c, W, MixTheta, dW_, dMixTheta_, &ar_opt); if (Error) ar_opt = (FLOAT)1.0;
    } else
    if (acceleration_t_ == acc_fixed) {
        ar_opt = ar_;
    } else {
        ar_opt = (FLOAT)1.0;
    }

    Error = UpdateMixtureParameters(c, W, MixTheta, dW_, dMixTheta_, ar_opt);

E0:
    if (TmpMeanVec != NULL) free(TmpMeanVec);
    if (TmpCovVec != NULL) free(TmpCovVec);

    return Error;
} // MaximizationStep


