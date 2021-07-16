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

	k_ = 0;

	Y_ = NULL;

    cmax_ = 0;

    TOL_ = (FLOAT)0.0;
    
    am_ = (FLOAT)0.0;

    max_iter_ = 0;

	K_ = 0;

    strategy_ = strategy_none;

	variant_ = varEM;

    accel_ = acc_fixed;

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

	if (Y_) {
		for (i = 0; i < length_pdf_ + 1; i++) {
			if (Y_[i]) free(Y_[i]);
		}

		free(Y_);
	}
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
	                  int                  K,             // Number of bins for histogram EM algorithm.
                      EmStrategyType_e     strategy,      // EM strategy utilization.
                      EmVariantType_e      variant,       // Type of EM variant algorithm.
                      EmAccelerationType_e accel)         // Type of acceleration of standard EM algorithm.
{
    int i, j, Error = 0;
    
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

	Y_ = (FLOAT**)malloc((length_pdf_ + 1) * sizeof(FLOAT*));

	Error = NULL == Y_; if (Error) goto E0;

	for (i = 0; i < length_pdf_ + 1; i++) {
		Y_[i] = (FLOAT*)malloc(n_ * sizeof(FLOAT));

		Error = NULL == Y_[i]; if (Error) goto E0;
	}
    
    TOL_ = TOL;

    am_ = am;

    max_iter_ = max_iter;

	K_ = K;

	if (K_ > 0) {
		Error = Transform(Y);
	}
	else {
		for (i = 0; i < n_; i++) {
			for (j = 0; j < length_pdf_; j++) {
				Y_[j][i] = Y[j][i];
			}

			Y_[length_pdf_][i] = (FLOAT)1.0;
		}

		k_ = n_;
	}

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

int Emmix::Transform(FLOAT **Y)
{
	int i, j, l, Error = 0;

	FLOAT *y0 = NULL, *ymin = NULL, *ymax = NULL, *h = NULL;

	y0 = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

	Error = NULL == y0; if (Error) goto E0;

	ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

	Error = NULL == ymin; if (Error) goto E0;

	for (i = 0; i < length_pdf_; i++) {
		ymin[i] = Y[i][0];

		for (j = 1; j < n_; j++) {
			if (Y[i][j] < ymin[i]) ymin[i] = Y[i][j];
		}
	}

	ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

	Error = NULL == ymax; if (Error) goto E0;

	for (i = 0; i < length_pdf_; i++) {
		ymax[i] = Y[i][0];

		for (j = 1; j < n_; j++) {
			if (Y[i][j] > ymax[i]) ymax[i] = Y[i][j];
		}
	}

	h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

	Error = NULL == h; if (Error) goto E0;

	for (j = 0; j < length_pdf_; j++) {
		h[j] = (ymax[j] - ymin[j]) / K_;

		y0[j] = ymin[j] + (FLOAT)0.5 * h[j];
	}

	k_ = 0;

	for (i = 0; i < n_; i++) {
		for (j = 0; j < length_pdf_; j++) {
			l = (int)floor((Y[j][i] - y0[j]) / h[j] + (FLOAT)0.5);

			Y_[j][k_] = y0[j] + l * h[j];

			if (Y_[j][k_] < ymin[j]) {
				Y_[j][k_] += h[j];
			}
			else
				if (Y_[j][k_] > ymax[j]) {
					Y_[j][k_] -= h[j];
				}
		}

		for (j = 0; j < k_; j++) {
			for (l = 0; l < length_pdf_; l++) if ((FLOAT)fabs(Y_[l][j] - Y_[l][k_]) > (FLOAT)0.5 * h[l]) goto S0;

			Y_[length_pdf_][j] += (FLOAT)1.0; goto S1;
		S0:;
		}

		Y_[length_pdf_][k_] = (FLOAT)1.0; k_++;
	S1:;
	}

E0:	if (h) free(h);

	if (ymax) free(ymax);

	if (ymin) free(ymin);

	if (y0) free(y0);

	return Error;
} // Transform

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

    for (i = 0; i < k_; i++) {
        Error = MixtureDist(i, Y_, c, W, MixTheta, &MixDist);

        if (Error) goto E0;

        if (MixDist > FLOAT_MIN) {
            *LogL += Y_[length_pdf_][i] * (FLOAT)log(MixDist);
        }
        else {
            *LogL += Y_[length_pdf_][i] * (FLOAT)log(FLOAT_MIN);
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

    for (i = 0; i < k_; i++) {
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

    for (i = 0; i < k_; i++) {
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
		case pfTNormal:
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
			y = -(Y[i][j] - CmpTheta->Theta_[0][i]) / CmpTheta->Theta_[1][i];

			*CmpDist += y - (FLOAT)exp(y) - (FLOAT)log(CmpTheta->Theta_[1][i]);

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
	FLOAT dM, dC, A[5], T[2];
    FLOAT W, am_opt = (FLOAT)1.0;
    FLOAT *M = NULL, *C = NULL;
    int   i, j, k, l, Error = 0;

    M = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = M == NULL; if (Error) goto E0;

    C = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    Error = C == NULL; if (Error) goto E0;

    for (l = 0; l < c_; l++) {
        W = (FLOAT)0.0;

        for (j = 0; j < k_; j++) {
            W += Y_[length_pdf_][j] * P_[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta_[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < k_; j++) {
                    M[i] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
                }

                M[i] = M[i] / (W + FLOAT_MIN);

				dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];

                break;
			case pfTNormal:
				break;
            case pfLognormal:
				for (j = 0; j < k_; j++)  if (Y_[i][j] > FLOAT_MIN) {
					M[i] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)log(Y_[i][j]);
				}

				M[i] = M[i] / (W + FLOAT_MIN);

				dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];

				break;
            case pfWeibull:
				M[i] = MixTheta_[l]->Theta_[1][i];

				j = 1; Error = 1;
				while ((j <= ItMax) && Error) {
					memset(&A, 0, 5 * sizeof(FLOAT));

					for (k = 0; k < k_; k++) if (Y_[i][k] > FLOAT_MIN) {
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

					if (IsNan(dM) || IsInf(dM)) {
						Error = 1; goto E0;
					}

					if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = 0;

					j++;
				}

				if (Error) goto E0;

				dMixTheta_[l]->Theta_[1][i] = M[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfGamma:
				M[i] = MixTheta_[l]->Theta_[1][i];

				memset(&A, 0, 4 * sizeof(FLOAT));

				for (j = 0; j < k_; j++) if (Y_[i][j] > FLOAT_MIN) {
					A[1] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
					A[2] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)log(Y_[i][j]);
				}

				A[0] = W + FLOAT_MIN;

				A[3] = (FLOAT)log(A[1] / A[0]) - A[2] / A[0];

				j = 1; Error = 1;
				while ((j <= ItMax) && Error) {
					if (Digamma(M[i], &T[0]) || Digamma(M[i] + Eps, &T[1])) goto E0;

					dM = (T[0] + A[3] - (FLOAT)log(M[i])) / ((T[1] - T[0]) / Eps - (FLOAT)1.0 / M[i]);

					M[i] -= dM;

					if (IsNan(dM) || IsInf(dM)) {
						Error = 1; goto E0;
					}

					if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = 0;

					j++;
				}

				if (Error) goto E0;

				dMixTheta_[l]->Theta_[1][i] = M[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfGumbel:
				M[i] = MixTheta_[l]->Theta_[1][i];

				j = 1; Error = 1;
				while ((j <= ItMax) && Error) {
					memset(&A, 0, 5 * sizeof(FLOAT));

					for (k = 0; k < k_; k++) {
						T[0] = (FLOAT)exp(-Y_[i][k] / M[i]);

						A[1] += Y_[length_pdf_][k] * P_[l][k] * Y_[i][k];
						A[2] += Y_[length_pdf_][k] * P_[l][k] * T[0];
						A[3] += Y_[length_pdf_][k] * P_[l][k] * Y_[i][k] * T[0];
						A[4] += Y_[length_pdf_][k] * P_[l][k] * Y_[i][k] * Y_[i][k] * T[0];
					}

					A[0] = W + FLOAT_MIN;

					T[0] = A[0] / A[2]; T[1] = A[3] / A[2];

					dM = (M[i] * A[0] + T[0] * A[3] - A[1]) / (A[0] + T[0] * (A[4] - T[1] * A[3]) / M[i] / M[i]);

					M[i] -= dM;

					if (IsNan(dM) || IsInf(dM)) {
						Error = 1; goto E0;
					}

					if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = 0;

					j++;
				}

				if (Error) goto E0;

				dMixTheta_[l]->Theta_[1][i] = M[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfvonMises:
				memset(&A, 0, 3 * sizeof(FLOAT));

				for (j = 0; j < k_; j++) {
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
					Error = 1; goto E0;
				}

				dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfBinomial:
				dMixTheta_[l]->Theta_[0][i] = (FLOAT)0.0;

                break;
            case pfPoisson:
                for (j = 0; j < k_; j++) {
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

                break;
            }
        }

        memset(C, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta_[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < k_; j++) {
                    C[i] += Y_[length_pdf_][j] * P_[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
                }

                C[i] = (FLOAT)sqrt(C[i] / (W + FLOAT_MIN));

				dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];

                break;
			case pfTNormal:
				break;
            case pfLognormal:
				for (j = 0; j < k_; j++) if (Y_[i][j] > FLOAT_MIN)  {
					C[i] += Y_[length_pdf_][j] * P_[l][j] * ((FLOAT)log(Y_[i][j]) - M[i]) * ((FLOAT)log(Y_[i][j]) - M[i]);
				}

				C[i] = (FLOAT)sqrt(C[i] / (W + FLOAT_MIN));

				dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];

				break;
            case pfWeibull:
				memset(&A, 0, 2 * sizeof(FLOAT));

				for (j = 0; j < k_; j++) if (Y_[i][j] > FLOAT_MIN) {
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

				for (j = 0; j < k_; j++) if (Y_[i][j] > FLOAT_MIN) {
					A[1] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
				}

				A[0] = W + FLOAT_MIN;

				C[i] = A[1] / A[0] / M[i];

				dMixTheta_[l]->Theta_[0][i] = C[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfGumbel:
				memset(&A, 0, 2 * sizeof(FLOAT));

				for (j = 0; j < k_; j++) {
					T[0] = (FLOAT)exp(-Y_[i][j] / M[i]);

					A[1] += Y_[length_pdf_][j] * P_[l][j] * T[0];
				}

				A[0] = W + FLOAT_MIN;

				T[0] = A[0] / A[1];

				C[i] = M[i] * (FLOAT)log(T[0]);

				dMixTheta_[l]->Theta_[0][i] = C[i] - MixTheta_[l]->Theta_[0][i];

                break;
            case pfvonMises:
				memset(&A, 0, 3 * sizeof(FLOAT));

				for (j = 0; j < k_; j++) {
					A[1] += Y_[length_pdf_][j] * P_[l][j] * (FLOAT)cos(Y_[i][j] - M[i]);
				}

				A[0] = W + FLOAT_MIN; A[2] = A[1] / A[0];

				C[i] = MixTheta_[l]->Theta_[1][i];

				j = 1; Error = 1;
				while ((j <= ItMax) && Error) {
					A[0] = BesselI0(C[i]); A[1] = BesselI1(C[i]);

					dC = (A[1] - A[2] * A[0]) / (A[0] - (A[2] + (FLOAT)1.0 / C[i]) * A[1]);

					if (IsNan(dC) || IsInf(dC)) {
						Error = 1; goto E0;
					}

					C[i] -= dC;

					if ((FLOAT)fabs(dC) < Max(Eps * (FLOAT)fabs(C[i]), Eps)) Error = 0;

					j++;
				}

				dMixTheta_[l]->Theta_[1][i] = C[i] - MixTheta_[l]->Theta_[1][i];

                break;
            case pfBinomial:
                for (j = 0; j < k_; j++) {
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

                break;
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

E0: if (C) free(C); 
	
	if (M) free(M);

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

        for (j = 0; j < k_; j++) {
            W += Y_[length_pdf_][j] * P_[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            for (j = 0; j < k_; j++) {
                M[i] += Y_[length_pdf_][j] * P_[l][j] * Y_[i][j];
            }

            M[i] = M[i] / (W + FLOAT_MIN);

            dMixTheta_[l]->Theta_[0][i] = M[i] - MixTheta_[l]->Theta_[0][i];
        }

        memset(C, 0, length_pdf_ * length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            p = i * length_pdf_ + i;

            for (j = 0; j < k_; j++) {
                C[p] += Y_[length_pdf_][j] * P_[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
            }

            dMixTheta_[l]->Theta_[1][p] = C[p] / (W + FLOAT_MIN) - MixTheta_[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii;

                for (j = 0; j < k_; j++) {
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
    
E0: if (C) free(C); 
	
	if (M) free(M);

    return Error;
} // MaximizationStep