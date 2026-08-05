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
    aam_ = (FLOAT)0.0;
    strategy_ = strategy_none;
    variant_ = varEM;
    accel_ = acc_fixed;
 /// Panic Branislav
    merge_ = merge_none;
 /// End
    toltype_ = likelihood_normalised;
    loglest_ = likelihood_standard;
    n_iter_ = 0;
    c_ = 0;
    W_ = NULL;
    MixTheta_ = NULL;
    dW_ = NULL;
    dMixTheta_ = NULL;
    TW_ = NULL;
    TMixTheta_ = NULL;
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

    if (MixTheta_) {
        for (i = 0; i < cmax_; i++) {
            if (MixTheta_[i]) delete MixTheta_[i];
        }

        delete[] MixTheta_;
    }

    if (TMixTheta_) {
        for (i = 0; i < cmax_; i++) {
            if (TMixTheta_[i]) delete TMixTheta_[i];
        }

        delete [] TMixTheta_;
    }
    
    if (W_) free(W_);

    if (dW_) free(dW_);

    if (TW_) free(TW_);

    if (Y_) {
        for (i = 0; i < length_pdf_ + 1; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_);
    }
} // ~Emmix

// Emmix initialize. Returns 0 on success, 1 otherwise.

INT Emmix::Initialize(INT                  n,                 // Number of observations.
                      INT                  nr,                // Number of rows.
                      INT                  nc,                // Number of columns.
                      FLOAT                **Y,               // Dataset.
                      INT                  cmax,              // Maximum number of components. 
                      INT                  length_pdf,        // Length of pdf.
                      INT                  length_Theta,      // Length of Theta.
                      INT                  *length_theta,     // Length of Theta[i].
                      FLOAT                TOL,               // Tolerance for EM algorithm.
                      FLOAT                am,                // Acceleration multiplier for EM algorithm.
                      INT                  max_iter,          // Maximum number of iterations of EM algorithm.
                      INT                  K,                 // Number of bins for histogram EM algorithm.
                      EmStrategyType_e     strategy,          // EM strategy utilization.
                      EmVariantType_e      variant,           // Type of EM variant algorithm.
                      EmAccelerationType_e accel,             // Type of acceleration of standard EM algorithm.
                      EmAccelParamType_e   accel_eq,          // Equation for calculation of optimal acceleration parameter.
                      EmLikelihoodEstimateType_e loglest,     // Likelihood estimation (normal vs aitken).
                      EmConvergenceType_e  toltype)           // convergence for log likelihood (absolute vs normalised).
{
    INT i, j, Error = E_OK;
    
    n_ = n;
    nr_ = nr;
    nc_ = nc;
    cmax_ = cmax;
    length_pdf_ = length_pdf;
    length_Theta_ = length_Theta;

    length_theta_ = (INT*)malloc(length_Theta_ * sizeof(INT));

    E_CHECK(NULL == length_theta_, E_MEM);

    for (i = 0; i < length_Theta_; i++) {
        length_theta_[i] = (INT)labs(length_theta[i]);
    }

    Y_ = (FLOAT**)malloc((length_pdf_ + 1) * sizeof(FLOAT*));

    E_CHECK(NULL == Y_, E_MEM);

    for (i = 0; i < length_pdf_ + 1; i++) {
        Y_[i] = (FLOAT*)malloc(nr_ * sizeof(FLOAT));

        E_CHECK(NULL == Y_[i], E_MEM);
    }
    
    TOL_ = TOL;
    am_ = am;
    max_iter_ = max_iter;
    K_ = K;

    if (nc_ == length_pdf_) {
        if (K_ > 0) {
            Error = Transform(Y);

            E_CHECK(Error != E_OK, Error);
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
        E_CHECK(1, E_ARG);
    }

    strategy_ = strategy;
    variant_ = variant;
    accel_ = accel;
    accel_eq_ = accel_eq;
    loglest_ = loglest;
    toltype_ = toltype;

    W_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    E_CHECK(NULL == W_, E_MEM);

    MixTheta_ = new CompnentDistribution* [(unsigned INT)cmax_];

    E_CHECK(NULL == MixTheta_, E_MEM);

    for (i = 0; i < cmax_; i++) {
        MixTheta_[i] = new CompnentDistribution(this);

        E_CHECK(NULL == MixTheta_[i], E_MEM);

        Error = MixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        E_CHECK(Error != E_OK, Error);
    }
    
    dW_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));
    
    E_CHECK(NULL == dW_, E_MEM);

    dMixTheta_ = new CompnentDistribution* [(unsigned INT)cmax_];

    E_CHECK(NULL == dMixTheta_, E_MEM);

    for (i = 0; i < cmax_; i++) {
        dMixTheta_[i] = new CompnentDistribution(this);
        
        E_CHECK(NULL == dMixTheta_[i], E_MEM);
        
        Error = dMixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);
        
        E_CHECK(Error != E_OK, Error);
    }

    P_ = (FLOAT**)malloc(cmax_ * sizeof(FLOAT*));

    E_CHECK(NULL == P_, E_MEM);

    TMixTheta_ = new CompnentDistribution* [(unsigned INT)cmax_];

    E_CHECK(NULL == TMixTheta_, E_MEM);

    for (i = 0; i < cmax_; i++) {
        TMixTheta_[i] = new CompnentDistribution(this);
        
        E_CHECK(NULL == TMixTheta_[i], E_MEM);
        
        Error = TMixTheta_[i]->Realloc(length_pdf_, length_Theta_, length_theta_);
        
        E_CHECK(Error != E_OK, Error);
    }


    TW_ = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));
    
    E_CHECK(NULL == TW_, E_MEM);

    if (nc_ == length_pdf_) {
        for (i = 0; i < cmax_; i++) {
            P_[i] = (FLOAT*)malloc(n_ * sizeof(FLOAT));

            E_CHECK(NULL == P_[i], E_MEM);
        }
    }
    else
    if (nc_ == length_pdf_ + 1) {
        for (i = 0; i < cmax_; i++) {
            P_[i] = (FLOAT*)malloc(nr_ * sizeof(FLOAT));

            E_CHECK(NULL == P_[i], E_MEM);
        }
    }

EEXIT:
    
    E_RETURN(Error);
} // Initialize

INT Emmix::Transform(FLOAT **Y)
{
    FLOAT *h = NULL, *y0 = NULL, *ymax = NULL, *ymin = NULL;
    INT   i, j, l, Error = E_OK;

    y0 = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == y0, E_MEM);

    ymin = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == ymin, E_MEM);

    for (i = 0; i < length_pdf_; i++) {
        ymin[i] = Y[i][0];

        for (j = 1; j < n_; j++) {
            if (Y[i][j] < ymin[i]) ymin[i] = Y[i][j];
        }
    }

    ymax = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == ymax, E_MEM);

    for (i = 0; i < length_pdf_; i++) {
        ymax[i] = Y[i][0];

        for (j = 1; j < n_; j++) {
            if (Y[i][j] > ymax[i]) ymax[i] = Y[i][j];
        }
    }

    h = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == h, E_MEM);

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

    E_RETURN(Error);
} // Transform

// Returns mixture p.d.f..

INT Emmix::MixturePdf(INT                  j,           // Indey of observation.  
                       FLOAT                **Y,        // Pointer to the input array [y0,...,yd-1,...]
                       INT                  c,          // Number of components.
                       FLOAT                *W,         // Component weights.
                       CompnentDistribution **MixTheta, // Mixture parameters.
                       FLOAT                *MixPdf)    // Mixture probability density.
{
    FLOAT CmpPdf;
    INT   i, Error = E_OK;

    *MixPdf = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = LogComponentPdf(j, Y, MixTheta[i], &CmpPdf);

        E_CHECK(Error != E_OK, Error);

        *MixPdf += W[i] * (FLOAT)exp(CmpPdf);
    }

EEXIT:

    E_RETURN(Error);
} // MixturePdf

// Calculates the log likelihood value of current mixture model parameters. 

INT Emmix::LogLikelihood(INT                  c,          // Number of components.
                         FLOAT                *W,         // Component weights.
                         CompnentDistribution **MixTheta, // Mixture parameters.    
                         FLOAT                *LogL)      // Value of log likelihood.
{
    FLOAT MixPdf;
    INT   i, Error = E_OK;

    *LogL = (FLOAT)0.0;

    for (i = 0; i < nr_; i++) {
        Error = MixturePdf(i, Y_, c, W, MixTheta, &MixPdf);

        E_CHECK(Error != E_OK, Error);

        if (MixPdf > FLOAT_MIN) {
            *LogL += Y_[length_pdf_][i] * (FLOAT)log(MixPdf);
        }
        else {
            *LogL += Y_[length_pdf_][i] * (FLOAT)log(FLOAT_MIN);
        }
    }

EEXIT:

    E_RETURN(Error);
} // LogLikelihood

// Calculates posterior probabilities P for all components in mixture model (expectation step of EM algorithm). 

INT Emmix::ExpectationStep(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P)
{
    FLOAT CmpPdf, *CmpPdfArr = NULL, PostProb;
    INT   i, j, Error = E_OK;

    CmpPdfArr = (FLOAT*)malloc(c_ * sizeof(FLOAT));

    E_CHECK(NULL == CmpPdfArr, E_MEM);

    for (i = 0; i < nr_; i++) {
        PostProb = (FLOAT)0.0;

        for (j = 0; j < c; j++) {
            Error = LogComponentPdf(i, Y_, MixTheta[j], &CmpPdf);

            E_CHECK(Error != E_OK, Error);

            CmpPdf = (FLOAT)exp(CmpPdf);

            CmpPdfArr[j] = W[j] * CmpPdf;

            PostProb += CmpPdfArr[j];
        }

        for (j = 0; j < c_; j++) {
            P[j][i] = CmpPdfArr[j] / (PostProb + FLOAT_MIN); // can be replaced with multiplication so there is only n divisions.
        }
    }

EEXIT:

    if (CmpPdfArr) free(CmpPdfArr);

    E_RETURN(Error);
} // ExpectationStep

// Performs hard clustering (k-means like variant of E-step).

INT Emmix::ConditionalStep(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P)
{
    FLOAT TmpVal, CmpPdf;
    INT   i, j, MaxPos, Error = E_OK;

    for (i = 0; i < nr_; i++) {
        MaxPos = -1; TmpVal = FLOAT_MIN; 

        for (j = 0; j < c; j++) {
            Error = LogComponentPdf(i, Y_, MixTheta[j], &CmpPdf);

            E_CHECK(Error != E_OK, Error);

            CmpPdf = (FLOAT)exp(CmpPdf);

            CmpPdf *= W[j];

            if (CmpPdf > TmpVal) {
                MaxPos = j; TmpVal = CmpPdf;
            }

            P[j][i] = (FLOAT)0.0;
        }

        if (MaxPos < 0) {
             Error = E_NO_SOLUTION;

             goto EEXIT;
        } 

        P[MaxPos][i] = (FLOAT)1.0;
    }
EEXIT:
    E_RETURN(Error);
} // ConditionalStep

// Draws an partition described by n categorical distributions Pj1...Pjd.

INT Emmix::StochasticStep(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT **P)
{
    FLOAT CmpPdf, *CmpPdfArr = NULL;
    INT   i, j, MaxPos, Seed, Error = E_OK;

    CmpPdfArr = (FLOAT*)malloc((c_ + 1) * sizeof(FLOAT));

    memset(CmpPdfArr, 0, (c_ + 1) * sizeof(FLOAT));

    E_CHECK(NULL == CmpPdfArr, E_MEM);

    for (i = 0; i < nr_; i++) {

        for (j = 0; j < c; j++) {
            Error = LogComponentPdf(i, Y_, MixTheta[j], &CmpPdf);

            E_CHECK(Error != E_OK, Error);

            CmpPdf = (FLOAT)exp(CmpPdf);

            CmpPdf *= W[j];

            CmpPdfArr[j + 1] = CmpPdfArr[j] + CmpPdf;

            P[j][i] = (FLOAT)0.0;
        }

        Seed = -i;

        MaxPos = Choice(&Seed, CmpPdfArr, c);

        P[MaxPos][i] = (FLOAT)1.0;
    }
EEXIT:

    if (CmpPdfArr) free(CmpPdfArr);

    E_RETURN(Error);
} // StochasticStep

// Performs golden ration search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment.

INT Emmix::GoldenRatioSearch(INT *c, 
                             FLOAT *W, 
                             CompnentDistribution **MixTheta, 
                             FLOAT *dW, 
                             CompnentDistribution **dMixTheta, 
                             FLOAT *TW, 
                             CompnentDistribution **TMixTheta, 
                             FLOAT *am_opt) // Optimal acceleration rate.
{
    FLOAT                arLowerBracket = MIN_SEARCH_RANGE, arUpdateLower = (FLOAT)0.0, arUpdateUpper = (FLOAT)0.0, arUpperBracket = MAX_SEARCH_RANGE;
    FLOAT                LogLLower = (FLOAT)0.0, LogLUpper = (FLOAT)0.0; 
    INT                  i, j, error, Error = E_OK;

    for (i = 0; i < *c; i++) {
        TW[i] = W[i];

        Error = TMixTheta[i]->Memmove(MixTheta[i]);

        E_CHECK(Error != E_OK, Error);
    }

    i = 1; Error = E_CON;

    while ((i <= ITER_NBR_GOLDEN_SEARCH) && (Error != E_OK)) {
        arUpdateLower = (FLOAT)(arUpperBracket - (arUpperBracket - arLowerBracket) * GoldR);

        arUpdateUpper = (FLOAT)(arLowerBracket + (arUpperBracket - arLowerBracket) * GoldR);

        if ((FLOAT)fabs(arUpdateUpper - arUpdateLower) < INC_LINE_SEARCH) Error = E_OK;
        
        if (((FLOAT)fabs(LogLUpper - LogLLower) < TOL_ && i > 1) || (FLOAT)fabs(arUpdateUpper - arUpdateLower) < GOLDEN_SEARCH_TOL)  Error = E_OK;

        error = UpdateMixParams(c, TW, TMixTheta, dW, dMixTheta, arUpdateLower); 
        
        E_CHECK(error != E_OK, error);

        error = LogLikelihood(*c, TW, TMixTheta, &LogLLower); 
        
        E_CHECK(error != E_OK, error);

        for (j = 0; j < c_; j++) {
            TW[j] = W[j];

            error = TMixTheta[j]->Memmove(MixTheta[j]); 
            
            E_CHECK(error != E_OK, error);
        }

        error = UpdateMixParams(c, TW, TMixTheta, dW, dMixTheta, arUpdateUpper); 
        
        E_CHECK(error != E_OK, error);

        error = LogLikelihood(*c, TW, TMixTheta, &LogLUpper); 
        
        E_CHECK(error != E_OK, error);

        for (j = 0; j < c_; j++) {
            TW[j] = W[j];

            error = TMixTheta[j]->Memmove(MixTheta[j]); 
            
            E_CHECK(error != E_OK, error);
        }
        // Here seems to be error: it should be > instead of <
        if (LogLLower > LogLUpper) {
            arUpperBracket = arUpdateUpper;
        }
        else {
            arLowerBracket = arUpdateLower;
        }

        i++;
    }

    *am_opt = (FLOAT)((arUpperBracket + arLowerBracket) / (FLOAT)2.0);

EEXIT:

    E_RETURN(Error);
} // GoldenRatioSearch

// Performs line search from minimum value (hardcoded to 1) to maximum value (hardcoded to 1.9) for acceleration constant of mixture parameter update increment.

INT Emmix::LineSearch(INT *c, 
                      FLOAT *W, 
                      CompnentDistribution **MixTheta, 
                      FLOAT *dW, 
                      CompnentDistribution **dMixTheta, 
                      FLOAT *TW, 
                      CompnentDistribution **TMixTheta, 
                      FLOAT *am_opt) // Return value for optimal acceleration rate.
{

    FLOAT am = (FLOAT)1.0, LogL = (FLOAT)0.0, LogLUpdate = (FLOAT)0.0;//, *W = NULL;
    INT   i, j, Error = E_OK;

    for (i = 0; i < c_; i++) {
        TW[i] = W[i];

        Error = TMixTheta[i]->Memmove(MixTheta[i]); 
        
        E_CHECK(Error != E_OK, Error);
    }

    Error = UpdateMixParams(c, TW, TMixTheta, dW, dMixTheta, am); 
    
    E_CHECK(Error != E_OK, Error);

    Error = LogLikelihood(*c, TW, TMixTheta, &LogL); 
    
    E_CHECK(Error != E_OK, Error);

    *am_opt = am;

    for (j = 0; j < c_; j++) {
        TW[j] = W[j];

        Error = TMixTheta[j]->Memmove(MixTheta[j]); 
        
        E_CHECK(Error != E_OK, Error);
    }
    
    for (i = 0; i < ITER_NBR_LINE_SEARCH; i++) {
        am += INC_LINE_SEARCH; 

        Error = UpdateMixParams(c, TW, TMixTheta, dW, dMixTheta, am); 
        
        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(*c, TW, TMixTheta, &LogLUpdate); 
        
        E_CHECK(Error != E_OK, Error);

        for (j = 0; j < c_; j++) {
            TW[j] = W[j];

            Error = TMixTheta[j]->Memmove(MixTheta[j]); 
            
            E_CHECK(Error != E_OK, Error);
        }

        if (LogLUpdate > LogL) {
            LogL = LogLUpdate; *am_opt = am;
        }
    }

EEXIT:

    E_RETURN(Error);
} // LineSearch

INT Emmix::Resids(INT *c, 
                  FLOAT *W, 
                  CompnentDistribution **MixTheta,
                  FLOAT *RW, 
                  CompnentDistribution **RMixTheta, 
                  FLOAT *VW, 
                  CompnentDistribution **VMixTheta,
                  FLOAT **P) 
{
    INT i, l, Error = E_OK;
    
    Error = F(*c, W, MixTheta, RW, RMixTheta, P);

    E_CHECK(Error != E_OK, Error);

    Error = UpdateMixParams(c, W, MixTheta, RW, RMixTheta, (FLOAT)1.0);

    E_CHECK(Error != E_OK, Error);

    Error = ExpectationStep(*c, W, MixTheta, P);

    E_CHECK(Error != E_OK, Error);

    Error = F(*c, W, MixTheta, VW, VMixTheta, P);

    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < *c; l++) {
        VW[l] -= RW[l];
        for (i = 0; i < length_pdf_; i++) {
            switch(VMixTheta[l]->pdf_[i]) {
            case pfNormal:
                VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];

                VMixTheta[l]->Theta_[1][i] -= RMixTheta[l]->Theta_[1][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];
                
                VMixTheta[l]->Theta_[1][i] -= RMixTheta[l]->Theta_[1][i];

                break;
            case pfWeibull:
                VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];
                
                VMixTheta[l]->Theta_[1][i] -= RMixTheta[l]->Theta_[1][i];

                break;
            case pfGamma:
                VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];
                
                VMixTheta[l]->Theta_[1][i] -= RMixTheta[l]->Theta_[1][i];

                break;
            case pfGumbel:
                VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];
                
                VMixTheta[l]->Theta_[1][i] -= RMixTheta[l]->Theta_[1][i];

                break;
            case pfvonMises:
                VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];
                
                VMixTheta[l]->Theta_[1][i] -= RMixTheta[l]->Theta_[1][i];

                break;
            case pfBinomial:
                VMixTheta[l]->Theta_[1][i] -= RMixTheta[l]->Theta_[1][i];

                break;
            case pfPoisson:
                VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];

                break;
            case pfDirac: case pfUniform:
                break;
            }
        }
    }

EEXIT:

    E_RETURN(Error);
}

INT Emmix::OptSLen1(INT *c, 
                    FLOAT *W, 
                    CompnentDistribution **MixTheta,
                    FLOAT *RW, 
                    CompnentDistribution **RMixTheta, 
                    FLOAT *VW, 
                    CompnentDistribution **VMixTheta,
                    FLOAT **P,
                    FLOAT *am_opt) 
{
    INT i, l, Error = E_OK;

    FLOAT A = (FLOAT)0.0, B = (FLOAT)0.0;

    Error = Resids(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P);
    
    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < c_; l++) {
        A += RW[l] * VW[l];
        B += VW[l] * VW[l];
        for (i = 0; i < length_pdf_; i++) {
            switch(VMixTheta[l]->pdf_[i]) {
            case pfNormal:
                A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfWeibull:
                A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfGamma:
                A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfGumbel:
                A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfvonMises:
                A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfBinomial:
                A += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfPoisson:
                A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                break;
            case pfDirac: case pfUniform:
                break;
            }
        }
    }

    *am_opt = Max(A / (B + Eps), (FLOAT)1.0);

EEXIT:

    E_RETURN(Error);
}
INT Emmix::OptSLen2(INT *c, 
                    FLOAT *W, 
                    CompnentDistribution **MixTheta, 
                    FLOAT *RW, 
                    CompnentDistribution **RMixTheta, 
                    FLOAT *VW, 
                    CompnentDistribution **VMixTheta,
                    FLOAT **P,
                    FLOAT *am_opt) 
{
    INT i, l, Error = E_OK;

    FLOAT A = (FLOAT)0.0, B = (FLOAT)0.0;

    Error = Resids(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P);
    
    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < c_; l++) {
        A += RW[l] * RW[l];
        B += RW[l] * VW[l];
        for (i = 0; i < length_pdf_; i++) {
            switch(VMixTheta[l]->pdf_[i]) {
            case pfNormal:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfWeibull:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfGamma:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfGumbel:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfvonMises:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfBinomial:
                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += RMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfPoisson:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                break;
            case pfDirac: case pfUniform:
                break;
            }
        }
    }

    *am_opt = Max(A / (B + Eps), (FLOAT)1.0);

EEXIT:

    E_RETURN(Error);

}

INT Emmix::OptSLen3(INT *c, 
                    FLOAT *W, 
                    CompnentDistribution **MixTheta, 
                    FLOAT *RW, 
                    CompnentDistribution **RMixTheta, 
                    FLOAT *VW, 
                    CompnentDistribution **VMixTheta,
                    FLOAT **P,
                    FLOAT *am_opt) 
{
    INT i, l, Error = E_OK;

    FLOAT A = (FLOAT)0.0, B = (FLOAT)0.0;

    Error = Resids(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P);
    
    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < c_; l++) {
        A += RW[l] * RW[l];
        B += VW[l] * VW[l];
        for (i = 0; i < length_pdf_; i++) {
            switch(VMixTheta[l]->pdf_[i]) {
            case pfNormal:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfWeibull:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfGamma:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfGumbel:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfvonMises:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfBinomial:
                A += RMixTheta[l]->Theta_[1][i] * RMixTheta[l]->Theta_[1][i];

                B += VMixTheta[l]->Theta_[1][i] * VMixTheta[l]->Theta_[1][i];

                break;
            case pfPoisson:
                A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];

                B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

                break;
            case pfDirac: case pfUniform:
                break;
            }
        }
    }

    *am_opt = Max(sqrt(A) / (sqrt(B) + Eps), (FLOAT)1.0);

EEXIT:

    E_RETURN(Error);

}
INT Emmix::OptSLen(INT *c, 
                   FLOAT *W, 
                   CompnentDistribution **MixTheta, 
                   FLOAT *RW, 
                   CompnentDistribution **RMixTheta, 
                   FLOAT *VW, 
                   CompnentDistribution **VMixTheta,
                   FLOAT **P,
                   FLOAT *am_opt) {
    INT Error = E_OK;

    switch (accel_eq_) {
    case acc_param_euclidean_squared:
        Error = OptSLen1(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P, am_opt);
        if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); *am_opt = (FLOAT)1.0;
        }
        break;
    case acc_param_euclidean_weighted:
        Error = OptSLen2(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P, am_opt);
        if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); *am_opt = (FLOAT)1.0;
        }
        break;
    case acc_param_euclidean_semi_weighted:
        Error = OptSLen3(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P, am_opt);
        if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); *am_opt = (FLOAT)1.0;
        }
        break;
    case acc_param_none: default:
        *am_opt = (FLOAT)1.0;

        break;
    }

    E_RETURN(Error);
}

INT Emmix::CheckSLenConv(INT *c, 
                         FLOAT *W, 
                         CompnentDistribution **MixTheta, 
                         FLOAT *RW, 
                         CompnentDistribution **RMixTheta, 
                         FLOAT *VW, 
                         CompnentDistribution **VMixTheta, 
                         FLOAT *am_opt) 
{
    // Init structs
    FLOAT *_W = NULL, logL, logLAcc, am;
    CompnentDistribution **_MixTheta = NULL;
    INT i,  Error = E_OK; 

    _W = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    E_CHECK(NULL == _W, E_MEM);

    _MixTheta = new CompnentDistribution* [(unsigned INT)cmax_];

    E_CHECK(NULL == _MixTheta, E_MEM);

    for (i = 0; i < cmax_; i++) {
        _MixTheta[i] = new CompnentDistribution(this);

        E_CHECK(NULL == _MixTheta[i], E_MEM);

        Error = _MixTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        E_CHECK(Error != E_OK, Error);
    }
    
    // Copy the structs and floats

    for (i = 0; i < *c; i++) {
        _W[i] = W[i];

        Error = _MixTheta[i]->Memmove(MixTheta[i]);

        E_CHECK(Error != E_OK, Error);
    }

    // Estimate the logl

    Error = LogLikelihood(*c, W, MixTheta, &logL);

    E_CHECK(Error != E_OK, Error);

    if (accel_ == acc_square) {
        Error = UpdateMixParamsSqr(c, _W, _MixTheta, RW, RMixTheta, VW, VMixTheta, *am_opt);

        E_CHECK(Error != E_OK, Error);

    } else 
    if (accel_ == acc_stem) {
        Error = UpdateMixParamsSt(c, _W, _MixTheta, RW, RMixTheta, VW, VMixTheta, *am_opt);

        E_CHECK(Error != E_OK, Error);

    } else { 
        Error = E_ARG;

        E_CHECK(Error != E_OK, Error);
    }

    Error = LogLikelihood(*c, _W, _MixTheta, &logLAcc);

    E_CHECK(Error != E_OK, Error);

    // Check if logL increases if not take half step until 1

    while ((logLAcc < logL) && fabs(logLAcc - logL) > TOL_) {

        am = *am_opt;

        *am_opt = (*am_opt + 1) / (FLOAT)2.0;

        if (fabs(*am_opt - am) < Eps) {
            Error = E_CON;

            E_CHECK(Error != E_OK, Error);
        }

        for (i = 0; i < *c; i++) {
            _W[i] = W[i];

            Error = _MixTheta[i]->Memmove(MixTheta[i]);

            E_CHECK(Error != E_OK, Error);
        }

        if (accel_ == acc_square) {
            Error = UpdateMixParamsSqr(c, _W, _MixTheta, RW, RMixTheta, VW, VMixTheta, *am_opt);

            E_CHECK(Error != E_OK, Error);

        } else 
        if (accel_ == acc_stem) {
            Error = UpdateMixParamsSt(c, _W, _MixTheta, RW, RMixTheta, VW, VMixTheta, *am_opt);

            E_CHECK(Error != E_OK, Error);
        } else { 
            Error = E_ARG;

            E_CHECK(Error != E_OK, Error);
        }

        Error = LogLikelihood(*c, _W, _MixTheta, &logLAcc);

        E_CHECK(Error != E_OK, Error);
    }



EEXIT:
    if (_MixTheta) {
        for (i = 0; i < cmax_; i++) {
            if (_MixTheta[i]) delete _MixTheta[i];
        }

        delete [] _MixTheta;
    }
    
    if (_W) free(_W);
    
    E_RETURN(Error);
}


INT Emmix::CheckSLenConvFast(INT *c, 
                         FLOAT *W, 
                         CompnentDistribution **MixTheta, 
                         FLOAT *RW, 
                         CompnentDistribution **RMixTheta, 
                         FLOAT *VW, 
                         CompnentDistribution **VMixTheta, 
                         FLOAT *am_opt) 
{
    // Init structs
    FLOAT *_W = NULL, logL, logLAcc;
    CompnentDistribution **_MixTheta = NULL;
    INT i,  Error = E_OK; 

    _W = (FLOAT*)malloc(cmax_ * sizeof(FLOAT));

    E_CHECK(NULL == _W, E_MEM);

    _MixTheta = new CompnentDistribution* [(unsigned INT)cmax_];

    E_CHECK(NULL == _MixTheta, E_MEM);

    for (i = 0; i < cmax_; i++) {
        _MixTheta[i] = new CompnentDistribution(this);

        E_CHECK(NULL == _MixTheta[i], E_MEM);

        Error = _MixTheta[i]->Realloc(length_pdf_, length_Theta_, length_theta_);

        E_CHECK(Error != E_OK, Error);
    }
    
    // Copy the structs and floats

    for (i = 0; i < *c; i++) {
        _W[i] = W[i];

        Error = _MixTheta[i]->Memmove(MixTheta[i]);

        E_CHECK(Error != E_OK, Error);
    }

    // Estimate the logl

    Error = LogLikelihood(*c, W, MixTheta, &logL);

    E_CHECK(Error != E_OK, Error);

    if (accel_ == acc_square) {
        Error = UpdateMixParamsSqr(c, _W, _MixTheta, RW, RMixTheta, VW, VMixTheta, *am_opt);

        E_CHECK(Error != E_OK, Error);

    } else 
    if (accel_ == acc_stem) {
        Error = UpdateMixParamsSt(c, _W, _MixTheta, RW, RMixTheta, VW, VMixTheta, *am_opt);

        E_CHECK(Error != E_OK, Error);

    } else { 
        Error = E_ARG;

        E_CHECK(Error != E_OK, Error);
    }

    Error = LogLikelihood(*c, _W, _MixTheta, &logLAcc);

    E_CHECK(Error != E_OK, Error);

    // Check if logL increases if not take half step until 1

    if (logLAcc < logL) {
        *am_opt = (FLOAT)1.0;
        Error = E_CON;
        goto EEXIT;
    }

EEXIT:
    if (_MixTheta) {
        for (i = 0; i < cmax_; i++) {
            if (_MixTheta[i]) delete _MixTheta[i];
        }

        delete [] _MixTheta;
    }
    
    if (_W) free(_W);
    
    E_RETURN(Error);
}

INT Emmix::EstimateAlphaNaive(INT i, FLOAT am, FLOAT *am_opt) 
{
    INT Error = E_OK;

    switch (accel_eq_) {
    case acc_param_growth_lin:
        *am_opt = (FLOAT)1.0 + (FLOAT)i * (am - (FLOAT)1.0) / (FLOAT)max_iter_;

        break;
    case acc_param_growth_exp:
        *am_opt = exp((FLOAT)i * log(am) / (FLOAT) max_iter_);

        break;
    case acc_param_decay_lin:
        *am_opt = am + (FLOAT)i * (- am + (FLOAT)1.0) / (FLOAT)max_iter_;

        break;
    case acc_param_decay_exp:
        *am_opt = am * exp((FLOAT)i * log((FLOAT)1.0 / am) / (FLOAT) max_iter_);

        break;
    case acc_param_none: case acc_param_euclidean_semi_weighted: case acc_param_euclidean_squared: case acc_param_euclidean_weighted:
        Error = E_ARG;
    }

    return Error;
}

INT Emmix::AitkenAcceleratedLogL(FLOAT *LogL, FLOAT *LogLPrev)
{
    INT Error = E_OK; FLOAT a;

    //LogL_[1] = LogL_[0];

    //LogL_[1] = LogL_[2];

    //LogL_[2] = *LogL;

    if (fabs(LogLPrev[1] - LogLPrev[0]) < FLOAT_MIN || fabs(LogLPrev[2] - LogLPrev[1]) < FLOAT_MIN) {
        Error = E_CON;

        goto EEXIT;
    }

    a = (LogLPrev[2] - LogLPrev[1]) / (LogLPrev[1] - LogLPrev[0]);

    a = (FLOAT)1.0 - a;

    if (fabs(a) < FLOAT_MIN) {
        Error = E_CON;

        goto EEXIT;
    }

    *LogL = LogLPrev[1] + (FLOAT)1.0 / a * (LogLPrev[2] - LogLPrev[1]);


EEXIT:

    return Error;


}

INT Emmix::Converged(FLOAT LogLNew, FLOAT LogLOld, INT n) {
    INT Error = E_CON;

    if (loglest_ == likelihood_aitken_bohning) {
        AitkenAcceleratedLogL(&LogLOld, LogL_);

        LogL_[1] = LogL_[0];

        LogL_[1] = LogL_[2];

        LogL_[2] = LogLNew;

        AitkenAcceleratedLogL(&LogLNew, LogL_);

    } else 
    if (loglest_ == likelihood_aitken_lindsay) {
        LogL_[1] = LogL_[0];

        LogL_[1] = LogL_[2];

        LogL_[2] = LogLNew;

        LogLOld = LogLNew;

        if (AitkenAcceleratedLogL(&LogLNew, LogL_) == E_CON) LogLOld = LogL_[1];
    } else
    if (loglest_ == likelihood_aitken_nicholas) {
        LogL_[1] = LogL_[0];

        LogL_[1] = LogL_[2];

        LogL_[2] = LogLNew;

        AitkenAcceleratedLogL(&LogLNew, LogL_);
    } else {

        LogL_[1] = LogL_[0];

        LogL_[1] = LogL_[2];

        LogL_[2] = LogLNew;
    }
    

    switch (toltype_) {
    case likelihood_standard:
        if ((FLOAT)fabs(LogLNew - LogLOld) <= TOL_) Error = E_OK;
        break;
    case likelihood_normalised:
        if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)n <= TOL_) Error = E_OK;
        break;
    case likelihood_percentage:
        if ((FLOAT)fabs((LogLNew - LogLOld) / Max(LogLNew, (FLOAT)1.0))  <= TOL_) Error = E_OK;
        break;
    default:
        if ((FLOAT)fabs(LogLNew - LogLOld) <= TOL_) Error = E_OK;
        break;
    }

    return Error;
}

// Performs the standard EM algorithm (reference).

INT Emmix::EM()
{
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;
    INT   i, Error = E_OK;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    E_CHECK(Error != E_OK, Error);

    LogL_[2] = LogL_[1] = LogL_[0] = LogLOld;

    n_iter_ = 0;

    for (i = 0; i < max_iter_; i++) {
        Error = ExpectationStep(c_, W_, MixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = MaximizationStep(&c_, W_, MixTheta_, dW_, dMixTheta_, TW_, TMixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != E_OK, Error);

        if (Converged(LogLNew, LogLOld, n_) == E_OK) break;

        /*if (loglest_ == likelihood_aitken) AitkenAcceleratedLogL(&LogLNew);

        if ((FLOAT)fabs(LogLNew - LogLOld) <= TOL_) break;

        if (loglest_ == likelihood_aitken) LogLNew = LogL_[2];*/

        //if ((FLOAT)fabs(LogLNew - LogLOld) <= TOL_) break;

        LogLOld = LogLNew;
    }

    //n_iter_ = i;

    aam_ /= (FLOAT)i;

EEXIT:

    E_RETURN(Error);
} // EM

// Performs the Expectation-Conditional-Maximization algorithm (reference).

INT Emmix::ECM()
{
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;
    INT   i, Error = E_OK;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    E_CHECK(Error != E_OK, Error);

    LogL_[2] = LogL_[1] = LogL_[0] = LogLOld;

    n_iter_ = 0;

    for (i = 0; i < max_iter_; i++) {
        Error = ConditionalStep(c_, W_, MixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = MaximizationStep(&c_, W_, MixTheta_, dW_, dMixTheta_, TW_, TMixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != E_OK, Error);

        //if (loglest_ == likelihood_aitken) AitkenAcceleratedLogL(&LogLNew);

        //if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        if (Converged(LogLNew, LogLOld, n_) == E_OK) break; 

        LogLOld = LogLNew;
    }

    //n_iter_ = i;

    aam_ /= (FLOAT)i;

EEXIT:

    E_RETURN(Error);
} // ECM

INT Emmix::SEM()
{
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0;
    INT   i, Error = E_OK;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    E_CHECK(Error != E_OK, Error);

    LogL_[2] = LogL_[1] = LogL_[0] = LogLOld;

    n_iter_ = 0;

    for (i = 0; i < max_iter_; i++) {
        Error = StochasticStep(c_, W_, MixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = MaximizationStep(&c_, W_, MixTheta_, dW_, dMixTheta_, TW_, TMixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != E_OK, Error);

        //if (loglest_ == likelihood_aitken) AitkenAcceleratedLogL(&LogLNew);

        //if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        if (Converged(LogLNew, LogLOld, n_) == E_OK) break;

        LogLOld = LogLNew;
    }

    //n_iter_ = i;

    aam_ /= (FLOAT)i;

EEXIT:

    E_RETURN(Error);
}

INT Emmix::ECM_EM()
{
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0, am = am_;
    INT   i, max_iter = max_iter_,  Error = E_OK;

    EmAccelerationType_e accel = accel_; accel_ = acc_fixed; am_ = (FLOAT)1.0;

    EmLikelihoodEstimateType_e loglest = loglest_; loglest_ = likelihood_standard;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    E_CHECK(Error != E_OK, Error);

    LogL_[2] = LogL_[1] = LogL_[0] = LogLOld;

    n_iter_ = 0;

    for (i = 0; i < max_iter; i++) {
        Error = ConditionalStep(c_, W_, MixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = MaximizationStep(&c_, W_, MixTheta_, dW_, dMixTheta_, TW_, TMixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != E_OK, Error);

        //if (loglest_ == likelihood_aitken) AitkenAcceleratedLogL(&LogLNew);

        //if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        if (Converged(LogLNew, LogLOld, n_) == E_OK) break; 

        LogLOld = LogLNew;
    }

    max_iter = max_iter - i; aam_ = (FLOAT)0.0; accel_ = accel; am_ = am; loglest_ = loglest;

    for (i = 0; i < max_iter; i++) {
        Error = ExpectationStep(c_, W_, MixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = MaximizationStep(&c_, W_, MixTheta_, dW_, dMixTheta_, TW_, TMixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != E_OK, Error);

        //if (loglest_ == likelihood_aitken) AitkenAcceleratedLogL(&LogLNew);

        //if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        //if (loglest_ == likelihood_aitken) LogLNew = LogL_[2];

        if (Converged(LogLNew, LogLOld, n_) == E_OK) break;

        LogLOld = LogLNew;
    }

    //n_iter_ = i;

    aam_ /= (FLOAT)i;

EEXIT:

    E_RETURN(Error);
}

INT Emmix::SEM_EM()
{
    FLOAT LogLOld = (FLOAT)0.0, LogLNew = (FLOAT)0.0, am = am_;
    INT   i, max_iter = max_iter_, Error = E_OK;

    EmAccelerationType_e accel = accel_; accel_ = acc_fixed; am_ = (FLOAT)1.0;

    EmLikelihoodEstimateType_e loglest = loglest_; loglest_ = likelihood_standard;

    Error = LogLikelihood(c_, W_, MixTheta_, &LogLOld);

    E_CHECK(Error != E_OK, Error);

    LogL_[2] = LogL_[1] = LogL_[0] = LogLOld;

    n_iter_ = 0;

    for (i = 0; i < max_iter; i++) {
        Error = StochasticStep(c_, W_, MixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = MaximizationStep(&c_, W_, MixTheta_, dW_, dMixTheta_, TW_, TMixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != E_OK, Error);

        //if (loglest_ == likelihood_aitken) AitkenAcceleratedLogL(&LogLNew);

        //if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        if (Converged(LogLNew, LogLOld, n_) == E_OK) break; 

        LogLOld = LogLNew;
    }

    max_iter = max_iter_ - i; aam_ = (FLOAT)0.0; accel_ = accel; am_ = am; loglest_ = loglest;

    for (i = 0; i<max_iter; i++) {
        Error = ExpectationStep(c_, W_, MixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = MaximizationStep(&c_, W_, MixTheta_, dW_, dMixTheta_, TW_, TMixTheta_, P_);

        E_CHECK(Error != E_OK, Error);

        Error = LogLikelihood(c_, W_, MixTheta_, &LogLNew);

        E_CHECK(Error != E_OK, Error);

        //if (loglest_ == likelihood_aitken) AitkenAcceleratedLogL(&LogLNew);

        //if ((FLOAT)fabs(LogLNew - LogLOld) / (FLOAT)fabs(LogLNew) <= TOL_) break;

        //if (loglest_ == likelihood_aitken) LogLNew = LogL_[2];

        if (Converged(LogLNew, LogLOld, n_) == E_OK) break;

        LogLOld = LogLNew;
    }

    //n_iter_ = i;

    aam_ /= (FLOAT)i;

EEXIT:

    E_RETURN(Error);
}

// Runs the EM algorithm or its variant.

INT Emmix::Run(INT *c, FLOAT *W, CompnentDistribution **MixTheta)
{
    INT i, Error = E_OK;

    c_ = *c;

    for (i = 0; i < c_; i++) {
        W_[i] = W[i];

        Error = MixTheta_[i]->Memmove(MixTheta[i]);

        E_CHECK(Error != E_OK, Error);
    }

    switch (variant_) {
    case varEM:
        Error = EM();

        E_CHECK(Error != E_OK, Error);

        break;
    case varECM:
        Error = ECM();

        E_CHECK(Error != E_OK, Error);

        break;
    case varSEM:
        Error = SEM();

        E_CHECK(Error != E_OK, Error);

        break;
    case varECM_EM:
        Error = ECM_EM();

        E_CHECK(Error != E_OK, Error);

        break;
    case varSEM_EM:
        Error = SEM_EM();

        E_CHECK(Error != E_OK, Error);

        break;
    }

    *c = c_;

    for (i = 0; i < *c; i++) {
        W[i] = W_[i];

        Error = MixTheta[i]->Memmove(MixTheta_[i]);

        E_CHECK(Error != E_OK, Error);
    }

EEXIT:

    E_RETURN(Error);
} // Run

// Returns logarithm of component p.d.f..

INT Emmix::LogComponentPdf(INT                  j,         // Indey of observation.  
                           FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                           CompnentDistribution *CmpTheta, // Component parameters.
                           FLOAT                *CmpPdf)   // Logarithm of component probability density.
{
    FLOAT p, Theta, y, ypb;
    INT   i, k, n, Error = E_OK;

    *CmpPdf = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        switch (CmpTheta->pdf_[i]) {
        case pfNormal:
            y = (Y[i][j] - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

            *CmpPdf += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]);

            break;
        case pfTNormal:
            break;
        case pfLognormal:
            if (Y[i][j] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i][j]) - CmpTheta->Theta_[0][i]) / (Sqrt2 * CmpTheta->Theta_[1][i]); y *= y;

                *CmpPdf += -y - LogSqrtPi2 - (FLOAT)log(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpPdf = -FLOAT_MAX;
            }

            break;
        case pfWeibull:
            if (Y[i][j] > FLOAT_MIN) {
                ypb = (FLOAT)exp(CmpTheta->Theta_[1][i] * (FLOAT)log(Y[i][j] / CmpTheta->Theta_[0][i]));

                *CmpPdf += (FLOAT)log(CmpTheta->Theta_[1][i]) + (FLOAT)log(ypb) - ypb - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpPdf = -FLOAT_MAX;
            }

            break;
        case pfGamma:
            if (Y[i][j] > FLOAT_MIN) {
                ypb = Y[i][j] / CmpTheta->Theta_[0][i];

                *CmpPdf += CmpTheta->Theta_[1][i] * (FLOAT)log(ypb) - ypb - Gammaln(CmpTheta->Theta_[1][i]) - (FLOAT)log(Y[i][j]);
            }
            else {
                *CmpPdf = -FLOAT_MAX;
            }

            break;
        case pfGumbel:
            ypb = CmpTheta->Theta_[2][i] * (Y[i][j] - CmpTheta->Theta_[0][i]) / CmpTheta->Theta_[1][i];

            *CmpPdf += ypb - (FLOAT)exp(ypb) - (FLOAT)log(CmpTheta->Theta_[1][i]);

            break;
        case pfvonMises:
            if ((Y[i][j] < (FLOAT)0.0) || (Y[i][j] > Pi2)) {
                *CmpPdf = -FLOAT_MAX;
            }
            else {
                *CmpPdf += CmpTheta->Theta_[1][i] * (FLOAT)cos(Y[i][j] - CmpTheta->Theta_[0][i]) - LogPi2 - (FLOAT)log(BesselI0(CmpTheta->Theta_[1][i]));
            }

            break;
        case pfBinomial:
            k = (INT)Y[i][j]; n = (INT)CmpTheta->Theta_[0][i]; p = CmpTheta->Theta_[1][i];

            if (k < 0) {
                *CmpPdf = -FLOAT_MAX;
            }
            else
            if (k == 0)
                *CmpPdf += n * (FLOAT)log((FLOAT)1.0 - p);
            else
            if (k == n)
                *CmpPdf += n * (FLOAT)log(p);
            else
            if (k > n) {
                *CmpPdf = -FLOAT_MAX;
            }
            else
               *CmpPdf += Gammaln(n + (FLOAT)1.0) - Gammaln(k + (FLOAT)1.0) - Gammaln(n - k + (FLOAT)1.0) +
                   k * (FLOAT)log(p) + (n - k) * (FLOAT)log((FLOAT)1.0 - p);

            break;
        case pfPoisson:
            k = (INT)Y[i][j]; Theta = CmpTheta->Theta_[0][i];

            *CmpPdf += k * (FLOAT)log(Theta) - Theta - Gammaln(k + (FLOAT)1.0);

            break;
        case pfDirac:
            if ((FLOAT)fabs(Y[i][j] - CmpTheta->Theta_[0][i]) > FLOAT_MIN) {
                *CmpPdf = -FLOAT_MAX;
            }
            else {
                *CmpPdf += (FLOAT)0.0;
            }

            break;
        case pfUniform:
            if ((Y[i][j] > CmpTheta->Theta_[1][i]) || (Y[i][j] < CmpTheta->Theta_[0][i])) {
                *CmpPdf = -FLOAT_MAX;
            }
            else {
                *CmpPdf -= (FLOAT)log(CmpTheta->Theta_[1][i] - CmpTheta->Theta_[0][i]);
            }
        }
    }

    E_RETURN(Error);
} // LogComponentPdf

// Updates mixture model parameters with appropriate increment.

INT Emmix::UpdateMixParams(INT                  *c,          // Number of components. 
                           FLOAT                *W,          // Mixture model weight values.
                           CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                           FLOAT                *dW,         // Update increment of mixture model weights.
                           CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                           FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, j, l, Error = E_OK;

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

INT Emmix::UpdateMixParamsSt(INT                  *c,          // Number of components. 
                              FLOAT                *W,          // Mixture model weight values.
                              CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                              FLOAT                *RW,         // Update increment of mixture model weights.
                              CompnentDistribution **RMixTheta, // Update increment of mixture model distribution parameter values.
                              FLOAT                *VW,         // Update increment of mixture model weights.
                              CompnentDistribution **VMixTheta, // Update increment of mixture model distribution parameter values.
                              FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, j, l, Error = E_OK;

    for (l = 0; l < *c; l++) {
        W[l] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RW[l] + am * am * VW[l];

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

/// Panic Branislav
        if (W[l] < FLOAT_MIN) {
            switch (merge_) {
            case merge_naive:
                for (j = l; j < *c - 1; j++) {
                    RW[j] = RW[j + 1];

                    VW[j] = VW[j + 1];

                    W[j] = W[j + 1];

                    for (i = 0; i < length_pdf_; i++) {
                        MixTheta[j]->Theta_[0][i] = MixTheta[j + 1]->Theta_[0][i];

                        RMixTheta[j]->Theta_[0][i] = RMixTheta[j + 1]->Theta_[0][i];

                        VMixTheta[j]->Theta_[0][i] = VMixTheta[j + 1]->Theta_[0][i];

                        MixTheta[j]->Theta_[1][i] = MixTheta[j + 1]->Theta_[1][i];

                        RMixTheta[j]->Theta_[1][i] = RMixTheta[j + 1]->Theta_[1][i];

                        VMixTheta[j]->Theta_[1][i] = VMixTheta[j + 1]->Theta_[1][i];

                        MixTheta[j]->Theta_[2][i] = MixTheta[j + 1]->Theta_[2][i];

                        RMixTheta[j]->Theta_[2][i] = RMixTheta[j + 1]->Theta_[2][i];

                        VMixTheta[j]->Theta_[2][i] = VMixTheta[j + 1]->Theta_[2][i];
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
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfWeibull:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[0][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[0][i] = Eps;
                }

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfGamma:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[0][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[0][i] = Eps;
                }

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfGumbel:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfvonMises:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfBinomial:
                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < (FLOAT)0.0) {
                    MixTheta[l]->Theta_[1][i] = (FLOAT)0.0;
                }
                else
                if (MixTheta[l]->Theta_[1][i] > (FLOAT)1.0) {
                    MixTheta[l]->Theta_[1][i] = (FLOAT)1.0;
                }

                break;
            case pfPoisson:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

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

INT Emmix::UpdateMixParamsSqr(INT                  *c,          // Number of components. 
                              FLOAT                *W,          // Mixture model weight values.
                              CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                              FLOAT                *RW,         // Update increment of mixture model weights.
                              CompnentDistribution **RMixTheta, // Update increment of mixture model distribution parameter values.
                              FLOAT                *VW,         // Update increment of mixture model weights.
                              CompnentDistribution **VMixTheta, // Update increment of mixture model distribution parameter values.
                              FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, j, l, Error = E_OK;

    for (l = 0; l < *c; l++) {
        W[l] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RW[l] + am * am * VW[l];

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

/// Panic Branislav
        if (W[l] < FLOAT_MIN) {
            switch (merge_) {
            case merge_naive:
                for (j = l; j < *c - 1; j++) {
                    RW[j] = RW[j + 1];

                    VW[j] = VW[j + 1];

                    W[j] = W[j + 1];

                    for (i = 0; i < length_pdf_; i++) {
                        MixTheta[j]->Theta_[0][i] = MixTheta[j + 1]->Theta_[0][i];

                        RMixTheta[j]->Theta_[0][i] = RMixTheta[j + 1]->Theta_[0][i];

                        VMixTheta[j]->Theta_[0][i] = VMixTheta[j + 1]->Theta_[0][i];

                        MixTheta[j]->Theta_[1][i] = MixTheta[j + 1]->Theta_[1][i];

                        RMixTheta[j]->Theta_[1][i] = RMixTheta[j + 1]->Theta_[1][i];

                        VMixTheta[j]->Theta_[1][i] = VMixTheta[j + 1]->Theta_[1][i];

                        MixTheta[j]->Theta_[2][i] = MixTheta[j + 1]->Theta_[2][i];

                        RMixTheta[j]->Theta_[2][i] = RMixTheta[j + 1]->Theta_[2][i];

                        VMixTheta[j]->Theta_[2][i] = VMixTheta[j + 1]->Theta_[2][i];
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
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfWeibull:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[0][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[0][i] = Eps;
                }

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfGamma:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[0][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[0][i] = Eps;
                }

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfGumbel:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfvonMises:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < Eps) {
                    W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][i] = Eps;
                }

                break;
            case pfBinomial:
                MixTheta[l]->Theta_[1][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][i] + am * am * VMixTheta[l]->Theta_[1][i];

                if (MixTheta[l]->Theta_[1][i] < (FLOAT)0.0) {
                    MixTheta[l]->Theta_[1][i] = (FLOAT)0.0;
                }
                else
                if (MixTheta[l]->Theta_[1][i] > (FLOAT)1.0) {
                    MixTheta[l]->Theta_[1][i] = (FLOAT)1.0;
                }

                break;
            case pfPoisson:
                MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

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

INT Emmix::F(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT **P) {
    FLOAT A[5], *C = NULL, dC, dM, *M = NULL, T[2];
    INT   i, j, k, l, error, Error = E_OK;

    M = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == M, E_MEM);

    C = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == C, E_MEM);

    for (l = 0; l < c; l++) {
        dW[l] = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) {
            dW[l] += Y_[length_pdf_][j] * P[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < nr_; j++) {
                    if (Y_[length_pdf_][j] > FLOAT_MIN) {
                        M[i] += Y_[length_pdf_][j] * P[l][j] * Y_[i][j];
                    }
                }

                M[i] = M[i] / (dW[l] + FLOAT_MIN);

                dMixTheta[l]->Theta_[0][i] = M[i] - MixTheta[l]->Theta_[0][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                for (j = 0; j < nr_; j++) {
                    if ((Y_[length_pdf_][j] > FLOAT_MIN) && (Y_[i][j] > FLOAT_MIN)) {
                        M[i] += Y_[length_pdf_][j] * P[l][j] * (FLOAT)log(Y_[i][j]);
                    }
                }

                M[i] = M[i] / (dW[l] + FLOAT_MIN);

                dMixTheta[l]->Theta_[0][i] = M[i] - MixTheta[l]->Theta_[0][i];

                break;
            case pfWeibull:
                M[i] = MixTheta[l]->Theta_[1][i];

                j = 1; Error = E_CON;

                while ((j <= ItMax) && (Error != E_OK)) {
                    memset(&A, 0, 5 * sizeof(FLOAT));

                    for (k = 0; k < nr_; k++) {
                        if ((Y_[length_pdf_][k] > FLOAT_MIN) && (Y_[i][k] > FLOAT_MIN)) {
                            T[0] = (FLOAT)log(Y_[i][k]);
                            T[1] = (FLOAT)exp(T[0] * M[i]);

                            A[1] += Y_[length_pdf_][k] * P[l][k] * T[0];
                            A[2] += Y_[length_pdf_][k] * P[l][k] * T[1];
                            A[3] += Y_[length_pdf_][k] * P[l][k] * T[1] * T[0];
                            A[4] += Y_[length_pdf_][k] * P[l][k] * T[1] * T[0] * T[0];
                        }
                    }

                    A[0] = dW[l] + FLOAT_MIN;

                    T[0] = A[0] / A[2];

                    dM = (A[0] / M[i] + A[1] - T[0] * A[3]) / (T[0] * (A[3] * A[3] / A[2] - A[4]) - A[0] / (M[i] * M[i]));

                    M[i] -= dM;

                    E_CHECK(IsNan(dM) || IsInf(dM), E_CON);

                    if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = E_OK;

                    j++;
                }

                dMixTheta[l]->Theta_[1][i] = M[i] - MixTheta[l]->Theta_[1][i];

                break;
            case pfGamma:
                M[i] = MixTheta[l]->Theta_[1][i];

                memset(&A, 0, 4 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    if ((Y_[length_pdf_][j] > FLOAT_MIN) && (Y_[i][j] > FLOAT_MIN)) {
                        A[1] += Y_[length_pdf_][j] * P[l][j] * Y_[i][j];
                        A[2] += Y_[length_pdf_][j] * P[l][j] * (FLOAT)log(Y_[i][j]);
                    }
                }

                A[0] = dW[l] + FLOAT_MIN;

                A[3] = (FLOAT)log(A[1] / A[0]) - A[2] / A[0];

                j = 1; Error = E_CON;

                while ((j <= ItMax) && (Error != E_OK)) {
                    error = Digamma(M[i], &T[0]);

                    E_CHECK(error != E_OK, error);
                    
                    error = Digamma(M[i] + Eps, &T[1]);

                    E_CHECK(error != E_OK, error);

                    T[1] = (T[1] - T[0]) / Eps;

                    dM = (T[0] + A[3] - (FLOAT)log(M[i]) ) / (T[1] - (FLOAT)1.0 / M[i]);

                    M[i] -= dM;

                    E_CHECK(IsNan(dM) || IsInf(dM), E_CON);

                    if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = E_OK;

                    j++;
                }

                dMixTheta[l]->Theta_[1][i] = M[i] - MixTheta[l]->Theta_[1][i];

                break;
            case pfGumbel:
                M[i] = MixTheta[l]->Theta_[1][i];

                j = 1; Error = E_CON;

                while ((j <= ItMax) && (Error != E_OK)) {
                    memset(&A, 0, 5 * sizeof(FLOAT));

                    for (k = 0; k < nr_; k++) {
                        if (Y_[length_pdf_][k] > FLOAT_MIN) {
                            T[0] = (FLOAT)exp(MixTheta[l]->Theta_[2][i] * Y_[i][k] / M[i]);

                            A[1] += Y_[length_pdf_][k] * P[l][k] * Y_[i][k];
                            A[2] += Y_[length_pdf_][k] * P[l][k] * T[0];
                            A[3] += Y_[length_pdf_][k] * P[l][k] * Y_[i][k] * T[0];
                            A[4] += Y_[length_pdf_][k] * P[l][k] * Y_[i][k] * Y_[i][k] * T[0];
                        }
                    }

                    A[0] = dW[l] + FLOAT_MIN;

                    T[0] = A[0] / A[2]; T[1] = A[3] / A[2];

                    dM = (M[i] * A[0] + MixTheta[l]->Theta_[2][i] * (A[1] - T[0] * A[3])) / (A[0] + T[0] * (A[4] - T[1] * A[3]) / (M[i] * M[i]));

                    M[i] -= dM;

                    E_CHECK(IsNan(dM) || IsInf(dM), E_CON);

                    if ((FLOAT)fabs(dM) < Max(Eps * (FLOAT)fabs(M[i]), Eps)) Error = E_OK;

                    j++;
                }

                dMixTheta[l]->Theta_[1][i] = M[i] - MixTheta[l]->Theta_[1][i];

                break;
            case pfvonMises:
                memset(&A, 0, 3 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    if (Y_[length_pdf_][j] > FLOAT_MIN) {
                        A[0] += Y_[length_pdf_][j] * P[l][j] * (FLOAT)cos(Y_[i][j]);
                        A[1] += Y_[length_pdf_][j] * P[l][j] * (FLOAT)sin(Y_[i][j]);
                    }
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
                    E_CHECK(1, E_CON);
                }

                dMixTheta[l]->Theta_[0][i] = M[i] - MixTheta[l]->Theta_[0][i];

                break;
            case pfBinomial:
                dMixTheta[l]->Theta_[0][i] = (FLOAT)0.0;

                break;
            case pfPoisson:
                for (j = 0; j < nr_; j++) {
                    M[i] += Y_[length_pdf_][j] * P[l][j] * Y_[i][j];
                }

                M[i] = M[i] / (dW[l] + FLOAT_MIN);

                dMixTheta[l]->Theta_[0][i] = M[i] - MixTheta[l]->Theta_[0][i];

                break;
            case pfDirac:
                dMixTheta[l]->Theta_[0][i] = (FLOAT)0.0;

                break;
            case pfUniform:
                dMixTheta[l]->Theta_[0][i] = (FLOAT)0.0;
            }
        }

        memset(C, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            switch (MixTheta[l]->pdf_[i]) {
            case pfNormal:
                for (j = 0; j < nr_; j++) {
                    if (Y_[length_pdf_][j] > FLOAT_MIN) {
                        C[i] += Y_[length_pdf_][j] * P[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
                    }
                }

                C[i] = (FLOAT)sqrt(C[i] / (dW[l] + FLOAT_MIN));

                E_CHECK(C[i] < FLOAT_MIN || IsNan(C[i]), E_CON);

                dMixTheta[l]->Theta_[1][i] = C[i] - MixTheta[l]->Theta_[1][i];

                break;
            case pfTNormal:
                break;
            case pfLognormal:
                for (j = 0; j < nr_; j++) {
                    if ((Y_[length_pdf_][j] > FLOAT_MIN) && (Y_[i][j] > FLOAT_MIN)) {
                        A[0] = (FLOAT)log(Y_[i][j]);

                        C[i] += Y_[length_pdf_][j] * P[l][j] * A[0] * A[0];
                    }
                }

                C[i] = (FLOAT)sqrt(C[i] / (dW[l] + FLOAT_MIN) - M[i] * M[i]);

                E_CHECK(C[i] < FLOAT_MIN || IsNan(C[i]), E_CON);

                dMixTheta[l]->Theta_[1][i] = C[i] - MixTheta[l]->Theta_[1][i];

                break;
            case pfWeibull:
                memset(&A, 0, 2 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) if (Y_[i][j] > FLOAT_MIN) {
                    if ((Y_[length_pdf_][j] > FLOAT_MIN) && (Y_[i][j] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y_[i][j]);
                        T[1] = (FLOAT)exp(T[0] * M[i]);

                        A[1] += Y_[length_pdf_][j] * P[l][j] * T[1];
                    }
                }

                A[0] = dW[l] + FLOAT_MIN;

                T[0] = A[1] / A[0];

                C[i] = (FLOAT)exp((FLOAT)log(T[0]) / M[i]);

                dMixTheta[l]->Theta_[0][i] = C[i] - MixTheta[l]->Theta_[0][i];

                break;
            case pfGamma:
                memset(&A, 0, 2 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    if ((Y_[length_pdf_][j] > FLOAT_MIN) && (Y_[i][j] > FLOAT_MIN)) {
                        A[1] += Y_[length_pdf_][j] * P[l][j] * Y_[i][j];
                    }
                }

                A[0] = dW[l] + FLOAT_MIN;

                C[i] = A[1] / A[0] / M[i];

                dMixTheta[l]->Theta_[0][i] = C[i] - MixTheta[l]->Theta_[0][i];

                break;
            case pfGumbel:
                memset(&A, 0, 2 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    if (Y_[length_pdf_][j] > FLOAT_MIN) {
                        T[0] = (FLOAT)exp(MixTheta[l]->Theta_[2][i] * Y_[i][j] / M[i]);

                        A[1] += Y_[length_pdf_][j] * P[l][j] * T[0];
                    }
                }

                A[0] = dW[l] + FLOAT_MIN;

                T[0] = A[1] / A[0];

                C[i] = MixTheta[l]->Theta_[2][i] * M[i] * (FLOAT)log(T[0]);

                dMixTheta[l]->Theta_[0][i] = C[i] - MixTheta[l]->Theta_[0][i];

                break;
            case pfvonMises:
                memset(&A, 0, 3 * sizeof(FLOAT));

                for (j = 0; j < nr_; j++) {
                    A[1] += Y_[length_pdf_][j] * P[l][j] * (FLOAT)cos(Y_[i][j] - M[i]);
                }

                A[0] = dW[l] + FLOAT_MIN; A[2] = A[1] / A[0];

                C[i] = MixTheta[l]->Theta_[1][i];

                j = 1; Error = E_CON;

                while ((j <= ItMax) && (Error != E_OK)) {
                    A[0] = BesselI0(C[i]); A[1] = BesselI1(C[i]);

                    dC = (A[1] - A[2] * A[0]) / (A[0] - (A[2] + (FLOAT)1.0 / C[i]) * A[1]);

                    E_CHECK(IsNan(dC) || IsInf(dC), E_CON);

                    C[i] -= dC;

                    if ((FLOAT)fabs(dC) < Max(Eps * (FLOAT)fabs(C[i]), Eps)) Error = E_OK;

                    j++;
                }

                dMixTheta[l]->Theta_[1][i] = C[i] - MixTheta[l]->Theta_[1][i];

                break;
            case pfBinomial:
                for (j = 0; j < nr_; j++) {
                    C[i] += Y_[length_pdf_][j] * P[l][j] * Y_[i][j];
                }

                C[i] = C[i] / (dW[l] + FLOAT_MIN) / MixTheta[l]->Theta_[0][i];

                dMixTheta[l]->Theta_[1][i] = C[i] - MixTheta[l]->Theta_[1][i];

                break;
            case pfPoisson:
                dMixTheta[l]->Theta_[1][i] = (FLOAT)0.0;

                break;
            case pfDirac:
                dMixTheta[l]->Theta_[1][i] = (FLOAT)0.0;

                break;
            case pfUniform:
                dMixTheta[l]->Theta_[1][i] = (FLOAT)0.0;
            }
        }

        dW[l] = dW[l] / n_ - W[l];
    }


EEXIT:

    if (C) free(C); 
    
    if (M) free(M);

    n_iter_++;

    E_RETURN(Error);
}

INT Emmix::MaximizationStep(INT *c,
                            FLOAT *W,
                            CompnentDistribution **MixTheta,
                            FLOAT *dW,
                            CompnentDistribution **dMixTheta,
                            FLOAT *TW,
                            CompnentDistribution **TMixTheta,
                            FLOAT **P) 
{
    INT Error; FLOAT am_opt = (FLOAT)1.0;

    Error = F(*c, W, MixTheta, dW, dMixTheta, P);

    E_CHECK(Error != E_OK, Error);

    switch (accel_) {
    case acc_golden:
        Error = GoldenRatioSearch(c, W, MixTheta, dW, dMixTheta, TW, TMixTheta, &am_opt);

        if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); am_opt = (FLOAT)1.0;
        }
        
        break;
    case acc_line:
        Error = LineSearch(c, W, MixTheta, dW, dMixTheta, TW, TMixTheta, &am_opt);

        if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); am_opt = (FLOAT)1.0;
        }

        break;
    case acc_fixed:
        am_opt = am_;
        
        break;
    case acc_naive:
        Error = EstimateAlphaNaive(n_iter_, am_, &am_opt);

        if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); am_opt = (FLOAT)1.0;
        }

        break;
    case acc_stem: case acc_square:
        Error = OptSLen(c, W, MixTheta, dW, dMixTheta, TW, TMixTheta, P, &am_opt);

         if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); am_opt = (FLOAT)1.0;
            goto SKIP_CHECK;
        }

        Error = CheckSLenConv(c, W, MixTheta, dW, dMixTheta, TW, TMixTheta, &am_opt);

        if (Error != E_OK) {
            W_CHECK(Error != E_OK, 1); am_opt = (FLOAT)1.0;
        }
SKIP_CHECK:
        break;
    default:
        am_opt = (FLOAT)1.0;
        break;
    }
   
    if (accel_ == acc_square) {
        Error = UpdateMixParamsSqr(c, W, MixTheta, dW, dMixTheta, TW, TMixTheta, am_opt);

    } else 
    if (accel_ == acc_stem) {
        Error = UpdateMixParamsSt(c, W, MixTheta, dW, dMixTheta, TW, TMixTheta, am_opt);

    } else { 
        Error = UpdateMixParams(c, W, MixTheta, dW, dMixTheta, am_opt);

    }

    E_CHECK(Error != E_OK, Error);

    aam_ += am_opt;

EEXIT:
    
    E_RETURN(Error);
} // MaximizationStep

// Returns logarithm of component p.d.f..

INT Emmvnorm::LogComponentPdf(INT                  j,         // Indey of observation.  
                              FLOAT                **Y,       // Pointer to the input array [y0,...,yd-1,...]
                              CompnentDistribution *CmpTheta, // Component parameters.
                              FLOAT                *CmpPdf)   // Logarithm of component probability density.
{
    FLOAT y, yi, yk;
    INT   i, k, Error = E_OK;

    y = (FLOAT)0.0;

    for (i = 0; i < CmpTheta->length_pdf_; i++) {
        yi = Y[i][j] - CmpTheta->Theta_[0][i]; y += (FLOAT)0.5 * CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + i] * yi * yi;

        for (k = i + 1; k < CmpTheta->length_pdf_; k++) {
            yk = Y[k][j] - CmpTheta->Theta_[0][k]; y += CmpTheta->Theta_[2][i * CmpTheta->length_pdf_ + k] * yi * yk;
        }
    }
    
    *CmpPdf = -y - CmpTheta->length_pdf_ * LogSqrtPi2 - (FLOAT)0.5 * CmpTheta->Theta_[3][0];

    E_RETURN(Error);
} // LogComponentPdf

// Updates mixture model parameters with appropriate increment.

INT Emmvnorm::UpdateMixParams(INT                  *c,          // Number of components. 
                              FLOAT                *W,          // Mixture model weight values.
                              CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                              FLOAT                *dW,         // Update increment of mixture model weights.
                              CompnentDistribution **dMixTheta, // Update increment of mixture model distribution parameter values.
                              FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, ii, j, l, p, q, Error = E_OK;
    
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

        E_CHECK(Error != E_OK, Error);
S1:;
    }

EEXIT:

    E_RETURN(Error);
} // UpdateMixtureParameters

INT Emmvnorm::UpdateMixParamsSt(INT                  *c,          // Number of components. 
                                 FLOAT                *W,          // Mixture model weight values.
                                 CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                 FLOAT                *RW,         // Update increment of mixture model weights.
                                 CompnentDistribution **RMixTheta, // Update increment of mixture model distribution parameter values.
                                 FLOAT                *VW,         // 
                                 CompnentDistribution **VMixTheta, //
                                 FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, ii, j, l, p, q, Error = E_OK;
    
    for (l = 0; l < *c; l++) {
        W[l] += am * ( RW[l] + VW[l]);

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

/// Panic Branislav
        if (W[l] < FLOAT_MIN) {
            switch (merge_) {
            case merge_naive:
                for (j = l; j < *c - 1; j++) {
                    RW[j] = RW[j + 1];
                    
                    VW[j] = VW[j + 1];

                    W[j] = W[j + 1];

                    for (i = 0; i < length_pdf_; i++) {
                        MixTheta[j]->Theta_[0][i] = MixTheta[j + 1]->Theta_[0][i];

                        RMixTheta[j]->Theta_[0][i] = RMixTheta[j + 1]->Theta_[0][i];

                        VMixTheta[j]->Theta_[0][i] = VMixTheta[j + 1]->Theta_[0][i];

                        p = i * length_pdf_ + i;

                        MixTheta[j]->Theta_[1][p] = MixTheta[j + 1]->Theta_[1][p];

                        RMixTheta[j]->Theta_[1][p] = RMixTheta[j + 1]->Theta_[1][p];

                        VMixTheta[j]->Theta_[1][p] = VMixTheta[j + 1]->Theta_[1][p];

                        for (ii = 0; ii < i; ii++) {
                            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                            MixTheta[j]->Theta_[1][p] = MixTheta[j + 1]->Theta_[1][p];

                            RMixTheta[j]->Theta_[1][p] = RMixTheta[j + 1]->Theta_[1][p];

                            VMixTheta[j]->Theta_[1][p] = VMixTheta[j + 1]->Theta_[1][p];

                            MixTheta[j]->Theta_[1][q] = MixTheta[j + 1]->Theta_[1][q];

                            RMixTheta[j]->Theta_[1][q] = RMixTheta[j + 1]->Theta_[1][q];

                            VMixTheta[j]->Theta_[1][q] = VMixTheta[j + 1]->Theta_[1][q];
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
            MixTheta[l]->Theta_[0][i] += am * (RMixTheta[l]->Theta_[0][i] + VMixTheta[l]->Theta_[0][i]);

            p = i * length_pdf_ + i;

            MixTheta[l]->Theta_[1][p] += am * (RMixTheta[l]->Theta_[1][p] + VMixTheta[l]->Theta_[1][p]);

            if (MixTheta[l]->Theta_[1][p] < Eps) {
                W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][p] = Eps;
            }

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                MixTheta[l]->Theta_[1][p] += am * (RMixTheta[l]->Theta_[1][p] + VMixTheta[l]->Theta_[1][p]);

                MixTheta[l]->Theta_[1][q] = MixTheta[l]->Theta_[1][p];
            }
        }

        Error = Cholinvdet(length_pdf_, MixTheta[l]->Theta_[1], MixTheta[l]->Theta_[2], MixTheta[l]->Theta_[3]);

        E_CHECK(Error != E_OK, Error);
S1:;
    }

EEXIT:

    E_RETURN(Error);
} // UpdateMixtureParameters

INT Emmvnorm::UpdateMixParamsSqr(INT                  *c,          // Number of components. 
                                 FLOAT                *W,          // Mixture model weight values.
                                 CompnentDistribution **MixTheta,  // Mixture model distribution parameter values.
                                 FLOAT                *RW,         // Update increment of mixture model weights.
                                 CompnentDistribution **RMixTheta, // Update increment of mixture model distribution parameter values.
                                 FLOAT                *VW,         // 
                                 CompnentDistribution **VMixTheta, //
                                 FLOAT                am)          // Acceleration multiplier for EM algorithm.
{
    INT i, ii, j, l, p, q, Error = E_OK;
    
    for (l = 0; l < *c; l++) {
        W[l] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RW[l] + am * am * VW[l];

        if (W[l] < (FLOAT)0.0) W[l] = (FLOAT)0.0;

/// Panic Branislav
        if (W[l] < FLOAT_MIN) {
            switch (merge_) {
            case merge_naive:
                for (j = l; j < *c - 1; j++) {
                    RW[j] = RW[j + 1];
                    
                    VW[j] = VW[j + 1];

                    W[j] = W[j + 1];

                    for (i = 0; i < length_pdf_; i++) {
                        MixTheta[j]->Theta_[0][i] = MixTheta[j + 1]->Theta_[0][i];

                        RMixTheta[j]->Theta_[0][i] = RMixTheta[j + 1]->Theta_[0][i];

                        VMixTheta[j]->Theta_[0][i] = VMixTheta[j + 1]->Theta_[0][i];

                        p = i * length_pdf_ + i;

                        MixTheta[j]->Theta_[1][p] = MixTheta[j + 1]->Theta_[1][p];

                        RMixTheta[j]->Theta_[1][p] = RMixTheta[j + 1]->Theta_[1][p];

                        VMixTheta[j]->Theta_[1][p] = VMixTheta[j + 1]->Theta_[1][p];

                        for (ii = 0; ii < i; ii++) {
                            p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                            MixTheta[j]->Theta_[1][p] = MixTheta[j + 1]->Theta_[1][p];

                            RMixTheta[j]->Theta_[1][p] = RMixTheta[j + 1]->Theta_[1][p];

                            VMixTheta[j]->Theta_[1][p] = VMixTheta[j + 1]->Theta_[1][p];

                            MixTheta[j]->Theta_[1][q] = MixTheta[j + 1]->Theta_[1][q];

                            RMixTheta[j]->Theta_[1][q] = RMixTheta[j + 1]->Theta_[1][q];

                            VMixTheta[j]->Theta_[1][q] = VMixTheta[j + 1]->Theta_[1][q];
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
            MixTheta[l]->Theta_[0][i] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[0][i] + am * am * VMixTheta[l]->Theta_[0][i];

            p = i * length_pdf_ + i;

            MixTheta[l]->Theta_[1][p] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][p] + am * am * VMixTheta[l]->Theta_[1][p];

            if (MixTheta[l]->Theta_[1][p] < Eps) {
                W[l] = (FLOAT)0.0; MixTheta[l]->Theta_[1][p] = Eps;
            }

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                MixTheta[l]->Theta_[1][p] += ((FLOAT)2.0 * am - (FLOAT)1.0) * RMixTheta[l]->Theta_[1][p] + am * am * VMixTheta[l]->Theta_[1][p];

                MixTheta[l]->Theta_[1][q] = MixTheta[l]->Theta_[1][p];
            }
        }

        Error = Cholinvdet(length_pdf_, MixTheta[l]->Theta_[1], MixTheta[l]->Theta_[2], MixTheta[l]->Theta_[3]);

        E_CHECK(Error != E_OK, Error);
S1:;
    }

EEXIT:

    E_RETURN(Error);
} // UpdateMixtureParameters

INT Emmvnorm::Resids(INT *c, 
                     FLOAT *W, 
                     CompnentDistribution **MixTheta, 
                     FLOAT *RW, 
                     CompnentDistribution **RMixTheta, 
                     FLOAT *VW, 
                     CompnentDistribution **VMixTheta,
                     FLOAT **P) 
{
    INT i, ii, l, p, q, Error = E_OK;
    
    Error = F(*c, W, MixTheta, RW, RMixTheta, P);

    E_CHECK(Error != E_OK, Error);

    Error = UpdateMixParams(c, W, MixTheta, RW, RMixTheta, (FLOAT)1.0);

    E_CHECK(Error != E_OK, Error);

    Error = ExpectationStep(*c, W, MixTheta, P);

    E_CHECK(Error != E_OK, Error);

    Error = F(*c, W, MixTheta, VW, VMixTheta, P);

    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < *c; l++) {
        VW[l] -= RW[l];
        for (i = 0; i < length_pdf_; i++) {
            VMixTheta[l]->Theta_[0][i] -= RMixTheta[l]->Theta_[0][i];

            p = i * length_pdf_ + i;

            VMixTheta[l]->Theta_[1][p] -= RMixTheta[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                VMixTheta[l]->Theta_[1][p] -= RMixTheta[l]->Theta_[1][p];

                VMixTheta[l]->Theta_[1][q] = VMixTheta[l]->Theta_[1][p];
            }   
        }
    }

EEXIT:

    E_RETURN(Error);

}

INT Emmvnorm::OptSLen1(INT *c, 
                       FLOAT *W, 
                       CompnentDistribution **MixTheta, 
                       FLOAT *RW, 
                       CompnentDistribution **RMixTheta, 
                       FLOAT *VW, 
                       CompnentDistribution **VMixTheta,
                       FLOAT **P,
                       FLOAT *am_opt) 
{
    INT i, ii, l, p, q, Error = E_OK;

    FLOAT A = (FLOAT)0.0, B = (FLOAT)0.0;

    Error = Resids(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P);
    
    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < c_; l++) {
        A += RW[l] * VW[l];
        B += VW[l] * VW[l];
        for (i = 0; i < length_pdf_; i++) {
            A += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];
            
            B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

            p = i * length_pdf_ + i;

            A += RMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];
            
            B += VMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                A += RMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];
            
                B += VMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];

                A += RMixTheta[l]->Theta_[1][q] * VMixTheta[l]->Theta_[1][q];
            
                B += VMixTheta[l]->Theta_[1][q] * VMixTheta[l]->Theta_[1][q];
            }  
            
        }
    }

    *am_opt = Max(- A / (B + Eps), (FLOAT)1.0);

EEXIT:

    E_RETURN(Error);

}
INT Emmvnorm::OptSLen2(INT *c, 
                       FLOAT *W, 
                       CompnentDistribution **MixTheta, 
                       FLOAT *RW, 
                       CompnentDistribution **RMixTheta, 
                       FLOAT *VW, 
                       CompnentDistribution **VMixTheta,
                       FLOAT **P,
                       FLOAT *am_opt) 
{
    INT i, ii, l, p, q, Error = E_OK;

    FLOAT A = (FLOAT)0.0, B = (FLOAT)0.0;

    Error = Resids(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P);
    
    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < c_; l++) {
        A += RW[l] * RW[l];
        B += RW[l] * VW[l];
        for (i = 0; i < length_pdf_; i++) {
            A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];
            
            B += RMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

            p = i * length_pdf_ + i;

            A += RMixTheta[l]->Theta_[1][p] * RMixTheta[l]->Theta_[1][p];
            
            B += RMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                A += RMixTheta[l]->Theta_[1][p] * RMixTheta[l]->Theta_[1][p];
            
                B += RMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];

                A += RMixTheta[l]->Theta_[1][q] * RMixTheta[l]->Theta_[1][q];
            
                B += RMixTheta[l]->Theta_[1][q] * VMixTheta[l]->Theta_[1][q];
            }  
            
        }
    }

    *am_opt = Max(- A / (B + Eps), (FLOAT)1.0);

EEXIT:

    E_RETURN(Error);
}
INT Emmvnorm::OptSLen3(INT *c, 
                       FLOAT *W, 
                       CompnentDistribution **MixTheta, 
                       FLOAT *RW, 
                       CompnentDistribution **RMixTheta, 
                       FLOAT *VW, 
                       CompnentDistribution **VMixTheta,
                       FLOAT **P,
                       FLOAT *am_opt) 
{
    INT i, ii, l, p, q, Error = E_OK;

    FLOAT A = (FLOAT)0.0, B = (FLOAT)0.0;

    Error = Resids(c, W, MixTheta, RW, RMixTheta, VW, VMixTheta, P);
    
    E_CHECK(Error != E_OK, Error);

    for (l = 0; l < c_; l++) {
        A += RW[l] * RW[l];
        B += VW[l] * VW[l];
        for (i = 0; i < length_pdf_; i++) {
            A += RMixTheta[l]->Theta_[0][i] * RMixTheta[l]->Theta_[0][i];
            
            B += VMixTheta[l]->Theta_[0][i] * VMixTheta[l]->Theta_[0][i];

            p = i * length_pdf_ + i;

            A += RMixTheta[l]->Theta_[1][p] * RMixTheta[l]->Theta_[1][p];
            
            B += VMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii; q = ii * length_pdf_ + i;

                A += RMixTheta[l]->Theta_[1][p] * RMixTheta[l]->Theta_[1][p];
            
                B += VMixTheta[l]->Theta_[1][p] * VMixTheta[l]->Theta_[1][p];

                A += RMixTheta[l]->Theta_[1][q] * RMixTheta[l]->Theta_[1][q];
            
                B += VMixTheta[l]->Theta_[1][q] * VMixTheta[l]->Theta_[1][q];
            }  
            
        }
    }

    A = sqrt(A);

    B = sqrt(B);

    if (fabs(A) < FLOAT_MIN || fabs(B) < FLOAT_MIN) Error = E_NO_SOLUTION; 

    *am_opt = Max(A / B, (FLOAT)1.0);

EEXIT:

    E_RETURN(Error);

}


// Maximization step of the EM algoritm.

INT Emmvnorm::F(INT c, FLOAT *W, CompnentDistribution **MixTheta, FLOAT *dW, CompnentDistribution **dMixTheta, FLOAT **P) {
    FLOAT *C = NULL, *M = NULL;
    INT   i, ii, j, l, p, q, Error = E_OK;
     
    M = (FLOAT*)malloc(length_pdf_ * sizeof(FLOAT));

    E_CHECK(NULL == M, E_MEM);
    
    C = (FLOAT*)malloc(length_pdf_ * length_pdf_ * sizeof(FLOAT));
    
    E_CHECK(NULL == C, E_MEM);
    
    for (l = 0; l < c; l++) {
        dW[l] = (FLOAT)0.0;

        for (j = 0; j < nr_; j++) {
            dW[l] += Y_[length_pdf_][j] * P[l][j];
        }

        memset(M, 0, length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            for (j = 0; j < nr_; j++) {
                if (Y_[length_pdf_][j] > FLOAT_MIN) {
                    M[i] += Y_[length_pdf_][j] * P[l][j] * Y_[i][j];
                }
            }

            M[i] = M[i] / (dW[l] + FLOAT_MIN);

            dMixTheta[l]->Theta_[0][i] = M[i] - MixTheta[l]->Theta_[0][i];
        }

        memset(C, 0, length_pdf_ * length_pdf_ * sizeof(FLOAT));

        for (i = 0; i < length_pdf_; i++) {
            p = i * length_pdf_ + i;

            for (j = 0; j < nr_; j++) {
                if (Y_[length_pdf_][j] > FLOAT_MIN) {
                    C[p] += Y_[length_pdf_][j] * P[l][j] * (Y_[i][j] - M[i]) * (Y_[i][j] - M[i]);
                }
            }

            dMixTheta[l]->Theta_[1][p] = C[p] / (dW[l] + FLOAT_MIN) - MixTheta[l]->Theta_[1][p];

            for (ii = 0; ii < i; ii++) {
                p = i * length_pdf_ + ii;

                for (j = 0; j < nr_; j++) {
                    if (Y_[length_pdf_][j] > FLOAT_MIN) {
                        C[p] += Y_[length_pdf_][j] * P[l][j] * (Y_[i][j] - M[i]) * (Y_[ii][j] - M[ii]);
                    }
                }

                dMixTheta[l]->Theta_[1][p] = C[p] / (dW[l] + FLOAT_MIN) - MixTheta[l]->Theta_[1][p];

                q = ii * length_pdf_ + i;

                dMixTheta[l]->Theta_[1][q] = dMixTheta[l]->Theta_[1][p];
            }
        }

        dW[l] = dW[l] / n_ - W[l];
    }

EEXIT:
    
    if (C) free(C); 
    
    if (M) free(M);

    n_iter_++;

    E_RETURN(Error);

} 