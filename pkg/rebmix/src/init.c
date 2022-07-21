#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void RRNGMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RREBMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RdensKNearestNeighbourXY(void *, void *, void *, void *, void *, void *, void *, void *);
extern void RdensKDEXY(void *, void *, void *, void *, void *, void *, void *);
extern void RdensHistogramXY(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RdensKXY(void *, void *, void *, void *, void *, void *, void *, void *);

extern void RdensKNearestNeighbourX(void *, void *, void *, void *, void *, void *);
extern void RdensKDEX(void *, void *, void *, void *, void *);
extern void RdensHistogramX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RdensKX(void *, void *, void *, void *, void *, void *);

extern void RCLSMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLRMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RPreprocessingKNNMIX(void *, void *, void *, void *, void *, void *, void *);
extern void RPreprocessingKDEMIX(void *, void *, void *, void *, void *, void *);
extern void RPreprocessingHMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RPreprocessingKMIX(void *, void *, void *, void *, void *);

extern void RInformationCriterionKNNMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionKDEMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionHMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionKMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionMIX(void *, void *, void *, void *, void *, void *,  void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCombineComponentsMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/// Panic Branislav
extern void REMMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *);
/// End

extern void RRNGMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RREBMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLSMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLRMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RPreprocessingKNNMVNORM(void *, void *, void *, void *, void *, void *, void *);
extern void RPreprocessingKDEMVNORM(void *, void *, void *, void *, void *, void *);
extern void RPreprocessingHMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RInformationCriterionKNNMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionKDEMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionHMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionKMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCombineComponentsMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/// Panic Branislav
extern void REMMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *);
/// End

extern void RvonMisesPdf(void *, void *, void *, void *, void *);
extern void RvonMisesCdf(void *, void *, void *, void *, void *);

extern void RGumbelPdf(void *, void *, void *, void *, void *, void *);
extern void RGumbelCdf(void *, void *, void *, void *, void *, void *);

extern void Roptbins(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Rbins(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RTvtNormalPdf(void *, void *, void *, void *, void *, void *);
extern void RMvtNormalPdf(void *, void *, void *, void *, void *, void *);

extern void Rfhistogram(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Rchistogram(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CMethods[] = {
    {"RRNGMIX", (DL_FUNC) &RRNGMIX, 13},
    {"RREBMIX", (DL_FUNC) &RREBMIX, 55},
    {"RdensKNearestNeighbourXY", (DL_FUNC) &RdensKNearestNeighbourXY, 8},
    {"RdensKDEXY", (DL_FUNC) &RdensKDEXY, 7},
    {"RdensHistogramXY", (DL_FUNC) &RdensHistogramXY, 16},
    {"RdensKXY", (DL_FUNC)&RdensKXY, 8},
    {"RdensKNearestNeighbourX", (DL_FUNC) &RdensKNearestNeighbourX, 6},
    {"RdensKDEX", (DL_FUNC) &RdensKDEX, 5},
    {"RdensHistogramX", (DL_FUNC) &RdensHistogramX, 10},
    {"RdensKX", (DL_FUNC)&RdensKX, 6},
    {"RCLSMIX", (DL_FUNC) &RCLSMIX, 14},
    {"RCLRMIX", (DL_FUNC) &RCLRMIX, 11},
    {"RPreprocessingKNNMIX", (DL_FUNC) &RPreprocessingKNNMIX, 7},
    {"RPreprocessingKDEMIX", (DL_FUNC) &RPreprocessingKDEMIX, 6},
    {"RPreprocessingHMIX", (DL_FUNC) &RPreprocessingHMIX, 10},
    {"RInformationCriterionKNNMIX", (DL_FUNC) &RInformationCriterionKNNMIX, 20},
    {"RInformationCriterionKDEMIX", (DL_FUNC) &RInformationCriterionKDEMIX, 19},
    {"RInformationCriterionHMIX", (DL_FUNC) &RInformationCriterionHMIX, 23},
    {"RInformationCriterionKMIX", (DL_FUNC)&RInformationCriterionKMIX, 19},
    {"RInformationCriterionMIX", (DL_FUNC)&RInformationCriterionMIX, 18},
    {"RCombineComponentsMIX", (DL_FUNC) &RCombineComponentsMIX, 17},
    {"RRNGMVNORM", (DL_FUNC) &RRNGMVNORM, 12},
    {"RREBMVNORM", (DL_FUNC) &RREBMVNORM, 54},
    {"RCLSMVNORM", (DL_FUNC) &RCLSMVNORM, 13},
    {"RCLRMVNORM", (DL_FUNC) &RCLRMVNORM, 10},
    {"RPreprocessingKNNMVNORM", (DL_FUNC) &RPreprocessingKNNMVNORM, 7},
    {"RPreprocessingKDEMVNORM", (DL_FUNC) &RPreprocessingKDEMVNORM, 6},
    {"RPreprocessingHMVNORM", (DL_FUNC) &RPreprocessingHMVNORM, 10},
    {"RPreprocessingKMIX", (DL_FUNC)&RPreprocessingKMIX, 5},
    {"RInformationCriterionKNNMVNORM", (DL_FUNC) &RInformationCriterionKNNMVNORM, 17},
    {"RInformationCriterionKDEMVNORM", (DL_FUNC) &RInformationCriterionKDEMVNORM, 16},
    {"RInformationCriterionHMVNORM", (DL_FUNC) &RInformationCriterionHMVNORM, 20},
    {"RInformationCriterionKMVNORM", (DL_FUNC)&RInformationCriterionKMVNORM, 16},
    {"RInformationCriterionMVNORM", (DL_FUNC)&RInformationCriterionMVNORM, 15},
    {"RCombineComponentsMVNORM", (DL_FUNC) &RCombineComponentsMVNORM, 17},
    {"RvonMisesPdf", (DL_FUNC) &RvonMisesPdf, 5},
    {"RvonMisesCdf", (DL_FUNC) &RvonMisesCdf, 5},
    {"RGumbelPdf", (DL_FUNC) &RGumbelPdf, 6},
    {"RGumbelCdf", (DL_FUNC) &RGumbelCdf, 6},
    {"Roptbins", (DL_FUNC) &Roptbins, 12},
    {"Rbins", (DL_FUNC) &Rbins, 11},
    {"RTvtNormalPdf", (DL_FUNC) &RTvtNormalPdf, 6},
    {"RMvtNormalPdf", (DL_FUNC) &RMvtNormalPdf, 6},
/// Panic Branislav
    {"REMMIX", (DL_FUNC)&REMMIX, 21},
    {"REMMVNORM", (DL_FUNC)&REMMVNORM, 20},
/// End
    {"Rfhistogram", (DL_FUNC)&Rfhistogram, 10},
    {"Rchistogram", (DL_FUNC)&Rchistogram, 10},
    {NULL, NULL, 0}
};

void R_init_rebmix(DllInfo *dll)
{
    R_registerRoutines(dll, CMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
