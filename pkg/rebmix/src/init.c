#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void RRNGMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RREBMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
	                void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
					void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
					void *, void *, void *, void *);

extern void RdensKNearestNeighbourXY(void *, void *, void *, void *, void *, void *, void *, void *);
extern void RdensParzenWindowXY(void *, void *, void *, void *, void *, void *, void *);
extern void RdensHistogramXY(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RdensKNearestNeighbourX(void *, void *, void *, void *, void *, void *);
extern void RdensParzenWindowX(void *, void *, void *, void *, void *);
extern void RdensHistogramX(void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLSMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLRMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RPreprocessingKNNMIX(void *, void *, void *, void *, void *, void *, void *);
extern void RPreprocessingPWMIX(void *, void *, void *, void *, void *, void *);
extern void RPreprocessingHMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RInformationCriterionKNNMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionPWMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionHMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCombineComponentsKNNMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RCombineComponentsPWMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RCombineComponentsHMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RRNGMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RREBMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
	                   void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
					   void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
					   void *, void *, void *, void *);

extern void RCLSMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLRMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RPreprocessingKNNMVNORM(void *, void *, void *, void *, void *, void *, void *);
extern void RPreprocessingPWMVNORM(void *, void *, void *, void *, void *, void *);
extern void RPreprocessingHMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RInformationCriterionKNNMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionPWMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionHMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCombineComponentsKNNMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RCombineComponentsPWMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RCombineComponentsHMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RvonMisesPdf(void *, void *, void *, void *, void *);
extern void RvonMisesCdf(void *, void *, void *, void *, void *);

static const R_CMethodDef CMethods[] = {
	{"RRNGMIX", (DL_FUNC) &RRNGMIX, 13},
	{"RREBMIX", (DL_FUNC) &RREBMIX, 43},
	{"RdensKNearestNeighbourXY", (DL_FUNC) &RdensKNearestNeighbourXY, 8},
	{"RdensParzenWindowXY", (DL_FUNC) &RdensParzenWindowXY, 7},
	{"RdensHistogramXY", (DL_FUNC) &RdensHistogramXY, 12},
	{"RdensKNearestNeighbourX", (DL_FUNC) &RdensKNearestNeighbourX, 6},
	{"RdensParzenWindowX", (DL_FUNC) &RdensParzenWindowX, 5},
	{"RdensHistogramX", (DL_FUNC) &RdensHistogramX, 8},
	{"RCLSMIX", (DL_FUNC) &RCLSMIX, 13},
	{"RCLRMIX", (DL_FUNC) &RCLRMIX, 10},
	{"RPreprocessingKNNMIX", (DL_FUNC) &RPreprocessingKNNMIX, 7},
	{"RPreprocessingPWMIX", (DL_FUNC) &RPreprocessingPWMIX, 6},
	{"RPreprocessingHMIX", (DL_FUNC) &RPreprocessingHMIX, 10},
	{"RInformationCriterionKNNMIX", (DL_FUNC) &RInformationCriterionKNNMIX, 17},
	{"RInformationCriterionPWMIX", (DL_FUNC) &RInformationCriterionPWMIX, 16},
	{"RInformationCriterionHMIX", (DL_FUNC) &RInformationCriterionHMIX, 18},
	{"RCombineComponentsKNNMIX", (DL_FUNC) &RCombineComponentsKNNMIX, 16},
	{"RCombineComponentsPWMIX", (DL_FUNC) &RCombineComponentsPWMIX, 15},
	{"RCombineComponentsHMIX", (DL_FUNC) &RCombineComponentsHMIX, 17},
	{"RRNGMVNORM", (DL_FUNC) &RRNGMVNORM, 12},
	{"RREBMVNORM", (DL_FUNC) &RREBMVNORM, 43},
	{"RCLSMVNORM", (DL_FUNC) &RCLSMVNORM, 13},
	{"RCLRMVNORM", (DL_FUNC) &RCLRMVNORM, 10},
	{"RPreprocessingKNNMVNORM", (DL_FUNC) &RPreprocessingKNNMVNORM, 7},
	{"RPreprocessingPWMVNORM", (DL_FUNC) &RPreprocessingPWMVNORM, 6},
	{"RPreprocessingHMVNORM", (DL_FUNC) &RPreprocessingHMVNORM, 10},
	{"RInformationCriterionKNNMVNORM", (DL_FUNC) &RInformationCriterionKNNMVNORM, 17},
	{"RInformationCriterionPWMVNORM", (DL_FUNC) &RInformationCriterionPWMVNORM, 16},
	{"RInformationCriterionHMVNORM", (DL_FUNC) &RInformationCriterionHMVNORM, 18},
	{"RCombineComponentsKNNMVNORM", (DL_FUNC) &RCombineComponentsKNNMVNORM, 16},
	{"RCombineComponentsPWMVNORM", (DL_FUNC) &RCombineComponentsPWMVNORM, 15},
	{"RCombineComponentsHMVNORM", (DL_FUNC) &RCombineComponentsHMVNORM, 17},
	{"RvonMisesPdf", (DL_FUNC) &RvonMisesPdf, 5},
	{"RvonMisesCdf", (DL_FUNC) &RvonMisesCdf, 5},
    {NULL, NULL, 0}
};

void R_init_rebmix(DllInfo *dll)
{
	R_registerRoutines(dll, CMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}