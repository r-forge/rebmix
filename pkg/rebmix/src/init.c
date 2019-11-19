#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void RRNGMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RREBMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                    void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RdensKNearestNeighbourXY(void *, void *, void *, void *, void *, void *, void *, void *);
extern void RdensKDEXY(void *, void *, void *, void *, void *, void *, void *);
extern void RdensHistogramXY(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RdensKNearestNeighbourX(void *, void *, void *, void *, void *, void *);
extern void RdensKDEX(void *, void *, void *, void *, void *);
extern void RdensHistogramX(void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLSMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLRMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RPreprocessingKNNMIX(void *, void *, void *, void *, void *, void *, void *);
extern void RPreprocessingKDEMIX(void *, void *, void *, void *, void *, void *);
extern void RPreprocessingHMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RInformationCriterionKNNMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionKDEMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionHMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCombineComponentsMIX(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RRNGMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RREBMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,
                       void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLSMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCLRMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RPreprocessingKNNMVNORM(void *, void *, void *, void *, void *, void *, void *);
extern void RPreprocessingKDEMVNORM(void *, void *, void *, void *, void *, void *);
extern void RPreprocessingHMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RInformationCriterionKNNMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionKDEMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RInformationCriterionHMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RCombineComponentsMVNORM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern void RvonMisesPdf(void *, void *, void *, void *, void *);
extern void RvonMisesCdf(void *, void *, void *, void *, void *);

static const R_CMethodDef CMethods[] = {
	{"RRNGMIX", (DL_FUNC) &RRNGMIX, 13},
	{"RREBMIX", (DL_FUNC) &RREBMIX, 52},
	{"RdensKNearestNeighbourXY", (DL_FUNC) &RdensKNearestNeighbourXY, 8},
	{"RdensKDEXY", (DL_FUNC) &RdensKDEXY, 7},
	{"RdensHistogramXY", (DL_FUNC) &RdensHistogramXY, 12},
	{"RdensKNearestNeighbourX", (DL_FUNC) &RdensKNearestNeighbourX, 6},
	{"RdensKDEX", (DL_FUNC) &RdensKDEX, 5},
	{"RdensHistogramX", (DL_FUNC) &RdensHistogramX, 8},
	{"RCLSMIX", (DL_FUNC) &RCLSMIX, 13},
	{"RCLRMIX", (DL_FUNC) &RCLRMIX, 10},
	{"RPreprocessingKNNMIX", (DL_FUNC) &RPreprocessingKNNMIX, 7},
	{"RPreprocessingKDEMIX", (DL_FUNC) &RPreprocessingKDEMIX, 6},
	{"RPreprocessingHMIX", (DL_FUNC) &RPreprocessingHMIX, 10},
	{"RInformationCriterionKNNMIX", (DL_FUNC) &RInformationCriterionKNNMIX, 17},
	{"RInformationCriterionKDEMIX", (DL_FUNC) &RInformationCriterionKDEMIX, 16},
	{"RInformationCriterionHMIX", (DL_FUNC) &RInformationCriterionHMIX, 18},
	{"RCombineComponentsMIX", (DL_FUNC) &RCombineComponentsMIX, 15},
	{"RRNGMVNORM", (DL_FUNC) &RRNGMVNORM, 12},
	{"RREBMVNORM", (DL_FUNC) &RREBMVNORM, 52},
	{"RCLSMVNORM", (DL_FUNC) &RCLSMVNORM, 13},
	{"RCLRMVNORM", (DL_FUNC) &RCLRMVNORM, 10},
	{"RPreprocessingKNNMVNORM", (DL_FUNC) &RPreprocessingKNNMVNORM, 7},
	{"RPreprocessingKDEMVNORM", (DL_FUNC) &RPreprocessingKDEMVNORM, 6},
	{"RPreprocessingHMVNORM", (DL_FUNC) &RPreprocessingHMVNORM, 10},
	{"RInformationCriterionKNNMVNORM", (DL_FUNC) &RInformationCriterionKNNMVNORM, 17},
	{"RInformationCriterionKDEMVNORM", (DL_FUNC) &RInformationCriterionKDEMVNORM, 16},
	{"RInformationCriterionHMVNORM", (DL_FUNC) &RInformationCriterionHMVNORM, 18},
	{"RCombineComponentsMVNORM", (DL_FUNC) &RCombineComponentsMVNORM, 15},
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
