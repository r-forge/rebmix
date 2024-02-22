#include "rngmixf.h"

#include <math.h>

#if (_MAINTAIN_SWITCH)
#include <ctype.h>
#include <stdio.h>
#endif

static INT   NDevISet = 0;
static FLOAT NDevVSet = (FLOAT)0.0;
static INT   LDevISet = 0;
static FLOAT LDevVSet = (FLOAT)0.0;
static FLOAT Bn = -(FLOAT)1.0, Bp = -(FLOAT)1.0, Be, Bg, Bplog, Bpc, Bpclog;
static FLOAT PTheta = -(FLOAT)1.0, Pg, Psq, PalTheta;

// Rngmix constructor.

Rngmix::Rngmix()
{
    curr_ = NULL;
    o_ = 0;
    open_ = NULL;
    save_ = NULL;
    IDum_ = 0;
    c_ = 0;
    IniTheta_ = NULL;
    n_ = 0;
    Y_ = NULL;
    N_ = NULL;
    MixTheta_ = NULL;
    Z_ = NULL;
} // Rngmix

// Rngmix destructor.

Rngmix::~Rngmix()
{
    INT i;

    if (Z_) free(Z_);

    if (MixTheta_) {
        for (i = 0; i < c_; i++) {
            if (MixTheta_[i]) delete MixTheta_[i];
        }

        delete[] MixTheta_;
    }

    if (N_) free(N_);

    if (Y_) {
        for (i = 0; i < length_pdf_; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_);
    }

    if (IniTheta_) delete IniTheta_;

    if (save_) free(save_);

    if (open_) {
        for (i = 0; i < o_; i++) {
            if (open_[i]) free(open_[i]);
        }

        free(open_);
    }
} // ~Rngmix

#if (_MAINTAIN_SWITCH)
// Writes data file.

INT Rngmix::WriteDataFile()
{
    FILE *fp = NULL;
    INT  i, j, Error = EOK;

    fp = fopen(curr_, "w");

    E_CHECK(NULL == fp, ERngmixWriteDataFile);

    for (i = 0; i < n_; i++) {
        fprintf(fp, "%E", Y_[0][i]);

        for (j = 1; j < length_pdf_; j++) fprintf(fp, "\t%E", Y_[j][i]);

        fprintf(fp, "\t%d\n", Z_[i]);
    }

EEXIT:

    if (fp) fclose(fp);

    if (Z_) free(Z_);

    Z_ = NULL;

    if (Y_) {
        for (i = 0; i < length_pdf_; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_); Y_ = NULL;
    }

    E_RETURN(Error);
} // WriteDataFile
#endif

#if (_MAINTAIN_SWITCH)
// Writes parameters file:

INT Rngmix::WriteParameterFile()
{
    char line[65536];
    char path[FILENAME_MAX];
    char ext[FILENAME_MAX];
    char *pchar = NULL;
    FILE *fp = NULL;
    INT  Error = EOK;

    strcpy(path, save_);

    pchar = strrchr(path, '.');

    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }

    sprintf(line, "%s%s%s", path, "_1", ext);

    fp = fopen(line, "w");

    E_CHECK(NULL == fp, ERngmixWriteParameterFile);

    fprintf(fp, "%s\n", "rseed");

    fprintf(fp, "%d\n", IDum_);

EEXIT:

    if (fp) fclose(fp);

    E_RETURN(Error);
} // WriteParameterFile
#endif

INT Rngmix::InvComponentDist(CompnentDistribution *CmpDist, INT j, FLOAT **Y)
{
    FLOAT C[8], y, p;
    INT   i, k, Error = EOK;

    for (i = 0; i < length_pdf_; i++) {
        switch (CmpDist->pdf_[i]) {
        case pfNormal:
            if (NDevISet == 0) {
                do {
                    C[0] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;
                    C[1] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;

                    C[2] = C[0] * C[0] + C[1] * C[1];
                } while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

                C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * (FLOAT)log(C[2]) / C[2]);

                y = C[3] * C[0];

                NDevISet = 1; NDevVSet = C[3] * C[1];
            }
            else {
                y = NDevVSet; NDevISet = 0;
            }

            Y[i][j] = CmpDist->Theta_[1][i] * y + CmpDist->Theta_[0][i];

            break;
        case pfTNormal:
            break;
        case pfLognormal:
            if (LDevISet == 0) {
                do {
                    C[0] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;
                    C[1] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;

                    C[2] = C[0] * C[0] + C[1] * C[1];
                } while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

                C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * (FLOAT)log(C[2]) / C[2]);

                y = C[3] * C[0];

                LDevISet = 1; LDevVSet = C[3] * C[1];
            }
            else {
                y = LDevVSet; LDevISet = 0;
            }

            Y[i][j] = (FLOAT)exp(CmpDist->Theta_[1][i] * y + CmpDist->Theta_[0][i]);

            break;
        case pfWeibull:
            Y[i][j] = CmpDist->Theta_[0][i] * (FLOAT)exp((FLOAT)log((FLOAT)log((FLOAT)1.0 / Ran1(&IDum_))) / CmpDist->Theta_[1][i]);

            break;
        case pfGamma:
            Error = GammaInv(Ran1(&IDum_), CmpDist->Theta_[0][i], CmpDist->Theta_[1][i], &y);

            E_CHECK(Error != EOK, Error);

            Y[i][j] = y;

            break;
        case pfGumbel:
            if (CmpDist->Theta_[2][i] > Eps) {
                Y[i][j] = CmpDist->Theta_[0][i] + CmpDist->Theta_[1][i] * (FLOAT)log((FLOAT)log((FLOAT)1.0 / ((FLOAT)1.0 - Ran1(&IDum_))));
            }
            else {
                Y[i][j] = CmpDist->Theta_[0][i] - CmpDist->Theta_[1][i] * (FLOAT)log((FLOAT)log((FLOAT)1.0 / Ran1(&IDum_)));
            }

            break;
        case pfvonMises:
            CmpDist->Theta_[0][i] -= Pi2 * INT(CmpDist->Theta_[0][i] / Pi2);

            Error = vonMisesInv(Ran1(&IDum_), CmpDist->Theta_[0][i], CmpDist->Theta_[1][i], &Y[i][j]);

            E_CHECK(Error != EOK, Error);

            break;
        case pfBinomial:
            if (CmpDist->Theta_[1][i] < (FLOAT)0.5) {
                p = CmpDist->Theta_[1][i];
            }
            else {
                p = (FLOAT)1.0 - CmpDist->Theta_[1][i];
            }

            C[0] = CmpDist->Theta_[0][i] * p;
            if ((INT)CmpDist->Theta_[0][i] < 25) {
                Y[i][j] = (FLOAT)0.0;

                for (k = 0; k < (INT)CmpDist->Theta_[0][i]; k++) {
                    if (Ran1(&IDum_) < p) ++Y[i][j];
                }
            }
            else
            if (C[0] < (FLOAT)1.0) {
                C[1] = (FLOAT)exp(-C[0]); C[2] = (FLOAT)1.0;

                for (k = 0; k < (INT)CmpDist->Theta_[0][i]; k++) {
                    C[2] *= Ran1(&IDum_); if (C[2] < C[1]) break;
                }

                if (k > (INT)CmpDist->Theta_[0][i]) {
                    Y[i][j] = CmpDist->Theta_[0][i];
                }
                else {
                    Y[i][j] = k;
                }
            }
            else {
                if (CmpDist->Theta_[0][i] != Bn) {
                    Be = CmpDist->Theta_[0][i];
                    Bg = Gammaln(Be + (FLOAT)1.0);
                    Bn = CmpDist->Theta_[0][i];
                }

                if (p != Bp) {
                    Bpc = (FLOAT)1.0 - p;
                    Bplog = (FLOAT)log(p);
                    Bpclog = (FLOAT)log(Bpc);
                    Bp = p;
                }

                C[3] = (FLOAT)sqrt((FLOAT)2.0 * C[0] * Bpc);

                do {
                    do {
                        C[4] = Pi * Ran1(&IDum_);

                        C[5] = (FLOAT)tan(C[4]);

                        C[6] = C[3] * C[5] + C[0];
                    } while ((C[6] < (FLOAT)0.0) || (C[6] >= Be + (FLOAT)1.0));


                    C[6] = (FLOAT)floor(C[6]);

                    C[7] = (FLOAT)1.2 * C[3] * ((FLOAT)1.0 + C[5] * C[5]) *
                           (FLOAT)exp(Bg - Gammaln(C[6] + (FLOAT)1.0) -
                           Gammaln(Be - C[6] + (FLOAT)1.0) +
                           C[6] * Bplog + (Be - C[6]) * Bpclog);

                } while (Ran1(&IDum_) > C[7]);

                Y[i][j] = C[6];
            }

            if (p != CmpDist->Theta_[1][i]) {
                Y[i][j] = CmpDist->Theta_[0][i] - Y[i][j];
            }

            break;
        case pfPoisson:
            if (CmpDist->Theta_[0][i] < (FLOAT)12.0) {
                if (CmpDist->Theta_[0][i] != PTheta) {
                    PTheta = CmpDist->Theta_[0][i];

                    Pg = (FLOAT)exp(-CmpDist->Theta_[0][i]);
                }

                C[0] = -(FLOAT)1.0; C[1] = (FLOAT)1.0;

                do {
                    ++C[0]; C[1] *= Ran1(&IDum_);
                } while (C[1] > Pg);
            }
            else {
                if (CmpDist->Theta_[0][i] != PTheta) {
                    PTheta = CmpDist->Theta_[0][i];

                    Psq = (FLOAT)sqrt((FLOAT)2.0 * CmpDist->Theta_[0][i]);

                    PalTheta = (FLOAT)log(CmpDist->Theta_[0][i]);

                    Pg = CmpDist->Theta_[0][i] * PalTheta - Gammaln(CmpDist->Theta_[0][i] + (FLOAT)1.0);
                }

                do {
                    do {
                        C[2] = (FLOAT)tan(Pi * Ran1(&IDum_));

                        C[0] = Psq * C[2] + CmpDist->Theta_[0][i];
                    } while (C[0] < (FLOAT)0.0);

                    C[0] = (FLOAT)floor(C[0]);

                    C[1] = (FLOAT)0.9 * ((FLOAT)1.0 + C[2] * C[2]) * (FLOAT)exp(C[0] * PalTheta - Gammaln(C[0] + (FLOAT)1.0) - Pg);
                } while (Ran1(&IDum_) > C[1]);
            }

            Y[i][j] = C[0];

            break;
        case pfDirac:
            Y[i][j] = CmpDist->Theta_[0][i];

            break;
        case pfUniform:
            Y[i][j] = CmpDist->Theta_[0][i] + Ran1(&IDum_) * (CmpDist->Theta_[1][i] - CmpDist->Theta_[0][i]);

            break;
        default:;
        }
    }

EEXIT:

    E_RETURN(Error);
} // InvComponentDist

// Returns random sample of independent observations.

INT Rngmix::RNGMIX()
{
    INT i, j, k, Error = EOK;

    n_ = 0; for (i = 0; i < c_; i++) n_ += N_[i];

    Y_ = (FLOAT**)malloc(length_pdf_ * sizeof(FLOAT*));

    E_CHECK(NULL == Y_, ERngmixRNGMIX);

    for (i = 0; i < length_pdf_; i++) {
        Y_[i] = (FLOAT*)malloc(n_ * sizeof(FLOAT));

        E_CHECK(NULL == Y_[i], ERngmixRNGMIX);
    }

    Z_ = (INT*)malloc(n_ * sizeof(INT));

    E_CHECK(NULL == Z_, ERngmixRNGMIX);

    k = 0;

    for (i = 0; i < c_; i++) {
        Trigger_ = 1;

        for (j = 0; j < N_[i]; j++) {
            Z_[k] = i + 1;

            Error = InvComponentDist(MixTheta_[i], k, Y_);

            E_CHECK(Error != EOK, Error);

            k++;
        }
    }

EEXIT:

    E_RETURN(Error);
} // RNGMIX

#if (_MAINTAIN_SWITCH)
// Runs template file.

INT Rngmix::RunTemplateFile(char *file)
{
    char                 line[65536], ident[65536], list[65536];
    char                 *pchar = NULL;
    FILE                 *fp = NULL;
    CompnentDistribution **MixTheta = NULL;
    INT                  i, imin, imax, isI, j, k, Error = EOK;

    fp = fopen(file, "r");

    E_CHECK(NULL == fp, ERngmixRunTemplateFile);

    printf("RNGMIX Version 2.16.0\n");

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n");

        pchar = strtok(pchar, "=");

        if (!pchar) goto S0;

        j = 0;

        for (i = 0; i < (INT)strlen(pchar); i++) {
            if (pchar[i] != ' ') {
                ident[j] = (char)toupper(pchar[i]); j++;
            }
        }

        ident[j] = '\0';

        j = 0; list[j] = '\0'; imin = 0; imax = 0;

        while((pchar = strtok(NULL, ",")) != NULL) {
            if (!strcmp(ident, "DATASET") || !strcmp(ident, "SAVE")) {
                for (i = 0; pchar[i] != '\0'; i++) {
                    list[j] = pchar[i]; j++;

                    if (pchar[i] != ' ') {
                        imax = i; if (!imin) imin = i;
                    }
                }

                j = imax + 1 - imin;

                for (i = 0; i < j; i++) {
                    list[i] = list[imin + i];
                }
            }
            else {
                for (i = 0; pchar[i] != '\0'; i++) if (pchar[i] != '[' && pchar[i] != ']') {
                    if (pchar[i] == ' ') {
                    }
                    else {
                        list[j] = (char)toupper(pchar[i]); j++;
                    }
                }
            }

            list[j] = '\t'; j++;
        }

        if (!j) goto S0; else list[j - 1] = '\0';

        pchar = strtok(list, "\t");

        if (!strcmp(ident, "RUN")) {
            Error = WriteParameterFile();

            E_CHECK(Error != EOK, Error);

            for (k = 0; k < o_; k++) {
                curr_ = open_[k];

                printf("Dataset = %s\n", curr_);

                Error = RNGMIX();

                E_CHECK(Error != EOK, Error);

                Error = WriteDataFile();

                E_CHECK(Error != EOK, Error);

                IDum_--;
            }
        }
        else
        if (!strcmp(ident, "DATASET")) {
            open_ = (char**)realloc(open_, (o_ + 1) * sizeof(char*));

            E_CHECK(NULL == open_, ERngmixRunTemplateFile);

            open_[o_] = (char*)malloc((strlen(pchar) + 1) * sizeof(char));

            E_CHECK(NULL == open_[o_], ERngmixRunTemplateFile);

            strcpy(open_[o_], pchar); o_++;
        }
        else
        if (!strcmp(ident, "RSEED")) {
            IDum_ = isI = (INT)atol(pchar);

            E_CHECK(isI >= 0, ERngmixRunTemplateFile);
        }
        else
        if (!strcmp(ident, "LENGTHPDF")) {
            length_pdf_ = isI = (INT)atol(pchar);

            E_CHECK(isI < 1, ERngmixRunTemplateFile);
        }
        else
        if (!strcmp(ident, "LENGTHTHETA")) {
            i = 0;

            while (pchar) {
                length_theta_ = (INT*)realloc(length_theta_, (i + 1) * sizeof(INT));

                E_CHECK(NULL == length_theta_, ERngmixRunTemplateFile);

                length_theta_[i] = isI = (INT)atol(pchar);

                E_CHECK(isI == 0, ERngmixRunTemplateFile);

                pchar = strtok(NULL, "\t"); ++i;
            }

            length_Theta_ = i;

            IniTheta_ = new CompnentDistribution(this);

            E_CHECK(NULL == IniTheta_, ERngmixRunTemplateFile);

            Error = IniTheta_->Realloc(length_pdf_, length_Theta_, length_theta_);

            E_CHECK(Error != EOK, Error);
        }
        else
        if (!strcmp(ident, "PDF")) {
            i = 0;

            while (pchar) {
                if (!strcmp(pchar, "NORMAL"))
                    IniTheta_->pdf_[i] = pfNormal;
                else
                if (!strcmp(pchar, "LOGNORMAL"))
                    IniTheta_->pdf_[i] = pfLognormal;
                else
                if (!strcmp(pchar, "WEIBULL"))
                    IniTheta_->pdf_[i] = pfWeibull;
                else
                if (!strcmp(pchar, "GAMMA"))
                    IniTheta_->pdf_[i] = pfGamma;
                else
                if (!strcmp(pchar, "GUMBEL"))
                    IniTheta_->pdf_[i] = pfGumbel;
                else
                if (!strcmp(pchar, "VONMISES"))
                    IniTheta_->pdf_[i] = pfvonMises;
                else
                if (!strcmp(pchar, "BINOMIAL"))
                    IniTheta_->pdf_[i] = pfBinomial;
                else
                if (!strcmp(pchar, "POISSON"))
                    IniTheta_->pdf_[i] = pfPoisson;
                else
                if (!strcmp(pchar, "DIRAC"))
                    IniTheta_->pdf_[i] = pfDirac;
                else
                if (!strcmp(pchar, "UNIFORM"))
                    IniTheta_->pdf_[i] = pfUniform;
                else {
                    E_CHECK(1, ERngmixRunTemplateFile);
                }

                pchar = strtok(NULL, "\t"); ++i;
            }

            E_CHECK((length_pdf_ > 0) && (length_pdf_ != i), ERngmixRunTemplateFile);
        }
        else
        if (!strcmp(ident, "NTHETA")) {
            N_ = (INT*)realloc(N_, (c_ + 1) * sizeof(INT));

            E_CHECK(NULL == N_, ERngmixRunTemplateFile);

            N_[c_] = isI = (INT)atol(pchar);

            E_CHECK(isI < 1, ERngmixRunTemplateFile);

            MixTheta = new CompnentDistribution* [(unsigned INT)(c_ + 1)];

            E_CHECK(NULL == MixTheta, ERngmixRunTemplateFile);

            MixTheta[c_] = new CompnentDistribution(this);

            E_CHECK(NULL == MixTheta[c_], ERngmixRunTemplateFile);

            for (i = 0; i < c_; i++) {
                MixTheta[i] = MixTheta_[i];
            }

            if (MixTheta_) delete[] MixTheta_;

            MixTheta_ = MixTheta;

            Error = MixTheta_[c_]->Realloc(length_pdf_, length_Theta_, length_theta_);

            E_CHECK(Error != EOK, Error);

            for (i = 0; i < length_pdf_; i++) {
                MixTheta_[c_]->pdf_[i] = IniTheta_->pdf_[i];
            }

            for (i = 0; i < length_Theta_; i++) if (IniTheta_->Theta_[i]) {
                for (j = 0; j < length_theta_[i]; j++) {
                    pchar = strtok(NULL, "\t");

                    MixTheta_[c_]->Theta_[i][j] = (FLOAT)atof(pchar);
                }
            }

            c_++;
        }
        else
        if (!strcmp(ident, "SAVE")) {
            save_ = (char*)realloc(save_, (strlen(pchar) + 1) * sizeof(char));

            E_CHECK(NULL == save_, ERngmixRunTemplateFile);

            strcpy(save_, pchar);
        }
    }

EEXIT:

    if (fp) fclose(fp);

    E_RETURN(Error);
} // RunTemplateFile
#endif
