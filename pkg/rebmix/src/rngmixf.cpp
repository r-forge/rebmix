#include <math.h>
#include <stdio.h>
#include <ctype.h>

#include "base.h"
#include "rngmixf.h"

static int   NDevISet = 0;
static FLOAT NDevVSet = (FLOAT)0.0;
static int   LDevISet = 0;
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
    d_ = 0;
    c_ = 0;
    IniTheta_ = NULL;
    n_ = 0;
    Dataset_ = NULL;
    Y_ = NULL;
    N_ = NULL;
    MixTheta_ = NULL;
} // Rngmix

// Rngmix destructor.

Rngmix::~Rngmix()
{
    int i;

    if (MixTheta_) {
        for (i = 0; i < c_; i++) {
            if (MixTheta_[i]) delete MixTheta_[i];
        }

        delete MixTheta_;
    }

    if (N_) free(N_);

    if (Y_) {
        for (i = 0; i < n_; i++) {
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

// Writes data file.

int Rngmix::WriteDataFile()
{
    int  i, j;
    FILE *fp = NULL;
    int  Error = 0;

    if ((fp = fopen(curr_, "w")) == NULL) {
        Error = 1; goto E0;
    }

    for (i = 0; i < n_; i++) {
        fprintf(fp, "%E", Y_[i][0]); 
        
        for (j = 1; j < d_; j++) fprintf(fp, "\t%E", Y_[i][j]); 
        
        fprintf(fp, "\n");
    }

E0: if (fp) fclose(fp);

    if (Y_) {
        for (i = 0; i < n_; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_); Y_ = NULL;
    }

    return Error;
} // WriteDataFile

// Writes parameters file:

int Rngmix::WriteParameterFile()
{
    char path[FILENAME_MAX];
    char ext[FILENAME_MAX];
    char *pchar = NULL;
    FILE *fp = NULL;
    int  Error = 0;

    strcpy(path, save_); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(path, "%s%s%s", path, "_1", ext);

    if ((fp = fopen(path, "w")) == NULL) {
        Error = 1; goto E0;
    }

    fprintf(fp, "%s\n", "rseed");

    fprintf(fp, "%d\n", IDum_);

E0: if (fp) fclose(fp);

    return Error;
} // WriteParameterFile

int Rngmix::InvComponentDist(CompnentDistribution *CmpDist, FLOAT *Y)
{
    FLOAT C[8];
    FLOAT y, p;
    int   i, j;
    int   Error = 0;

    for (i = 0; i < length_pdf_; i++) {
        switch (CmpDist->pdf_[i]) {
        case pfNormal:
            if (NDevISet == 0) {
                do {
                    C[0] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;
                    C[1] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;

                    C[2] = C[0] * C[0] + C[1] * C[1];
                } while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

                C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * log(C[2]) / C[2]);

                y = C[3] * C[0];

                NDevISet = 1; NDevVSet = C[3] * C[1];
            }
            else {
                y = NDevVSet; NDevISet = 0;
            }

            Y[i] = CmpDist->Theta_[1][i] * y + CmpDist->Theta_[0][i];

            break;
        case pfLognormal:
            if (LDevISet == 0) {
                do {
                    C[0] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;
                    C[1] = (FLOAT)2.0 * Ran1(&IDum_) - (FLOAT)1.0;

                    C[2] = C[0] * C[0] + C[1] * C[1];
                } while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

                C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * log(C[2]) / C[2]);

                y = C[3] * C[0];

                LDevISet = 1; LDevVSet = C[3] * C[1];
            }
            else {
                y = LDevVSet; LDevISet = 0;
            }

            Y[i] = (FLOAT)exp(CmpDist->Theta_[1][i] * y + CmpDist->Theta_[0][i]);

            break;
        case pfWeibull:
            Y[i] = CmpDist->Theta_[0][i] * (FLOAT)exp(log(log((FLOAT)1.0 / Ran1(&IDum_))) / CmpDist->Theta_[1][i]);

            break;
        case pfGamma:
            Error = GammaInv(Ran1(&IDum_), CmpDist->Theta_[0][i], CmpDist->Theta_[1][i], &y);

            if (Error) goto E0;

            Y[i] = y;

            break;
        case pfBinomial:
            if (CmpDist->Theta_[1][i] < (FLOAT)0.5) {
                p = CmpDist->Theta_[1][i];
            }
            else {
                p = (FLOAT)1.0 - CmpDist->Theta_[1][i];
            }

            C[0] = CmpDist->Theta_[0][i] * p;
            if ((int)CmpDist->Theta_[0][i] < 25) {
                Y[i] = (FLOAT)0.0;

                for (j = 0; j < (int)CmpDist->Theta_[0][i]; j++) {
                    if (Ran1(&IDum_) < p) ++Y[i];
                }
            }
            else
            if (C[0] < (FLOAT)1.0) {
                C[1] = (FLOAT)exp(-C[0]); C[2] = (FLOAT)1.0;

                for (j = 0; j < (int)CmpDist->Theta_[0][i]; j++) {
                    C[2] *= Ran1(&IDum_); if (C[2] < C[1]) break;
                }

                if (j > (int)CmpDist->Theta_[0][i]) {
                    Y[i] = CmpDist->Theta_[0][i];
                }
                else {
                    Y[i] = j;
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

                Y[i] = C[6];
            }

            if (p != CmpDist->Theta_[1][i]) {
                Y[i] = CmpDist->Theta_[0][i] - Y[i];
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

            Y[i] = C[0];

            break;
        case pfDirac:
            Y[i] = CmpDist->Theta_[0][i];

            break;
        default:;
        }
    }

E0: return Error;
} // InvComponentDist

// Returns random sample of independent observations.

int Rngmix::RNGMIX()      
{
    int i, j, k;
    int Error = 0;
    
    n_ = 0; for (i = 0; i < c_; i++) n_ += N_[i];
    
    Y_ = (FLOAT**)malloc(n_ * sizeof(FLOAT*));

    Error = NULL == Y_; if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        Y_[i] = (FLOAT*)malloc(d_ * sizeof(FLOAT));

        Error = NULL == Y_[i]; if (Error) goto E0;
    }

    k = 0;

    for (i = 0; i < c_; i++) {
        Trigger_ = 1;

        for (j = 0; j < N_[i]; j++) {
            Error = InvComponentDist(MixTheta_[i], Y_[k]);

            if (Error) goto E0;

            k++;
        }
    }

E0: return Error;
} // RNGMIX

// Runs template file.

int Rngmix::RunTemplateFile(char *file)
{
    int                  i, imin, imax, j, k, isI;
    char                 line[65536], ident[65536], list[65536];
    char                 *pchar = NULL;
    FILE                 *fp = NULL;
    CompnentDistribution **MixTheta = NULL;
    int                  Error = 0;

    if ((fp = fopen(file, "r")) == NULL) {
        Error = 1; goto E0;
    }

    #if (_REBMIXEXE)
    printf("RNGMIX Version 2.7.3\n");
    #endif

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n"); 
        
        pchar = strtok(pchar, "=");

        if (!pchar) goto S0;

        j = 0;

        for (i = 0; i < (int)strlen(pchar); i++) {
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

            if (Error) goto E0;

            for (k = 0; k < o_; k++) {
                curr_ = open_[k]; 

                #if (_REBMIXEXE)
                printf("Dataset = %s\n", curr_);
                #endif

                Error = RNGMIX();

                if (Error) goto E0;

                Error = WriteDataFile();

                if (Error) goto E0;

                IDum_--;
            }
        }
        else
        if (!strcmp(ident, "DATASET")) {
            open_ = (char**)realloc(open_, (o_ + 1) * sizeof(char*));

            Error = NULL == open_; if (Error) goto E0;

            open_[o_] = (char*)malloc((strlen(pchar) + 1) * sizeof(char));

            Error = NULL == open_[o_]; if (Error) goto E0;

            strcpy(open_[o_], pchar); o_++;
        }
        else
        if (!strcmp(ident, "RSEED")) {
            IDum_ = isI = (int)atol(pchar);

            Error = isI >= 0; if (Error) goto E0;
        }
        else
        if (!strcmp(ident, "D")) {
            d_ = isI = (int)atol(pchar);

            Error = isI < 1; if (Error) goto E0;
        }
        else
        if (!strcmp(ident, "LENGTHPDF")) {
            length_pdf_ = isI = (int)atol(pchar);

            Error = isI < 1; if (Error) goto E0;
        } 
        else
        if (!strcmp(ident, "LENGTHTHETA")) {
            i = 0;

            while (pchar) {
                length_theta_ = (int*)realloc(length_theta_, (i + 1) * sizeof(int));

                Error = NULL == length_theta_; if (Error) goto E0;

                length_theta_[i] = isI = (int)atol(pchar);

                Error = isI == 0; if (Error) goto E0;

                pchar = strtok(NULL, "\t"); ++i;
            }

            length_Theta_ = i;

            IniTheta_ = new CompnentDistribution(this);

            Error = NULL == IniTheta_; if (Error) goto E0;

            Error = IniTheta_->Realloc(length_pdf_, length_Theta_, length_theta_);

            if (Error) goto E0;
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
                if (!strcmp(pchar, "BINOMIAL"))
                    IniTheta_->pdf_[i] = pfBinomial;
                else
                if (!strcmp(pchar, "POISSON"))
                    IniTheta_->pdf_[i] = pfPoisson;
                else
                if (!strcmp(pchar, "DIRAC"))
                    IniTheta_->pdf_[i] = pfDirac;
                else
                if (!strcmp(pchar, "MVNORMAL"))
                    IniTheta_->pdf_[i] = pfNormal;
                else {
                    Error = 1; goto E0;
                }

                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((length_pdf_ > 0) && (length_pdf_ != i)) {
                Error = 1; goto E0;
            }
        }
        else
        if (!strcmp(ident, "NTHETA")) {
            N_ = (int*)realloc(N_, (c_ + 1) * sizeof(int));

            Error = NULL == N_; if (Error) goto E0;

            N_[c_] = isI = (int)atol(pchar);

            Error = isI < 1; if (Error) goto E0;

            MixTheta = new CompnentDistribution* [c_ + 1];

            Error = NULL == MixTheta; if (Error) goto E0;

            MixTheta[c_] = new CompnentDistribution(this);

            Error = NULL == MixTheta[c_]; if (Error) goto E0;

            for (i = 0; i < c_; i++) {
                MixTheta[i] = MixTheta_[i];
            }

            if (MixTheta_) delete MixTheta_;

            MixTheta_ = MixTheta;

            Error = MixTheta_[c_]->Realloc(length_pdf_, length_Theta_, length_theta_);

            if (Error) goto E0;

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

            Error = NULL == save_; if (Error) goto E0;

            strcpy(save_, pchar);
        }
    }

E0: if (fp) fclose(fp);

    return Error;
} // RunTemplateFile
