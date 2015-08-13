#include <math.h>
#include <stdio.h>
#include <ctype.h>

#include "base.h"
#include "rngmixf.h"

int   NDevISet = 0;
FLOAT NDevVSet = (FLOAT)0.0;
int   LDevISet = 0;
FLOAT LDevVSet = (FLOAT)0.0;
FLOAT Bn = -(FLOAT)1.0, Bp = -(FLOAT)1.0, Be, Bg, Bplog, Bpc, Bpclog;
FLOAT PTheta = -(FLOAT)1.0, Pg, Psq, PalTheta;
long  IY = 0;
long  IV[NTAB];

// Minimal random number generator of Park and Miller with Bays-Durham shuffle and added
// safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
// values). Call with IDum a negative integer to initialize; thereafter do not alter IDum
// between successive deviates in a sequence.RNMX should approximate the largest floating
// value that is less than 1. See http://www.nrbook.com/a/bookcpdf/c7-1.pdf

FLOAT Ran1(int *IDum)
{
    int   j, k;
    FLOAT Tmp;

    if (*IDum <= 0 || !IY) {
        *IDum = (-(*IDum) < 1) ? 1 : -(*IDum);

        for (j = NTAB + 7; j >= 0; j--) {
            k = *IDum / IQ;

            *IDum = IA * (*IDum - k * IQ) - IR * k;

            if (*IDum < 0) *IDum += IM;

            if (j < NTAB) IV[j] = *IDum;
        }

        IY = IV[0];
    }

    k = *IDum / IQ; *IDum = IA * (*IDum - k * IQ) - IR * k;

    if (*IDum < 0) *IDum += IM;

    j = IY / NDIV; IY = IV[j]; IV[j] = *IDum;

    if ((Tmp = AM * IY) > RNMX) return (RNMX); else return (Tmp);
} // Ran1

// Constructor.

Rngmix::Rngmix()
{
    curr_ = NULL;
    o_ = 0;
    open_ = NULL;
    save_ = NULL;
    IDum_ = 0;
    d_ = 0;
    c_ = 0;
    n_ = 0;
    Dataset_ = NULL;
    Y_ = NULL;
    N_ = NULL;
    Theta_ = NULL;
} // Rngmix

// Rebmix destructor.

Rngmix::~Rngmix()
{
    int i;

    if (Theta_) {
        for (i = 0; i < c_; i++) {
            if (Theta_[i].Theta2) free(Theta_[i].Theta2);
            if (Theta_[i].Theta1) free(Theta_[i].Theta1);
            if (Theta_[i].pdf) free(Theta_[i].pdf);
        }

        free(Theta_);
    }

    if (N_) free(N_);

    if (Y_) {
        for (i = 0; i < n_; i++) {
            if (Y_[i]) free(Y_[i]);
        }

        free(Y_);
    }

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

    return (Error);
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

    return (Error);
} // WriteParameterFile

int Rngmix::InvComponentDist(ComponentDistributionType *CmpTheta, FLOAT *Y)
{
    FLOAT C[8];
    FLOAT y, p;
    int   i, j;
    int   Error = 0;

    for (i = 0; i < d_; i++) {
        switch (CmpTheta->pdf[i]) {
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

            Y[i] = CmpTheta->Theta2[i] * y + CmpTheta->Theta1[i];

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

            Y[i] = (FLOAT)exp(CmpTheta->Theta2[i] * y + CmpTheta->Theta1[i]);

            break;
        case pfWeibull:
            Y[i] = CmpTheta->Theta1[i] * (FLOAT)exp(log(log((FLOAT)1.0 / Ran1(&IDum_))) / CmpTheta->Theta2[i]);

            break;
        case pfGamma:
            Error = GammaInv(Ran1(&IDum_), CmpTheta->Theta1[i], CmpTheta->Theta2[i], &y);

            if (Error) goto E0;

            Y[i] = y;

            break;
        case pfBinomial:
            if (CmpTheta->Theta2[i] < (FLOAT)0.5) {
                p = CmpTheta->Theta2[i];
            }
            else {
                p = (FLOAT)1.0 - CmpTheta->Theta2[i];
            }

            C[0] = CmpTheta->Theta1[i] * p;
            if ((int)CmpTheta->Theta1[i] < 25) {
                Y[i] = (FLOAT)0.0;

                for (j = 0; j < (int)CmpTheta->Theta1[i]; j++) {
                    if (Ran1(&IDum_) < p) ++Y[i];
                }
            }
            else
            if (C[0] < (FLOAT)1.0) {
                C[1] = (FLOAT)exp(-C[0]); C[2] = (FLOAT)1.0;

                for (j = 0; j < (int)CmpTheta->Theta1[i]; j++) {
                    C[2] *= Ran1(&IDum_); if (C[2] < C[1]) break;
                }

                if (j > (int)CmpTheta->Theta1[i]) {
                    Y[i] = CmpTheta->Theta1[i];
                }
                else {
                    Y[i] = j;
                }
            }
            else {
                if (CmpTheta->Theta1[i] != Bn) {
                    Be = CmpTheta->Theta1[i];
                    Bg = Gammaln(Be + (FLOAT)1.0);
                    Bn = CmpTheta->Theta1[i];
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

            if (p != CmpTheta->Theta2[i]) {
                Y[i] = CmpTheta->Theta1[i] - Y[i];
            }

            break;
        case pfPoisson:
            if (CmpTheta->Theta1[i] < (FLOAT)12.0) {
                if (CmpTheta->Theta1[i] != PTheta) {
                    PTheta = CmpTheta->Theta1[i];

                    Pg = (FLOAT)exp(-CmpTheta->Theta1[i]);
                }

                C[0] = -(FLOAT)1.0; C[1] = (FLOAT)1.0;

                do {
                    ++C[0]; C[1] *= Ran1(&IDum_);
                } while (C[1] > Pg);
            }
            else {
                if (CmpTheta->Theta1[i] != PTheta) {
                    PTheta = CmpTheta->Theta1[i];

                    Psq = (FLOAT)sqrt((FLOAT)2.0 * CmpTheta->Theta1[i]);

                    PalTheta = (FLOAT)log(CmpTheta->Theta1[i]);

                    Pg = CmpTheta->Theta1[i] * PalTheta - Gammaln(CmpTheta->Theta1[i] + (FLOAT)1.0);
                }

                do {
                    do {
                        C[2] = (FLOAT)tan(Pi * Ran1(&IDum_));

                        C[0] = Psq * C[2] + CmpTheta->Theta1[i];
                    } while (C[0] < (FLOAT)0.0);

                    C[0] = (FLOAT)floor(C[0]);

                    C[1] = (FLOAT)0.9 * ((FLOAT)1.0 + C[2] * C[2]) * (FLOAT)exp(C[0] * PalTheta - Gammaln(C[0] + (FLOAT)1.0) - Pg);
                } while (Ran1(&IDum_) > C[1]);
            }

            Y[i] = C[0];

            break;
        case pfDirac:
            Y[i] = CmpTheta->Theta1[i];
        }
    }

E0: return (Error);
} //InvComponentDist

// Returns random sample of independent observations.

int Rngmix::RNGMIX()      
{
    int i, j, l;
    int Error = 0;
    
    n_ = 0; for (i = 0; i < c_; i++) n_ += N_[i];
    
    Y_ = (FLOAT**)malloc(n_ * sizeof(FLOAT*));

    Error = NULL == Y_; if (Error) goto E0;

    for (i = 0; i < n_; i++) {
        Y_[i] = (FLOAT*)malloc(d_ * sizeof(FLOAT));

        Error = NULL == Y_[i]; if (Error) goto E0;
    }

    l = 0;
    for (i = 0; i < c_; i++) {
        for (j = 0; j < N_[i]; j++) {
            Error = InvComponentDist(&Theta_[i], Y_[l]);

            if (Error) goto E0;

            l++;
        }
    }

E0: return (Error);
} // RNGMIX

// Runs template file.

int Rngmix::RunTemplateFile(char *file)
{
    int   i, imin, imax, j, k, isI;
    FLOAT isF;
    char  line[65536], ident[65536], list[65536];
    char  *pchar = NULL;
    FILE  *fp = NULL;
    int   Error = 0;

    if ((fp = fopen(file, "r")) == NULL) {
        Error = 1; goto E0;
    }

    #if (_REBMIXEXE)
    printf("RNGMIX Version 2.7.2\n");
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
        } else
        if (!strcmp(ident, "NTHETA")) {
            N_ = (int*)realloc(N_, (c_ + 1) * sizeof(int));

            Error = NULL == N_; if (Error) goto E0;

            N_[c_] = isI = (int)atol(pchar);

            Error = isI < 1; if (Error) goto E0;

            pchar = strtok(NULL, "\t"); 

            while (pchar) {
                Theta_ = (ComponentDistributionType*)realloc(Theta_, (c_ + 1) * sizeof(ComponentDistributionType));

                Error = NULL == Theta_; if (Error) goto E0;

                d_ = 0; memset(&Theta_[c_], 0, sizeof(ComponentDistributionType));

                while (pchar) {
                    Theta_[c_].pdf = (ParametricFamilyType_e*)realloc(Theta_[c_].pdf, (d_ + 1) * sizeof(ParametricFamilyType_e));

                    Error = NULL == Theta_[c_].pdf; if (Error) goto E0;

                    Theta_[c_].Theta1 = (FLOAT*)realloc(Theta_[c_].Theta1, (d_ + 1) * sizeof(FLOAT));

                    Error = NULL == Theta_[c_].Theta1; if (Error) goto E0;

                    Theta_[c_].Theta2 = (FLOAT*)realloc(Theta_[c_].Theta2, (d_ + 1) * sizeof(FLOAT));

                    Error = NULL == Theta_[c_].Theta2; if (Error) goto E0;

                    if (!strcmp(pchar, "NORMAL")) {
                        Theta_[c_].pdf[d_] = pfNormal;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta1[d_] = (FLOAT)atof(pchar);

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta2[d_] = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "LOGNORMAL")) {
                        Theta_[c_].pdf[d_] = pfLognormal;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta1[d_] = (FLOAT)atof(pchar);

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta2[d_] = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "WEIBULL")) {
                        Theta_[c_].pdf[d_] = pfWeibull;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta1[d_] = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta2[d_] = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "GAMMA")) {
                        Theta_[c_].pdf[d_] = pfGamma;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta1[d_] = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta2[d_] = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "BINOMIAL")) {
                        Theta_[c_].pdf[d_] = pfBinomial;

                        pchar = strtok(NULL, "\t");

                        Theta_[c_].Theta1[d_] = isF = (FLOAT)floor(atof(pchar) + (FLOAT)0.5);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta2[d_] = isF = (FLOAT)atof(pchar);

                        Error = (isF < (FLOAT)0.0) || (isF > (FLOAT)1.0) ; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "POISSON")) {
                        Theta_[c_].pdf[d_] = pfPoisson;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta1[d_] = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "DIRAC")) {
                        Theta_[c_].pdf[d_] = pfDirac;

                        pchar = strtok(NULL, "\t"); 

                        Theta_[c_].Theta1[d_] = isF = (FLOAT)atof(pchar);
                    }
                    else {
                        Error = 1; goto E0;
                    }

                    d_++;

                    pchar = strtok(NULL, "\t");
                }

                c_++;
            }
        } else
        if (!strcmp(ident, "SAVE")) {
            save_ = (char*)realloc(save_, (strlen(pchar) + 1) * sizeof(char));

            Error = NULL == save_; if (Error) goto E0;

            strcpy(save_, pchar);
        }
    }

E0: if (fp) fclose(fp);

    return (Error);
} // RunTemplateFile
