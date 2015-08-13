#ifndef RNGMIXF_H_INCLUDED
#define RNGMIXF_H_INCLUDED

#include <stdlib.h>
#include <string.h>

#include "base.h"

class Rngmix {
    // Methods.
    int WriteDataFile();
    int WriteParameterFile();
public:
    // Members.
    char                      *curr_;         // Path to the currently open data file.
    int                       o_;             // Number of paths.
    char                      **open_;        // Paths to open data files.
    char                      *save_;         // Path to the save data file.
    int                       IDum_;          // Random seed.
    int                       d_;             // Number of independent random variables.
    int                       c_;             // Number of components.
    int                       n_;             // Number of observations.
    char                      *Dataset_;      // Dataset name.
    FLOAT                     **Y_;           // Dataset.
    int                       *N_;            // Numbers of observations.
    ComponentDistributionType *Theta_;        // Component parameters.
    // Constructor.
    Rngmix();
    // Destructor.
    ~Rngmix();
    // Methods.
    virtual int InvComponentDist(ComponentDistributionType *CmpTheta, FLOAT *Y);
    int RNGMIX();
    int RunTemplateFile(char *file);
}; // Rngmix

#endif