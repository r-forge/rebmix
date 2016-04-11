#ifndef RNGMIXF_H_INCLUDED
#define RNGMIXF_H_INCLUDED

#include <stdlib.h>
#include <string.h>

#include "rebmixf.h"
#include "base.h"

class Rngmix : public Base {
    // Methods.
    int WriteDataFile();
    int WriteParameterFile();
public:
    // Members.
    char                 *curr_;         // Path to the currently open data file.
    int                  o_;             // Number of paths.
    char                 **open_;        // Paths to open data files.
    char                 *save_;         // Path to the save data file.
    int                  IDum_;          // Random seed.
    int                  c_;             // Number of components.
    CompnentDistribution *IniTheta_;     // Initial component parameters.
    int                  n_;             // Number of observations.
    char                 *Dataset_;      // Dataset name.
    FLOAT                **Y_;           // Dataset.
    int                  *N_;            // Numbers of observations.
    CompnentDistribution **MixTheta_;    // Mixture parameters.
    int                  *Z_;            // Component membership.
    // Constructor.
    Rngmix();
    // Destructor.
    virtual ~Rngmix();
    // Methods.
    virtual int InvComponentDist(CompnentDistribution *CmpDist, FLOAT *Y);
    int RNGMIX();
    int RunTemplateFile(char *file);
}; // Rngmix

#endif