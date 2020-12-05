#ifndef RNGMIXF_H_INCLUDED
#define RNGMIXF_H_INCLUDED

#include "base.h"
#include "rebmixf.h"

class Rngmix : public Base {
    // Methods.
    #if (_MAINTAIN_SWITCH)
    int WriteDataFile();
    int WriteParameterFile();
    #endif
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
    FLOAT                **Y_;           // Dataset.
    int                  *N_;            // Numbers of observations.
    CompnentDistribution **MixTheta_;    // Mixture parameters.
    int                  *Z_;            // Component membership.
    // Constructor.
    Rngmix();
    // Destructor.
    virtual ~Rngmix();
    // Methods.
    virtual int InvComponentDist(CompnentDistribution *CmpDist, int j, FLOAT **Y);
    int RNGMIX();
    #if (_MAINTAIN_SWITCH)
    int RunTemplateFile(char *file);
    #endif
}; // Rngmix

#endif