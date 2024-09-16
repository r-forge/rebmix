#ifndef RNGMIXF_H_INCLUDED
#define RNGMIXF_H_INCLUDED

#include "base.h"
#include "rebmixf.h"

class Rngmix : public Base {
    // Methods.
    #if (_MAINTAIN_SWITCH)
    INT WriteDataFile();
    INT WriteParameterFile();
    #endif
public:
    // Members.
    char                 *curr_;         // Path to the currently open data file.
    INT                  o_;             // Number of paths.
    char                 **open_;        // Paths to open data files.
    char                 *save_;         // Path to the save data file.
    INT                  IDum_;          // Random seed.
    INT                  c_;             // Number of components.
    CompnentDistribution *IniTheta_;     // Initial component parameters.
    INT                  n_;             // Number of observations.
    FLOAT                **Y_;           // Dataset.
    INT                  *N_;            // Numbers of observations.
    CompnentDistribution **MixTheta_;    // Mixture parameters.
    INT                  *Z_;            // Component membership.
    // Constructor.
    Rngmix();
    // Destructor.
    virtual ~Rngmix();
    // Methods.
    virtual INT ComponentInv(CompnentDistribution *CmpCdf, INT j, FLOAT **Y);
    INT RNGMIX();
    #if (_MAINTAIN_SWITCH)
    INT RunTemplateFile(char *file);
    #endif
}; // Rngmix

#endif
