#include <stdio.h>

#include "rngmixf.h"
#include "rebmixf.h"

int main(int argc, char* argv[])
{
    int Error = 0;

    if (argc != 3) goto E0;

    if (!strcmp(argv[2], "RNGMIX")) {
        Error = RunRNGMIXTemplateFile(argv[1]);

        if (Error) goto E0;
    }
    else
    if (!strcmp(argv[2], "REBMIX")) {
        Error = RunREBMIXTemplateFile(argv[1]);

        if (Error) goto E0;
    }

E0:
    #if (_REBMIXEXE)
    printf("\n%s%d\n", "Error: ", Error);
    #endif
    
    return (Error);
} /* main */
