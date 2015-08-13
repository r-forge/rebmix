#include <stdio.h>

#include "base.h"
#include "rngmixf.h"
#include "rebmixf.h"

#if (_MEMORY_LEAK_SWITCH)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

int main(int argc, char* argv[])
{
    #if (_MEMORY_LEAK_SWITCH)
    _CrtMemState s1, s2, s3;

    _CrtMemCheckpoint(&s1);
    #endif

    Rngmix *rngmix = NULL;
    Rebmix *rebmix = NULL;
    int    Error = 0;

    if (argc != 3) goto E0;

    if (!strcmp(argv[2], "RNGMIX")) {
        rngmix = new Rngmix;

        Error = rngmix->RunTemplateFile(argv[1]);

        if (Error) goto E0;
    }
    else
    if (!strcmp(argv[2], "REBMIX")) {
        rebmix = new Rebmix;
 
        Error = rebmix->RunTemplateFile(argv[1]);

        if (Error) goto E0;
    }

E0: if (rngmix) delete(rngmix);
    if (rebmix) delete(rebmix);

    #if (_MEMORY_LEAK_SWITCH)
    _CrtMemCheckpoint(&s2);

    if (_CrtMemDifference(&s3, &s1, &s2)) _CrtMemDumpStatistics(&s3);
    #endif
    
    #if (_REBMIXEXE)
    printf("\n%s%d\n", "Error: ", Error);
    #endif

    return (Error);
} // main
