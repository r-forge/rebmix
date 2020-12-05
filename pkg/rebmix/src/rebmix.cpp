#include "rngmvnormf.h"
#include "rebmvnormf.h"

#if (_MEMORY_LEAK_SWITCH)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

#if (_MAINTAIN_SWITCH)
#include <stdio.h>

int main(int argc, char* argv[])
{
    #if (_MEMORY_LEAK_SWITCH)
    _CrtMemState s1, s2, s3;

    _CrtMemCheckpoint(&s1);
    #endif

    Rngmix    *rngmix = NULL;
    Rebmix    *rebmix = NULL;
    Rngmvnorm *rngmvnorm = NULL;
    Rebmvnorm *rebmvnorm = NULL;
    int       Error = 0;

    if (argc != 3) goto E0;

    if (!strcmp(argv[2], "RNGMIX")) {
        rngmix = new Rngmix;

        Error = NULL == rngmix; if (Error) goto E0;

        Error = rngmix->RunTemplateFile(argv[1]);

        if (Error) goto E0;
    }
    else
    if (!strcmp(argv[2], "REBMIX")) {
        rebmix = new Rebmix;

        Error = NULL == rebmix; if (Error) goto E0;

        Error = rebmix->RunTemplateFile(argv[1]);

        if (Error) goto E0;
    }
    else
    if (!strcmp(argv[2], "REBMIX")) {
        rebmix = new Rebmix;

        Error = NULL == rebmix; if (Error) goto E0;

        Error = rebmix->RunTemplateFile(argv[1]);

        if (Error) goto E0;
    }
    else
    if (!strcmp(argv[2], "REBMVNORM")) {
        rebmvnorm = new Rebmvnorm;

        Error = NULL == rebmvnorm; if (Error) goto E0;

        Error = rebmvnorm->RunTemplateFile(argv[1]);

        if (Error) goto E0;
    }

E0: if (rngmix) delete rngmix;
    if (rebmix) delete rebmix;

    if (rngmvnorm) delete rngmvnorm;
    if (rebmvnorm) delete rebmvnorm;

    #if (_MEMORY_LEAK_SWITCH)
    _CrtMemCheckpoint(&s2);

    if (_CrtMemDifference(&s3, &s1, &s2)) _CrtMemDumpStatistics(&s3);
    #endif

    printf("\n%s%d\n", "Error: ", Error);

    return Error;
} // main
#else
int main()
{
    return 0;
}
#endif
