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
    INT       Error = EOK;

    if (argc != 3) goto EEXIT;

    if (!strcmp(argv[2], "RNGMIX")) {
        rngmix = new Rngmix;

        E_CHECK(NULL == rngmix, Emain);

        Error = rngmix->RunTemplateFile(argv[1]);

        E_CHECK(Error != EOK, Error);
    }
    else
    if (!strcmp(argv[2], "REBMIX")) {
        rebmix = new Rebmix;

        E_CHECK(NULL == rebmix, Emain);

        Error = rebmix->RunTemplateFile(argv[1]);

        E_CHECK(Error != EOK, Error);
    }
    else
    if (!strcmp(argv[2], "REBMVNORM")) {
        rebmvnorm = new Rebmvnorm;

        E_CHECK(NULL == rebmvnorm, Emain);

        Error = rebmvnorm->RunTemplateFile(argv[1]);

        E_CHECK(Error != EOK, Error);
    }

EEXIT: 
    
    if (rngmix) delete rngmix;
    if (rebmix) delete rebmix;

    if (rngmvnorm) delete rngmvnorm;
    if (rebmvnorm) delete rebmvnorm;

    #if (_MEMORY_LEAK_SWITCH)
    _CrtMemCheckpoint(&s2);

    if (_CrtMemDifference(&s3, &s1, &s2)) _CrtMemDumpStatistics(&s3);
    #endif

    printf("\n%s%d\n", "Error: ", Error);

    E_RETURN(Error);
} // main
#else
int main()
{
    return 0;
}
#endif
