#include <stdio.h>
#include "par.h"

int main(int argc, char *argv[])
{
    struct parList pars = PAR_DEFAULT;

    if(argc != 2)
    {
        printf("usage: hydro1c <parameter file>\n");
        return 0;
    }

    printf("\nWelcome to Hydro-1D!\n\n");

    read_pars(&pars, argv[1]);
    print_pars(&pars);

    return 0;
}
