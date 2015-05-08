#include <stdio.h>
#include "grid.h"
#include "par.h"

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        printf("usage: hydro1c <parameter file>\n");
        return 0;
    }

    printf("\nWelcome to Hydro-1D!\n\n");

    struct parList pars = PAR_DEFAULT;
    struct grid grid = GRID_DEFAULT;

    read_pars(&pars, argv[1]);
    print_pars(&pars);

    make_grid(&grid, &pars);
    print_grid(&grid, "grid.txt");
    free_grid(&grid);

    return 0;
}
