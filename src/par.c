#include <stdio.h>
#include <string.h>
#include "par.h"

int readvar(char filename[], char key[], int vtype, void *ptr)
{
    FILE *f = fopen(filename, "r");

    char line[256];
    char word[256];
    int found = 0;

    while(fgets(line,256,f) != NULL)
    {
        sscanf(line, "%s ", word);
        if(strcmp(word,key) == 0)
        {
            found = 1;
            break;
        }
    }
    fclose(f);
    if(!found)
        return 0;

    char *sval = line + strlen(key) + strspn(line+strlen(key)," \t:=");

    if(vtype == VAR_DBL)
    {
        double val;
        sscanf(sval, "%lf", &val);
        *((double *)ptr) = val;
    }
    else if(vtype == VAR_INT)
    {
        int val;
        sscanf(sval, "%d", &val);
        *((int *)ptr) = val;
    }
    else if(vtype == VAR_LON)
    {
        long val;
        sscanf(sval, "%ld", &val);
        *((long *)ptr) = val;
    }
    else
    {
        strcpy((char *) ptr, sval);
    }

    return 0;
}

void read_pars(struct parList *theParList, char filename[])
{
    readvar(filename, "Hydro", VAR_INT, &(theParList->hydro));
    readvar(filename, "EOS",   VAR_INT, &(theParList->eos));
    readvar(filename, "Cool",  VAR_INT, &(theParList->cool));
    readvar(filename, "Timestep",  VAR_INT, &(theParList->step));
    readvar(filename, "Nx",     VAR_INT, &(theParList->nx));
    readvar(filename, "Nghost",     VAR_INT, &(theParList->nghost));
    readvar(filename, "Ncons",     VAR_INT, &(theParList->nc));
    readvar(filename, "Ntot",     VAR_INT, &(theParList->nq));
    readvar(filename, "Xmin",     VAR_DBL, &(theParList->xmin));
    readvar(filename, "Xmax",     VAR_DBL, &(theParList->xmax));
}

void print_pars(struct parList *theParList)
{
    printf("===Input Parameters===\n");
    printf("Hydro: %d\n", theParList->hydro);
    printf("EOS: %d\n", theParList->eos);
    printf("Cool: %d\n", theParList->cool);
    printf("Timestep: %d\n", theParList->step);
    printf("Nx: %d\n", theParList->nx);
    printf("Nghost: %d\n", theParList->nghost);
    printf("Ncons: %d\n", theParList->nc);
    printf("Ntot: %d\n", theParList->nq);
    printf("Xmin: %g\n", theParList->xmin);
    printf("Xmax: %g\n", theParList->xmax);
    printf("\n");
}
