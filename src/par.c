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
    readvar(filename, "Recon",  VAR_INT, &(theParList->recon));
    readvar(filename, "Riemann",  VAR_INT, &(theParList->riemann));
    readvar(filename, "Timestep",  VAR_INT, &(theParList->step));
    readvar(filename, "Movement",  VAR_INT, &(theParList->movement));
    readvar(filename, "Frame",  VAR_INT, &(theParList->frame));
    readvar(filename, "BCInner",  VAR_INT, &(theParList->bcInner));
    readvar(filename, "BCOuter",  VAR_INT, &(theParList->bcOuter));
    readvar(filename, "Nx",     VAR_INT, &(theParList->nx));
    readvar(filename, "Nghost",     VAR_INT, &(theParList->nghost));
    readvar(filename, "Ncons",     VAR_INT, &(theParList->nc));
    readvar(filename, "Npass",     VAR_INT, &(theParList->np));
    readvar(filename, "Tmin",     VAR_DBL, &(theParList->tmin));
    readvar(filename, "Tmax",     VAR_DBL, &(theParList->tmax));
    readvar(filename, "Xmin",     VAR_DBL, &(theParList->xmin));
    readvar(filename, "Xmax",     VAR_DBL, &(theParList->xmax));
    readvar(filename, "PLM",     VAR_DBL, &(theParList->plm));
    readvar(filename, "CFL",     VAR_DBL, &(theParList->cfl));
    readvar(filename, "GammaLaw",     VAR_DBL, &(theParList->gammalaw));
    readvar(filename, "M",     VAR_DBL, &(theParList->M));
    readvar(filename, "NumCheckpoints",     VAR_INT, &(theParList->nChkpt));
    readvar(filename, "Init",         VAR_INT, &(theParList->init));
    readvar(filename, "InitPar1",     VAR_DBL, &(theParList->initPar1));
    readvar(filename, "InitPar2",     VAR_DBL, &(theParList->initPar2));
    readvar(filename, "InitPar3",     VAR_DBL, &(theParList->initPar3));
    readvar(filename, "InitPar4",     VAR_DBL, &(theParList->initPar4));
    readvar(filename, "InitPar5",     VAR_DBL, &(theParList->initPar5));
    readvar(filename, "InitPar6",     VAR_DBL, &(theParList->initPar6));
    readvar(filename, "InitPar7",     VAR_DBL, &(theParList->initPar7));
    readvar(filename, "InitPar8",     VAR_DBL, &(theParList->initPar8));
}

void print_pars(struct parList *theParList)
{
    printf("===Input Parameters===\n");
    printf("Hydro: %d\n", theParList->hydro);
    printf("EOS: %d\n", theParList->eos);
    printf("Cool: %d\n", theParList->cool);
    printf("Recon: %d\n", theParList->recon);
    printf("Riemann: %d\n", theParList->riemann);
    printf("Timestep: %d\n", theParList->step);
    printf("Movement: %d\n", theParList->movement);
    printf("Frame: %d\n", theParList->frame);
    printf("BCInner: %d\n", theParList->bcInner);
    printf("BCOuter: %d\n", theParList->bcOuter);
    printf("Nx: %d\n", theParList->nx);
    printf("Nghost: %d\n", theParList->nghost);
    printf("Ncons: %d\n", theParList->nc);
    printf("Npass: %d\n", theParList->np);
    printf("Tmin: %g\n", theParList->tmin);
    printf("Tmax: %g\n", theParList->tmax);
    printf("Xmin: %g\n", theParList->xmin);
    printf("Xmax: %g\n", theParList->xmax);
    printf("PLM: %g\n", theParList->plm);
    printf("CFL: %g\n", theParList->cfl);
    printf("GammaLaw: %g\n", theParList->gammalaw);
    printf("M: %g\n", theParList->M);
    printf("NumCheckpoints: %d\n", theParList->nChkpt);
    printf("Init: %d\n", theParList->init);
    printf("InitPar1: %g\n", theParList->initPar1);
    printf("InitPar2: %g\n", theParList->initPar2);
    printf("InitPar3: %g\n", theParList->initPar3);
    printf("InitPar4: %g\n", theParList->initPar4);
    printf("InitPar5: %g\n", theParList->initPar5);
    printf("InitPar6: %g\n", theParList->initPar6);
    printf("InitPar7: %g\n", theParList->initPar7);
    printf("InitPar8: %g\n", theParList->initPar8);
    printf("\n");
}
