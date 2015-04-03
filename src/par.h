#ifndef ACCRETOR_PAR
#define ACCRETOR_PAR

enum{VAR_DBL, VAR_INT, VAR_LON};

struct parList
{
    int hydro;
    int eos;
    int cool;
    int step;
    int nx;
    int nghost;
    int nc;
    int nq;
    double xmin;
    double xmax;
};

const static struct parList PAR_DEFAULT = {
    .hydro = 0,
    .eos = 0,
    .cool = 0,
    .step = 0,
    .nx = 1,
    .nghost = 1,
    .nc = 3,
    .nq = 0,
    .xmin = 0.0,
    .xmax = 1.0
};

int readvar(char filename[], char key[], int vtype, void *ptr);
void read_pars(struct parList *theParList, char filename[]);

#endif
