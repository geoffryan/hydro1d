#ifndef HYDRO1D_BOUNDARY
#define HYDRO1D_BOUNDARY

void (*bc_inner)(struct grid *, struct parList *);
void (*bc_outer)(struct grid *, struct parList *);

int set_bc(struct parList *pars);
void bc_fixed_inner(struct grid *g, struct parList *pars);
void bc_fixed_outer(struct grid *g, struct parList *pars);
void bc_outflow_inner(struct grid *g, struct parList *pars);
void bc_outflow_outer(struct grid *g, struct parList *pars);
void bc_geodesic_inner(struct grid *g, struct parList *pars);
void bc_geodesic_outer(struct grid *g, struct parList *pars);

#endif
