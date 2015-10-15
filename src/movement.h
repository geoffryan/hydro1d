#ifndef HYDRO1D_MOVEMENT
#define HYDRO1D_MOVEMENT

void (*calc_grid_movement)(struct grid *, struct parList *);

int set_movement(struct parList *pars);

void move_none(struct grid *g, struct parList *pars);
void move_local(struct grid *g, struct parList *pars);

#endif
