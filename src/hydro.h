#ifndef HYDRO1D_HYDRO
#define HYDRO1D_HYDRO

enum{RHO, PPP, VX1, VX2};
enum{DDD, TAU, SX1, SX2};

void (*prim2cons)(double *prim, double *cons, double x, double dV,
                    struct parList *pars);
void (*cons2prim)(double *cons, double *prim, double x, double dV,
                    struct parList *pars);
void (*flux)(double *prim, double *F, double x, struct parList *pars);
void (*add_source)(double *prim, double *cons, double x, double dVdt, 
                    struct parList *pars);
void (*wave_speeds)(double *prim1, double *prim2, double *sL, double *sR,
                    double *sC, double x, struct parList *pars);
double (*mindt)(double *prim, double x, double dx, struct parList *pars);

int set_hydro(struct parList *pars);

void prim2cons_newt_cart(double *prim, double *cons, double x, double dV,
                            struct parList *pars);
void cons2prim_newt_cart(double *cons, double *prim, double x, double dV,
                            struct parList *pars);
void flux_newt_cart(double *prim, double *F, double x, struct parList *pars);
void add_source_newt_cart(double *prim, double *cons, double x, double dVdt, 
                            struct parList *pars);
void wave_speeds_newt_cart(double *prim1, double *prim2, double *sL, 
                            double *sR, double *sC, double x, 
                            struct parList *pars);
double mindt_newt_cart(double *prim, double x, double dx, 
                        struct parList *pars);

void prim2cons_newt_cyl(double *prim, double *cons, double x, double dV,
                            struct parList *pars);
void cons2prim_newt_cyl(double *cons, double *prim, double x, double dV,
                            struct parList *pars);
void flux_newt_cyl(double *prim, double *F, double x, struct parList *pars);
void add_source_newt_cyl(double *prim, double *cons, double x, double dVdt, 
                            struct parList *pars);
void wave_speeds_newt_cyl(double *prim1, double *prim2, double *sL, 
                            double *sR, double *sC, double x, 
                            struct parList *pars);
double mindt_newt_cyl(double *prim, double x, double dx, 
                        struct parList *pars);

void prim2cons_newt_sph(double *prim, double *cons, double x, double dV,
                            struct parList *pars);
void cons2prim_newt_sph(double *cons, double *prim, double x, double dV,
                            struct parList *pars);
void flux_newt_sph(double *prim, double *F, double x, struct parList *pars);
void add_source_newt_sph(double *prim, double *cons, double x, double dVdt, 
                            struct parList *pars);
void wave_speeds_newt_sph(double *prim1, double *prim2, double *sL, 
                            double *sR, double *sC, double x, 
                            struct parList *pars);
double mindt_newt_sph(double *prim, double x, double dx, 
                        struct parList *pars);

void prim2cons_rel_cart(double *prim, double *cons, double x, double dV,
                            struct parList *pars);
void cons2prim_rel_cart(double *cons, double *prim, double x, double dV,
                            struct parList *pars);
void flux_rel_cart(double *prim, double *F, double x, struct parList *pars);
void add_source_rel_cart(double *prim, double *cons, double x, double dVdt, 
                            struct parList *pars);
void wave_speeds_rel_cart(double *prim1, double *prim2, double *sL, 
                            double *sR, double *sC, double x, 
                            struct parList *pars);
double mindt_rel_cart(double *prim, double x, double dx, 
                        struct parList *pars);

void prim2cons_rel_metric(double *prim, double *cons, double x, double dV,
                            struct parList *pars);
void cons2prim_rel_metric(double *cons, double *prim, double x, double dV,
                            struct parList *pars);
void flux_rel_metric(double *prim, double *F, double x, struct parList *pars);
void add_source_rel_metric(double *prim, double *cons, double x, double dVdt, 
                            struct parList *pars);
void wave_speeds_rel_metric(double *prim1, double *prim2, double *sL, 
                            double *sR, double *sC, double x, 
                            struct parList *pars);
double mindt_rel_metric(double *prim, double x, double dx, 
                        struct parList *pars);

void (*frame_U)(double, double, double *);
void (*frame_dU)(double, double, double *);
void frame_normal_U(double r, double M, double U[]);
void frame_normal_dU(double r, double M, double dU[]);
void frame_geo_schw_U(double r, double M, double U[]);
void frame_geo_schw_dU(double r, double M, double dU[]);

#endif
