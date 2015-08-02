#ifndef HYDRO1D_METRIC_SCHW_SC
#define HYDRO1D_METRIC_SCHW_SC

/*
 * Metric definitions for Schwarzschild metric in Cylindrical Schwarzschild 
 * coordinates.
 */

#define G400 (-1.0+2.0*M/r)
#define G40R (0.0)
#define G40P (0.0)
#define G4RR (1.0/(1.0-2.0*M/r))
#define G4PP (r*r)
#define G4RP (0.0)
#define IG400 (1.0/(-1.0+2.0*M/r))
#define IG40R (0.0)
#define IG40P (0.0)
#define IG4RR (1.0-2.0*M/r)
#define IG4PP (1.0/(r*r))
#define IG4RP (0.0)
#define DG400 (-2.0*M/(r*r))
#define DG40R (0.0)
#define DG40P (0.0)
#define DG4RR (-2.0*M/(r*r) / ((1.0-2.0*M/r)*(1.0-2.0*M/r)))
#define DG4PP (2.0*r)
#define DG4RP (0.0)
#define IG3RR (1.0-2.0*M/r)
#define IG3PP (1.0/(r*r))
#define IG3RP (0.0)
#define J3 (r/sqrt(1-2*M/r))
#define J4 (r)
#define AL (sqrt(1.0-2.0*M/r))

#endif
