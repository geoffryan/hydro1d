#ifndef HYDRO1D_METRIC_SCHW_KS
#define HYDRO1D_METRIC_SCHW_KS

/*
 * Metric definitions for Schwarzschild metric in cylindrical Kerr-Schild 
 * coordinates.
 */

#define G400 (-1.0+2.0*M/r)
#define G40R (2.0*M/r)
#define G40P (0.0)
#define G4RR (1.0+2.0*M/r)
#define G4PP (r*r)
#define G4RP (0.0)
#define IG400 (-1.0-2.0*M/r)
#define IG40R (2.0*M/r)
#define IG40P (0.0)
#define IG4RR (1.0-2.0*M/r)
#define IG4PP (1.0/(r*r))
#define IG4RP (0.0)
#define DG400 (-2.0*M/(r*r))
#define DG40R (-2.0*M/(r*r))
#define DG40P (0.0)
#define DG4RR (-2.0*M/(r*r))
#define DG4PP (2.0*r)
#define DG4RP (0.0)
#define IG3RR (r/(2.0*M+r))
#define IG3PP (1.0/(r*r))
#define IG3RP (0.0)
#define J3 (sqrt(r*r+2.0*M*r))
#define J4 (r)
#define AL (1.0/sqrt(1.0+2.0*M/r))

#endif
