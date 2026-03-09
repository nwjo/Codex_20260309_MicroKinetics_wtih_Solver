#ifndef STIFF_SOLVER_H
#define STIFF_SOLVER_H

#include "chem_jacobian.h"

typedef struct {
  int initialized;
} SolverGlobal;

void solver_global_init(SolverGlobal *g);
void solver_prepare_problem(const ChemInputs *in, const ChemState *s);
int solver_integrate_surface_state(const ChemInputs *in, real dt, ChemState *s);
void solver_free(SolverGlobal *g);

#endif
