#include "udf.h"
#include <string.h>
#include "stiff_solver.h"
#include "chem_constants.h"
#include "chem_state.h"

static SolverGlobal g_solver;

void solver_global_init(SolverGlobal *g) { g->initialized = 1; g_solver = *g; }
void solver_prepare_problem(const ChemInputs *in, const ChemState *s) { (void)in; (void)s; }

int solver_integrate_surface_state(const ChemInputs *in, real dt, ChemState *s)
{
  int i;
  real rhs[CHEM_N_SPECIES];
  ChemState old = *s;

  chem_build_rhs(in, &old, rhs);
  for (i=0;i<CHEM_N_SPECIES;i++) s->theta[i] = old.theta[i] + dt*rhs[i];
  chem_state_project_constraints(s);

  s->dt_last = dt;
  s->err_last = 0.0;
  s->solver_status = 0;
  return 0;
}

void solver_free(SolverGlobal *g) { g->initialized = 0; g_solver = *g; }
