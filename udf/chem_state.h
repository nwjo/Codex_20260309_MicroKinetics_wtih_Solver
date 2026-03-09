#ifndef CHEM_STATE_H
#define CHEM_STATE_H

#include "chem_indices.h"

typedef struct {
  real theta[CHEM_N_SPECIES];
  real theta_prev[CHEM_N_SPECIES];
  real rates_intrinsic[CHEM_N_REACTIONS];
  real rates_effective[CHEM_N_REACTIONS];
  real rates_export[CHEM_N_REACTIONS];
  real eta_ref;
  real dt_last;
  real err_last;
  int solver_status;
} ChemState;

void chem_state_set_default(ChemState *s);
void chem_state_project_constraints(ChemState *s);

#endif
