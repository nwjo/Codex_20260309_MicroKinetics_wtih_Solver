#ifndef CHEM_RATES_H
#define CHEM_RATES_H

#include "chem_state.h"

typedef struct {
  real T_wall;
  real rho_cell;
  real T_cell;
  real yi[CHEM_N_SPECIES];
} ChemInputs;

void chem_compute_intrinsic_rates(const ChemInputs *in, ChemState *s);
void chem_apply_effectiveness_and_export(ChemState *s);
int chem_reaction_id_from_name(const char *name);

#endif
