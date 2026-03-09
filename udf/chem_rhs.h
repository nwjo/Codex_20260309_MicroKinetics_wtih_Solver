#ifndef CHEM_RHS_H
#define CHEM_RHS_H

#include "chem_rates.h"

void chem_build_rhs(const ChemInputs *in, const ChemState *state, real rhs[CHEM_N_SPECIES]);

#endif
