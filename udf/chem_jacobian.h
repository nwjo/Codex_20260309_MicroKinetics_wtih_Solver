#ifndef CHEM_JACOBIAN_H
#define CHEM_JACOBIAN_H

#include "chem_rhs.h"

void chem_build_jacobian(const ChemInputs *in, const ChemState *state,
                         real jac[CHEM_N_SPECIES][CHEM_N_SPECIES]);

#endif
