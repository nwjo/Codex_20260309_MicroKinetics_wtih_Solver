#ifndef LEGACY61_BRIDGE_H
#define LEGACY61_BRIDGE_H

#include "udf.h"
#include "chem_state.h"
#include "chem_rates.h"

/*
 * Evaluate all reaction-1..reaction-61 using legacy mechanism formulas
 * internalized in udf/chem_internalized_61rxn.c
 * and store migrated state arrays.
 */
void legacy61_eval_all_rates(face_t f, Thread *tf, cell_t c0, Thread *t0,
                             const ChemInputs *in, ChemState *s);

#endif
