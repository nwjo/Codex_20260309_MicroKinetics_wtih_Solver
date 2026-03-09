#ifndef CHEMISTRY_CACHE_H
#define CHEMISTRY_CACHE_H

#include "chem_state.h"

enum {
  UDM_BASE = 0,
  UDM_THETA_CUR_BASE = 0,
  UDM_THETA_PREV_BASE = UDM_THETA_CUR_BASE + CHEM_N_SPECIES,
  UDM_RATE_EXPORT_BASE = UDM_THETA_PREV_BASE + CHEM_N_SPECIES,
  UDM_DIAG_STATUS = UDM_RATE_EXPORT_BASE + CHEM_N_REACTIONS,
  UDM_DIAG_DT,
  UDM_DIAG_ERR,
  UDM_REQUIRED_COUNT
};

void chemistry_cache_read_cell(cell_t c, Thread *tc, ChemState *s);
void chemistry_cache_write_cell(cell_t c, Thread *tc, const ChemState *s);

#endif
