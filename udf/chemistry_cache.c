#include "udf.h"
#include "chemistry_cache.h"

void chemistry_cache_read_cell(cell_t c, Thread *tc, ChemState *s)
{
  int i;
  for (i=0;i<CHEM_N_SPECIES;i++) {
    s->theta[i] = C_UDMI(c, tc, UDM_THETA_CUR_BASE + i);
    s->theta_prev[i] = C_UDMI(c, tc, UDM_THETA_PREV_BASE + i);
  }
  for (i=0;i<CHEM_N_REACTIONS;i++) s->rates_export[i] = C_UDMI(c, tc, UDM_RATE_EXPORT_BASE + i);
  s->solver_status = (int)C_UDMI(c, tc, UDM_DIAG_STATUS);
  s->dt_last = C_UDMI(c, tc, UDM_DIAG_DT);
  s->err_last = C_UDMI(c, tc, UDM_DIAG_ERR);
}

void chemistry_cache_write_cell(cell_t c, Thread *tc, const ChemState *s)
{
  int i;
  for (i=0;i<CHEM_N_SPECIES;i++) {
    C_UDMI(c, tc, UDM_THETA_CUR_BASE + i) = s->theta[i];
    C_UDMI(c, tc, UDM_THETA_PREV_BASE + i) = s->theta_prev[i];
  }
  for (i=0;i<CHEM_N_REACTIONS;i++) C_UDMI(c, tc, UDM_RATE_EXPORT_BASE + i) = s->rates_export[i];
  C_UDMI(c, tc, UDM_DIAG_STATUS) = (real)s->solver_status;
  C_UDMI(c, tc, UDM_DIAG_DT) = s->dt_last;
  C_UDMI(c, tc, UDM_DIAG_ERR) = s->err_last;
}
