#include "udf.h"
#include "chem_state.h"
#include "chem_constants.h"
#include "chem_indices.h"

#define MAX2(a,b) ((a)>(b)?(a):(b))

void chem_state_set_default(ChemState *s)
{
  int i;
  for (i=0;i<CHEM_N_SPECIES;i++) {
    s->theta[i] = 0.0;
    s->theta_prev[i] = 0.0;
  }
  for (i=0;i<CHEM_N_REACTIONS;i++) {
    s->rates_intrinsic[i] = 0.0;
    s->rates_effective[i] = 0.0;
    s->rates_export[i] = 0.0;
  }
  s->theta[IDX_Pt_Vac] = 1.0;
  s->theta[IDX_Rh_Vac] = 1.0;
  s->eta_ref = 1.0;
  s->dt_last = 0.0;
  s->err_last = 0.0;
  s->solver_status = 0;
}

void chem_state_project_constraints(ChemState *s)
{
#if CHEM_ENFORCE_SITE_BALANCE_PROJECTION
  real sum_pt = 0.0, sum_rh = 0.0;
  int i;

  for (i = IDX_Pt_Vac; i <= IDX_CH3CO_Pt; ++i) {
    s->theta[i] = MAX2(0.0, s->theta[i]);
    sum_pt += s->theta[i];
  }
  for (i = IDX_Rh_Vac; i <= IDX_N_Rh; ++i) {
    s->theta[i] = MAX2(0.0, s->theta[i]);
    sum_rh += s->theta[i];
  }

  if (sum_pt > CHEM_EPS) {
    for (i = IDX_Pt_Vac; i <= IDX_CH3CO_Pt; ++i) s->theta[i] /= sum_pt;
  }
  if (sum_rh > CHEM_EPS) {
    for (i = IDX_Rh_Vac; i <= IDX_N_Rh; ++i) s->theta[i] /= sum_rh;
  }
#endif
}
