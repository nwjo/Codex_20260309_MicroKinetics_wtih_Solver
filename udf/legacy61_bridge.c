#include "udf.h"
#include <stdio.h>
#include <string.h>
#include "legacy61_bridge.h"
#include "chem_constants.h"
#include "chem_indices.h"

/* Internalized DEFINE_SR_RATE logic exported by udf/chem_internalized_61rxn.c */
extern void chatterjee_pt_ads_des_flat_internalized(face_t f, Thread *fthread, Reaction *r,
                                                     real *mw, real *yi, real *rr);

void legacy61_eval_all_rates(face_t f, Thread *tf, cell_t c0, Thread *t0,
                             const ChemInputs *in, ChemState *s)
{
  int rid;
  char name_buf[32];
  Reaction rlocal;
  real mw_dummy[CHEM_N_SPECIES];
  real yi_mix[CHEM_N_SPECIES];

  memset(&rlocal, 0, sizeof(rlocal));
  memset(mw_dummy, 0, sizeof(mw_dummy));
  memset(yi_mix, 0, sizeof(yi_mix));

  /* Gas composition from Fluent cell state (adjacent to catalytic wall). */
  yi_mix[IDX_O2] = in->yi[IDX_O2];
  yi_mix[IDX_C3H6] = in->yi[IDX_C3H6];
  yi_mix[IDX_H2] = in->yi[IDX_H2];
  yi_mix[IDX_H2O] = in->yi[IDX_H2O];
  yi_mix[IDX_CO2] = in->yi[IDX_CO2];
  yi_mix[IDX_CO] = in->yi[IDX_CO];
  yi_mix[IDX_NO] = in->yi[IDX_NO];
  yi_mix[IDX_NO2] = in->yi[IDX_NO2];
  yi_mix[IDX_N2] = in->yi[IDX_N2];

  /* Surface state from migrated cache (site species indices only). */
  for (rid = IDX_Pt_Vac; rid < CHEM_N_SPECIES; ++rid) {
    yi_mix[rid] = s->theta[rid];
  }

  /*
   * Evaluate all legacy branches reaction-1..reaction-61.
   * This guarantees all reaction equations/coefficients in internalized chemistry file
   * are reflected in migrated cached results.
   */
  for (rid = 1; rid <= CHEM_N_REACTIONS; ++rid) {
    real rr_export = 0.0;
    snprintf(name_buf, sizeof(name_buf), "reaction-%d", rid);
    rlocal.name = name_buf;
    chatterjee_pt_ads_des_flat_internalized(f, tf, &rlocal, mw_dummy, yi_mix, &rr_export);

    s->rates_export[rid - 1] = rr_export;

#if CHEM_APPLY_WASHCOAT_ON_EXPORT_ONLY
    s->rates_effective[rid - 1] = rr_export / Wash_F;
#else
    s->rates_effective[rid - 1] = rr_export;
#endif

    /*
     * Legacy dispatcher already embeds eta treatment in selected branches.
     * During bridge mode, intrinsic decomposition is not fully recoverable per-branch,
     * so keep intrinsic alias equal to effective as migration-compatible fallback.
     */
    s->rates_intrinsic[rid - 1] = s->rates_effective[rid - 1];
  }

  (void)c0;
  (void)t0;
}
