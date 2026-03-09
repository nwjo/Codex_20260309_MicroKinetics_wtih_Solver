#include "udf.h"
#include <string.h>
#include <stdio.h>
#include "chem_constants.h"
#include "chem_state.h"
#include "chem_rates.h"
#include "chem_rhs.h"
#include "stiff_solver.h"
#include "chemistry_cache.h"
#include "legacy61_bridge.h"

static SolverGlobal g;


DEFINE_INIT(chemistry_init, domain)
{
  Thread *t;
  cell_t c;
  ChemState s;

  solver_global_init(&g);
  chem_state_set_default(&s);

  thread_loop_c(t, domain) {
    begin_c_loop(c, t) {
      chemistry_cache_write_cell(c, t, &s);
    } end_c_loop(c, t)
  }
}

DEFINE_ADJUST(chemistry_adjust, domain)
{
  Thread *tf;
  face_t f;
  const real dt = RP_Get_Real("physical-time-step");

  thread_loop_f(tf, domain) {
    if (!BOUNDARY_FACE_THREAD_P(tf)) continue;

    begin_f_loop(f, tf) {
      cell_t c0 = F_C0(f, tf);
      Thread *t0 = THREAD_T0(tf);
      ChemState s;
      ChemInputs in;
      int i;

      chemistry_cache_read_cell(c0, t0, &s);
      for (i=0;i<CHEM_N_SPECIES;i++) s.theta_prev[i] = s.theta[i];

      memset(&in, 0, sizeof(in));
      in.T_wall = F_T(f, tf);
      in.rho_cell = C_R(c0, t0);
      in.T_cell = C_T(c0, t0);
      in.yi[IDX_O2] = C_YI(c0, t0, IDX_O2);
      in.yi[IDX_C3H6] = C_YI(c0, t0, IDX_C3H6);
      in.yi[IDX_H2] = C_YI(c0, t0, IDX_H2);
      in.yi[IDX_H2O] = C_YI(c0, t0, IDX_H2O);
      in.yi[IDX_CO2] = C_YI(c0, t0, IDX_CO2);
      in.yi[IDX_CO] = C_YI(c0, t0, IDX_CO);
      in.yi[IDX_NO] = C_YI(c0, t0, IDX_NO);
      in.yi[IDX_NO2] = C_YI(c0, t0, IDX_NO2);
      in.yi[IDX_N2] = C_YI(c0, t0, IDX_N2);

      solver_prepare_problem(&in, &s);
      solver_integrate_surface_state(&in, dt, &s);

#if CHEM_USE_LEGACY_61_RATE_BRIDGE
      legacy61_eval_all_rates(f, tf, c0, t0, &in, &s);
#else
      chem_compute_intrinsic_rates(&in, &s);
      chem_apply_effectiveness_and_export(&s);
#endif

      chemistry_cache_write_cell(c0, t0, &s);
    } end_f_loop(f, tf)
  }
}

/* Compatibility wrapper: returns cached rates only. No reintegration. */
DEFINE_SR_RATE(chatterjee_cached_rates, f, fthread, r, mw, yi, rr)
{
  const cell_t c0 = F_C0(f, fthread);
  Thread *t0 = THREAD_T0(fthread);
  const int rid = chem_reaction_id_from_name(r->name);
  (void)mw; (void)yi;

  if (rid >= 1 && rid <= CHEM_N_REACTIONS) {
    *rr = C_UDMI(c0, t0, UDM_RATE_EXPORT_BASE + (rid-1));
  } else {
    *rr = 0.0;
  }
}
