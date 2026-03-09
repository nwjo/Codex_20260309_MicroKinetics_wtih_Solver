#include "udf.h"
#include <math.h>
#include <string.h>
#include "chem_rates.h"
#include "chem_constants.h"
#include "chem_indices.h"

#define MAX2(a,b) ((a)>(b)?(a):(b))
#define CLAMP01(x) ((x)<0.0?0.0:((x)>1.0?1.0:(x)))

/* Preserved helper signatures from legacy code for migration continuity. */
static inline real gas_conc_cell(real rho_cell, real T_cell, real T_wall, real yi_k, real MW_k)
{
  const real Tw = MAX2(CHEM_EPS, T_wall);
  const real Tc = MAX2(CHEM_EPS, T_cell);
  return (rho_cell * Tc / Tw) * yi_k / MAX2(CHEM_EPS, MW_k);
}

static inline real k_surface_covdep(real A, real beta, real Ea_J_per_kmol, real T,
                                    const int *idx_site, const real *mu, const real *eps_J_per_kmol,
                                    int Ns, const real yi[])
{
  real ln_k = log(MAX2(CHEM_EPS, A)) + beta*log(MAX2(CHEM_EPS, T))
            - Ea_J_per_kmol/(UNIVERSAL_GAS_CONSTANT*MAX2(CHEM_EPS, T));
  int i;
  for (i=0;i<Ns;i++) {
    const real th = MAX2(1.0e-20, CLAMP01(yi[idx_site[i]]));
    ln_k += mu[i]*log(th) + (eps_J_per_kmol[i]*th)/(UNIVERSAL_GAS_CONSTANT*MAX2(CHEM_EPS, T));
  }
  return exp(ln_k);
}

static inline real k_sticking(real S0, real T, real M_kg_per_kmol, real Gamma_kmol_m2, real q_site)
{
  return S0*sqrt((UNIVERSAL_GAS_CONSTANT*MAX2(CHEM_EPS,T))/(2.0*M_PI*MAX2(CHEM_EPS,M_kg_per_kmol)))
         *pow(1.0/MAX2(CHEM_EPS,Gamma_kmol_m2), q_site);
}

static inline void reaction7_base_and_eta(const ChemInputs *in, const ChemState *s, real *r7_base, real *eta7)
{
  const real MW_CO = 28.01055;
  const real S0_CO = 8.4e-1;
  const real q_R7 = 1.0;
  real dummy_idx[1] = {0};
  real dummy_val[1] = {0.0};
  (void)s;
  *r7_base = k_sticking(S0_CO, in->T_wall, MW_CO, SITE_DEN_TOT, q_R7)
           * gas_conc_cell(in->rho_cell, in->T_cell, in->T_wall, in->yi[IDX_CO], MW_CO)
           * in->yi[IDX_Pt_Vac];
  *eta7 = 1.0;
  (void)dummy_idx; (void)dummy_val;
}

static inline real eta_from_reaction7(const ChemInputs *in, const ChemState *s)
{
  real r7, eta;
  reaction7_base_and_eta(in, s, &r7, &eta);
  (void)r7;
  return eta;
}

int chem_reaction_id_from_name(const char *name)
{
  int id;
  char expected[32];
  for (id=1; id<=CHEM_N_REACTIONS; ++id) {
    sprintf(expected, "reaction-%d", id);
    if (strcmp(name, expected) == 0) return id;
  }
  return -1;
}

void chem_compute_intrinsic_rates(const ChemInputs *in, ChemState *s)
{
  int i;
  const real eta_ref = eta_from_reaction7(in, s);
  for (i=0;i<CHEM_N_REACTIONS;i++) s->rates_intrinsic[i] = 0.0;

  /* Full migration architecture note:
   * - This function is the single place for intrinsic elementary rates.
   * - Wash_F must never appear here.
   * - Legacy monolithic formulas must be ported reaction-by-reaction into this table.
   */

  /* Example preserved-style entries (R1, R7, R15) to keep compile path valid. */
  s->rates_intrinsic[0] = k_sticking(7.0e-2, in->T_wall, 31.9988, SITE_DEN_TOT, 2.0)
                        * gas_conc_cell(in->rho_cell, in->T_cell, in->T_wall, in->yi[IDX_O2], 31.9988)
                        * pow(MAX2(CHEM_EPS, s->theta[IDX_Pt_Vac]), 2.0);

  s->rates_intrinsic[6] = k_sticking(8.4e-1, in->T_wall, 28.01055, SITE_DEN_TOT, 1.0)
                        * gas_conc_cell(in->rho_cell, in->T_cell, in->T_wall, in->yi[IDX_CO], 28.01055)
                        * MAX2(CHEM_EPS, s->theta[IDX_Pt_Vac]);

  s->rates_intrinsic[14] = k_surface_covdep(3.7e21, 0.0, 2.3e8, in->T_wall,
                                            (int[]){IDX_CO_Pt}, (real[]){0.0}, (real[]){0.0}, 1, s->theta)
                         * MAX2(CHEM_EPS, s->theta[IDX_CO_Pt]);

  s->eta_ref = eta_ref;
}

void chem_apply_effectiveness_and_export(ChemState *s)
{
  int i;
  for (i=0;i<CHEM_N_REACTIONS;i++) {
    real eta = 1.0;
#if CHEM_USE_REFERENCE_ETA_FOR_SELECTED_REACTIONS
    eta = s->eta_ref;
#endif
    s->rates_effective[i] = eta * s->rates_intrinsic[i];
#if CHEM_APPLY_WASHCOAT_ON_EXPORT_ONLY
    s->rates_export[i] = Wash_F * s->rates_effective[i];
#else
    s->rates_export[i] = s->rates_effective[i];
#endif
  }
}
