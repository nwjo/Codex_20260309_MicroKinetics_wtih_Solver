#include "udf.h"
#include <math.h>
#include <string.h>

/*
 * Homogeneous gas-phase net reaction-rate UDF
 * Macro prototype required by Fluent manual:
 * void name(cell_t c, Thread *t, Particle *particle,
 *           double *pressure, double *temp, double *yi,
 *           double *rr, double *jac)
 *
 * [REACTION_MODEL_DESCRIPTION]
 * - Placeholder model implemented here: single irreversible global reaction
 *   O2 + CO -> CO2 with Arrhenius kinetics.
 * - Replace with project mechanism and stoichiometry as needed.
 *
 * [SPECIES_INDEX_MAP]
 * IMPORTANT: These indices must match Fluent species order in the material panel.
 */
#define NSPEC 9         /* [NSPEC] Replace with actual Fluent species count */

enum {
  SPEC_O2  = 0,
  SPEC_CO  = 1,
  SPEC_CO2 = 2
  /* Add remaining species indices in Fluent order */
};

/* [KINETIC_PARAMETERS] */
#define R_UNIVERSAL 8.31446261815324        /* J/(mol*K) */
#define A_PREEXP    2.50e7                  /* placeholder units for concentration-based rate */
#define E_ACT       8.00e4                  /* J/mol */
#define MIN_POS     1.0e-30

/* Placeholder molecular weights [kg/mol], align with Fluent species order */
static void get_mw(double mw[NSPEC])
{
  int k;
  for (k = 0; k < NSPEC; ++k) mw[k] = 0.028; /* placeholder default */
  mw[SPEC_O2]  = 0.0319988;
  mw[SPEC_CO]  = 0.02801055;
  mw[SPEC_CO2] = 0.04400995;
}

/* Convert Yi -> Xi, Ctot and Ck using pressure[0], temp[0], yi[] from Fluent */
static void state_from_fluent_inputs(const double p,
                                     const double T,
                                     const double yi[NSPEC],
                                     const double mw[NSPEC],
                                     double xi[NSPEC],
                                     double ck[NSPEC],
                                     double *mw_mix,
                                     double *ctot)
{
  int k;
  double sum_y_over_mw = 0.0;

  for (k = 0; k < NSPEC; ++k)
    sum_y_over_mw += yi[k] / MAX(mw[k], MIN_POS);

  *mw_mix = 1.0 / MAX(sum_y_over_mw, MIN_POS);
  *ctot = p / MAX(R_UNIVERSAL * T, MIN_POS); /* mol/m^3 */

  for (k = 0; k < NSPEC; ++k) {
    xi[k] = yi[k] * (*mw_mix) / MAX(mw[k], MIN_POS);
    ck[k] = xi[k] * (*ctot);
  }
}

/*
 * [JACOBIAN_LAYOUT_IF_KNOWN]
 * This implementation assumes dense row-major storage for an NSPEC x NSPEC Jacobian:
 *   jac[i*NSPEC + j] = d(rr_i)/d(yi_j)
 * If your Fluent setup expects a different layout/variable basis, replace accordingly.
 */
static void clear_outputs(double rr[NSPEC], double jac[NSPEC*NSPEC])
{
  int i;
  for (i = 0; i < NSPEC; ++i) rr[i] = 0.0;
  for (i = 0; i < NSPEC*NSPEC; ++i) jac[i] = 0.0;
}

DEFINE_NET_REACTION_RATE(homogeneous_net_rates, c, t, particle, pressure, temp, yi, rr, jac)
{
  int i, j;
  double p = pressure[0];
  double T = temp[0];
  double mw[NSPEC];
  double xi[NSPEC], ck[NSPEC], mw_mix, ctot;
  double kf, rop;  /* forward rate constant, rate-of-progress */
  double dck_dyj;

  (void)c;
  (void)t;
  (void)particle;

  /* Input argument meaning:
   * c,t      : cell and thread where source is evaluated
   * particle : DPM particle pointer (unused here)
   * pressure : local absolute pressure array from Fluent (pressure[0])
   * temp     : local temperature array from Fluent (temp[0])
   * yi       : species mass-fraction array from Fluent (must match species order)
   * rr       : output net molar rate array for each species [mol/m^3/s]
   * jac      : output Jacobian terms used by Fluent linearization
   */

  clear_outputs(rr, jac);
  get_mw(mw);
  state_from_fluent_inputs(p, T, yi, mw, xi, ck, &mw_mix, &ctot);

  /* Conservative Arrhenius model with concentration basis.
   * rop = kf * C_O2 * C_CO
   */
  kf = A_PREEXP * exp(-E_ACT / MAX(R_UNIVERSAL * T, MIN_POS));
  rop = kf * ck[SPEC_O2] * ck[SPEC_CO];

  /* rr assembly:
   * rr[k] is net molar production/destruction rate of species k.
   * For O2 + CO -> CO2:
   *   nu_O2  = -1
   *   nu_CO  = -1
   *   nu_CO2 = +1
   */
  rr[SPEC_O2]  -= rop;
  rr[SPEC_CO]  -= rop;
  rr[SPEC_CO2] += rop;

  /* jac assembly (d(rr_i)/d(yi_j)):
   * ck[s] = yi[s]*mw_mix/mw[s]*ctot
   * neglecting d(mw_mix)/d(yi_j) for conservative robustness placeholder.
   * Replace with full derivatives for production use.
   */
  for (j = 0; j < NSPEC; ++j) {
    double drop_dyj = 0.0;

    dck_dyj = ((j == SPEC_O2) ? 1.0 : 0.0) * (mw_mix / MAX(mw[SPEC_O2], MIN_POS)) * ctot;
    drop_dyj += kf * dck_dyj * ck[SPEC_CO];

    dck_dyj = ((j == SPEC_CO) ? 1.0 : 0.0) * (mw_mix / MAX(mw[SPEC_CO], MIN_POS)) * ctot;
    drop_dyj += kf * ck[SPEC_O2] * dck_dyj;

    jac[SPEC_O2*NSPEC + j]  -= drop_dyj;
    jac[SPEC_CO*NSPEC + j]  -= drop_dyj;
    jac[SPEC_CO2*NSPEC + j] += drop_dyj;
  }

  /* Keep unused-loop variable warning-free across compilers configured by Fluent */
  i = 0;
  (void)i;
}
