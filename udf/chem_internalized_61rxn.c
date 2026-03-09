/*
 * Based on Chatterjee et al., 2001
 * Rewritten by NWJO, KATECH
 * Date: 2025-10-29 (Refactored Version)
 * coverage-dependent k^(s) applied
 * Update: Added full surface reaction set (R1-R61) as continuous else-if chain
 * Notes:
 *  - gas_conc_cell uses species molecular weight for ideal-gas concentration.
 *  - surface intermediates indices expanded for hydrocarbon oxidation network.
 *  - Thiele modulus and effectiveness factor calculation based on power-law diffusion.
 *  - Reaction-7 (CO adsorption) used as reference for Thiele modulus and eta.
 *  - UDM slots added for reaction rate diagnostics.
 *  - Adsorption, desorption, and surface reaction rate constants separated for clarity.
 *  - Pre-exponential factors calculated from sticking coefficients.
 *  - Added one-time logging for adsorption reaction rates.
 *  - Added guards for macros (MAX, M_PI, CLAMP01, STREQ).
 *  - Refactored coverage-dependent rate constant calculation into a helper function.
 *  - Reaction rate constants are based on Chatterjee et al., 2001 formulations.
 *    - Note that ANSYS Fluent using different formulation for reaction rate constants.
 *    - ANSYS Fluent use negative exponents for coverage dependent activation energy.
*/

#include "udf.h"
#include <math.h>
#include <string.h>

/* ---------- Print Gates ---------- */
/* Using an array for cleaner management */
#define NUM_REACTIONS 65
static int print_gate[NUM_REACTIONS];

DEFINE_ADJUST(r1_reset_print_gate, d)
{
    /* Reset all print gates */
    int i;
    for (i = 0; i < NUM_REACTIONS; i++) print_gate[i] = 0;
}

/* ---------- UDM slots ---------- */
enum {
  UDM_RCO_NET   = 0,        /* r_CO,net (without base, eta)         [kmol/m^2/s] */
  UDM_RCO_NETP  = 1,        /* r_CO,net with Cco+(1+eps)            [kmol/m^2/s] */
  UDM_DRDC      = 2,        /* dr/dCco (coverage fixed)             [m/s]        */
  UDM_KAPP      = 3,        /* k_app,CO = a_w * dr/dCco             [1/s]        */
  UDM_PHI       = 4,        /* Thiele modulus phi                   [-]         */
  UDM_ETA       = 5         /* eta = tanh(phi)/phi                  [-]         */
};

/* guards */
#ifndef MAX
#define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef CLAMP01
#define CLAMP01(x) ((x) < 0.0 ? 0.0 : ((x) > 1.0 ? 1.0 : (x)))
#endif

#ifndef STREQ
#define STREQ(a,b) (strcmp((a),(b))==0)
#endif

#define EPS 1.0E-30

/*** species index mapping (Combined) ***/
#define IDX_O2 0
#define IDX_C3H6 1
#define IDX_H2 2
#define IDX_H2O 3
#define IDX_CO2 4
#define IDX_CO 5
#define IDX_NO 6
#define IDX_NO2 7
#define IDX_N2 8
#define IDX_Pt_Vac 9
#define IDX_O_Pt 10
#define IDX_C3H6_Pt 11
#define IDX_H_Pt 12
#define IDX_H2O_Pt 13
#define IDX_CO2_Pt 14
#define IDX_CO_Pt 15
#define IDX_NO_Pt 16
#define IDX_N_Pt 17
#define IDX_C3H5_Pt 18
#define IDX_C2H3_Pt 19
#define IDX_CH2_Pt 20
#define IDX_CH3_Pt 21
#define IDX_OH_Pt 22
#define IDX_CH_Pt 23
#define IDX_C_Pt 24
#define IDX_CC2H5_Pt 25
#define IDX_CH3CO_Pt 26
#define IDX_Rh_Vac 27
#define IDX_O_Rh 28
#define IDX_CO_Rh 29
#define IDX_NO_Rh 30
#define IDX_N_Rh 31

/*** Constants ***/
#define SITE_DEN_TOT 2.72E-8                        /* kmol/m2 */
#define Pt_Frac 0.75
#define SITE_DEN_Pt (SITE_DEN_TOT*Pt_Frac)
#define SITE_DEN_Rh (SITE_DEN_TOT*(1.0-Pt_Frac))

// Molecular Weight [kg/kmol]
#define MW_O2   31.9988
#define MW_C3H6 42.081
#define MW_H2   2.01594
#define MW_H2O  18.01534
#define MW_CO2  44.00995
#define MW_CO   28.01055
#define MW_NO   30.0061

// Sticking Coefficients S0
#define S0_O2       7.0E-2      // R-1
#define S0_C3H6     9.8E-1      // R-2
#define S0_C3H6_O   5.00E-2     // R-3
#define S0_H2       4.60E-2     // R-4
#define S0_H2O      7.50E-1     // R-5
#define S0_CO2      5.00E-3     // R-6
#define S0_CO       8.40E-1     // R-7
#define S0_NO       8.50E-1     // R-48
#define S1_O2       1.00E-2     // R-53
#define S1_CO       5.00E-1     // R-54
#define S1_NO       5.00E-1     // R-55

// Vacant Site Requirement
#define q_R1    2.0
#define q_R2    2.0
#define q_R3    2.0
#define q_R4    2.0
#define q_R5    1.0
#define q_R6    1.0
#define q_R7    1.0
#define q_R48   1.0
#define q_R53   2.0
#define q_R54   1.0
#define q_R55   1.0

// Washcoat Factor
#define Wash_F 70.0

/* Thiele Constants */
#define D_EFF_FIXED (3.40776e-7)
#define DELTA_W     (2.5e-5)
#define A_V         (2600)

#define NU_GAS_R7   1.0

/* ===== Helpers ===== */

/*
 * Calculate k_surface with coverage dependency.
 * Parameters:
 *  - A: pre-exponential factor.
 *  - beta: temperature exponent for Arrhenius form.
 *  - Ea_J_per_kmol: activation energy [J/kmol].
 *  - T: wall temperature [K].
 *  - idx_site: indices of coverage-dependent species in yi[].
 *  - mu: logarithmic coverage exponents for each site species.
 *  - eps_J_per_kmol: coverage-dependent energy terms [J/kmol].
 *  - Ns: number of coverage-dependent species.
 *  - yi: local species mass fractions (gas + surface).
 */
static inline real k_surface_covdep(real A, real beta, real Ea_J_per_kmol, real T,
                                    const int *idx_site, const real *mu, const real *eps_J_per_kmol,
                                    int Ns, const real yi[])
{
    const real Tpos  = MAX(EPS, T);
    const real invRT = 1.0 / (UNIVERSAL_GAS_CONSTANT * Tpos);
    real ln_k = log(MAX(EPS, A)) + beta*log(Tpos) - Ea_J_per_kmol*invRT;

    for (int i = 0; i < Ns; ++i) {
        const int  si   = idx_site ? idx_site[i] : 0;
        const real th   = MAX(1.0e-20, CLAMP01(yi[si]));
        const real mui  = mu  ? mu[i]  : 0.0;
        const real epsi = eps_J_per_kmol ? eps_J_per_kmol[i] : 0.0;
        ln_k += mui * log(th) + (epsi * th) * invRT;
    }
    return exp(ln_k);
}

/*
 * Calculate gas-phase molar concentration at the wall [kmol/m^3] from the wall mass fraction.
 * Manual-consistent conversion for DEFINE_SR_RATE() (see Fluent UDF Manual Example for SR_RATE):
 *   rho_w = rho_cell * (T_cell / T_w)
 *   Cg,k  = rho_w * Yk,w / MW_k
 * Parameters:
 *  - c0, t0 : adjacent cell and its thread (t->t0)
 *  - Tw     : wall temperature [K] (e.g., F_T(f,t))
 *  - yi_k   : wall mass fraction Yk,w [-] (from yi[] argument)
 *  - MW_k   : species molecular weight [kg/kmol]
 * Notes:
 *  - Uses absolute pressure (static + operating) and local temperature.
 *  - Uses species-specific gas constant (UNIVERSAL_GAS_CONSTANT / MW_k).
 */
static inline real gas_conc_cell(cell_t c0, Thread *t0, real Tw, real yi_k, real MW_k)
{
    const real T_cell = MAX(EPS, C_T(c0, t0));
    const real T_w = MAX(EPS, Tw);

    /* Scale cell density to the wall temperature at (approximately) constant pressure */
    const real rho_w  = C_R(c0, t0) * T_cell / T_w;
    /* kmol/m^3 */
    return (rho_w * yi_k) / MAX(EPS, MW_k);
}

/*
 * Thiele modulus and effectiveness factor calculation.
 * Parameters:
 *  - rpp: intrinsic surface rate [kmol/m2-s].
 *  - Cg: gas-phase concentration [kmol/m3].
 *  - nu_g: gas stoichiometric coefficient for diffusion scaling.
 *  - phi: (output) Thiele modulus.
 *  - eta: (output) effectiveness factor.
 */
static inline void thiele_eta_from_powerlaw(real rpp, real Cg, real nu_g, real *phi, real *eta)
{
    real kv_lin = 0.0;
    if (nu_g > 0.0 && Cg > EPS) kv_lin = A_V * rpp * (nu_g / Cg);
    const real ph = (kv_lin > 0.0) ? (DELTA_W * sqrt(kv_lin / D_EFF_FIXED)) : 0.0;
    const real et = (ph < 1.0e-6) ? (1.0 - ph*ph/3.0) : tanh(ph)/ph;
    *phi = ph;
    *eta = et;
}

/*
 * Sticking coefficient to pre-exponential factor A.
 * Parameters:
 *  - T: temperature [K].
 *  - S0: sticking coefficient.
 *  - M_kg_per_kmol: molecular weight [kg/kmol].
 *  - Gamma_kmol_m2: site density [kmol/m2].
 *  - q_site: number of occupied sites (reaction order in sites).
 */
static inline real k_sticking(real S0, real T, real M_kg_per_kmol, real Gamma_kmol_m2, real q_site)
{
    const real Tpos  = MAX(EPS, T);
    const real root = sqrt((UNIVERSAL_GAS_CONSTANT * Tpos) / MAX(EPS, 2.0 * M_PI * M_kg_per_kmol));
    const real invG = pow(1.0 / MAX(EPS, Gamma_kmol_m2), q_site);
    return S0 * root * invG;
}

/* ======================================================================= */
/* Reaction Constants Definitions */
/* ======================================================================= */

/* Define reusable parameters for coverage dependency */
static const int  idx_null[] = {0};
static const real val_null[] = {0.0};

/* Adsorption Parameters */
static const int  idx_site_r3[] = { IDX_Pt_Vac };
static const real mu_r3[]       = { -0.9 };
static const int  idx_site_r4[] = { IDX_Pt_Vac };
static const real mu_r4[]       = { -1.0 };
static const int  idx_site_r53[] = { IDX_Rh_Vac };
static const real mu_r53[]       = { -1.0 };

/* Pre-exponential Factors (Calculated from Sticking Coefficients)
 - R1: O2 + Pt(s) + Pt(s) -> O(s) + O(s)        
 - R2: C3H6 + Pt(s) + Pt(s) -> C3H6(s)         
 - R3: C3H6 + O(s) + Pt(s) -> C3H5(s) + OH(s) 
 - R4: H2 + Pt(s) + Pt(s) -> H(s) + H(s) 
 - R5: H2O + Pt(s) -> H2O(s) 
 - R6: CO2 + Pt(s) -> CO2(s) 
 - R7: CO + Pt(s) -> CO(s) 
 - R48: NO + Pt(s) -> NO(s)
 - R53: O2 + Rh(s) + Rh(s) -> O(Rh) + O(Rh) 
 - R54: CO + Rh(s) -> CO(Rh) 
 - R55: NO + Rh(s) -> NO(Rh) 
 - A1_k = 1.081678e+15 (based on Pt site density)
 - A2_k = 1.320535e+16 (based on Pt site density)
 - A3_k = 6.737423e+14 (based on Pt site density)
 - A4_k = 2.831952e+15 (based on Pt site density)
 - A5_k = 3.150917e+08 (based on Pt site density)
 - A6_k = 1.343976e+06 (based on Pt site density)
 - A7_k = 2.830189e+08 (based on Pt site density)
 - A48_k = 2.767013e+08 (based on Pt site density)
 - A53_k = 1.390729e+15 (based on Rh site density)
 - A54_k = 5.053909e+08 (based on Rh site density)
 - A55_k = 4.882963e+08 (based on Rh site density)                    
 */

/* --- DESORPTION CONSTANTS --- */
/* R8: O(s) + O(s) -> O2 + Pt(s) + Pt(s) */
static const int  idx_site_r8[] = {IDX_O_Pt};
static const real eps_r8[]      = {9.0E7};
#define NS_R8 1
#define A8_k 3.7E20
#define B8_beta 0.0
#define Ea8_Jpm 2.322E8
// CHECK

/* R9: C3H6(s) -> C3H6 + Pt(s) + Pt(s) */
#define NS_R9 0
#define A9_k 1E13
#define B9_beta 0.0
#define Ea9_Jpm 7.27E7
// CHECK

/* R10: C3H5(s) + OH(s) -> C3H6 + O(s) + Pt(s) */
#define NS_R10 0
#define A10_k 3.7E20
#define B10_beta 0.0
#define Ea10_Jpm 3.1E7
// CHECK

/* R11: H(s) + H(s) -> H2 + Pt(s) + Pt(s) */
static const int  idx_site_r11[] = {IDX_H_Pt};
static const real eps_r11[]      = {6.0E6};
#define NS_R11 1
#define A11_k 3.7E20
#define B11_beta 0.0
#define Ea11_Jpm 6.74E7
// CHECK

/* R12: H2O(s) -> H2O + Pt(s) */
#define NS_R12 0
#define A12_k 1.0E13
#define B12_beta 0.0
#define Ea12_Jpm 4.03E7
// CHECK

/* R13: CO(s) -> CO + Pt(s) */
static const int  idx_site_r13[] = {IDX_CO_Pt};
static const real eps_r13[]      = {3.3E7};
#define NS_R13 1
#define A13_k 1.0E13
#define B13_beta 0.0
#define Ea13_Jpm 1.364E8
// CHECK

/* R14: CO2(s) -> CO2 + Pt(s) */
#define NS_R14 0
#define A14_k 1.0E13
#define B14_beta 0.0
#define Ea14_Jpm 2.71E7
// CHECK

/* R49: NO(s) -> NO + Pt(s) */
#define A49_k 1.0E16
#define B49_beta 0.0
#define Ea49_Jpm 1.4E8
// CHECK

/* R50: N(s) + N(s) -> N2 + Pt(s) + Pt(s) */
static const int  idx_site_r50[] = {IDX_CO_Pt};
static const real eps_r50[]      = {7.5E7};
#define NS_R50 1
#define A50_k 3.7E20
#define B50_beta 0.0
#define Ea50_Jpm 1.139E8
// CHECK

/* R56: O(Rh) + O(Rh) -> O2 + Rh(s) + Rh(s) */
#define A56_k 3.0E20
#define B56_beta 0.0
#define Ea56_Jpm 2.933E8
// CHECK

/* R57: CO(Rh) -> CO + Rh(s) */
static const int  idx_site_r57[] = {IDX_CO_Rh, IDX_N_Rh};
static const real eps_r57[]      = {1.88E7, 4.19E7};
#define NS_R57 2
#define A57_k 1.0E14
#define B57_beta 0.0
#define Ea57_Jpm 1.323E8
// CHECK

/* R58: NO(Rh) -> NO + Rh(s) */
#define A58_k 5.0E13
#define B58_beta 0.0
#define Ea58_Jpm 1.089E8
// CHECK

/* R59: N(Rh) + N(Rh) -> N2 + Rh(s) + Rh(s) */
static const int  idx_site_r59[] = {IDX_N_Rh};
static const real eps_r59[]      = {1.67E7};
#define NS_R59 1
#define A59_k 1.11E18
#define B59_beta 0.0
#define Ea59_Jpm 1.369E8
// CHECK

/* --- SURFACE REACTION CONSTANTS --- */
/* R15: C3H5(s) + 5O(s) -> 5OH(s) + 3C(s) */
#define A15_k 3.7E16
#define Ea15_Jpm 9.5E7
// CHECK

/* R16: C3H6(s) + H(s) -> CC2H5(s) */
#define A16_k 1.0E13                            // Fixed. Reaction order=1.0, convert factor=1.0(no change)
#define Ea16_Jpm 7.54E7
// CHECK

/* R17: CC2H5(s) + H(s) -> C3H6(s) */
#define A17_k 3.7E20
#define Ea17_Jpm 4.88E7
// CHECK

/* R18: CC2H5(s) + Pt(s) -> C2H3(s) + CH2(s) */
#define A18_k 3.7E20
#define Ea18_Jpm 1.082E8
// CHECK

/* R19: C2H3(s) + CH2(s) -> Pt(s) + CC2H5(s) */
#define A19_k 3.7E20
#define Ea19_Jpm 3.2E6
// CHECK

/* R20: C2H3(s) + Pt(s) -> CH3(s) + C(s) */
#define A20_k 3.7E20
#define Ea20_Jpm 4.6E7
// CHECK

/* R21: CH3(s) + C(s) -> C2H3(s) + Pt(s) */
#define A21_k 3.7E20
#define Ea21_Jpm 4.69E7
// CHECK

/* R22: CH3(s) + Pt(s) -> CH2(s) + H(s) */
#define A22_k 1.26E21
#define Ea22_Jpm 7.04E7
// CHECK

/* R23: CH2(s) + H(s) -> CH3(s) + Pt(s) */
#define A23_k 3.09E21
#define Ea23_Jpm 0.0
// CHECK

/* R24: CH2(s) + Pt(s) -> CH(s) + H(s) */
#define A24_k 7.0E21
#define Ea24_Jpm 5.92E7
// CHECK

/* R25: CH(s) + H(s) -> CH2(s) + Pt(s) */
#define A25_k 3.09E21
#define Ea25_Jpm 0.0
// CHECK

/* R26: CH(s) + Pt(s) -> C(s) + H(s) */
#define A26_k 3.09E21
#define Ea26_Jpm 0.0
// CHECK

/* R27: C(s) + H(s) -> CH(s) + Pt(s) */
#define A27_k 1.25E21
#define Ea27_Jpm 1.38E8
// CHECK

/* R28: C2H3(s) + O(s) -> Pt(s) + CH3CO(s) */
#define A28_k 3.7E18
#define Ea28_Jpm 6.23E7
// CHECK

/* R29: CH3CO(s) + Pt(s) -> C2H3(s) + O(s) */
static const int  idx_site_r29[] = {IDX_O_Pt};
static const real eps_r29[]      = {-4.5E7};
#define NS_R29 1
#define A29_k 3.7E20
#define Ea29_Jpm 1.967E8
// CHECK

/* R30: CH3(s) + CO(s) -> Pt(s) + CH3CO(s) */
#define A30_k 3.7E20
#define Ea30_Jpm 8.29E7
// CHECK

/* R31: CH3CO(s) + Pt(s) -> CH3(s) + CO(s) */
#define A31_k 3.7E20
#define Ea31_Jpm 0.0
// CHECK

/* R32: CH3(s) + O(s) -> CH2(s) + OH(s) */
#define A32_k 3.7E20
#define Ea32_Jpm 3.66E7
// CHECK

/* R33: CH2(s) + OH(s) -> CH3(s) + O(s) */
#define A33_k 3.7E20
#define Ea33_Jpm 2.51E7
// CHECK

/* R34: CH2(s) + O(s) -> CH(s) + OH(s) */
#define A34_k 3.7E20
#define Ea34_Jpm 2.51E7
// CHECK

/* R35: CH(s) + OH(s) -> CH2(s) + O(s) */
#define A35_k 3.7E20
#define Ea35_Jpm 2.52E7
// CHECK

/* R36: CH(s) + O(s) -> C(s) + OH(s) */
#define A36_k 3.7E20
#define Ea36_Jpm 2.51E7
// CHECK

/* R37: C(s) + OH(s) -> CH(s) + O(s) */
#define A37_k 3.7E20
#define Ea37_Jpm 2.248E8
// CHECK

/* R38: O(s) + H(s) -> OH(s) + Pt(s) */
#define A38_k 3.7E20
#define Ea38_Jpm 1.15E7
// CHECK

/* R39: OH(s) + Pt(s) -> O(s) + H(s) */
#define A39_k 5.77E21
#define Ea39_Jpm 7.49E7
// CHECK

/* R40: H(s) + OH(s) -> H2O(s) + Pt(s) */
#define A40_k 3.7E20
#define Ea40_Jpm 1.74E7
// CHECK

/* R41: H2O(s) + Pt(s) -> H(s) + OH(s) */
#define A41_k 3.66E20
#define Ea41_Jpm 7.36E7
// CHECK

/* R42: OH(s) + OH(s) -> H2O(s) + O(s) */
#define A42_k 3.7E20
#define Ea42_Jpm 4.82E7
// CHECK

/* R43: H2O(s) + O(s) -> OH(s) + OH(s) */
#define A43_k 2.35E19
#define Ea43_Jpm 4.1E7
// CHECK

/* R44: CO(s) + O(s) -> CO2(s) + Pt(s) */
static const int  idx_site_r44[] = {IDX_CO_Pt, IDX_NO_Pt};
static const real eps_r44[]      = {3.3E7, -9.0E7};
#define NS_R44 2
#define A44_k 3.7E19
#define Ea44_Jpm 1.08E8
// CHECK

/* R45: CO2(s) + Pt(s) -> CO(s) + O(s) */
static const int  idx_site_r45[] = {IDX_O_Pt};
static const real eps_r45[]      = {-4.5E7};
#define NS_R45 1
#define A45_k 3.7E20
#define Ea45_Jpm 1.651E8
// CHECK

/* R46: C(s) + O(s) -> CO(s) + Pt(s) */
static const int  idx_site_r46[] = {IDX_CO_Pt};
static const real eps_r46[]      = {-3.3E7};
#define NS_R46 1
#define A46_k 3.7E20
#define Ea46_Jpm 0.0
// CHECK

/* R47: CO(s) + Pt(s) -> C(s) + O(s) */
static const int  idx_site_r47[] = {IDX_O_Pt};
static const real eps_r47[]      = {-4.5E7};
#define NS_R47 1
#define A47_k 3.7E20
#define Ea47_Jpm 2.185E8
// CHECK

/* R51: NO(s) + Pt(s) -> N(s) + O(s) */
static const int  idx_site_r51[] = {IDX_CO_Pt};
static const real eps_r51[]      = {-3.0E6};
#define NS_R51 1
#define A51_k 5.0E19
#define Ea51_Jpm 1.078E8
// CHECK

/* R52: N(s) + O(s) -> NO(s) + Pt(s) */
static const int  idx_site_r52[] = {IDX_O_Pt};
static const real eps_r52[]      = {4.5E7};
#define NS_R52 1
#define A52_k 3.7E20
#define Ea52_Jpm 1.281E8
// CHECK

/* R60: CO(Rh) + O(Rh) -> CO2 + Rh(s) + Rh(s) */
#define A60_k 3.7E19
#define Ea60_Jpm 5.99E7
// CHECK

/* R61: NO(Rh) + Rh(s) -> N(Rh) + O(Rh) */
#define A61_k 2.22E21
#define Ea61_Jpm 7.95E7
// CHECK


/* Adsorption rate constant. This wrapper exists to keep adsorption logic
 * visually separated from desorption / surface reactions.
 */
static inline real k_surface_ads(real A, real beta, real Ea_J_per_kmol, real T,
                                 const int *idx_site, const real *mu,
                                 const real *eps_J_per_kmol, int Ns,
                                 const real yi[])
{
    return k_surface_covdep(A, beta, Ea_J_per_kmol, T,
                            idx_site, mu, eps_J_per_kmol, Ns, yi);
}


/* ======================================================================= */
/* Helper: Reaction 7 Base Rate & Eta Calculation (The Reference)          */
/* Parameters:                                                             */
/*  - c0, t0: cell and thread handles.                                     */
/*  - Tw: wall temperature [K].                                            */
/*  - yi: local species mass fractions.                                    */
/*  - Cv_R7: available Pt site concentration for reaction 7 [kmol/m2].     */
/*  - rate7_base: (output) intrinsic base rate r7'' [kmol/m2-s].            */
/*  - phi7: (output) Thiele modulus for reaction 7.                        */
/*  - eta7: (output) effectiveness factor for reaction 7.                  */
/* ======================================================================= */
static inline void reaction7_base_and_eta(cell_t c0, Thread *t0,
                                          real Tw, const real *yi,
                                          real Cv_R7,
                                          real *rate7_base,
                                          real *phi7, real *eta7)
{
    /* Based on Reaction-7 (CO adsorption) */
    const real k7  = k_sticking(S0_CO, Tw, MW_CO, SITE_DEN_TOT, q_R7);
    const real Cg7 = gas_conc_cell(c0, t0, Tw, yi[IDX_CO], MW_CO);
    const real r7  = k7 * Cg7 * Cv_R7;      /* base rate r7'' */

    real phi, eta;
    thiele_eta_from_powerlaw(r7, Cg7, NU_GAS_R7, &phi, &eta);

    if (rate7_base) *rate7_base = r7;
    if (phi7)       *phi7       = phi;
    if (eta7)       *eta7       = eta;
}

/* ======================================================================= */
/* Adsorption-only helpers (explicit, Reaction-8 style)                     */
/* ======================================================================= */





/* Desorption rate constant wrapper (for readability) */
static inline real k_surface_des(real A, real beta, real Ea_J_per_kmol, real T,
                                 const int *idx_site, const real *mu,
                                 const real *eps_J_per_kmol, int Ns,
                                 const real yi[])
{
    return k_surface_covdep(A, beta, Ea_J_per_kmol, T,
                            idx_site, mu, eps_J_per_kmol, Ns, yi);
}

/* Surface-reaction rate constant wrapper (for readability) */
static inline real k_surface_sr(real A, real beta, real Ea_J_per_kmol, real T,
                                const int *idx_site, const real *mu,
                                const real *eps_J_per_kmol, int Ns,
                                const real yi[])
{
    return k_surface_covdep(A, beta, Ea_J_per_kmol, T,
                            idx_site, mu, eps_J_per_kmol, Ns, yi);
}

/* Convert a surface coverage theta_i to a surface concentration [kmol/m2]
 * on Pt or Rh site pools.
 */
static inline real Cs_Pt(const real yi[], int idx_theta)
{
    return SITE_DEN_Pt * CLAMP01(yi[idx_theta]);
}
static inline real Cs_Rh(const real yi[], int idx_theta)
{
    return SITE_DEN_Rh * CLAMP01(yi[idx_theta]);
}
/* Return effectiveness factor based on the reference (Reaction-7) Thiele model.
 * Optionally returns the reference base rate and phi.
 */
static inline real eta_from_reaction7(cell_t c0, Thread *t0,
                                      real Tw, const real *yi,
                                      real Cv_R7,
                                      real *r7_base_out,
                                      real *phi7_out)
{
    real r7_base, phi7, eta7;
    reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7, &r7_base, &phi7, &eta7);

    if (r7_base_out) *r7_base_out = r7_base;
    if (phi7_out)    *phi7_out    = phi7;

    /* Optionally override computed effectiveness factor (force 1.0 if desired).
       To enable diffusion-limited correction use the computed `eta7` from above.
       To force fast-diffusion assumption, assign eta7 = 1.0. */

    //eta7 = 1.0; /* no diffusion limitation (override). */

    return eta7;
}

/* UDM write helper for the reference reaction (Reaction-7). */
static inline void write_udm_r7(cell_t c0, Thread *t0,
                                real r7_base, real Cg7,
                                real phi7, real eta7)
{
    if (Cg7 > EPS) {
        C_UDMI(c0,t0, UDM_RCO_NET)    = r7_base;
        C_UDMI(c0,t0, UDM_RCO_NETP)   = 0;
        C_UDMI(c0,t0, UDM_DRDC)       = r7_base / Cg7;
        C_UDMI(c0,t0, UDM_KAPP)       = A_V * r7_base / Cg7;
        C_UDMI(c0,t0, UDM_PHI)        = phi7;
        C_UDMI(c0,t0, UDM_ETA)        = eta7;
    } else {
        C_UDMI(c0,t0, UDM_RCO_NET)    = 0;
        C_UDMI(c0,t0, UDM_RCO_NETP)   = 0;
        C_UDMI(c0,t0, UDM_DRDC)       = 0;
        C_UDMI(c0,t0, UDM_KAPP)       = 0;
        C_UDMI(c0,t0, UDM_PHI)        = 0;
        C_UDMI(c0,t0, UDM_ETA)        = 0;
    }
}

/* One-time logger for adsorption reaction rates. */
static inline void log_ads_rate_once(int r_id_num, const char* r_short_name,
                                     real rate_base, real eta)
{
    if (r_id_num < NUM_REACTIONS && print_gate[r_id_num] == 0) {
        #if RP_HOST
            Message("\n[%s] rate_base = %.6e | rate_wash = %.6e | rate_eta = %.6e (eta=%.3e)\n",
                    r_short_name, rate_base, rate_base * Wash_F,
                    rate_base * Wash_F * eta, eta);
        #endif
        #if RP_NODE
            if (myid == 0) {
                Message0("\n[%s] rate_base = %.6e | rate_wash = %.6e | rate_eta = %.6e (eta=%.3e)\n",
                         r_short_name, rate_base, rate_base * Wash_F,
                         rate_base * Wash_F * eta, eta);
            }
        #endif
        print_gate[r_id_num] = 1;
    }
}


/**************************************************/
/* ======================================================================= */
/* Main UDF                                                                */
/* ======================================================================= */
DEFINE_SR_RATE(chatterjee_pt_ads_des_flat_internalized, f, fthread, r, mw, yi, rr)
{
    cell_t  c0  = F_C0(f, fthread);
    Thread *t0  = THREAD_T0(fthread);
    const real Tw  = F_T(f, fthread);

    /* Vacant Site Concentrations */
    const real theta_pt_vac = CLAMP01(yi[IDX_Pt_Vac]);
    const real theta_rh_vac = CLAMP01(yi[IDX_Rh_Vac]);
    const real theta_O_Pt = CLAMP01(yi[IDX_O_Pt]);

    /* Pre-calculate site terms */
    const real term_pt = MAX(EPS, SITE_DEN_Pt * theta_pt_vac);
    const real term_rh = MAX(EPS, SITE_DEN_Rh * theta_rh_vac);
    const real term_O_Pt = MAX(EPS, SITE_DEN_Pt * theta_O_Pt);

    /* Cv for Reaction 7 (Ref) */
    const real Cv_R7_val = term_pt;

    // *rr = 0.0; // if set = 0.0, no reaction will occur

     /* Precompute reference Reaction-7 (CO adsorption) effectiveness and base rate
         so we don't recompute it in every reaction branch. */
     real r7_base_ref = 0.0;
     real phi7_ref = 0.0;
     const real eta_ref = eta_from_reaction7(c0, t0, Tw, yi, Cv_R7_val, &r7_base_ref, &phi7_ref);

    /* --- REACTION DISPATCHER (Continuous else if chain) --- */

/* ======================================================================= */
/* Adsorption                                                              */
/* ======================================================================= */
    if (STREQ(r->name, "reaction-1")) {
        /* R1: O2 + Pt(s) + Pt(s) -> O(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_ads  = k_sticking(S0_O2, Tw, MW_O2, SITE_DEN_TOT, q_R1);
        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_O2], MW_O2);
        const real Cv_ads = term_pt * term_pt;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(1, "r1", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-2")) {
        /* R2: C3H6 + Pt(s) + Pt(s) -> C3H6(s) */
        const real eta7 = eta_ref;
        const real k_ads  = k_sticking(S0_C3H6, Tw, MW_C3H6, SITE_DEN_TOT, q_R2);
        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_C3H6], MW_C3H6);
        const real Cv_ads = term_pt * term_pt;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(2, "r2", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-3")) {
        /* R3: C3H6 + O(s) + Pt(s) -> C3H5(s) + OH(s), cov-dependency */
        const real eta7 = eta_ref;
        const real k0   = k_sticking(S0_C3H6_O, Tw, MW_C3H6, SITE_DEN_TOT, q_R3); 
        const real cov  = pow(MAX(1.0E-20, theta_pt_vac), -0.9); /* mu_r3 = -0.9*/
        const real k_ads  = k0*cov;

        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_C3H6], MW_C3H6);
        const real Cv_ads = term_pt * term_O_Pt;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(3, "r3", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-4")) {
        /* R4: H2 + Pt(s) + Pt(s) -> H(s) + H(s), cov-dependency */
        const real eta7 = eta_ref;
        const real k0   = k_sticking(S0_H2, Tw, MW_H2, SITE_DEN_TOT, q_R4); 
        const real cov  = pow(MAX(1.0E-20, theta_pt_vac), -1.0); /* mu_r4 = -1.0*/
        const real k_ads  = k0 * cov;

        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_H2], MW_H2);
        const real Cv_ads = term_pt * term_pt;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(4, "r4", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-5")) {
        /* R5: H2O + Pt(s) -> H2O(s) */
        const real eta7 = eta_ref;
        const real k_ads  = k_sticking(S0_H2O, Tw, MW_H2O, SITE_DEN_TOT, q_R5);
        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_H2O], MW_H2O);
        const real Cv_ads = term_pt;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(5, "r5", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-6")) {
        /* R6: CO2 + Pt(s) -> CO2(s) */
        const real eta7 = eta_ref;
        const real k_ads  = k_sticking(S0_CO2, Tw, MW_CO2, SITE_DEN_TOT, q_R6);
        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_CO2], MW_CO2);
        const real Cv_ads = term_pt;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(6, "r6", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-7")) {
        /* R7: CO + Pt(s) -> CO(s) (reference for eta) */
        const real eta7 = eta_ref;
        const real Cg7  = gas_conc_cell(c0, t0, Tw, yi[IDX_CO], MW_CO);
        write_udm_r7(c0, t0, r7_base_ref, Cg7, phi7_ref, eta_ref);
        *rr = r7_base_ref * Wash_F * eta_ref;
        log_ads_rate_once(7, "r7", r7_base_ref, eta_ref);
    }
    else if (STREQ(r->name, "reaction-48")) {
        /* R48: NO + Pt(s) -> NO(s) */
        const real eta7 = eta_ref;
        const real k_ads  = k_sticking(S0_NO, Tw, MW_NO, SITE_DEN_TOT, q_R48);
        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_NO], MW_NO);
        const real Cv_ads = term_pt;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(48, "r48", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-53")) {
        /* R53: O2 + Rh(s) + Rh(s) -> O(Rh) + O(Rh), cov-dependency */
        const real eta7 = eta_ref;
        const real k0   = k_sticking(S1_O2, Tw, MW_O2, SITE_DEN_TOT, q_R53); 
        const real cov  = pow(MAX(1.0E-20, theta_rh_vac), -1.0); /* mu_r53 = -1.0*/
        const real k_ads  = k0 * cov;

        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_O2], MW_O2);
        const real Cv_ads = term_rh * term_rh;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(53, "r53", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-54")) {
        /* R54: CO + Rh(s) -> CO(Rh) */
        const real eta7 = eta_ref;
        const real k_ads  = k_sticking(S1_CO, Tw, MW_CO, SITE_DEN_TOT, q_R54);
        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_CO], MW_CO);
        const real Cv_ads = term_rh;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(54, "r54", rate_base, eta7);
    }
    else if (STREQ(r->name, "reaction-55")) {
        /* R55: NO + Rh(s) -> NO(Rh) */
        const real eta7 = eta_ref;
        const real k_ads  = k_sticking(S1_NO, Tw, MW_NO, SITE_DEN_TOT, q_R55);
        const real Cg     = gas_conc_cell(c0, t0, Tw, yi[IDX_NO], MW_NO);
        const real Cv_ads = term_rh;
        const real rate_base = k_ads * Cg * Cv_ads;

        *rr = rate_base * Wash_F * eta7;
        log_ads_rate_once(55, "r55", rate_base, eta7);
    }

/* ======================================================================= */
/* Desorption                                                              */
/* ======================================================================= */
    else if (STREQ(r->name, "reaction-8")) {

        /* R8: 2O(s) -> O2 (Pt), cov-dependency */
        const real eta7   = eta_ref;
        const real k_des  = k_surface_des(A8_k, B8_beta, Ea8_Jpm, Tw, idx_site_r8, val_null, eps_r8, NS_R8, yi);
        const real C_O    = Cs_Pt(yi, IDX_O_Pt);
        const real Cv_des = C_O * C_O;
        const real rate_base = k_des * Cv_des;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-9")) {

        /* R9: C3H6(s) -> C3H6 (Pt) */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A9_k, B9_beta, Ea9_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real C_C3H6 = Cs_Pt(yi, IDX_C3H6_Pt);
        const real rate_base = k_des * C_C3H6;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-10")) {
        /* R10: C3H5(s) + OH(s) -> C3H6 + O(s) + Pt(s)*/
        const real eta7 = eta_ref;
        const real k10 = k_surface_covdep(A10_k, B10_beta, Ea10_Jpm, Tw, idx_null, val_null, val_null, NS_R10, yi);
        const real thC3H5 = CLAMP01(yi[IDX_C3H5_Pt]);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k10 * thC3H5 * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-11")) {

        /* R11: 2H(s) -> H2 (Pt), cov-dependency */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A11_k, B11_beta, Ea11_Jpm, Tw, idx_site_r11, val_null, eps_r11, NS_R11, yi);
        const real C_H   = Cs_Pt(yi, IDX_H_Pt);
        const real Cv_des = C_H * C_H;
        const real rate_base = k_des * Cv_des;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-12")) {

        /* R12: H2O(s) -> H2O (Pt) */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A12_k, B12_beta, Ea12_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real C_H2O = Cs_Pt(yi, IDX_H2O_Pt);
        const real rate_base = k_des * C_H2O;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-13")) {

        /* R13: CO(s) -> CO (Pt), cov-dependency */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A13_k, B13_beta, Ea13_Jpm, Tw, idx_site_r13, val_null, eps_r13, NS_R13, yi);
        const real C_CO  = Cs_Pt(yi, IDX_CO_Pt);
        const real rate_base = k_des * C_CO;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-14")) {

        /* R14: CO2(s) -> CO2 (Pt) */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A14_k, B14_beta, Ea14_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real C_CO2 = Cs_Pt(yi, IDX_CO2_Pt);
        const real rate_base = k_des * C_CO2;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-49")) {

        /* R49: NO(s) -> NO (Pt) */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A49_k, B49_beta, Ea49_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real C_NO  = Cs_Pt(yi, IDX_NO_Pt);
        const real rate_base = k_des * C_NO;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-50")) {

        /* R50: 2N(s) -> N2 (Pt), cov-dependency */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A50_k, B50_beta, Ea50_Jpm, Tw, idx_site_r50, val_null, eps_r50, NS_R50, yi);
        const real C_N   = Cs_Pt(yi, IDX_N_Pt);
        const real Cv_des = C_N * C_N;
        const real rate_base = k_des * Cv_des;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-56")) {

        /* R56: 2O(Rh) -> O2 */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A56_k, B56_beta, Ea56_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real C_O_rh = Cs_Rh(yi, IDX_O_Rh);
        const real Cv_des = C_O_rh * C_O_rh;
        const real rate_base = k_des * Cv_des;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-57")) {

        /* R57: CO(Rh) -> CO, cov-dependency */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A57_k, B57_beta, Ea57_Jpm, Tw, idx_site_r57, val_null, eps_r57, NS_R57, yi);
        const real C_CO_rh = Cs_Rh(yi, IDX_CO_Rh);
        const real rate_base = k_des * C_CO_rh;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-58")) {

        /* R58: NO(Rh) -> NO */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A58_k, B58_beta, Ea58_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real C_NO_rh = Cs_Rh(yi, IDX_NO_Rh);
        const real rate_base = k_des * C_NO_rh;

        *rr = rate_base * Wash_F * eta7;
    }
    else if (STREQ(r->name, "reaction-59")) {

        /* R59: 2N(Rh) -> N2, cov-dependency */
        const real eta7  = eta_ref;
        const real k_des = k_surface_des(A59_k, B59_beta, Ea59_Jpm, Tw, idx_site_r59, val_null, eps_r59, NS_R59, yi);
        const real C_N_rh = Cs_Rh(yi, IDX_N_Rh);
        const real Cv_des = C_N_rh * C_N_rh;
        const real rate_base = k_des * Cv_des;

        *rr = rate_base * Wash_F * eta7;
}

/* ======================================================================= */
/* Surface Reactions (LH / ER)                                             */
/* - Pt oxidation / NOx network: R15-R52                                   */
/* - Rh network: R60-R61 (Rh adsorption/desorption handled above)          */
/* ======================================================================= */
    else if (STREQ(r->name, "reaction-15")) {

        /* R15: C3H5(s) + 5O(s) -> 5OH(s) + 3C(s) */
        const real eta7 = eta_ref;
        const real k_sr  = k_surface_sr(A15_k, 0.0, Ea15_Jpm, Tw, idx_null, val_null, val_null, 0, yi);

        const real C_C3H5 = Cs_Pt(yi, IDX_C3H5_Pt);         /* Cs_Pt = SITE_DEN_Pt * yi[] */
        const real C_O    = MAX(EPS, Cs_Pt(yi, IDX_O_Pt));  /* legacy MAX(EPS, ...) for pow(.,5) */
        const real Cv_ads = 1.0;              /* term_pt * term_pt */
        const real rate_base = k_sr * C_C3H5 * pow(C_O, 5.0) * Cv_ads;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-16")) {

        /* R16: C3H6(s) -> H(s) + CC2H5(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A16_k, 0.0, Ea16_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_C3H6 = Cs_Pt(yi, IDX_C3H6_Pt);

        const real rate_base = k_sr * C_C3H6;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-17")) {

        /* R17: CC2H5(s) + H(s) -> C3H6(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A17_k, 0.0, Ea17_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CC2H5 = Cs_Pt(yi, IDX_CC2H5_Pt);
        const real C_H = Cs_Pt(yi, IDX_H_Pt);

        const real rate_base = k_sr * C_CC2H5 * C_H;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-18")) {

        /* R18: CC2H5(s) + Pt(s) -> C2H3(s) + CH2(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A18_k, 0.0, Ea18_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CC2H5 = Cs_Pt(yi, IDX_CC2H5_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));
        // SAME
        // term_pt = SITE_DEN_Pt * theta_pt_vac{=yi[IDX_Pt_Vac]} 
        // C_vac   = Cs_Pt(yi, IDX_Pt_Vac){=SITE_DEN_Pt * yi[IDX_Pt_Vac]}

        const real rate_base = k_sr * C_CC2H5 * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-19")) {

        /* R19: C2H3(s) + CH2(s) -> Pt(s) + CC2H5(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A19_k, 0.0, Ea19_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_C2H3 = Cs_Pt(yi, IDX_C2H3_Pt);
        const real C_CH2 = Cs_Pt(yi, IDX_CH2_Pt);

        const real rate_base = k_sr * C_C2H3 * C_CH2;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-20")) {

        /* R20: C2H3(s) + Pt(s) -> CH3(s) + C(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A20_k, 0.0, Ea20_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_C2H3 = Cs_Pt(yi, IDX_C2H3_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_C2H3 * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-21")) {

        /* R21: CH3(s) + C(s) -> C2H3(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A21_k, 0.0, Ea21_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH3 = Cs_Pt(yi, IDX_CH3_Pt);
        const real C_C = Cs_Pt(yi, IDX_C_Pt);

        const real rate_base = k_sr * C_CH3 * C_C;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-22")) {

        /* R22: CH3(s) + Pt(s) -> CH2(s) + H(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A22_k, 0.0, Ea22_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH3 = Cs_Pt(yi, IDX_CH3_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_CH3 * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-23")) {

        /* R23: CH2(s) + H(s) -> CH3(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A23_k, 0.0, Ea23_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH2 = Cs_Pt(yi, IDX_CH2_Pt);
        const real C_H = Cs_Pt(yi, IDX_H_Pt);

        const real rate_base = k_sr * C_CH2 * C_H;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-24")) {

        /* R24: CH2(s) + Pt(s) -> CH(s) + H(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A24_k, 0.0, Ea24_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH2 = Cs_Pt(yi, IDX_CH2_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_CH2 * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-25")) {

        /* R25: CH(s) + H(s) -> CH2(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A25_k, 0.0, Ea25_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH = Cs_Pt(yi, IDX_CH_Pt);
        const real C_H = Cs_Pt(yi, IDX_H_Pt);

        const real rate_base = k_sr * C_CH * C_H;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-26")) {

        /* R26: CH(s) + Pt(s) -> C(s) + H(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A26_k, 0.0, Ea26_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH = Cs_Pt(yi, IDX_CH_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_CH * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-27")) {

        /* R27: C(s) + H(s) -> CH(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A27_k, 0.0, Ea27_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_C = Cs_Pt(yi, IDX_C_Pt);
        const real C_H = Cs_Pt(yi, IDX_H_Pt);

        const real rate_base = k_sr * C_C * C_H;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-28")) {

        /* R28: C2H3(s) + O(s) -> Pt(s) + CH3CO(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A28_k, 0.0, Ea28_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_C2H3 = Cs_Pt(yi, IDX_C2H3_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_C2H3 * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-29")) {

        /* R29: CH3CO(s) + Pt(s) -> C2H3(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A29_k, 0.0, Ea29_Jpm, Tw,
                                       idx_site_r29, val_null, eps_r29, NS_R29, yi);

        const real C_CH3CO = Cs_Pt(yi, IDX_CH3CO_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_CH3CO * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-30")) {

        /* R30: CH3(s) + CO(s) -> Pt(s) + CH3CO(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A30_k, 0.0, Ea30_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH3 = Cs_Pt(yi, IDX_CH3_Pt);
        const real C_CO = Cs_Pt(yi, IDX_CO_Pt);

        const real rate_base = k_sr * C_CH3 * C_CO;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-31")) {

        /* R31: CH3CO(s) + Pt(s) -> CH3(s) + CO(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A31_k, 0.0, Ea31_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH3CO = Cs_Pt(yi, IDX_CH3CO_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_CH3CO * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-32")) {

        /* R32: CH3(s) + O(s) -> CH2(s) + OH(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A32_k, 0.0, Ea32_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH3 = Cs_Pt(yi, IDX_CH3_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_CH3 * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-33")) {

        /* R33: CH2(s) + OH(s) -> CH3(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A33_k, 0.0, Ea33_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH2 = Cs_Pt(yi, IDX_CH2_Pt);
        const real C_OH = Cs_Pt(yi, IDX_OH_Pt);

        const real rate_base = k_sr * C_CH2 * C_OH;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-34")) {

        /* R34: CH2(s) + O(s) -> CH(s) + OH(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A34_k, 0.0, Ea34_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH2 = Cs_Pt(yi, IDX_CH2_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_CH2 * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-35")) {

        /* R35: CH(s) + OH(s) -> CH2(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A35_k, 0.0, Ea35_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH = Cs_Pt(yi, IDX_CH_Pt);
        const real C_OH = Cs_Pt(yi, IDX_OH_Pt);

        const real rate_base = k_sr * C_CH * C_OH;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-36")) {

        /* R36: CH(s) + O(s) -> C(s) + OH(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A36_k, 0.0, Ea36_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_CH = Cs_Pt(yi, IDX_CH_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_CH * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-37")) {

        /* R37: C(s) + OH(s) -> CH(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A37_k, 0.0, Ea37_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_C = Cs_Pt(yi, IDX_C_Pt);
        const real C_OH = Cs_Pt(yi, IDX_OH_Pt);

        const real rate_base = k_sr * C_C * C_OH;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-38")) {

        /* R38: O(s) + H(s) -> OH(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A38_k, 0.0, Ea38_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_O = Cs_Pt(yi, IDX_O_Pt);
        const real C_H = Cs_Pt(yi, IDX_H_Pt);

        const real rate_base = k_sr * C_O * C_H;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-39")) {

        /* R39: OH(s) + Pt(s) -> O(s) + H(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A39_k, 0.0, Ea39_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_OH = Cs_Pt(yi, IDX_OH_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_OH * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-40")) {

        /* R40: H(s) + OH(s) -> H2O(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A40_k, 0.0, Ea40_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_H = Cs_Pt(yi, IDX_H_Pt);
        const real C_OH = Cs_Pt(yi, IDX_OH_Pt);

        const real rate_base = k_sr * C_H * C_OH;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-41")) {

        /* R41: H2O(s) + Pt(s) -> H(s) + OH(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A41_k, 0.0, Ea41_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_H2O = Cs_Pt(yi, IDX_H2O_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_H2O * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-42")) {

        /* R42: OH(s) + OH(s) -> H2O(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A42_k, 0.0, Ea42_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_OH = Cs_Pt(yi, IDX_OH_Pt);

        const real rate_base = k_sr * C_OH * C_OH;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-43")) {

        /* R43: H2O(s) + O(s) -> OH(s) + OH(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A43_k, 0.0, Ea43_Jpm, Tw,
                                       idx_null, val_null, val_null, 0, yi);

        const real C_H2O = Cs_Pt(yi, IDX_H2O_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_H2O * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-44")) {

        /* R44: CO(s) + O(s) -> CO2(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A44_k, 0.0, Ea44_Jpm, Tw,
                                       idx_site_r44, val_null, eps_r44, NS_R44, yi);

        const real C_CO = Cs_Pt(yi, IDX_CO_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_CO * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-45")) {

        /* R45: CO2(s) + Pt(s) -> CO(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A45_k, 0.0, Ea45_Jpm, Tw,
                                       idx_site_r45, val_null, eps_r45, NS_R45, yi);

        const real C_CO2 = Cs_Pt(yi, IDX_CO2_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_CO2 * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-46")) {

        /* R46: C(s) + O(s) -> CO(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A46_k, 0.0, Ea46_Jpm, Tw,
                                       idx_site_r46, val_null, eps_r46, NS_R46, yi);

        const real C_C = Cs_Pt(yi, IDX_C_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_C * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-47")) {

        /* R47: CO(s) + Pt(s) -> C(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A47_k, 0.0, Ea47_Jpm, Tw,
                                       idx_site_r47, val_null, eps_r47, NS_R47, yi);

        const real C_CO = Cs_Pt(yi, IDX_CO_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_CO * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-51")) {

        /* R51: NO(s) + Pt(s) -> N(s) + O(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A51_k, 0.0, Ea51_Jpm, Tw,
                                       idx_site_r51, val_null, eps_r51, NS_R51, yi);

        const real C_NO = Cs_Pt(yi, IDX_NO_Pt);
        const real C_vac = MAX(EPS, Cs_Pt(yi, IDX_Pt_Vac));

        const real rate_base = k_sr * C_NO * C_vac;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-52")) {

        /* R52: N(s) + O(s) -> NO(s) + Pt(s) */
        const real eta7 = eta_ref;
        const real k_sr = k_surface_sr(A52_k, 0.0, Ea52_Jpm, Tw,
                                       idx_site_r52, val_null, eps_r52, NS_R52, yi);

        const real C_N = Cs_Pt(yi, IDX_N_Pt);
        const real C_O = Cs_Pt(yi, IDX_O_Pt);

        const real rate_base = k_sr * C_N * C_O;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-60")) {

        /* R60: CO(Rh) + O(Rh) -> CO2 + Rh(s) + Rh(s) */
        const real eta7 = eta_ref;
        const real k_sr  = k_surface_sr(A60_k, 0.0, Ea60_Jpm, Tw, idx_null, val_null, val_null, 0, yi);

        const real C_CO_rh = Cs_Rh(yi, IDX_CO_Rh);
        const real C_O_rh  = Cs_Rh(yi, IDX_O_Rh);

        const real rate_base = k_sr * C_CO_rh * C_O_rh;

        *rr = rate_base * Wash_F;
        // CHECK
    }
    else if (STREQ(r->name, "reaction-61")) {

        /* R61: NO(Rh) + Rh(s) -> N(Rh) + O(Rh) */
        const real eta7 = eta_ref;
        const real k_sr  = k_surface_sr(A61_k, 0.0, Ea61_Jpm, Tw, idx_null, val_null, val_null, 0, yi);

        const real C_NO_rh  = Cs_Rh(yi, IDX_NO_Rh);
        const real C_vac_rh = Cs_Rh(yi, IDX_Rh_Vac);

        const real rate_base = k_sr * C_NO_rh * C_vac_rh;

        *rr = rate_base * Wash_F;
        // CHECK
    }
}
