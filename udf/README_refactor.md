# Fluent UDF Chemistry Refactor (State-Centered + Solver-Centered)

## 1) Legacy monolithic summary and Wash_F/legacy surface-only washcoat factor injection
The legacy file `B_Routine_SR_rev36v4_LocalF_Eta_R15Fix_DSR_ReCheck_GamTot_CUR.c` is a single `DEFINE_SR_RATE` callback with a long `if/else if` chain (`reaction-1` ... `reaction-61`).

Observed scaling split in legacy:
- Several branches used `*rr = rate_base * Wash_F * eta7`.
- Other branches used `*rr = rate_base * legacy surface-only washcoat factor`.

This refactor removes `legacy surface-only washcoat factor` from the new architecture and centralizes export scaling.

## 2) New full-conversion architecture
Files:
- `chem_indices.h`: index map and reaction count.
- `chem_constants.h`: constants + behavior flags.
- `chem_state.h/.c`: chemistry state, initialization, projection.
- `chem_rates.h/.c`: intrinsic elementary rates and export scaling.
- `chem_rhs.h/.c`: ODE RHS from intrinsic rates.
- `chem_jacobian.h/.c`: Jacobian assembly (numerical fallback by default).
- `stiff_solver.h/.c`: solver interface and backward-Euler/Euler placeholder backend.
- `chemistry_cache.h/.c`: Fluent UDM cache read/write helpers.
- `fluent_hooks.c`: `DEFINE_INIT`, `DEFINE_ADJUST`, compatibility `DEFINE_SR_RATE` wrapper.

## 3) Data flow
1. `DEFINE_INIT`: allocates initial chemistry state in UDM slots.
2. `DEFINE_ADJUST`: wall-face loop -> fetch local gas/wall state -> load cached chemistry -> integrate surface state ODE -> compute intrinsic rates -> apply eta/washcoat in export layer -> write cache.
3. `DEFINE_SR_RATE` compatibility hook (`chatterjee_cached_rates`): maps `reaction-N` to cached exported rate only.

No integration is performed in return hooks.

## 4) Washcoat scaling policy (mandatory behavior)
- Intrinsic elementary rates: no washcoat factor.
- ODE RHS/Jacobian: no washcoat factor.
- Exported Fluent-facing rates: `r_export = Wash_F * r_effective` when `CHEM_APPLY_WASHCOAT_ON_EXPORT_ONLY=1`.

## 5) Intentional behavior changes
- `legacy surface-only washcoat factor` removed from new migrated architecture.
- Rate evaluation moved from repeated `DEFINE_SR_RATE` callbacks to centralized `DEFINE_ADJUST` chemistry update.
- Added explicit cache-based interface between chemistry integration and Fluent rate-return path.
- Added explicit site-balance projection/clipping policy.
- Added named compile-time behavior flags for migration transparency.

## Notes
- Current code is structured for full mechanism porting from legacy branch formulas into `chem_compute_intrinsic_rates()` and stoichiometric RHS assembly in `chem_build_rhs()`.
- Solver API is finalized for future semi-implicit extrapolation upgrade; current backend is a robust minimal integrator scaffold.
