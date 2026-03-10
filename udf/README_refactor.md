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
- `legacy61_bridge.h/.c`: migration bridge evaluating `reaction-1..reaction-61` through an internalized chemistry module and caching them.
- `chem_internalized_61rxn.c`: internalized copy of all 61 reaction equations/rate coefficients (no runtime dependency on repository legacy file).

## 3) Data flow
1. `DEFINE_INIT`: allocates initial chemistry state in UDM slots.
2. `DEFINE_ADJUST`: wall-face loop -> fetch local gas/wall state -> load cached chemistry -> integrate surface state ODE -> compute intrinsic rates -> apply eta/washcoat in export layer -> write cache.
3. `DEFINE_NET_REACTION_RATE` and/or `DEFINE_SOURCE` hooks consume cached chemistry outputs for returned rates/sources.

`DEFINE_SR_RATE` is kept out of active hooks in this migrated path (reference-only semantics). No integration is performed in return hooks.

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


- `CHEM_USE_LEGACY_61_RATE_BRIDGE` (default ON): imports all 61 reaction-rate formulas through the validated legacy dispatcher once per face update in `DEFINE_ADJUST`, then caches them for return hooks.

## Notes
- All 61 reactions are now reflected in the new flow through the legacy-bridge path (`CHEM_USE_LEGACY_61_RATE_BRIDGE=1`) during migration hardening. Native table-based intrinsic-rate porting remains available via `chem_compute_intrinsic_rates()` when the bridge is disabled.
- Solver API is finalized for future semi-implicit extrapolation upgrade; current backend is a robust minimal integrator scaffold.

## Homogeneous gas-phase option (DEFINE_NET_REACTION_RATE)
- Gas species order for net-rate UDF is explicitly fixed to 9 species: O2, C3H6, H2, H2O, CO2, CO, NO, NO2, N2 (must match Fluent material species order).
- Added `net_reaction_rate_udf.c` implementing `DEFINE_NET_REACTION_RATE(homogeneous_net_rates, ...)`.
- This path computes species net molar rates (`rr[]`) and Jacobian (`jac[]`) directly from `pressure[0]`, `temp[0]`, and `yi[]`.
- If your workflow returns source terms from net gas-phase reaction rates via `DEFINE_SOURCE`, a `DEFINE_SR_RATE` hook is not required (it may remain only as reference/compatibility).


## 61 reactions reflection policy
- In bridge mode (`CHEM_USE_LEGACY_61_RATE_BRIDGE=1`), all 61 reaction rates are evaluated by calling the internalized dispatcher (`chatterjee_pt_ads_des_flat_internalized`) for `reaction-1` to `reaction-61` each chemistry update.
- Therefore all equations and kinetic coefficients already internalized in `udf/chem_internalized_61rxn.c` are reflected directly in cached migrated rates.


## Hooking policy
- `DEFINE_SR_RATE` and `DEFINE_NET_REACTION_RATE` are not used together in active hook registration.
- In this refactor, `DEFINE_SR_RATE` is not provided as an active Fluent hook and must not be hooked.
- Runtime return path should use `DEFINE_NET_REACTION_RATE` and/or `DEFINE_SOURCE` only.
