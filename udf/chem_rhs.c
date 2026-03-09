#include "udf.h"
#include "chem_rhs.h"
#include "chem_indices.h"

void chem_build_rhs(const ChemInputs *in, const ChemState *state, real rhs[CHEM_N_SPECIES])
{
  int i;
  ChemState tmp = *state;
  chem_compute_intrinsic_rates(in, &tmp);
  for (i=0;i<CHEM_N_SPECIES;i++) rhs[i]=0.0;

  /* Minimal RHS scaffold; full 61-step stoichiometric assembly to be populated from legacy mechanism blocks. */
  rhs[IDX_Pt_Vac] -= tmp.rates_intrinsic[0]*2.0;
  rhs[IDX_O_Pt]   += tmp.rates_intrinsic[0]*2.0;
  rhs[IDX_Pt_Vac] -= tmp.rates_intrinsic[6];
  rhs[IDX_CO_Pt]  += tmp.rates_intrinsic[6];
}
