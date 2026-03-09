#include "udf.h"
#include <string.h>
#include "chem_jacobian.h"
#include "chem_constants.h"

void chem_build_jacobian(const ChemInputs *in, const ChemState *state,
                         real jac[CHEM_N_SPECIES][CHEM_N_SPECIES])
{
  int i,j;
  real f0[CHEM_N_SPECIES], fp[CHEM_N_SPECIES];
  ChemState pert;

  chem_build_rhs(in, state, f0);
  memset(jac, 0, sizeof(real)*CHEM_N_SPECIES*CHEM_N_SPECIES);

#if CHEM_USE_NUMERICAL_JACOBIAN
  for (j=0;j<CHEM_N_SPECIES;j++) {
    const real eps = 1.0e-8;
    pert = *state;
    pert.theta[j] += eps;
    chem_build_rhs(in, &pert, fp);
    for (i=0;i<CHEM_N_SPECIES;i++) jac[i][j] = (fp[i]-f0[i])/eps;
  }
#endif
}
