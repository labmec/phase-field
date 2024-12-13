#ifndef TPZPHASEFIELDANALYSIS_H
#define TPZPHASEFIELDANALYSIS_H
#include "TPZLinearAnalysis.h"

class TPZPhaseFieldAnalysis : public TPZLinearAnalysis {
 public:
  /** @name Constructors */
  /** @{ */
  /** @brief Create an empty TPZPhaseFieldAnalysis object */
  TPZPhaseFieldAnalysis();

  /** @brief Create an TPZPhaseFieldAnalysis object from one mesh pointer */
  TPZPhaseFieldAnalysis(TPZCompMesh *mesh, const RenumType &renumtype = RenumType::EDefault, std::ostream &out = std::cout);
  /** @brief Create an TPZPhaseFieldAnalysis object from one mesh auto pointer object */
  TPZPhaseFieldAnalysis(TPZAutoPointer<TPZCompMesh> mesh, const RenumType &renumtype = RenumType::EDefault, std::ostream &out = std::cout);
  /** @} */



};

#endif