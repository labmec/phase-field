#include "TPZElasticityPhaseField.h"
#include "pzaxestools.h"



TPZElasticityPhaseField::TPZElasticityPhaseField() : TPZRegisterClassId(&TPZElasticityPhaseField::ClassId),
                                                     TPZMatCombinedSpacesT<STATE>(),
                                                     TPZElasticity2D() {}

TPZElasticityPhaseField::TPZElasticityPhaseField(int id, int dim) : TPZRegisterClassId(&TPZElasticityPhaseField::ClassId),
                                                                    TPZElasticity2D(id) {}

int TPZElasticityPhaseField::VariableIndex(const std::string &name) const {
  if(!strcmp("PhaseField",name.c_str()))     return 50;
  if(!strcmp("Phasefield",name.c_str()))     return 50;
  if(!strcmp("PF",name.c_str()))     return 50;
  return TPZElasticity2D::VariableIndex(name);
}

int TPZElasticityPhaseField::NSolutionVariables(int var) const {
  switch(var) {
		case 50:
			return 1;
  } 
  return TPZElasticity2D::NSolutionVariables(var);
}

int TPZElasticityPhaseField::ClassId() const {
  return Hash("TPZElasticityPhaseField") ^ TPZElasticity2D::ClassId() << 1;
}

TPZMaterial *TPZElasticityPhaseField::NewMaterial() const {
  return new TPZElasticityPhaseField(*this);
}

void TPZElasticityPhaseField::Print(std::ostream &out) const {
  out << "Material Name: " << this->Name() << "\n";
  out << "Material Id: " << TPZElasticity2D::Id() << "\n";
  out << "Dimension: " << TPZElasticity2D::Dimension() << "\n\n";
}

void TPZElasticityPhaseField::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                         REAL weight, TPZFMatrix<STATE> &ek,
                                         TPZFMatrix<STATE> &ef) {
  
  // Call the contribute method of father class TPZElasticity2D
  TPZFMatrix<STATE> temp_ek(ek.Rows(), ek.Cols());
  temp_ek.Zero();
  TPZElasticity2D::Contribute(datavec[0], weight, temp_ek, ef);
  
  // Get the value of the phase field material
  const STATE pf = datavec[1].sol[0][0];

  // Multiply temp_ek by cte and add to ek
  const REAL g = fPhaseFieldMat->DegradationFunc(pf);
  temp_ek *= g;
  ek += temp_ek;

}

void TPZElasticityPhaseField::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

  /// Multiply val2 by fTime  
  if (bc.Type() == -1) { // Dirichlet incremental    
    TPZVec<STATE> bcVal2cp = bc.Val2(), bcMem = bc.Val2();
    for (auto &val : bcVal2cp) {
      val *= fTime;
    }
    bc.SetType(0);
    bc.SetVal2(bcVal2cp);
    TPZElasticity2D::ContributeBC(datavec[0], weight, ek, ef, bc);
    bc.SetVal2(bcMem);
    bc.SetType(-1);
  }
  else {
    // Call the contributeBC of father class TPZElasticity2D
    TPZElasticity2D::ContributeBC(datavec[0], weight, ek, ef, bc);    
  }
}

void TPZElasticityPhaseField::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                                       int var, TPZVec<STATE> &sol) {
  if (var == 50) {
    sol[0] = datavec[1].sol[0][0];
    return;
  }
  TPZElasticity2D::Solution(datavec[0], var, sol);
}

void TPZElasticityPhaseField::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
  DebugStop(); // not implemented
}

void TPZElasticityPhaseField::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
  int nref = datavec.size();
  for (int i = 0; i < nref; i++) {
    datavec[i].SetAllRequirements(false);
    datavec[i].fNeedsNeighborSol = false;
    datavec[i].fNeedsNeighborCenter = false;
    datavec[i].fNeedsNormal = false;
    datavec[i].fNeedsHSize = false;
    datavec[i].fNeedsSol = true;
  }
}

void TPZElasticityPhaseField::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
  // default is no specific data requirements
  int nref = datavec.size();
  for (int iref = 0; iref < nref; iref++) {
    datavec[iref].SetAllRequirements(false);
    datavec[iref].fNeedsSol = true;
  }
  datavec[0].fNeedsNormal = true;
  if (type == 50) {
    DebugStop();  // What is this?
    for (int iref = 0; iref < nref; iref++) {
      datavec[iref].fNeedsSol = false;
    }
  }
}