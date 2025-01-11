#include "TPZPhaseField.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"

const REAL TPZPhaseField::eta = 1.0e-8;

const bool oldWay = false;

TPZPhaseField::TPZPhaseField() : TPZRegisterClassId(&TPZPhaseField::ClassId),
                                 TBase(),
                                 fDim(-1) {}

TPZPhaseField::TPZPhaseField(int id, int dim, STATE Gc, STATE l0, STATE c0) : TPZRegisterClassId(&TPZPhaseField::ClassId),
                                                                              TBase(id),
                                                                              fDim(dim),
                                                                              fGc(Gc),
                                                                              fl0(l0),
                                                                              fc0(c0) {
}

/**
         copy constructor
 */
TPZPhaseField::TPZPhaseField(const TPZPhaseField &copy) : TBase(copy), fDim(copy.fDim) {
}
/**
         copy constructor
 */
TPZPhaseField &TPZPhaseField::operator=(const TPZPhaseField &copy) {
  TBase::operator=(copy);
  fDim = copy.fDim;
  return *this;
}

void TPZPhaseField::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef) {
  const STATE Gc = fGc;
  const STATE l0 = fl0;
  const STATE c0 = fc0;

  // Setting the phis
  TPZFMatrix<REAL> &phi = datavec[0].phi;
  TPZFMatrix<REAL> &dphiaxes = datavec[0].dphix;
  TPZFNMatrix<9, REAL> dphi(3, dphiaxes.Cols());
  TPZAxesTools<REAL>::Axes2XYZ(dphiaxes, dphi, datavec[0].axes);

  int phr = phi.Rows();

  int nactive = 0;
  for (const auto &i : datavec) {
    if (i.fActiveApproxSpace) {
      nactive++;
    }
  }
  const STATE sigmaDotEps = fElasMat->CalculateSigmaDotEps(datavec[0]);
  const STATE elasEnergy = sigmaDotEps * 0.5;
  
  if(oldWay){
    for (int i = 0; i < phr; i++) {
      ef(i,0) += weight * ( (Gc/ l0 ) * phi(i) );
      for (int j = 0; j < phr; j++) {
        ek(i, j) += weight * (Gc * l0 * (dphi(0,i) * dphi(0,j) + dphi(1,i) * dphi(1,j) ) + (sigmaDotEps  + Gc/l0) * phi(i) * phi(j) );
      }
    }
  } 
  else {
    for (int i = 0; i < phr; i++) {
      ef(i,0) += (1.) * weight * elasEnergy * phi(i);
      for (int j = 0; j < phr; j++) {
        ek(i, j) += weight * (Gc * l0 / c0 * (dphi(0,i) * dphi(0,j) + dphi(1,i) * dphi(1,j) ) + (Gc/(c0*l0) + elasEnergy ) * phi(j) * phi(i) );
      }
    }        
  }
}

void TPZPhaseField::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                 TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
  DebugStop();  // Should not be called unless imposing the phase-field as a dirichlet bc (not that common).
}

void TPZPhaseField::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) {
  solOut.Resize(this->NSolutionVariables(var));
  solOut.Fill(0.);
  TPZManVector<STATE, 10> SolPF;

  // SolQ = datavec[0].sol[0];
  SolPF = datavec[1].sol[0];
  if (SolPF.size() == 0) SolPF.Resize(1, 0.);

  if (var == 1) {  // phase field
    solOut[0] = datavec[0].sol[0][0];
    return;
  }
}

void TPZPhaseField::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
  // Not implemented
  DebugStop();
}

int TPZPhaseField::VariableIndex(const std::string &name) const {
  if (!strcmp("Phasefield", name.c_str())) return 1;
  if (!strcmp("PhaseField", name.c_str())) return 1;
  if (!strcmp("PF", name.c_str())) return 1;
  DebugStop();
  return -1;
}

int TPZPhaseField::NSolutionVariables(int var) const {
  if (var == 1) return 1;
  DebugStop();
  return -1;
}

void TPZPhaseField::SetDimension(int dim) {
  if (dim > 3 || dim < 1) DebugStop();
  if (dim == 3) DebugStop();  // Never tested with 3d... if decide to do it, please proceed with caution
  fDim = dim;
}

int TPZPhaseField::ClassId() const {
  return Hash("TPZPhaseField") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZPhaseField::NewMaterial() const {
  return new TPZPhaseField(*this);
}

void TPZPhaseField::Print(std::ostream &out) const {
  out << "Material Name: " << this->Name() << "\n";
  out << "Material Id: " << this->Id() << "\n";
  out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZPhaseField::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
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

void TPZPhaseField::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
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

const STATE TPZPhaseField::DegradationFunc(const STATE pf) const {
  int type = 1;
  if(oldWay) type = 0;
  if(type == 0){
    return pf*pf + eta;
  }
  else{
    return (1. - pf)*(1. - pf) + eta;
  }
}