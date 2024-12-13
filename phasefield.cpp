#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <Elasticity/TPZElasticity2D.h>
#include <TPZGmshReader.h>
#include <TPZLinearAnalysis.h>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZNullMaterial.h>
#include <TPZSSpStructMatrix.h>
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>

#include <iostream>

#include "TPZAnalyticSolution.h"
#include "TPZCompElH1.h"
#include "TPZElementMatrixT.h"
#include "TPZGenGrid2D.h"
#include "TPZGeoMeshTools.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzlog.h"
#include "pzshapequad.h"
#include "pzvisualmatrix.h"
#include "tpzchangeel.h"
#include "TPZPhaseFieldAnalysis.h"
#include "TPZPhaseField.h"
#include "TPZElasticityPhaseField.h"

enum EMatid { ENone,
              EDomain,
              EDomainFrac,
              EDispX,
              EDispY,
              EDispXY,
              EPtDispX,
              EPtDispY,
              EPtDispXY,
              EFixedX,
              EFixedY,
              EFixedXY,
              EPtFixedX,
              EPtFixedY,
              EPtFixedXY,
              EForceX,
              EForceY,
              EForceXY,
              EPtForceX,
              EPtForceY,
              EPtForceXY };
const int global_nthread = 8;

REAL pseudotime = 0;
// const REAL maxdisp = 0.1;
const REAL maxdisp = 0.05;
auto applied_disp = [](const TPZVec<REAL> &coord, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal) {
  rhsVal[0] = 0.0;
  rhsVal[1] = - maxdisp * pseudotime * 1.e11; // Bignumber
};

TPZGeoMesh* CreateGMesh(int ndivx, int ndivy);
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
void SetPointBC(TPZGeoMesh* gr, TPZVec<REAL>& x, int bc);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity2DAnalytic* elas);
void SolveProblemDirect(TPZLinearAnalysis& an, TPZCompMesh* cmesh);
void PrintResults(TPZLinearAnalysis& an, TPZCompMesh* cmesh);

void GetSolVec(TPZInterpolationSpace* intel, TPZFMatrix<STATE>& u);
void SolveNormalH1Elatisticity(TPZGeoMesh* gmesh, const int pord);

TPZCompMesh* CreateElasticityAtomicMesh(TPZGeoMesh* gmesh, const int pord);
TPZCompMesh* CreatePhaseFieldAtomicMesh(TPZGeoMesh* gmesh, const int pord);

TPZMultiphysicsCompMesh* CreateElasticityMultiphysicsMesh(TPZManVector<TPZCompMesh*, 2>& mesh_vec, const int pord, const REAL E, const REAL nu);
TPZMultiphysicsCompMesh* CreatePhaseFieldMultiphysicsMesh(TPZManVector<TPZCompMesh*, 2>& mesh_vec, const int pord, const REAL Gc, const REAL l0);

void SolveStaggeredProblem(TPZLinearAnalysis& anElas, TPZLinearAnalysis& anPF);
void SolveIncrementalProblem(TPZLinearAnalysis& anElas, TPZLinearAnalysis& anPF, const int ntimesteps);

void SetPFAndElasMaterialPointer(TPZMultiphysicsCompMesh* mp_cmeshElas, TPZMultiphysicsCompMesh* mp_cmeshPF);
void TransferFromMultiBothMeshes(TPZMultiphysicsCompMesh* mp_cmeshElas, TPZMultiphysicsCompMesh* mp_cmeshPF);


std::string plotfile = "post_pfelas";
int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  const int whichprob = 1;

  int pord = 1;
  REAL E = 20.8e3, nu = 0.3;
  REAL Gc = 0.5, l0 = 0.03;  

  if(whichprob == 1){
    E = 20.e6;
    nu = 0.3;
    Gc = 1.0;
    l0 = 0.06;
  }

  bool isReadFromGmsh = true;
  TPZGeoMesh* gmesh = nullptr;
  if (isReadFromGmsh) {
    if (whichprob == 0)
      gmesh = ReadMeshFromGmsh("../gmsh_meshes/3-pt-bending-disp-surf-topmat.msh");
    else if (whichprob == 1)
      gmesh = ReadMeshFromGmsh("../gmsh_meshes/bittencourt.msh");
  } else {
    int ndivx = 25, ndivy = 50;
    gmesh = CreateGMesh(ndivx, ndivy);
  }

  std::ofstream out("gmesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

  const bool justSolveNormalElasticityAndQuit = false;
  if (justSolveNormalElasticityAndQuit) {
    plotfile = "postprocess";
    SolveNormalH1Elatisticity(gmesh, pord);
    return 0;
  }

  TPZCompMesh* cmeshElas = CreateElasticityAtomicMesh(gmesh, pord);
  TPZCompMesh* cmeshPF = CreatePhaseFieldAtomicMesh(gmesh, pord);
  

  TPZManVector<TPZCompMesh*, 2> mesh_vec(2);
  mesh_vec[0] = cmeshElas;
  mesh_vec[1] = cmeshPF;
  // TPZMultiphysicsCompMesh* mpmesh = CreateMultiphysicsMesh(mesh_vec, pord);
  TPZMultiphysicsCompMesh* mp_cmeshElas = CreateElasticityMultiphysicsMesh(mesh_vec, pord, E, nu);
  TPZMultiphysicsCompMesh* mp_cmeshPF = CreatePhaseFieldMultiphysicsMesh(mesh_vec, pord, Gc, nu);
  SetPFAndElasMaterialPointer(mp_cmeshElas, mp_cmeshPF);
  // matelas->SetPhaseFieldMaterial(matpf);
  // matpf->SetElasticityMaterial(matelas);


  // Solve problems in staggered manner
  TPZLinearAnalysis anElas(mp_cmeshElas);
  TPZSSpStructMatrix<STATE> matskl_elas(mp_cmeshElas);
  matskl_elas.SetNumThreads(global_nthread);
  anElas.SetStructuralMatrix(matskl_elas);

  TPZStepSolver<STATE> stepElas;
  stepElas.SetDirect(ECholesky);  // ELU //ECholesky // ELDLt
  anElas.SetSolver(stepElas);


  TPZLinearAnalysis anPF(mp_cmeshPF);
  TPZSSpStructMatrix<STATE> matskl_PF(mp_cmeshPF);
  matskl_PF.SetNumThreads(global_nthread);
  anPF.SetStructuralMatrix(matskl_PF);

  TPZStepSolver<STATE> stepPF;
  stepPF.SetDirect(ECholesky);  // ELU //ECholesky // ELDLt
  anPF.SetSolver(stepPF);

  const int ntimesteps = 100;
  SolveIncrementalProblem(anElas,anPF,ntimesteps);

  
  delete mp_cmeshElas;
  delete mp_cmeshPF;
  delete cmeshElas;
  delete cmeshPF;
  delete gmesh;

  std::cout << "--------- Simulation finished ---------" << std::endl;
}

void SolveIncrementalProblem(TPZLinearAnalysis& anElas, TPZLinearAnalysis& anPF, const int ntimesteps) {
  TPZElasticityPhaseField* elaspffrac = dynamic_cast<TPZElasticityPhaseField*>(anElas.Mesh()->FindMaterial(EDomainFrac));
  if (!elaspffrac) {
    DebugStop();
  }
  TPZElasticityPhaseField* elaspf = dynamic_cast<TPZElasticityPhaseField*>(anElas.Mesh()->FindMaterial(EDomain));
  if (!elaspf) {
    DebugStop();
  }

  // Initialize solution vector of anPF with 1 and anElas with 0
  TPZFMatrix<STATE> UPF = anPF.Solution();
  UPF = 1.;
  anPF.LoadSolution(UPF);
  // Update atomic meshes solution in multiphysics cmesh of anPF
  TPZMultiphysicsCompMesh* mp_cmeshPF = dynamic_cast<TPZMultiphysicsCompMesh*>(anPF.Mesh());
  if(!mp_cmeshPF) DebugStop();
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mp_cmeshPF->MeshVector(), mp_cmeshPF);
  
  TPZFMatrix<STATE> UElas = anElas.Solution();
  UElas.Zero();
  anElas.LoadSolution(UElas);
  // Update atomic meshes solution in multiphysics cmesh of anElas
  TPZMultiphysicsCompMesh* mp_cmeshElas = dynamic_cast<TPZMultiphysicsCompMesh*>(anElas.Mesh());
  if(!mp_cmeshElas) DebugStop();
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mp_cmeshElas->MeshVector(), mp_cmeshElas);

  for (int t = 0; t < ntimesteps; t++) {
    REAL time = (REAL)(t+1) / ntimesteps;
    std::cout << "******************** Time Step " << t << " | Pseudo time = " << std::fixed << std::setprecision(6) << time << " ********************" << std::endl;
    elaspffrac->SetTime(time);
    elaspf->SetTime(time);
    pseudotime = time;
    SolveStaggeredProblem(anElas, anPF);
    std::cout << "--------- PostProcess ---------" << std::endl;
    // TransferFromMultiBothMeshes(mp_cmeshElas, mp_cmeshPF);
    PrintResults(anElas, mp_cmeshElas);
  }
}

void SolveStaggeredProblem(TPZLinearAnalysis& anElas, TPZLinearAnalysis& anPF) {
  const int max_iterations = 50;
  const REAL tol = 5e-5, tolNormUPF = 1.e-10;
  int iteration = 0;
  REAL prevNormUPF = std::numeric_limits<REAL>::max(), currentNormUPF = std::numeric_limits<REAL>::max();
  REAL resElas = std::numeric_limits<REAL>::max();

  TPZMultiphysicsCompMesh* mp_cmeshPF = dynamic_cast<TPZMultiphysicsCompMesh*>(anPF.Mesh());
  if(!mp_cmeshPF) DebugStop();
  TPZMultiphysicsCompMesh* mp_cmeshElas = dynamic_cast<TPZMultiphysicsCompMesh*>(anElas.Mesh());
  if(!mp_cmeshElas) DebugStop();

  TPZSimpleTimer time_stag("Staggered Iteration");

  for (; iteration < max_iterations ; iteration++) {
    std::cout << "------ Staggered Iteration " << iteration << " ------" << std::endl;

    anElas.Assemble();

    // Right now, I will consider as converged if, after computing the updated phasefield, the norm of the residual in the elasticity problem is less than tol
    // Computer K*U - F
    auto matKU = anElas.MatrixSolver<STATE>().Matrix();
    TPZFMatrix<STATE> resElasVec = anElas.Rhs();
    TPZFMatrix<STATE> &UElas = anElas.Solution();
    TPZFMatrix<STATE> KU;
    matKU->Multiply(UElas, KU);
    TPZFMatrix<STATE> resElasComp = resElasVec - KU;
    resElas = Norm(resElasComp);
    std::cout << "Residual Elasticity Norm: " << std::scientific << std::setprecision(2) << resElas << std::endl;

    if(resElas < tol && iteration != 0) {
      std::cout << "Staggered scheme converged in " << iteration << " iterations." << std::endl;
      break;
    }

    anElas.Solve();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mp_cmeshElas->MeshVector(), mp_cmeshElas);
    
    anPF.Assemble();   
    anPF.Solve();
    currentNormUPF = Norm(anPF.Solution());
    const REAL varUPFNorm = fabs(currentNormUPF - prevNormUPF);
    if (iteration != 0)
      std::cout << "Variation on the norm of the phase field solution: " << std::scientific << std::setprecision(2) << varUPFNorm << std::endl;
    if(varUPFNorm < tolNormUPF) {
      std::cout << "Solution is not changing anymore. Stopping the staggered iterations and considering converged at iteration " << iteration << std::endl;
      break;
    }
    prevNormUPF = currentNormUPF;
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mp_cmeshElas->MeshVector(), mp_cmeshElas);
  }

  if (iteration == max_iterations) {
    std::cout << "WARNING! Maximum number of staggered iterations. Considering as converged and continuing" << std::endl;

  }

  std::cout << "Staggered iteration time: " << time_stag.ReturnTimeDouble() / 1000.0 << " seconds." << std::endl;

}

TPZMultiphysicsCompMesh* CreateElasticityMultiphysicsMesh(TPZManVector<TPZCompMesh*, 2>& mesh_vec, const int pord, const REAL E, const REAL nu) {
  TPZGeoMesh* gmesh = mesh_vec[0]->Reference();
  TPZMultiphysicsCompMesh* mp_cmesh = new TPZMultiphysicsCompMesh(gmesh);
  mp_cmesh->SetDefaultOrder(pord);
  mp_cmesh->SetDimModel(gmesh->Dimension());
  mp_cmesh->SetAllCreateFunctionsMultiphysicElem();

  // Create the TPZElasticityPhaseField material
  TPZElasticityPhaseField* mat = new TPZElasticityPhaseField(EDomain, gmesh->Dimension());
  mat->SetElasticity(E, nu);
  mat->SetPlaneStrain();
  mp_cmesh->InsertMaterialObject(mat);

  // Create the TPZElasticityPhaseField material for the fracture area
  TPZElasticityPhaseField* matfracArea = new TPZElasticityPhaseField(EDomainFrac, gmesh->Dimension());
  matfracArea->SetElasticity(E, nu);
  matfracArea->SetPlaneStrain();
  mp_cmesh->InsertMaterialObject(matfracArea);

  // Create the boundary conditions analogous to the H1 Cmesh
  const int diri = 0, neu = 1, mixed = 2, normaltrac = 4, diriincremental = -1;
  TPZFMatrix<STATE> val1(2, 2, 0.);
  TPZManVector<STATE> val2(2, 0.);
  auto* BCCondFixed1 = mat->CreateBC(mat, EPtFixedXY, diri, val1, val2);
  mp_cmesh->InsertMaterialObject(BCCondFixed1);

  val1(1, 1) = mat->BigNumber();
  auto* BCCondSym = mat->CreateBC(mat, EPtFixedY, mixed, val1, val2);
  mp_cmesh->InsertMaterialObject(BCCondSym);
  
  val1.Zero();
  val2[1] = -0.05;
  // auto* BCCondPoint = mat->CreateBC(mat, EPtDispY, diriincremental, val1, val2);
  // mp_cmesh->InsertMaterialObject(BCCondPoint);
  val1(1,1) = mat->BigNumber();
  auto* BCCondSurf = mat->CreateBC(mat, EDispY, mixed, val1, val2);
  BCCondSurf->SetForcingFunctionBC(applied_disp,2);
  mp_cmesh->InsertMaterialObject(BCCondSurf);

  val1(1,1) = mat->BigNumber();
  auto* BCCondPoint = mat->CreateBC(mat, EPtDispY, mixed, val1, val2);
  BCCondPoint->SetForcingFunctionBC(applied_disp,2);
  mp_cmesh->InsertMaterialObject(BCCondPoint);


  mp_cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
  TPZManVector<int, 2> active(2, 1);
  active[1] = 0;
  mp_cmesh->BuildMultiphysicsSpace(active, mesh_vec);

  return mp_cmesh;
}

TPZMultiphysicsCompMesh* CreatePhaseFieldMultiphysicsMesh(TPZManVector<TPZCompMesh*, 2>& mesh_vec, const int pord, const REAL Gc, const REAL l0) {
  TPZGeoMesh* gmesh = mesh_vec[0]->Reference();
  TPZMultiphysicsCompMesh* mp_cmesh = new TPZMultiphysicsCompMesh(gmesh);
  mp_cmesh->SetDefaultOrder(pord);
  mp_cmesh->SetDimModel(gmesh->Dimension());
  mp_cmesh->SetAllCreateFunctionsMultiphysicElem();

  // Create the TPZPhaseField material for the region without fracture
  const REAL c0 = -1.0; // not being used yet
  const REAL GcBig = Gc * 1000.;
  TPZPhaseField* mat = new TPZPhaseField(EDomain, gmesh->Dimension(), GcBig, l0, c0);
  mp_cmesh->InsertMaterialObject(mat);

  // Create the TPZPhaseField material for the fracture area
  TPZPhaseField* matfracArea = new TPZPhaseField(EDomainFrac, gmesh->Dimension(), Gc, l0, c0);
  mp_cmesh->InsertMaterialObject(matfracArea);

  // No BCs in Phase field because Neumann zero on all boundary

  // Build the multiphysic mesh
  mp_cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
  TPZManVector<int, 2> active(2, 1);
  active[0] = 0;
  mp_cmesh->BuildMultiphysicsSpace(active, mesh_vec);

  return mp_cmesh;
}

TPZCompMesh* CreateElasticityAtomicMesh(TPZGeoMesh* gmesh, const int pord) {
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
  const int dim = gmesh->Dimension();
  cmesh->SetDimModel(dim);
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  auto mat = new TPZNullMaterial(EDomain);
  mat->SetNStateVariables(2);
  cmesh->InsertMaterialObject(mat);
  auto* matfracArea = new TPZNullMaterial(EDomainFrac);
  matfracArea->SetNStateVariables(2);
  cmesh->InsertMaterialObject(matfracArea);

  TPZFMatrix<STATE> val1(2, 2, 0.);
  TPZManVector<STATE> val2(2, 0.);

  const int diri = 0, neu = 1, mixed = 2, normaltrac = 4, diriincremental = -1;;

  // BCs
  auto* BCCondFixed1 = mat->CreateBC(mat, EPtFixedXY, diri, val1, val2);
  cmesh->InsertMaterialObject(BCCondFixed1);

  val1(1, 1) = mat->BigNumber();
  auto* BCCondSym = mat->CreateBC(mat, EPtFixedY, mixed, val1, val2);
  cmesh->InsertMaterialObject(BCCondSym);

  val1.Zero();
  val2[1] = -0.05;
  auto* BCCondPoint = mat->CreateBC(mat, EPtDispY, diri, val1, val2);
  cmesh->InsertMaterialObject(BCCondPoint);
  auto* BCCondSurf = mat->CreateBC(mat, EDispY, mixed, val1, val2);
  cmesh->InsertMaterialObject(BCCondSurf);
  
  cmesh->AutoBuild();
  return cmesh;
}

TPZCompMesh* CreatePhaseFieldAtomicMesh(TPZGeoMesh* gmesh, const int pord) {
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
  const int dim = gmesh->Dimension();
  cmesh->SetDimModel(dim);
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  auto mat = new TPZNullMaterial(EDomain);
  cmesh->InsertMaterialObject(mat);
  auto* matfracArea = new TPZNullMaterial(EDomainFrac);
  cmesh->InsertMaterialObject(matfracArea);

  // No BCs in Phase field because Neumann zero on all boundary

  cmesh->AutoBuild();
  return cmesh;
}

void SolveNormalH1Elatisticity(TPZGeoMesh* gmesh, const int pord) {
  TElasticity2DAnalytic* elas = new TElasticity2DAnalytic;
  elas->gE = 20.8e3;
  elas->gPoisson = 0.3;
  elas->fProblemType = TElasticity2DAnalytic::EStretchx;
  const REAL Gc = 0.5;
  const REAL l0 = 0.03;
  const REAL GcBig = Gc * 1000.;
  TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh, pord, elas);

  TPZLinearAnalysis an(cmeshH1);

  SolveProblemDirect(an, cmeshH1);

  std::cout << "--------- PostProcess ---------" << std::endl;
  PrintResults(an, cmeshH1);

  delete cmeshH1;
}

TPZGeoMesh* CreateGMesh(int ndivx, int ndivy) {
  TPZGeoMesh* gmesh = new TPZGeoMesh;

  MMeshType meshType = MMeshType::EQuadrilateral;
  int dim = 2;
  TPZManVector<REAL, 3> minX = {0, 0, 0};
  TPZManVector<REAL, 3> maxX = {100, 200, 0};
  int nMats = 2 * dim + 1;

  constexpr bool createBoundEls{true};
  TPZVec<int> matIds(nMats, ENone);
  matIds[0] = EDomain;

  TPZManVector<int, 2> ndivvec = {ndivx, ndivy};
  gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX, matIds, ndivvec, meshType, createBoundEls);

  TPZManVector<REAL, 2> xfixed1 = {0., 0., 0.}, xfixed2 = {0., 200., 0.}, xforce = {100., 100., 0.};
  SetPointBC(gmesh, xfixed1, EPtDispXY);
  SetPointBC(gmesh, xfixed2, EPtDispXY);
  SetPointBC(gmesh, xforce, EPtForceY);

  return gmesh;
}

TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity2DAnalytic* elas) {
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
  const int dim = gmesh->Dimension();
  cmesh->SetDimModel(dim);
  cmesh->SetDefaultOrder(pord);
  cmesh->SetAllCreateFunctionsContinuous();

  const STATE E = elas->gE, nu = elas->gPoisson;
  TPZManVector<STATE> force = {0, 0, 0};
  TPZElasticity2D* mat = new TPZElasticity2D(EDomain, E, nu, 0., 0., true);
  //    mat->SetExactSol(elas->ExactSolution(), 2);
  //    mat->SetForcingFunction(elas->ForceFunc(), 4);
  cmesh->InsertMaterialObject(mat);
  TPZElasticity2D* matfrac = new TPZElasticity2D(EDomainFrac, E, nu, 0., 0., true);
  cmesh->InsertMaterialObject(matfrac);

  TPZFMatrix<STATE> val1(2, 2, 0.);
  TPZManVector<STATE> val2(2, 0.);

  const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;

  // BCs
  auto* BCCondFixed1 = mat->CreateBC(mat, EPtFixedXY, diri, val1, val2);
  cmesh->InsertMaterialObject(BCCondFixed1);

  val1(1, 1) = mat->BigNumber();
  auto* BCCondSym = mat->CreateBC(mat, EPtFixedY, mixed, val1, val2);
  cmesh->InsertMaterialObject(BCCondSym);

  val1.Zero();
  val2[1] = -0.05;
  auto* BCCondPoint = mat->CreateBC(mat, EPtDispY, diri, val1, val2);
  cmesh->InsertMaterialObject(BCCondPoint);
  auto* BCCondSurf = mat->CreateBC(mat, EDispY, diri, val1, val2);
  cmesh->InsertMaterialObject(BCCondSurf);


  cmesh->AutoBuild();
  return cmesh;
}

void SolveProblemDirect(TPZLinearAnalysis& an, TPZCompMesh* cmesh) {
  //    TPZSkylineStructMatrix<STATE> matskl(cmesh);
  TPZSSpStructMatrix<STATE> matskl(cmesh);
  matskl.SetNumThreads(global_nthread);
  an.SetStructuralMatrix(matskl);

  TPZStepSolver<STATE> step;
  step.SetDirect(ECholesky);  // ELU //ECholesky // ELDLt
  an.SetSolver(step);

  std::cout << "--------- Assemble ---------" << std::endl;
  TPZSimpleTimer time_ass;
  std::cout << "NElements = " << an.Mesh()->NElements() << std::endl;
  std::cout << "NEquations = " << an.Mesh()->NEquations() << std::endl;
  an.Assemble();
  std::cout << "Total time = " << time_ass.ReturnTimeDouble() / 1000. << " s" << std::endl;

  std::cout << "--------- Solve ---------" << std::endl;
  TPZSimpleTimer time_sol;
  an.Solve();
  std::cout << "Total time = " << time_sol.ReturnTimeDouble() / 1000. << " s" << std::endl;

  return;
}

void PrintResults(TPZLinearAnalysis& an, TPZCompMesh* cmesh) {
  std::cout << "--------- Post Process ---------" << std::endl;
  TPZSimpleTimer postProc("Post processing time");  
  constexpr int vtkRes{0};
  // TPZVec<std::string> fields = {
  //     "Displacement",
  //     "Stress"};

  TPZVec<std::string> fields = {
      "Displacement",
      "Stress",
      "PhaseField"};
  static auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.SetNThreads(0);
  vtk.Do();
  std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;

  return;
}

void SetPointBC(TPZGeoMesh* gr, TPZVec<REAL>& x, int bc) {
  // look for an element/corner node whose distance is close to start
  TPZGeoNode* gn1 = gr->FindNode(x);
  int64_t iel;
  int64_t nelem = gr->ElementVec().NElements();
  TPZGeoEl* gel;
  for (iel = 0; iel < nelem; iel++) {
    gel = gr->ElementVec()[iel];
    if (!gel) continue;
    int nc = gel->NCornerNodes();
    int c;
    for (c = 0; c < nc; c++) {
      TPZGeoNode* gn = gel->NodePtr(c);
      if (gn == gn1) {
        break;
      }
    }
    if (c < nc) {
      TPZGeoElBC(gel, c, bc);
      return;
    }
  }
}

void GetSolVec(TPZInterpolationSpace* intel, TPZFMatrix<STATE>& u) {
  const int nstate = intel->Material()->NStateVariables();
  const int ncon = intel->NConnects();
  TPZBlock& block = intel->Mesh()->Block();
  TPZFMatrix<STATE>& MeshSol = intel->Mesh()->Solution();
  const int64_t numbersol = MeshSol.Cols();
  if (numbersol != 1) DebugStop();  // I did not think about this case yet, but it can be done

  int64_t iv = 0;
  for (int in = 0; in < ncon; in++) {
    TPZConnect* df = &intel->Connect(in);
    const int64_t dfseq = df->SequenceNumber();
    const int dfvar = block.Size(dfseq);
    const int64_t pos = block.Position(dfseq);
    for (int jn = 0; jn < dfvar; jn++) {
      u(iv++, 0) = MeshSol(pos + jn, 0);
    }
  }
}

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name) {
  TPZGeoMesh* gmesh;
  gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(20);
    stringtoint[2]["dom"] = EDomain;
    stringtoint[2]["domfrac"] = EDomainFrac;

    stringtoint[1]["dispx"] = EDispX;
    stringtoint[1]["dispy"] = EDispY;
    stringtoint[1]["dispxy"] = EDispXY;
    stringtoint[0]["ptdispx"] = EPtDispX;
    stringtoint[0]["ptdispy"] = EPtDispY;
    stringtoint[0]["ptdispxy"] = EPtDispXY;

    stringtoint[1]["fixedx"] = EFixedX;
    stringtoint[1]["fixedy"] = EFixedY;
    stringtoint[1]["fixedxy"] = EFixedXY;
    stringtoint[0]["ptfixedx"] = EPtFixedX;
    stringtoint[0]["ptfixedy"] = EPtFixedY;
    stringtoint[0]["ptfixedxy"] = EPtFixedXY;

    stringtoint[1]["forcex"] = EForceX;
    stringtoint[1]["forcey"] = EForceY;
    stringtoint[1]["forcexy"] = EForceXY;
    stringtoint[0]["ptforcex"] = EPtForceX;
    stringtoint[0]["ptforcey"] = EPtForceY;
    stringtoint[0]["ptforcexy"] = EPtForceXY;
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}

void SetPFAndElasMaterialPointer(TPZMultiphysicsCompMesh* mp_cmeshElas, TPZMultiphysicsCompMesh* mp_cmeshPF) {
  TPZElasticityPhaseField* matelasPF = dynamic_cast<TPZElasticityPhaseField*>(mp_cmeshElas->FindMaterial(EDomain));
  TPZPhaseField* matpfPF = dynamic_cast<TPZPhaseField*>(mp_cmeshPF->FindMaterial(EDomain));  
  if(!matelasPF || !matpfPF) {
    DebugStop();
  }
  matelasPF->SetPhaseFieldMaterial(matpfPF);
  matpfPF->SetElasticityMaterial(matelasPF);

  TPZElasticityPhaseField* matelasPFFrac = dynamic_cast<TPZElasticityPhaseField*>(mp_cmeshElas->FindMaterial(EDomainFrac));
  TPZPhaseField* matpfPFFrac = dynamic_cast<TPZPhaseField*>(mp_cmeshPF->FindMaterial(EDomainFrac));  
  if(!matelasPFFrac || !matpfPFFrac) {
    DebugStop();
  }
  matelasPFFrac->SetPhaseFieldMaterial(matpfPFFrac);
  matpfPFFrac->SetElasticityMaterial(matelasPFFrac);
}

void TransferFromMultiBothMeshes(TPZMultiphysicsCompMesh* mp_cmeshElas, TPZMultiphysicsCompMesh* mp_cmeshPF) {

}