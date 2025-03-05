#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for SMConstants.C_vev0, SMConstants.C_MassTop, SMConstants.C_g
#include <algorithm> // for max, copy
#include <iomanip>
#include <iostream> // for operator<<, endl, basic_o...
#include <memory>   // for allocator_traits<>::value...
#include <stddef.h> // for std::size_t

#include <BSMPT/models/ClassPotentialRxSM_OS_lagparams.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_RxSM_OS_lagparams::Class_Potential_RxSM_OS_lagparams(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::RXSM_OS_LAGPARAMS;

  nPar = 5;   // number of parameters in the tree-Level Lagrangian AFTER using
               // tadpole equations
  nParCT = 10; // number of parameters in the counterterm potential

  nVEV = 2; // number of VEVs to minimize the potential

  NHiggs = 5; // number of scalar d.o.f.

  NGauge = 4; // number of gauge fields

  NLepton = 9; // number of lepton fields

  NQuarks = 12; // number of quark fields

  VevOrder.resize(nVEV);
  VevOrder[0] = 2; // wh
  VevOrder[1] = 4; // ws

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_RxSM_OS_lagparams::~Class_Potential_RxSM_OS_lagparams()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_RxSM_OS_lagparams::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmh");
  labels.push_back("dms");
  labels.push_back("dLh");
  labels.push_back("dLs");
  labels.push_back("dLhs");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Potential_RxSM_OS_lagparams::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");     // Label for the critical temperature
  labels.push_back("v_c");     // Label for the critical vev
  labels.push_back("v_c/T_c"); // Label for xi_c
  // out += "VEV order";
  labels.push_back("wh(T_c)");
  labels.push_back("ws(T_c)");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 */
std::vector<std::string>
Class_Potential_RxSM_OS_lagparams::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  // mass basis, you can identify here your particles
  particles.push_back("h_1");
  particles.push_back("h_2");
  particles.push_back("h_3");
  particles.push_back("h_4");
  particles.push_back("h_5");

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back("Tree_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CT_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CW_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_RxSM_OS_lagparams::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order"
  labels.push_back("wh");
  labels.push_back("ws");

  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_RxSM_OS_lagparams::ReadAndSet(const std::string &linestr,
                                             std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 5; k++)
  {
    ss >> tmp;
    if (k == 1)
      par[0] = tmp; // Lh
    if (k == 2)
      par[1] = tmp; // Ls
    if (k == 3)
      par[2] = tmp; // Lhs
    if (k == 4)
      par[3] = tmp; // vh
    if (k == 5)
      par[4] = tmp; // vs

  }

  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_RxSM_OS_lagparams::set_gen(const std::vector<double> &par)
{

  Lh = par[0]; 
  Ls = par[1]; 
  Lhs = par[2]; 
  vh = par[3]; 
  vs = par[4]; 

  mh = (Lh*pow(vh,2))/6. + Lhs*pow(vs,2); 
  ms = Lhs*pow(vh,2) + (Ls*pow(vs,2))/6.; 

  scale = SMConstants.C_vev0; // renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // set the vector vevTreeMin. vevTree will then be set by the
  // function MinimizeOrderVEV
  vevTreeMin[0] = vh; // wh
  vevTreeMin[1] = vs; // ws

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_Potential_RxSM_OS_lagparams::set_CT_Pot_Par(const std::vector<double> &par)
{
  dmh = par[0];
  dms = par[1];
  dLh = par[2];
  dLs = par[3];
  dLhs = par[4];
  dT1 = par[5];
  dT2 = par[6];
  dT3 = par[7];
  dT4 = par[8];
  dT5 = par[9];

  // assign the non-zero entries
  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;

  Curvature_Higgs_CT_L2[0][0] = -dmh;
  Curvature_Higgs_CT_L2[1][1] = -dmh;
  Curvature_Higgs_CT_L2[2][2] = -dmh;
  Curvature_Higgs_CT_L2[3][3] = -dmh;
  Curvature_Higgs_CT_L2[4][4] = -dms;

  Curvature_Higgs_CT_L4[0][0][0][0] = dLh;
  Curvature_Higgs_CT_L4[0][0][1][1] = dLh/3.;
  Curvature_Higgs_CT_L4[0][0][2][2] = dLh/3.;
  Curvature_Higgs_CT_L4[0][0][3][3] = dLh/3.;
  Curvature_Higgs_CT_L4[0][0][4][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[0][1][0][1] = dLh/3.;
  Curvature_Higgs_CT_L4[0][1][1][0] = dLh/3.;
  Curvature_Higgs_CT_L4[0][2][0][2] = dLh/3.;
  Curvature_Higgs_CT_L4[0][2][2][0] = dLh/3.;
  Curvature_Higgs_CT_L4[0][3][0][3] = dLh/3.;
  Curvature_Higgs_CT_L4[0][3][3][0] = dLh/3.;
  Curvature_Higgs_CT_L4[0][4][0][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[0][4][4][0] = 2*dLhs;
  Curvature_Higgs_CT_L4[1][0][0][1] = dLh/3.;
  Curvature_Higgs_CT_L4[1][0][1][0] = dLh/3.;
  Curvature_Higgs_CT_L4[1][1][0][0] = dLh/3.;
  Curvature_Higgs_CT_L4[1][1][1][1] = dLh;
  Curvature_Higgs_CT_L4[1][1][2][2] = dLh/3.;
  Curvature_Higgs_CT_L4[1][1][3][3] = dLh/3.;
  Curvature_Higgs_CT_L4[1][1][4][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[1][2][1][2] = dLh/3.;
  Curvature_Higgs_CT_L4[1][2][2][1] = dLh/3.;
  Curvature_Higgs_CT_L4[1][3][1][3] = dLh/3.;
  Curvature_Higgs_CT_L4[1][3][3][1] = dLh/3.;
  Curvature_Higgs_CT_L4[1][4][1][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[1][4][4][1] = 2*dLhs;
  Curvature_Higgs_CT_L4[2][0][0][2] = dLh/3.;
  Curvature_Higgs_CT_L4[2][0][2][0] = dLh/3.;
  Curvature_Higgs_CT_L4[2][1][1][2] = dLh/3.;
  Curvature_Higgs_CT_L4[2][1][2][1] = dLh/3.;
  Curvature_Higgs_CT_L4[2][2][0][0] = dLh/3.;
  Curvature_Higgs_CT_L4[2][2][1][1] = dLh/3.;
  Curvature_Higgs_CT_L4[2][2][2][2] = dLh;
  Curvature_Higgs_CT_L4[2][2][3][3] = dLh/3.;
  Curvature_Higgs_CT_L4[2][2][4][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[2][3][2][3] = dLh/3.;
  Curvature_Higgs_CT_L4[2][3][3][2] = dLh/3.;
  Curvature_Higgs_CT_L4[2][4][2][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[2][4][4][2] = 2*dLhs;
  Curvature_Higgs_CT_L4[3][0][0][3] = dLh/3.;
  Curvature_Higgs_CT_L4[3][0][3][0] = dLh/3.;
  Curvature_Higgs_CT_L4[3][1][1][3] = dLh/3.;
  Curvature_Higgs_CT_L4[3][1][3][1] = dLh/3.;
  Curvature_Higgs_CT_L4[3][2][2][3] = dLh/3.;
  Curvature_Higgs_CT_L4[3][2][3][2] = dLh/3.;
  Curvature_Higgs_CT_L4[3][3][0][0] = dLh/3.;
  Curvature_Higgs_CT_L4[3][3][1][1] = dLh/3.;
  Curvature_Higgs_CT_L4[3][3][2][2] = dLh/3.;
  Curvature_Higgs_CT_L4[3][3][3][3] = dLh;
  Curvature_Higgs_CT_L4[3][3][4][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[3][4][3][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[3][4][4][3] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][0][0][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][0][4][0] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][1][1][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][1][4][1] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][2][2][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][2][4][2] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][3][3][4] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][3][4][3] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][4][0][0] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][4][1][1] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][4][2][2] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][4][3][3] = 2*dLhs;
  Curvature_Higgs_CT_L4[4][4][4][4] = dLs;


}

/**
 * console output of all parameters
 */
void Class_Potential_RxSM_OS_lagparams::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "Model = " << Model << "\n";

  ss << "\nThe input parameters are : \n";
  ss << "Lh = " << Lh << "\n";
  ss << "Ls = " << Ls << "\n";
  ss << "Lhs = " << Lhs << "\n";
  ss << "vh = " << vh << "\n";
  ss << "vs = " << vs << "\n";

  ss << "\nThe parameters are : \n";
  ss << "mh = " << mh << "\n";
  ss << "ms = " << ms << "\n";
  ss << "Lh = " << Lh << "\n";
  ss << "Ls = " << Ls << "\n";
  ss << "Lhs = " << Lhs << "\n";

  ss << "\nThe counterterm parameters are : \n";
  ss << "dmh = " << dmh << "\n";
  ss << "dms = " << dms << "\n";
  ss << "dLh = " << dLh << "\n";
  ss << "dLs = " << dLs << "\n";
  ss << "dLhs = " << dLhs << "\n";
  ss << "dT1 = " << dT1 << "\n";
  ss << "dT2 = " << dT2 << "\n";
  ss << "dT3 = " << dT3 << "\n";
  ss << "dT4 = " << dT4 << "\n";
  ss << "dT5 = " << dT5 << "\n";

  ss << "\nThe scale is given by mu = " << scale << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Potential_RxSM_OS_lagparams::calc_CT() const
{
  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsdone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // formulae for the counterterm scheme
  parCT.push_back(-0.5*(vh*(-3*HesseWeinberg(0,0) + HesseWeinberg(2,2)) + vs*HesseWeinberg(2,4))/vh); //dmh;
  parCT.push_back(-0.5*(vh*HesseWeinberg(2,4) + vs*HesseWeinberg(4,4) - 3*NablaWeinberg(4))/vs); //dms;
  parCT.push_back((3*(HesseWeinberg(0,0) - HesseWeinberg(2,2)))/pow(vh,2)); //dLh;
  parCT.push_back((-3*vs*HesseWeinberg(4,4) + 3*NablaWeinberg(4))/pow(vs,3)); //dLs;
  parCT.push_back(-0.5*HesseWeinberg(2,4)/(vh*vs)); //dLhs;
  parCT.push_back(-NablaWeinberg(0)); //dT1;
  parCT.push_back(-NablaWeinberg(1)); //dT2;
  parCT.push_back(vh*HesseWeinberg(0,0) - NablaWeinberg(2)); //dT3;
  parCT.push_back(-NablaWeinberg(3)); //dT4;
  parCT.push_back(0); //dT5;

  return parCT;
}

// mass basis triple couplings
void Class_Potential_RxSM_OS_lagparams::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  // new rotation matrix with
  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  std::vector<double> HiggsOrder(NHiggs);

  // example for keeping the mass order
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void Class_Potential_RxSM_OS_lagparams::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // assign the non-zero entries
  Curvature_Higgs_L2[0][0] = -mh;
  Curvature_Higgs_L2[1][1] = -mh;
  Curvature_Higgs_L2[2][2] = -mh;
  Curvature_Higgs_L2[3][3] = -mh;
  Curvature_Higgs_L2[4][4] = -ms;

  Curvature_Higgs_L4[0][0][0][0] = Lh;
  Curvature_Higgs_L4[0][0][1][1] = Lh/3.;
  Curvature_Higgs_L4[0][0][2][2] = Lh/3.;
  Curvature_Higgs_L4[0][0][3][3] = Lh/3.;
  Curvature_Higgs_L4[0][0][4][4] = 2*Lhs;
  Curvature_Higgs_L4[0][1][0][1] = Lh/3.;
  Curvature_Higgs_L4[0][1][1][0] = Lh/3.;
  Curvature_Higgs_L4[0][2][0][2] = Lh/3.;
  Curvature_Higgs_L4[0][2][2][0] = Lh/3.;
  Curvature_Higgs_L4[0][3][0][3] = Lh/3.;
  Curvature_Higgs_L4[0][3][3][0] = Lh/3.;
  Curvature_Higgs_L4[0][4][0][4] = 2*Lhs;
  Curvature_Higgs_L4[0][4][4][0] = 2*Lhs;
  Curvature_Higgs_L4[1][0][0][1] = Lh/3.;
  Curvature_Higgs_L4[1][0][1][0] = Lh/3.;
  Curvature_Higgs_L4[1][1][0][0] = Lh/3.;
  Curvature_Higgs_L4[1][1][1][1] = Lh;
  Curvature_Higgs_L4[1][1][2][2] = Lh/3.;
  Curvature_Higgs_L4[1][1][3][3] = Lh/3.;
  Curvature_Higgs_L4[1][1][4][4] = 2*Lhs;
  Curvature_Higgs_L4[1][2][1][2] = Lh/3.;
  Curvature_Higgs_L4[1][2][2][1] = Lh/3.;
  Curvature_Higgs_L4[1][3][1][3] = Lh/3.;
  Curvature_Higgs_L4[1][3][3][1] = Lh/3.;
  Curvature_Higgs_L4[1][4][1][4] = 2*Lhs;
  Curvature_Higgs_L4[1][4][4][1] = 2*Lhs;
  Curvature_Higgs_L4[2][0][0][2] = Lh/3.;
  Curvature_Higgs_L4[2][0][2][0] = Lh/3.;
  Curvature_Higgs_L4[2][1][1][2] = Lh/3.;
  Curvature_Higgs_L4[2][1][2][1] = Lh/3.;
  Curvature_Higgs_L4[2][2][0][0] = Lh/3.;
  Curvature_Higgs_L4[2][2][1][1] = Lh/3.;
  Curvature_Higgs_L4[2][2][2][2] = Lh;
  Curvature_Higgs_L4[2][2][3][3] = Lh/3.;
  Curvature_Higgs_L4[2][2][4][4] = 2*Lhs;
  Curvature_Higgs_L4[2][3][2][3] = Lh/3.;
  Curvature_Higgs_L4[2][3][3][2] = Lh/3.;
  Curvature_Higgs_L4[2][4][2][4] = 2*Lhs;
  Curvature_Higgs_L4[2][4][4][2] = 2*Lhs;
  Curvature_Higgs_L4[3][0][0][3] = Lh/3.;
  Curvature_Higgs_L4[3][0][3][0] = Lh/3.;
  Curvature_Higgs_L4[3][1][1][3] = Lh/3.;
  Curvature_Higgs_L4[3][1][3][1] = Lh/3.;
  Curvature_Higgs_L4[3][2][2][3] = Lh/3.;
  Curvature_Higgs_L4[3][2][3][2] = Lh/3.;
  Curvature_Higgs_L4[3][3][0][0] = Lh/3.;
  Curvature_Higgs_L4[3][3][1][1] = Lh/3.;
  Curvature_Higgs_L4[3][3][2][2] = Lh/3.;
  Curvature_Higgs_L4[3][3][3][3] = Lh;
  Curvature_Higgs_L4[3][3][4][4] = 2*Lhs;
  Curvature_Higgs_L4[3][4][3][4] = 2*Lhs;
  Curvature_Higgs_L4[3][4][4][3] = 2*Lhs;
  Curvature_Higgs_L4[4][0][0][4] = 2*Lhs;
  Curvature_Higgs_L4[4][0][4][0] = 2*Lhs;
  Curvature_Higgs_L4[4][1][1][4] = 2*Lhs;
  Curvature_Higgs_L4[4][1][4][1] = 2*Lhs;
  Curvature_Higgs_L4[4][2][2][4] = 2*Lhs;
  Curvature_Higgs_L4[4][2][4][2] = 2*Lhs;
  Curvature_Higgs_L4[4][3][3][4] = 2*Lhs;
  Curvature_Higgs_L4[4][3][4][3] = 2*Lhs;
  Curvature_Higgs_L4[4][4][0][0] = 2*Lhs;
  Curvature_Higgs_L4[4][4][1][1] = 2*Lhs;
  Curvature_Higgs_L4[4][4][2][2] = 2*Lhs;
  Curvature_Higgs_L4[4][4][3][3] = 2*Lhs;
  Curvature_Higgs_L4[4][4][4][4] = Ls;

  Curvature_Gauge_G2H2[0][0][0][0] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[0][0][1][1] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[0][0][2][2] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[0][0][3][3] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[0][3][0][2] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[0][3][1][3] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[0][3][2][0] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[0][3][3][1] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[1][1][0][0] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[1][1][1][1] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[1][1][2][2] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[1][1][3][3] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[1][3][0][3] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[1][3][1][2] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][2][1] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][3][0] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[2][2][0][0] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[2][2][1][1] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[2][2][2][2] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[2][2][3][3] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[2][3][0][0] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[2][3][1][1] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[2][3][2][2] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][3][3] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][0][0][2] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][0][1][3] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][0][2][0] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][0][3][1] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][1][0][3] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][1][1][2] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][2][1] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][3][0] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][2][0][0] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][2][1][1] = SMConstants.C_g*SMConstants.C_gs;
  Curvature_Gauge_G2H2[3][2][2][2] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][3][3] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][3][0][0] = pow(SMConstants.C_gs,2);
  Curvature_Gauge_G2H2[3][3][1][1] = pow(SMConstants.C_gs,2);
  Curvature_Gauge_G2H2[3][3][2][2] = pow(SMConstants.C_gs,2);
  Curvature_Gauge_G2H2[3][3][3][3] = pow(SMConstants.C_gs,2);

  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = SMConstants.C_Vud;
  V12 = SMConstants.C_Vus;
  V13 = SMConstants.C_Vub;
  V21 = SMConstants.C_Vcd;
  V22 = SMConstants.C_Vcs;
  V23 = SMConstants.C_Vcb;
  V31 = SMConstants.C_Vtd;
  V32 = SMConstants.C_Vts;
  V33 = SMConstants.C_Vtb;

  Curvature_Lepton_F2H1[0][1][2] = SMConstants.C_MassElectron/vh;
  Curvature_Lepton_F2H1[0][1][3] = (II*SMConstants.C_MassElectron)/vh;
  Curvature_Lepton_F2H1[1][0][2] = SMConstants.C_MassElectron/vh;
  Curvature_Lepton_F2H1[1][0][3] = (II*SMConstants.C_MassElectron)/vh;
  Curvature_Lepton_F2H1[1][6][0] = SMConstants.C_MassElectron/vh;
  Curvature_Lepton_F2H1[1][6][1] = (II*SMConstants.C_MassElectron)/vh;
  Curvature_Lepton_F2H1[2][3][2] = SMConstants.C_MassMu/vh;
  Curvature_Lepton_F2H1[2][3][3] = (II*SMConstants.C_MassMu)/vh;
  Curvature_Lepton_F2H1[3][2][2] = SMConstants.C_MassMu/vh;
  Curvature_Lepton_F2H1[3][2][3] = (II*SMConstants.C_MassMu)/vh;
  Curvature_Lepton_F2H1[3][7][0] = SMConstants.C_MassMu/vh;
  Curvature_Lepton_F2H1[3][7][1] = (II*SMConstants.C_MassMu)/vh;
  Curvature_Lepton_F2H1[4][5][2] = SMConstants.C_MassTau/vh;
  Curvature_Lepton_F2H1[4][5][3] = (II*SMConstants.C_MassTau)/vh;
  Curvature_Lepton_F2H1[5][4][2] = SMConstants.C_MassTau/vh;
  Curvature_Lepton_F2H1[5][4][3] = (II*SMConstants.C_MassTau)/vh;
  Curvature_Lepton_F2H1[5][8][0] = SMConstants.C_MassTau/vh;
  Curvature_Lepton_F2H1[5][8][1] = (II*SMConstants.C_MassTau)/vh;
  Curvature_Lepton_F2H1[6][1][0] = SMConstants.C_MassElectron/vh;
  Curvature_Lepton_F2H1[6][1][1] = (II*SMConstants.C_MassElectron)/vh;
  Curvature_Lepton_F2H1[7][3][0] = SMConstants.C_MassMu/vh;
  Curvature_Lepton_F2H1[7][3][1] = (II*SMConstants.C_MassMu)/vh;
  Curvature_Lepton_F2H1[8][5][0] = SMConstants.C_MassTau/vh;
  Curvature_Lepton_F2H1[8][5][1] = (II*SMConstants.C_MassTau)/vh;

  Curvature_Quark_F2H1[0][6][2] = SMConstants.C_MassUp/vh;
  Curvature_Quark_F2H1[0][6][3] = (-II*SMConstants.C_MassUp)/vh;
  Curvature_Quark_F2H1[0][9][0] = -((SMConstants.C_MassUp*conj(V11))/vh);
  Curvature_Quark_F2H1[0][9][1] = (II*SMConstants.C_MassUp*conj(V11))/vh;
  Curvature_Quark_F2H1[0][10][0] = -((SMConstants.C_MassUp*conj(V12))/vh);
  Curvature_Quark_F2H1[0][10][1] = (II*SMConstants.C_MassUp*conj(V12))/vh;
  Curvature_Quark_F2H1[0][11][0] = -((SMConstants.C_MassUp*conj(V13))/vh);
  Curvature_Quark_F2H1[0][11][1] = (II*SMConstants.C_MassUp*conj(V13))/vh;
  Curvature_Quark_F2H1[1][7][2] = SMConstants.C_MassCharm/vh;
  Curvature_Quark_F2H1[1][7][3] = (-II*SMConstants.C_MassCharm)/vh;
  Curvature_Quark_F2H1[1][9][0] = -((SMConstants.C_MassCharm*conj(V21))/vh);
  Curvature_Quark_F2H1[1][9][1] = (II*SMConstants.C_MassCharm*conj(V21))/vh;
  Curvature_Quark_F2H1[1][10][0] = -((SMConstants.C_MassCharm*conj(V22))/vh);
  Curvature_Quark_F2H1[1][10][1] = (II*SMConstants.C_MassCharm*conj(V22))/vh;
  Curvature_Quark_F2H1[1][11][0] = -((SMConstants.C_MassCharm*conj(V23))/vh);
  Curvature_Quark_F2H1[1][11][1] = (II*SMConstants.C_MassCharm*conj(V23))/vh;
  Curvature_Quark_F2H1[2][8][2] = SMConstants.C_MassTop/vh;
  Curvature_Quark_F2H1[2][8][3] = (-II*SMConstants.C_MassTop)/vh;
  Curvature_Quark_F2H1[2][9][0] = -((SMConstants.C_MassTop*conj(V31))/vh);
  Curvature_Quark_F2H1[2][9][1] = (II*SMConstants.C_MassTop*conj(V31))/vh;
  Curvature_Quark_F2H1[2][10][0] = -((SMConstants.C_MassTop*conj(V32))/vh);
  Curvature_Quark_F2H1[2][10][1] = (II*SMConstants.C_MassTop*conj(V32))/vh;
  Curvature_Quark_F2H1[2][11][0] = -((SMConstants.C_MassTop*conj(V33))/vh);
  Curvature_Quark_F2H1[2][11][1] = (II*SMConstants.C_MassTop*conj(V33))/vh;
  Curvature_Quark_F2H1[3][6][0] = (SMConstants.C_MassDown*V11)/vh;
  Curvature_Quark_F2H1[3][6][1] = (II*SMConstants.C_MassDown*V11)/vh;
  Curvature_Quark_F2H1[3][7][0] = (SMConstants.C_MassDown*V21)/vh;
  Curvature_Quark_F2H1[3][7][1] = (II*SMConstants.C_MassDown*V21)/vh;
  Curvature_Quark_F2H1[3][8][0] = (SMConstants.C_MassDown*V31)/vh;
  Curvature_Quark_F2H1[3][8][1] = (II*SMConstants.C_MassDown*V31)/vh;
  Curvature_Quark_F2H1[3][9][2] = SMConstants.C_MassDown/vh;
  Curvature_Quark_F2H1[3][9][3] = (II*SMConstants.C_MassDown)/vh;
  Curvature_Quark_F2H1[4][6][0] = (SMConstants.C_MassStrange*V12)/vh;
  Curvature_Quark_F2H1[4][6][1] = (II*SMConstants.C_MassStrange*V12)/vh;
  Curvature_Quark_F2H1[4][7][0] = (SMConstants.C_MassStrange*V22)/vh;
  Curvature_Quark_F2H1[4][7][1] = (II*SMConstants.C_MassStrange*V22)/vh;
  Curvature_Quark_F2H1[4][8][0] = (SMConstants.C_MassStrange*V32)/vh;
  Curvature_Quark_F2H1[4][8][1] = (II*SMConstants.C_MassStrange*V32)/vh;
  Curvature_Quark_F2H1[4][10][2] = SMConstants.C_MassStrange/vh;
  Curvature_Quark_F2H1[4][10][3] = (II*SMConstants.C_MassStrange)/vh;
  Curvature_Quark_F2H1[5][6][0] = (SMConstants.C_MassBottom*V13)/vh;
  Curvature_Quark_F2H1[5][6][1] = (II*SMConstants.C_MassBottom*V13)/vh;
  Curvature_Quark_F2H1[5][7][0] = (SMConstants.C_MassBottom*V23)/vh;
  Curvature_Quark_F2H1[5][7][1] = (II*SMConstants.C_MassBottom*V23)/vh;
  Curvature_Quark_F2H1[5][8][0] = (SMConstants.C_MassBottom*V33)/vh;
  Curvature_Quark_F2H1[5][8][1] = (II*SMConstants.C_MassBottom*V33)/vh;
  Curvature_Quark_F2H1[5][11][2] = SMConstants.C_MassBottom/vh;
  Curvature_Quark_F2H1[5][11][3] = (II*SMConstants.C_MassBottom)/vh;
  Curvature_Quark_F2H1[6][0][2] = SMConstants.C_MassUp/vh;
  Curvature_Quark_F2H1[6][0][3] = (-II*SMConstants.C_MassUp)/vh;
  Curvature_Quark_F2H1[6][3][0] = (SMConstants.C_MassDown*V11)/vh;
  Curvature_Quark_F2H1[6][3][1] = (II*SMConstants.C_MassDown*V11)/vh;
  Curvature_Quark_F2H1[6][4][0] = (SMConstants.C_MassStrange*V12)/vh;
  Curvature_Quark_F2H1[6][4][1] = (II*SMConstants.C_MassStrange*V12)/vh;
  Curvature_Quark_F2H1[6][5][0] = (SMConstants.C_MassBottom*V13)/vh;
  Curvature_Quark_F2H1[6][5][1] = (II*SMConstants.C_MassBottom*V13)/vh;
  Curvature_Quark_F2H1[7][1][2] = SMConstants.C_MassCharm/vh;
  Curvature_Quark_F2H1[7][1][3] = (-II*SMConstants.C_MassCharm)/vh;
  Curvature_Quark_F2H1[7][3][0] = (SMConstants.C_MassDown*V21)/vh;
  Curvature_Quark_F2H1[7][3][1] = (II*SMConstants.C_MassDown*V21)/vh;
  Curvature_Quark_F2H1[7][4][0] = (SMConstants.C_MassStrange*V22)/vh;
  Curvature_Quark_F2H1[7][4][1] = (II*SMConstants.C_MassStrange*V22)/vh;
  Curvature_Quark_F2H1[7][5][0] = (SMConstants.C_MassBottom*V23)/vh;
  Curvature_Quark_F2H1[7][5][1] = (II*SMConstants.C_MassBottom*V23)/vh;
  Curvature_Quark_F2H1[8][2][2] = SMConstants.C_MassTop/vh;
  Curvature_Quark_F2H1[8][2][3] = (-II*SMConstants.C_MassTop)/vh;
  Curvature_Quark_F2H1[8][3][0] = (SMConstants.C_MassDown*V31)/vh;
  Curvature_Quark_F2H1[8][3][1] = (II*SMConstants.C_MassDown*V31)/vh;
  Curvature_Quark_F2H1[8][4][0] = (SMConstants.C_MassStrange*V32)/vh;
  Curvature_Quark_F2H1[8][4][1] = (II*SMConstants.C_MassStrange*V32)/vh;
  Curvature_Quark_F2H1[8][5][0] = (SMConstants.C_MassBottom*V33)/vh;
  Curvature_Quark_F2H1[8][5][1] = (II*SMConstants.C_MassBottom*V33)/vh;
  Curvature_Quark_F2H1[9][0][0] = -((SMConstants.C_MassUp*conj(V11))/vh);
  Curvature_Quark_F2H1[9][0][1] = (II*SMConstants.C_MassUp*conj(V11))/vh;
  Curvature_Quark_F2H1[9][1][0] = -((SMConstants.C_MassCharm*conj(V21))/vh);
  Curvature_Quark_F2H1[9][1][1] = (II*SMConstants.C_MassCharm*conj(V21))/vh;
  Curvature_Quark_F2H1[9][2][0] = -((SMConstants.C_MassTop*conj(V31))/vh);
  Curvature_Quark_F2H1[9][2][1] = (II*SMConstants.C_MassTop*conj(V31))/vh;
  Curvature_Quark_F2H1[9][3][2] = SMConstants.C_MassDown/vh;
  Curvature_Quark_F2H1[9][3][3] = (II*SMConstants.C_MassDown)/vh;
  Curvature_Quark_F2H1[10][0][0] = -((SMConstants.C_MassUp*conj(V12))/vh);
  Curvature_Quark_F2H1[10][0][1] = (II*SMConstants.C_MassUp*conj(V12))/vh;
  Curvature_Quark_F2H1[10][1][0] = -((SMConstants.C_MassCharm*conj(V22))/vh);
  Curvature_Quark_F2H1[10][1][1] = (II*SMConstants.C_MassCharm*conj(V22))/vh;
  Curvature_Quark_F2H1[10][2][0] = -((SMConstants.C_MassTop*conj(V32))/vh);
  Curvature_Quark_F2H1[10][2][1] = (II*SMConstants.C_MassTop*conj(V32))/vh;
  Curvature_Quark_F2H1[10][4][2] = SMConstants.C_MassStrange/vh;
  Curvature_Quark_F2H1[10][4][3] = (II*SMConstants.C_MassStrange)/vh;
  Curvature_Quark_F2H1[11][0][0] = -((SMConstants.C_MassUp*conj(V13))/vh);
  Curvature_Quark_F2H1[11][0][1] = (II*SMConstants.C_MassUp*conj(V13))/vh;
  Curvature_Quark_F2H1[11][1][0] = -((SMConstants.C_MassCharm*conj(V23))/vh);
  Curvature_Quark_F2H1[11][1][1] = (II*SMConstants.C_MassCharm*conj(V23))/vh;
  Curvature_Quark_F2H1[11][2][0] = -((SMConstants.C_MassTop*conj(V33))/vh);
  Curvature_Quark_F2H1[11][2][1] = (II*SMConstants.C_MassTop*conj(V33))/vh;
  Curvature_Quark_F2H1[11][5][2] = SMConstants.C_MassBottom/vh;
  Curvature_Quark_F2H1[11][5][3] = (II*SMConstants.C_MassBottom)/vh;

}

bool Class_Potential_RxSM_OS_lagparams::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_RxSM_OS_lagparams::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */
  return false;
}
double
Class_Potential_RxSM_OS_lagparams::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (Lh*pow(v[2],4) - 12*ms*pow(v[4],2) + Ls*pow(v[4],4) - 12*pow(v[2],2)*(mh - Lhs*pow(v[4],2)))/24.;
  
  return res;
}

double Class_Potential_RxSM_OS_lagparams::VCounterSimplified(
    const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = dT3*v[2] + dT5*v[4] + (dLh*pow(v[2],4) - 12*dms*pow(v[4],2) + dLs*pow(v[4],4) - 12*pow(v[2],2)*(dmh - dLhs*pow(v[4],2)))/24.;
  
  return res;
}

void Class_Potential_RxSM_OS_lagparams::Debugging(const std::vector<double> &input,
                                            std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
