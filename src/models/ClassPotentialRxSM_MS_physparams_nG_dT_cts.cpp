#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for SMConstants.C_vev0, SMConstants.C_MassTop, SMConstants.C_g
#include <algorithm> // for max, copy
#include <iomanip>
#include <iostream> // for operator<<, endl, basic_o...
#include <memory>   // for allocator_traits<>::value...
#include <stddef.h> // for std::size_t

#include <BSMPT/models/ClassPotentialRxSM_MS_physparams_nG_dT_cts.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_RxSM_MS_physparams_nG_dT_cts::Class_Potential_RxSM_MS_physparams_nG_dT_cts(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::RXSM_MS_PHYSPARAMS_NG_DT_CTS;

  nPar = 9;   // number of parameters in the tree-Level Lagrangian AFTER using
               // tadpole equations
  nParCT = 10; // number of parameters in the counterterm potential

  nVEV = 2; // number of VEVs to minimize the potential

  NHiggs = 2; // number of scalar d.o.f.

  NGauge = 4; // number of gauge fields

  NLepton = 9; // number of lepton fields

  NQuarks = 12; // number of quark fields

  VevOrder.resize(nVEV);
  VevOrder[0] = 0; // wh
  VevOrder[1] = 1; // ws

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_RxSM_MS_physparams_nG_dT_cts::~Class_Potential_RxSM_MS_physparams_nG_dT_cts()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_RxSM_MS_physparams_nG_dT_cts::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dmh");
  labels.push_back("dms");
  labels.push_back("dLh");
  labels.push_back("dLs");
  labels.push_back("dLhs");
  labels.push_back("ddLhOS");
  labels.push_back("ddLsOS");
  labels.push_back("ddLhsOS");
  labels.push_back("dT1");
  labels.push_back("dT2");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Potential_RxSM_MS_physparams_nG_dT_cts::addLegendTemp() const
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
Class_Potential_RxSM_MS_physparams_nG_dT_cts::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  // mass basis, you can identify here your particles
  particles.push_back("h_1");
  particles.push_back("h_2");

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
std::vector<std::string> Class_Potential_RxSM_MS_physparams_nG_dT_cts::addLegendVEV() const
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
void Class_Potential_RxSM_MS_physparams_nG_dT_cts::ReadAndSet(const std::string &linestr,
                                             std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 9; k++)
  {
    ss >> tmp;
    if (k == 1)
      par[0] = tmp; // m1
    if (k == 2)
      par[1] = tmp; // m2
    if (k == 3)
      par[2] = tmp; // vh
    if (k == 4)
      par[3] = tmp; // vs
    if (k == 5)
      par[4] = tmp; // a
    if (k == 6)
      par[5] = tmp; // dLhOS
    if (k == 7)
      par[6] = tmp; // dLsOS
    if (k == 8)
      par[7] = tmp; // dLhsOS
    if (k == 9)
      par[8] = tmp; // mu

  }

  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_RxSM_MS_physparams_nG_dT_cts::set_gen(const std::vector<double> &par)
{

  m1 = par[0]; 
  m2 = par[1]; 
  vh = par[2]; 
  vs = par[3]; 
  a = par[4]; 
  dLhOS = par[5]; 
  dLsOS = par[6]; 
  dLhsOS = par[7]; 

  mh = (dLhOS*pow(vh,2) + 6*dLhsOS*pow(vs,2) + 3*pow(m1,2)*pow(cos(a),2) - (3*(pow(m1,2) - pow(m2,2))*vs*cos(a)*sin(a))/vh + 3*pow(m2,2)*pow(sin(a),2))/6.; 
  ms = (6*dLhsOS*pow(vh,2) + dLsOS*pow(vs,2) + 3*pow(m2,2)*pow(cos(a),2) - (3*(pow(m1,2) - pow(m2,2))*vh*cos(a)*sin(a))/vs + 3*pow(m1,2)*pow(sin(a),2))/6.; 
  Lh = (3*(pow(m1,2)*pow(cos(a),2) + pow(m2,2)*pow(sin(a),2)))/pow(vh,2); 
  Ls = (3*(pow(m2,2)*pow(cos(a),2) + pow(m1,2)*pow(sin(a),2)))/pow(vs,2); 
  Lhs = -0.5*((pow(m1,2) - pow(m2,2))*cos(a)*sin(a))/(vh*vs); 

  scale = par[8]; // renormalisation scale is given as input

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
void Class_Potential_RxSM_MS_physparams_nG_dT_cts::set_CT_Pot_Par(const std::vector<double> &par)
{
  dmh = par[0];
  dms = par[1];
  dLh = par[2];
  dLs = par[3];
  dLhs = par[4];
  ddLhOS = par[5];
  ddLsOS = par[6];
  ddLhsOS = par[7];
  dT1 = par[8];
  dT2 = par[9];

  // assign the non-zero entries
  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;




}

/**
 * console output of all parameters
 */
void Class_Potential_RxSM_MS_physparams_nG_dT_cts::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "Model = " << Model << "\n";

  ss << "\nThe input parameters are : \n";
  ss << "m1 = " << m1 << "\n";
  ss << "m2 = " << m2 << "\n";
  ss << "vh = " << vh << "\n";
  ss << "vs = " << vs << "\n";
  ss << "a = " << a << "\n";
  ss << "dLhOS = " << dLhOS << "\n";
  ss << "dLsOS = " << dLsOS << "\n";
  ss << "dLhsOS = " << dLhsOS << "\n";

  ss << "\nThe parameters are : \n";
  ss << "mh = " << mh << "\n";
  ss << "ms = " << ms << "\n";
  ss << "Lh = " << Lh << "\n";
  ss << "Ls = " << Ls << "\n";
  ss << "Lhs = " << Lhs << "\n";
  ss << "dLhOS = " << dLhOS << "\n";
  ss << "dLsOS = " << dLsOS << "\n";
  ss << "dLhsOS = " << dLhsOS << "\n";

  ss << "\nThe counterterm parameters are : \n";
  ss << "dmh = " << dmh << "\n";
  ss << "dms = " << dms << "\n";
  ss << "dLh = " << dLh << "\n";
  ss << "dLs = " << dLs << "\n";
  ss << "dLhs = " << dLhs << "\n";
  ss << "ddLhOS = " << ddLhOS << "\n";
  ss << "ddLsOS = " << ddLsOS << "\n";
  ss << "ddLhsOS = " << ddLhsOS << "\n";
  ss << "dT1 = " << dT1 << "\n";
  ss << "dT2 = " << dT2 << "\n";

  ss << "\nThe scale is given by mu = " << scale << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Potential_RxSM_MS_physparams_nG_dT_cts::calc_CT() const
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
  parCT.push_back(0); //dmh;
  parCT.push_back(0); //dms;
  parCT.push_back(0); //dLh;
  parCT.push_back(0); //dLs;
  parCT.push_back(0); //dLhs;
  parCT.push_back(0); //ddLhOS;
  parCT.push_back(0); //ddLsOS;
  parCT.push_back(0); //ddLhsOS;
  parCT.push_back(-NablaWeinberg(0)); //dT1;
  parCT.push_back(-NablaWeinberg(1)); //dT2;

  return parCT;
}

// mass basis triple couplings
void Class_Potential_RxSM_MS_physparams_nG_dT_cts::TripleHiggsCouplings()
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

void Class_Potential_RxSM_MS_physparams_nG_dT_cts::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // assign the non-zero entries
  Curvature_Higgs_L2[0][0] = -mh;
  Curvature_Higgs_L2[1][1] = -ms;

  Curvature_Higgs_L4[0][0][0][0] = dLhOS + Lh;
  Curvature_Higgs_L4[0][0][1][1] = 2*(dLhsOS + Lhs);
  Curvature_Higgs_L4[0][1][0][1] = 2*(dLhsOS + Lhs);
  Curvature_Higgs_L4[0][1][1][0] = 2*(dLhsOS + Lhs);
  Curvature_Higgs_L4[1][0][0][1] = 2*(dLhsOS + Lhs);
  Curvature_Higgs_L4[1][0][1][0] = 2*(dLhsOS + Lhs);
  Curvature_Higgs_L4[1][1][0][0] = 2*(dLhsOS + Lhs);
  Curvature_Higgs_L4[1][1][1][1] = dLsOS + Ls;

  Curvature_Gauge_G2H2[0][0][0][0] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[1][1][0][0] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[2][2][0][0] = pow(SMConstants.C_g,2);
  Curvature_Gauge_G2H2[2][3][0][0] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][0][0] = -(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][3][0][0] = pow(SMConstants.C_gs,2);

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

  Curvature_Lepton_F2H1[1][6][0] = SMConstants.C_MassElectron/vh;
  Curvature_Lepton_F2H1[3][7][0] = SMConstants.C_MassMu/vh;
  Curvature_Lepton_F2H1[5][8][0] = SMConstants.C_MassTau/vh;
  Curvature_Lepton_F2H1[6][1][0] = SMConstants.C_MassElectron/vh;
  Curvature_Lepton_F2H1[7][3][0] = SMConstants.C_MassMu/vh;
  Curvature_Lepton_F2H1[8][5][0] = SMConstants.C_MassTau/vh;

  Curvature_Quark_F2H1[0][6][0] = -(SMConstants.C_MassUp/vh);
  Curvature_Quark_F2H1[1][7][0] = -(SMConstants.C_MassCharm/vh);
  Curvature_Quark_F2H1[2][8][0] = -(SMConstants.C_MassTop/vh);
  Curvature_Quark_F2H1[3][9][0] = SMConstants.C_MassDown/vh;
  Curvature_Quark_F2H1[4][10][0] = SMConstants.C_MassStrange/vh;
  Curvature_Quark_F2H1[5][11][0] = SMConstants.C_MassBottom/vh;
  Curvature_Quark_F2H1[6][0][0] = -(SMConstants.C_MassUp/vh);
  Curvature_Quark_F2H1[7][1][0] = -(SMConstants.C_MassCharm/vh);
  Curvature_Quark_F2H1[8][2][0] = -(SMConstants.C_MassTop/vh);
  Curvature_Quark_F2H1[9][3][0] = SMConstants.C_MassDown/vh;
  Curvature_Quark_F2H1[10][4][0] = SMConstants.C_MassStrange/vh;
  Curvature_Quark_F2H1[11][5][0] = SMConstants.C_MassBottom/vh;

}

bool Class_Potential_RxSM_MS_physparams_nG_dT_cts::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_RxSM_MS_physparams_nG_dT_cts::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */
  return false;
}
double
Class_Potential_RxSM_MS_physparams_nG_dT_cts::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (pow(v[0],4)*(dLhOS + Lh) + pow(v[1],4)*(dLsOS + Ls) + 12*pow(v[0],2)*(pow(v[1],2)*(dLhsOS + Lhs) - mh) - 12*pow(v[1],2)*ms)/24.;
  
  return res;
}

double Class_Potential_RxSM_MS_physparams_nG_dT_cts::VCounterSimplified(
    const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = v[0]*dT1 + v[1]*dT2;
  
  return res;
}

void Class_Potential_RxSM_MS_physparams_nG_dT_cts::Debugging(const std::vector<double> &input,
                                            std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
