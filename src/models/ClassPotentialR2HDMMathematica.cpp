#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for SMConstants.C_vev0, SMConstants.C_MassTop, SMConstants.C_g
#include <algorithm> // for max, copy
#include <iomanip>
#include <iostream> // for operator<<, endl, basic_o...
#include <memory>   // for allocator_traits<>::value...
#include <stddef.h> // for std::size_t

#include <BSMPT/models/ClassPotentialR2HDMMathematica.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_R2HDMMathematica::Class_Potential_R2HDMMathematica(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::R2HDMMATHEMATICA;

  nPar = 8;   // number of parameters in the tree-Level Lagrangian AFTER using
               // tadpole equations
  nParCT = 16; // number of parameters in the counterterm potential

  nVEV = 4; // number of VEVs to minimize the potential

  NHiggs = 8; // number of scalar d.o.f.

  NGauge = 4; // number of gauge fields

  NLepton = 9; // number of lepton fields

  NQuarks = 12; // number of quark fields

  VevOrder.resize(nVEV);
  VevOrder[0] = 2; // wcb
  VevOrder[1] = 4; // w1
  VevOrder[2] = 6; // w2
  VevOrder[3] = 7; // wcp

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_R2HDMMathematica::~Class_Potential_R2HDMMathematica()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_R2HDMMathematica::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("dm11Sq");
  labels.push_back("dm22Sq");
  labels.push_back("dm12Sq");
  labels.push_back("dL1");
  labels.push_back("dL2");
  labels.push_back("dL3");
  labels.push_back("dL4");
  labels.push_back("dL5");
  labels.push_back("dT1");
  labels.push_back("dT2");
  labels.push_back("dT3");
  labels.push_back("dT4");
  labels.push_back("dT5");
  labels.push_back("dT6");
  labels.push_back("dT7");
  labels.push_back("dT8");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Potential_R2HDMMathematica::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");     // Label for the critical temperature
  labels.push_back("v_c");     // Label for the critical vev
  labels.push_back("v_c/T_c"); // Label for xi_c
  // out += "VEV order";
  labels.push_back("wcb(T_c)");
  labels.push_back("w1(T_c)");
  labels.push_back("w2(T_c)");
  labels.push_back("wcp(T_c)");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 */
std::vector<std::string>
Class_Potential_R2HDMMathematica::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  // mass basis, you can identify here your particles
  particles.push_back("h_1");
  particles.push_back("h_2");
  particles.push_back("h_3");
  particles.push_back("h_4");
  particles.push_back("h_5");
  particles.push_back("h_6");
  particles.push_back("h_7");
  particles.push_back("h_8");

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
std::vector<std::string> Class_Potential_R2HDMMathematica::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order"
  labels.push_back("wcb");
  labels.push_back("w1");
  labels.push_back("w2");
  labels.push_back("wcp");

  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_R2HDMMathematica::ReadAndSet(const std::string &linestr,
                                             std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 8; k++)
  {
    ss >> tmp;
    if (k == 1)
      par[0] = tmp; // type
    if (k == 2)
      par[1] = tmp; // L1
    if (k == 3)
      par[2] = tmp; // L2
    if (k == 4)
      par[3] = tmp; // L3
    if (k == 5)
      par[4] = tmp; // L4
    if (k == 6)
      par[5] = tmp; // L5
    if (k == 7)
      par[6] = tmp; // m12Sq
    if (k == 8)
      par[7] = tmp; // tanbeta

  }

  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_R2HDMMathematica::set_gen(const std::vector<double> &par)
{

  type = par[0]; 
  L1 = par[1]; 
  L2 = par[2]; 
  L3 = par[3]; 
  L4 = par[4]; 
  L5 = par[5]; 
  m12Sq = par[6]; 
  tanbeta = par[7]; 

  m11Sq = (2*m12Sq*(tanbeta + pow(tanbeta,3)) - (L1 + (L3 + L4 + L5)*pow(tanbeta,2))*pow(SMConstants.C_vev0,2))/(2.*(1 + pow(tanbeta,2))); 
  m22Sq = m12Sq/tanbeta - ((L3 + L4 + L5 + L2*pow(tanbeta,2))*pow(SMConstants.C_vev0,2))/(2.*(1 + pow(tanbeta,2))); 
  v1 = SMConstants.C_vev0/sqrt(1 + pow(tanbeta,2)); 
  v2 = (tanbeta*SMConstants.C_vev0)/sqrt(1 + pow(tanbeta,2)); 

  scale = SMConstants.C_vev0; // renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // set the vector vevTreeMin. vevTree will then be set by the
  // function MinimizeOrderVEV
  vevTreeMin[0] = 0; // wcb
  vevTreeMin[1] = v1; // w1
  vevTreeMin[2] = v2; // w2
  vevTreeMin[3] = 0; // wcp

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_Potential_R2HDMMathematica::set_CT_Pot_Par(const std::vector<double> &par)
{
  dm11Sq = par[0];
  dm22Sq = par[1];
  dm12Sq = par[2];
  dL1 = par[3];
  dL2 = par[4];
  dL3 = par[5];
  dL4 = par[6];
  dL5 = par[7];
  dT1 = par[8];
  dT2 = par[9];
  dT3 = par[10];
  dT4 = par[11];
  dT5 = par[12];
  dT6 = par[13];
  dT7 = par[14];
  dT8 = par[15];

  // assign the non-zero entries
  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;
  Curvature_Higgs_CT_L1[5] = dT6;
  Curvature_Higgs_CT_L1[6] = dT7;
  Curvature_Higgs_CT_L1[7] = dT8;

  Curvature_Higgs_CT_L2[0][0] = dm11Sq;
  Curvature_Higgs_CT_L2[0][2] = -dm12Sq;
  Curvature_Higgs_CT_L2[1][1] = dm11Sq;
  Curvature_Higgs_CT_L2[1][3] = -dm12Sq;
  Curvature_Higgs_CT_L2[2][0] = -dm12Sq;
  Curvature_Higgs_CT_L2[2][2] = dm22Sq;
  Curvature_Higgs_CT_L2[3][1] = -dm12Sq;
  Curvature_Higgs_CT_L2[3][3] = dm22Sq;
  Curvature_Higgs_CT_L2[4][4] = dm11Sq;
  Curvature_Higgs_CT_L2[4][6] = -dm12Sq;
  Curvature_Higgs_CT_L2[5][5] = dm11Sq;
  Curvature_Higgs_CT_L2[5][7] = -dm12Sq;
  Curvature_Higgs_CT_L2[6][4] = -dm12Sq;
  Curvature_Higgs_CT_L2[6][6] = dm22Sq;
  Curvature_Higgs_CT_L2[7][5] = -dm12Sq;
  Curvature_Higgs_CT_L2[7][7] = dm22Sq;

  Curvature_Higgs_CT_L4[0][0][0][0] = 3*dL1;
  Curvature_Higgs_CT_L4[0][0][1][1] = dL1;
  Curvature_Higgs_CT_L4[0][0][2][2] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[0][0][3][3] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[0][0][4][4] = dL1;
  Curvature_Higgs_CT_L4[0][0][5][5] = dL1;
  Curvature_Higgs_CT_L4[0][0][6][6] = dL3;
  Curvature_Higgs_CT_L4[0][0][7][7] = dL3;
  Curvature_Higgs_CT_L4[0][1][0][1] = dL1;
  Curvature_Higgs_CT_L4[0][1][1][0] = dL1;
  Curvature_Higgs_CT_L4[0][1][2][3] = dL5;
  Curvature_Higgs_CT_L4[0][1][3][2] = dL5;
  Curvature_Higgs_CT_L4[0][2][0][2] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[0][2][1][3] = dL5;
  Curvature_Higgs_CT_L4[0][2][2][0] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[0][2][3][1] = dL5;
  Curvature_Higgs_CT_L4[0][2][4][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][2][5][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][2][6][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][2][7][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][3][0][3] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[0][3][1][2] = dL5;
  Curvature_Higgs_CT_L4[0][3][2][1] = dL5;
  Curvature_Higgs_CT_L4[0][3][3][0] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[0][3][4][7] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][3][5][6] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][3][6][5] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][3][7][4] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][4][0][4] = dL1;
  Curvature_Higgs_CT_L4[0][4][2][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][4][3][7] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][4][4][0] = dL1;
  Curvature_Higgs_CT_L4[0][4][6][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][4][7][3] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][5][0][5] = dL1;
  Curvature_Higgs_CT_L4[0][5][2][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][5][3][6] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][5][5][0] = dL1;
  Curvature_Higgs_CT_L4[0][5][6][3] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][5][7][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][6][0][6] = dL3;
  Curvature_Higgs_CT_L4[0][6][2][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][6][3][5] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][6][4][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][6][5][3] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][6][6][0] = dL3;
  Curvature_Higgs_CT_L4[0][7][0][7] = dL3;
  Curvature_Higgs_CT_L4[0][7][2][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][7][3][4] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][7][4][3] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][7][5][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[0][7][7][0] = dL3;
  Curvature_Higgs_CT_L4[1][0][0][1] = dL1;
  Curvature_Higgs_CT_L4[1][0][1][0] = dL1;
  Curvature_Higgs_CT_L4[1][0][2][3] = dL5;
  Curvature_Higgs_CT_L4[1][0][3][2] = dL5;
  Curvature_Higgs_CT_L4[1][1][0][0] = dL1;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3*dL1;
  Curvature_Higgs_CT_L4[1][1][2][2] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[1][1][3][3] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[1][1][4][4] = dL1;
  Curvature_Higgs_CT_L4[1][1][5][5] = dL1;
  Curvature_Higgs_CT_L4[1][1][6][6] = dL3;
  Curvature_Higgs_CT_L4[1][1][7][7] = dL3;
  Curvature_Higgs_CT_L4[1][2][0][3] = dL5;
  Curvature_Higgs_CT_L4[1][2][1][2] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[1][2][2][1] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[1][2][3][0] = dL5;
  Curvature_Higgs_CT_L4[1][2][4][7] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][2][5][6] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][2][6][5] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][2][7][4] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][3][0][2] = dL5;
  Curvature_Higgs_CT_L4[1][3][1][3] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[1][3][2][0] = dL5;
  Curvature_Higgs_CT_L4[1][3][3][1] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[1][3][4][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][3][5][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][3][6][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][3][7][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][4][1][4] = dL1;
  Curvature_Higgs_CT_L4[1][4][2][7] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][4][3][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][4][4][1] = dL1;
  Curvature_Higgs_CT_L4[1][4][6][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][4][7][2] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][5][1][5] = dL1;
  Curvature_Higgs_CT_L4[1][5][2][6] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][5][3][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][5][5][1] = dL1;
  Curvature_Higgs_CT_L4[1][5][6][2] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][5][7][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][6][1][6] = dL3;
  Curvature_Higgs_CT_L4[1][6][2][5] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][6][3][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][6][4][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][6][5][2] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][6][6][1] = dL3;
  Curvature_Higgs_CT_L4[1][7][1][7] = dL3;
  Curvature_Higgs_CT_L4[1][7][2][4] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][7][3][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][7][4][2] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][7][5][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[1][7][7][1] = dL3;
  Curvature_Higgs_CT_L4[2][0][0][2] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[2][0][1][3] = dL5;
  Curvature_Higgs_CT_L4[2][0][2][0] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[2][0][3][1] = dL5;
  Curvature_Higgs_CT_L4[2][0][4][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][0][5][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][0][6][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][0][7][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][1][0][3] = dL5;
  Curvature_Higgs_CT_L4[2][1][1][2] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[2][1][2][1] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[2][1][3][0] = dL5;
  Curvature_Higgs_CT_L4[2][1][4][7] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][1][5][6] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][1][6][5] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][1][7][4] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][2][0][0] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[2][2][1][1] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3*dL2;
  Curvature_Higgs_CT_L4[2][2][3][3] = dL2;
  Curvature_Higgs_CT_L4[2][2][4][4] = dL3;
  Curvature_Higgs_CT_L4[2][2][5][5] = dL3;
  Curvature_Higgs_CT_L4[2][2][6][6] = dL2;
  Curvature_Higgs_CT_L4[2][2][7][7] = dL2;
  Curvature_Higgs_CT_L4[2][3][0][1] = dL5;
  Curvature_Higgs_CT_L4[2][3][1][0] = dL5;
  Curvature_Higgs_CT_L4[2][3][2][3] = dL2;
  Curvature_Higgs_CT_L4[2][3][3][2] = dL2;
  Curvature_Higgs_CT_L4[2][4][0][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][4][1][7] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][4][2][4] = dL3;
  Curvature_Higgs_CT_L4[2][4][4][2] = dL3;
  Curvature_Higgs_CT_L4[2][4][6][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][4][7][1] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][5][0][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][5][1][6] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][5][2][5] = dL3;
  Curvature_Higgs_CT_L4[2][5][5][2] = dL3;
  Curvature_Higgs_CT_L4[2][5][6][1] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][5][7][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][6][0][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][6][1][5] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][6][2][6] = dL2;
  Curvature_Higgs_CT_L4[2][6][4][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][6][5][1] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][6][6][2] = dL2;
  Curvature_Higgs_CT_L4[2][7][0][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][7][1][4] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][7][2][7] = dL2;
  Curvature_Higgs_CT_L4[2][7][4][1] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][7][5][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[2][7][7][2] = dL2;
  Curvature_Higgs_CT_L4[3][0][0][3] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[3][0][1][2] = dL5;
  Curvature_Higgs_CT_L4[3][0][2][1] = dL5;
  Curvature_Higgs_CT_L4[3][0][3][0] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[3][0][4][7] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][0][5][6] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][0][6][5] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][0][7][4] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][1][0][2] = dL5;
  Curvature_Higgs_CT_L4[3][1][1][3] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[3][1][2][0] = dL5;
  Curvature_Higgs_CT_L4[3][1][3][1] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[3][1][4][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][1][5][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][1][6][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][1][7][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][2][0][1] = dL5;
  Curvature_Higgs_CT_L4[3][2][1][0] = dL5;
  Curvature_Higgs_CT_L4[3][2][2][3] = dL2;
  Curvature_Higgs_CT_L4[3][2][3][2] = dL2;
  Curvature_Higgs_CT_L4[3][3][0][0] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[3][3][1][1] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[3][3][2][2] = dL2;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3*dL2;
  Curvature_Higgs_CT_L4[3][3][4][4] = dL3;
  Curvature_Higgs_CT_L4[3][3][5][5] = dL3;
  Curvature_Higgs_CT_L4[3][3][6][6] = dL2;
  Curvature_Higgs_CT_L4[3][3][7][7] = dL2;
  Curvature_Higgs_CT_L4[3][4][0][7] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][4][1][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][4][3][4] = dL3;
  Curvature_Higgs_CT_L4[3][4][4][3] = dL3;
  Curvature_Higgs_CT_L4[3][4][6][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][4][7][0] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][5][0][6] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][5][1][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][5][3][5] = dL3;
  Curvature_Higgs_CT_L4[3][5][5][3] = dL3;
  Curvature_Higgs_CT_L4[3][5][6][0] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][5][7][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][6][0][5] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][6][1][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][6][3][6] = dL2;
  Curvature_Higgs_CT_L4[3][6][4][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][6][5][0] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][6][6][3] = dL2;
  Curvature_Higgs_CT_L4[3][7][0][4] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][7][1][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][7][3][7] = dL2;
  Curvature_Higgs_CT_L4[3][7][4][0] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][7][5][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[3][7][7][3] = dL2;
  Curvature_Higgs_CT_L4[4][0][0][4] = dL1;
  Curvature_Higgs_CT_L4[4][0][2][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][0][3][7] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][0][4][0] = dL1;
  Curvature_Higgs_CT_L4[4][0][6][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][0][7][3] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][1][1][4] = dL1;
  Curvature_Higgs_CT_L4[4][1][2][7] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][1][3][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][1][4][1] = dL1;
  Curvature_Higgs_CT_L4[4][1][6][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][1][7][2] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][2][0][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][2][1][7] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][2][2][4] = dL3;
  Curvature_Higgs_CT_L4[4][2][4][2] = dL3;
  Curvature_Higgs_CT_L4[4][2][6][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][2][7][1] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][3][0][7] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][3][1][6] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][3][3][4] = dL3;
  Curvature_Higgs_CT_L4[4][3][4][3] = dL3;
  Curvature_Higgs_CT_L4[4][3][6][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][3][7][0] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][4][0][0] = dL1;
  Curvature_Higgs_CT_L4[4][4][1][1] = dL1;
  Curvature_Higgs_CT_L4[4][4][2][2] = dL3;
  Curvature_Higgs_CT_L4[4][4][3][3] = dL3;
  Curvature_Higgs_CT_L4[4][4][4][4] = 3*dL1;
  Curvature_Higgs_CT_L4[4][4][5][5] = dL1;
  Curvature_Higgs_CT_L4[4][4][6][6] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[4][4][7][7] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[4][5][4][5] = dL1;
  Curvature_Higgs_CT_L4[4][5][5][4] = dL1;
  Curvature_Higgs_CT_L4[4][5][6][7] = dL5;
  Curvature_Higgs_CT_L4[4][5][7][6] = dL5;
  Curvature_Higgs_CT_L4[4][6][0][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][6][1][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][6][2][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][6][3][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][6][4][6] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[4][6][5][7] = dL5;
  Curvature_Higgs_CT_L4[4][6][6][4] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[4][6][7][5] = dL5;
  Curvature_Higgs_CT_L4[4][7][0][3] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][7][1][2] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][7][2][1] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][7][3][0] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[4][7][4][7] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[4][7][5][6] = dL5;
  Curvature_Higgs_CT_L4[4][7][6][5] = dL5;
  Curvature_Higgs_CT_L4[4][7][7][4] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[5][0][0][5] = dL1;
  Curvature_Higgs_CT_L4[5][0][2][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][0][3][6] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][0][5][0] = dL1;
  Curvature_Higgs_CT_L4[5][0][6][3] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][0][7][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][1][1][5] = dL1;
  Curvature_Higgs_CT_L4[5][1][2][6] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][1][3][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][1][5][1] = dL1;
  Curvature_Higgs_CT_L4[5][1][6][2] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][1][7][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][2][0][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][2][1][6] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][2][2][5] = dL3;
  Curvature_Higgs_CT_L4[5][2][5][2] = dL3;
  Curvature_Higgs_CT_L4[5][2][6][1] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][2][7][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][3][0][6] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][3][1][7] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][3][3][5] = dL3;
  Curvature_Higgs_CT_L4[5][3][5][3] = dL3;
  Curvature_Higgs_CT_L4[5][3][6][0] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][3][7][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][4][4][5] = dL1;
  Curvature_Higgs_CT_L4[5][4][5][4] = dL1;
  Curvature_Higgs_CT_L4[5][4][6][7] = dL5;
  Curvature_Higgs_CT_L4[5][4][7][6] = dL5;
  Curvature_Higgs_CT_L4[5][5][0][0] = dL1;
  Curvature_Higgs_CT_L4[5][5][1][1] = dL1;
  Curvature_Higgs_CT_L4[5][5][2][2] = dL3;
  Curvature_Higgs_CT_L4[5][5][3][3] = dL3;
  Curvature_Higgs_CT_L4[5][5][4][4] = dL1;
  Curvature_Higgs_CT_L4[5][5][5][5] = 3*dL1;
  Curvature_Higgs_CT_L4[5][5][6][6] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[5][5][7][7] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[5][6][0][3] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][6][1][2] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][6][2][1] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][6][3][0] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][6][4][7] = dL5;
  Curvature_Higgs_CT_L4[5][6][5][6] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[5][6][6][5] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[5][6][7][4] = dL5;
  Curvature_Higgs_CT_L4[5][7][0][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][7][1][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][7][2][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][7][3][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[5][7][4][6] = dL5;
  Curvature_Higgs_CT_L4[5][7][5][7] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[5][7][6][4] = dL5;
  Curvature_Higgs_CT_L4[5][7][7][5] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[6][0][0][6] = dL3;
  Curvature_Higgs_CT_L4[6][0][2][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][0][3][5] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][0][4][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][0][5][3] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][0][6][0] = dL3;
  Curvature_Higgs_CT_L4[6][1][1][6] = dL3;
  Curvature_Higgs_CT_L4[6][1][2][5] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][1][3][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][1][4][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][1][5][2] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][1][6][1] = dL3;
  Curvature_Higgs_CT_L4[6][2][0][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][2][1][5] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][2][2][6] = dL2;
  Curvature_Higgs_CT_L4[6][2][4][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][2][5][1] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][2][6][2] = dL2;
  Curvature_Higgs_CT_L4[6][3][0][5] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][3][1][4] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][3][3][6] = dL2;
  Curvature_Higgs_CT_L4[6][3][4][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][3][5][0] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][3][6][3] = dL2;
  Curvature_Higgs_CT_L4[6][4][0][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][4][1][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][4][2][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][4][3][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][4][4][6] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[6][4][5][7] = dL5;
  Curvature_Higgs_CT_L4[6][4][6][4] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[6][4][7][5] = dL5;
  Curvature_Higgs_CT_L4[6][5][0][3] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][5][1][2] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][5][2][1] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][5][3][0] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[6][5][4][7] = dL5;
  Curvature_Higgs_CT_L4[6][5][5][6] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[6][5][6][5] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[6][5][7][4] = dL5;
  Curvature_Higgs_CT_L4[6][6][0][0] = dL3;
  Curvature_Higgs_CT_L4[6][6][1][1] = dL3;
  Curvature_Higgs_CT_L4[6][6][2][2] = dL2;
  Curvature_Higgs_CT_L4[6][6][3][3] = dL2;
  Curvature_Higgs_CT_L4[6][6][4][4] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[6][6][5][5] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[6][6][6][6] = 3*dL2;
  Curvature_Higgs_CT_L4[6][6][7][7] = dL2;
  Curvature_Higgs_CT_L4[6][7][4][5] = dL5;
  Curvature_Higgs_CT_L4[6][7][5][4] = dL5;
  Curvature_Higgs_CT_L4[6][7][6][7] = dL2;
  Curvature_Higgs_CT_L4[6][7][7][6] = dL2;
  Curvature_Higgs_CT_L4[7][0][0][7] = dL3;
  Curvature_Higgs_CT_L4[7][0][2][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][0][3][4] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][0][4][3] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][0][5][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][0][7][0] = dL3;
  Curvature_Higgs_CT_L4[7][1][1][7] = dL3;
  Curvature_Higgs_CT_L4[7][1][2][4] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][1][3][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][1][4][2] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][1][5][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][1][7][1] = dL3;
  Curvature_Higgs_CT_L4[7][2][0][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][2][1][4] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][2][2][7] = dL2;
  Curvature_Higgs_CT_L4[7][2][4][1] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][2][5][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][2][7][2] = dL2;
  Curvature_Higgs_CT_L4[7][3][0][4] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][3][1][5] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][3][3][7] = dL2;
  Curvature_Higgs_CT_L4[7][3][4][0] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][3][5][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][3][7][3] = dL2;
  Curvature_Higgs_CT_L4[7][4][0][3] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][4][1][2] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][4][2][1] = (-4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][4][3][0] = (4*dL4 - 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][4][4][7] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[7][4][5][6] = dL5;
  Curvature_Higgs_CT_L4[7][4][6][5] = dL5;
  Curvature_Higgs_CT_L4[7][4][7][4] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[7][5][0][2] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][5][1][3] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][5][2][0] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][5][3][1] = (4*dL4 + 4*dL5)/8.;
  Curvature_Higgs_CT_L4[7][5][4][6] = dL5;
  Curvature_Higgs_CT_L4[7][5][5][7] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[7][5][6][4] = dL5;
  Curvature_Higgs_CT_L4[7][5][7][5] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[7][6][4][5] = dL5;
  Curvature_Higgs_CT_L4[7][6][5][4] = dL5;
  Curvature_Higgs_CT_L4[7][6][6][7] = dL2;
  Curvature_Higgs_CT_L4[7][6][7][6] = dL2;
  Curvature_Higgs_CT_L4[7][7][0][0] = dL3;
  Curvature_Higgs_CT_L4[7][7][1][1] = dL3;
  Curvature_Higgs_CT_L4[7][7][2][2] = dL2;
  Curvature_Higgs_CT_L4[7][7][3][3] = dL2;
  Curvature_Higgs_CT_L4[7][7][4][4] = (8*dL3 + 8*dL4 - 8*dL5)/8.;
  Curvature_Higgs_CT_L4[7][7][5][5] = (8*dL3 + 8*dL4 + 8*dL5)/8.;
  Curvature_Higgs_CT_L4[7][7][6][6] = dL2;
  Curvature_Higgs_CT_L4[7][7][7][7] = 3*dL2;


}

/**
 * console output of all parameters
 */
void Class_Potential_R2HDMMathematica::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "Model = " << Model << "\n";

  ss << "\nThe input parameters are : \n";
  ss << "type = " << type << "\n";
  ss << "L1 = " << L1 << "\n";
  ss << "L2 = " << L2 << "\n";
  ss << "L3 = " << L3 << "\n";
  ss << "L4 = " << L4 << "\n";
  ss << "L5 = " << L5 << "\n";
  ss << "m12Sq = " << m12Sq << "\n";
  ss << "tanbeta = " << tanbeta << "\n";

  ss << "\nThe parameters are : \n";
  ss << "m11Sq = " << m11Sq << "\n";
  ss << "m22Sq = " << m22Sq << "\n";
  ss << "m12Sq = " << m12Sq << "\n";
  ss << "L1 = " << L1 << "\n";
  ss << "L2 = " << L2 << "\n";
  ss << "L3 = " << L3 << "\n";
  ss << "L4 = " << L4 << "\n";
  ss << "L5 = " << L5 << "\n";

  ss << "\nThe counterterm parameters are : \n";
  ss << "dm11Sq = " << dm11Sq << "\n";
  ss << "dm22Sq = " << dm22Sq << "\n";
  ss << "dm12Sq = " << dm12Sq << "\n";
  ss << "dL1 = " << dL1 << "\n";
  ss << "dL2 = " << dL2 << "\n";
  ss << "dL3 = " << dL3 << "\n";
  ss << "dL4 = " << dL4 << "\n";
  ss << "dL5 = " << dL5 << "\n";
  ss << "dT1 = " << dT1 << "\n";
  ss << "dT2 = " << dT2 << "\n";
  ss << "dT3 = " << dT3 << "\n";
  ss << "dT4 = " << dT4 << "\n";
  ss << "dT5 = " << dT5 << "\n";
  ss << "dT6 = " << dT6 << "\n";
  ss << "dT7 = " << dT7 << "\n";
  ss << "dT8 = " << dT8 << "\n";

  ss << "\nThe scale is given by mu = " << scale << " GeV \n";

  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Potential_R2HDMMathematica::calc_CT() const
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
  parCT.push_back((v1*(-3*HesseWeinberg(0,0) + HesseWeinberg(4,4)) + v2*(HesseWeinberg(0,2) + HesseWeinberg(4,6) - 2*HesseWeinberg(5,7)))/(2.*v1)); //dm11Sq;
  parCT.push_back((v1*(HesseWeinberg(0,2) + HesseWeinberg(4,6) - 2*HesseWeinberg(5,7)) + v2*(-3*HesseWeinberg(2,2) + HesseWeinberg(6,6)))/(2.*v2)); //dm22Sq;
  parCT.push_back(2*HesseWeinberg(0,2) - HesseWeinberg(5,7)); //dm12Sq;
  parCT.push_back((v1*HesseWeinberg(0,0) - v2*HesseWeinberg(0,2) - v1*HesseWeinberg(4,4) + v2*HesseWeinberg(5,7))/pow(v1,3)); //dL1;
  parCT.push_back((v1*(-HesseWeinberg(0,2) + HesseWeinberg(5,7)) + v2*(HesseWeinberg(2,2) - HesseWeinberg(6,6)))/pow(v2,3)); //dL2;
  parCT.push_back((-HesseWeinberg(4,6) + HesseWeinberg(5,7))/(v1*v2)); //dL3;
  parCT.push_back(0); //dL4;
  parCT.push_back((2*(HesseWeinberg(0,2) - HesseWeinberg(5,7)))/(v1*v2)); //dL5;
  parCT.push_back(-NablaWeinberg(0)); //dT1;
  parCT.push_back(-NablaWeinberg(1)); //dT2;
  parCT.push_back(-NablaWeinberg(2)); //dT3;
  parCT.push_back(-NablaWeinberg(3)); //dT4;
  parCT.push_back(v1*HesseWeinberg(0,0) + v2*HesseWeinberg(0,2) - NablaWeinberg(4)); //dT5;
  parCT.push_back(-NablaWeinberg(5)); //dT6;
  parCT.push_back(v1*HesseWeinberg(0,2) + v2*HesseWeinberg(2,2) - NablaWeinberg(6)); //dT7;
  parCT.push_back(-NablaWeinberg(7)); //dT8;

  return parCT;
}

// mass basis triple couplings
void Class_Potential_R2HDMMathematica::TripleHiggsCouplings()
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

void Class_Potential_R2HDMMathematica::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // assign the non-zero entries
  Curvature_Higgs_L2[0][0] = m11Sq;
  Curvature_Higgs_L2[0][2] = -m12Sq;
  Curvature_Higgs_L2[1][1] = m11Sq;
  Curvature_Higgs_L2[1][3] = -m12Sq;
  Curvature_Higgs_L2[2][0] = -m12Sq;
  Curvature_Higgs_L2[2][2] = m22Sq;
  Curvature_Higgs_L2[3][1] = -m12Sq;
  Curvature_Higgs_L2[3][3] = m22Sq;
  Curvature_Higgs_L2[4][4] = m11Sq;
  Curvature_Higgs_L2[4][6] = -m12Sq;
  Curvature_Higgs_L2[5][5] = m11Sq;
  Curvature_Higgs_L2[5][7] = -m12Sq;
  Curvature_Higgs_L2[6][4] = -m12Sq;
  Curvature_Higgs_L2[6][6] = m22Sq;
  Curvature_Higgs_L2[7][5] = -m12Sq;
  Curvature_Higgs_L2[7][7] = m22Sq;

  Curvature_Higgs_L4[0][0][0][0] = 3*L1;
  Curvature_Higgs_L4[0][0][1][1] = L1;
  Curvature_Higgs_L4[0][0][2][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][5][5] = L1;
  Curvature_Higgs_L4[0][0][6][6] = L3;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][1][0][1] = L1;
  Curvature_Higgs_L4[0][1][1][0] = L1;
  Curvature_Higgs_L4[0][1][2][3] = L5;
  Curvature_Higgs_L4[0][1][3][2] = L5;
  Curvature_Higgs_L4[0][2][0][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][2][1][3] = L5;
  Curvature_Higgs_L4[0][2][2][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][2][3][1] = L5;
  Curvature_Higgs_L4[0][2][4][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][2][5][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][2][6][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][2][7][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][3][0][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][3][1][2] = L5;
  Curvature_Higgs_L4[0][3][2][1] = L5;
  Curvature_Higgs_L4[0][3][3][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][3][4][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][3][5][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][3][6][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][3][7][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][4][0][4] = L1;
  Curvature_Higgs_L4[0][4][2][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][4][3][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][4][4][0] = L1;
  Curvature_Higgs_L4[0][4][6][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][4][7][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][5][0][5] = L1;
  Curvature_Higgs_L4[0][5][2][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][5][3][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][5][5][0] = L1;
  Curvature_Higgs_L4[0][5][6][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][5][7][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][6][0][6] = L3;
  Curvature_Higgs_L4[0][6][2][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][6][3][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][6][4][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][6][5][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][6][6][0] = L3;
  Curvature_Higgs_L4[0][7][0][7] = L3;
  Curvature_Higgs_L4[0][7][2][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][7][3][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][7][4][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][7][5][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][7][7][0] = L3;
  Curvature_Higgs_L4[1][0][0][1] = L1;
  Curvature_Higgs_L4[1][0][1][0] = L1;
  Curvature_Higgs_L4[1][0][2][3] = L5;
  Curvature_Higgs_L4[1][0][3][2] = L5;
  Curvature_Higgs_L4[1][1][0][0] = L1;
  Curvature_Higgs_L4[1][1][1][1] = 3*L1;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][1][3][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][1][4][4] = L1;
  Curvature_Higgs_L4[1][1][5][5] = L1;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L3;
  Curvature_Higgs_L4[1][2][0][3] = L5;
  Curvature_Higgs_L4[1][2][1][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][2][2][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][2][3][0] = L5;
  Curvature_Higgs_L4[1][2][4][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][2][5][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][2][6][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][2][7][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][3][0][2] = L5;
  Curvature_Higgs_L4[1][3][1][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][3][2][0] = L5;
  Curvature_Higgs_L4[1][3][3][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][3][4][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][3][5][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][3][6][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][3][7][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][4][1][4] = L1;
  Curvature_Higgs_L4[1][4][2][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][4][3][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][4][4][1] = L1;
  Curvature_Higgs_L4[1][4][6][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][4][7][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][5][1][5] = L1;
  Curvature_Higgs_L4[1][5][2][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][5][3][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][5][5][1] = L1;
  Curvature_Higgs_L4[1][5][6][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][5][7][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][6][1][6] = L3;
  Curvature_Higgs_L4[1][6][2][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][6][3][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][6][4][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][6][5][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][6][6][1] = L3;
  Curvature_Higgs_L4[1][7][1][7] = L3;
  Curvature_Higgs_L4[1][7][2][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][7][3][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][7][4][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][7][5][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][7][7][1] = L3;
  Curvature_Higgs_L4[2][0][0][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][0][1][3] = L5;
  Curvature_Higgs_L4[2][0][2][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][0][3][1] = L5;
  Curvature_Higgs_L4[2][0][4][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][0][5][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][0][6][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][0][7][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][1][0][3] = L5;
  Curvature_Higgs_L4[2][1][1][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][1][2][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][1][3][0] = L5;
  Curvature_Higgs_L4[2][1][4][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][1][5][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][1][6][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][1][7][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][2][0][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][2][1][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][2][2][2] = 3*L2;
  Curvature_Higgs_L4[2][2][3][3] = L2;
  Curvature_Higgs_L4[2][2][4][4] = L3;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][6][6] = L2;
  Curvature_Higgs_L4[2][2][7][7] = L2;
  Curvature_Higgs_L4[2][3][0][1] = L5;
  Curvature_Higgs_L4[2][3][1][0] = L5;
  Curvature_Higgs_L4[2][3][2][3] = L2;
  Curvature_Higgs_L4[2][3][3][2] = L2;
  Curvature_Higgs_L4[2][4][0][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][4][1][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][4][2][4] = L3;
  Curvature_Higgs_L4[2][4][4][2] = L3;
  Curvature_Higgs_L4[2][4][6][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][4][7][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][5][0][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][5][1][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][5][2][5] = L3;
  Curvature_Higgs_L4[2][5][5][2] = L3;
  Curvature_Higgs_L4[2][5][6][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][5][7][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][6][0][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][6][1][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][6][2][6] = L2;
  Curvature_Higgs_L4[2][6][4][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][6][5][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][6][6][2] = L2;
  Curvature_Higgs_L4[2][7][0][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][7][1][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][7][2][7] = L2;
  Curvature_Higgs_L4[2][7][4][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][7][5][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][7][7][2] = L2;
  Curvature_Higgs_L4[3][0][0][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][0][1][2] = L5;
  Curvature_Higgs_L4[3][0][2][1] = L5;
  Curvature_Higgs_L4[3][0][3][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][0][4][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][0][5][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][0][6][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][0][7][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][1][0][2] = L5;
  Curvature_Higgs_L4[3][1][1][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][1][2][0] = L5;
  Curvature_Higgs_L4[3][1][3][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][1][4][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][1][5][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][1][6][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][1][7][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][2][0][1] = L5;
  Curvature_Higgs_L4[3][2][1][0] = L5;
  Curvature_Higgs_L4[3][2][2][3] = L2;
  Curvature_Higgs_L4[3][2][3][2] = L2;
  Curvature_Higgs_L4[3][3][0][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][3][1][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][3][2][2] = L2;
  Curvature_Higgs_L4[3][3][3][3] = 3*L2;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][5][5] = L3;
  Curvature_Higgs_L4[3][3][6][6] = L2;
  Curvature_Higgs_L4[3][3][7][7] = L2;
  Curvature_Higgs_L4[3][4][0][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][4][1][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][4][3][4] = L3;
  Curvature_Higgs_L4[3][4][4][3] = L3;
  Curvature_Higgs_L4[3][4][6][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][4][7][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][5][0][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][5][1][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][5][3][5] = L3;
  Curvature_Higgs_L4[3][5][5][3] = L3;
  Curvature_Higgs_L4[3][5][6][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][5][7][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][6][0][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][6][1][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][6][3][6] = L2;
  Curvature_Higgs_L4[3][6][4][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][6][5][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][6][6][3] = L2;
  Curvature_Higgs_L4[3][7][0][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][7][1][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][7][3][7] = L2;
  Curvature_Higgs_L4[3][7][4][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][7][5][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][7][7][3] = L2;
  Curvature_Higgs_L4[4][0][0][4] = L1;
  Curvature_Higgs_L4[4][0][2][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][0][3][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][0][4][0] = L1;
  Curvature_Higgs_L4[4][0][6][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][0][7][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][1][1][4] = L1;
  Curvature_Higgs_L4[4][1][2][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][1][3][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][1][4][1] = L1;
  Curvature_Higgs_L4[4][1][6][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][1][7][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][2][0][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][2][1][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][2][2][4] = L3;
  Curvature_Higgs_L4[4][2][4][2] = L3;
  Curvature_Higgs_L4[4][2][6][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][2][7][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][3][0][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][3][1][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][3][3][4] = L3;
  Curvature_Higgs_L4[4][3][4][3] = L3;
  Curvature_Higgs_L4[4][3][6][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][3][7][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][4][0][0] = L1;
  Curvature_Higgs_L4[4][4][1][1] = L1;
  Curvature_Higgs_L4[4][4][2][2] = L3;
  Curvature_Higgs_L4[4][4][3][3] = L3;
  Curvature_Higgs_L4[4][4][4][4] = 3*L1;
  Curvature_Higgs_L4[4][4][5][5] = L1;
  Curvature_Higgs_L4[4][4][6][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][5][4][5] = L1;
  Curvature_Higgs_L4[4][5][5][4] = L1;
  Curvature_Higgs_L4[4][5][6][7] = L5;
  Curvature_Higgs_L4[4][5][7][6] = L5;
  Curvature_Higgs_L4[4][6][0][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][6][1][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][6][2][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][6][3][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][6][4][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][6][5][7] = L5;
  Curvature_Higgs_L4[4][6][6][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][6][7][5] = L5;
  Curvature_Higgs_L4[4][7][0][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][7][1][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][7][2][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][7][3][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][7][4][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][7][5][6] = L5;
  Curvature_Higgs_L4[4][7][6][5] = L5;
  Curvature_Higgs_L4[4][7][7][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][0][0][5] = L1;
  Curvature_Higgs_L4[5][0][2][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][0][3][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][0][5][0] = L1;
  Curvature_Higgs_L4[5][0][6][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][0][7][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][1][1][5] = L1;
  Curvature_Higgs_L4[5][1][2][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][1][3][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][1][5][1] = L1;
  Curvature_Higgs_L4[5][1][6][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][1][7][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][2][0][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][2][1][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][2][2][5] = L3;
  Curvature_Higgs_L4[5][2][5][2] = L3;
  Curvature_Higgs_L4[5][2][6][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][2][7][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][3][0][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][3][1][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][3][3][5] = L3;
  Curvature_Higgs_L4[5][3][5][3] = L3;
  Curvature_Higgs_L4[5][3][6][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][3][7][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][4][4][5] = L1;
  Curvature_Higgs_L4[5][4][5][4] = L1;
  Curvature_Higgs_L4[5][4][6][7] = L5;
  Curvature_Higgs_L4[5][4][7][6] = L5;
  Curvature_Higgs_L4[5][5][0][0] = L1;
  Curvature_Higgs_L4[5][5][1][1] = L1;
  Curvature_Higgs_L4[5][5][2][2] = L3;
  Curvature_Higgs_L4[5][5][3][3] = L3;
  Curvature_Higgs_L4[5][5][4][4] = L1;
  Curvature_Higgs_L4[5][5][5][5] = 3*L1;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][5][7][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][6][0][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][6][1][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][6][2][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][6][3][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][6][4][7] = L5;
  Curvature_Higgs_L4[5][6][5][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][6][6][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][6][7][4] = L5;
  Curvature_Higgs_L4[5][7][0][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][7][1][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][7][2][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][7][3][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][7][4][6] = L5;
  Curvature_Higgs_L4[5][7][5][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][7][6][4] = L5;
  Curvature_Higgs_L4[5][7][7][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][0][0][6] = L3;
  Curvature_Higgs_L4[6][0][2][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][0][3][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][0][4][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][0][5][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][0][6][0] = L3;
  Curvature_Higgs_L4[6][1][1][6] = L3;
  Curvature_Higgs_L4[6][1][2][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][1][3][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][1][4][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][1][5][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][1][6][1] = L3;
  Curvature_Higgs_L4[6][2][0][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][2][1][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][2][2][6] = L2;
  Curvature_Higgs_L4[6][2][4][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][2][5][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][2][6][2] = L2;
  Curvature_Higgs_L4[6][3][0][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][3][1][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][3][3][6] = L2;
  Curvature_Higgs_L4[6][3][4][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][3][5][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][3][6][3] = L2;
  Curvature_Higgs_L4[6][4][0][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][4][1][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][4][2][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][4][3][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][4][4][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][4][5][7] = L5;
  Curvature_Higgs_L4[6][4][6][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][4][7][5] = L5;
  Curvature_Higgs_L4[6][5][0][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][5][1][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][5][2][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][5][3][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][5][4][7] = L5;
  Curvature_Higgs_L4[6][5][5][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][5][6][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][5][7][4] = L5;
  Curvature_Higgs_L4[6][6][0][0] = L3;
  Curvature_Higgs_L4[6][6][1][1] = L3;
  Curvature_Higgs_L4[6][6][2][2] = L2;
  Curvature_Higgs_L4[6][6][3][3] = L2;
  Curvature_Higgs_L4[6][6][4][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][6][5][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][6][6][6] = 3*L2;
  Curvature_Higgs_L4[6][6][7][7] = L2;
  Curvature_Higgs_L4[6][7][4][5] = L5;
  Curvature_Higgs_L4[6][7][5][4] = L5;
  Curvature_Higgs_L4[6][7][6][7] = L2;
  Curvature_Higgs_L4[6][7][7][6] = L2;
  Curvature_Higgs_L4[7][0][0][7] = L3;
  Curvature_Higgs_L4[7][0][2][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][0][3][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][0][4][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][0][5][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][0][7][0] = L3;
  Curvature_Higgs_L4[7][1][1][7] = L3;
  Curvature_Higgs_L4[7][1][2][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][1][3][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][1][4][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][1][5][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][1][7][1] = L3;
  Curvature_Higgs_L4[7][2][0][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][2][1][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][2][2][7] = L2;
  Curvature_Higgs_L4[7][2][4][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][2][5][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][2][7][2] = L2;
  Curvature_Higgs_L4[7][3][0][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][3][1][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][3][3][7] = L2;
  Curvature_Higgs_L4[7][3][4][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][3][5][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][3][7][3] = L2;
  Curvature_Higgs_L4[7][4][0][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][4][1][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][4][2][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][4][3][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][4][4][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][4][5][6] = L5;
  Curvature_Higgs_L4[7][4][6][5] = L5;
  Curvature_Higgs_L4[7][4][7][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][5][0][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][5][1][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][5][2][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][5][3][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][5][4][6] = L5;
  Curvature_Higgs_L4[7][5][5][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][5][6][4] = L5;
  Curvature_Higgs_L4[7][5][7][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][6][4][5] = L5;
  Curvature_Higgs_L4[7][6][5][4] = L5;
  Curvature_Higgs_L4[7][6][6][7] = L2;
  Curvature_Higgs_L4[7][6][7][6] = L2;
  Curvature_Higgs_L4[7][7][0][0] = L3;
  Curvature_Higgs_L4[7][7][1][1] = L3;
  Curvature_Higgs_L4[7][7][2][2] = L2;
  Curvature_Higgs_L4[7][7][3][3] = L2;
  Curvature_Higgs_L4[7][7][4][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][7][5][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][7][6][6] = L2;
  Curvature_Higgs_L4[7][7][7][7] = 3*L2;

  Curvature_Gauge_G2H2[0][0][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][4][4] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][5][5] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][6][6] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][7][7] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][3][0][4] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][1][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][2][6] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][3][7] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][4][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][5][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][6][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][7][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][1][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][4][4] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][5][5] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][6][6] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][7][7] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][3][0][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][1][4] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][2][7] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][3][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][4][1] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][5][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][6][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][7][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][2][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][4][4] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][5][5] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][6][6] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][7][7] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][3][0][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][1][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][2][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][3][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][4][4] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][5][5] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][6][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][7][7] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][0][0][4] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][1][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][2][6] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][3][7] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][4][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][5][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][6][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][7][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][0][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][1][4] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][2][7] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][3][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][4][1] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][5][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][6][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][7][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][0][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][1][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][2][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][3][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][4][4] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][5][5] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][6][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][7][7] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][3][0][0] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][1][1] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][2][2] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][3][3] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][4][4] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][5][5] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][6][6] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][7][7] = pow(SMConstants.C_gs,2)/2.;

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

  Curvature_Lepton_F2H1[0][1][6] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[0][1][7] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[1][0][6] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[1][0][7] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[1][6][2] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[1][6][3] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[2][3][6] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[2][3][7] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[3][2][6] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[3][2][7] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[3][7][2] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[3][7][3] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[4][5][6] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[4][5][7] = (II*SMConstants.C_MassTau)/v2;
  Curvature_Lepton_F2H1[5][4][6] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[5][4][7] = (II*SMConstants.C_MassTau)/v2;
  Curvature_Lepton_F2H1[5][8][2] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[5][8][3] = (II*SMConstants.C_MassTau)/v2;
  Curvature_Lepton_F2H1[6][1][2] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[6][1][3] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[7][3][2] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[7][3][3] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[8][5][2] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[8][5][3] = (II*SMConstants.C_MassTau)/v2;

  Curvature_Quark_F2H1[0][6][6] = SMConstants.C_MassUp/v2;
  Curvature_Quark_F2H1[0][6][7] = (-II*SMConstants.C_MassUp)/v2;
  Curvature_Quark_F2H1[0][9][2] = -((SMConstants.C_MassUp*conj(V11))/v2);
  Curvature_Quark_F2H1[0][9][3] = (II*SMConstants.C_MassUp*conj(V11))/v2;
  Curvature_Quark_F2H1[0][10][2] = -((SMConstants.C_MassUp*conj(V12))/v2);
  Curvature_Quark_F2H1[0][10][3] = (II*SMConstants.C_MassUp*conj(V12))/v2;
  Curvature_Quark_F2H1[0][11][2] = -((SMConstants.C_MassUp*conj(V13))/v2);
  Curvature_Quark_F2H1[0][11][3] = (II*SMConstants.C_MassUp*conj(V13))/v2;
  Curvature_Quark_F2H1[1][7][6] = SMConstants.C_MassCharm/v2;
  Curvature_Quark_F2H1[1][7][7] = (-II*SMConstants.C_MassCharm)/v2;
  Curvature_Quark_F2H1[1][9][2] = -((SMConstants.C_MassCharm*conj(V21))/v2);
  Curvature_Quark_F2H1[1][9][3] = (II*SMConstants.C_MassCharm*conj(V21))/v2;
  Curvature_Quark_F2H1[1][10][2] = -((SMConstants.C_MassCharm*conj(V22))/v2);
  Curvature_Quark_F2H1[1][10][3] = (II*SMConstants.C_MassCharm*conj(V22))/v2;
  Curvature_Quark_F2H1[1][11][2] = -((SMConstants.C_MassCharm*conj(V23))/v2);
  Curvature_Quark_F2H1[1][11][3] = (II*SMConstants.C_MassCharm*conj(V23))/v2;
  Curvature_Quark_F2H1[2][8][6] = SMConstants.C_MassTop/v2;
  Curvature_Quark_F2H1[2][8][7] = (-II*SMConstants.C_MassTop)/v2;
  Curvature_Quark_F2H1[2][9][2] = -((SMConstants.C_MassTop*conj(V31))/v2);
  Curvature_Quark_F2H1[2][9][3] = (II*SMConstants.C_MassTop*conj(V31))/v2;
  Curvature_Quark_F2H1[2][10][2] = -((SMConstants.C_MassTop*conj(V32))/v2);
  Curvature_Quark_F2H1[2][10][3] = (II*SMConstants.C_MassTop*conj(V32))/v2;
  Curvature_Quark_F2H1[2][11][2] = -((SMConstants.C_MassTop*conj(V33))/v2);
  Curvature_Quark_F2H1[2][11][3] = (II*SMConstants.C_MassTop*conj(V33))/v2;
  Curvature_Quark_F2H1[3][6][2] = (SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[3][6][3] = (II*SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[3][7][2] = (SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[3][7][3] = (II*SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[3][8][2] = (SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[3][8][3] = (II*SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[3][9][6] = SMConstants.C_MassDown/v2;
  Curvature_Quark_F2H1[3][9][7] = (II*SMConstants.C_MassDown)/v2;
  Curvature_Quark_F2H1[4][6][2] = (SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[4][6][3] = (II*SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[4][7][2] = (SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[4][7][3] = (II*SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[4][8][2] = (SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[4][8][3] = (II*SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[4][10][6] = SMConstants.C_MassStrange/v2;
  Curvature_Quark_F2H1[4][10][7] = (II*SMConstants.C_MassStrange)/v2;
  Curvature_Quark_F2H1[5][6][2] = (SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[5][6][3] = (II*SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[5][7][2] = (SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[5][7][3] = (II*SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[5][8][2] = (SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[5][8][3] = (II*SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[5][11][6] = SMConstants.C_MassBottom/v2;
  Curvature_Quark_F2H1[5][11][7] = (II*SMConstants.C_MassBottom)/v2;
  Curvature_Quark_F2H1[6][0][6] = SMConstants.C_MassUp/v2;
  Curvature_Quark_F2H1[6][0][7] = (-II*SMConstants.C_MassUp)/v2;
  Curvature_Quark_F2H1[6][3][2] = (SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[6][3][3] = (II*SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[6][4][2] = (SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[6][4][3] = (II*SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[6][5][2] = (SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[6][5][3] = (II*SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[7][1][6] = SMConstants.C_MassCharm/v2;
  Curvature_Quark_F2H1[7][1][7] = (-II*SMConstants.C_MassCharm)/v2;
  Curvature_Quark_F2H1[7][3][2] = (SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[7][3][3] = (II*SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[7][4][2] = (SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[7][4][3] = (II*SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[7][5][2] = (SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[7][5][3] = (II*SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[8][2][6] = SMConstants.C_MassTop/v2;
  Curvature_Quark_F2H1[8][2][7] = (-II*SMConstants.C_MassTop)/v2;
  Curvature_Quark_F2H1[8][3][2] = (SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[8][3][3] = (II*SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[8][4][2] = (SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[8][4][3] = (II*SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[8][5][2] = (SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[8][5][3] = (II*SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[9][0][2] = -((SMConstants.C_MassUp*conj(V11))/v2);
  Curvature_Quark_F2H1[9][0][3] = (II*SMConstants.C_MassUp*conj(V11))/v2;
  Curvature_Quark_F2H1[9][1][2] = -((SMConstants.C_MassCharm*conj(V21))/v2);
  Curvature_Quark_F2H1[9][1][3] = (II*SMConstants.C_MassCharm*conj(V21))/v2;
  Curvature_Quark_F2H1[9][2][2] = -((SMConstants.C_MassTop*conj(V31))/v2);
  Curvature_Quark_F2H1[9][2][3] = (II*SMConstants.C_MassTop*conj(V31))/v2;
  Curvature_Quark_F2H1[9][3][6] = SMConstants.C_MassDown/v2;
  Curvature_Quark_F2H1[9][3][7] = (II*SMConstants.C_MassDown)/v2;
  Curvature_Quark_F2H1[10][0][2] = -((SMConstants.C_MassUp*conj(V12))/v2);
  Curvature_Quark_F2H1[10][0][3] = (II*SMConstants.C_MassUp*conj(V12))/v2;
  Curvature_Quark_F2H1[10][1][2] = -((SMConstants.C_MassCharm*conj(V22))/v2);
  Curvature_Quark_F2H1[10][1][3] = (II*SMConstants.C_MassCharm*conj(V22))/v2;
  Curvature_Quark_F2H1[10][2][2] = -((SMConstants.C_MassTop*conj(V32))/v2);
  Curvature_Quark_F2H1[10][2][3] = (II*SMConstants.C_MassTop*conj(V32))/v2;
  Curvature_Quark_F2H1[10][4][6] = SMConstants.C_MassStrange/v2;
  Curvature_Quark_F2H1[10][4][7] = (II*SMConstants.C_MassStrange)/v2;
  Curvature_Quark_F2H1[11][0][2] = -((SMConstants.C_MassUp*conj(V13))/v2);
  Curvature_Quark_F2H1[11][0][3] = (II*SMConstants.C_MassUp*conj(V13))/v2;
  Curvature_Quark_F2H1[11][1][2] = -((SMConstants.C_MassCharm*conj(V23))/v2);
  Curvature_Quark_F2H1[11][1][3] = (II*SMConstants.C_MassCharm*conj(V23))/v2;
  Curvature_Quark_F2H1[11][2][2] = -((SMConstants.C_MassTop*conj(V33))/v2);
  Curvature_Quark_F2H1[11][2][3] = (II*SMConstants.C_MassTop*conj(V33))/v2;
  Curvature_Quark_F2H1[11][5][6] = SMConstants.C_MassBottom/v2;
  Curvature_Quark_F2H1[11][5][7] = (II*SMConstants.C_MassBottom)/v2;

}

bool Class_Potential_R2HDMMathematica::CalculateDebyeSimplified()
{
  return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */
}

bool Class_Potential_R2HDMMathematica::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */
  return false;
}
double
Class_Potential_R2HDMMathematica::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (pow(v[4],4)*L1 + 2*pow(v[4],2)*(pow(v[2],2)*L3 + pow(v[7],2)*L3 + pow(v[7],2)*L4 - pow(v[7],2)*L5 + pow(v[6],2)*(L3 + L4 + L5) + 2*m11Sq) - 8*v[4]*v[6]*m12Sq + (pow(v[2],2) + pow(v[6],2) + pow(v[7],2))*(pow(v[2],2)*L2 + pow(v[6],2)*L2 + pow(v[7],2)*L2 + 4*m22Sq))/8.;
  
  return res;
}

double Class_Potential_R2HDMMathematica::VCounterSimplified(
    const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (pow(v[4],4)*dL1 + 2*pow(v[4],2)*(pow(v[2],2)*dL3 + pow(v[7],2)*dL3 + pow(v[7],2)*dL4 - pow(v[7],2)*dL5 + pow(v[6],2)*(dL3 + dL4 + dL5) + 2*dm11Sq) - 8*v[4]*v[6]*dm12Sq + (pow(v[2],2) + pow(v[6],2) + pow(v[7],2))*(pow(v[2],2)*dL2 + pow(v[6],2)*dL2 + pow(v[7],2)*dL2 + 4*dm22Sq))/8. + v[2]*dT3 + v[4]*dT5 + v[6]*dT7 + v[7]*dT8;
  
  return res;
}

void Class_Potential_R2HDMMathematica::Debugging(const std::vector<double> &input,
                                            std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
