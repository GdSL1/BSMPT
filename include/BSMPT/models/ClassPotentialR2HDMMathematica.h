#ifndef SRC_CLASSPOTENTIALR2HDMMATHEMATICA_H_
#define SRC_CLASSPOTENTIALR2HDMMATHEMATICA_H_

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{
class Class_Potential_R2HDMMathematica : public Class_Potential_Origin
{
public:
  Class_Potential_R2HDMMathematica(const ISMConstants &smConstants);
  virtual ~Class_Potential_R2HDMMathematica();

  // Initialize input parameters
  double type = 0;
  double L1 = 0;
  double L2 = 0;
  double L3 = 0;
  double L4 = 0;
  double L5 = 0;
  double m12Sq = 0;
  double tanbeta = 0;

  // Initialize dependent parameters
  double m11Sq = 0;
  double m22Sq = 0;
  double v1 = 0;
  double v2 = 0;

  // Initialize counter terms
  double dm11Sq = 0;
  double dm22Sq = 0;
  double dm12Sq = 0;
  double dL1 = 0;
  double dL2 = 0;
  double dL3 = 0;
  double dL4 = 0;
  double dL5 = 0;
  double dT1 = 0;
  double dT2 = 0;
  double dT3 = 0;
  double dT4 = 0;
  double dT5 = 0;
  double dT6 = 0;
  double dT7 = 0;
  double dT8 = 0;


  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
#endif /* SRC_R2HDMMATHEMATICA_H_ */
