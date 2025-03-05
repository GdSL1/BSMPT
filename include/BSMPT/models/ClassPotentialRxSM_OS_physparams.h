#ifndef SRC_CLASSPOTENTIALRXSM_OS_PHYSPARAMS_H_
#define SRC_CLASSPOTENTIALRXSM_OS_PHYSPARAMS_H_

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{
class Class_Potential_RxSM_OS_physparams : public Class_Potential_Origin
{
public:
  Class_Potential_RxSM_OS_physparams(const ISMConstants &smConstants);
  virtual ~Class_Potential_RxSM_OS_physparams();

  // Initialize input parameters
  double m1 = 0;
  double m2 = 0;
  double a = 0;
  double vh = 0;
  double vs = 0;

  // Initialize dependent parameters
  double mh = 0;
  double ms = 0;
  double Lh = 0;
  double Ls = 0;
  double Lhs = 0;

  // Initialize counter terms
  double dmh = 0;
  double dms = 0;
  double dLh = 0;
  double dLs = 0;
  double dLhs = 0;
  double dT1 = 0;
  double dT2 = 0;
  double dT3 = 0;
  double dT4 = 0;
  double dT5 = 0;


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
#endif /* SRC_RXSM_OS_PHYSPARAMS_H_ */
