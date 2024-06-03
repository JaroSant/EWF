
#ifndef __EA__WrightFisher__
#define __EA__WrightFisher__

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace boost::multiprecision;

#include "Polynomial.h"
#include "myHelpers.h"

/// typedef boost::multiprecision::cpp_dec_float_100 double100; /// Can switch
/// to this for increased precision, but makes computations much slower!!
typedef double double100;

class WrightFisher {
 public:
  WrightFisher(vector<double> thetaP_in, bool non_neut_in, double100 sigma_in,
               int selectionsetup_in, double dom_in, int SelPolyDeg_in,
               vector<double> selCoefs_in)
      : thetaP(thetaP_in),
        non_neutral(non_neut_in),
        sigma(sigma_in),
        SelectionSetup(selectionsetup_in),
        dominanceParameter(dom_in),
        SelPolyDeg(SelPolyDeg_in),
        selectionCoeffs(selCoefs_in) {
    ThetaSetter();
    SelectionSetter();
    PhiSetter();
    std::cout << R"(
                           
     +++           )|(     
    (o o)         (o o)    
ooO--(_)--Ooo-ooO--(_)--Ooo
    _______       ________
   / ____/ |     / / ____/
  / __/  | | /| / / /_    
 / /___  | |/ |/ / __/    
/_____/  |__/|__/_/       
                          
)" << '\n';
    std::cout << "You've instantiated a WrightFisher class with the following "
                 "parameters:"
              << std::endl;
    if (!thetaP.empty()) {
      std::cout << "Theta vector: (" << thetaP.front() << ", " << thetaP.back()
                << ")" << std::endl;
    } else {
      std::cout << "Theta vector: (0.0, 0.0)" << std::endl;
    }
    if (non_neutral) {
      if (SelectionSetup == 0) {
        std::cout << "Genic selection: " << sigma << std::endl;
      } else if (SelectionSetup == 1) {
        std::cout << "Diploid selection with:" << std::endl;
        std::cout << "Sigma: " << sigma << std::endl;
        std::cout << "Dominance parameter: " << dominanceParameter << std::endl;
      } else {
        std::cout << "Polynomial selection with degree " << SelPolyDeg
                  << " with entries:" << std::endl;
        for (vector<double>::iterator sc_it = selectionCoeffs.begin();
             sc_it != selectionCoeffs.end(); sc_it++) {
          std::cout << *sc_it << ", ";
        }
        std::cout << "." << std::endl;
      }
    } else {
      std::cout << "No selection." << std::endl;
    }
  }

  double100 phiMin, phiMax, AtildeMax;

  /// HELPER FUNCTIONS
  void ThetaSetter();
  void ThetaResetter();
  void SelectionSetter();
  void PhiSetter();
  vector<double100> get_Theta();
  double100 Phitilde(double100 y);
  vector<double100> PhitildeMinMaxRange();
  double100 Atilde(double100 x);
  double100 Atildeplus();
  pair<double, double> GriffithsParas(double100 t);
  double100 computeLogBeta(int m, int k);
  double100 NormCDF(double100 x, double100 m, double100 v);
  double100 DiscretisedNormCDF(int m, double100 t);
  double100 LogBinomialCoefficientCalculator(int n, int k);
  double100 UnconditionedDiffusionDensity(double100 x, double100 y, double100 t,
                                          const Options &o);
  double100 DiffusionDensityApproximationDenom(double100 x, double100 t,
                                               const Options &o);
  double100 DiffusionDensityApproximation(double100 x, double100 y, double100 t,
                                          const Options &o);
  double100 QmApprox(int m, double100 t, const Options &o);
  double100 UnconditionedBridgeDensity(double100 x, double100 z, double100 y,
                                       double100 s, double100 t,
                                       const Options &o);
  double100 BridgeDenominatorApprox(double100 x, double100 z, double100 t,
                                    const Options &o);
  double100 DiffusionBridgeDensityApproximation(double100 x, double100 z,
                                                double100 y, double100 t,
                                                double100 s, const Options &o);
  double100 BridgeDenom(double100 x, double100 z, double100 y, double100 s,
                        double100 t, const Options &o);
  double100 ComputeDensity1(double100 x, double100 z, double100 y, double100 s,
                            double100 t, const Options &o);
  double100 ComputeDensity2(double100 x, double100 z, double100 y, double100 s,
                            double100 t, const Options &o);
  double100 ComputeDensity3(double100 x, double100 z, double100 y, double100 s,
                            double100 t, const Options &o);
  double100 ComputeDensity4(double100 x, double100 z, double100 y, double100 s,
                            double100 t, const Options &o);
  double100 BridgeDensity(double100 x, double100 z, double100 y, double100 s,
                          double100 t, const Options &o);
  double100 LogSumExp(vector<double100> &vecProb, double100 maxProb);
  double100 CustomGammaRatio(double100 a, double100 b);

  /// DIFFUSION SIMULATION - NEUTRAL PATHS

  pair<int, int> DrawAncestralProcess(double100 t, const Options &o,
                                      boost::random::mt19937 &gen);
  pair<int, int> DrawAncestralProcessConditionalZero(
      double100 t, const Options &o, boost::random::mt19937 &gen);
  pair<pair<int, int>, int> DrawAncestralProcessConditionalInterior(
      double100 t, double100 x, const Options &o, boost::random::mt19937 &gen);
  int DrawAncestralProcessG1984(double100 t, boost::random::mt19937 &gen);
  int DrawSizebiasedAncestralProcess(int d, double100 t,
                                     boost::random::mt19937 &gen);
  pair<double100, int> DrawEndpoint(double100 x, double100 t1, double100 t2,
                                    const Options &o,
                                    boost::random::mt19937 &gen);
  pair<double100, int> DrawUnconditionedDiffusion(double100 x, double100 t,
                                                  const Options &o,
                                                  boost::random::mt19937 &gen);

  /// DIFFUSION SIMULATION - NON-NEUTRAL PATHS

  vector<vector<double100>> NonNeutralDraw(double100 x, double100 t1,
                                           double100 t2, bool Absorption,
                                           const Options &o,
                                           boost::random::mt19937 &gen);
  pair<double100, int> NonNeutralDrawEndpoint(double100 x, double100 t1,
                                              double100 t2, bool Absorption,
                                              const Options &o,
                                              boost::random::mt19937 &gen);

  /// BRIDGE SIMULATION - NEUTRAL PATHS

  vector<int> DrawBridgePMFUnconditional(double100 x, double100 z, double100 s,
                                         double100 t, const Options &o,
                                         boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFUnconditionalOneQApprox(double100 x, double100 z,
                                                   double100 s, double100 t,
                                                   const Options &o,
                                                   boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFUnconditionalApprox(double100 x, double100 z,
                                               double100 s, double100 t,
                                               const Options &o,
                                               boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFSameMutation(double100 x, double100 s, double100 t,
                                        const Options &o,
                                        boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFSameMutationOneQApprox(double100 x, double100 s,
                                                  double100 t, const Options &o,
                                                  boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFSameMutationApprox(double100 x, double100 s,
                                              double100 t,
                                              boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDifferentMutation(double100 s, double100 t,
                                             double100 x, const Options &o,
                                             boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDifferentMutationOneQApprox(
      double100 s, double100 t, double100 x, const Options &o,
      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDifferentMutationApprox(double100 s, double100 t,
                                                   double100 x,
                                                   const Options &o,
                                                   boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFInteriorMutation(double100 x, double100 z,
                                            double100 s, double100 t,
                                            const Options &o,
                                            boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFInteriorMutationOneQApprox(
      double100 x, double100 z, double100 s, double100 t, const Options &o,
      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFInteriorMutationApprox(double100 x, double100 z,
                                                  double100 s, double100 t,
                                                  const Options &o,
                                                  boost::random::mt19937 &gen);
  vector<int> DrawBridgePMF(double100 x, double100 z, double100 s, double100 t,
                            const Options &o, boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFOneQApprox(double100 x, double100 z, double100 s,
                                      double100 t, const Options &o,
                                      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFG1984(double100 x, double100 z, double100 s,
                                 double100 t, const Options &o,
                                 boost::random::mt19937 &gen);
  double100 mkModeFinder_Evaluator(int m, int k, double100 x, double100 z,
                                   double100 s, double100 t, const Options &o);
  double100 mkjModeFinder_Evaluator(int m, int k, int j, double100 x,
                                    double100 z, double100 s, double100 t,
                                    const Options &o);
  double100 mkljModeFinder_Evaluator(int m, int k, int l, int j, double100 x,
                                     double100 z, double100 s, double100 t,
                                     const Options &o);
  double100 mklModeFinder_Evaluator(int m, int k, int l, double100 x,
                                    double100 z, double100 s, double100 t,
                                    const Options &o);
  vector<int> mkModeFinder(double100 x, double100 z, double100 s, double100 t,
                           const Options &o);
  vector<int> mkjModeFinder(double100 x, double100 z, double100 s, double100 t,
                            const Options &o);
  vector<int> mkljModeFinder(double100 x, double100 z, double100 s, double100 t,
                             const Options &o);
  vector<int> mklModeFinder(double100 x, double100 z, double100 s, double100 t,
                            const Options &o);
  pair<double100, int> DrawBridgepoint(double100 x, double100 z, double100 t1,
                                       double100 t2, double100 s,
                                       const Options &o,
                                       boost::random::mt19937 &gen);
  pair<double100, int> DrawUnconditionedBridge(double100 x, double100 z,
                                               double100 t1, double100 t2,
                                               double100 s, const Options &o,
                                               boost::random::mt19937 &gen);

  /// BRIDGE SIMULATION - NON-NEUTRAL PATHS

  vector<vector<double100>> NonNeutralDrawBridge(double100 x, double100 t1,
                                                 double100 t2, double100 z,
                                                 bool Absorption,
                                                 const Options &o,
                                                 boost::random::mt19937 &gen);
  pair<double100, int> NonNeutralDrawBridgepoint(
      double100 x, double100 t1, double100 t2, double100 z, double100 testT,
      bool Absorption, const Options &o, boost::random::mt19937 &gen);

  /// SIMULATION RUNNER FUNCTIONS

  void DiffusionRunner(int nSim, double100 x, double100 startT, double100 endT,
                       bool Absorption, string &Filename,
                       double100 diffusion_threshold,
                       double100 bridge_threshold);
  void BridgeDiffusionRunner(int nSim, double100 x, double100 z,
                             double100 startT, double100 endT,
                             double100 sampleT, bool Absorption,
                             string &Filename, double100 diffusion_threshold,
                             double100 bridge_threshold);
  void DiffusionDensityCalculator(int meshSize, double100 x, double100 startT,
                                  double100 endT, bool Absorption,
                                  string &Filename,
                                  double100 diffusion_threshold,
                                  double100 bridge_threshold);
  void BridgeDiffusionDensityCalculator(int meshSize, double100 x, double100 z,
                                        double100 startT, double100 endT,
                                        double100 sampleT, bool Absorption,
                                        string &Filename,
                                        double100 diffusion_threshold,
                                        double100 bridge_threshold);

 private:
  /// WRIGHT-FISHER PROPERTIES

  vector<double> thetaP;
  bool non_neutral;
  double theta, sigma;
  int SelectionSetup;
  double dominanceParameter;
  int SelPolyDeg;
  vector<double> selectionCoeffs;
  Polynomial SelectionFunction, PhiFunction, AtildeFunction;
  int thetaIndex;
  vector<vector<double100>> akm;

  /// HELPER FUNCTIONS

  template <typename T>
  T Getlogakm(int k, int m);
  int radiate_from_mode(int index, const double100 t) const;
  void increment_on_mk(vector<int> &mk, const double100 s,
                       const double100 t) const;
  double100 Getd(vector<double100> &d, int i, double100 x, double100 z,
                 double100 t);
  double100 Getd2(vector<double100> &d, int i, double100 x, double100 t);
  double100 GetdBridgeSame(vector<double100> &d, int i, double100 x,
                           double100 t);
  double100 GetdBridgeInterior(vector<double100> &d, int i, double100 x,
                               double100 z, double100 t);
  double100 GetdBridgeUnconditional(vector<double100> &d, int i, double100 x,
                                    double100 z, double100 t);
  double100 computeA(int m, int k, int l, int j, double100 x, double100 z);
  double100 computeAUnconditional(int m, int k, int l, double100 x,
                                  double100 z);
  int computeC(int m, pair<vector<int>, double100> &C);
  int computeE(pair<vector<int>, double100> &C);
};

template <typename T>
T WrightFisher::Getlogakm(int k, int m) {
  assert(k >= m);
  if (m > static_cast<int>(akm.size()) - 1) akm.resize(m + 1);
  if ((k - m > static_cast<int>(akm[m].size()) - 1)) {
    int oldsize = static_cast<int>(akm[m].size());
    akm[m].resize(k - m + 1, 0);

    for (int i = oldsize + m; i <= k; ++i) {
      if (i == 0)
        akm[m][i - m] = 0.0;
      else {
        T a = log(theta + 2.0 * i - 1);
        for (int j = 2; j <= i; ++j) {
          a += log(theta + m + j - 2.0);
          if (j <= i - m)
            a -= log(static_cast<T>(j));  /// We need (n-m)! in the denominator
          if (j <= m)
            a -= log(static_cast<T>(j));  /// We need m! in the denominator
        }
        akm[m][i - m] = a;  /// So akm[m] = (a{m,m}, a{m+1,m} a{m+2,m}, ...)
      }
    }
  }
  return akm[m][k - m];
}

#endif
