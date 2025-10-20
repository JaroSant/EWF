
#ifndef __EA__WrightFisher__
#define __EA__WrightFisher__

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

using namespace boost::multiprecision;

#include "Polynomial.h"
#include "myHelpers.h"

class WrightFisher {
 public:
  WrightFisher(vector<double> changepoints, vector<vector<double>> thetaP_in,
               bool non_neut_in, vector<double> sigma_in, int selectionsetup_in,
               double dom_in, int SelPolyDeg_in, vector<double> selCoefs_in)
      : changepts(std::move(changepoints)),
        thetaP(std::move(thetaP_in)),
        non_neutral(std::move(non_neut_in)),
        sigma(std::move(sigma_in)),
        SelectionSetup(std::move(selectionsetup_in)),
        dominanceParameter(std::move(dom_in)),
        SelPolyDeg(std::move(SelPolyDeg_in)),
        selectionCoeffs(std::move(selCoefs_in)) {
    size_t N = thetaP.size();
    theta.assign(N, 0.0);
    SelectionFunction.assign(N, Polynomial());
    PhiFunction.assign(N, Polynomial());
    AtildeFunction.assign(N, Polynomial());

    for (size_t N_i = 0; N_i < N; ++N_i) {
      ThetaSetter(N_i);
      SelectionSetter(N_i);
      // PhiSetter(N_i);
    }
  }

  vector<double> phiMin, phiMax, AtildeMax, AtildeMin;

  /// HELPER FUNCTIONS
  void ThetaSetter(size_t N_i);
  void ThetaResetter(size_t N_i);
  void SelectionSetter(size_t N_i);
  void PhiSetter(size_t N_i);
  vector<vector<double>> get_Theta();
  double Phitilde(size_t N_i, double y);
  vector<double> PhitildeMinMaxRange(size_t N_i);
  double Atilde(size_t N_i, double x);
  double Atildeplus(size_t N_i);
  double Atildeminus(size_t N_i);
  size_t first_greater_index(double s);
  size_t first_ge_index(double s);
  size_t getIndex(double s);
  pair<double, double> GriffithsParas(size_t N_i, double t);
  double computeLogBeta(int m, int k);
  double NormCDF(double x, double m, double v);
  double DiscretisedNormCDF(size_t N_i, int m, double t);
  double LogBinomialCoefficientCalculator(int n, int k);
  double UnconditionedDiffusionDensity(size_t N_i, double x, double y, double t,
                                       const Options &o);
  double DiffusionDensityApproximationDenom(size_t N_i, double x, double t,
                                            const Options &o);
  double DiffusionDensityApproximation(size_t N_i, double x, double y, double t,
                                       const Options &o);
  double QmApprox(size_t N_i, int m, double t, const Options &o);
  double UnconditionedBridgeDensity(size_t N_i, double x, double z, double y,
                                    double s, double t, const Options &o);
  double BridgeDenom(size_t N_i, double x, double z, double y, double s,
                     double t, const Options &o);
  double ComputeDensity1(size_t N_i, double x, double z, double y, double s,
                         double t, const Options &o);
  double ComputeDensity2(size_t N_i, double x, double z, double y, double s,
                         double t, const Options &o);
  double ComputeDensity3(size_t N_i, double x, double z, double y, double s,
                         double t, const Options &o);
  double ComputeDensity4(size_t N_i, double x, double z, double y, double s,
                         double t, const Options &o);
  double ComputeDensityDiffTheta(size_t N_i, double x, double z, double y,
                                 double s, double t, const Options &o);
  double ComputeDensityDiffThetaBoundaries(size_t N_i, double x, double z,
                                           double y, double s, double t,
                                           const Options &o);
  double ComputeDensityDiffThetaInterior(size_t N_i, double x, double z,
                                         double y, double s, double t,
                                         const Options &o);
  double BridgeDensity(size_t N_i, double x, double z, double y, double s,
                       double t, const Options &o);
  double BridgeDiffThetaDensity(size_t N_i, double x, double z, double y,
                                double s, double t, const Options &o);
  double LogSumExp(vector<double> &vecProb, double maxProb);

  /// DIFFUSION SIMULATION - NEUTRAL PATHS

  pair<int, int> DrawAncestralProcess(size_t N_i, double t, const Options &o,
                                      boost::random::mt19937 &gen);
  pair<int, int> DrawAncestralProcessConditionalZero(
      size_t N_i, double t, const Options &o, boost::random::mt19937 &gen);
  pair<pair<int, int>, int> DrawAncestralProcessConditionalInterior(
      size_t N_i, double t, double x, const Options &o,
      boost::random::mt19937 &gen);
  int DrawAncestralProcessG1984(size_t N_i, double t,
                                boost::random::mt19937 &gen);
  int DrawSizebiasedAncestralProcess(size_t N_i, int d, double t,
                                     boost::random::mt19937 &gen);
  pair<double, int> DrawEndpoint(size_t N_i, double x, double t1, double t2,
                                 const Options &o, boost::random::mt19937 &gen);
  pair<double, int> DrawUnconditionedDiffusion(size_t N_i, double x, double t,
                                               const Options &o,
                                               boost::random::mt19937 &gen);

  /// DIFFUSION SIMULATION - NON-NEUTRAL PATHS

  vector<vector<double>> NonNeutralDraw(size_t N_i, double x, double t1,
                                        double t2, bool Absorption,
                                        const Options &o,
                                        boost::random::mt19937 &gen);
  pair<double, int> NonNeutralDrawEndpoint(size_t N_i, double x, double t1,
                                           double t2, bool Absorption,
                                           const Options &o,
                                           boost::random::mt19937 &gen);

  /// BRIDGE SIMULATION - NEUTRAL PATHS

  vector<int> DrawBridgePMFUnconditional(size_t N_i, double x, double z,
                                         double s, double t, const Options &o,
                                         boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFUnconditionalOneQApprox(size_t N_i, double x,
                                                   double z, double s, double t,
                                                   const Options &o,
                                                   boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFUnconditionalApprox(size_t N_i, double x, double z,
                                               double s, double t,
                                               const Options &o,
                                               boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFSameMutation(size_t N_i, double x, double s,
                                        double t, const Options &o,
                                        boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFSameMutationOneQApprox(size_t N_i, double x,
                                                  double s, double t,
                                                  const Options &o,
                                                  boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFSameMutationApprox(size_t N_i, double x, double s,
                                              double t,
                                              boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDifferentMutation(size_t N_i, double s, double t,
                                             double x, const Options &o,
                                             boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDifferentMutationOneQApprox(
      size_t N_i, double s, double t, double x, const Options &o,
      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDifferentMutationApprox(size_t N_i, double s,
                                                   double t, double x,
                                                   const Options &o,
                                                   boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFInteriorMutation(size_t N_i, double x, double z,
                                            double s, double t,
                                            const Options &o,
                                            boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFInteriorMutationOneQApprox(
      size_t N_i, double x, double z, double s, double t, const Options &o,
      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFInteriorMutationApprox(size_t N_i, double x,
                                                  double z, double s, double t,
                                                  const Options &o,
                                                  boost::random::mt19937 &gen);
  vector<int> DrawBridgePMF(size_t N_i, double x, double z, double s, double t,
                            const Options &o, boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFOneQApprox(size_t N_i, double x, double z, double s,
                                      double t, const Options &o,
                                      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFG1984(size_t N_i, double x, double z, double s,
                                 double t, const Options &o,
                                 boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffTheta(size_t N_i, double x, double z, double s,
                                     double t, const Options &o,
                                     boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaOneQApprox(size_t N_i, double x, double z,
                                               double s, double t,
                                               const Options &o,
                                               boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaApprox(size_t N_i, double x, double z,
                                           double s, double t, const Options &o,
                                           boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaBoundaries(size_t N_i, double x, double z,
                                               double s, double t,
                                               const Options &o,
                                               boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaBoundariesOneQApprox(
      size_t N_i, double x, double z, double s, double t, const Options &o,
      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaBoundariesApprox(
      size_t N_i, double x, double z, double s, double t, const Options &o,
      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaInterior(size_t N_i, double x, double z,
                                             double s, double t,
                                             const Options &o,
                                             boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaInteriorOneQApprox(
      size_t N_i, double x, double z, double s, double t, const Options &o,
      boost::random::mt19937 &gen);
  vector<int> DrawBridgePMFDiffThetaInteriorApprox(size_t N_i, double x,
                                                   double z, double s, double t,
                                                   const Options &o,
                                                   boost::random::mt19937 &gen);
  double mkModeFinder_Evaluator(bool isDiffTheta, size_t N_i, int m, int k,
                                double x, double z, double s, double t,
                                const Options &o);
  double mkjModeFinder_Evaluator(bool isDiffTheta, size_t N_i, int m, int k,
                                 int j, double x, double z, double s, double t,
                                 const Options &o);
  double mkjDensityModeFinder_Evaluator(bool isDiffTheta, size_t N_i, int m,
                                        int k, int j, double x, double z,
                                        double y, double s, double t,
                                        const Options &o);
  double mkljModeFinder_Evaluator(bool isDiffTheta, size_t N_i, int m, int k,
                                  int l, int j, double x, double z, double s,
                                  double t, const Options &o);
  double mklModeFinder_Evaluator(bool isDiffTheta, size_t N_i, int m, int k,
                                 int l, double x, double z, double s, double t,
                                 const Options &o);
  vector<int> mkModeFinder(bool isDiffTheta, size_t N_i, double x, double z,
                           double s, double t, const Options &o);
  vector<int> mkjModeFinder(bool isDiffTheta, size_t N_i, double x, double z,
                            double s, double t, const Options &o);
  vector<int> mkjDensityModeFinder(bool isDiffTheta, size_t N_i, double x,
                                   double z, double y, double s, double t,
                                   const Options &o);
  vector<int> mkljModeFinder(bool isDiffTheta, size_t N_i, double x, double z,
                             double s, double t, const Options &o);
  vector<int> mklModeFinder(bool isDiffTheta, size_t N_i, double x, double z,
                            double s, double t, const Options &o);
  pair<double, int> DrawBridgepoint(size_t N_i, double x, double z, double t1,
                                    double t2, double s, const Options &o,
                                    boost::random::mt19937 &gen);
  pair<double, int> DrawBridgepointDiffTheta(size_t N_i, double x, double z,
                                             double t1, double t2, double s,
                                             const Options &o,
                                             boost::random::mt19937 &gen);
  pair<double, int> DrawUnconditionedBridge(size_t N_i, double x, double z,
                                            double t1, double t2, double s,
                                            const Options &o,
                                            boost::random::mt19937 &gen);

  /// BRIDGE SIMULATION - NON-NEUTRAL PATHS

  vector<vector<double>> NonNeutralDrawBridge(size_t N_i, double x, double t1,
                                              double t2, double z,
                                              bool Absorption, const Options &o,
                                              boost::random::mt19937 &gen);
  pair<double, vector<vector<double>>> NonNeutralDrawBridgeDiffTheta(
      size_t N_i, double x, double t1, double t2, double z, double s,
      const Options &o, boost::random::mt19937 &gen);
  pair<double, int> NonNeutralDrawBridgepoint(size_t N_i, double x, double t1,
                                              double t2, double z, double testT,
                                              bool Absorption, const Options &o,
                                              boost::random::mt19937 &gen);
  pair<double, int> NonNeutralDrawBridgepointDiffTheta(
      size_t N_i, double x, double t1, double t2, double z, double testT,
      const Options &o, boost::random::mt19937 &gen);

  /// SIMULATION RUNNER FUNCTIONS

  void DiffusionRunner(int nSim, double x, double startT, double endT,
                       bool Absorption, string &Filename,
                       double diffusion_threshold, double bridge_threshold);
  void DiffusionRunnerVector(int nSim, vector<double> x, double startT,
                             double endT, bool Absorption, string &Filename,
                             double diffusion_threshold,
                             double bridge_threshold);
  void DiffusionTrajectoryVector(int nSim, double x, vector<double> times,
                                 bool Absorption, string &Filename,
                                 double diffusion_threshold,
                                 double bridge_threshold);
  void BridgeDiffusionRunner(int nSim, double x, double z, double startT,
                             double endT, double sampleT, bool Absorption,
                             string &Filename, bool verbose,
                             double diffusion_threshold,
                             double bridge_threshold);
  void DiffusionDensityCalculator(int meshSize, double x, double startT,
                                  double endT, bool Absorption,
                                  string &Filename, bool verbose,
                                  double diffusion_threshold,
                                  double bridge_threshold);
  void BridgeDiffusionDensityCalculator(int meshSize, double x, double z,
                                        double startT, double endT,
                                        double sampleT, bool Absorption,
                                        string &Filename, bool verbose,
                                        double diffusion_threshold,
                                        double bridge_threshold);
  void DrawBridgeDiffTheta(string Filename, int nSim, size_t N_i, double x,
                           double z, double startT, double endT, double sampleT,
                           double diffusion_threshold, double bridge_threshold);
  void BridgeDiffusionDiffThetaDensityCalculator(string Filename, size_t N_i,
                                                 int meshSize, double x,
                                                 double z, double startT,
                                                 double endT, double sampleT,
                                                 double diffusion_threshold,
                                                 double bridge_threshold);

 private:
  /// WRIGHT-FISHER PROPERTIES

  vector<vector<double>> thetaP;
  bool non_neutral;
  vector<double> theta, sigma;
  int SelectionSetup;
  double dominanceParameter;
  int SelPolyDeg;
  vector<double> selectionCoeffs;
  vector<Polynomial> SelectionFunction, PhiFunction, AtildeFunction;
  int thetaIndex;
  vector<vector<vector<double>>> akm;
  boost::random::mt19937 WF_gen;
  vector<double> factorials;
  vector<vector<double>> lg_theta1, lg_theta2, lg_theta;
  vector<double> changepts;

  /// HELPER FUNCTIONS

  template <typename T>
  T Getlogakm(size_t theta_index, int k, int m);
  int radiate_from_mode(size_t N_i, int index, const double t) const;
  void increment_on_mk(bool isDiffTheta, size_t N_i, vector<int> &mk,
                       const double s, const double t) const;
  vector<int> radiate2d_mode(bool isDiffTheta, size_t N_i, size_t idx, double s,
                             double t) const;
  double Getd(size_t N_i, vector<double> &d, int i, double x, double z,
              double t);
  double GetdDiffTheta(size_t N_i, vector<double> &d, int i, double x, double z,
                       double s, double t, const Options &o);
  double GetdDiffThetaBoundaries(size_t N_i, vector<double> &d, int i, double x,
                                 double z, double s, double t,
                                 const Options &o);
  double GetdDiffThetaOneQApprox(size_t N_i, vector<double> &d, int i, double x,
                                 double z, double s, double t,
                                 const Options &o);
  double Getd2(size_t N_i, vector<double> &d, int i, double x, double t);
  double GetdBridgeSame(size_t N_i, vector<double> &d, int i, double x,
                        double t);
  double GetdBridgeInterior(size_t N_i, vector<double> &d, int i, double x,
                            double z, double t);
  double GetdBridgeUnconditional(size_t N_i, vector<double> &d, int i, double x,
                                 double z, double t);
  double computeA(size_t N_i, int m, int k, int l, int j, double x, double z);
  double computeAUnconditional(size_t N_i, int m, int k, int l, double x,
                               double z);
  double computeADiffTheta(size_t N_i, int m, int k, int l, int j, double x,
                           double z);
  vector<vector<double>> computeADiffTheta_new(size_t N_i, int m, int k,
                                               double x, double z);
  void precomputeADiffTheta(size_t N_i, int m, int k);
  void precomputeA(size_t N_i_s, size_t N_i_t, int m, int k);
  double computeAGammaLambda(size_t N_i, int gamma, int lambda, double x,
                             double z, double s, double t, const Options &o);
  double calculate_expectation(size_t N_i, int m, int k, double x, double z);
  double calculate_expectation_qApprox(size_t N_i, int k, double x, double z,
                                       double s, double t, const Options &o);
  int computeC(size_t N_i, int m, pair<vector<int>, double> &C);
  int computeF(size_t N_i, int m, pair<vector<int>, double> &C);
  int computeE(size_t N_i, pair<vector<int>, double> &C);
  int computeG(size_t N_i, double x, double z, double s, double t,
               const Options &o);
  int computeGBoundaries(size_t N_i, double x, double z, double s, double t,
                         const Options &o);
  int computeGZInterior(size_t N_i, double x, double z, double s, double t,
                        const Options &o);
  int computeGXInterior(size_t N_i, double x, double z, double s, double t,
                        const Options &o);
};

template <typename T>
T WrightFisher::Getlogakm(size_t theta_index, int k, int m) {
  double theta_s = theta[theta_index];
  if (k < m) {
    return static_cast<T>(0);
  }
  // Ensure outer dimension covers theta_index
  if (theta_index >= static_cast<int>(akm.size())) {
    akm.resize(theta_index + 1);
  }
  auto &akm_theta = akm[theta_index];
  // Ensure second dimension covers m
  if (static_cast<int>(akm_theta.size()) <= m) {
    akm_theta.resize(m + 1);
  }
  // For each row m0 = 0..m, ensure full columns 0..k and populate
  for (int m0 = 0; m0 <= m; ++m0) {
    auto &row = akm_theta[m0];
    int needed = k + 1;  // indices 0..k inclusive
    if (static_cast<int>(row.size()) < needed) {
      int old_size = static_cast<int>(row.size());
      row.resize(needed);
      // fill new entries
      for (int idx = old_size; idx < needed; ++idx) {
        if (idx < m0) {
          // entries k < m -> 0
          row[idx] = static_cast<T>(0);
        } else {
          int i = idx;
          T a;
          if (i == 0) {
            a = static_cast<T>(0);
          } else {
            a = std::log(theta_s + 2.0 * i - 1);
            for (int j = 2; j <= i; ++j) {
              a += std::log(theta_s + m0 + j - 2.0);
              if (j <= i - m0) a -= std::log(static_cast<T>(j));
              if (j <= m0) a -= std::log(static_cast<T>(j));
            }
          }
          row[idx] = a;
          // std::cout << "Setting akm[" << theta_index << "][" << m0 << "]["
          //           << idx << "] = " << a << std::endl;
        }
      }
    }
  }
  return akm_theta[m][k];
}

#endif