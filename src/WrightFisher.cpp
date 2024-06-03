#include "WrightFisher.h"

#include <time.h>

#include <algorithm>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <fstream>
#include <map>

#include "Polynomial.h"

/// HELPER FUNCTIONS

void WrightFisher::ThetaSetter()  /// Set theta depending on what thetaP is
                                  /// entered
{
  if (!thetaP.empty())  /// When thetaP is not empty, set theta to be the L1
                        /// norm of thetaP
  {
    theta = accumulate(thetaP.begin(), thetaP.end(), 0.0);
  } else {
    theta = 0.0;  /// Otherwise thetaP is empty and we set theta to be 0.0
  }
}

void WrightFisher::ThetaResetter()  /// For conditioned bridge ONLY - set thetaP
                                    /// empty to (2,2), or thetaP =
                                    /// (0,theta)/(theta,0) to
                                    /// (2,theta)/(theta,2)
{
  if (thetaP.empty())  /// If no mutation present, we need to set both thetaP
                       /// entries to 2
  {
    thetaP.push_back(2.0);
    thetaP.push_back(2.0);

    theta = 4.0;
  } else if (!(thetaP[0] > 0.0) ||
             !(thetaP[1] > 0.0))  /// If one sided mutation, we need to set
                                  /// corresponding parameter to 2
  {
    if (!(thetaP[0] > 0.0)) {
      thetaP[0] = 2.0;
    }
    if (!(thetaP[1] > 0.0)) {
      thetaP[1] = 2.0;
    }

    theta = accumulate(thetaP.begin(), thetaP.end(), 0.0);
  }
}

void WrightFisher::SelectionSetter()  /// Setting up selection mechanisms
{
  if (non_neutral)  /// If false, we don't need to run anything
  {
    Polynomial GenicSelection(-1.0, 1.0,
                              0.0);  /// Set up genic selection (which will
    /// always be present) as Polynomial class
    double100 AlphaC1 = (thetaP.empty() ? 0.0 : 0.5 * thetaP[0]),
              AlphaC2 = (thetaP.empty() ? 0.0 : -0.5 * theta);
    Polynomial Alpha(AlphaC2, AlphaC1);

    if (SelectionSetup == 0)  /// Genic Selection - sigma*x*(1-x)
    {
      SelPolyDeg = 2;
      SelectionFunction.Copy(0.5 * sigma * GenicSelection);
      Polynomial temp_Phi = 0.5 * (2.0 * 0.5 * sigma * Alpha +
                                   pow(0.5 * sigma, 2.0) * GenicSelection);
      Polynomial temp_A(0.5 * sigma, 0.0);  //= sigma*;
      PhiFunction.Copy(temp_Phi);
      AtildeFunction.Copy(temp_A);
    } else if (SelectionSetup ==
               1)  /// Diploid Selection - sigma*x*(1-x)*(h+x*(1-2*h))
    {
      SelPolyDeg = 4;
      Polynomial Eta(1.0 - 2.0 * dominanceParameter, dominanceParameter);
      SelectionFunction.Copy(0.5 * sigma * GenicSelection * Eta);
      Polynomial temp_Phi =
          0.5 * (2.0 * 0.5 * sigma * Alpha * Eta +
                 pow(0.5 * sigma, 2.0) * GenicSelection * Eta * Eta +
                 0.5 * sigma * GenicSelection * Eta.Derivative());
      Polynomial temp_A(0.25 * sigma * (1.0 - 2.0 * dominanceParameter),
                        0.5 * sigma * dominanceParameter, 0.0);
      PhiFunction.Copy(temp_Phi);
      AtildeFunction.Copy(temp_A);
    } else if (SelectionSetup ==
               2)  /// General polynomial selection - sigma*x*(1-x)*eta(x)
    {
      SelPolyDeg = 2 * selectionCoeffs.size();
      double *coeff_ptr = &selectionCoeffs[0];
      Polynomial Eta(coeff_ptr, SelPolyDeg);
      Polynomial temp_Phi =
          0.5 * (2.0 * 0.5 * sigma * Alpha * Eta +
                 pow(0.5 * sigma, 2.0) * GenicSelection * Eta * Eta +
                 0.5 * sigma * GenicSelection * Eta.Derivative());
      std::cout << temp_Phi.EvaluateReal(0.5) << endl;
      vector<double> ACoeffs(SelPolyDeg + 1);
      for (int i = 0; i <= SelPolyDeg; i++) {
        ACoeffs[i] =
            (0.5 * sigma * selectionCoeffs[i]) /
            (static_cast<double>(
                i + 1));  /// Scaling polynomials coefficients appropriately
      }
      coeff_ptr = &ACoeffs[0];
      Polynomial temp_A(coeff_ptr, SelPolyDeg);
      SelectionFunction.Copy(0.5 * sigma * GenicSelection * Eta);
      PhiFunction.Copy(temp_Phi);
      AtildeFunction.Copy(temp_A);
    } else {
      std::cout << "Please enter a valid selection setup!"
                << endl;  /// If values not within {0,1,2} complain!
    }
  }

  PhiSetter();  /// Run function to set phi values
}

void WrightFisher::PhiSetter()  /// Compute max and min of Phi & Atilde
                                /// functions
{
  if (non_neutral)  /// If false, we don't need to bother
  {
    if (SelectionSetup == 0)  /// Genic selection case
    {
      vector<double> GenicSelMinMax = PhitildeMinMaxRange();
      phiMin = GenicSelMinMax[0];
      phiMax = GenicSelMinMax[1];
      AtildeMax = 0.5 * sigma;
    } else  /// Otherwise we set Phi & Atilde fns, then use the
            /// PolynomialRootFinder to compute their extrema
    {
      Polynomial temp1 = PhiFunction.Derivative(),
                 temp2 = AtildeFunction.Derivative();
      vector<double> validRoots, Aroots;

      std::vector<double> real_vector, imag_vector, r_Avec, i_Avec;
      int degree = temp1.Degree(), Adegree = temp2.Degree();
      real_vector.resize(degree), imag_vector.resize(degree),
          r_Avec.resize(Adegree), i_Avec.resize(Adegree);

      double *real_vector_ptr = &real_vector[0];
      double *imag_vector_ptr = &imag_vector[0];
      double *r_Avec_ptr = &r_Avec[0];
      double *i_Avec_ptr = &i_Avec[0];
      int root_count = 0, root_Acount = 0;

      if (temp1.FindRoots(real_vector_ptr, imag_vector_ptr, &root_count) ==
          PolynomialRootFinder::SUCCESS) {
        int i = 0;

        for (i = 0; i < root_count; ++i) {
          if (!(imag_vector_ptr[i] > 0.0 && imag_vector_ptr[i] < 0.0) &&
              real_vector_ptr[i] >= 0.0 && real_vector_ptr[i] <= 1.0) {
            validRoots.push_back(real_vector_ptr[i]);
          }
        }
      }

      if (temp2.FindRoots(r_Avec_ptr, i_Avec_ptr, &root_Acount) ==
          PolynomialRootFinder::SUCCESS) {
        int i = 0;

        for (i = 0; i < root_Acount; ++i) {
          if (!(i_Avec_ptr[i] > 0.0 && i_Avec_ptr[i] < 0.0) &&
              r_Avec_ptr[i] >= 0.0 && r_Avec_ptr[i] <= 1.0) {
            Aroots.push_back(r_Avec_ptr[i]);
          }
        }
      }

      double potMax = max(PhiFunction.EvaluateReal(0.0),
                          PhiFunction.EvaluateReal(1.0)),
             potMin = min(PhiFunction.EvaluateReal(0.0),
                          PhiFunction.EvaluateReal(1.0));
      double AMax = max(AtildeFunction.EvaluateReal(0.0),
                        AtildeFunction.EvaluateReal(1.0));

      for (vector<double>::iterator vRit = validRoots.begin();
           vRit != validRoots.end(); vRit++) {
        potMax = max(potMax, PhiFunction.EvaluateReal(*vRit));
        potMin = min(potMin, PhiFunction.EvaluateReal(*vRit));
      }

      phiMax = potMax;
      phiMin = potMin;

      for (vector<double>::iterator vRit = Aroots.begin(); vRit != Aroots.end();
           vRit++) {
        AMax = max(AMax, AtildeFunction.EvaluateReal(*vRit));
      }

      AtildeMax = AMax;
    }
  }
}

vector<double100> WrightFisher::get_Theta() { return thetaP; }

double100 WrightFisher::Phitilde(
    double100 y)  /// Returns value of the quadratic function phitilde(y)
{
  assert((y >= 0.0) && (y <= 1.0));
  if (SelectionSetup == 0)  /// For genic selection
  {
    if (sigma == 0.0) {
      return 0.0;
    } else {
      return ((0.25) * sigma *
              (-(0.5) * sigma * y * y + ((0.5 * sigma) - theta) * y +
               (thetaP.empty() ? 0.0 : thetaP[0])));
    }
  } else  /// Otherwise use Polynomial class to evaluate
  {
    return PhiFunction.EvaluateReal(y);
  }
}

vector<double100>
WrightFisher::PhitildeMinMaxRange()  /// Returns the minimum, maximum and range
                                     /// of phitilde respectively for genic
                                     /// selection case
{
  double100 phiargmin,
      phiargmax = max(min(0.5 - (theta / sigma), 1.0),
                      0.0);  /// Find out max value of phitilde by plugging in
  /// derivative(phitilde)=0

  if (!(phiargmax > 0.0))  /// Figure out if max is at either boundary and set
                           /// min to be other
  {
    phiargmin = 1.0;
  } else if (!(phiargmax < 1.0)) {
    phiargmin = 0.0;
  } else  /// Otherwise we need to choose the smaller value at each endpoints as
          /// the min
  {
    double100 phizero = Phitilde(0.0), phione = Phitilde(1.0);

    if (min(phizero, phione) == phizero) {
      phiargmin = 0.0;
    } else {
      phiargmin = 1.0;
    }
  }

  vector<double100> MinMaxRange;

  MinMaxRange.push_back(Phitilde(phiargmin));
  MinMaxRange.push_back(Phitilde(phiargmax));
  MinMaxRange.push_back(MinMaxRange[1] - MinMaxRange[0]);

  return MinMaxRange;
}

double100 WrightFisher::Atilde(double100 x)  /// Returns the value of Atilde(x)
{
  assert((x >= 0.0) && (x <= 1.0));
  if (SelectionSetup == 0)  /// Genic selection case
  {
    return 0.5 * (sigma * x);
  } else  /// Otherwise use Polynomial class
  {
    return AtildeFunction.EvaluateReal(x);
  }
}

double100 WrightFisher::Atildeplus()  /// Returns the max of Atilde
{
  return AtildeMax;
}

pair<double, double> WrightFisher::GriffithsParas(
    double100 t)  /// Compute parameters for Griffiths approximation
{
  assert(t > 0.0);
  double beta = (theta - 1.0) * static_cast<double>(t) / 2.0;
  double eta = (abs(beta) <= 2.5e-5 ? 1.0 : beta / (exp(beta) - 1.0));
  double mu = 2 * (eta / static_cast<double>(t));
  double sigma =
      (abs(beta) <= 2.5e-5
           ? 2.0 / (3.0 * static_cast<double>(t))
           : 2.0 * eta / static_cast<double>(t) * pow(eta + beta, 2.0) *
                 (1.0 + eta / (eta + beta) - 2.0 * eta) / pow(beta, 2.0));
  assert(sigma > 0.0);
  return make_pair(mu, sigma);
}

int WrightFisher::radiate_from_mode(int index, const double100 t)
    const  /// Moving from mode as opposed to starting from m = 0 and working
           /// upwards
{
  double beta = (theta - 1.0) * static_cast<double>(t) / 2.0;
  double eta = (abs(beta) <= 2.5e-5 ? 1.0 : beta / (exp(beta) - 1.0));
  int mmode = static_cast<int>(round(2 * eta / static_cast<double>(t))),
      threshold = (thetaP.empty() ? 1 : 0), adjindex = index + threshold;
  if (adjindex > 2 * (mmode)-threshold) {
    return adjindex;
  } else {
    return mmode + (index % 2 == 0 ? 1 : -1) * ((index + 1) / 2);
  }
}

void WrightFisher::increment_on_mk(vector<int> &mk, const double100 s,
                                   const double100 t)
    const  /// Incrementing (m,k) in bridge sampler using bijective fn
{
  int &m_index = mk[2], &k_index = mk[3];
  --m_index;
  ++k_index;
  if (m_index < 0) {
    m_index = m_index + k_index + 1;
    k_index = 0;
  }
  mk[0] = radiate_from_mode(m_index, s);
  mk[1] = radiate_from_mode(k_index, t - s);
}

double100 WrightFisher::Getd(
    vector<double100> &d, int i, double100 x, double100 z,
    double100 t)  /// Compute contributions to denom for DrawBridgePMF cases
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset, c2 = m - j;
      boost::math::binomial_distribution<double100> B(c2, x);
      double100 Expected_Dirichlet = 0.0;
      for (int l = 0; l <= c2; ++l) {
        double para1 = (thetaP.empty() ? static_cast<double>(l)
                                       : static_cast<double>(thetaP[0] + l)),
               para2 =
                   (thetaP.empty() ? static_cast<double>(c2 - l)
                                   : static_cast<double>(thetaP[1] + c2 - l));

        boost::math::beta_distribution<double100> D(para1, para2);
        Expected_Dirichlet += pdf(B, l) * pdf(D, z);
      }

      d[i] += exp(Getlogakm<double100>(c1, c2) +
                  static_cast<double100>(-c1 * (c1 + theta - 1) * t / 2.0)) *
              Expected_Dirichlet;
      assert(exp(Getlogakm<double100>(c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double100 WrightFisher::Getd2(
    vector<double100> &d, int i, double100 x,
    double100 t)  /// Compute contributions to denom for
                  /// AncestralProcessConditional function
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;
    int bumper = (thetaP.empty()) ? 2 : 1;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + bumper + j + offset, c2 = m + bumper - j;
      double100 binomialremainder =
          (bumper == 2) ? (1.0 - pow(1.0 - x, static_cast<double100>(c2)) -
                           pow(x, static_cast<double100>(c2)))
                        : (1.0 - pow(1.0 - x, static_cast<double100>(c2)));
      boost::math::binomial_distribution<double100> B(
          c2, x);  // Could make these distributions static for the duration of
      // this x and z
      d[i] += exp(Getlogakm<double100>(c1, c2) +
                  static_cast<double100>(-c1 * (c1 + theta - 1) * t / 2.0)) *
              binomialremainder;
      assert(exp(Getlogakm<double100>(c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double100 WrightFisher::GetdBridgeSame(
    vector<double100> &d, int i, double100 x,
    double100
        t)  /// Compute contributions to denom for DrawBridgePMFSame function
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset, c2 = m - j;
      double100 addon =
          (thetaP[0] > 0.0 && thetaP[1] > 0.0
               ? (!(x > 0.0) ? (log(boost::math::tgamma_ratio(
                                    static_cast<double100>(theta + c2),
                                    static_cast<double100>(thetaP[1] + c2))) -
                                log(boost::math::tgamma(thetaP[0])))
                             : (log(boost::math::tgamma_ratio(
                                    static_cast<double100>(theta + c2),
                                    static_cast<double100>(thetaP[0] + c2))) -
                                log(boost::math::tgamma(thetaP[1]))))
               : (log(boost::math::tgamma_ratio(
                      static_cast<double100>(theta + c2),
                      static_cast<double100>(c2))) -
                  log(boost::math::tgamma(theta))));
      d[i] +=
          exp(Getlogakm<double100>(c1, c2) +
              static_cast<double100>(-c1 * (c1 + theta - 1) * t / 2.0) + addon);
      assert(exp(Getlogakm<double100>(c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double100 WrightFisher::GetdBridgeInterior(
    vector<double100> &d, int i, double100 x, double100 z,
    double100 t)  /// Compute contributions to denom for DrawBridgePMFInterior
                  /// function
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset, c2 = m - j;
      double100 para1, para2;

      if (!(x > 0.0)) {
        para1 = thetaP[0];
        para2 = thetaP[1] + c2;
      } else {
        para1 = thetaP[0] + c2;
        para2 = thetaP[1];
      }

      boost::math::beta_distribution<double100> BETA(para1, para2);
      double100 zcontribution = log(pdf(BETA, z));
      d[i] += exp(Getlogakm<double100>(c1, c2) +
                  static_cast<double100>(-c1 * (c1 + theta - 1) * t / 2.0) +
                  zcontribution);
      assert(exp(Getlogakm<double100>(c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double100 WrightFisher::GetdBridgeUnconditional(
    vector<double100> &d, int i, double100 x, double100 z,
    double100 t)  /// Compute contributions to denom for
                  /// DrawBridgePMFUnconditioned function
{
  int thetaDep = (thetaP.empty() ? 1 : 0);
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset + thetaDep, c2 = m - j + thetaDep;
      double100 xcontribution =
          (!(z > 0.0) ? static_cast<double100>(c2) * log(1.0 - x)
                      : static_cast<double100>(c2) * log(x));
      d[i] += exp(Getlogakm<double100>(c1, c2) +
                  static_cast<double100>(-c1 * (c1 + theta - 1) * t / 2.0) +
                  xcontribution);
      assert(exp(Getlogakm<double100>(c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double100 WrightFisher::computeA(
    int m, int k, int l, int j, double100 x,
    double100 z)  /// Compute weights for bridge diffusion decomposition
{
  assert((m >= 0) && (k >= 0) && (l >= 0) && (j >= 0) && (m >= l) && (k >= j) &&
         (x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0));
  double theta1 = thetaP[0], theta2 = thetaP[1];
  boost::math::binomial_distribution<double100> BIN(m, x);
  boost::math::beta_distribution<double100> BETA(theta1 + j, theta2 + k - j);

  double100 betaBinPdfs;

  if (!(x > 0.0) || !(x < 1.0)) {
    betaBinPdfs = 1.0;  /// If x = 0.0 or 1.0, then there is no contribution
  } else {
    betaBinPdfs = pdf(BIN, l);
  }

  if (!(z > 0.0) || !(z < 1.0))  /// Same as above but for z
  {
  } else {
    betaBinPdfs *= pdf(BETA, z);
  }

  double100 log_bin = LogBinomialCoefficientCalculator(k, j);
  double100 log_beta_1 =
      log(boost::math::beta<double100>(theta1 + l + j, theta2 + m - l + k - j));
  double100 log_beta_2 =
      log(boost::math::beta<double100>(theta1 + l, theta2 + m - l));

  return betaBinPdfs * exp(log_bin + log_beta_1 - log_beta_2);
}

double100 WrightFisher::computeAUnconditional(
    int m, int k, int l, double100 x,
    double100
        z)  /// Compute weights for unconditioned bridge diffusion decomposition
{
  assert((m >= 0) && (k >= 0) && (l >= 0) && (m >= l) && (x > 0.0) &&
         (x < 1.0) && (!(z > 0.0) || !(z < 1.0)));
  boost::math::binomial_distribution<double100> BIN(m, x);

  if ((l == 0) && (thetaP.empty() || !(thetaP[0] > 0.0))) {
    return pdf(BIN, l);
  } else if ((l == m) && (thetaP.empty() || !(thetaP[1] > 0.0))) {
    return pdf(BIN, l);
  }

  if (!(z > 0.0)) {
    return pdf(BIN, l) *
           exp(boost::math::lgamma(static_cast<double100>(theta + m - l + k)) +
               boost::math::lgamma(static_cast<double100>(theta + m)) -
               boost::math::lgamma(static_cast<double100>(theta + m + k)) -
               boost::math::lgamma(static_cast<double100>(theta + m - l)));
  } else if (!(z < 1.0)) {
    return pdf(BIN, l) *
           exp(boost::math::lgamma(static_cast<double100>(theta + l + k)) +
               boost::math::lgamma(static_cast<double100>(theta + m)) -
               boost::math::lgamma(static_cast<double100>(theta + l)) -
               boost::math::lgamma(static_cast<double100>(theta + k + m)));
  } else {
    cerr << "z was not 0 or 1! Returning 0.";
    return 0.0;
  }
}

int WrightFisher::computeC(
    int m, pair<vector<int>, double100> &C)  /// Compute the quantity C_m
{
  assert(m >= 0);
  if (m > static_cast<int>(C.first.size() - 1) || C.first[m] < 0) {
    C.first.resize(m + 1, -1);
    int i = 0;
    double100 bimnext =
                  exp(Getlogakm<double100>(i + m + 1, m) -
                      (i + m + 1) * (i + m + 1 + theta - 1) * C.second / 2.0),
              bim = exp(Getlogakm<double100>(i + m, m) -
                        (i + m) * (i + m + theta - 1) * C.second / 2.0);
    while (bimnext > bim) {
      ++i;
      bim = bimnext;
      bimnext = exp(Getlogakm<double100>(i + m + 1, m) -
                    (i + m + 1) * (i + m + 1 + theta - 1) * C.second / 2.0);
    }
    C.first[m] = i;
  }
  assert(C.first[m] >= 0);
  return C.first[m];
}

int WrightFisher::computeE(
    pair<vector<int>, double100> &C)  /// Compute the quantity E
{
  int next_constraint = computeC(0, C), curr_row = 0;
  int diag_index = next_constraint + (next_constraint % 2);
  bool Efound = false;
  while (!Efound) {
    ++curr_row;
    next_constraint = computeC(curr_row, C) + curr_row;
    if (diag_index - curr_row < next_constraint) {
      diag_index = next_constraint + curr_row;
      diag_index += (diag_index % 2);
    }
    if (curr_row == diag_index / 2) Efound = true;
  }
  return curr_row;
}

double100 WrightFisher::NormCDF(
    double100 x, double100 m,
    double100 v)  /// CDF for N(m,v) - v is the variance
{
  return 0.5 * erfc(-(x - m) / sqrt(2 * v));
}

double100 WrightFisher::DiscretisedNormCDF(
    int m,
    double100
        t)  /// CDF for discretised normal, binning mass into nearest integer
{
  double100 returnval;
  int threshold = (thetaP.empty()                           ? 2
                   : (thetaP[0] == 0 || !(thetaP[1] > 0.0)) ? 1
                                                            : 0);
  double beta = (theta - 1.0) * static_cast<double>(t) / 2.0;
  double eta = (abs(beta) <= 2.5e-5 ? 1.0 : beta / (exp(beta) - 1.0));
  double v =
      (abs(beta) <= 2.5e-5 ? 2.0 / (3.0 * static_cast<double>(t))
                           : 2.0 * (eta / static_cast<double>(t)) *
                                 pow((eta + beta) / beta, 2.0) *
                                 (1.0 + (eta / (eta + beta)) - 2.0 * eta));
  assert(v > 0.0);

  if (m == threshold) {
    returnval = NormCDF(static_cast<double>(threshold) + 0.5,
                        2.0 * (eta / static_cast<double>(t)), v);
  } else {
    returnval = NormCDF(static_cast<double>(m) + 0.5,
                        2.0 * eta / static_cast<double>(t), v) -
                NormCDF(static_cast<double>(m) - 0.5,
                        2.0 * eta / static_cast<double>(t), v);
  }

  return returnval;
}

double100 WrightFisher::LogBinomialCoefficientCalculator(
    int n, int k)  /// Calculate usual binomial, with approximations kicking in
                   /// for large n and k
{
  assert(n >= 0 && k >= 0 && k <= n);
  if (n <= 1000) {
    return log(boost::math::binomial_coefficient<double100>(n, k));
  } else {  // Compute log approximation using Stirling's formula
    return static_cast<double100>(n) * log(static_cast<double100>(n)) -
           static_cast<double100>(k) * log(static_cast<double100>(k)) -
           static_cast<double100>(n - k) * log(static_cast<double100>(n - k)) +
           0.5 * (log(static_cast<double100>(n)) -
                  log(static_cast<double100>(k)) -
                  log(static_cast<double100>(n - k)) -
                  log(2 * boost::math::constants::pi<double100>()));
  }
}

double100 WrightFisher::UnconditionedDiffusionDensity(
    double100 x, double100 y, double100 t,
    const Options &o)  /// Compute truncation to diffusion transition density
                       /// where diffusion can be absorbed at any time
{
  assert((x > 0.0) && (x < 1.0) && (y >= 0.0) && (y <= 1.0) && (t > 0.0));

  int thetaDependent =
      (thetaP.empty() ? 1 : 0);  /// Check what thetaP configuration we have
  double100 density = 0.0, density_inc, threshold = 1.0e-12;
  int mMode = static_cast<int>(floor(GriffithsParas(t).first)), m = mMode,
      mFlip = 1, mU = 0, mD = 0;
  bool mSwitch = false, mUpSwitch = false, mDownSwitch = false;

  while (!mSwitch) {
    if (thetaP.empty()) {
      double100 addon =
          (!(y > 0.0) ? pow(1.0 - x, static_cast<double100>(m))
                      : (!(y < 1.0) ? pow(x, static_cast<double100>(m)) : 0.0));
      for (int l = 1; l <= m - 1; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<> BETA(static_cast<double100>(l),
                                              static_cast<double100>(m - l));
        addon += pdf(BIN, l) * pdf(BETA, y);
      }
      density_inc = QmApprox(m, t, o) * addon;
    } else if (!(thetaP[0] > 0.0)) {
      double100 addon =
          (!(y > 0.0) ? pow(1.0 - x, static_cast<double100>(m)) : 0.0);
      for (int l = 1; l <= m; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<> BETA(
            static_cast<double100>(l), static_cast<double100>(theta + m - l));
        addon += pdf(BIN, l) * pdf(BETA, y);
      }
      density_inc = QmApprox(m, t, o) * addon;
    } else {
      double100 addon = (!(y < 1.0) ? pow(x, static_cast<double100>(m)) : 0.0);
      for (int l = 0; l <= m - 1; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<> BETA(static_cast<double100>(theta + l),
                                              static_cast<double100>(m - l));
        addon += pdf(BIN, l) * pdf(BETA, y);
      }
      density_inc = QmApprox(m, t, o) * addon;
    }

    density += density_inc;

    if (!(mDownSwitch))  /// Switching mechanism
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch =
            ((density_inc < threshold) || (mMode - mD - 1 < thetaDependent));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = (density_inc < threshold);
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= thetaDependent)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }

  return density;
}

double100 WrightFisher::DiffusionDensityApproximationDenom(
    double100 x, double100 t,
    const Options
        &o)  /// Compute denominator for truncation to diffusion transition
             /// density when diffusion conditioned on non-absorption
{
  assert((x >= 0.0) && (x <= 1.0) && (t > 0.0));
  int thetaDependent =
      (thetaP.empty() ? 2
                      : ((!(thetaP[0] > 0.0) || !(thetaP[1] > 0.0)) ? 1 : 0));
  double100 denom;

  if ((thetaDependent == 1) ||
      (thetaDependent == 2))  /// Check whether we have a zero mutation entry
  {
    double100 denom_inc = 1.0;
    denom = 0.0;
    int dMode = static_cast<int>(ceil(GriffithsParas(t).first)), d = dMode,
        Dflip = 1, Djm = 0, Djp = 0;

    while (denom_inc >
           0.0)  /// As long as increments are positive, keep computing
    {
      if (!(x > 0.0) || !(x < 1.0)) {
        denom_inc = QmApprox(d, t, o) * static_cast<double100>(d);
      } else {
        denom_inc =
            QmApprox(d, t, o) * (1.0 - pow(x, static_cast<double100>(d)) -
                                 pow(1.0 - x, static_cast<double100>(d)));
      }

      denom += denom_inc;

      if (Dflip == -1 &&
          (dMode - Djm - 1 >= thetaDependent))  /// Mechanism to explore either
                                                /// side around the mode
      {
        Djm++;
        d = dMode - Djm;
      } else {
        Djp++;
        d = dMode + Djp;
      }
      assert(d >= thetaDependent);
      Dflip *= -1;
    }
  } else {
    return 1.0;
  }

  return denom;
}

double100 WrightFisher::DiffusionDensityApproximation(
    double100 x, double100 y, double100 t,
    const Options &o)  /// Compute truncation to diffusion transition density
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (t > 0.0));
  int thetaDependent =
      (thetaP.empty() ? 2
                      : ((thetaP[0] > 0.0 && thetaP[1] > 0.0)
                             ? 0
                             : 1));  /// Check what thetaP configuration we have
  double100 density = 0.0, density_inc,
            denom = DiffusionDensityApproximationDenom(x, t, o),
            threshold = max(1.0e-12 * denom, 1.0e-50);
  int mMode = static_cast<int>(floor(GriffithsParas(t).first)), m = mMode,
      mFlip = 1, mU = 0, mD = 0;
  bool mSwitch = false, mUpSwitch = false, mDownSwitch = false;

  while (!mSwitch) {
    if (thetaP.empty()) {
      if (!(x > 0.0)) {
        density_inc = QmApprox(m, t, o) * static_cast<double100>(m * (m - 1)) *
                      pow(1.0 - y, static_cast<double100>(m - 2));
      } else if (!(x < 1.0)) {
        density_inc = QmApprox(m, t, o) * static_cast<double100>(m * (m - 1)) *
                      pow(y, static_cast<double100>(m - 2));
      } else {
        double100 addon = 0.0;
        for (int l = 1; l != m; l++) {
          boost::math::binomial_distribution<> BIN(m, x);
          boost::math::beta_distribution<double100> BETA(
              static_cast<double100>(l), static_cast<double100>(m - l));
          if (l < 1 && y == 0.0) {
            addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
          } else if (m - l < 1 && y == 1.0) {
            addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
          } else {
            addon += pdf(BIN, l) * pdf(BETA, y);
          }
        }
        density_inc = QmApprox(m, t, o) * addon;
      }
    } else if (!(thetaP[0] > 0.0))  //|| !( thetaP[1] > 0.0 ) )
    {
      /*bool thetaSwitch = false;
      if ( !( thetaP[1] > 0.0 ) )
      {
          iter_swap(thetaP.begin(),thetaP.begin()+1);
          thetaSwitch = true;
      }*/

      if (!(x > 0.0)) {
        density_inc = QmApprox(m, t, o) *
                      static_cast<double100>(m * (theta + m - 1)) *
                      pow(1.0 - y, static_cast<double100>(theta + m - 2));
      } else {
        double100 addon = 0.0;
        for (int l = 1; l != m + 1; l++) {
          boost::math::binomial_distribution<> BIN(m, x);
          boost::math::beta_distribution<double100> BETA(
              static_cast<double100>(l), static_cast<double100>(theta + m - l));
          if (l < 1 && y == 0.0) {
            addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
          } else if (theta + m - l < 1.0 && y == 1.0) {
            addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
          } else {
            addon += pdf(BIN, l) * pdf(BETA, y);
          }
        }
        density_inc = QmApprox(m, t, o) * addon;
      }

      /*if ( thetaSwitch )
      {
          iter_swap(thetaP.begin(),thetaP.begin()+1);
      }*/
    } else if (!(thetaP[1] > 0.0)) {
      if (!(x < 1.0)) {
        density_inc = QmApprox(m, t, o) *
                      static_cast<double100>(m * (theta + m - 1)) *
                      pow(y, static_cast<double100>(theta + m - 2));
      } else {
        double100 addon = 0.0;
        for (int l = 0; l != m; l++) {
          boost::math::binomial_distribution<> BIN(m, x);
          boost::math::beta_distribution<double100> BETA(
              static_cast<double100>(theta + l), static_cast<double100>(m - l));
          if (theta + l < 1.0 && y == 0.0) {
            addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
          } else if (m - l < 1 && y == 1.0) {
            addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
          } else {
            addon += pdf(BIN, l) * pdf(BETA, y);
          }
        }
        density_inc = QmApprox(m, t, o) * addon;
      }
    } else {
      double100 addon = 0.0;
      for (int l = 0; l != m + 1; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<double100> BETA(
            static_cast<double100>(thetaP[0] + l),
            static_cast<double100>(thetaP[1] + m - l));
        if (thetaP[0] + l < 1.0 && y == 0.0) {
          addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
        } else if (thetaP[1] + m - l < 1.0 && y == 1.0) {
          addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
        } else {
          addon += pdf(BIN, l) * pdf(BETA, y);
        }
      }
      density_inc = QmApprox(m, t, o) * addon;
    }

    density += density_inc / denom;

    if (!(mDownSwitch))  /// Switching mechanism
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch =
            ((density_inc < threshold) || (mMode - mD - 1 < thetaDependent));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = (density_inc < threshold);
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= thetaDependent)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }

  return density;
}

double100 WrightFisher::QmApprox(
    int m, double100 t,
    const Options &o)  /// Compute an approximation to q_m(t)
{
  assert((m >= 0) && (t > 0.0));
  double100 qm = 0.0, qmold = -1.0;

  if (t <= o.g1984threshold)  /// If time increment is too small, use
                              /// discretised Gaussian
  {
    qm = DiscretisedNormCDF(m, t);
  } else {
    int mkIndex = m;
    while (abs(qm - qmold) > 1.0e-12 || qm < 0.0 ||
           qm > 1.0)  /// If increments are big enough, keep going
    {
      qmold = qm;
      qm += pow(-1.0, mkIndex - m) *
            exp(Getlogakm<double100>(mkIndex, m) +
                static_cast<double100>(-(mkIndex) * (mkIndex + theta - 1) *
                                       (t) / 2.0));
      mkIndex++;

      if (!(qm > qmold) && !(qm < qmold) &&
          (qm < 0.0 || qm > 1.0))  /// We have lost precision, so use
                                   /// discretised normal approximation
      {
        return DiscretisedNormCDF(m, t);
      }
    }
  }
  assert(qm >= 0.0);

  return qm;
}

double100 WrightFisher::UnconditionedBridgeDensity(
    double100 x, double100 z, double100 y, double100 s, double100 t,
    const Options
        &o)  /// Compute an approximation to the transition density of the
             /// diffusion bridge when absorption can happen at any time
{
  assert((x > 0.0) && (x < 1.0) && (y >= 0.0) && (y <= 1.0) && (s > 0.0) &&
         (s < t));

  if (thetaP.empty() &&
      ((!(z > 0.0) && !(y < 1.0)) || (!(z < 1.0) && !(y > 0.0)))) {
    return 0.0;
  }

  double100 eC = 0.0, denom_inc = 1.0;
  int dmode = static_cast<int>(ceil(GriffithsParas(t).first)), d = dmode,
      Dflip = 1, Djm = 0, Djp = 0, mkdLower = (thetaP.empty() ? 1 : 0);

  while (denom_inc >
         0.0)  /// As long as increments are positive, keep computing
  {
    if (!(z > 0.0))  /// z = 0
    {
      denom_inc = QmApprox(d, t, o) * pow(1.0 - x, static_cast<double100>(d));
    } else  /// z = 1
    {
      denom_inc = QmApprox(d, t, o) * pow(x, static_cast<double100>(d));
    }

    eC += denom_inc;

    if (Dflip == -1 &&
        (dmode - Djm - 1 >=
         mkdLower))  /// Mechanism to explore either side around the mode
    {
      Djm++;
      d = dmode - Djm;
    } else {
      Djp++;
      d = dmode + Djp;
    }
    assert(d >= 0);
    Dflip *= -1;
  }

  if ((!(y > 0.0) && ((thetaP.empty() || !(thetaP[0] > 0.0)))) ||
      (!(y < 1.0) && ((thetaP.empty() || !(thetaP[1] > 0.0))))) {
    int mMode = GriffithsParas(s).first,
        kMode = GriffithsParas(t - s).first;  /// Use these together eC to get
    /// estimate of suitable threshold
    double100 constcontr =
        ((!(y > 0.0) && ((thetaP.empty() || !(thetaP[0] > 0.0))))
             ? static_cast<double100>(mMode) * log(1.0 - x)
             : static_cast<double100>(mMode) * log(x));
    double100 density = 0.0,
              threshold =
                  max(exp(log(max(1.0e-300, QmApprox(mMode, s, o))) +
                          log(max(1.0e-300, QmApprox(kMode, t - s, o))) +
                          constcontr - log(eC)) *
                          1.0e-20,
                      1.0e-50);

    int m = mMode, mFlip = 1, mD = 0, mU = 0;
    bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;

    while (!mSwitch) {
      int k = kMode, kFlip = 1, kD = 0, kU = 0;
      bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;

      while (!kSwitch) {
        double100 num_inc;
        if (!(y > 0.0) && ((thetaP.empty() || !(thetaP[0] > 0.0)))) {
          num_inc =
              exp(log(max(1.0e-300, QmApprox(m, s, o))) +
                  log(max(1.0e-300, QmApprox(k, t - s, o))) +
                  static_cast<double100>(m) * log(1.0 - x) -
                  log(eC));  /// Putting all the separate contributions together
        } else {
          num_inc =
              exp(log(max(1.0e-300, QmApprox(m, s, o))) +
                  log(max(1.0e-300, QmApprox(k, t - s, o))) +
                  static_cast<double100>(m) * log(x) -
                  log(eC));  /// Putting all the separate contributions together
        }

        density += num_inc;

        if (!(kDownSwitch))  /// Switching mechanism for k
        {
          if (sgn(k - kMode) <= 0) {
            kDownSwitch =
                ((num_inc < threshold) || (kMode - kD - 1 < mkdLower));
          }
        }

        if (!(kUpSwitch)) {
          if (sgn(k - kMode) >= 0) {
            kUpSwitch = (num_inc < threshold);
          }
        }

        kSwitch = (kDownSwitch && kUpSwitch);

        if (!kSwitch) {
          if (kFlip == 1) {
            kU++;
            k = kMode + kU;
            kFlip *= (kDownSwitch ? 1 : -1);
          } else if ((kFlip == -1) && (kMode - kD - 1 >= mkdLower)) {
            kD++;
            k = kMode - kD;
            kFlip *= (kUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(mDownSwitch))  /// Switching mechanism for m
      {
        if (sgn(m - mMode) <= 0) {
          mDownSwitch =
              (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < mkdLower));
        }
      }

      if (!(mUpSwitch)) {
        if (sgn(m - mMode) >= 0) {
          mUpSwitch = ((kU == 0) && (kD == 0));
        }
      }

      mSwitch = (mDownSwitch && mUpSwitch);

      if (!mSwitch) {
        if (mFlip == 1) {
          mU++;
          m = mMode + mU;
          mFlip *= (mDownSwitch ? 1 : -1);
        } else if ((mFlip == -1) && (mMode - mD - 1 >= mkdLower)) {
          mD++;
          m = mMode - mD;
          mFlip *= (mUpSwitch ? 1 : -1);
        }
      }
    }
    assert(density >= 0.0);

    return density;
  } else if (thetaP.empty() || !(thetaP[0] > 0.0)) {
    vector<int> modeGuess = mklModeFinder(
        x, z, s, t, o);  /// Compute mode over (m,k,l) to use as start points
    int mMode = modeGuess[0], kMode = modeGuess[1],
        lMode = modeGuess[2];  /// Use these estimates together with eC as a
    /// gauge for suitable threshold
    double100 ycontr, xcontr;
    if (((lMode == 0) && (!(z > 0.0))) || ((lMode == mMode) && (!(z < 1.0)))) {
      xcontr =
          static_cast<double100>(mMode) * (!(z > 0.0) ? log(1.0 - x) : log(x));
      ycontr = 0.0;
    } else {
      boost::math::binomial_distribution<> BIN(mMode, x);
      xcontr = log(pdf(BIN, lMode));
      double100 p1 = static_cast<double100>(lMode);
      double100 p2 = static_cast<double100>(theta + mMode - lMode);
      boost::math::beta_distribution<double100> BETA(p1, p2);
      ycontr = log(pdf(BETA, y)) + static_cast<double100>(kMode) *
                                       (!(z > 0.0) ? log(1.0 - y) : log(y));
    }
    double100 density = 0.0,
              threshold =
                  max(exp(log(max(1.0e-300, QmApprox(mMode, s, o))) +
                          log(max(1.0e-300, QmApprox(kMode, t - s, o))) +
                          xcontr + ycontr - log(eC)) *
                          1.0e-20,
                      1.0e-50);

    int m = mMode, mFlip = 1, mD = 0, mU = 0;
    bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
    double100 mContr_D = boost::math::lgamma(static_cast<double100>(theta + m)),
              mContr_U = mContr_D, mContr;
    double100 constContr = -log(y) - log(1.0 - y);
    /// Allows to avoid calling beta functions, and thus faster
    while (!mSwitch) {
      double100 qm = max(1.0e-300, QmApprox(m, s, o));
      if (m != mMode)  /// Increment m contributions accordingly
      {
        if (mU > mD) {
          mContr_U += log(static_cast<double100>(theta + (m - 1)));
          mContr = log(qm) + mContr_U;
        } else {
          mContr_D -= log(static_cast<double100>(theta + (m + 1) - 1));
          mContr = log(qm) + mContr_D;
        }
      } else {
        mContr = log(qm) + mContr_U;
      }

      int k = kMode, kFlip = 1, kD = 0, kU = 0;
      bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
      double100 kContr_D = static_cast<double100>(k) *
                           (!(z > 0.0) ? log(1.0 - y) : log(y)),
                kContr_U = kContr_D, kContr;  /// Calculate k contributions

      while (!kSwitch) {
        double100 qk = max(1.0e-300, QmApprox(k, t - s, o));
        if (k != kMode) {
          if (kU > kD) {
            kContr_U += (!(z > 0.0) ? log(1.0 - y) : log(y));
            kContr = log(qk) + kContr_U;
          } else {
            kContr_D -= (!(z > 0.0) ? log(1.0 - y) : log(y));
            kContr = log(qk) + kContr_D;
          }
        } else {
          kContr = log(qk) + kContr_U;
        }

        int lLower =
                (thetaP.empty() ? 1
                                : ((!(thetaP[0] > 0.0) && !(z > 0.0)) ? 1 : 0)),
            lUpper = (thetaP.empty()
                          ? m - 1
                          : ((!(thetaP[1] > 0.0) && !(z < 1.0)) ? m - 1 : m));
        int lFlip = 1, newlMode = min(lMode, lUpper), l = newlMode, lU = 0,
            lD = 0;  /// Need to ensure l <= m as m changes!
        bool lSwitch = false, lDownSwitch = false,
             lUpSwitch = false;  /// Compute l contributions
        boost::math::binomial_distribution<double100> BIN(m, x);
        double100 lContr_D =
            log(pdf(BIN, l)) + static_cast<double100>(l) * log(y) +
            static_cast<double100>(m - l) * log(1.0 - y) -
            boost::math::lgamma(static_cast<double100>(l)) -
            boost::math::lgamma(static_cast<double100>(theta + m - l));
        double100 lContr_U = lContr_D, lContr;

        while (!lSwitch) {
          if (l != newlMode) {
            if (lU > lD) {
              lContr_U += log(static_cast<double100>(m - (l - 1))) -
                          log(static_cast<double100>((l - 1) + 1)) + log(x) -
                          log(1.0 - x) + log(y) - log(1.0 - y) +
                          log(static_cast<double100>(theta + m - (l - 1) - 1)) -
                          log(static_cast<double100>((l - 1)));
              lContr = lContr_U;
            } else {
              lContr_D += log(static_cast<double100>(l + 1)) -
                          log(static_cast<double100>(m - (l + 1) + 1)) +
                          log(1.0 - x) - log(x) + log(1.0 - y) - log(y) +
                          log(static_cast<double100>((l + 1) - 1)) -
                          log(static_cast<double100>(theta + m - (l + 1)));
              lContr = lContr_D;
            }
          } else {
            lContr = lContr_U;
          }
          double100 density_inc =
              exp(constContr + mContr + kContr + lContr -
                  log(eC));  /// Putting all separate contributions together
          density += density_inc;

          if (!(lDownSwitch))  /// Switching mechanism for l
          {
            if (sgn(l - newlMode) <= 0) {
              lDownSwitch =
                  ((density_inc < threshold) || (newlMode - lD - 1) < lLower);
            }
          }

          if (!(lUpSwitch)) {
            if (sgn(l - newlMode) >= 0) {
              lUpSwitch =
                  ((density_inc < threshold) || (newlMode + lU + 1) > lUpper);
            }
          }

          lSwitch = (lDownSwitch && lUpSwitch);

          if (!lSwitch) {
            if ((lFlip == 1 && (newlMode + lU + 1 <= lUpper)) ||
                (lDownSwitch && !(lUpSwitch))) {
              lU++;
              l = newlMode + lU;
              lFlip *= (lDownSwitch ? 1 : -1);
            } else if ((lFlip == -1 && (newlMode - lD - 1 >= lLower)) ||
                       (lUpSwitch && !(lDownSwitch))) {
              lD++;
              l = newlMode - lD;
              lFlip *= (lUpSwitch ? 1 : -1);
            }
          }
        }

        if (!(kDownSwitch))  /// Switching mechanism for k
        {
          if (sgn(k - kMode) <= 0) {
            kDownSwitch =
                (((lU == 0) && (lD == 0)) || (kMode - kD - 1 < mkdLower));
          }
        }

        if (!(kUpSwitch)) {
          if (sgn(k - kMode) >= 0) {
            kUpSwitch = ((lU == 0) && (lD == 0));
          }
        }

        kSwitch = (kDownSwitch && kUpSwitch);

        if (!kSwitch) {
          if (kFlip == 1) {
            kU++;
            k = kMode + kU;
            kFlip *= (kDownSwitch ? 1 : -1);
          } else if ((kFlip == -1) && (kMode - kD - 1 >= mkdLower)) {
            kD++;
            k = kMode - kD;
            kFlip *= (kUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(mDownSwitch))  /// Switching mechanism for m
      {
        if (sgn(m - mMode) <= 0) {
          mDownSwitch =
              (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < mkdLower));
        }
      }

      if (!(mUpSwitch)) {
        if (sgn(m - mMode) >= 0) {
          mUpSwitch = ((kU == 0) && (kD == 0));
        }
      }

      mSwitch = (mDownSwitch && mUpSwitch);

      if (!mSwitch) {
        if (mFlip == 1) {
          mU++;
          m = mMode + mU;
          mFlip *= (mDownSwitch ? 1 : -1);
        } else if ((mFlip == -1) && (mMode - mD - 1 >= mkdLower)) {
          mD++;
          m = mMode - mD;
          mFlip *= (mUpSwitch ? 1 : -1);
        }
      }
    }
    assert(density >= 0.0);

    return density;
  } else {
    vector<int> modeGuess = mklModeFinder(
        x, z, s, t, o);  /// Compute mode over (m,k,l) to use as start points
    int mMode = modeGuess[0], kMode = modeGuess[1],
        lMode = modeGuess[2];  /// Use these estimates together with eC as a
    /// gauge for suitable threshold
    double100 xcontr, ycontr;
    if (lMode == 0) {
      xcontr = static_cast<double100>(mMode) * log(x);
      ycontr = 0.0;
    } else {
      boost::math::binomial_distribution<> BIN(mMode, x);
      xcontr = log(pdf(BIN, lMode));
      double100 p1 = static_cast<double100>(theta + lMode),
                p2 = static_cast<double100>(mMode - lMode);
      boost::math::beta_distribution<double100> B1(p1, p2);
      ycontr = log(pdf(B1, y)) + static_cast<double100>(kMode) * log(y);
    }
    double100 density = 0.0,
              threshold =
                  max(exp(log(max(1.0e-300, QmApprox(mMode, s, o))) +
                          log(max(1.0e-300, QmApprox(kMode, t - s, o))) +
                          xcontr + ycontr - log(eC)) *
                          1.0e-20,
                      1.0e-50);

    int m = mMode, mFlip = 1, mD = 0, mU = 0;
    bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
    double100 mContr_D = boost::math::lgamma(static_cast<double100>(theta + m)),
              mContr_U = mContr_D, mContr;
    double100 constContr = -log(y) - log(1.0 - y);
    /// Allows to avoid calling beta functions, and thus faster
    while (!mSwitch) {
      double100 qm = max(1.0e-300, QmApprox(m, s, o));
      if (m != mMode)  /// Increment m contributions accordingly
      {
        if (mU > mD) {
          mContr_U += log(static_cast<double100>(theta + (m - 1)));
          mContr = log(qm) + mContr_U;
        } else {
          mContr_D -= log(static_cast<double100>(theta + (m + 1) - 1));
          mContr = log(qm) + mContr_D;
        }
      } else {
        mContr = log(qm) + mContr_U;
      }

      int k = kMode, kFlip = 1, kD = 0, kU = 0;
      bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
      double100 kContr_D = static_cast<double100>(k) * log(y),
                kContr_U = kContr_D, kContr;  /// Calculate k contributions

      while (!kSwitch) {
        double100 qk = max(1.0e-300, QmApprox(k, t - s, o));
        if (k != kMode) {
          if (kU > kD) {
            kContr_U += log(y);
            kContr = log(qk) + kContr_U;
          } else {
            kContr_D -= log(y);
            kContr = log(qk) + kContr_D;
          }
        } else {
          kContr = log(qk) + kContr_U;
        }

        int lFlip = 1, newlMode = min(lMode, m), l = newlMode, lU = 0,
            lD = 0;  /// Need to ensure l <= m as m changes!
        int lLower =
                (thetaP.empty() ? 1
                                : ((!(thetaP[0] > 0.0) && !(z > 0.0)) ? 1 : 0)),
            lUpper = (thetaP.empty()
                          ? m - 1
                          : ((!(thetaP[1] > 0.0) && !(z < 1.0)) ? m - 1 : m));
        bool lSwitch = false, lDownSwitch = false,
             lUpSwitch = false;  /// Compute l contributions
        boost::math::binomial_distribution<double100> BIN(m, x);
        double100 lContr_D =
            log(pdf(BIN, l)) + static_cast<double100>(l) * log(y) +
            static_cast<double100>(m - l) * log(1.0 - y) -
            boost::math::lgamma(static_cast<double100>(theta + l)) -
            boost::math::lgamma(static_cast<double100>(m - l));
        double100 lContr_U = lContr_D, lContr;

        while (!lSwitch) {
          if (l != newlMode) {
            if (lU > lD) {
              lContr_U += log(static_cast<double100>(m - (l - 1))) -
                          log(static_cast<double100>((l - 1) + 1)) + log(x) -
                          log(1.0 - x) + log(y) - log(1.0 - y) +
                          log(static_cast<double100>(m - (l - 1) - 1)) -
                          log(static_cast<double100>(theta + (l - 1)));
              lContr = lContr_U;
            } else {
              lContr_D += log(static_cast<double100>(l + 1)) -
                          log(static_cast<double100>(m - (l + 1) + 1)) +
                          log(1.0 - x) - log(x) + log(1.0 - y) - log(y) +
                          log(static_cast<double100>(theta + (l + 1) - 1)) -
                          log(static_cast<double100>(m - (l + 1)));
              lContr = lContr_D;
            }
          } else {
            lContr = lContr_U;
          }
          double100 density_inc =
              exp(constContr + mContr + kContr + lContr -
                  log(eC));  /// Putting all separate contributions together
          density += density_inc;

          if (!(lDownSwitch))  /// Switching mechanism for l
          {
            if (sgn(l - newlMode) <= 0) {
              lDownSwitch =
                  ((density_inc < threshold) || (newlMode - lD - 1) < lLower);
            }
          }

          if (!(lUpSwitch)) {
            if (sgn(l - newlMode) >= 0) {
              lUpSwitch =
                  ((density_inc < threshold) || (newlMode + lU + 1) > lUpper);
            }
          }

          lSwitch = (lDownSwitch && lUpSwitch);

          if (!lSwitch) {
            if ((lFlip == 1 && (newlMode + lU + 1 <= lUpper)) ||
                (lDownSwitch && !(lUpSwitch))) {
              lU++;
              l = newlMode + lU;
              lFlip *= (lDownSwitch ? 1 : -1);
            } else if ((lFlip == -1 && (newlMode - lD - 1 >= lLower)) ||
                       (lUpSwitch && !(lDownSwitch))) {
              lD++;
              l = newlMode - lD;
              lFlip *= (lUpSwitch ? 1 : -1);
            }
          }
        }

        if (!(kDownSwitch))  /// Switching mechanism for k
        {
          if (sgn(k - kMode) <= 0) {
            kDownSwitch =
                (((lU == 0) && (lD == 0)) || (kMode - kD - 1 < mkdLower));
          }
        }

        if (!(kUpSwitch)) {
          if (sgn(k - kMode) >= 0) {
            kUpSwitch = ((lU == 0) && (lD == 0));
          }
        }

        kSwitch = (kDownSwitch && kUpSwitch);

        if (!kSwitch) {
          if (kFlip == 1) {
            kU++;
            k = kMode + kU;
            kFlip *= (kDownSwitch ? 1 : -1);
          } else if ((kFlip == -1) && (kMode - kD - 1 >= mkdLower)) {
            kD++;
            k = kMode - kD;
            kFlip *= (kUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(mDownSwitch))  /// Switching mechanism for m
      {
        if (sgn(m - mMode) <= 0) {
          mDownSwitch =
              (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < mkdLower));
        }
      }

      if (!(mUpSwitch)) {
        if (sgn(m - mMode) >= 0) {
          mUpSwitch = ((kU == 0) && (kD == 0));
        }
      }

      mSwitch = (mDownSwitch && mUpSwitch);

      if (!mSwitch) {
        if (mFlip == 1) {
          mU++;
          m = mMode + mU;
          mFlip *= (mDownSwitch ? 1 : -1);
        } else if ((mFlip == -1) && (mMode - mD - 1 >= mkdLower)) {
          mD++;
          m = mMode - mD;
          mFlip *= (mUpSwitch ? 1 : -1);
        }
      }
    }
    assert(density >= 0.0);

    return density;
  }
}

double100 WrightFisher::BridgeDenom(
    double100 x, double100 z, double100 y, double100 s, double100 t,
    const Options &o)  /// Compute an approximation to the denominator for the
                       /// transition density of the diffusion bridge
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  double100 denom = 0.0;
  double100 inc = 1.0;
  int dmode = static_cast<int>(ceil(GriffithsParas(t).first)), d = dmode,
      Dflip = 1, Djm = 0, Djp = 0;

  while (inc > 0.0)  /// As long as increments are positive, keep computing
  {
    if ((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0))  /// x,z in (0,1)
    {
      double100 betabin = 0.0;
      for (int f = 0; f != d + 1; f++) {
        boost::math::binomial_distribution<double100> BIN(d, x);
        double100 para1 = static_cast<double100>(thetaP[0] + f),
                  para2 = static_cast<double100>(thetaP[1] + d - f);
        boost::math::beta_distribution<double100> BETA(para1, para2);

        betabin += pdf(BIN, f) * pdf(BETA, z);
      }

      inc = QmApprox(d, t, o) * betabin;
    } else if ((z > 0.0) && (z < 1.0) &&
               (!(x > 0.0) || !(x < 1.0)))  /// x in {0,1}, z in (0,1)
    {
      double100 para1 = (!(x > 0.0) ? thetaP[0]
                                    : static_cast<double100>(thetaP[0] + d)),
                para2 = (!(x > 0.0) ? static_cast<double100>(thetaP[1] + d)
                                    : thetaP[1]);
      boost::math::beta_distribution<double100> BETAZ(para1, para2);
      inc = QmApprox(d, t, o) * pdf(BETAZ, z);
    } else if ((x > 0.0) && (x < 1.0) &&
               (!(z > 0.0) || !(z < 1.0)))  /// x in (0,1), z in {0,1}
    {
      double100 para1 = (!(z > 0.0) ? thetaP[0]
                                    : static_cast<double100>(thetaP[0] + d)),
                para2 = (!(z > 0.0) ? static_cast<double100>(thetaP[1] + d)
                                    : thetaP[1]);
      double100 xcon = (!(z > 0.0) ? pow(1.0 - x, static_cast<double100>(d))
                                   : pow(x, static_cast<double100>(d)));
      inc = QmApprox(d, t, o) * xcon *
            (1.0 / boost::math::beta<double100>(para1, para2));
    } else if (!(x < z) && !(x > z))  /// x,z in {0,1} and x=z
    {
      double100 para1 = (!(x > 0.0) ? thetaP[0]
                                    : static_cast<double100>(thetaP[0] + d)),
                para2 = (!(x > 0.0) ? static_cast<double100>(thetaP[1] + d)
                                    : thetaP[1]);
      inc = QmApprox(d, t, o) / boost::math::beta<double100>(para1, para2);
    } else  /// x,z in {0,1} and x!=z
    {
      return QmApprox(0, t, o) /
             boost::math::beta<double100>(thetaP[0], thetaP[1]);
    }

    denom += inc;

    if (Dflip == -1 &&
        (dmode - Djm > 0))  /// Mechanism to explore either side around the mode
    {
      Djm++;
      d = dmode - Djm;
    } else {
      Djp++;
      d = dmode + Djp;
    }
    assert(d >= 0);
    Dflip *= -1;
  }

  return denom;
}

double100 WrightFisher::ComputeDensity1(
    double100 x, double100 z, double100 y, double100 s, double100 t,
    const Options &o)  /// Compute bridge density when x,z in (0,1) (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }

  vector<int> modeGuess = mkljModeFinder(
      x, z, s, t, o);  /// Compute mode over (m,k,l,j) to use as start points
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2],
      jMode = modeGuess[3];  /// Use these estimates together with eC as a gauge
  /// for suitable threshold
  double100 qmmode = max(1.0e-300, QmApprox(mMode, s, o)),
            qkmode = max(1.0e-300, QmApprox(kMode, t - s, o));
  boost::math::binomial_distribution<> BINx(mMode, x), BINy(kMode, y);
  boost::math::beta_distribution<double100> BETAz(
      static_cast<double100>(thetaP[0] + jMode),
      static_cast<double100>(thetaP[1] + kMode - jMode)),
      BETAy(static_cast<double100>(thetaP[0] + lMode),
            static_cast<double100>(thetaP[1] + mMode - lMode));
  double100 density = 0.0, eC = BridgeDenom(x, z, y, s, t, o),
            threshold =
                max(exp(log(qmmode) + log(qkmode) + log(pdf(BINx, lMode)) +
                        log(pdf(BINy, jMode)) + log(pdf(BETAz, z)) +
                        log(pdf(BETAy, y)) - log(eC)) *
                        1.0e-6,
                    1.0e-50);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double100 mContr_D = boost::math::lgamma(
                static_cast<double100>(thetaP[0] + thetaP[1] + m)),
            mContr_U = mContr_D,
            mContr;  /// Compute contributions depending only on m
  /// Allows to avoid calling beta functions, and thus faster
  while (!mSwitch) {
    double100 qm = QmApprox(m, s, o);
    if (m != mMode)  /// Increment m contributions accordingly
    {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m - 1)));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D -=
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m + 1) - 1));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double100 kContr_D = boost::math::lgamma(
                  static_cast<double100>(thetaP[0] + thetaP[1] + k)),
              kContr_U = kContr_D, kContr;  /// Calculate k contributions

    while (!kSwitch) {
      double100 qk = QmApprox(k, t - s, o);
      if (k != kMode) {
        if (kU > kD) {
          kContr_U +=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k - 1)));
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D -=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k + 1) - 1));
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      int lFlip = 1, lU = 0, lD = 0, newlMode = min(lMode, m),
          l = newlMode;  /// Need to ensure l <= m since m changes!
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;
      boost::math::binomial_distribution<double100> BINL(
          m, x);  /// Contributions from l
      double100 lContr_D =
          (log(pdf(BINL, l)) -
           boost::math::lgamma(static_cast<double100>(thetaP[0] + l)) -
           boost::math::lgamma(static_cast<double100>(thetaP[1] + m - l)) +
           static_cast<double100>(l) * log(y) +
           static_cast<double100>(m - l) * log(1.0 - y));
      double100 lContr_U = lContr_D, lContr;

      while (!lSwitch) {
        assert((l >= 0) && (l <= m));
        if (l != newlMode) {
          if (lU > lD) {
            lContr_U +=
                log(static_cast<double100>(m - (l - 1))) -
                log(static_cast<double100>((l - 1) + 1)) + log(x) -
                log(1.0 - x) +
                log(static_cast<double100>(thetaP[1] + m - (l - 1) - 1)) -
                log(static_cast<double100>(thetaP[0] + (l - 1))) + log(y) -
                log(1.0 - y);
            lContr = lContr_U;
          } else {
            lContr_D += log(static_cast<double100>(l + 1)) -
                        log(static_cast<double100>(m - (l + 1) + 1)) +
                        log(1.0 - x) - log(x) +
                        log(static_cast<double100>(thetaP[0] + (l + 1) - 1)) -
                        log(static_cast<double100>(thetaP[1] + m - (l + 1))) +
                        log(1.0 - y) - log(y);
            lContr = lContr_D;
          }
        } else {
          lContr = lContr_U;
        }

        int jFlip = 1, jU = 0, jD = 0, newjMode = min(jMode, k),
            j = newjMode;  /// Need to ensure j <= k as k changes!
        bool jSwitch = false, jDownSwitch = false,
             jUpSwitch = false;  /// Compute j contributions

        double100 jContr_D =
            LogBinomialCoefficientCalculator(k, j) -
            boost::math::lgamma(static_cast<double100>(thetaP[0] + j)) -
            boost::math::lgamma(static_cast<double100>(thetaP[1] + k - j)) +
            static_cast<double100>(thetaP[0] + j - 1) * log(z) +
            static_cast<double100>(thetaP[1] + k - j - 1) * log(1.0 - z) +
            static_cast<double100>(thetaP[0] + j - 1) * log(y) +
            static_cast<double100>(thetaP[1] + k - j - 1) * log(1.0 - y);
        double100 jContr_U = jContr_D, jContr;

        while (!jSwitch) {
          if (j != newjMode) {
            if (jU > jD) {
              jContr_U +=
                  log(static_cast<double100>(k - (j - 1))) -
                  log(static_cast<double100>((j - 1) + 1)) + log(z) -
                  log(1.0 - z) +
                  log(static_cast<double100>(thetaP[1] + k - (j - 1) - 1)) -
                  log(static_cast<double100>(thetaP[0] + (j - 1))) + log(y) -
                  log(1.0 - y);
              jContr = jContr_U;
            } else {
              jContr_D += log(static_cast<double100>(j + 1)) -
                          log(static_cast<double100>(k - (j + 1) + 1)) +
                          log(1.0 - z) - log(z) +
                          log(static_cast<double100>(thetaP[0] + (j + 1) - 1)) -
                          log(static_cast<double100>(thetaP[1] + k - (j + 1))) +
                          log(1.0 - y) - log(y);
              jContr = jContr_D;
            }
          } else {
            jContr = jContr_U;
          }

          double100 density_inc =
              exp(mContr + kContr + lContr + jContr -
                  log(eC));  /// Put all separate contributions together
          density += density_inc;

          if (!(jDownSwitch))  /// Check whether we can still move downwards
          {
            if (sgn(j - newjMode) <= 0) {
              jDownSwitch =
                  ((density_inc < threshold) || (newjMode - jD - 1) < 0);
            }
          }

          if (!(jUpSwitch))  /// Check whether we can still move downwards
          {
            if (sgn(j - newjMode) >= 0) {
              jUpSwitch =
                  ((density_inc < threshold) || (newjMode + jU + 1) > k);
            }
          }

          jSwitch = (jDownSwitch &&
                     jUpSwitch);  /// Decide if we can move out and change l

          if (!jSwitch)  /// If we cannot, we need to move upwards or downwards
          {
            if ((jFlip == 1 && (newjMode + jU + 1 <= k)) ||
                (jDownSwitch && !(jUpSwitch))) {
              jU++;
              j = newjMode + jU;
              jFlip *= (jDownSwitch ? 1 : -1);
            } else if ((jFlip == -1 && (newjMode - jD - 1 >= 0)) ||
                       (jUpSwitch && !(jDownSwitch))) {
              jD++;
              j = newjMode - jD;
              jFlip *= (jUpSwitch ? 1 : -1);
            }
          }
        }

        if (!(lDownSwitch))  /// Same procedure as for j contributions, but now
                             /// we don't need to worry about increments being
                             /// below threshold
        {
          if (sgn(l - newlMode) <= 0) {
            lDownSwitch = (((jU == 0) && (jD == 0)) || (newlMode - lD - 1) < 0);
          }
        }

        if (!(lUpSwitch)) {
          if (sgn(l - newlMode) >= 0) {
            lUpSwitch = (((jU == 0) && (jD == 0)) || (newlMode + lU + 1) > m);
          }
        }

        lSwitch = (lDownSwitch && lUpSwitch);

        if (!lSwitch) {
          if ((lFlip == 1 && (newlMode + lU + 1 <= m)) ||
              (lDownSwitch && !(lUpSwitch))) {
            lU++;
            l = newlMode + lU;
            lFlip *= (lDownSwitch ? 1 : -1);
          } else if ((lFlip == -1 && (newlMode - lD - 1 >= 0)) ||
                     (lUpSwitch && !(lDownSwitch))) {
            lD++;
            l = newlMode - lD;
            lFlip *= (lUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(kDownSwitch))  /// Same switching procedure as for l
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch = (((lU == 0) && (lD == 0)) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = ((lU == 0) && (lD == 0));
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1) && (kMode - kD - 1 >= 0)) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Same switching procedure as for k
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < 0));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= 0)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }
  assert(density >= 0.0);

  return density;
}

double100 WrightFisher::ComputeDensity2(
    double100 x, double100 z, double100 y, double100 s, double100 t,
    const Options
        &o)  /// Compute bridge density when x in {0,1},z in (0,1) (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }

  vector<int> modeGuess = mkjModeFinder(
      x, z, s, t, o);  /// Compute mode over (m,k,j) to use as start points
  int mMode = modeGuess[0], kMode = modeGuess[1],
      jMode = modeGuess[2];  /// Use these estimates together with eC as a gauge
  /// for suitable threshold
  double100 qmmode = max(1.0e-300, QmApprox(mMode, s, o)),
            qkmode = max(1.0e-300, QmApprox(kMode, t - s, o));
  double100 p1 = static_cast<double100>(!(x > 0.0) ? thetaP[0]
                                                   : thetaP[0] + mMode),
            p2 = static_cast<double100>(!(x < 1.0) ? thetaP[1] + mMode
                                                   : thetaP[1]);
  boost::math::binomial_distribution<> BINy(kMode, y);
  boost::math::beta_distribution<double100> BETAz(
      static_cast<double100>(thetaP[0] + jMode),
      static_cast<double100>(thetaP[1] + kMode - jMode)),
      BETAy(p1, p2);
  double100 density = 0.0, eC = BridgeDenom(x, z, y, s, t, o),
            threshold =
                max(exp(log(qmmode) + log(qkmode) + log(pdf(BINy, jMode)) +
                        log(pdf(BETAz, z)) + log(pdf(BETAy, y)) - log(eC)) *
                        1.0e-6,
                    1.0e-50);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double100 mContr_D =
                boost::math::lgamma(
                    static_cast<double100>(thetaP[0] + thetaP[1] + m)) +
                (!(x > 0.0) ? -boost::math::lgamma(
                                  static_cast<double100>(thetaP[1] + m)) +
                                  static_cast<double100>(m) * log(1.0 - y)
                            : -boost::math::lgamma(
                                  static_cast<double100>(thetaP[0] + m)) +
                                  static_cast<double100>(m) * log(y)),
            mContr_U = mContr_D, mContr;
  double100 constContr =
      (static_cast<double100>(thetaP[0] - 1.0)) * log(z) +
      (static_cast<double100>(thetaP[1] - 1.0)) * log(1.0 - z) +
      (static_cast<double100>(thetaP[0] - 1.0)) * log(y) +
      (static_cast<double100>(thetaP[1] - 1.0)) * log(1.0 - y) +
      (!(x > 0.0) ? -boost::math::lgamma(static_cast<double100>(thetaP[0]))
                  : -boost::math::lgamma(static_cast<double100>(thetaP[1])));
  /// Allows to avoid calling beta functions, and thus faster
  while (!mSwitch) {
    double100 qm = QmApprox(m, s, o);
    if (m != mMode)  /// Increment m contributions accordingly
    {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m - 1))) +
            (!(x > 0.0)
                 ? -log(static_cast<double100>(thetaP[1] + (m - 1))) +
                       log(1.0 - y)
                 : -log(static_cast<double100>(thetaP[0] + (m - 1))) + log(y));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D +=
            -log(static_cast<double100>(thetaP[0] + thetaP[1] + (m + 1) - 1)) +
            (!(x > 0.0) ? log(static_cast<double100>(thetaP[1] + (m + 1) - 1)) -
                              log(1.0 - y)
                        : log(static_cast<double100>(thetaP[0] + (m + 1) - 1)) -
                              log(y));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double100 kContr_D = boost::math::lgamma(
                  static_cast<double100>(thetaP[0] + thetaP[1] + k)),
              kContr_U = kContr_D, kContr;  /// Calculate k contributions

    while (!kSwitch) {
      double100 qk = QmApprox(k, t - s, o);
      if (k != kMode) {
        if (kU > kD) {
          kContr_U +=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k - 1)));
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D -=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k + 1) - 1));
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      int jFlip = 1, newjMode = min(jMode, k), j = newjMode, jU = 0,
          jD = 0;  /// Need to ensure j <= k as k changes!
      bool jSwitch = false, jDownSwitch = false,
           jUpSwitch = false;  /// Compute j contributions

      double100 jContr_D =
          LogBinomialCoefficientCalculator(k, j) -
          boost::math::lgamma(static_cast<double100>(thetaP[0] + j)) -
          boost::math::lgamma(static_cast<double100>(thetaP[1] + k - j)) +
          static_cast<double100>(thetaP[0] + j) * log(z) +
          static_cast<double100>(thetaP[1] + k - j) * log(1.0 - z) +
          static_cast<double100>(j) * log(y) +
          static_cast<double100>(k - j) * log(1.0 - y);
      double100 jContr_U = jContr_D, jContr;

      while (!jSwitch) {
        if (j != newjMode) {
          if (jU > jD) {
            jContr_U +=
                log(static_cast<double100>(k - (j - 1))) -
                log(static_cast<double100>((j - 1) + 1)) + log(z) -
                log(1.0 - z) +
                log(static_cast<double100>(thetaP[1] + k - (j - 1) - 1)) -
                log(static_cast<double100>(thetaP[0] + (j - 1))) + log(y) -
                log(1.0 - y);
            jContr = jContr_U;
          } else {
            jContr_D += log(static_cast<double100>(j + 1)) -
                        log(static_cast<double100>(k - (j + 1) + 1)) +
                        log(1.0 - z) - log(z) +
                        log(static_cast<double100>(thetaP[0] + (j + 1) - 1)) -
                        log(static_cast<double100>(thetaP[1] + k - (j + 1))) +
                        log(1.0 - y) - log(y);
            jContr = jContr_D;
          }
        } else {
          jContr = jContr_U;
        }
        double100 density_inc =
            exp(constContr + mContr + kContr + jContr -
                log(eC));  /// Putting all separate contributions together
        density += density_inc;

        if (!(jDownSwitch))  /// Switching mechnism for j, exactly as in
                             /// ComputeDensity1 function
        {
          if (sgn(j - newjMode) <= 0) {
            jDownSwitch =
                ((density_inc < threshold) || (newjMode - jD - 1) < 0);
          }
        }

        if (!(jUpSwitch)) {
          if (sgn(j - newjMode) >= 0) {
            jUpSwitch = ((density_inc < threshold) || (newjMode + jU + 1) > k);
          }
        }

        jSwitch = (jDownSwitch && jUpSwitch);

        if (!jSwitch) {
          if ((jFlip == 1 && (newjMode + jU + 1 <= k)) ||
              (jDownSwitch && !(jUpSwitch))) {
            jU++;
            j = newjMode + jU;
            jFlip *= (jDownSwitch ? 1 : -1);
          } else if ((jFlip == -1 && (newjMode - jD - 1 >= 0)) ||
                     (jUpSwitch && !(jDownSwitch))) {
            jD++;
            j = newjMode - jD;
            jFlip *= (jUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch = (((jU == 0) && (jD == 0)) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = ((jU == 0) && (jD == 0));
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1) && (kMode - kD - 1 >= 0)) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Switching mechanism for m
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < 0));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= 0)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }
  assert(density >= 0.0);

  return density;
}

double100 WrightFisher::ComputeDensity3(
    double100 x, double100 z, double100 y, double100 s, double100 t,
    const Options
        &o)  /// Compute bridge density when x,z in {0,1}, x=z (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }

  vector<int> modeGuess = mkModeFinder(x, z, s, t, o);  /// Find mode over (m,k)
  int mMode = modeGuess[0],
      kMode = modeGuess[1];  /// Use these together eC to get estimate of
  /// suitable threshold
  double100 qmmode = max(1.0e-300, QmApprox(mMode, s, o)),
            qkmode = max(1.0e-300, QmApprox(kMode, t - s, o));
  double100 th1 = static_cast<double100>(!(x > 0.0) ? thetaP[0] : thetaP[1]),
            th2 = theta - th1;
  double100 gammaratios =
      -boost::math::lgamma(th1) +
      boost::math::lgamma(static_cast<double100>(theta + mMode)) -
      boost::math::lgamma(static_cast<double100>(th2 + mMode)) +
      boost::math::lgamma(static_cast<double100>(theta + kMode)) -
      boost::math::lgamma(static_cast<double100>(th2 + kMode)) +
      boost::math::lgamma(static_cast<double100>(th2 + mMode + kMode)) -
      boost::math::lgamma(static_cast<double100>(theta + mMode + kMode));
  boost::math::beta_distribution<double100> BETA(
      static_cast<double100>(
          !(x > 0.0 ? thetaP[0] : thetaP[0] + mMode + kMode)),
      static_cast<double100>(!(x > 0.0) ? thetaP[1] + mMode + kMode
                                        : thetaP[1]));
  double100 density = 0.0, eC = BridgeDenom(x, z, y, s, t, o),
            threshold = max(exp(log(qmmode) + log(qkmode) + gammaratios +
                                log(pdf(BETA, y)) - log(eC)) *
                                1.0e-6,
                            1.0e-50);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double100 constContr =
      static_cast<double100>(thetaP[0] - 1.0) * log(y) +
      static_cast<double100>(thetaP[1] - 1.0) * log(1.0 - y) +
      (!(x > 0.0)
           ? -2.0 * boost::math::lgamma(static_cast<double100>(thetaP[0]))
           : -2.0 * boost::math::lgamma(static_cast<double100>(thetaP[1])));
  double100 mContr_D =
      boost::math::lgamma(static_cast<double100>(thetaP[0] + thetaP[1] + m)) +
      (!(x > 0.0)
           ? -boost::math::lgamma(static_cast<double100>(thetaP[1] + m)) +
                 static_cast<double100>(m) * log(1.0 - y)
           : -boost::math::lgamma(static_cast<double100>(thetaP[0] + m)) +
                 static_cast<double100>(m) * log(y));
  double100 mContr_U = mContr_D, mContr;  /// Contributions from m

  while (!mSwitch) {
    double100 qm = QmApprox(m, s, o);
    if (m != mMode) {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m - 1))) +
            (!(x > 0.0)
                 ? -log(static_cast<double100>(thetaP[1] + (m - 1))) +
                       log(1.0 - y)
                 : -log(static_cast<double100>(thetaP[0] + (m - 1))) + log(y));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D +=
            (!(x > 0.0) ? log(static_cast<double100>(thetaP[1] + (m + 1) - 1)) -
                              log(1.0 - y)
                        : log(static_cast<double100>(thetaP[0] + (m + 1) - 1)) -
                              log(y)) -
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m + 1) - 1));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double100 kContr_D =
        (boost::math::lgamma(
             static_cast<double100>(thetaP[0] + thetaP[1] + k)) +
         (!(x > 0.0)
              ? -boost::math::lgamma(static_cast<double100>(thetaP[1] + k)) +
                    static_cast<double100>(k) * log(1.0 - y)
              : -boost::math::lgamma(static_cast<double100>(thetaP[0] + k)) +
                    static_cast<double100>(k) * log(y)));
    double100 kContr_U = kContr_D, kContr;  /// Contributions from k

    while (!kSwitch) {
      double100 qk = QmApprox(k, t - s, o);
      if (k != kMode) {
        if (kU > kD) {
          kContr_U +=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k - 1))) +
              (!(x > 0.0)
                   ? log(1.0 - y) -
                         log(static_cast<double100>(thetaP[1] + (k - 1)))
                   : log(y) - log(static_cast<double100>(thetaP[0] + (k - 1))));
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D +=
              -log(
                  static_cast<double100>(thetaP[0] + thetaP[1] + (k + 1) - 1)) +
              (!(x > 0.0)
                   ? log(static_cast<double100>(thetaP[1] + (k + 1) - 1)) -
                         log(1.0 - y)
                   : log(static_cast<double100>(thetaP[0] + (k + 1) - 1)) -
                         log(y));
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      double100 density_inc =
          exp(constContr + mContr + kContr -
              log(eC));  /// Putting all the separate contributions together
      density += density_inc;

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch = ((density_inc < threshold) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = (density_inc < threshold);
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1) && (kMode - kD - 1 >= 0)) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Switching mechanism for m
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < 0));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= 0)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }
  assert(density >= 0.0);

  return density;
}

double100 WrightFisher::ComputeDensity4(
    double100 x, double100 z, double100 y, double100 s, double100 t,
    const Options
        &o)  /// Compute bridge density when x,z in {0,1}, x!=z (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }

  vector<int> modeGuess = mkModeFinder(x, z, s, t, o);  /// Find mode over (m,k)
  int mMode = modeGuess[0],
      kMode = modeGuess[1];  /// Use together with eC to get suitable threshold
  double100 qmmode = max(1.0e-300, QmApprox(mMode, s, o)),
            qkmode = max(1.0e-300, QmApprox(kMode, t - s, o));
  double100 gammaratios =
      -boost::math::lgamma(static_cast<double100>(thetaP[0])) -
      boost::math::lgamma(static_cast<double100>(thetaP[1])) +
      boost::math::lgamma(static_cast<double100>(theta + mMode)) +
      boost::math::lgamma(static_cast<double100>(theta + kMode)) -
      boost::math::lgamma(static_cast<double100>(theta + mMode + kMode));
  boost::math::beta_distribution<double100> BETA(
      static_cast<double100>(
          !(x > 0.0 ? thetaP[0] + kMode : thetaP[0] + mMode)),
      static_cast<double100>(!(x > 0.0) ? thetaP[1] + mMode
                                        : thetaP[1] + mMode));
  double100 density = 0.0, eC = BridgeDenom(x, z, y, s, t, o),
            threshold = max(exp(log(qmmode) + log(qkmode) + gammaratios +
                                log(pdf(BETA, y)) - log(eC)) *
                                1.0e-6,
                            1.0e-50);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double100 constContr =
      static_cast<double100>(thetaP[0] - 1.0) * log(y) +
      static_cast<double100>(thetaP[1] - 1.0) * log(1.0 - y) -
      (boost::math::lgamma(static_cast<double100>(thetaP[0])) +
       boost::math::lgamma(static_cast<double100>(thetaP[1])));
  double100 mContr_D =
      boost::math::lgamma(static_cast<double100>(theta + m)) +
      (!(x > 0.0)
           ? -boost::math::lgamma(static_cast<double100>(thetaP[1] + m)) +
                 static_cast<double100>(m) * log(1.0 - y)
           : -boost::math::lgamma(static_cast<double100>(thetaP[0] + m)) +
                 static_cast<double100>(m) * log(y));
  double100 mContr_U = mContr_D, mContr;  /// Contributions from m

  while (!mSwitch) {
    double100 qm = QmApprox(m, s, o);
    if (m != mMode) {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double100>(theta + (m - 1))) +
            (!(x > 0.0)
                 ? -log(static_cast<double100>(thetaP[1] + (m - 1))) +
                       log(1.0 - y)
                 : -log(static_cast<double100>(thetaP[0] + (m - 1))) + log(y));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D +=
            (!(x > 0.0) ? log(static_cast<double100>(thetaP[1] + (m + 1) - 1)) -
                              log(1.0 - y)
                        : log(static_cast<double100>(thetaP[0] + (m + 1) - 1)) -
                              log(y)) -
            log(static_cast<double100>(theta + (m + 1) - 1));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double100 kContr_D =
        (boost::math::lgamma(static_cast<double100>(theta + k)) +
         (!(x > 0.0)
              ? -boost::math::lgamma(static_cast<double100>(thetaP[0] + k)) +
                    static_cast<double100>(k) * log(y)
              : -boost::math::lgamma(static_cast<double100>(thetaP[1] + k)) +
                    static_cast<double100>(k) * log(1.0 - y)));
    double100 kContr_U = kContr_D, kContr;  /// Contributions from k

    while (!kSwitch) {
      double100 qk = QmApprox(k, t - s, o);
      if (k != kMode) {
        if (kU > kD) {
          kContr_U +=
              log(static_cast<double100>(theta + (k - 1))) +
              (!(x > 0.0)
                   ? log(y) - log(static_cast<double100>(thetaP[0] + (k - 1)))
                   : log(1.0 - y) -
                         log(static_cast<double100>(thetaP[1] + (k - 1))));
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D +=
              -log(static_cast<double100>(theta + (k + 1) - 1)) +
              (!(x > 0.0)
                   ? log(static_cast<double100>(thetaP[0] + (k + 1) - 1)) -
                         log(y)
                   : log(static_cast<double100>(thetaP[1] + (k + 1) - 1)) -
                         log(1.0 - y));
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      double100 density_inc =
          exp(constContr + mContr + kContr -
              log(eC));  /// Putting all separate contributions together
      density += density_inc;

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch = ((density_inc < threshold) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = (density_inc < threshold);
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1) && (kMode - kD - 1 >= 0)) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Switching mechanism for m
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < 0));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= 0)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }
  assert(density >= 0.0);

  return density;
}

double100 WrightFisher::BridgeDensity(
    double100 x, double100 z, double100 y, double100 s, double100 t,
    const Options &o)  /// Function to determine which case from the above
                       /// computeDensity to invoke
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));

  if ((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0))  /// x,z in (0,1)
  {
    return ComputeDensity1(x, z, y, s, t, o);
  } else if ((x > 0.0) && (x < 1.0))  /// x in (0,1), z in {0,1}
  {
    return ComputeDensity2(z, x, y, t - s, t,
                           o);  /// Equivalent to a reverse bridge going from z
                                /// to y in time t-s and ending at x at time t
  } else if ((z > 0.0) && (z < 1.0))  /// x in {0,1}, z in (0,1)
  {
    return ComputeDensity2(x, z, y, s, t, o);
  } else if (x == z)  /// x,z in {0,1}, x=z
  {
    return ComputeDensity3(x, z, y, s, t, o);
  } else  /// x,z in {0,1}, x!=z
  {
    return ComputeDensity4(x, z, y, s, t, o);
  }
}

double100 WrightFisher::LogSumExp(
    vector<double100> &vecProb,
    double100 maxProb)  /// Routine to perform log-sum-exp trick
{
  double100 sumexp = 0.0;
  for (vector<double100>::iterator vPi = vecProb.begin(); vPi != vecProb.end();
       vPi++) {
    sumexp += exp(*vPi - maxProb);
  }
  for (vector<double100>::iterator vPi = vecProb.begin(); vPi != vecProb.end();
       vPi++) {
    *vPi -= maxProb + log(sumexp);
    *vPi = exp(*vPi);
  }

  return exp(maxProb + log(sumexp));
}

/// DIFFUSION SIMULATION - NEUTRAL PATHS

pair<int, int> WrightFisher::DrawAncestralProcess(
    double100 t, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of the Ancestral
                                  /// Process using upper and lower sums
{
  assert(t > 0.0);
  /// Pre-computing all necessary quantities
  int mup = -1, mdown = 0, mmode = GriffithsParas(t).first, m = -1;
  bool m_found = false;

  /// Setting up all the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen), currsumU = 0.0, currsumL = 0.0;
  map<int, double100> dL, dU;
  map<int, int> n_used;

  int threshold = ((thetaP.empty()) ? 1 : 0);

  while (!m_found) {
    if (mup > mdown &&
        mmode - mdown >
            threshold)  /// Checking if down moves still possible - ensures we
                        /// do not get m smaller than it should be allowed!
    {
      ++mdown;
      m = mmode - mdown;
    } else  /// Up moves can keep happening
    {
      ++mup;
      m = mmode + mup;
    }
    assert(m >= 0);

    int n = m - 2;
    bool coefficients_converging = false;
    dU.insert(make_pair(m, 0.0));
    dL.insert(make_pair(m, 0.0));

    while (n < 1 || !coefficients_converging ||
           ((dU[m]) > 1.0 &&
            dL[m] < 0.0))  /// Checking if upper and lower sums converging
    {
      n += 2;
      double100 newcoefficientU =
          (((n - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double100>(n, m) +
              static_cast<double100>(-n * (n + theta - 1) * t / 2.0));
      double100 newcoefficientL =
          (((n + 1 - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double100>(n + 1, m) +
              static_cast<double100>(-(n + 1) * ((n + 1) + theta - 1) * t /
                                     2.0));

      if (o.debug > 2) {
        cerr << setprecision(0) << n << " newcoeffU " << setprecision(50)
             << newcoefficientU << "\n"
             << n + 1 << " newcoeffL " << newcoefficientL << endl;
      }

      dU[m] = dL[m] + newcoefficientU;
      dL[m] = dU[m] + newcoefficientL;

      coefficients_converging = (newcoefficientU + newcoefficientL >= 0.0);

      if (o.debug > 2) {
        cerr << "\nk   ";
        for (int k = mmode - mdown; k <= mmode + mup; ++k) cerr << k << " ";
        cerr << "\ndL ";
        for (int k = mmode - mdown; k <= mmode + mup; ++k) cerr << dL[k] << " ";
        cerr << "\ndU ";
        for (int k = mmode - mdown; k <= mmode + mup; ++k) cerr << dU[k] << " ";
        cerr << endl;
      }
    }

    if (dL[m] > 1.0 ||
        dU[m] < 0.0)  /// If we fall out of bounds, use Gaussian approximations
    {
      cerr << "Numerical error detected: dL[" << m << "] = " << dL[m] << ", dU["
           << m << "] = " << dU[m]
           << ". Resorting to G1984 approximation (t = " << t << ") ..."
           << endl;
      return make_pair(DrawAncestralProcessG1984(t, gen), -1);
    }

    currsumL += dL[m];
    currsumU += dU[m];

    if (currsumL > currsumU) {
      cerr << "Error: currsumL = " << currsumL << " > " << currsumU
           << " = currsumU." << endl;
      exit(1);
    }

    n_used.insert(make_pair(m, n + 1));

    bool decision_on_m_made = (currsumL > u || currsumU < u);

    while (
        !decision_on_m_made)  /// We need to refine our upper and lower bounds
    {
      double100 currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      for (int k = mmode - mdown; k <= mmode + mup; ++k) {
        ++n_used[k];
        dU[k] = dL[k] +
                (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                    exp(Getlogakm<double100>(n_used[k], k) +
                        static_cast<double100>(
                            -n_used[k] * (n_used[k] + theta - 1) * t / 2.0));

        ++n_used[k];
        dL[k] = dU[k] +
                (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                    exp(Getlogakm<double100>(n_used[k], k) +
                        static_cast<double100>(
                            -n_used[k] * (n_used[k] + theta - 1) * t / 2.0));

        currsumL += dL[k];
        currsumU += dU[k];

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU." << endl;
          exit(1);
        }
      }

      if (currsumLold > currsumL) {
        // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
        std::cout << "Error: currsumLold = " << currsumLold << " > " << currsumL
                  << " = currsumL." << endl;
        // exit(1);
      }
      if (currsumUold < currsumU) {
        // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
        std::cout << "Error: currsumUold = " << currsumUold << " < " << currsumU
                  << " = currsumU." << endl;
        // exit(1);
      }

      decision_on_m_made = (currsumL > u || currsumU < u);
    }

    m_found = (currsumL > u);
  }

  int coeffcount = 0;
  for (int k = mmode - mdown; k <= mmode + mup; ++k)
    coeffcount += (n_used[k] - k + 1);

  if (o.debug > 2) {
    cerr << "d_m(t): Returned m = " << m << "\n";
    cerr << "m =\t\t\t";
    for (int k = mmode - mdown; k <= mmode + mup; ++k) cerr << k << "\t";
    cerr << "\nn_used =\t";
    for (int k = mmode - mdown; k <= mmode + mup; ++k)
      cerr << n_used[k] << "\t";

    cerr << "Coeffcount = " << coeffcount << endl;
  }
  return make_pair(m, coeffcount);
}

pair<int, int> WrightFisher::DrawAncestralProcessConditionalZero(
    double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of the Ancestral Process when conditioning
               /// the diffusion on non-absorption but started from boundary
{
  assert(t > 0.0);
  /// Pre-computing all necessary quantities
  int mup = -1, mdown = 0,
      mmode = static_cast<int>(round(GriffithsParas(t).first)), m = -1,
      eCindex_computed = 0, q1index = 0;
  bool m_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen), currsumU = 0.0, currsumL = 0.0, eCL = 0.0, eCU,
            q1L = 0.0, q1U = 0.0;
  double100 eCnew = (thetaP.empty()
                         ? 1.0
                         : static_cast<double100>((theta + 1.0)) *
                               exp(-t * static_cast<double100>(theta) * 0.5)),
            consForD = (thetaP.empty()) ? 3.0 : 5.0;
  map<int, double100> eAL, eAU, dL, dU;
  map<int, int> n_used, v_used;

  int threshold = ((thetaP.empty()) ? 2 : 1);

  while (!m_found) {
    if (mup > mdown &&
        mmode - mdown >
            threshold)  /// Checking whether down moves still allow - makes sure
                        /// m does not go below allowed values
    {
      ++mdown;
      m = mmode - mdown;
    } else  /// Otherwise move upwards, which is always allowed
    {
      ++mup;
      m = max(mmode + mup, threshold);
    }
    assert(m >= threshold);

    int n = m - 2, v = 0, n1 = -1;
    bool coefficients_converging = false;
    dU.insert(make_pair(m, 0.0));
    dL.insert(make_pair(m, 0.0));
    eAU.insert(make_pair(m, 0.0));
    eAL.insert(make_pair(m, 0.0));
    /// Computing quantities beyond which upper and lower bounds converge
    /// theoretically
    pair<vector<int>, double100> Ct;
    Ct.second = t;
    int F1 = computeC(m, Ct),
        F2 = static_cast<int>(
            ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
        F3 = static_cast<int>(ceil((log(consForD) / t) - (theta / 2.0)));
    while ((theta + static_cast<double>(2 * F2 + 1)) *
               exp(-(static_cast<double>(2 * F2) + theta) *
                   static_cast<double>(t) / 2.0) >=
           1 - o.eps)
      ++F2;
    int F = max(F1, F2) + 1;

    while (n < F || v < F3 || !coefficients_converging ||
           ((dU[m]) > 1.0 &&
            dL[m] < 0.0))  /// Checking if upper and lower sums converging
    {
      n += 2;
      n1 += 2;
      v++;
      double100 newcoefficientU =
          (((n - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double100>(n, m) +
              static_cast<double100>(-n * (n + theta - 1) * t / 2.0));
      double100 newcoefficientL =
          (((n + 1 - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double100>(n + 1, m) +
              static_cast<double100>(-(n + 1) * ((n + 1) + theta - 1) * t /
                                     2.0));

      if ((thetaP.empty()) &&
          (2 * (n1 + 1) >
           q1index))  /// If mutation vector is empty, we need to subtract q_1
                      /// from normalising constant (because m can only be >= 2)
      {
        double100 newcoefficientUq1 =
                      (((n1 - 1) % 2 != 0) ? -1.0 : 1.0) *
                      exp(Getlogakm<double100>(n1, 1) +
                          static_cast<double100>(-n1 * (n1 + theta - 1) * t /
                                                 2.0)),
                  newcoefficientLq1 =
                      ((n1 % 2 != 0) ? -1.0 : 1.0) *
                      exp(Getlogakm<double100>(n1 + 1, 1) +
                          static_cast<double100>(
                              -(n1 + 1) * ((n1 + 1) + theta - 1) * t / 2.0));
        q1U = q1L + newcoefficientUq1;
        q1L = q1U + newcoefficientLq1;
        q1index += 2;
      }

      if ((v + 1) > eCindex_computed)  /// Computing the denominator by using
                                       /// explicit expression for falling
                                       /// factorial moments of ancestral
                                       /// process and subtracting q_1
      {
        eCL += eCnew;
        eCnew = exp(static_cast<double100>(-(v + 1) * (v + theta) * t / 2.0) +
                    log(static_cast<double100>(2 * (v + 1) + theta - 1)));
        eCU = eCL + (eCnew / (1.0 - consForD * exp(-static_cast<double100>(
                                                       v + 1 + (theta / 2.0)) *
                                                   t)));
        eCindex_computed += 1;
      }

      if (o.debug > 2) {
        cerr << setprecision(0) << n << " newcoeffU " << setprecision(50)
             << newcoefficientU << "\n"
             << n + 1 << " newcoeffL " << newcoefficientL << endl;
      }

      eAU[m] = eAL[m] + newcoefficientU;
      eAL[m] = eAU[m] + newcoefficientL;

      coefficients_converging = (newcoefficientU + newcoefficientL >= 0.0);

      if (o.debug > 2) {
        cerr << "\nk   ";
        for (int k = mmode - mdown; k <= mmode + mup; ++k) cerr << k << " ";
        cerr << "\ndL ";
        for (int k = mmode - mdown; k <= mmode + mup; ++k) cerr << dL[k] << " ";
        cerr << "\ndU ";
        for (int k = mmode - mdown; k <= mmode + mup; ++k) cerr << dU[k] << " ";
        cerr << endl;
      }
    }

    if (eAL[m] < 0.0 ||
        eCL < 0.0)  /// If we fall out of bounds, use Gaussian approximations
    {
      cerr << "Numerical error detected: eAL[" << m << "] = " << eAL[m]
           << ", eAU[" << m << "] = " << eAU[m] << ", eCL = " << eCL
           << ", eCU = " << eCU
           << ". Resorting to G1984 approximation (t = " << t << ") ..."
           << endl;
      return make_pair(DrawAncestralProcessG1984(t, gen), -1);
    }

    dU[m] =
        ((!(eAU[m] > 0.0) || eCL < 0.0 || eCL < q1U)
             ? static_cast<double100>(nan(""))
             : exp(log(static_cast<double100>(m)) + log(eAU[m])) / (eCL - q1U));
    dL[m] =
        ((!(eAL[m] > 0.0) || eCU < q1L)
             ? static_cast<double100>(0.0)
             : exp(log(static_cast<double100>(m)) + log(eAL[m])) / (eCU - q1L));

    currsumL += dL[m];
    currsumU += dU[m];

    if (currsumL > currsumU) {
      cerr << "Error: currsumL = " << currsumL << " > " << currsumU
           << " = currsumU." << endl;
      exit(1);
    }

    n_used.insert(make_pair(m, n + 1));
    v_used.insert(make_pair(m, v));

    bool decision_on_m_made = (currsumL > u || currsumU < u);

    while (!decision_on_m_made)  /// Refine upper and lower bounds
    {
      double100 currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      for (int k = max(mmode, threshold) - mdown;
           k <= max(mmode, threshold) + mup; ++k) {
        ++n_used[k];
        ++v_used[k];  /// Adding on more contributions to sharpen bounds

        if (thetaP.empty() && ((n_used[k] + 1) > q1index)) {
          double100 newcoefficientUq1 =
                        (((n_used[k] - 1) % 2 != 0) ? -1.0 : 1.0) *
                        exp(Getlogakm<double100>(n_used[k], 1) +
                            static_cast<double100>(-n_used[k] *
                                                   (n_used[k] + theta - 1) * t /
                                                   2.0)),
                    newcoefficientLq1 =
                        (((n_used[k]) % 2 != 0) ? -1.0 : 1.0) *
                        exp(Getlogakm<double100>(n_used[k] + 1, 1) +
                            static_cast<double100>(
                                -(n_used[k] + 1) *
                                ((n_used[k] + 1) + theta - 1) * t / 2.0));
          q1U = q1L + newcoefficientUq1;
          q1L = q1U + newcoefficientLq1;
        }

        if ((v_used[k] + 1) > eCindex_computed) {
          eCL += eCnew;
          eCnew =
              exp(static_cast<double100>(-(v_used[k] + 2) *
                                         (v_used[k] + theta + 1) * t / 2.0) +
                  log(static_cast<double100>(2 * (v_used[k] + 2) + theta - 1)));
          eCU = eCL + (eCnew /
                       (1.0 - consForD * exp(-static_cast<double100>(
                                                 v_used[k] + 2 + theta / 2.0) *
                                             t)));
          eCindex_computed += 1;
        }

        eAU[k] += (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                  exp(Getlogakm<double100>(n_used[k], k) +
                      static_cast<double100>(
                          -n_used[k] * (n_used[k] + theta - 1) * t / 2.0));
        ++n_used[k];
        eAL[k] += (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                  exp(Getlogakm<double100>(n_used[k], k) +
                      static_cast<double100>(
                          -n_used[k] * (n_used[k] + theta - 1) * t / 2.0));

        dU[k] = ((eCL < 0.0 || !(eAU[k] > 0.0) || eCL < q1U)
                     ? static_cast<double100>(nan(""))
                     : exp(log(static_cast<double100>(k)) + log(eAU[k])) /
                           (eCL - q1U));
        dL[k] = ((eAL[k] < 0.0 || eCU < q1L)
                     ? static_cast<double100>(0.0)
                     : exp(log(static_cast<double100>(k)) + log(eAL[k])) /
                           (eCU - q1L));

        currsumL += dL[k];
        currsumU += dU[k];

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU." << endl;
          exit(1);
        }
      }

      if (currsumLold > currsumL) {
        // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
        std::cout << "Error: currsumLold = " << currsumLold << " > " << currsumL
                  << " = currsumL." << endl;
        // exit(1);
      }
      if (currsumUold < currsumU) {
        // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
        std::cout << "Error: currsumUold = " << currsumUold << " < " << currsumU
                  << " = currsumU." << endl;
        // exit(1);
      }

      decision_on_m_made = (currsumL > u || currsumU < u);
    }

    m_found = (currsumL > u);
  }

  int coeffcount = 0;
  for (int k = max(mmode, threshold) - mdown; k <= max(mmode, threshold) + mup;
       ++k)
    coeffcount += (n_used[k] - k + 1);

  if (o.debug > 2) {
    cerr << "d_m(t): Returned m = " << m << "\n";
    cerr << "m =\t\t\t";
    for (int k = max(mmode, threshold) - mdown;
         k <= max(mmode, threshold) + mup; ++k)
      cerr << k << "\t";
    cerr << "\nn_used =\t";
    for (int k = max(mmode, threshold) - mdown;
         k <= max(mmode, threshold) + mup; ++k)
      cerr << n_used[k] << "\t";

    cerr << "Coeffcount = " << coeffcount << endl;
  }
  return make_pair(m, coeffcount);
}

pair<pair<int, int>, int> WrightFisher::DrawAncestralProcessConditionalInterior(
    double100 t, double100 x, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of the Ancestral Process when conditioning
               /// the diffusion on non-absorption but started from inside (0,1)
{
  assert((x > 0.0) && (x <= 1.0) && (t > 0.0));
  /// Pre-computing all necessary quantities
  int mup = -1, mdown = 0,
      mmode = static_cast<int>(round(GriffithsParas(t).first)), m = -1,
      eCindex_computed = 0, lstore;
  bool ml_found = false;

  /// Setting up all the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU = Getd2(eCvec, 0, x, t);
  map<int, double100> eAL, eAU, dL, dU;

  map<int, int> n_used, v_used;

  int threshold = ((thetaP.empty()) ? 2 : 1);

  while (!ml_found) {
    if (mup > mdown &&
        mmode - mdown >
            threshold)  /// Checking whether down moves are still possible
    {
      ++mdown;
      m = mmode - mdown;
    } else  /// If not, upward moves are always allowed
    {
      ++mup;
      m = max(mmode + mup, threshold);
    }

    int n = m - 2;
    bool coefficients_converging = false;
    int v = -1;
    dU.insert(make_pair(m, 0.0));
    dL.insert(make_pair(m, 0.0));
    eAU.insert(make_pair(m, 0.0));
    eAL.insert(make_pair(m, 0.0));
    /// Compute F ensuring that theoretically upper and lower bounds converge
    pair<vector<int>, double100> Ct;
    Ct.second = t;
    int F1 = computeC(m, Ct),
        F2 = static_cast<int>(
            ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
        F3 = static_cast<int>(ceil((log(3.0) / t) - (theta / 2.0)));
    while ((theta + static_cast<double>(2 * F2 + 1)) *
               exp(-(static_cast<double>(2 * F2) + theta) *
                   static_cast<double>(t) / 2.0) >=
           1 - o.eps)
      ++F2;
    int F = max(F1, max(F2, F3)) + 1;

    while (2 * v < F || n < F || !coefficients_converging ||
           ((dU[m]) > 1.0 &&
            dL[m] < 0.0))  /// Checking if upper and lower sums converging
    {
      n += 2;
      v++;  /// Adding on new contributions for numerator
      double100 newcoefficientU =
                    (((n - m) % 2 != 0) ? -1.0 : 1.0) *
                    exp(Getlogakm<double100>(n, m) +
                        static_cast<double100>(-n * (n + theta - 1) * t / 2.0)),
                newcoefficientL =
                    (((n + 1 - m) % 2 != 0) ? -1.0 : 1.0) *
                    exp(Getlogakm<double100>(n + 1, m) +
                        static_cast<double100>(
                            -(n + 1) * ((n + 1) + theta - 1) * t / 2.0));

      if (2 * (v + 1) > eCindex_computed)  /// Computing denominator
      {
        assert(2 * v == eCindex_computed);
        eCL = eCU - Getd2(eCvec, 2 * v + 1, x, t);
        eCU = eCL + Getd2(eCvec, 2 * v + 2, x, t);
        eCindex_computed += 2;
      }

      if (o.debug > 2) {
        cerr << setprecision(0) << n << " newcoeffU " << setprecision(50)
             << newcoefficientU << "\n"
             << n + 1 << " newcoeffL " << newcoefficientL << endl;
      }

      eAU[m] = eAL[m] + newcoefficientU;
      eAL[m] = eAU[m] + newcoefficientL;

      coefficients_converging = (newcoefficientU + newcoefficientL >= 0.0);

      if (o.debug > 2) {
        cerr << "\nk   ";
        for (int k = max(mmode, threshold) - mdown;
             k <= max(mmode, threshold) + mup; ++k)
          cerr << k << " ";
        cerr << "\ndL ";
        for (int k = max(mmode, threshold) - mdown;
             k <= max(mmode, threshold) + mup; ++k)
          cerr << dL[k] << " ";
        cerr << "\ndU ";
        for (int k = max(mmode, threshold) - mdown;
             k <= max(mmode, threshold) + mup; ++k)
          cerr << dU[k] << " ";
        cerr << endl;
      }
    }

    if (eAL[m] < 0.0 ||
        eCL < 0.0)  /// If we fall out of bounds, use Gaussian approximations
    {
      cerr << "Numerical error detected: eAL[" << m << "] = " << eAL[m]
           << ", eAU[" << m << "] = " << eAU[m] << ", eCL = " << eCL
           << ", eCU = " << eCU
           << ". Resorting to G1984 approximation (t = " << t << ") ..."
           << endl;
      m = DrawAncestralProcessG1984(t, gen);
      vector<double100> binomialBins;
      boost::math::binomial_distribution<double100> BINOMIAL(m, x);
      int mlimit = (threshold == 2) ? m - 1 : m;
      for (int k = 1; k <= mlimit; ++k) {
        binomialBins.push_back(pdf(BINOMIAL, k));
      }
      boost::random::discrete_distribution<> CATEGORICAL(binomialBins.begin(),
                                                         binomialBins.end());
      int l = CATEGORICAL(gen) + 1;

      return make_pair(make_pair(m, l), -1);
    }

    dU[m] = (eCL < 0.0 ? static_cast<double100>(nan("")) : eAU[m] / eCL);
    dL[m] = (eAL[m] < 0.0 ? static_cast<double100>(0.0) : eAL[m] / eCU);

    int llimit = (threshold == 2) ? m - 1 : m;
    boost::math::binomial_distribution<double100> BIN(m, x);
    for (int l = 1; l <= llimit && !ml_found;
         ++l)  /// Adding on the corresponding binomial terms (minus the edge
               /// cases i.e. l=0 and potentially l=m (depending on theta))
    {
      double100 currsummaddon =
          pdf(BIN, l);  //( ( !(x < 1.0) || !(x > 0.0) ) ? 1.0 : pdf(BIN,l) );
      currsumL += currsummaddon * dL[m];
      currsumU += currsummaddon * dU[m];

      bool decision_on_ml_made = (currsumL > u || currsumU < u);

      if (currsumL > currsumU) {
        cerr << "Error: currsumL = " << currsumL << " > " << currsumU
             << " = currsumU (n,m,l) = (" << n << "," << m << "," << l << ")."
             << endl;
        exit(1);
      }

      while (!decision_on_ml_made)  /// Refine upper and lower bounds
      {
        double100 currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const map<int, double100> dUold(dU), dLold(dL);
        for (int k = max(mmode, threshold) - mdown;
             k <= max(mmode, threshold) + mup; ++k) {
          ++v_used[k];
          ++n_used[k];  /// Adding on new contributions to sharpen bounds
          double100 newcoefficientU =
                        exp(Getlogakm<double100>(k + 2 * n_used[k], k) +
                            static_cast<double100>(
                                -(k + 2 * v_used[k]) *
                                (k + 2 * n_used[k] + theta - 1) * t / 2.0)),
                    newcoefficientL =
                        exp(Getlogakm<double100>(k + 2 * n_used[k] + 1, k) +
                            static_cast<double100>(
                                -(k + 2 * n_used[k] + 1) *
                                (k + 2 * n_used[k] + 1 + theta - 1) * t / 2.0));

          eAU[k] = eAL[k] + newcoefficientU;
          eAL[k] = eAU[k] - newcoefficientL;

          if (2 * v_used[k] + 2 > eCindex_computed) {
            assert(2 * v_used[k] == eCindex_computed);
            eCL = eCU - Getd2(eCvec, 2 * v_used[k] + 1, x, t);
            eCU = eCL + Getd2(eCvec, 2 * v_used[k] + 2, x, t);
            eCindex_computed += 2;
          }

          dU[k] = (eCL < 0.0 ? static_cast<double100>(nan("")) : eAU[k] / eCL);
          dL[k] = (eAL[k] < 0.0 ? static_cast<double100>(0.0) : eAL[k] / eCU);

          int l1limit = (threshold == 2) ? k - 1 : k;
          boost::math::binomial_distribution<double100> BIN2(k, x);
          for (int l1 = 1; l1 <= l1limit; ++l1) {
            double100 addon =
                ((!(x < 1.0) || !(x > 0.0)) ? 1.0 : pdf(BIN2, l1));
            currsumL += addon * dL[k];
            currsumU += addon * dU[k];
          }

          if (currsumL > currsumU) {
            cerr << "Error: currsumL = " << currsumL << " > " << currsumU
                 << " = currsumU (n,m,l) = (" << n << "," << m << "," << l
                 << ")." << endl;
            exit(1);
          }
        }

        if (currsumLold > currsumL) {
          // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,l) = (" << n << "," << m
                    << "," << l << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
          std::cout << "Error: currsumUold = " << currsumUold << " < "
                    << currsumU << " = currsumU (n,m,l) = (" << n << "," << m
                    << "," << l << ")." << endl;
          // exit(1);
        }

        decision_on_ml_made = (currsumL > u || currsumU < u);
      }

      ml_found = (currsumL > u);
      { lstore = l; }
    }
  }

  int coeffcount = 0;
  for (int k = max(mmode, threshold) - mdown; k <= max(mmode, threshold) + mup;
       ++k)
    coeffcount += (n_used[k] - k + 1);

  if (o.debug > 2) {
    cerr << "d_m(t): Returned m = " << m << "\n";
    cerr << "m =\t\t\t";
    for (int k = max(mmode, threshold) - mdown;
         k <= max(mmode, threshold) + mup; ++k)
      cerr << k << "\t";
    cerr << "\nn_used =\t";
    for (int k = max(mmode, threshold) - mdown;
         k <= max(mmode, threshold) + mup; ++k)
      cerr << n_used[k] << "\t";

    cerr << "Coeffcount = " << coeffcount << endl;
  }
  return make_pair(make_pair(m, lstore), coeffcount);
}

int WrightFisher::DrawAncestralProcessG1984(
    double100 t, boost::random::mt19937
                     &gen)  /// Draws from the law of the Ancestral Process when
                            /// time increment falls below threshold
{
  assert(t > 0.0);
  int threshold = (thetaP.empty()                               ? 2
                   : (!(thetaP[0] > 0.0) || !(thetaP[1] > 0.0)) ? 1
                                                                : 0);
  /// Pre-compute all necessary quantities
  double mean = GriffithsParas(t).first, v = GriffithsParas(t).second;
  assert(v > 0.0);
  boost::random::normal_distribution<double100> NORMAL(mean, sqrt(v));
  return max(static_cast<int>(round(NORMAL(gen))),
             threshold);  // Return discretised Gaussian approximation
}

int WrightFisher::DrawSizebiasedAncestralProcess(
    int d, double100 t,
    boost::random::mt19937
        &gen)  /// Draws from the size-biased distribution of the Ancestral
               /// process when the time increment falls below threshold
{
  assert((d > 0) && (t > 0.0));
  /// Pre-compute all necessary quantities
  double mean = GriffithsParas(t).first, v = GriffithsParas(t).second;
  assert(v > 0.0);
  /// Normalisation obtained from integral_{0}^{inf} of x*Gaussian pdf(x) dx
  double normalisation =
      0.5 * sqrt(v) *
      (exp(-pow(mean, 2.0) * 0.5 * (1.0 / v)) *
           sqrt(2.0 / boost::math::constants::pi<double>()) +
       (mean / sqrt(v)) * (1.0 + boost::math::erf(mean / (sqrt(2.0)))));

  bool decision = false;
  int flip = 1, mlimit = thetaP.empty() ? 2 : 1,
      mode = static_cast<int>(round(mean)), m = mode, jp = 0,
      jm = 0;  /// Starting m from mode to speed up routine
  double mlim = static_cast<double>(mlimit) + 0.5, summqm = 0.0;

  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);

  while (!decision) {
    if (m ==
        mlimit)  /// Computing correct border (i.e. m = 1 or 2) contribution
    {
      summqm += exp(static_cast<double>(d) * log(static_cast<double>(mlim))) *
                0.5 *
                erfc(-(static_cast<double>(mlim) - mean) / (sqrt(2.0 * v)));
    } else  /// Add on discretised Gaussian contributions
    {
      summqm +=
          exp(static_cast<double>(d) * log(static_cast<double>(m))) * 0.5 *
          (erfc(-(static_cast<double>(m) + 0.5 - mean) / (sqrt(2.0 * v))) -
           erfc(-(static_cast<double>(m) - 0.5 - mean) / (sqrt(2.0 * v))));
    }
    if (summqm / normalisation > u)  /// Return last m
    {
      decision = true;
    } else {
      if (flip == -1 &&
          (mode - jm > 2))  /// Mechanism to explore either side around mode
      {
        jm++;
        m = mode - jm;
      } else {
        jp++;
        m = mode + jp;
      }
      flip *= -1;
    }
  }
  assert(m >= mlimit);

  return m;
}

pair<double100, int> WrightFisher::DrawEndpoint(
    double100 x, double100 t1, double100 t2, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of a Wright-Fisher
                                  /// diffusion conditioned on non-absorption
{
  // std::cout << "x is " << x << ", t1 is " << t1 << ", t2 is " << t2
  //           << std::endl;
  assert((x >= 0.0) && (x <= 1.0) && (t1 < t2));
  double100 y;
  int m, coeffcount;
  pair<int, int> m_and_coeffcount;

  if (thetaP.empty())  /// No mutation - condition on avoiding either boundary
  {
    if (!(x >
          0.0))  /// Diffusion conditioned on non-absorption but started from 0
    {
      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount = DrawAncestralProcessConditionalZero(t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(1.0, 1.0),
          GAMMA2(static_cast<double>(m) - 1.0, 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    } else if (!(x < 1.0))  /// Diffusion conditioned on non-absorption but
                            /// started from 1
    {
      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount = DrawAncestralProcessConditionalZero(t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(static_cast<double>(m) - 1.0,
                                                 1.0),
          GAMMA2(1.0, 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    } else  /// Diffusion conditioned on non-absorption but started from the
            /// interior (0,1)
    {
      int l;

      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawAncestralProcessG1984(t2 - t1, gen);
        coeffcount = -1;
        vector<double100> binomialBins;
        boost::math::binomial_distribution<double100> BINOMIAL(
            m, static_cast<double>(x));
        for (int k = 1; k <= m - 1; ++k) {
          binomialBins.push_back(pdf(BINOMIAL, k));
        }
        boost::random::discrete_distribution<> CATEGORICAL(binomialBins.begin(),
                                                           binomialBins.end());
        l = CATEGORICAL(gen) + 1;
      } else  /// Otherwise use alternating series technique
      {
        pair<pair<int, int>, int> storer =
            DrawAncestralProcessConditionalInterior(t2 - t1, x, o, gen);
        m = (storer.first).first;
        l = (storer.first).second;
        coeffcount = storer.second;
      }

      boost::random::gamma_distribution<> GAMMA1(static_cast<double>(l), 1.0),
          GAMMA2(static_cast<double>(m - l), 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    }
  } else if ((!(thetaP[0] > 0.0)) ||
             (!(thetaP[1] > 0.0)))  /// One sided mutation - condition on
                                    /// avoiding one boundary only
  {
    if ((!(thetaP[0] > 0.0)) &&
        (!(x >
           0.0)))  /// Diffusion conditioned on avoiding 0 but started from 0
    {
      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount = DrawAncestralProcessConditionalZero(t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(1.0, 1.0),
          GAMMA2(thetaP[1] + static_cast<double>(m) - 1.0, 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    } else if ((!(thetaP[1] > 0.0)) &&
               (!(x < 1.0)))  /// Diffusion conditioned on avoiding 1 but
                              /// started from 1 - can translate to diffusion
                              /// with mutation entries swapped avoiding 0
                              /// started from 0 by symmetry
    {
      iter_swap(thetaP.begin(), thetaP.begin() + 1);  /// Swap mutation entries

      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount = DrawAncestralProcessConditionalZero(t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(1.0, 1.0),
          GAMMA2(thetaP[1] + static_cast<double>(m) - 1.0, 1.0);

      double100 y1 = -1.0;
      while (!(0.0 < y1 && y1 < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y1 = G1 / (G1 + G2);
      }

      y = 1.0 - y1;

      iter_swap(thetaP.begin(),
                thetaP.begin() + 1);  /// Swap back mutation entries
    } else if (!(thetaP[1] >
                 0.0))  /// Diffusion conditioned on avoiding 1 but started from
                        /// x in [0,1) - can translate to diffusion with
                        /// mutation entries swapped avoiding 0 started from 1-x
                        /// in interior (0,1] by symmetry
    {
      iter_swap(thetaP.begin(), thetaP.begin() + 1);  /// Swap mutation entries

      int l;

      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawAncestralProcessG1984(t2 - t1, gen);
        coeffcount = -1;
        vector<double100> binomialBins;
        boost::math::binomial_distribution<double100> BINOMIAL(m, 1.0 - x);
        for (int k = 1; k <= m; ++k) {
          binomialBins.push_back(pdf(BINOMIAL, k));
        }
        boost::random::discrete_distribution<> CATEGORICAL(binomialBins.begin(),
                                                           binomialBins.end());
        l = CATEGORICAL(gen) + 1;
      } else  /// Otherwise use alternating series technique
      {
        pair<pair<int, int>, int> storer =
            DrawAncestralProcessConditionalInterior(t2 - t1, 1.0 - x, o, gen);
        m = (storer.first).first;
        l = (storer.first).second;
        coeffcount = storer.second;
      }

      boost::random::gamma_distribution<> GAMMA1(
          thetaP[0] + static_cast<double>(l), 1.0),
          GAMMA2(thetaP[1] + static_cast<double>(m - l), 1.0);

      double100 y1 = -1.0;
      while (!(0.0 < y1 && y1 < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y1 = G1 / (G1 + G2);
      }
      y = 1.0 - y1;

      iter_swap(thetaP.begin(),
                thetaP.begin() + 1);  /// Swap mutation entries back
    } else  /// Diffusion conditioned on avoiding 0 but started from x in (0,1]
    {
      int l;

      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawAncestralProcessG1984(t2 - t1, gen);
        coeffcount = -1;
        vector<double100> binomialBins;
        boost::math::binomial_distribution<> BINOMIAL(m, x);
        for (int j = 1; j <= m; j++) {
          binomialBins.push_back(pdf(BINOMIAL, j));
        }
        boost::random::discrete_distribution<> CAT(binomialBins.begin(),
                                                   binomialBins.end());
        l = CAT(gen) + 1;
      } else  /// Otherwise use alternating series technique
      {
        pair<pair<int, int>, int> storer =
            DrawAncestralProcessConditionalInterior(t2 - t1, x, o, gen);
        m = (storer.first).first;
        l = (storer.first).second;
        coeffcount = storer.second;
      }

      boost::random::gamma_distribution<> GAMMA1(
          thetaP[0] + static_cast<double100>(l), 1.0),
          GAMMA2(thetaP[1] + static_cast<double100>(m - l), 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    }
  } else  /// Strictly positive entries in mutation vector, so unconditioned
          /// diffusion
  {
    if (t2 - t1 <= o.g1984threshold)  /// Time increment smaller than threshold,
                                      /// use Gaussian approximations
    {
      m = DrawAncestralProcessG1984(t2 - t1, gen);
      coeffcount = -1;
    } else  /// Otherwise use alternating series technique
    {
      m_and_coeffcount = DrawAncestralProcess(t2 - t1, o, gen);
      m = m_and_coeffcount.first;
      coeffcount = m_and_coeffcount.second;
    }

    boost::random::binomial_distribution<int> BINOMIAL(m,
                                                       static_cast<double>(x));
    int l = BINOMIAL(gen);

    boost::random::gamma_distribution<> GAMMA1(
        thetaP[0] + static_cast<double100>(l), 1.0),
        GAMMA2(thetaP[1] + static_cast<double100>(m - l), 1.0);

    y = -1.0;
    while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
    {
      double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }
  }

  return make_pair(y, coeffcount);
}

pair<double100, int> WrightFisher::DrawUnconditionedDiffusion(
    double100 x, double100 t, const Options &o,
    boost::random::mt19937 &
        gen)  /// Draws from the law of an unconditioned Wright-Fisher diffusion
{
  assert((x >= 0.0) && (x <= 1.0) && (t > 0.0));

  if (!(x > 0.0) &&
      (thetaP.empty() ||
       (!(thetaP[0] > 0.0))))  /// If we start at 0 and there is no
                               /// mutation, diffusion stays there
  {
    return make_pair(0.0, -1);
  } else if (!(x < 1.0) &&
             (thetaP.empty() ||
              (!(thetaP[1] >
                 0.0))))  /// Similarly if started from 1 with no mutation
  {
    return make_pair(1.0, -1);
  } else  /// Otherwise we first draw M ~ {q_m}
  {
    int m, coefcount;
    if (t <= o.g1984threshold) {
      m = DrawAncestralProcessG1984(t, gen);
      coefcount = -1;
    } else {
      pair<int, int> temp = DrawAncestralProcess(t, o, gen);
      m = temp.first;
      coefcount = temp.second;
    }

    double100 para1, para2;
    boost::random::binomial_distribution<> BIN(m, x);  /// Draw L ~ Bin(m,x)
    int l = BIN(gen);

    if (thetaP.empty())  /// Sort out cases based on what mutation parameters is
    {
      if (l == 0) {
        return make_pair(0.0, coefcount);
      } else if (l == m) {
        return make_pair(1.0, coefcount);
      } else {
        para1 = static_cast<double100>(l);
        para2 = static_cast<double100>(m - l);
      }
    } else if (!(thetaP[0] > 0.0)) {
      if (l == 0) {
        return make_pair(0.0, coefcount);
      } else {
        para1 = static_cast<double100>(l);
        para2 = static_cast<double100>(thetaP[1] + m - l);
      }
    } else if (!(thetaP[1] > 0.0)) {
      if (l == m) {
        return make_pair(1.0, coefcount);
      } else {
        para1 = static_cast<double100>(thetaP[0] + l);
        para2 = static_cast<double100>(m - l);
      }
    } else {
      para1 = static_cast<double100>(thetaP[0] + l);
      para2 = static_cast<double100>(thetaP[0] + m - l);
    }

    boost::random::gamma_distribution<double100> GAMMA1(para1, 1.0),
        GAMMA2(para2, 1.0);
    double100 y = -1.0;
    while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
    {
      double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }

    return make_pair(y, coefcount);
  }
}

/// DIFFUSION SIMULATION - NON-NEUTRAL PATHS

vector<vector<double100>> WrightFisher::NonNeutralDraw(
    double100 x, double100 t1, double100 t2, bool Absorption, const Options &o,
    boost::random::mt19937 &gen)  /// Draws of paths from non-neutral WF
                                  /// diffusion started from x at time t1
{
  assert((x >= 0.0) && (x <= 1.0) && (t1 < t2));
  bool accept = false;
  vector<double100> paras{phiMin, phiMax, phiMax - phiMin};
  vector<vector<double100>> ptr;
  double100 kapmean = paras[2] * (t2 - t1);  /// Rate of Poisson point process

  boost::random::poisson_distribution<int> kap(static_cast<double>(kapmean));

  boost::random::uniform_01<double100> unift, unifm,
      unifU;  /// Set up uniform generators for points over [t1,t2] *
  /// [0,phiMax-phiMin] * [0,1]
  int rcount = 0;
  while (!accept)  /// Until all skeleton points get accepted, keep going
  {
    int kappa = kap(gen);  /// Generate kappa ~ Pois
    double100 u = unifU(gen);
    double100 small_offset = 1.0e-14;

    vector<double100> path, times(kappa), marks(kappa), rejcount;
    auto gent = [&t1, &t2, &unift, &small_offset,
                 &gen]()  /// time stamps ~ Unif [t1,t2]
    { return (t1 + small_offset + ((t2 - t1) * unift(gen))); };
    auto genm = [&paras, &unifm, &gen]()  /// marks ~ Unif [0,phiMax-phiMin]
    { return (paras[2] * unifm(gen)); };
    std::generate(begin(times), end(times), gent);
    std::generate(begin(marks), end(marks), genm);
    sortVectorsAscending(times, times,
                         marks);  /// Sort vectors according to timestamps

    times.push_back(t2);
    marks.push_back(
        u);  /// Add on end point - to be checked differently to skeleton points

    for (vector<double100>::iterator itt = times.begin(), itm = marks.begin();
         itt != times.end(); itt++, itm++) {
      if (kappa == 0)  /// No skeleton points -> check end point directly
      {
        if (Absorption) {
          path.push_back(DrawUnconditionedDiffusion(x, t2 - t1, o, gen)
                             .first);  /// Generate endpoint
        } else {
          path.push_back(
              DrawEndpoint(x, t1, t2, o, gen).first);  /// Generate endpoint
        }

        if (exp(Atilde(path.back()) - Atildeplus()) <
            *itm)  /// Test if generated point is good; if not generate new
                   /// Poisson point process
        {
          rcount++;
          break;
        }

        rejcount.push_back(rcount);
        ptr.push_back(path);
        ptr.push_back(times);
        ptr.push_back(rejcount);
        accept = true;
      } else  /// kappa > 0 - generate skeleton points and endpoint
      {
        if (itt == times.begin()) {
          if (Absorption) {
            path.push_back(
                DrawUnconditionedDiffusion(x, (*itt) - t1, o, gen)
                    .first);  /// Generate skeleton points sequentially
          } else {
            path.push_back(
                DrawEndpoint(x, t1, *itt, o, gen)
                    .first);  /// Generate skeleton points sequentially
          }

          if (Phitilde(path.back()) - paras[0] >
              *itm)  /// Test generated point is OK, otherwise can stop and
                     /// generate a new Poisson point process
          {
            rcount++;
            break;
          }

        } else if (*itt != t2)  /// Separate first time stamp and rest to ensure
                                /// correct sequential sampling
        {
          if (Absorption) {
            path.push_back(DrawUnconditionedDiffusion(
                               path.back(), (*itt) - *(itt - 1), o, gen)
                               .first);
          } else {
            path.push_back(
                DrawEndpoint(path.back(), *(itt - 1), *itt, o, gen).first);
          }

          if (Phitilde(path.back()) - paras[0] > *itm) {
            rcount++;
            break;
          }

        } else  /// Endpoint draw
        {
          if (Absorption) {
            path.push_back(DrawUnconditionedDiffusion(
                               path.back(), (*itt) - *(itt - 1), o, gen)
                               .first);
          } else {
            path.push_back(
                DrawEndpoint(path.back(), *(itt - 1), *itt, o, gen).first);
          }

          if (exp(Atilde(path.back()) - Atildeplus()) <
              *itm)  /// Check corresponding endpoint condition
          {
            rcount++;
            break;
          }

          rejcount.push_back(rcount);
          ptr.push_back(path);
          ptr.push_back(times);
          ptr.push_back(rejcount);
          accept = true;
        }
      }
    }
  }

  return ptr;
}

pair<double100, int> WrightFisher::NonNeutralDrawEndpoint(
    double100 x, double100 t1, double100 t2, bool Absorption, const Options &o,
    boost::random::mt19937
        &gen)  /// Invoke NonNeutralDraw to generate a whole path, but retains
               /// only the endpoint at time t2
{
  vector<vector<double100>> ptr = NonNeutralDraw(x, t1, t2, Absorption, o, gen);
  return make_pair(ptr[0].back(), -1);
}

/// BRIDGE SIMULATION - NEUTRAL PATHS

vector<int> WrightFisher::DrawBridgePMFUnconditional(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting within (0,1), ending
               /// at some boundary point (absorption can happen at any time)
               /// with time increments being large enough
{
  assert((x > 0.0) && (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (s > 0.0) &&
         (t > 0.0));
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkl_found = false;

  /// Setting up necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(GetdBridgeUnconditional(eCvec, 0, x, z, t));
  vector<int> v_used;
  double100 Amkl;
  int n = -1, Fmkl = 0, eCindex_computed = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double100> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkl_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    eBL.push_back(0.0);
    eBU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(m, Cs), F2 = computeC(k, Cts);
    if (o.debug > 2)
      cerr << "(F1,F2,F3,F4,F5) = (" << F1 << "," << F2 << "," << F3 << ","
           << F4 << "," << F5 << ")" << endl;
    Fmkl = max(max(F1, F2), F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmkl << endl;
    }

    while (2 * v < Fmkl || eAU[n] > 1.0 || eBU[n] > 1.0) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(m + 2 * v, m) +
              static_cast<double100>(-(m + 2 * v) * (m + 2 * v + theta - 1) *
                                     s / 2.0));
      double100 newcoefficientL =
          exp(Getlogakm<double100>(m + 2 * v + 1, m) +
              static_cast<double100>(-(m + 2 * v + 1) *
                                     (m + 2 * v + 1 + theta - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double100>(k + 2 * v, k) +
              static_cast<double100>(-(k + 2 * v) * (k + 2 * v + theta - 1) *
                                     (t - s) / 2.0));
      newcoefficientL = exp(Getlogakm<double100>(k + 2 * v + 1, k) +
                            static_cast<double100>(-(k + 2 * v + 1) *
                                                   (k + 2 * v + 1 + theta - 1) *
                                                   (t - s) / 2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeUnconditional(eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeUnconditional(eCvec, 2 * v + 2, x, z, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] && eBU[n] == eBL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmkj.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmkj = " << Fmkl
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eBU[n] - eBL[n] << ", " << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFUnconditionalApprox(x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[n]) + log(eBL[n]) - log(eCU)));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    /// Add on l contributions
    int llower =
            (((thetaP.empty() || !(thetaP[0] > 0.0)) && !(z > 0.0)) ? 0 : 1),
        lupper =
            (((thetaP.empty() || !(thetaP[1] > 0.0)) && !(z < 1.0)) ? m
                                                                    : m - 1);
    for (int l = llower; l <= lupper && !mkl_found; ++l) {
      Amkl = computeAUnconditional(m, k, l, x, z);
      if (o.debug > 2)
        cerr << "Adding to currsums with A(n,m,k,l) = A(" << n << "," << m
             << "," << k << "," << l << ") = " << Amkl << endl;

      if (Amkl * dL[n] > 1.0 || Amkl * dU[n] < 0.0 || eAU[n] < 0.0 ||
          eBU[n] < 0.0 || eAL[n] > 1.0 || eBL[n] > 1.0 || eCU < 0.0) {
        cerr << "Numerical error detected: (m,k,l) = " << m << "," << k << ","
             << l << "), Amkj = " << Amkl << ", dL[" << n << "] = " << dL[n]
             << ", dU[" << n << "] = " << dU[n] << ", eAU[" << n
             << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
             << ",  eBU[" << n << "] = " << eBU[n] << ", eBL[" << n
             << "] = " << eBL[n] << ", eCU = " << eCU << ", eCL = " << eCL
             << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
             << "," << z << "," << s << "," << t << ") ..." << endl;
        return DrawBridgePMFUnconditionalApprox(x, z, s, t, o, gen);
      }

      currsumL += Amkl * dL[n];
      currsumU += Amkl * dU[n];
      bool decision_on_mkl_made = (currsumL > u || currsumU < u);

      if (currsumL > currsumU) {
        cerr << "Error: currsumL = " << currsumL << " > " << currsumU
             << " = currsumU (n,m,k,l) = (" << n << "," << m << "," << k << ","
             << l << ")." << endl;
        exit(1);
      }

      while (!decision_on_mkl_made)  /// Refine upper and lower bounds
      {
        double100 currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double100> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
          ++v_used[i];
          double100 newcoefficientU =
              exp(Getlogakm<double100>(mi + 2 * v_used[i], mi) +
                  static_cast<double100>(-(mi + 2 * v_used[i]) *
                                         (mi + 2 * v_used[i] + theta - 1) * s /
                                         2.0));
          double100 newcoefficientL =
              exp(Getlogakm<double100>(mi + 2 * v_used[i] + 1, mi) +
                  static_cast<double100>(-(mi + 2 * v_used[i] + 1) *
                                         (mi + 2 * v_used[i] + 1 + theta - 1) *
                                         s / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          newcoefficientU =
              exp(Getlogakm<double100>(ki + 2 * v_used[i], ki) +
                  static_cast<double100>(-(ki + 2 * v_used[i]) *
                                         (ki + 2 * v_used[i] + theta - 1) *
                                         (t - s) / 2.0));
          newcoefficientL =
              exp(Getlogakm<double100>(ki + 2 * v_used[i] + 1, ki) +
                  static_cast<double100>(-(ki + 2 * v_used[i] + 1) *
                                         (ki + 2 * v_used[i] + 1 + theta - 1) *
                                         (t - s) / 2.0));

          eBU[i] = eBL[i] + newcoefficientU;
          eBL[i] = eBU[i] - newcoefficientL;

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU -
                  GetdBridgeUnconditional(eCvec, 2 * v_used[i] + 1, x, z, t);
            eCU = eCL +
                  GetdBridgeUnconditional(eCvec, 2 * v_used[i] + 2, x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                             : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
          dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                       ? static_cast<double100>(0.0)
                       : exp(log(eAL[i]) + log(eBL[i]) - log(eCU)));

          if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
              static_cast<double>(eAU[i]) < 0.0 ||
              static_cast<double>(eBU[i]) < 0.0 ||
              static_cast<double>(eAL[i]) > 1.0 ||
              static_cast<double>(eBL[i]) > 1.0 ||
              static_cast<double>(eCU) < 0.0) {
            cerr << "Numerical error detected: (m,k,l) = " << mi << "," << ki
                 << ", *, *), dL[" << i << "] = " << dL[i] << ", dU[" << i
                 << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
                 << ", eAL[" << i << "] = " << eAL[i] << ",  eBU[" << i
                 << "] = " << eBU[i] << ", eBL[" << i << "] = " << eBL[i]
                 << ", eCU = " << eCU << ", eCL = " << eCL
                 << ". Resorting to G1984-style approximation (x,z,s,t) = ("
                 << x << "," << z << "," << s << "," << t << ") ..." << endl;
            return DrawBridgePMFUnconditionalApprox(x, z, s, t, o, gen);
          }

          int lilower =
                  (((thetaP.empty() || !(thetaP[0] > 0.0)) && !(z > 0.0)) ? 0
                                                                          : 1),
              liupper = (((thetaP.empty() || !(thetaP[1] > 0.0)) && !(z < 1.0))
                             ? mi
                             : mi - 1);
          for (int l2 = lilower; l2 <= liupper; ++l2) {
            Amkl = computeAUnconditional(mi, ki, l2, x, z);
            currsumL += Amkl * dL[i];
            currsumU += Amkl * dU[i];
            if (o.debug > 3)
              cerr << "Recomputing currsums with A(n,m,k,l) = A(" << i << ","
                   << mi << "," << ki << "," << l2 << ") = " << Amkl << endl;
          }
        }

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k,l) = (" << n << "," << m << "," << k
               << "," << l << ")." << endl;
          exit(1);
        }

        if (o.debug > 2) {
          cerr << "\ndL ";
          printVec(dL, cerr);
          cerr << "\ndU ";
          printVec(dU, cerr);
        }

        if (currsumLold > currsumL) {
          // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,l) = (" << n << "," << m
                    << "," << k << "," << l << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
          std::cout << "Error: currsumUold = " << currsumUold << " < "
                    << currsumU << " = currsumU (n,m,k,l) = (" << n << "," << m
                    << "," << k << "," << l << ")." << endl;
          // exit(1);
        }

        decision_on_mkl_made = (currsumL > u || currsumU < u);
      }

      mkl_found = (currsumL > u);

      if (mkl_found) {
        if (!(z > 0.0)) {
          curr_mk[n][2] = l;  /// Sets j = 0
          curr_mk[n][3] = 0;
        } else {
          curr_mk[n][2] = l;  /// Sets j = k
          curr_mk[n][3] = k;
        }
      }
    }
  }

  int coeffcount = 0;
  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,j: Returned (m,k,j) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << "," << curr_mk[n][3] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFUnconditionalOneQApprox(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting within (0,1), ending
               /// at some boundary point (absorption can happen at any time)
               /// with one time increment falling below the threshold
{
  assert((x > 0.0) && (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (t > s) &&
         (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double100 sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkl_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec, eAL, eAU, dL, dU, currsumStore;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(GetdBridgeUnconditional(eCvec, 0, x, z, t)), qApprox = 0.0,
            runningMax = 1.0e-300;
  vector<int> v_used;
  double100 Amkl;
  int n = -1, Fmkl = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to be
  /// needed
  double mmode = GriffithsParas(s).first, vm = GriffithsParas(s).second,
         kmode = GriffithsParas(t - s).first, vk = GriffithsParas(t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double100> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkl_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    int mork = (ind1 ? m : k), morkapprox = (mork == m ? k : m);
    if (m <= mlimU && m >= mlimL && k <= klimU && k >= klimL) counter++;
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(mork, Csorts);
    if (o.debug > 2)
      cerr << "(F1,F3,F4,F5) = (" << F1 << "," << F3 << "," << F4 << "," << F5
           << ")" << endl;
    Fmkl = max(F1, F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmkl << endl;
    }

    while (2 * v < Fmkl) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(mork + 2 * v, mork) +
              static_cast<double100>(-(mork + 2 * v) *
                                     (mork + 2 * v + theta - 1) * sorts / 2.0));
      double100 newcoefficientL = exp(
          Getlogakm<double100>(mork + 2 * v + 1, mork) +
          static_cast<double100>(-(mork + 2 * v + 1) *
                                 (mork + 2 * v + 1 + theta - 1) * sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeUnconditional(eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeUnconditional(eCvec, 2 * v + 2, x, z, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmkj.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmkj = " << Fmkl
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFUnconditionalApprox(x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double100>(0.0)
                          : exp(log(eAL[n]) + log(qApprox) - log(eCU)));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    /// Add on l contributions
    int llower =
            (((thetaP.empty() || !(thetaP[0] > 0.0)) && !(z > 0.0)) ? 0 : 1),
        lupper =
            (((thetaP.empty() || !(thetaP[1] > 0.0)) && !(z < 1.0)) ? m
                                                                    : m - 1);
    for (int l = llower; l <= lupper && !mkl_found; ++l) {
      Amkl = computeAUnconditional(m, k, l, x, z);
      if (o.debug > 2)
        cerr << "Adding to currsums with A(n,m,k,l) = A(" << n << "," << m
             << "," << k << "," << l << ") = " << Amkl << endl;

      if (Amkl * dL[n] > 1.0 || Amkl * dU[n] < 0.0 || eAU[n] < 0.0 ||
          eAL[n] > 1.0 || eCU < 0.0) {
        cerr << "Numerical error detected: (m,k,l) = " << m << "," << k << ","
             << l << "), Amkl = " << Amkl << ", dL[" << n << "] = " << dL[n]
             << ", dU[" << n << "] = " << dU[n] << ", eAU[" << n
             << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
             << ", eCU = " << eCU << ", eCL = " << eCL
             << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
             << "," << z << "," << s << "," << t << ") ..." << endl;
        return DrawBridgePMFUnconditionalApprox(x, z, s, t, o, gen);
      }

      currsumL += Amkl * dL[n];
      currsumU += Amkl * dU[n];

      currsumStore.push_back(log(Amkl) + log(dL[n]));
      runningMax = max(runningMax, currsumStore.back());
      mkljStore.resize(mkljStore.size() + 4, -1);
      mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
      mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
      mkljStore[2 + (4 * (currsumStore.size() - 1))] = l;
      mkljStore[3 + (4 * (currsumStore.size() - 1))] = (!(z > 0.0) ? 0 : k);

      bool decision_on_mkl_made = (currsumL > u || currsumU < u);

      if (currsumL > currsumU) {
        cerr << "Error: currsumL = " << currsumL << " > " << currsumU
             << " = currsumU (n,m,k,l) = (" << n << "," << m << "," << k << ","
             << l << ")." << endl;
        exit(1);
      }

      while (!decision_on_mkl_made)  /// Refine upper and lower bounds
      {
        double100 currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double100> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
              morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
          ++v_used[i];
          double100 newcoefficientU =
              exp(Getlogakm<double100>(morki + 2 * v_used[i], morki) +
                  static_cast<double100>(-(morki + 2 * v_used[i]) *
                                         (morki + 2 * v_used[i] + theta - 1) *
                                         sorts / 2.0));
          double100 newcoefficientL =
              exp(Getlogakm<double100>(morki + 2 * v_used[i] + 1, morki) +
                  static_cast<double100>(
                      -(morki + 2 * v_used[i] + 1) *
                      (morki + 2 * v_used[i] + 1 + theta - 1) * sorts / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          qApprox = DiscretisedNormCDF(morkapproxi, sortsapprox);

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU -
                  GetdBridgeUnconditional(eCvec, 2 * v_used[i] + 1, x, z, t);
            eCU = eCL +
                  GetdBridgeUnconditional(eCvec, 2 * v_used[i] + 2, x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                             : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
          dL[i] = (eAL[i] < 0.0 ? static_cast<double100>(0.0)
                                : exp(log(eAL[i]) + log(qApprox) - log(eCU)));

          if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
              static_cast<double>(eAU[i]) < 0.0 ||
              static_cast<double>(eAL[i]) > 1.0 ||
              static_cast<double>(eCU) < 0.0) {
            cerr << "Numerical error detected: (m,k,l,j) = " << mi << "," << ki
                 << ", *, *), dL[" << i << "] = " << dL[i] << ", dU[" << i
                 << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
                 << ", eAL[" << i << "] = " << eAL[i] << ", eCU = " << eCU
                 << ", eCL = " << eCL
                 << ". Resorting to G1984-style approximation (x,z,s,t) = ("
                 << x << "," << z << "," << s << "," << t << ") ..." << endl;
            return DrawBridgePMFUnconditionalApprox(x, z, s, t, o, gen);
          }

          int lilower =
                  (((thetaP.empty() || !(thetaP[0] > 0.0)) && !(z > 0.0)) ? 0
                                                                          : 1),
              liupper = (((thetaP.empty() || !(thetaP[1] > 0.0)) && !(z < 1.0))
                             ? mi
                             : mi - 1);
          for (int l2 = lilower; l2 <= liupper; ++l2) {
            Amkl = computeAUnconditional(mi, ki, l2, x, z);
            currsumL += Amkl * dL[i];
            currsumU += Amkl * dU[i];

            currsumStore[i] = log(Amkl) + log(dL[i]);
            runningMax = max(runningMax, currsumStore[i]);

            if (o.debug > 3)
              cerr << "Recomputing currsums with A(n,m,k,l) = A(" << i << ","
                   << mi << "," << ki << "," << l2 << ") = " << Amkl << endl;
          }
        }

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k,l) = (" << n << "," << m << "," << k
               << "," << l << ")." << endl;
          exit(1);
        }

        if (o.debug > 2) {
          cerr << "\ndL ";
          printVec(dL, cerr);
          cerr << "\ndU ";
          printVec(dU, cerr);
        }

        if (currsumLold > currsumL) {
          // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,l) = (" << n << "," << m
                    << "," << k << "," << l << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
          std::cout << "Error: currsumUold = " << currsumUold << " < "
                    << currsumU << " = currsumU (n,m,k,l) = (" << n << "," << m
                    << "," << k << "," << l << ")." << endl;
          // exit(1);
        }

        decision_on_mkl_made = (currsumL > u || currsumU < u);
      }

      mkl_found = (currsumL > u);

      if (mkl_found) {
        if (!(z > 0.0)) {
          curr_mk[n][2] = l;  /// Sets j = 0
          curr_mk[n][3] = 0;
        } else {
          curr_mk[n][2] = l;  /// Sets j = k
          curr_mk[n][3] = k;
        }
      }
    }

    if (counter ==
        totalpts)  /// Gaussian approximation leads to currsum summing to < 1.0,
                   /// so we renormalise and sample from the computed quantities
    {
      LogSumExp(currsumStore, runningMax);
      double100 sum = 0.0;
      int ind = 0;

      bool found = false;

      while (!found) {
        sum += currsumStore[ind];
        if (sum > u) {
          found = true;
        }
        if (ind == static_cast<int>(currsumStore.size() - 1)) {
          found = true;
        }
        ind++;
      }

      vector<int> returnvec;

      returnvec.push_back(mkljStore[0 + (4 * ind)]);
      returnvec.push_back(mkljStore[1 + (4 * ind)]);
      returnvec.push_back(mkljStore[2 + (4 * ind)]);
      returnvec.push_back(mkljStore[3 + (4 * ind)]);

      return returnvec;
    }
  }

  int coeffcount = 0;
  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);
  if (o.debug > 2) {
    cerr << "p_m,k,j: Returned (m,k,j) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << "," << curr_mk[n][3] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFUnconditionalApprox(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting within (0,1), ending
               /// at some boundary point (absorption can happen at any time)
               /// with both time increments falling below the threshold
{
  assert((x > 0.0) && (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (t > s) &&
         (s > 0.0));
  vector<int> returnvec;
  vector<double100> currsumStore;
  vector<int> mkljStore;

  /// Compute denominator
  double100 eC = 0.0, eCInc = 1.0, eCOldInc = 1.0;
  int Dflip = 1, Djm = 0, Djp = 0, mkLower = (thetaP.empty() ? 1 : 0);
  int dmode = static_cast<int>(ceil(GriffithsParas(t).first)), d = dmode;

  while (max(eCOldInc, eCInc) > 0.0 || (!(eC > 0.0))) {
    eCOldInc = eCInc;
    double100 xcont = (!(z > 0.0) ? static_cast<double100>(d) * log(1.0 - x)
                                  : static_cast<double100>(d) * log(x));

    eCInc = exp(log(QmApprox(d, t, o)) + xcont);
    eC += eCInc;

    if (Dflip == -1 &&
        (dmode - Djm - 1 > 0))  /// Mechanism to explore either side around mode
    {
      Djm++;
      d = dmode - Djm;
    } else {
      Djp++;
      d = dmode + Djp;
    }
    Dflip *= -1;
  }

  vector<int> modeGuess = mklModeFinder(
      x, z, s, t, o);  /// Get a guess to location of mode over (m,k,l)
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2];

  boost::random::uniform_01<double100>
      U01;  /// Use these guesses & eC to set a suitable threshold for
  /// subsequent computations
  double100 currsum = 0.0, u = U01(gen),
            threshold = exp(mklModeFinder_Evaluator(mMode, kMode, lMode, x, z,
                                                    s, t, o) -
                            log(eC)) *
                        1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0, kFlip = 1, kD = 0, kU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double100 mContr_D = boost::math::lgamma(static_cast<double100>(theta + m));
  double100 mContr_U = mContr_D, mContr,
            runningMax = -1.0e100;  /// Computing m contributions

  while (!mSwitch) {
    double100 qm = QmApprox(m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to very
                      /// small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    if (m != mMode) {
      if (mU > mD) {
        mContr_U += log(static_cast<double100>(theta + m - 1));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D -= log(static_cast<double100>(theta + m + 1));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode, j = (!(z > 0.0) ? 0 : k);
    kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Computing k contributions
    double100 kContr_D =
                  -boost::math::lgamma(static_cast<double100>(theta + m + k)),
              kContr_U = kContr_D, kContr;

    while (!kSwitch) {
      double100 qk = QmApprox(k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      if (k != kMode) {
        if (kU > kD) {
          kContr_U -= log(static_cast<double100>(theta + m + k - 1));
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D += log(static_cast<double100>(theta + m + k + 1));
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      int lFlip = 1, newlMode = min(m, lMode), l = newlMode, lU = 0,
          lD = 0;  /// Need to redefine lMode as m might be too small!
      int lLower = ((thetaP.empty() && !(z < 1.0)) ? 1 : 0),
          lUpper = ((thetaP.empty() && !(z > 0.0)) ? m : m - 1);
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;

      boost::math::binomial_distribution<> BIN(m, x);
      double100 lContr_D =
          log(pdf(BIN, l)) +
          (!(z > 0.0)
               ? boost::math::lgamma(
                     static_cast<double100>(theta + m - l + k)) -
                     boost::math::lgamma(static_cast<double100>(theta + m - l))
               : boost::math::lgamma(static_cast<double100>(theta + l + k)) -
                     boost::math::lgamma(static_cast<double100>(theta + l)));
      double100 lContr_U = lContr_D, lContr;

      while (!lSwitch) {
        if (l != newlMode) {
          if (lU > lD) {
            lContr_U +=
                log(static_cast<double100>(m - (l - 1))) -
                log(static_cast<double100>((l - 1) + 1)) + log(x) -
                log(1.0 - x) +
                (!(z > 0.0)
                     ? log(static_cast<double100>(theta + m - (l - 1) - 1)) -
                           log(static_cast<double100>(theta + m - (l - 1) + k -
                                                      1))
                     : log(static_cast<double100>(theta + (l - 1) + k)) -
                           log(static_cast<double100>(theta + (l - 1))));
            lContr = lContr_U;
          } else {
            lContr_D +=
                log(static_cast<double100>(l + 1)) -
                log(static_cast<double100>(m - (l + 1) + 1)) + log(1.0 - x) -
                log(x) +
                (!(z > 0.0)
                     ? log(static_cast<double100>(theta + m - (l + 1) + k)) -
                           log(static_cast<double100>(theta + m - (l + 1)))
                     : log(static_cast<double100>(theta + (l + 1) - 1)) -
                           log(static_cast<double100>(theta + (l + 1) + k -
                                                      1)));
            lContr = lContr_D;
          }
        } else {
          lContr = lContr_U;
        }
        double100 currsum_inc = mContr + kContr + lContr - log(eC);
        runningMax =
            max(currsum_inc, runningMax);  /// Running max of log probabilities
        /// for use in log-sum-exp trick
        currsum += exp(currsum_inc);
        currsumStore.push_back(currsum_inc);

        mkljStore.resize(mkljStore.size() + 4, -1);
        mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
        mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
        mkljStore[2 + (4 * (currsumStore.size() - 1))] = l;
        mkljStore[3 + (4 * (currsumStore.size() - 1))] = j;

        if (!(lDownSwitch))  /// Switching mechanism for l
        {
          if (sgn(l - newlMode) <= 0) {
            lDownSwitch = ((exp(currsum_inc) < threshold) ||
                           (newlMode - lD - 1) < lLower);
          }
        }

        if (!(lUpSwitch)) {
          if (sgn(l - newlMode) >= 0) {
            lUpSwitch = ((exp(currsum_inc) < threshold) ||
                         (newlMode + lU + 1) > lUpper);
          }
        }

        lSwitch = (lDownSwitch && lUpSwitch);

        if (!lSwitch) {
          if ((lFlip == 1 && (newlMode + lU + 1 <= lUpper)) ||
              (lDownSwitch && !(lUpSwitch))) {
            lU++;
            l = newlMode + lU;
            lFlip *= (lDownSwitch ? 1 : -1);
          } else if ((lFlip == -1 && (newlMode - lD - 1 >= lLower)) ||
                     (lUpSwitch && !(lDownSwitch))) {
            lD++;
            l = newlMode - lD;
            lFlip *= (lUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch =
              (((lU == 0) && (lD == 0)) || (kMode - kD - 1 < mkLower));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = ((lU == 0) && (lD == 0));
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1) && (kMode - kD - 1 >= mkLower)) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Switching mechanism for m
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < mkLower));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= mkLower)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }

  LogSumExp(currsumStore, runningMax);
  double100 sum = 0.0;
  int index, ind = 0;

  vector<int> indexing(currsumStore.size(), 0);
  for (int i = 0; i != static_cast<int>(indexing.size()); i++) {
    indexing[i] = i;
  }
  sort(indexing.begin(), indexing.end(), [&](const int &a, const int &b) {
    return (currsumStore[a] > currsumStore[b]);
  });

  bool found = false;

  while (!found) {
    sum += currsumStore[indexing[ind]];
    if (sum > u) {
      index = indexing[ind];
      found = true;
    }
    if (ind == static_cast<int>(currsumStore.size() - 1)) {
      index = indexing[ind];
      found = true;
    }
    ind++;
  }

  returnvec.push_back(mkljStore[0 + (4 * index)]);
  returnvec.push_back(mkljStore[1 + (4 * index)]);
  returnvec.push_back(mkljStore[2 + (4 * index)]);
  returnvec.push_back(mkljStore[3 + (4 * index)]);

  return returnvec;
}

vector<int> WrightFisher::DrawBridgePMFSameMutation(
    double100 x, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending at the
               /// same boundary point and time increments are large enough
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double100> U01;
  vector<double100> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double100 u = U01(gen), currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(GetdBridgeSame(eCvec, 0, x, t));
  vector<int> v_used;
  map<vector<int>, double100> Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0;

  /// Compute F ensuring theoretical convergence of lower and upper sums - F1
  /// and F2 depend on m and k but F3-5 do not. None depend on j.
  pair<vector<int>, double100> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    eBL.push_back(0.0);
    eBU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(m, Cs), F2 = computeC(k, Cts);
    if (o.debug > 2)
      cerr << "(F1,F2,F3,F4,F5) = (" << F1 << "," << F2 << "," << F3 << ","
           << F4 << "," << F5 << ")" << endl;
    Fmkj = max(max(F1, F2), F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << s << "," << t << ") = " << Fmkj << endl;
    }

    while (2 * v < Fmkj || eAU[n] > 1.0 || eBU[n] > 1.0) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(m + 2 * v, m) +
              static_cast<double100>(-(m + 2 * v) * (m + 2 * v + theta - 1) *
                                     s / 2.0));
      double100 newcoefficientL =
          exp(Getlogakm<double100>(m + 2 * v + 1, m) +
              static_cast<double100>(-(m + 2 * v + 1) *
                                     (m + 2 * v + 1 + theta - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double100>(k + 2 * v, k) +
              static_cast<double100>(-(k + 2 * v) * (k + 2 * v + theta - 1) *
                                     (t - s) / 2.0));
      newcoefficientL = exp(Getlogakm<double100>(k + 2 * v + 1, k) +
                            static_cast<double100>(-(k + 2 * v + 1) *
                                                   (k + 2 * v + 1 + theta - 1) *
                                                   (t - s) / 2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeSame(eCvec, 2 * v + 1, x, t);
        eCU = eCL + GetdBridgeSame(eCvec, 2 * v + 2, x, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] && eBU[n] == eBL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmkj.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmkj = " << Fmkj
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eBU[n] - eBL[n] << ", " << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFSameMutationApprox(x, s, t, gen);
        }

        break;
      }
    }

    double100 addon;  /// Compute the appropriate additional terms

    if ((thetaP[0] > 0.0 && thetaP[1] > 0.0) && !(x > 0.0)) {
      addon = log(boost::math::beta<double100>(thetaP[0], thetaP[1] + m + k)) -
              log(boost::math::beta<double100>(thetaP[0], thetaP[1] + m)) -
              log(boost::math::beta<double100>(thetaP[0], thetaP[1] + k));
    } else if (thetaP[0] > 0.0 && thetaP[1] > 0.0) {
      addon = log(boost::math::beta<double100>(thetaP[0] + m + k, thetaP[1])) -
              log(boost::math::beta<double100>(thetaP[0] + m, thetaP[1])) -
              log(boost::math::beta<double100>(thetaP[0] + k, thetaP[1]));
    } else {
      addon = log(boost::math::beta<double100>(m + k, theta)) -
              log(boost::math::beta<double100>(m, theta)) -
              log(boost::math::beta<double100>(k, theta));
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) + addon - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[n]) + log(eBL[n]) + addon - log(eCU)));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    currsumL += dL[n];
    currsumU += dU[n];

    bool decision_on_mk_made = (currsumL > u || currsumU < u);

    if (currsumL > currsumU) {
      cerr << "Error: currsumL = " << currsumL << " > " << currsumU
           << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
           << endl;
      exit(1);
    }

    while (!decision_on_mk_made)  /// Need to refine upper and lower bounds
    {
      double100 currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double100> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        ++v_used[i];
        double100 newcoefficientU =
                      exp(Getlogakm<double100>(mi + 2 * v_used[i], mi) +
                          static_cast<double100>(
                              -(mi + 2 * v_used[i]) *
                              (mi + 2 * v_used[i] + theta - 1) * s / 2.0)),
                  newcoefficientL =
                      exp(Getlogakm<double100>(mi + 2 * v_used[i] + 1, mi) +
                          static_cast<double100>(
                              -(mi + 2 * v_used[i] + 1) *
                              (mi + 2 * v_used[i] + 1 + theta - 1) * s / 2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        newcoefficientU =
            exp(Getlogakm<double100>(ki + 2 * v_used[i], ki) +
                static_cast<double100>(-(ki + 2 * v_used[i]) *
                                       (ki + 2 * v_used[i] + theta - 1) *
                                       (t - s) / 2.0));
        newcoefficientL =
            exp(Getlogakm<double100>(ki + 2 * v_used[i] + 1, ki) +
                static_cast<double100>(-(ki + 2 * v_used[i] + 1) *
                                       (ki + 2 * v_used[i] + 1 + theta - 1) *
                                       (t - s) / 2.0));

        eBU[i] = eBL[i] + newcoefficientU;
        eBL[i] = eBU[i] - newcoefficientL;

        if (2 * (v_used[i] + 1) > eCindex_computed) {
          assert(2 * v_used[i] == eCindex_computed);
          eCL = eCU - GetdBridgeSame(eCvec, 2 * v_used[i] + 1, x, t);
          eCU = eCL + GetdBridgeSame(eCvec, 2 * v_used[i] + 2, x, t);
          eCindex_computed += 2;
        }

        if ((thetaP[0] > 0.0 && thetaP[1] > 0.0) && !(x > 0.0)) {
          addon = log(boost::math::beta<double100>(thetaP[0],
                                                   thetaP[1] + mi + ki)) -
                  log(boost::math::beta<double100>(thetaP[0], thetaP[1] + mi)) -
                  log(boost::math::beta<double100>(thetaP[0], thetaP[1] + ki));
        } else if (thetaP[0] > 0.0 && thetaP[1] > 0.0) {
          addon = log(boost::math::beta<double100>(thetaP[0] + mi + ki,
                                                   thetaP[1])) -
                  log(boost::math::beta<double100>(thetaP[0] + mi, thetaP[1])) -
                  log(boost::math::beta<double100>(thetaP[0] + ki, thetaP[1]));
        } else {
          addon = log(boost::math::beta<double100>(mi + ki, theta)) -
                  log(boost::math::beta<double100>(mi, theta)) -
                  log(boost::math::beta<double100>(ki, theta));
        }

        dU[i] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                           : exp(log(eAU[i]) + log(eBU[i]) + addon - log(eCL)));
        dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                     ? static_cast<double100>(0.0)
                     : exp(log(eAL[i]) + log(eBL[i]) + addon - log(eCU)));

        if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
            static_cast<double>(eAU[i]) < 0.0 ||
            static_cast<double>(eBU[i]) < 0.0 ||
            static_cast<double>(eAL[i]) > 1.0 ||
            static_cast<double>(eBL[i]) > 1.0 ||
            static_cast<double>(eCU) < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << mi << "," << ki
               << "), dL[" << i << "] = " << dL[i] << ", dU[" << i
               << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
               << ", eAL[" << i << "] = " << eAL[i] << ",  eBU[" << i
               << "] = " << eBU[i] << ", eBL[" << i << "] = " << eBL[i]
               << ", eCU = " << eCU << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFSameMutationApprox(x, s, t, gen);
        }

        currsumL += dL[i];
        currsumU += dU[i];

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
               << endl;
          exit(1);
        }
      }

      if (o.debug > 2) {
        cerr << "\ndL ";
        printVec(dL, cerr);
        cerr << "\ndU ";
        printVec(dU, cerr);
      }

      if (currsumLold > currsumL) {
        // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
        std::cout << "Error: currsumLold = " << currsumLold << " > " << currsumL
                  << " = currsumL (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }
      if (currsumUold < currsumU) {
        // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
        std::cout << "Error: currsumUold = " << currsumUold << " < " << currsumU
                  << " = currsumU (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }

      decision_on_mk_made = (currsumL > u || currsumU < u);
    }

    mk_found = (currsumL > u);
  }

  if (!(x > 0.0)) {
    curr_mk[n][2] = 0;  /// Setting l & j to be 0
    curr_mk[n][3] = 0;
  } else {
    curr_mk[n][2] = curr_mk[n][0];  /// Setting l & j to be m & k respectively
    curr_mk[n][3] = curr_mk[n][1];
  }

  int coeffcount = 0;

  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,l,j: Returned (m,k) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFSameMutationOneQApprox(
    double100 x, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending at the
               /// same boundary point, but one of the time increments falls
               /// below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double100 sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec, eAL, eAU, dL, dU, currsumStore;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(GetdBridgeSame(eCvec, 0, x, t)), qApprox = 0.0,
            runningMax = 1.0e-300;
  vector<int> v_used;
  map<vector<int>, double100> Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0,
      threshold = (thetaP.empty()
                       ? 2
                       : ((!(thetaP[0] > 0.0) || !(thetaP[1] > 0.0)) ? 1 : 0));

  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to be
  /// needed
  double mmode = GriffithsParas(s).first, vm = GriffithsParas(s).second;
  double kmode = GriffithsParas(t - s).first, vk = GriffithsParas(t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(threshold, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(threshold, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double100> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    int mork = (ind1 ? m : k), morkapprox = (mork == m ? k : m);
    if (m <= mlimU && m >= mlimL && k <= klimU && k >= klimL)
      counter++;  /// Keep track of the points we have visited inside the
    /// specified (m,k)-square
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(mork, Csorts);
    if (o.debug > 2)
      cerr << "(F1,F3,F4,F5) = (" << F1 << "," << F3 << "," << F4 << "," << F5
           << ")" << endl;
    Fmkj = max(F1, F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << s << "," << t << ") = " << Fmkj << endl;
    }

    while (2 * v < Fmkj) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(mork + 2 * v, mork) +
              static_cast<double100>(-(mork + 2 * v) *
                                     (mork + 2 * v + theta - 1) * sorts / 2.0));
      double100 newcoefficientL = exp(
          Getlogakm<double100>(mork + 2 * v + 1, mork) +
          static_cast<double100>(-(mork + 2 * v + 1) *
                                 (mork + 2 * v + 1 + theta - 1) * sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeSame(eCvec, 2 * v + 1, x, t);
        eCU = eCL + GetdBridgeSame(eCvec, 2 * v + 2, x, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmkj.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmkj = " << Fmkj
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFSameMutationApprox(x, s, t, gen);
        }

        break;
      }
    }

    double100 addon;  /// Computing the appropriate additional contributions

    if ((thetaP[0] > 0.0 && thetaP[1] > 0.0) && !(x > 0.0)) {
      addon = log(boost::math::beta<double100>(thetaP[0], thetaP[1] + m + k)) -
              log(boost::math::beta<double100>(thetaP[0], thetaP[1] + m)) -
              log(boost::math::beta<double100>(thetaP[0], thetaP[1] + k));
    } else if (thetaP[0] > 0.0 && thetaP[1] > 0.0) {
      addon = log(boost::math::beta<double100>(thetaP[0] + m + k, thetaP[1])) -
              log(boost::math::beta<double100>(thetaP[0] + m, thetaP[1])) -
              log(boost::math::beta<double100>(thetaP[0] + k, thetaP[1]));
    } else {
      addon = log(boost::math::beta<double100>(m + k, theta)) -
              log(boost::math::beta<double100>(m, theta)) -
              log(boost::math::beta<double100>(k, theta));
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) + addon - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double100>(0.0)
                          : exp(log(eAL[n]) + log(qApprox) + addon - log(eCU)));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    currsumL += dL[n];
    currsumU += dU[n];

    currsumStore.push_back(log(dL[n]));
    runningMax = max(runningMax, currsumStore.back());
    mkljStore.resize(mkljStore.size() + 4, -1);
    mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
    mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
    mkljStore[2 + (4 * (currsumStore.size() - 1))] = (!(x > 0.0) ? 0 : m);
    mkljStore[3 + (4 * (currsumStore.size() - 1))] = (!(x > 0.0) ? 0 : k);

    if (counter == totalpts)  /// Gaussian approximation leads to currsum
                              /// summing to < 1.0, so we renormalise and sample
    {
      LogSumExp(currsumStore, runningMax);
      double100 sum = 0.0;
      int ind = 0;

      bool found = false;

      while (!found) {
        sum += currsumStore[ind];
        if (sum > u) {
          found = true;
        }
        if (ind == static_cast<int>(currsumStore.size() - 1)) {
          found = true;
        }
        ind++;
      }

      vector<int> returnvec;

      returnvec.push_back(mkljStore[0 + (4 * ind)]);
      returnvec.push_back(mkljStore[1 + (4 * ind)]);
      returnvec.push_back(mkljStore[2 + (4 * ind)]);
      returnvec.push_back(mkljStore[3 + (4 * ind)]);

      return returnvec;
    }

    bool decision_on_mk_made = (currsumL > u || currsumU < u);

    if (currsumL > currsumU) {
      cerr << "Error: currsumL = " << currsumL << " > " << currsumU
           << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
           << endl;
      exit(1);
    }

    while (!decision_on_mk_made)  /// Need to refine the upper and lower bounds
    {
      double100 currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double100> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        int morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
        ++v_used[i];
        double100 newcoefficientU =
            exp(Getlogakm<double100>(morki + 2 * v_used[i], morki) +
                static_cast<double100>(-(morki + 2 * v_used[i]) *
                                       (morki + 2 * v_used[i] + theta - 1) *
                                       sorts / 2.0));
        double100 newcoefficientL =
            exp(Getlogakm<double100>(morki + 2 * v_used[i] + 1, morki) +
                static_cast<double100>(-(morki + 2 * v_used[i] + 1) *
                                       (morki + 2 * v_used[i] + 1 + theta - 1) *
                                       sorts / 2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        qApprox = DiscretisedNormCDF(morkapproxi, sortsapprox);

        if (2 * (v_used[i] + 1) > eCindex_computed) {
          assert(2 * v_used[i] == eCindex_computed);
          eCL = eCU - GetdBridgeSame(eCvec, 2 * v_used[i] + 1, x, t);
          eCU = eCL + GetdBridgeSame(eCvec, 2 * v_used[i] + 2, x, t);
          eCindex_computed += 2;
        }

        if ((thetaP[0] > 0.0 && thetaP[1] > 0.0) && !(x > 0.0)) {
          addon = log(boost::math::beta<double100>(thetaP[0],
                                                   thetaP[1] + mi + ki)) -
                  log(boost::math::beta<double100>(thetaP[0], thetaP[1] + mi)) -
                  log(boost::math::beta<double100>(thetaP[0], thetaP[1] + ki));
        } else if (thetaP[0] > 0.0 && thetaP[1] > 0.0) {
          addon = log(boost::math::beta<double100>(thetaP[0] + mi + ki,
                                                   thetaP[1])) -
                  log(boost::math::beta<double100>(thetaP[0] + mi, thetaP[1])) -
                  log(boost::math::beta<double100>(thetaP[0] + ki, thetaP[1]));
        } else {
          addon = log(boost::math::beta<double100>(mi + ki, theta)) -
                  log(boost::math::beta<double100>(mi, theta)) -
                  log(boost::math::beta<double100>(ki, theta));
        }

        dU[i] =
            (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[i]) + log(qApprox) + addon - log(eCL)));
        dL[i] =
            (eAL[i] < 0.0 ? static_cast<double100>(0.0)
                          : exp(log(eAL[i]) + log(qApprox) + addon - log(eCU)));

        if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
            static_cast<double>(eAU[i]) < 0.0 ||
            static_cast<double>(eAL[i]) > 1.0 ||
            static_cast<double>(eCU) < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << mi << "," << ki
               << "), dL[" << i << "] = " << dL[i] << ", dU[" << i
               << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
               << ", eAL[" << i << "] = " << eAL[i] << ", eCU = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFSameMutationApprox(x, s, t, gen);
        }

        currsumL += dL[i];
        currsumU += dU[i];

        currsumStore[i] = log(dL[i]);
        runningMax = max(runningMax, currsumStore[i]);

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
               << endl;
          exit(1);
        }
      }

      if (o.debug > 2) {
        cerr << "\ndL ";
        printVec(dL, cerr);
        cerr << "\ndU ";
        printVec(dU, cerr);
      }

      if (currsumLold > currsumL) {
        // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
        std::cout << "Error: currsumLold = " << currsumLold << " > " << currsumL
                  << " = currsumL (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }
      if (currsumUold < currsumU) {
        // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
        std::cout << "Error: currsumUold = " << currsumUold << " < " << currsumU
                  << " = currsumU (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }

      decision_on_mk_made = (currsumL > u || currsumU < u);
    }

    mk_found = (currsumL > u);
  }

  if (!(x > 0.0)) {
    curr_mk[n][2] = 0;  /// Setting l & j to be 0
    curr_mk[n][3] = 0;
  } else {
    curr_mk[n][2] = curr_mk[n][0];  /// Setting l & j to be m & k respectively
    curr_mk[n][3] = curr_mk[n][1];
  }

  int coeffcount = 0;

  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,l,j: Returned (m,k) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFSameMutationApprox(
    double100 x, double100 s, double100 t,
    boost::random::mt19937 &
        gen)  /// Draws from the law of a bridge starting and ending at the same
              /// boundary point, but both time increments fall below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  vector<int> returnvec;
  bool accept = false;
  int m, k;

  /// Return (m,k) by running a rejection sampler
  double mean = GriffithsParas(t).first, v = GriffithsParas(t).second;
  int l = max(static_cast<int>(round(mean - 4.0 * sqrt(v))), 1);
  /// Precompute an upper bound for the rejection sampler
  double norm =
      boost::math::lgamma(static_cast<double100>(thetaP[1] + (2 * l))) +
      boost::math::lgamma(static_cast<double100>(thetaP[0] + thetaP[1] + l)) +
      boost::math::lgamma(static_cast<double100>(thetaP[0] + thetaP[1] + l)) -
      boost::math::lgamma(
          static_cast<double100>(thetaP[0] + thetaP[1] + (2 * l))) -
      boost::math::lgamma(static_cast<double100>(thetaP[1] + l)) -
      boost::math::lgamma(static_cast<double100>(thetaP[1] + l));

  while (!accept) {
    m = DrawSizebiasedAncestralProcess(2, s, gen),
    k = DrawSizebiasedAncestralProcess(
        2, t - s, gen);  /// Draw (m,k) from a size-biased Ancestral Process
    boost::random::uniform_01<>
        U01;  /// Compute alpha and run an accept/reject step
    double u = U01(gen),
           alpha =
               boost::math::lgamma(static_cast<double100>(thetaP[1] + m + k)) +
               boost::math::lgamma(
                   static_cast<double100>(thetaP[0] + thetaP[1] + m)) +
               boost::math::lgamma(
                   static_cast<double100>(thetaP[0] + thetaP[1] + k)) -
               boost::math::lgamma(
                   static_cast<double100>(thetaP[0] + thetaP[1] + m + k)) -
               boost::math::lgamma(static_cast<double100>(thetaP[1] + m)) -
               boost::math::lgamma(static_cast<double100>(thetaP[1] + k));

    if (exp(alpha - norm) > u) {
      accept = true;
    }
  }

  returnvec.push_back(m);
  returnvec.push_back(k);

  if (!(x > 0.0)) {
    returnvec.push_back(0);  /// Setting l & j to be m & k respectively
    returnvec.push_back(0);
  } else {
    returnvec.push_back(m);  /// Setting l & j to be m & k respectively
    returnvec.push_back(k);
  }

  return returnvec;
}

vector<int> WrightFisher::DrawBridgePMFDifferentMutation(
    double100 s, double100 t, double100 x, const Options &o,
    boost::random::mt19937 &
        gen)  /// Draws from the law of a bridge starting and ending at
              /// different boundary points and time increments are large enough
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eAL, eAU, eBL, eBU, eCU, eCL, dL, dU;
  double100 currsumU = 0.0, currsumL = 0.0;
  vector<int> v_used;
  int n = -1, Fmk = 0,
      denomqindex = ((thetaP[0] > 0.0 && thetaP[1] > 0.0) ? 0 : 1);

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double100> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F4 = static_cast<int>(
      ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0)));
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    eBL.push_back(0.0);
    eBU.push_back(0.0);
    eCL.push_back(0.0);
    eCU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(m, Cs), F2 = computeC(k, Cts),
        F3 = computeC(denomqindex, Ct);
    if (o.debug > 2)
      cerr << "(F1,F2,F3) = (" << F1 << "," << F2 << "," << F3 << ")" << endl;
    Fmk = max(max(max(F1, F2), F3), F4);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << s << "," << t << ") = " << Fmk << endl;
    }

    while (2 * v < Fmk || eAU[n] > 1.0 || eBU[n] > 1.0 || eCU[n] > 1.0 ||
           eAL[n] < 0.0 || eBL[n] < 0.0 || eCL[n] < 0.0) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(m + 2 * v, m) +
              static_cast<double100>(-(m + 2 * v) * (m + 2 * v + theta - 1) *
                                     s / 2.0));
      double100 newcoefficientL =
          exp(Getlogakm<double100>(m + 2 * v + 1, m) +
              static_cast<double100>(-(m + 2 * v + 1) *
                                     (m + 2 * v + 1 + theta - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double100>(k + 2 * v, k) +
              static_cast<double100>(-(k + 2 * v) * (k + 2 * v + theta - 1) *
                                     (t - s) / 2.0));
      newcoefficientL = exp(Getlogakm<double100>(k + 2 * v + 1, k) +
                            static_cast<double100>(-(k + 2 * v + 1) *
                                                   (k + 2 * v + 1 + theta - 1) *
                                                   (t - s) / 2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double100>(denomqindex + 2 * v, denomqindex) +
              static_cast<double100>(-(denomqindex + 2 * v) *
                                     (denomqindex + 2 * v + theta - 1) * (t) /
                                     2.0));  // Computing q_2
      newcoefficientL =
          exp(Getlogakm<double100>(denomqindex + 2 * v + 1, denomqindex) +
              static_cast<double100>(-(denomqindex + 2 * v + 1) *
                                     (denomqindex + 2 * v + 1 + theta - 1) *
                                     (t) / 2.0));

      eCU[n] = eCL[n] + newcoefficientU;
      eCL[n] = eCU[n] - newcoefficientL;

      if (eAU[n] == eAL[n] && eBU[n] == eBL[n] &&
          eCL[n] ==
              eCU[n])  /// ...then we have lost precision before reaching Fmk.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmk = " << Fmk
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eBU[n] - eBL[n] << ", " << eCU[n] - eCL[n] << endl;
        }

        if (eAL[n] < 0.0 || eCL[n] < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU[n]
               << ", eCL = " << eCL[n]
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFDifferentMutationApprox(s, t, x, o, gen);
        }

        break;
      }
    }

    /// Compute the corresponding additional contributions
    double100 addon =
        boost::math::lgamma(static_cast<double100>(theta + k)) +
        boost::math::lgamma(static_cast<double100>(theta + m)) -
        boost::math::lgamma(static_cast<double100>(theta + m + k)) -
        boost::math::lgamma(static_cast<double100>(thetaP[0])) -
        boost::math::lgamma(static_cast<double100>(thetaP[1]));

    dU[n] = (eCL[n] < 0.0
                 ? static_cast<double100>(nan(""))
                 : exp(log(eAU[n]) + log(eBU[n]) + addon -
                       log(eCL[n] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[n]) + log(eBL[n]) + addon -
                       log(eCU[n] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    currsumL += dL[n];
    currsumU += dU[n];

    bool decision_on_mk_made = (currsumL > u || currsumU < u);

    if (currsumL > currsumU) {
      cerr << "Error: currsumL = " << currsumL << " > " << currsumU
           << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
           << endl;
      exit(1);
    }

    while (!decision_on_mk_made)  /// Refine upper and lower bounds
    {
      double100 currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double100> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        ++v_used[i];
        double100 newcoefficientU =
                      exp(Getlogakm<double100>(mi + 2 * v_used[i], mi) +
                          static_cast<double100>(
                              -(mi + 2 * v_used[i]) *
                              (mi + 2 * v_used[i] + theta - 1) * s / 2.0)),
                  newcoefficientL =
                      exp(Getlogakm<double100>(mi + 2 * v_used[i] + 1, mi) +
                          static_cast<double100>(
                              -(mi + 2 * v_used[i] + 1) *
                              (mi + 2 * v_used[i] + 1 + theta - 1) * s / 2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        newcoefficientU =
            exp(Getlogakm<double100>(ki + 2 * v_used[i], ki) +
                static_cast<double100>(-(ki + 2 * v_used[i]) *
                                       (ki + 2 * v_used[i] + theta - 1) *
                                       (t - s) / 2.0));
        newcoefficientL =
            exp(Getlogakm<double100>(ki + 2 * v_used[i] + 1, ki) +
                static_cast<double100>(-(ki + 2 * v_used[i] + 1) *
                                       (ki + 2 * v_used[i] + 1 + theta - 1) *
                                       (t - s) / 2.0));

        eBU[i] = eBL[i] + newcoefficientU;
        eBL[i] = eBU[i] - newcoefficientL;

        newcoefficientU =
            exp(Getlogakm<double100>(denomqindex + 2 * v_used[i], denomqindex) +
                static_cast<double100>(
                    -(denomqindex + 2 * v_used[i]) *
                    (denomqindex + 2 * v_used[i] + theta - 1) * (t) / 2.0));
        newcoefficientL = exp(
            Getlogakm<double100>(denomqindex + 2 * v_used[i] + 1, denomqindex) +
            static_cast<double100>(
                -(denomqindex + 2 * v_used[i] + 1) *
                (denomqindex + 2 * v_used[i] + 1 + theta - 1) * (t) / 2.0));

        eCU[i] = eCL[i] + newcoefficientU;
        eCL[i] = eCU[i] - newcoefficientL;

        addon = boost::math::lgamma(static_cast<double100>(theta + ki)) +
                boost::math::lgamma(static_cast<double100>(theta + mi)) -
                boost::math::lgamma(static_cast<double100>(theta + mi + ki)) -
                boost::math::lgamma(static_cast<double100>(thetaP[0])) -
                boost::math::lgamma(static_cast<double100>(thetaP[1]));

        dU[i] =
            (eCL[i] < 0.0
                 ? static_cast<double100>(nan(""))
                 : exp(log(eAU[i]) + log(eBU[i]) + addon -
                       log(eCL[i] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));
        dL[i] =
            ((eAL[i] < 0.0 || eBL[i] < 0.0)
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[i]) + log(eBL[i]) + addon -
                       log(eCU[i] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));

        if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
            static_cast<double>(eAU[i]) < 0.0 ||
            static_cast<double>(eBU[i]) < 0.0 ||
            static_cast<double>(eAL[i]) > 1.0 ||
            static_cast<double>(eBL[i]) > 1.0 ||
            static_cast<double>(eCU[i]) < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << mi << "," << ki
               << "), dL[" << i << "] = " << dL[i] << ", dU[" << i
               << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
               << ", eAL[" << i << "] = " << eAL[i] << ",  eBU[" << i
               << "] = " << eBU[i] << ", eBL[" << i << "] = " << eBL[i]
               << ", eCU = " << eCU[i] << ", eCL = " << eCL[i]
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFDifferentMutationApprox(s, t, x, o, gen);
        }

        currsumL += dL[i];
        currsumU += dU[i];

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
               << endl;
          exit(1);
        }
      }

      if (o.debug > 2) {
        cerr << "\ndL ";
        printVec(dL, cerr);
        cerr << "\ndU ";
        printVec(dU, cerr);
      }

      if (currsumLold > currsumL) {
        // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
        std::cout << "Error: currsumLold = " << currsumLold << " > " << currsumL
                  << " = currsumL (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }
      if (currsumUold < currsumU) {
        // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
        std::cout << "Error: currsumUold = " << currsumUold << " < " << currsumU
                  << " = currsumU (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }

      decision_on_mk_made = (currsumL > u || currsumU < u);
    }

    mk_found = (currsumL > u);
  }

  if (!(x > 0.0)) {
    curr_mk[n][2] = 0;  /// Setting l & j to be 0 & k respectively
    curr_mk[n][3] = curr_mk[n][1];
  } else {
    curr_mk[n][2] = curr_mk[n][0];  /// Setting l & j to be m & 0 respectively
    curr_mk[n][3] = 0;
  }

  int coeffcount = 0;

  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,l,j: Returned (m,k) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFDifferentMutationOneQApprox(
    double100 s, double100 t, double100 x, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending at
               /// different boundary points, but one of the time increments
               /// falls below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double100 sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eAL, eAU, eCU, eCL, dL, dU, currsumStore;
  double100 currsumU = 0.0, currsumL = 0.0, qApprox = 0.0,
            runningMax = 1.0e-300;
  vector<int> v_used;
  int n = -1, Fmk = 0,
      denomqindex = ((thetaP[0] > 0.0 && thetaP[1] > 0.0) ? 0 : 1);
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to be
  /// needed
  double mmode = GriffithsParas(s).first, vm = GriffithsParas(s).second,
         kmode = GriffithsParas(t - s).first, vk = GriffithsParas(t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double100> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F4 = static_cast<int>(
      ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0)));
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    int mork = (ind1 ? m : k), morkapprox = (mork == m ? k : m);
    if (m <= mlimU && m >= mlimL && k <= klimU && k >= klimL) counter++;
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    eCL.push_back(0.0);
    eCU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(mork, Csorts), F3 = computeC(denomqindex, Ct);
    if (o.debug > 2) cerr << "(F1,F3) = (" << F1 << "," << F3 << ")" << endl;
    Fmk = max(max(F1, F3), F4);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << s << "," << t << ") = " << Fmk << endl;
    }

    while (2 * v < Fmk) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(mork + 2 * v, mork) +
              static_cast<double100>(-(mork + 2 * v) *
                                     (mork + 2 * v + theta - 1) * sorts / 2.0));
      double100 newcoefficientL = exp(
          Getlogakm<double100>(mork + 2 * v + 1, mork) +
          static_cast<double100>(-(mork + 2 * v + 1) *
                                 (mork + 2 * v + 1 + theta - 1) * sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(morkapprox, sortsapprox);

      newcoefficientU =
          exp(Getlogakm<double100>(denomqindex + 2 * v, denomqindex) +
              static_cast<double100>(-(denomqindex + 2 * v) *
                                     (denomqindex + 2 * v + theta - 1) * (t) /
                                     2.0));  // Computing q_2
      newcoefficientL =
          exp(Getlogakm<double100>(denomqindex + 2 * v + 1, denomqindex) +
              static_cast<double100>(-(denomqindex + 2 * v + 1) *
                                     (denomqindex + 2 * v + 1 + theta - 1) *
                                     (t) / 2.0));

      eCU[n] = eCL[n] + newcoefficientU;
      eCL[n] = eCU[n] - newcoefficientL;

      if (eAU[n] == eAL[n] &&
          eCL[n] ==
              eCU[n])  /// ...then we have lost precision before reaching Fmk.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmk = " << Fmk
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eCU[n] - eCL[n] << endl;
        }

        if (eAL[n] < 0.0 || eCL[n] < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU[n]
               << ", eCL = " << eCL[n]
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFDifferentMutationApprox(s, t, x, o, gen);
        }

        break;
      }
    }

    /// Compute the appropriate additional contributions
    double100 addon =
        boost::math::lgamma(static_cast<double100>(theta + k)) +
        boost::math::lgamma(static_cast<double100>(theta + m)) -
        boost::math::lgamma(static_cast<double100>(theta + m + k)) -
        boost::math::lgamma(static_cast<double100>(thetaP[0])) -
        boost::math::lgamma(static_cast<double100>(thetaP[1]));

    dU[n] = (eCL[n] < 0.0
                 ? static_cast<double100>(nan(""))
                 : exp(log(eAU[n]) + log(qApprox) + addon -
                       log(eCL[n] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));
    dL[n] = (eAL[n] < 0.0
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[n]) + log(qApprox) + addon -
                       log(eCU[n] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    currsumL += dL[n];
    currsumU += dU[n];

    currsumStore.push_back(log(dL[n]));
    runningMax = max(runningMax, currsumStore.back());
    mkljStore.resize(mkljStore.size() + 4, -1);
    mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
    mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
    mkljStore[2 + (4 * (currsumStore.size() - 1))] = (!(x > 0.0) ? 0 : m);
    mkljStore[3 + (4 * (currsumStore.size() - 1))] = (!(x > 0.0) ? k : 0);

    if (counter == totalpts)  /// Gaussian approximation leads to currsum
                              /// summing to < 1.0, so we renormalise and sample
    {
      LogSumExp(currsumStore, runningMax);
      double100 sum = 0.0;
      int ind = 0;

      bool found = false;

      while (!found) {
        sum += currsumStore[ind];
        if (sum > u) {
          found = true;
        }
        if (ind == static_cast<int>(currsumStore.size() - 1)) {
          found = true;
        }
        ind++;
      }

      vector<int> returnvec;

      returnvec.push_back(mkljStore[0 + (4 * ind)]);
      returnvec.push_back(mkljStore[1 + (4 * ind)]);
      returnvec.push_back(mkljStore[2 + (4 * ind)]);
      returnvec.push_back(mkljStore[3 + (4 * ind)]);

      return returnvec;
    }

    bool decision_on_mk_made = (currsumL > u || currsumU < u);

    if (currsumL > currsumU) {
      cerr << "Error: currsumL = " << currsumL << " > " << currsumU
           << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
           << endl;
      exit(1);
    }

    while (!decision_on_mk_made)  /// Refine upper and lower bounds
    {
      double100 currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double100> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1], morki = (ind1 ? mi : ki),
            morkapproxi = (morki == mi ? ki : mi);
        ++v_used[i];
        double100 newcoefficientU = exp(
                      Getlogakm<double100>(morki + 2 * v_used[i], morki) +
                      static_cast<double100>(
                          -(morki + 2 * v_used[i]) *
                          (morki + 2 * v_used[i] + theta - 1) * sorts / 2.0)),
                  newcoefficientL = exp(
                      Getlogakm<double100>(morki + 2 * v_used[i] + 1, morki) +
                      static_cast<double100>(
                          -(morki + 2 * v_used[i] + 1) *
                          (morki + 2 * v_used[i] + 1 + theta - 1) * sorts /
                          2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        qApprox = DiscretisedNormCDF(morkapproxi, sortsapprox);

        newcoefficientU =
            exp(Getlogakm<double100>(denomqindex + 2 * v_used[i], denomqindex) +
                static_cast<double100>(
                    -(denomqindex + 2 * v_used[i]) *
                    (denomqindex + 2 * v_used[i] + theta - 1) * (t) / 2.0));
        newcoefficientL = exp(
            Getlogakm<double100>(denomqindex + 2 * v_used[i] + 1, denomqindex) +
            static_cast<double100>(
                -(denomqindex + 2 * v_used[i] + 1) *
                (denomqindex + 2 * v_used[i] + 1 + theta - 1) * (t) / 2.0));

        eCU[i] = eCL[i] + newcoefficientU;
        eCL[i] = eCU[i] - newcoefficientL;

        addon = boost::math::lgamma(static_cast<double100>(theta + ki)) +
                boost::math::lgamma(static_cast<double100>(theta + mi)) -
                boost::math::lgamma(static_cast<double100>(theta + mi + ki)) -
                boost::math::lgamma(static_cast<double100>(thetaP[0])) -
                boost::math::lgamma(static_cast<double100>(thetaP[1]));

        dU[i] =
            (eCL[i] < 0.0
                 ? static_cast<double100>(nan(""))
                 : exp(log(eAU[i]) + log(qApprox) + addon -
                       log(eCL[i] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));
        dL[i] =
            (eAL[i] < 0.0
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[i]) + log(qApprox) + addon -
                       log(eCU[i] / boost::math::beta<double100>(
                                        static_cast<double100>(thetaP[0]),
                                        static_cast<double100>(thetaP[1])))));

        if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
            static_cast<double>(eAU[i]) < 0.0 ||
            static_cast<double>(eAL[i]) > 1.0 ||
            static_cast<double>(eCU[i]) < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << mi << "," << ki
               << "), dL[" << i << "] = " << dL[i] << ", dU[" << i
               << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
               << ", eAL[" << i << "] = " << eAL[i] << ", eCU = " << eCU[i]
               << ", eCL = " << eCL[i]
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFDifferentMutationApprox(s, t, x, o, gen);
        }

        currsumL += dL[i];
        currsumU += dU[i];

        currsumStore[i] = log(dL[i]);
        runningMax = max(runningMax, currsumStore[i]);

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k) = (" << n << "," << m << "," << k << ")."
               << endl;
          exit(1);
        }
      }

      if (o.debug > 2) {
        cerr << "\ndL ";
        printVec(dL, cerr);
        cerr << "\ndU ";
        printVec(dU, cerr);
      }

      if (currsumLold > currsumL) {
        // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
        std::cout << "Error: currsumLold = " << currsumLold << " > " << currsumL
                  << " = currsumL (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }
      if (currsumUold < currsumU) {
        // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
        std::cout << "Error: currsumUold = " << currsumUold << " < " << currsumU
                  << " = currsumU (n,m,k) = (" << n << "," << m << "," << k
                  << ")." << endl;
        // exit(1);
      }

      decision_on_mk_made = (currsumL > u || currsumU < u);
    }

    mk_found = (currsumL > u);
  }

  if (!(x > 0.0)) {
    curr_mk[n][2] = 0;  /// Setting l & j to be 0 & k respectively
    curr_mk[n][3] = curr_mk[n][1];
  } else {
    curr_mk[n][2] = curr_mk[n][0];  /// Setting l & j to be m & 0 respectively
    curr_mk[n][3] = 0;
  }

  int coeffcount = 0;

  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,l,j: Returned (m,k) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFDifferentMutationApprox(
    double100 s, double100 t, double100 x, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending at
               /// different boundary points, but both time increments fall
               /// below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  vector<int> returnvec;
  vector<double100> currsumStore;
  vector<int> mkStore;

  /// Compute denominator depending on value of time increment t
  double100 eC =
      (t <= o.g1984threshold ? DiscretisedNormCDF(0, t) : QmApprox(0, t, o)) /
      boost::math::beta<double100>(static_cast<double100>(thetaP[0]),
                                   static_cast<double100>(thetaP[1]));

  /// Compute a guess to the mode over (m,k)
  vector<int> modeGuess = mkModeFinder(x, 1.0 - x, s, t, o);
  int mMode = modeGuess[0],
      kMode = modeGuess[1];  /// Use this guess & eC to obtain a suitable
  /// threshold for subsequent calculations
  double100 constContr = -log(eC) - boost::math::lgamma(thetaP[0]) -
                         boost::math::lgamma(thetaP[1]);
  double100 currsum = 0.0,
            threshold =
                exp(mkModeFinder_Evaluator(mMode, kMode, x, 1.0 - x, s, t, o) -
                    constContr) *
                1.0e-4,
            runningMax = -1.0e100;

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;

  while (!mSwitch) {
    int kFlip = 1, k = kMode, kU = 0, kD = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;

    while (!kSwitch) {
      double100 currsum_inc =
          log(QmApprox(m, s, o)) + log(QmApprox(k, t - s, o)) +
          boost::math::lgamma(static_cast<double100>(theta + k)) +
          boost::math::lgamma(static_cast<double100>(theta + m)) -
          boost::math::lgamma(static_cast<double100>(theta + m + k)) -
          constContr;

      currsum += exp(currsum_inc);
      runningMax = max(currsum_inc, runningMax);  /// Storing maximum of log
      /// probabilities for log-sum-exp trick
      currsumStore.push_back(currsum_inc);

      mkStore.resize(mkStore.size() + 2, -1);
      mkStore[0 + (2 * (currsumStore.size() - 1))] = m;
      mkStore[1 + (2 * (currsumStore.size() - 1))] = k;

      if (!(kDownSwitch))  /// Mechanism to explore k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch =
              ((exp(currsum_inc) < threshold) || (kMode - kD - 1) < 0);
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = (exp(currsum_inc) < threshold);
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1 || (kDownSwitch && !(kUpSwitch))) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1 && (kMode - kD - 1 >= 0)) ||
                   (kUpSwitch && !(kDownSwitch))) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Mechanism to explore m
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < 0));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= 0)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }

  LogSumExp(currsumStore, runningMax);  /// Log-sum-exp trick to normalise
  /// vector of log probabilities
  double100 sum = 0.0;
  int index, ind = 0;
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);

  vector<int> indexing(currsumStore.size(), 0);
  for (int i = 0; i != static_cast<int>(indexing.size()); i++) {
    indexing[i] = i;
  }
  sort(indexing.begin(), indexing.end(),  /// Sort vector
       [&](const int &a, const int &b) {
         return (currsumStore[a] > currsumStore[b]);
       });

  bool found = false;

  while (!found) {
    sum += currsumStore[indexing[ind]];
    if (sum > u) {
      index =
          indexing[ind];  /// Figure out what the correct index for mkljStore is
      found = true;
    }
    if (ind == static_cast<int>(currsumStore.size()))  /// Ending condition
    {
      index = indexing[ind];
      found = true;
    }
    ind++;
  }

  if (!(x > 0.0))  /// Setting l = 0, j = k
  {
    returnvec.push_back(mkStore[0 + (2 * index)]);
    returnvec.push_back(mkStore[1 + (2 * index)]);
    returnvec.push_back(0);
    returnvec.push_back(mkStore[1 + (2 * index)]);
  } else  /// Setting l = m, j = 0
  {
    returnvec.push_back(mkStore[0 + (2 * index)]);
    returnvec.push_back(mkStore[1 + (2 * index)]);
    returnvec.push_back(mkStore[0 + (2 * index)]);
    returnvec.push_back(0);
  }

  return returnvec;
}

vector<int> WrightFisher::DrawBridgePMFInteriorMutation(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting at a boundary point
               /// but ending in the interior of (0,1) with both time increments
               /// large enough
{
  assert((!(x > 0.0) || !(x < 1.0)) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkj_found = false;

  /// Setting up necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(GetdBridgeInterior(eCvec, 0, x, z, t));
  vector<int> v_used;
  double100 Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double100> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
      F5 = (!(x > 0.0) ? static_cast<int>(floor(
                             (theta + 2.0) / (o.eps * (thetaP[1] + 1.0)) - 1)) +
                             1
                       : static_cast<int>(floor(
                             (theta + 2.0) / (o.eps * (thetaP[0] + 1.0)) - 1)) +
                             1);
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    eBL.push_back(0.0);
    eBU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(m, Cs), F2 = computeC(k, Cts);
    if (o.debug > 2)
      cerr << "(F1,F2,F3,F4,F5) = (" << F1 << "," << F2 << "," << F3 << ","
           << F4 << "," << F5 << ")" << endl;
    Fmkj = max(max(F1, F2), F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmkj << endl;
    }

    while (2 * v < Fmkj || eAU[n] > 1.0 || eBU[n] > 1.0) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(m + 2 * v, m) +
              static_cast<double100>(-(m + 2 * v) * (m + 2 * v + theta - 1) *
                                     s / 2.0));
      double100 newcoefficientL =
          exp(Getlogakm<double100>(m + 2 * v + 1, m) +
              static_cast<double100>(-(m + 2 * v + 1) *
                                     (m + 2 * v + 1 + theta - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double100>(k + 2 * v, k) +
              static_cast<double100>(-(k + 2 * v) * (k + 2 * v + theta - 1) *
                                     (t - s) / 2.0));
      newcoefficientL = exp(Getlogakm<double100>(k + 2 * v + 1, k) +
                            static_cast<double100>(-(k + 2 * v + 1) *
                                                   (k + 2 * v + 1 + theta - 1) *
                                                   (t - s) / 2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeInterior(eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeInterior(eCvec, 2 * v + 2, x, z, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] && eBU[n] == eBL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmkj.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmkj = " << Fmkj
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eBU[n] - eBL[n] << ", " << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFInteriorMutationApprox(x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[n]) + log(eBL[n]) - log(eCU)));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    /// Add on j contributions
    int lindex = (!(x < 1.0) ? m : 0);
    for (int j = 0; j <= k && !mkj_found; ++j) {
      Amkj = computeA(m, k, lindex, j, x, z);
      if (o.debug > 2)
        cerr << "Adding to currsums with A(n,m,k,j) = A(" << n << "," << m
             << "," << k << "," << j << ") = " << Amkj << endl;

      if (Amkj * dL[n] > 1.0 || Amkj * dU[n] < 0.0 || eAU[n] < 0.0 ||
          eBU[n] < 0.0 || eAL[n] > 1.0 || eBL[n] > 1.0 || eCU < 0.0) {
        cerr << "Numerical error detected: (m,k,j) = " << m << "," << k << ","
             << j << "), Amkj = " << Amkj << ", dL[" << n << "] = " << dL[n]
             << ", dU[" << n << "] = " << dU[n] << ", eAU[" << n
             << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
             << ",  eBU[" << n << "] = " << eBU[n] << ", eBL[" << n
             << "] = " << eBL[n] << ", eCU = " << eCU << ", eCL = " << eCL
             << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
             << "," << z << "," << s << "," << t << ") ..." << endl;
        return DrawBridgePMFInteriorMutationApprox(x, z, s, t, o, gen);
      }

      currsumL += Amkj * dL[n];
      currsumU += Amkj * dU[n];
      bool decision_on_mkj_made = (currsumL > u || currsumU < u);

      if (currsumL > currsumU) {
        cerr << "Error: currsumL = " << currsumL << " > " << currsumU
             << " = currsumU (n,m,k,j) = (" << n << "," << m << "," << k << ","
             << j << ")." << endl;
        exit(1);
      }

      while (!decision_on_mkj_made)  /// Refine upper and lower bounds
      {
        double100 currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double100> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
          ++v_used[i];
          double100 newcoefficientU =
              exp(Getlogakm<double100>(mi + 2 * v_used[i], mi) +
                  static_cast<double100>(-(mi + 2 * v_used[i]) *
                                         (mi + 2 * v_used[i] + theta - 1) * s /
                                         2.0));
          double100 newcoefficientL =
              exp(Getlogakm<double100>(mi + 2 * v_used[i] + 1, mi) +
                  static_cast<double100>(-(mi + 2 * v_used[i] + 1) *
                                         (mi + 2 * v_used[i] + 1 + theta - 1) *
                                         s / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          newcoefficientU =
              exp(Getlogakm<double100>(ki + 2 * v_used[i], ki) +
                  static_cast<double100>(-(ki + 2 * v_used[i]) *
                                         (ki + 2 * v_used[i] + theta - 1) *
                                         (t - s) / 2.0));
          newcoefficientL =
              exp(Getlogakm<double100>(ki + 2 * v_used[i] + 1, ki) +
                  static_cast<double100>(-(ki + 2 * v_used[i] + 1) *
                                         (ki + 2 * v_used[i] + 1 + theta - 1) *
                                         (t - s) / 2.0));

          eBU[i] = eBL[i] + newcoefficientU;
          eBL[i] = eBU[i] - newcoefficientL;

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU - GetdBridgeInterior(eCvec, 2 * v_used[i] + 1, x, z, t);
            eCU = eCL + GetdBridgeInterior(eCvec, 2 * v_used[i] + 2, x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                             : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
          dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                       ? static_cast<double100>(0.0)
                       : exp(log(eAL[i]) + log(eBL[i]) - log(eCU)));

          if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
              static_cast<double>(eAU[i]) < 0.0 ||
              static_cast<double>(eBU[i]) < 0.0 ||
              static_cast<double>(eAL[i]) > 1.0 ||
              static_cast<double>(eBL[i]) > 1.0 ||
              static_cast<double>(eCU) < 0.0) {
            cerr << "Numerical error detected: (m,k,j) = " << mi << "," << ki
                 << ", *, *), dL[" << i << "] = " << dL[i] << ", dU[" << i
                 << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
                 << ", eAL[" << i << "] = " << eAL[i] << ",  eBU[" << i
                 << "] = " << eBU[i] << ", eBL[" << i << "] = " << eBL[i]
                 << ", eCU = " << eCU << ", eCL = " << eCL
                 << ". Resorting to G1984-style approximation (x,z,s,t) = ("
                 << x << "," << z << "," << s << "," << t << ") ..." << endl;
            return DrawBridgePMFInteriorMutationApprox(x, z, s, t, o, gen);
          }

          int liindex = (!(x < 1.0) ? mi : 0);
          for (int j2 = 0; j2 <= ki; ++j2) {
            Amkj = computeA(mi, ki, liindex, j2, x, z);
            currsumL += Amkj * dL[i];
            currsumU += Amkj * dU[i];
            if (o.debug > 3)
              cerr << "Recomputing currsums with A(n,m,k,j) = A(" << i << ","
                   << mi << "," << ki << "," << j2 << ") = " << Amkj << endl;
          }
        }

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k,j) = (" << n << "," << m << "," << k
               << "," << j << ")." << endl;
          exit(1);
        }

        if (o.debug > 2) {
          cerr << "\ndL ";
          printVec(dL, cerr);
          cerr << "\ndU ";
          printVec(dU, cerr);
        }

        if (currsumLold > currsumL) {
          // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,j) = (" << n << "," << m
                    << "," << k << "," << j << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
          std::cout << "Error: currsumUold = " << currsumUold << " < "
                    << currsumU << " = currsumU (n,m,k,j) = (" << n << "," << m
                    << "," << k << "," << j << ")." << endl;
          // exit(1);
        }

        decision_on_mkj_made = (currsumL > u || currsumU < u);
      }

      mkj_found = (currsumL > u);

      if (mkj_found) {
        if (!(x > 0.0)) {
          curr_mk[n][2] = 0;  /// Sets l = 0
          curr_mk[n][3] = j;
        } else {
          curr_mk[n][2] = m;  /// Sets l = m
          curr_mk[n][3] = j;
        }
      }
    }
  }

  int coeffcount = 0;
  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,j: Returned (m,k,j) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << "," << curr_mk[n][3] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFInteriorMutationOneQApprox(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting at a boundary point
               /// but ending in the interior of (0,1) but one time increment is
               /// below the threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double100 sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkj_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec, eAL, eAU, dL, dU, currsumStore;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(GetdBridgeInterior(eCvec, 0, x, z, t)), qApprox = 0.0;
  vector<int> v_used;
  double100 Amkj, runningMax = 1.0e-300;
  int n = -1, Fmkj = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to be
  /// needed
  double mmode = GriffithsParas(s).first, vm = GriffithsParas(s).second,
         kmode = GriffithsParas(t - s).first, vk = GriffithsParas(t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double100> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    int mork = (ind1 ? m : k), morkapprox = (mork == m ? k : m);
    if (m <= mlimU && m >= mlimL && k <= klimU && k >= klimL) counter++;
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(mork, Csorts);
    if (o.debug > 2)
      cerr << "(F1,F3,F4,F5) = (" << F1 << "," << F3 << "," << F4 << "," << F5
           << ")" << endl;
    Fmkj = max(F1, F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmkj << endl;
    }

    while (2 * v < Fmkj) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(mork + 2 * v, mork) +
              static_cast<double100>(-(mork + 2 * v) *
                                     (mork + 2 * v + theta - 1) * sorts / 2.0));
      double100 newcoefficientL = exp(
          Getlogakm<double100>(mork + 2 * v + 1, mork) +
          static_cast<double100>(-(mork + 2 * v + 1) *
                                 (mork + 2 * v + 1 + theta - 1) * sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeInterior(eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeInterior(eCvec, 2 * v + 2, x, z, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmkj.
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmkj = " << Fmkj
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFInteriorMutationApprox(x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double100>(0.0)
                          : exp(log(eAL[n]) + log(qApprox) - log(eCU)));
    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    /// Add on j contributions
    int lindex = (!(x < 1.0) ? m : 0);
    for (int j = 0; j <= k && !mkj_found; ++j) {
      Amkj = computeA(m, k, lindex, j, x, z);
      if (o.debug > 2)
        cerr << "Adding to currsums with A(n,m,k,j) = A(" << n << "," << m
             << "," << k << "," << j << ") = " << Amkj << endl;

      if (Amkj * dL[n] > 1.0 || Amkj * dU[n] < 0.0 || eAU[n] < 0.0 ||
          eAL[n] > 1.0 || eCU < 0.0) {
        cerr << "Numerical error detected: (m,k,j) = " << m << "," << k << ","
             << j << "), Amkj = " << Amkj << ", dL[" << n << "] = " << dL[n]
             << ", dU[" << n << "] = " << dU[n] << ", eAU[" << n
             << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
             << ", eCU = " << eCU << ", eCL = " << eCL
             << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
             << "," << z << "," << s << "," << t << ") ..." << endl;
        return DrawBridgePMFInteriorMutationApprox(x, z, s, t, o, gen);
      }

      currsumL += Amkj * dL[n];
      currsumU += Amkj * dU[n];

      currsumStore.push_back(log(Amkj) + log(dL[n]));
      runningMax = max(runningMax, currsumStore.back());
      mkljStore.resize(mkljStore.size() + 4, -1);
      mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
      mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
      mkljStore[2 + (4 * (currsumStore.size() - 1))] = lindex;
      mkljStore[3 + (4 * (currsumStore.size() - 1))] = j;

      bool decision_on_mkj_made = (currsumL > u || currsumU < u);

      if (currsumL > currsumU) {
        cerr << "Error: currsumL = " << currsumL << " > " << currsumU
             << " = currsumU (n,m,k,j) = (" << n << "," << m << "," << k << ","
             << j << ")." << endl;
        exit(1);
      }

      while (!decision_on_mkj_made)  /// Refine upper and lower bounds
      {
        double100 currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double100> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
              morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
          ++v_used[i];
          double100 newcoefficientU =
              exp(Getlogakm<double100>(morki + 2 * v_used[i], morki) +
                  static_cast<double100>(-(morki + 2 * v_used[i]) *
                                         (morki + 2 * v_used[i] + theta - 1) *
                                         sorts / 2.0));
          double100 newcoefficientL =
              exp(Getlogakm<double100>(morki + 2 * v_used[i] + 1, morki) +
                  static_cast<double100>(
                      -(morki + 2 * v_used[i] + 1) *
                      (morki + 2 * v_used[i] + 1 + theta - 1) * sorts / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          qApprox = DiscretisedNormCDF(morkapproxi, sortsapprox);

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU - GetdBridgeInterior(eCvec, 2 * v_used[i] + 1, x, z, t);
            eCU = eCL + GetdBridgeInterior(eCvec, 2 * v_used[i] + 2, x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                             : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
          dL[i] = (eAL[i] < 0.0 ? static_cast<double100>(0.0)
                                : exp(log(eAL[i]) + log(qApprox) - log(eCU)));

          if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
              static_cast<double>(eAU[i]) < 0.0 ||
              static_cast<double>(eAL[i]) > 1.0 ||
              static_cast<double>(eCU) < 0.0) {
            cerr << "Numerical error detected: (m,k,l,j) = " << mi << "," << ki
                 << ", *, *), dL[" << i << "] = " << dL[i] << ", dU[" << i
                 << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
                 << ", eAL[" << i << "] = " << eAL[i] << ", eCU = " << eCU
                 << ", eCL = " << eCL
                 << ". Resorting to G1984-style approximation (x,z,s,t) = ("
                 << x << "," << z << "," << s << "," << t << ") ..." << endl;
            return DrawBridgePMFInteriorMutationApprox(x, z, s, t, o, gen);
          }

          int liindex = (!(x < 1.0) ? mi : 0);
          for (int j2 = 0; j2 <= ki; ++j2) {
            Amkj = computeA(mi, ki, liindex, j2, x, z);
            currsumL += Amkj * dL[i];
            currsumU += Amkj * dU[i];
            currsumStore[i] = log(Amkj) + log(dL[i]);
            runningMax = max(runningMax, currsumStore[i]);
            if (o.debug > 3)
              cerr << "Recomputing currsums with A(n,m,k,j) = A(" << i << ","
                   << mi << "," << ki << "," << j2 << ") = " << Amkj << endl;
          }
        }

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k,j) = (" << n << "," << m << "," << k
               << "," << j << ")." << endl;
          exit(1);
        }

        if (o.debug > 2) {
          cerr << "\ndL ";
          printVec(dL, cerr);
          cerr << "\ndU ";
          printVec(dU, cerr);
        }

        if (currsumLold > currsumL) {
          // cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,l,j) = (" << n << ","
                    << m << "," << k << "," << j << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
          std::cout << "Error: currsumUold = " << currsumUold << " < "
                    << currsumU << " = currsumU (n,m,k,l,j) = (" << n << ","
                    << m << "," << k << "," << j << ")." << endl;
          // exit(1);
        }

        decision_on_mkj_made = (currsumL > u || currsumU < u);
      }

      mkj_found = (currsumL > u);

      if (mkj_found) {
        if (!(x > 0.0)) {
          curr_mk[n][2] = 0;  /// Sets l = 0
          curr_mk[n][3] = j;
        } else {
          curr_mk[n][2] = m;  /// Sets l = m
          curr_mk[n][3] = j;
        }
      }
    }

    if (counter == totalpts)  /// Gaussian approximation leads to currsum
                              /// summing to < 1.0, so we renormalise and sample
    {
      LogSumExp(currsumStore, runningMax);
      double100 sum = 0.0;
      int ind = 0;

      bool found = false;

      while (!found) {
        sum += currsumStore[ind];
        if (sum > u) {
          found = true;
        }
        if (ind == static_cast<int>(currsumStore.size() - 1)) {
          found = true;
        }
        ind++;
      }

      vector<int> returnvec;

      returnvec.push_back(mkljStore[0 + (4 * ind)]);
      returnvec.push_back(mkljStore[1 + (4 * ind)]);
      returnvec.push_back(mkljStore[2 + (4 * ind)]);
      returnvec.push_back(mkljStore[3 + (4 * ind)]);

      return returnvec;
    }
  }

  int coeffcount = 0;
  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);
  if (o.debug > 2) {
    cerr << "p_m,k,j: Returned (m,k,j) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << "," << curr_mk[n][3] << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFInteriorMutationApprox(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting at a boundary point
               /// but ending in the interior of (0,1) with both time increments
               /// below the threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  vector<int> returnvec;
  vector<double100> currsumStore;
  vector<int> mkljStore;

  /// Compute denominator
  double100 eC = 0.0, eCInc = 1.0, eCOldInc = 1.0;
  int Dflip = 1, Djm = 0, Djp = 0;
  int dmode = static_cast<int>(ceil(GriffithsParas(t).first)), d = dmode;

  while (max(eCOldInc, eCInc) > 0.0 || (!(eC > 0.0))) {
    eCOldInc = eCInc;
    double100 para1 = (!(x > 0.0) ? static_cast<double100>(thetaP[0])
                                  : static_cast<double100>(thetaP[0] + d));
    double100 para2 = (!(x > 0.0) ? static_cast<double100>(thetaP[1] + d)
                                  : static_cast<double100>(thetaP[1]));
    double100 zcont = (!(x > 0.0) ? static_cast<double100>(d) * log(1.0 - z)
                                  : static_cast<double100>(d) * log(z));

    eCInc = exp(log(QmApprox(d, t, o)) + zcont -
                boost::math::lgamma(static_cast<double100>(para1)) -
                boost::math::lgamma(static_cast<double100>(para2)) +
                boost::math::lgamma(static_cast<double100>(para1 + para2)));
    eC += eCInc;

    if (Dflip == -1 &&
        (dmode - Djm - 1 > 0))  /// Mechanism to explore either side around mode
    {
      Djm++;
      d = dmode - Djm;
    } else {
      Djp++;
      d = dmode + Djp;
    }
    Dflip *= -1;
  }

  vector<int> modeGuess = mkjModeFinder(
      x, z, s, t, o);  /// Get a guess to location of mode over (m,k,j)
  int mMode = modeGuess[0], kMode = modeGuess[1], jMode = modeGuess[2];

  boost::random::uniform_01<double100>
      U01;  /// Use these guesses & eC to set a suitable threshold for
  /// subsequent computations
  double100 currsum = 0.0, u = U01(gen),
            threshold = exp(mkjModeFinder_Evaluator(mMode, kMode, jMode, x, z,
                                                    s, t, o) -
                            log(eC)) *
                        1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0, kFlip = 1, kD = 0, kU = 0, l;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double100 constContr =
      (static_cast<double100>(thetaP[0] - 1.0)) * log(z) +
      (static_cast<double100>(thetaP[1] - 1.0)) * log(1.0 - z) +
      (!(x > 0.0) ? -boost::math::lgamma(static_cast<double100>(thetaP[0]))
                  : -boost::math::lgamma(static_cast<double100>(thetaP[1])));
  double100 mContr_D =
      boost::math::lgamma(static_cast<double100>(thetaP[0] + thetaP[1] + m)) +
      (!(x > 0.0)
           ? -boost::math::lgamma(static_cast<double100>(thetaP[1] + m))
           : -boost::math::lgamma(static_cast<double100>(thetaP[0] + m)));
  double100 mContr_U = mContr_D, mContr,
            runningMax = -1.0e100;  /// Computing m contributions

  while (!mSwitch) {
    double100 qm = QmApprox(m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to very
                      /// small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    l = (!(x > 0.0) ? 0 : m);
    if (m != mMode) {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m - 1))) +
            (!(x > 0.0) ? -log(static_cast<double100>(thetaP[1] + (m - 1)))
                        : -log(static_cast<double100>(thetaP[0] + (m - 1))));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D +=
            -log(static_cast<double100>(thetaP[0] + thetaP[1] + (m + 1) - 1)) +
            (!(x > 0.0) ? log(static_cast<double100>(thetaP[1] + (m + 1) - 1))
                        : log(static_cast<double100>(thetaP[0] + (m + 1) - 1)));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode;
    kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Computing k contributions
    double100 kContr_D =
                  (boost::math::lgamma(
                       static_cast<double100>(thetaP[0] + thetaP[1] + k)) -
                   boost::math::lgamma(
                       static_cast<double100>(thetaP[0] + thetaP[1] + m + k))) +
                  static_cast<double100>(k) * log(1.0 - z),
              kContr_U = kContr_D, kContr;

    while (!kSwitch) {
      double100 qk = QmApprox(k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      if (k != kMode) {
        if (kU > kD) {
          kContr_U +=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k - 1))) -
              log(static_cast<double100>(thetaP[0] + thetaP[1] + m + (k - 1))) +
              log(1.0 - z);
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D +=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + m + (k + 1) -
                                         1)) -
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k + 1) - 1)) -
              log(1.0 - z);
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      int jFlip = 1, newjMode = min(k, jMode), j = newjMode, jU = 0,
          jD = 0;  /// Need to redefine jMode as k might be too small!
      bool jSwitch = false, jDownSwitch = false, jUpSwitch = false;

      double100 jContr_D =
          LogBinomialCoefficientCalculator(k, j) +
          (!(x > 0.0) ? boost::math::lgamma(
                            static_cast<double100>(thetaP[1] + m + k - j)) -
                            boost::math::lgamma(
                                static_cast<double100>(thetaP[1] + k - j))
                      : boost::math::lgamma(
                            static_cast<double100>(thetaP[0] + m + j)) -
                            boost::math::lgamma(
                                static_cast<double100>(thetaP[0] + j))) +
          static_cast<double100>(j) * log(z) +
          static_cast<double100>(-j) * log(1.0 - z);
      double100 jContr_U = jContr_D, jContr;

      while (!jSwitch) {
        if (j != newjMode) {
          if (jU > jD) {
            jContr_U +=
                log(static_cast<double100>(k - (j - 1))) -
                log(static_cast<double100>((j - 1) + 1)) + log(z) -
                log(1.0 - z) +
                (!(x > 0.0)
                     ? log(static_cast<double100>(thetaP[1] + k - (j - 1) -
                                                  1)) -
                           log(static_cast<double100>(thetaP[1] + m + k -
                                                      (j - 1) - 1))
                     : log(static_cast<double100>(thetaP[0] + m + (j - 1))) -
                           log(static_cast<double100>(thetaP[0] + (j - 1))));
            jContr = jContr_U;
          } else {
            jContr_D +=
                log(static_cast<double100>(j + 1)) -
                log(static_cast<double100>(k - (j + 1) + 1)) + log(1.0 - z) -
                log(z) +
                (!(x > 0.0)
                     ? log(static_cast<double100>(thetaP[1] + m + k -
                                                  (j + 1))) -
                           log(static_cast<double100>(thetaP[1] + k - (j + 1)))
                     : log(static_cast<double100>(thetaP[0] + (j + 1) - 1)) -
                           log(static_cast<double100>(thetaP[0] + m + (j + 1) -
                                                      1)));
            jContr = jContr_D;
          }
        } else {
          jContr = jContr_U;
        }
        double100 currsum_inc = constContr + mContr + kContr + jContr - log(eC);
        runningMax =
            max(currsum_inc, runningMax);  /// Running max of log probabilities
        /// for use in log-sum-exp trick
        currsum += exp(currsum_inc);
        currsumStore.push_back(currsum_inc);

        mkljStore.resize(mkljStore.size() + 4, -1);
        mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
        mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
        mkljStore[2 + (4 * (currsumStore.size() - 1))] = l;
        mkljStore[3 + (4 * (currsumStore.size() - 1))] = j;

        if (!(jDownSwitch))  /// Switching mechanism for j
        {
          if (sgn(j - newjMode) <= 0) {
            jDownSwitch =
                ((exp(currsum_inc) < threshold) || (newjMode - jD - 1) < 0);
          }
        }

        if (!(jUpSwitch)) {
          if (sgn(j - newjMode) >= 0) {
            jUpSwitch =
                ((exp(currsum_inc) < threshold) || (newjMode + jU + 1) > k);
          }
        }

        jSwitch = (jDownSwitch && jUpSwitch);

        if (!jSwitch) {
          if ((jFlip == 1 && (newjMode + jU + 1 <= k)) ||
              (jDownSwitch && !(jUpSwitch))) {
            jU++;
            j = newjMode + jU;
            jFlip *= (jDownSwitch ? 1 : -1);
          } else if ((jFlip == -1 && (newjMode - jD - 1 >= 0)) ||
                     (jUpSwitch && !(jDownSwitch))) {
            jD++;
            j = newjMode - jD;
            jFlip *= (jUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch = (((jU == 0) && (jD == 0)) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = ((jU == 0) && (jD == 0));
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1) && (kMode - kD - 1 >= 0)) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Switching mechanism for m
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < 0));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= 0)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }

  LogSumExp(currsumStore, runningMax);
  double100 sum = 0.0;
  int index, ind = 0;

  vector<int> indexing(currsumStore.size(), 0);
  for (int i = 0; i != static_cast<int>(indexing.size()); i++) {
    indexing[i] = i;
  }
  sort(indexing.begin(), indexing.end(), [&](const int &a, const int &b) {
    return (currsumStore[a] > currsumStore[b]);
  });

  bool found = false;

  while (!found) {
    sum += currsumStore[indexing[ind]];
    if (sum > u) {
      index = indexing[ind];
      found = true;
    }
    if (ind == static_cast<int>(currsumStore.size() - 1)) {
      index = indexing[ind];
      found = true;
    }
    ind++;
  }

  returnvec.push_back(mkljStore[0 + (4 * index)]);
  returnvec.push_back(mkljStore[1 + (4 * index)]);
  returnvec.push_back(mkljStore[2 + (4 * index)]);
  returnvec.push_back(mkljStore[3 + (4 * index)]);

  return returnvec;
}

vector<int> WrightFisher::DrawBridgePMF(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending in the
               /// interior of (0,1) with both time increments large enough
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(Getd(eCvec, 0, x, z, t));
  vector<int> v_used;
  double100 Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double100> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F5;
  if (thetaP.empty()) {
    F5 = static_cast<int>(
        ceil(((1.0 - x) / (1.0 - z)) + (x / z) * (1.0 + z) * o.eps));
  } else {
    F5 = static_cast<int>(ceil(
        2.0 *
        (max((theta / thetaP[1]) * (1.0 - static_cast<double>(z)),
             (1.0 + theta) / (1.0 - static_cast<double>(z))) *
             (1.0 - static_cast<double>(x)) +
         (1.0 + (1.0 / static_cast<double>(z))) *
             (max(((static_cast<double>(z) * theta) + 1.0) / thetaP[0], 1.0)) *
             static_cast<double>(x)) /
        o.eps));
  }
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0)));
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    eBL.push_back(0.0);
    eBU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(m, Cs), F2 = computeC(k, Cts);
    if (o.debug > 2)
      cerr << "(F1,F2,F3,F4,F5) = (" << F1 << "," << F2 << "," << F3 << ","
           << F4 << "," << F5 << ")" << endl;
    Fmklj = max(max(F1, F2), F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmklj << endl;
    }

    while (2 * v < Fmklj || eAU[n] > 1.0 || eBU[n] > 1.0) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(m + 2 * v, m) +
              static_cast<double100>(-(m + 2 * v) * (m + 2 * v + theta - 1) *
                                     s / 2.0));
      double100 newcoefficientL =
          exp(Getlogakm<double100>(m + 2 * v + 1, m) +
              static_cast<double100>(-(m + 2 * v + 1) *
                                     (m + 2 * v + 1 + theta - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double100>(k + 2 * v, k) +
              static_cast<double100>(-(k + 2 * v) * (k + 2 * v + theta - 1) *
                                     (t - s) / 2.0));
      newcoefficientL = exp(Getlogakm<double100>(k + 2 * v + 1, k) +
                            static_cast<double100>(-(k + 2 * v + 1) *
                                                   (k + 2 * v + 1 + theta - 1) *
                                                   (t - s) / 2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - Getd(eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + Getd(eCvec, 2 * v + 2, x, z, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] && eBU[n] == eBL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmklj
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmklj = " << Fmklj
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eBU[n] - eBL[n] << ", " << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFG1984(x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double100>(0.0)
                 : exp(log(eAL[n]) + log(eBL[n]) - log(eCU)));

    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    for (int l = 0; l <= m && !mklj_found; ++l) {
      for (int j = 0; j <= k && !mklj_found; ++j) {
        Amklj = computeA(m, k, l, j, x, z);
        if (o.debug > 2)
          cerr << "Adding to currsums with A(n,m,k,l,j) = A(" << n << "," << m
               << "," << k << "," << l << "," << j << ") = " << Amklj << endl;

        if (Amklj * dL[n] > 1.0 || Amklj * dU[n] < 0.0 || eAU[n] < 0.0 ||
            eBU[n] < 0.0 || eAL[n] > 1.0 || eBL[n] > 1.0 || eCU < 0.0) {
          cerr << "Numerical error detected: (m,k,l,j) = " << m << "," << k
               << "," << l << "," << j << "), Amklj = " << Amklj << ", dL[" << n
               << "] = " << dL[n] << ", dU[" << n << "] = " << dU[n] << ", eAU["
               << n << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
               << ",  eBU[" << n << "] = " << eBU[n] << ", eBL[" << n
               << "] = " << eBL[n] << ", eCU = " << eCU << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
               << "," << z << "," << s << "," << t << ") ..." << endl;
          return DrawBridgePMFG1984(x, z, s, t, o, gen);
        }

        currsumL += Amklj * dL[n];
        currsumU += Amklj * dU[n];
        bool decision_on_mklj_made = (currsumL > u || currsumU < u);

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k,l,j) = (" << n << "," << m << "," << k
               << "," << l << "," << j << ")." << endl;
          exit(1);
        }

        while (!decision_on_mklj_made)  /// Refine upper and lower bounds
        {
          double100 currsumLold = currsumL, currsumUold = currsumU;
          currsumL = 0.0;
          currsumU = 0.0;

          const vector<double100> dUold(dU), dLold(dL);
          for (int i = 0; i <= n; ++i) {
            int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
            ++v_used[i];
            double100 newcoefficientU =
                          exp(Getlogakm<double100>(mi + 2 * v_used[i], mi) +
                              static_cast<double100>(
                                  -(mi + 2 * v_used[i]) *
                                  (mi + 2 * v_used[i] + theta - 1) * s / 2.0)),
                      newcoefficientL = exp(
                          Getlogakm<double100>(mi + 2 * v_used[i] + 1, mi) +
                          static_cast<double100>(
                              -(mi + 2 * v_used[i] + 1) *
                              (mi + 2 * v_used[i] + 1 + theta - 1) * s / 2.0));

            eAU[i] = eAL[i] + newcoefficientU;
            eAL[i] = eAU[i] - newcoefficientL;

            newcoefficientU =
                exp(Getlogakm<double100>(ki + 2 * v_used[i], ki) +
                    static_cast<double100>(-(ki + 2 * v_used[i]) *
                                           (ki + 2 * v_used[i] + theta - 1) *
                                           (t - s) / 2.0));
            newcoefficientL =
                exp(Getlogakm<double100>(ki + 2 * v_used[i] + 1, ki) +
                    static_cast<double100>(
                        -(ki + 2 * v_used[i] + 1) *
                        (ki + 2 * v_used[i] + 1 + theta - 1) * (t - s) / 2.0));

            eBU[i] = eBL[i] + newcoefficientU;
            eBL[i] = eBU[i] - newcoefficientL;

            if (2 * v_used[i] + 2 > eCindex_computed) {
              assert(2 * v_used[i] == eCindex_computed);
              eCL = eCU - Getd(eCvec, 2 * v_used[i] + 1, x, z, t);
              eCU = eCL + Getd(eCvec, 2 * v_used[i] + 2, x, z, t);
              eCindex_computed += 2;
            }

            dU[i] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                               : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
            dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                         ? static_cast<double100>(0.0)
                         : exp(log(eAL[i]) + log(eBL[i]) - log(eCU)));
            if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
                static_cast<double>(eAU[i]) < 0.0 ||
                static_cast<double>(eBU[i]) < 0.0 ||
                static_cast<double>(eAL[i]) > 1.0 ||
                static_cast<double>(eBL[i]) > 1.0 ||
                static_cast<double>(eCU) < 0.0) {
              cerr << "Numerical error detected: (m,k,l,j) = " << mi << ","
                   << ki << ", *, *), dL[" << i << "] = " << dL[i] << ", dU["
                   << i << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
                   << ", eAL[" << i << "] = " << eAL[i] << ",  eBU[" << i
                   << "] = " << eBU[i] << ", eBL[" << i << "] = " << eBL[i]
                   << ", eCU = " << eCU << ", eCL = " << eCL
                   << ". Resorting to G1984-style approximation (x,z,s,t) = ("
                   << x << "," << z << "," << s << "," << t << ") ..." << endl;
              return DrawBridgePMFG1984(x, z, s, t, o, gen);
            }

            int l_upper = (i == n ? l : mi);
            for (int l2 = 0; l2 <= l_upper; ++l2) {
              int j_upper = ((i == n && l2 == l_upper) ? j : ki);
              for (int j2 = 0; j2 <= j_upper; ++j2) {
                Amklj = computeA(mi, ki, l2, j2, x, z);
                currsumL += Amklj * dL[i];
                currsumU += Amklj * dU[i];
                if (o.debug > 3)
                  cerr << "Recomputing currsums with A(n,m,k,l,j) = A(" << i
                       << "," << mi << "," << ki << "," << l2 << "," << j2
                       << ") = " << Amklj << endl;
              }
            }

            if (currsumL > currsumU) {
              cerr << "Error: currsumL = " << currsumL << " > " << currsumU
                   << " = currsumU (n,m,k,l,j) = (" << n << "," << m << "," << k
                   << "," << l << "," << j << ")." << endl;
              exit(1);
            }
          }

          if (o.debug > 2) {
            cerr << "\ndL ";
            printVec(dL, cerr);
            cerr << "\ndU ";
            printVec(dU, cerr);
          }

          if (currsumLold > currsumL) {
            // cerr << "Error: currsumLold = " << currsumLold << " > " <<
            // currsumL
            std::cout << "Error: currsumLold = " << currsumLold << " > "
                      << currsumL << " = currsumL (n,m,k,l,j) = (" << n << ","
                      << m << "," << k << "," << l << "," << j << ")." << endl;
            // exit(1);
          }
          if (currsumUold < currsumU) {
            // cerr << "Error: currsumUold = " << currsumUold << " < " <<
            // currsumU
            std::cout << "Error: currsumUold = " << currsumUold << " < "
                      << currsumU << " = currsumU (n,m,k,l,j) = (" << n << ","
                      << m << "," << k << "," << l << "," << j << ")." << endl;
            // exit(1);
          }

          decision_on_mklj_made = (currsumL > u || currsumU < u);
        }

        mklj_found = (currsumL > u);
        if (mklj_found) {
          curr_mk[n][2] = l;
          curr_mk[n][3] = j;
        }
      }
    }
  }

  int coeffcount = 0;
  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,l,j: Returned (m,k,l,j) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << "," << curr_mk[n][2] << "," << curr_mk[n][3]
         << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFOneQApprox(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of a bridge starting
                                  /// and ending in the interior of (0,1), but
                                  /// one time increment is below the threshold
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double100 sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double100> U01;
  double100 u = U01(gen);
  vector<double100> eCvec, eAL, eAU, dL, dU, currsumStore;
  double100 currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
            eCU(Getd(eCvec, 0, x, z, t)), qApprox = 0.0, runningMax = 1.0e-300;
  vector<int> v_used;
  double100 Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to be
  /// needed
  double mmode = GriffithsParas(s).first, vm = GriffithsParas(s).second,
         kmode = GriffithsParas(t - s).first, vk = GriffithsParas(t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double100> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F5;
  if (thetaP.empty()) {
    F5 = static_cast<int>(
        ceil(((1.0 - x) / (1.0 - z)) + (x / z) * (1.0 + z) * o.eps));
  } else {
    F5 = static_cast<int>(ceil(
        2.0 *
        (max((theta / thetaP[1]) * (1.0 - static_cast<double>(z)),
             (1.0 + theta) / (1.0 - static_cast<double>(z))) *
             (1.0 - static_cast<double>(x)) +
         (1.0 + (1.0 / static_cast<double>(z))) *
             (max(((static_cast<double>(z) * theta) + 1.0) / thetaP[0], 1.0)) *
             static_cast<double>(x)) /
        o.eps));
  }
  int F3 = computeE(Ct),
      F4 = static_cast<int>(
          ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta + 1.0) / 2.0)));
  while ((theta + 2 * F4 + 1) * exp(-(2 * F4 + theta) * t / 2.0) >= 1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(curr_mk.back(), s, t);
    int &m = curr_mk[n][0];
    int &k = curr_mk[n][1];
    int mork = (ind1 ? m : k), morkapprox = (mork == m ? k : m);
    if (m <= mlimU && m >= mlimL && k <= klimU && k >= klimL) counter++;
    if (o.debug > 2)
      cerr << "New n = " << n << ", (m,k) = (" << m << ", " << k << ")" << endl;

    /// Add the new n to the mix
    eAL.push_back(0.0);
    eAU.push_back(0.0);
    dL.push_back(0.0);
    dU.push_back(0.0);
    int F1 = computeC(mork, Csorts);
    if (o.debug > 2)
      cerr << "(F1,F3,F4,F5) = (" << F1 << "," << F3 << "," << F4 << "," << F5
           << ")" << endl;
    Fmklj = max(F1, F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmklj << endl;
    }

    while (2 * v < Fmklj) {
      ++v;
      double100 newcoefficientU =
          exp(Getlogakm<double100>(mork + 2 * v, mork) +
              static_cast<double100>(-(mork + 2 * v) *
                                     (mork + 2 * v + theta - 1) * sorts / 2.0));
      double100 newcoefficientL = exp(
          Getlogakm<double100>(mork + 2 * v + 1, mork) +
          static_cast<double100>(-(mork + 2 * v + 1) *
                                 (mork + 2 * v + 1 + theta - 1) * sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - Getd(eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + Getd(eCvec, 2 * v + 2, x, z, t);
        eCindex_computed += 2;
      }

      if (eAU[n] == eAL[n] &&
          eCL == eCU)  /// ...then we have lost precision before reaching Fmklj
      {
        if (o.debug > 2) {
          cerr << "Abandoning loop for n = " << n << ", Fmklj = " << Fmklj
               << " at v = " << v << endl;
          cerr << "Leftovers: " << setprecision(50) << eAU[n] - eAL[n] << ", "
               << eCU - eCL << endl;
        }

        if (eAL[n] < 0.0 || eCL < 0.0) {
          cerr << "Numerical error detected: (m,k) = " << m << "," << k
               << ", eAU[" << n << "] = " << eAU[n] << ", eAL[" << n
               << "] = " << eAL[n] << ",  eCU[" << n << "] = " << eCU
               << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (s,t) = (" << s
               << "," << t << ") ..." << endl;
          return DrawBridgePMFG1984(x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double100>(0.0)
                          : exp(log(eAL[n]) + log(qApprox) - log(eCU)));

    v_used.push_back(v);

    if (o.debug > 2) {
      cerr << "\nn   ";
      for (int k = 0; k <= n; ++k) cerr << k << " ";
      cerr << "\ndL ";
      for (int k = 0; k <= n; ++k) cerr << dL[k] << " ";
      cerr << "\ndU ";
      for (int k = 0; k <= n; ++k) cerr << dU[k] << " ";
      cerr << endl;
    }

    for (int l = 0; l <= m && !mklj_found; ++l) {
      for (int j = 0; j <= k && !mklj_found; ++j) {
        Amklj = computeA(m, k, l, j, x, z);
        if (o.debug > 2)
          cerr << "Adding to currsums with A(n,m,k,l,j) = A(" << n << "," << m
               << "," << k << "," << l << "," << j
               << ") = " << endl;  // Amklj[Akey] << endl;

        if (Amklj * dL[n] > 1.0 || Amklj * dU[n] < 0.0 || eAU[n] < 0.0 ||
            eAL[n] > 1.0 || eCU < 0.0) {
          cerr << "Numerical error detected: (m,k,l,j) = " << m << "," << k
               << "," << l << "," << j << "), Amklj = " << Amklj << ", dL[" << n
               << "] = " << dL[n] << ", dU[" << n << "] = " << dU[n] << ", eAU["
               << n << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
               << ", eCU = " << eCU << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
               << "," << z << "," << s << "," << t << ") ..." << endl;
          return DrawBridgePMFG1984(x, z, s, t, o, gen);
        }

        currsumL += Amklj * dL[n];
        currsumU += Amklj * dU[n];

        currsumStore.push_back(log(Amklj) + log(dL[n]));
        runningMax = max(runningMax, currsumStore.back());
        mkljStore.resize(mkljStore.size() + 4, -1);
        mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
        mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
        mkljStore[2 + (4 * (currsumStore.size() - 1))] = l;
        mkljStore[3 + (4 * (currsumStore.size() - 1))] = j;

        bool decision_on_mklj_made = (currsumL > u || currsumU < u);

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k,l,j) = (" << n << "," << m << "," << k
               << "," << l << "," << j << ")." << endl;
          exit(1);
        }

        while (!decision_on_mklj_made)  /// Refine upper and lower bounds
        {
          double100 currsumLold = currsumL, currsumUold = currsumU;
          currsumL = 0.0;
          currsumU = 0.0;

          const vector<double100> dUold(dU), dLold(dL);
          for (int i = 0; i <= n; ++i) {
            int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
                morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
            ;
            ++v_used[i];
            double100
                newcoefficientU =
                    exp(Getlogakm<double100>(morki + 2 * v_used[i], morki) +
                        static_cast<double100>(
                            -(morki + 2 * v_used[i]) *
                            (morki + 2 * v_used[i] + theta - 1) * sorts / 2.0)),
                newcoefficientL = exp(
                    Getlogakm<double100>(morki + 2 * v_used[i] + 1, morki) +
                    static_cast<double100>(
                        -(morki + 2 * v_used[i] + 1) *
                        (morki + 2 * v_used[i] + 1 + theta - 1) * sorts / 2.0));

            eAU[i] = eAL[i] + newcoefficientU;
            eAL[i] = eAU[i] - newcoefficientL;

            qApprox = DiscretisedNormCDF(morkapproxi, sortsapprox);

            if (2 * v_used[i] + 2 > eCindex_computed) {
              assert(2 * v_used[i] == eCindex_computed);
              eCL = eCU - Getd(eCvec, 2 * v_used[i] + 1, x, z, t);
              eCU = eCL + Getd(eCvec, 2 * v_used[i] + 2, x, z, t);
              eCindex_computed += 2;
            }

            dU[i] = (eCL < 0.0 ? static_cast<double100>(nan(""))
                               : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
            dL[i] = (eAL[i] < 0.0 ? static_cast<double100>(0.0)
                                  : exp(log(eAL[i]) + log(qApprox) - log(eCU)));
            if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
                static_cast<double>(eAU[i]) < 0.0 ||
                static_cast<double>(eAL[i]) > 1.0 ||
                static_cast<double>(eCU) < 0.0) {
              cerr << "Numerical error detected: (m,k,l,j) = " << mi << ","
                   << ki << ", *, *), dL[" << i << "] = " << dL[i] << ", dU["
                   << i << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
                   << ", eAL[" << i << "] = " << eAL[i] << ", eCU = " << eCU
                   << ", eCL = " << eCL
                   << ". Resorting to G1984-style approximation (x,z,s,t) = ("
                   << x << "," << z << "," << s << "," << t << ") ..." << endl;
              return DrawBridgePMFG1984(x, z, s, t, o, gen);
            }

            int l_upper = (i == n ? l : mi);
            for (int l2 = 0; l2 <= l_upper; ++l2) {
              int j_upper = ((i == n && l2 == l_upper) ? j : ki);
              for (int j2 = 0; j2 <= j_upper; ++j2) {
                Amklj = computeA(mi, ki, l2, j2, x, z);
                currsumL += Amklj * dL[i];
                currsumU += Amklj * dU[i];
                currsumStore[i] = log(Amklj) + log(dL[i]);
                runningMax = max(runningMax, currsumStore[i]);
                if (o.debug > 3)
                  cerr << "Recomputing currsums with A(n,m,k,l,j) = A(" << i
                       << "," << mi << "," << ki << "," << l2 << "," << j2
                       << ") = " << Amklj << endl;
              }
            }

            if (currsumL > currsumU) {
              cerr << "Error: currsumL = " << currsumL << " > " << currsumU
                   << " = currsumU (n,m,k,l,j) = (" << n << "," << m << "," << k
                   << "," << l << "," << j << ")." << endl;
              exit(1);
            }
          }

          if (o.debug > 2) {
            cerr << "\ndL ";
            printVec(dL, cerr);
            cerr << "\ndU ";
            printVec(dU, cerr);
          }

          if (currsumLold > currsumL) {
            // cerr << "Error: currsumLold = " << currsumLold << " > " <<
            // currsumL
            std::cout << "Error: currsumLold = " << currsumLold << " > "
                      << currsumL << " = currsumL (n,m,k,l,j) = (" << n << ","
                      << m << "," << k << "," << l << "," << j << ")." << endl;
            // exit(1);
          }
          if (currsumUold < currsumU) {
            // cerr << "Error: currsumUold = " << currsumUold << " < " <<
            // currsumU
            std::cout << "Error: currsumUold = " << currsumUold << " < "
                      << currsumU << " = currsumU (n,m,k,l,j) = (" << n << ","
                      << m << "," << k << "," << l << "," << j << ")." << endl;
            // exit(1);
          }

          decision_on_mklj_made = (currsumL > u || currsumU < u);
        }

        mklj_found = (currsumL > u);
        if (mklj_found) {
          curr_mk[n][2] = l;
          curr_mk[n][3] = j;
        }
      }
    }

    if (counter == totalpts)  /// Gaussian approximation leads to currsum
                              /// summing to < 1.0, so we renormalise and sample
    {
      LogSumExp(currsumStore, runningMax);
      double100 sum = 0.0;
      int ind = 0;

      bool found = false;

      while (!found) {
        sum += currsumStore[ind];
        if (sum > u) {
          found = true;
        }
        if (ind == static_cast<int>(currsumStore.size() - 1)) {
          found = true;
        }
        ind++;
      }

      vector<int> returnvec;

      returnvec.push_back(mkljStore[0 + (4 * ind)]);
      returnvec.push_back(mkljStore[1 + (4 * ind)]);
      returnvec.push_back(mkljStore[2 + (4 * ind)]);
      returnvec.push_back(mkljStore[3 + (4 * ind)]);

      return returnvec;
    }
  }

  int coeffcount = 0;
  for (int i = 0; i <= n; ++i) coeffcount += (v_used[i] + 1);
  curr_mk[n].push_back(coeffcount);
  curr_mk[n].push_back(0);

  if (o.debug > 2) {
    cerr << "p_m,k,l,j: Returned (m,k,l,j) = (" << curr_mk[n][0] << ","
         << curr_mk[n][1] << "," << curr_mk[n][2] << "," << curr_mk[n][3]
         << ")\n";
    cerr << "n =\t\t\t";
    for (int i = 0; i <= n; ++i) cerr << i << "\t";
    cerr << "\nv_used =\t";
    for (int i = 0; i <= n; ++i) cerr << v_used[i] << "\t";
    cerr << "Coeffcount = " << coeffcount << endl;
  }

  return curr_mk[n];
}

vector<int> WrightFisher::DrawBridgePMFG1984(
    double100 x, double100 z, double100 s, double100 t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending in the
               /// interior of (0,1), but both time increments are below the
               /// threshold
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  vector<int> returnvec;
  vector<double100> currsumStore;
  vector<int> mkljStore;
  bool mklj_found = false,
       earlyStop =
           (abs(x - z) <= 0.6);  /// earlyStop gauges whether we can use currsum
  /// as is, or whether we should compute all
  /// probabilities and then sample

  /// Compute denominator
  double100 eC = 0.0, eCInc = 1.0, eCOldInc = 1.0;
  int Dflip = 1, Djm = 0, Djp = 0;
  int dmode = static_cast<int>(ceil(GriffithsParas(t).first)), d = dmode;

  while (max(eCOldInc, eCInc) > 0.0 || (!(eC > 0.0))) {
    eCOldInc = eCInc;
    double100 addon = 0.0;
    for (int f = 0; f != d; f++) {
      double100 para1 = static_cast<double100>(thetaP[0] + f),
                para2 = static_cast<double100>(thetaP[1] + d - f);
      boost::math::binomial_distribution<double100> BIN(d, x);
      boost::math::beta_distribution<double100> BETA(para1, para2);

      addon += pdf(BIN, f) * pdf(BETA, z);
    }
    eCInc = QmApprox(d, t, o) * addon;
    eC += eCInc;

    if (Dflip == -1 &&
        (dmode - Djm - 1 > 0))  // Mechanism to explore either side around mode
    {
      Djm++;
      d = dmode - Djm;
    } else {
      Djp++;
      d = dmode + Djp;
    }
    Dflip *= -1;
  }

  vector<int> modeGuess =
      mkljModeFinder(x, z, s, t, o);  /// Get a guess on mode over (m,k,l,j)
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2],
      jMode = modeGuess[3];

  boost::random::uniform_01<double100>
      U01;  /// Use this guess & eC to compute a suitable threshold for
  /// subsequent computations
  double100 currsum = 0.0, u = U01(gen),
            threshold = exp(mkljModeFinder_Evaluator(mMode, kMode, lMode, jMode,
                                                     x, z, s, t, o) -
                            log(eC)) *
                        1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false,
       mUpSwitch = false;  /// Compute m contributions
  double100 mContr_D = boost::math::lgamma(
                static_cast<double100>(thetaP[0] + thetaP[1] + m)),
            mContr_U = mContr_D, mContr, runningMax = -1.0e100;

  while (!mSwitch) {
    double100 qm = QmApprox(m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to very
                      /// small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    if (m != mMode) {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m - 1)));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D -=
            log(static_cast<double100>(thetaP[0] + thetaP[1] + (m + 1) - 1));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Compute k contributions
    double100 kContr_D =
                  (boost::math::lgamma(
                       static_cast<double100>(thetaP[0] + thetaP[1] + k)) -
                   boost::math::lgamma(
                       static_cast<double100>(thetaP[0] + thetaP[1] + m + k))),
              kContr_U = kContr_D, kContr;

    while (!kSwitch) {
      double100 qk = QmApprox(k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      if (!(qk > 0.0)) {
        break;
      }
      if (k != kMode) {
        if (kU > kD) {
          kContr_U +=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k - 1))) -
              log(static_cast<double100>(thetaP[0] + thetaP[1] + m + (k - 1)));
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D +=
              log(static_cast<double100>(thetaP[0] + thetaP[1] + m + (k + 1) -
                                         1)) -
              log(static_cast<double100>(thetaP[0] + thetaP[1] + (k + 1) - 1));
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      int lFlip = 1, lU = 0, lD = 0, newlMode = min(lMode, m),
          l = newlMode;  /// Redefine lMode in case m is too small!
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;
      boost::math::binomial_distribution<double100> BINL(
          m, x);  /// Compute l contributions
      double100 lContr_D = (log(pdf(BINL, l)) -
                            boost::math::lgamma(
                                static_cast<double100>(thetaP[0] + l)) -
                            boost::math::lgamma(
                                static_cast<double100>(thetaP[1] + m - l))),
                lContr_U = lContr_D, lContr;

      while (!lSwitch) {
        assert((l >= 0) && (l <= m));
        if (l != newlMode) {
          if (lU > lD) {
            lContr_U +=
                log(static_cast<double100>(m - (l - 1))) -
                log(static_cast<double100>((l - 1) + 1)) + log(x) -
                log(1.0 - x) +
                log(static_cast<double100>(thetaP[1] + m - (l - 1) - 1)) -
                log(static_cast<double100>(thetaP[0] + (l - 1)));
            lContr = lContr_U;
          } else {
            lContr_D += log(static_cast<double100>(l + 1)) -
                        log(static_cast<double100>(m - (l + 1) + 1)) +
                        log(1.0 - x) - log(x) +
                        log(static_cast<double100>(thetaP[0] + (l + 1) - 1)) -
                        log(static_cast<double100>(thetaP[1] + m - (l + 1)));
            lContr = lContr_D;
          }
        } else {
          lContr = lContr_U;
        }

        int jFlip = 1, jU = 0, jD = 0, newjMode = min(jMode, k),
            j = newjMode;  /// Redefine jMode in case k is too small!
        bool jSwitch = false, jDownSwitch = false,
             jUpSwitch = false;  /// Compute j contributions

        double100 jContr_D =
            LogBinomialCoefficientCalculator(k, j) +
            boost::math::lgamma(static_cast<double100>(thetaP[0] + l + j)) -
            boost::math::lgamma(static_cast<double100>(thetaP[0] + j)) +
            boost::math::lgamma(
                static_cast<double100>(thetaP[1] + m - l + k - j)) -
            boost::math::lgamma(static_cast<double100>(thetaP[1] + k - j)) +
            static_cast<double100>(thetaP[0] + j - 1) * log(z) +
            static_cast<double100>(thetaP[1] + k - j - 1) * log(1.0 - z);
        double100 jContr_U = jContr_D, jContr;

        while (!jSwitch) {
          if (j != newjMode) {
            if (jU > jD) {
              jContr_U +=
                  log(static_cast<double100>(k - (j - 1))) -
                  log(static_cast<double100>((j - 1) + 1)) + log(z) -
                  log(1.0 - z) +
                  log(static_cast<double100>(thetaP[0] + l + (j - 1))) -
                  log(static_cast<double100>(thetaP[0] + (j - 1))) +
                  log(static_cast<double100>(thetaP[1] + k - (j - 1) - 1)) -
                  log(static_cast<double100>(thetaP[1] + m - l + k - (j - 1) -
                                             1));
              jContr = jContr_U;
            } else {
              jContr_D +=
                  log(static_cast<double100>(j + 1)) -
                  log(static_cast<double100>(k - (j + 1) + 1)) + log(1.0 - z) -
                  log(z) +
                  log(static_cast<double100>(thetaP[0] + (j + 1) - 1)) -
                  log(static_cast<double100>(thetaP[0] + l + (j + 1) - 1)) +
                  log(static_cast<double100>(thetaP[1] + m - l + k - (j + 1))) -
                  log(static_cast<double100>(thetaP[1] + k - (j + 1)));
              jContr = jContr_D;
            }
          } else {
            jContr = jContr_U;
          }
          double100 currsum_inc =
              exp(mContr + kContr + lContr + jContr - log(eC));
          runningMax = max(
              currsum_inc,
              runningMax);  /// Running max needed for log-sum-exp trick later
          currsum += currsum_inc;
          currsumStore.push_back(currsum);

          mkljStore.resize(mkljStore.size() + 4, -1);
          mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
          mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
          mkljStore[2 + (4 * (currsumStore.size() - 1))] = l;
          mkljStore[3 + (4 * (currsumStore.size() - 1))] = j;

          if ((currsum > u) && earlyStop)  /// if earlyStop is allowed, we can
                                           /// stop once currsum exceeds u
          {
            returnvec.push_back(m);
            returnvec.push_back(k);
            returnvec.push_back(l);
            returnvec.push_back(j);

            mklj_found = true;
            goto End;
          }

          if (!(jDownSwitch))  /// Switching mechanism for j
          {
            if (sgn(j - newjMode) <= 0) {
              jDownSwitch =
                  ((currsum_inc < threshold) || (newjMode - jD - 1) < 0);
            }
          }

          if (!(jUpSwitch)) {
            if (sgn(j - newjMode) >= 0) {
              jUpSwitch =
                  ((currsum_inc < threshold) || (newjMode + jU + 1) > k);
            }
          }

          jSwitch = (jDownSwitch && jUpSwitch);

          if (!jSwitch) {
            if ((jFlip == 1 && (newjMode + jU + 1 <= k)) ||
                (jDownSwitch && !(jUpSwitch))) {
              jU++;
              j = newjMode + jU;
              jFlip *= (jDownSwitch ? 1 : -1);
            } else if ((jFlip == -1 && (newjMode - jD - 1 >= 0)) ||
                       (jUpSwitch && !(jDownSwitch))) {
              jD++;
              j = newjMode - jD;
              jFlip *= (jUpSwitch ? 1 : -1);
            }
          }
        }

        if (!(lDownSwitch))  /// Switching mechanism for l
        {
          if (sgn(l - newlMode) <= 0) {
            lDownSwitch = (((jU == 0) && (jD == 0)) || (newlMode - lD - 1) < 0);
          }
        }

        if (!(lUpSwitch)) {
          if (sgn(l - newlMode) >= 0) {
            lUpSwitch = (((jU == 0) && (jD == 0)) || (newlMode + lU + 1) > m);
          }
        }

        lSwitch = (lDownSwitch && lUpSwitch);

        if (!lSwitch) {
          if ((lFlip == 1 && (newlMode + lU + 1 <= m)) ||
              (lDownSwitch && !(lUpSwitch))) {
            lU++;
            l = newlMode + lU;
            lFlip *= (lDownSwitch ? 1 : -1);
          } else if ((lFlip == -1 && (newlMode - lD - 1 >= 0)) ||
                     (lUpSwitch && !(lDownSwitch))) {
            lD++;
            l = newlMode - lD;
            lFlip *= (lUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch = (((lU == 0) && (lD == 0)) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = ((lU == 0) && (lD == 0));
        }
      }

      kSwitch = (kDownSwitch && kUpSwitch);

      if (!kSwitch) {
        if (kFlip == 1) {
          kU++;
          k = kMode + kU;
          kFlip *= (kDownSwitch ? 1 : -1);
        } else if ((kFlip == -1) && (kMode - kD - 1 >= 0)) {
          kD++;
          k = kMode - kD;
          kFlip *= (kUpSwitch ? 1 : -1);
        }
      }
    }

    if (!(mDownSwitch))  /// Switching mechanism for m
    {
      if (sgn(m - mMode) <= 0) {
        mDownSwitch = (((kU == 0) && (kD == 0)) || (mMode - mD - 1 < 0));
      }
    }

    if (!(mUpSwitch)) {
      if (sgn(m - mMode) >= 0) {
        mUpSwitch = ((kU == 0) && (kD == 0));
      }
    }

    mSwitch = (mDownSwitch && mUpSwitch);

    if (!mSwitch) {
      if (mFlip == 1) {
        mU++;
        m = mMode + mU;
        mFlip *= (mDownSwitch ? 1 : -1);
      } else if ((mFlip == -1) && (mMode - mD - 1 >= 0)) {
        mD++;
        m = mMode - mD;
        mFlip *= (mUpSwitch ? 1 : -1);
      }
    }
  }

End:

  if (!mklj_found)  /// Either earlyStop disable or even with it we still did
                    /// not get currsum > u - guards against cases when
                    /// earlyStop is not a good enough indicator of currsum
                    /// summing to < 1
  {
    LogSumExp(currsumStore,
              runningMax);  /// log-sum-exp trick to normalise and transform
    /// vector of log probabilities
    double100 sum = 0.0;
    int index, ind = 0;

    vector<int> indexing(currsumStore.size(), 0);
    for (int i = 0; i != static_cast<int>(indexing.size()); i++) {
      indexing[i] = i;
    }
    sort(indexing.begin(), indexing.end(),  /// Sort probabilities
         [&](const int &a, const int &b) {
           return (currsumStore[a] > currsumStore[b]);
         });

    bool found = false;

    while (!found) {
      sum += currsumStore[indexing[ind]];
      if (sum > u) {
        index = indexing[ind];  /// Return correct index to sorted probabilities
        found = true;
      }
      if (ind ==
          static_cast<int>(currsumStore.size() - 1))  /// Ending condition
      {
        index = indexing[ind];
        found = true;
      }
      ind++;
    }

    returnvec.push_back(mkljStore[0 + (4 * index)]);
    returnvec.push_back(mkljStore[1 + (4 * index)]);
    returnvec.push_back(mkljStore[2 + (4 * index)]);
    returnvec.push_back(mkljStore[3 + (4 * index)]);
  }

  return returnvec;
}

double100 WrightFisher::mkModeFinder_Evaluator(
    int m, int k, double100 x, double100 z, double100 s, double100 t,
    const Options &o)  /// Evaluation function for finding mode over (m,k)
{
  assert((m >= 0) && (k >= 0) && (x >= 0.0) && (x <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (s > 0.0) && (s < t));
  double100 qm = QmApprox(m, s, o), qk = QmApprox(k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If they
                    /// are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }

  if (x != z) {
    return log(qm) + log(qk) + boost::math::lgamma(theta + m) +
           boost::math::lgamma(theta + k) - boost::math::lgamma(theta + m + k);
  } else if (!(x > 0.0)) {
    return log(qm) + log(qk) + boost::math::lgamma(thetaP[1] + m + k) +
           boost::math::lgamma(theta + m) + boost::math::lgamma(theta + k) -
           boost::math::lgamma(thetaP[1] + m) -
           boost::math::lgamma(thetaP[1] + k) -
           boost::math::lgamma(theta + m + k);
  } else {
    return log(qm) + log(qk) + boost::math::lgamma(thetaP[0] + m + k) +
           boost::math::lgamma(theta + m) + boost::math::lgamma(theta + k) -
           boost::math::lgamma(thetaP[0] + m) -
           boost::math::lgamma(thetaP[0] + k) -
           boost::math::lgamma(theta + m + k);
  }
}

double100 WrightFisher::mkjModeFinder_Evaluator(
    int m, int k, int j, double100 x, double100 z, double100 s, double100 t,
    const Options &o)  /// Evaluation function for finding mode over (m,k,j)
{
  assert((m >= 0) && (k >= 0) && (j >= 0) && (j <= k) && (x >= 0.0) &&
         (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > 0.0) && (s < t));
  boost::math::binomial_distribution<> BIN(k, z);
  double100 qm = QmApprox(m, s, o), qk = QmApprox(k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If they
                    /// are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }
  return log(qm) + log(qk) + log(pdf(BIN, j)) +
         boost::math::lgamma(static_cast<double100>(theta + m)) +
         boost::math::lgamma(static_cast<double100>(theta + k)) -
         boost::math::lgamma(static_cast<double100>(theta + m + k)) +
         (!(x > 0.0)
              ? boost::math::lgamma(
                    static_cast<double100>(thetaP[1] + m + k - j)) -
                    boost::math::lgamma(static_cast<double100>(thetaP[1] + m)) -
                    boost::math::lgamma(
                        static_cast<double100>(thetaP[1] + k - j))
              : boost::math::lgamma(static_cast<double100>(thetaP[0] + m + j)) -
                    boost::math::lgamma(static_cast<double100>(thetaP[0] + m)) -
                    boost::math::lgamma(static_cast<double100>(thetaP[0] + j)));
}

double100 WrightFisher::mkljModeFinder_Evaluator(
    int m, int k, int l, int j, double100 x, double100 z, double100 s,
    double100 t,
    const Options &o)  /// Evaluation function for finding mode over (m,k,l,j)
{
  assert((m >= 0) && (k >= 0) && (j >= 0) && (j <= k) && (l >= 0) && (l <= m) &&
         (x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > 0.0) &&
         (s < t));
  boost::math::binomial_distribution<> BIN(m, x), BINZ(k, z);
  double100 qm = QmApprox(m, s, o), qk = QmApprox(k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If they
                    /// are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }
  return log(qm) + log(qk) + log(pdf(BIN, l)) + log(pdf(BINZ, j)) +
         boost::math::lgamma(static_cast<double100>(thetaP[0] + l + j)) +
         boost::math::lgamma(
             static_cast<double100>(thetaP[1] + m - l + k - j)) +
         boost::math::lgamma(static_cast<double100>(theta + m)) -
         boost::math::lgamma(static_cast<double100>(theta + m + k)) -
         boost::math::lgamma(static_cast<double100>(thetaP[0] + l)) -
         boost::math::lgamma(static_cast<double100>(thetaP[1] + m - l));
}

double100 WrightFisher::mklModeFinder_Evaluator(
    int m, int k, int l, double100 x, double100 z, double100 s, double100 t,
    const Options &o)  /// Evaluation function for finding mode over (m,k,l)
{
  assert((m >= 0) && (k >= 0) && (l >= 0) && (l <= m) && (x > 0.0) &&
         (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (s > 0.0) && (s < t));
  double100 qm = QmApprox(m, s, o), qk = QmApprox(k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If they
                    /// are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }
  boost::math::binomial_distribution<> BIN(m, x);
  return log(qm) + log(qk) + log(pdf(BIN, l)) +
         boost::math::lgamma(static_cast<double100>(theta + m)) -
         boost::math::lgamma(theta + m + k) +
         (!(z > 0.0)
              ? boost::math::lgamma(static_cast<double100>(theta + m - l + k)) -
                    boost::math::lgamma(static_cast<double100>(theta + m - l))
              : boost::math::lgamma(static_cast<double100>(theta + l + k)) -
                    boost::math::lgamma(static_cast<double100>(theta + l)));
}

vector<int> WrightFisher::mkModeFinder(
    double100 x, double100 z, double100 s, double100 t,
    const Options &o)  /// Routine for finding mode over (m,k)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(s).first)),
      k = static_cast<int>(floor(GriffithsParas(t - s).first));

  int m_ud, k_ud;  /// Start at mode from qm and qk
  double100 currMode_eval = mkModeFinder_Evaluator(m, k, x, z, s, t, o);

  bool stop = false;

  /// Iteratively increment m and k depending on whether function increases or
  /// decreases with proposed move
  while (!stop) {
    if (currMode_eval < mkModeFinder_Evaluator(m + 1, k, x, z, s, t, o)) {
      m_ud = 1;
    } else if ((m - 1 >= 0) && currMode_eval < mkModeFinder_Evaluator(
                                                   m - 1, k, x, z, s, t, o)) {
      m_ud = -1;
    } else {
      m_ud = 0;
    }

    m += m_ud;
    currMode_eval = mkModeFinder_Evaluator(m, k, x, z, s, t, o);

    if (currMode_eval < mkModeFinder_Evaluator(m, k + 1, x, z, s, t, o)) {
      k_ud = 1;
    } else if ((k - 1 >= 0) && currMode_eval < mkModeFinder_Evaluator(
                                                   m, k - 1, x, z, s, t, o)) {
      k_ud = -1;
    } else {
      k_ud = 0;
    }

    k += k_ud;
    currMode_eval = mkModeFinder_Evaluator(m, k, x, z, s, t, o);

    if (!(m_ud > 0) && !(m_ud < 0) && !(k_ud > 0) &&
        !(k_ud < 0))  /// Keep iterating until both m & k told to not change
    {
      stop = true;
      returnvec.push_back(m);
      returnvec.push_back(k);
    }
  }

  return returnvec;
}

vector<int> WrightFisher::mkjModeFinder(
    double100 x, double100 z, double100 s, double100 t,
    const Options &o)  /// Routine for finding mode over (m,k,j)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(s).first)),
      k = static_cast<int>(floor(GriffithsParas(t - s).first));
  int j = static_cast<int>(floor(static_cast<double100>(k) * z));

  int m_ud, k_ud, j_ud;  /// Starting from the modes of qm, qk and Bin(k,z)
  double100 currMode_eval = mkjModeFinder_Evaluator(m, k, j, x, z, s, t, o);

  bool stop = false;

  /// Iteratively increment m and (k,j) depending on whether function increases
  /// or decreases with proposed move - (k,j) updated jointly to make sure 0 <=
  /// j <= k at all times
  while (!stop) {
    if (currMode_eval < mkjModeFinder_Evaluator(m + 1, k, j, x, z, s, t, o)) {
      m_ud = 1;
    } else if ((m - 1 >= 0) &&
               currMode_eval <
                   mkjModeFinder_Evaluator(m - 1, k, j, x, z, s, t, o)) {
      m_ud = -1;
    } else {
      m_ud = 0;
    }

    m += m_ud;
    currMode_eval = mkjModeFinder_Evaluator(m, k, j, x, z, s, t, o);

    if (currMode_eval < mkjModeFinder_Evaluator(m, k + 1, j, x, z, s, t, o)) {
      k_ud = 1;
      if (currMode_eval <
          mkjModeFinder_Evaluator(m, k + 1, j + 1, x, z, s, t, o)) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval <
                  mkjModeFinder_Evaluator(m, k + 1, j - 1, x, z, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j <= k - 1) &&
               currMode_eval <
                   mkjModeFinder_Evaluator(m, k - 1, j, x, z, s, t, o)) {
      k_ud = -1;
      if ((j + 1 <= k - 1) &&
          (currMode_eval <
           mkjModeFinder_Evaluator(m, k - 1, j + 1, x, z, s, t, o))) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval <
                  mkjModeFinder_Evaluator(m, k - 1, j - 1, x, z, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j - 1 >= 0) &&
               (currMode_eval <
                mkjModeFinder_Evaluator(m, k - 1, j - 1, x, z, s, t, o))) {
      k_ud = -1;
      j_ud = -1;
    } else {
      k_ud = 0;
      if (k == 0) {
        j_ud = 0;
      } else {
        if ((j + 1 <= k) && (currMode_eval < mkjModeFinder_Evaluator(
                                                 m, k, j + 1, x, z, s, t, o))) {
          j_ud = 1;
        } else if ((j - 1 >= 0) &&
                   (currMode_eval <
                    mkjModeFinder_Evaluator(m, k, j - 1, x, z, s, t, o))) {
          j_ud = -1;
        } else {
          j_ud = 0;
        }
      }
    }

    k += k_ud;
    j += j_ud;
    currMode_eval = mkjModeFinder_Evaluator(m, k, j, x, z, s, t, o);

    /// Iterate until told to not update (m,k,j) anymore
    if (!(m_ud > 0) && !(m_ud < 0) && !(k_ud > 0) && !(k_ud < 0) &&
        !(j_ud > 0) && !(j_ud < 0)) {
      stop = true;
      returnvec.push_back(m);
      returnvec.push_back(k);
      returnvec.push_back(j);
    }
  }

  return returnvec;
}

vector<int> WrightFisher::mkljModeFinder(
    double100 x, double100 z, double100 s, double100 t,
    const Options &o)  /// Routine for finding mode over (m,k,l,j)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(s).first)),
      k = static_cast<int>(floor(GriffithsParas(t - s).first));
  int l = static_cast<int>(floor(static_cast<double100>(m) * x)),
      j = static_cast<int>(floor(static_cast<double100>(k) * z));

  int m_ud, k_ud, l_ud,
      j_ud;  /// Initialise at modes from qm, qk, Bin(m,x), Bin(k,z)
  double100 currMode_eval = mkljModeFinder_Evaluator(m, k, l, j, x, z, s, t, o);

  bool stop = false;

  /// Iteratively increment (m,l) and (k,j) depending on whether function
  /// increases or decreases with proposed move - (m,l) & (k,j) updated jointly
  /// to make sure 0 <= l <= m & 0 <= j <= k at all times
  while (!stop) {
    if (currMode_eval <
        mkljModeFinder_Evaluator(m + 1, k, l, j, x, z, s, t, o)) {
      m_ud = 1;
      if (currMode_eval <
          mkljModeFinder_Evaluator(m + 1, k, l + 1, j, x, z, s, t, o)) {
        l_ud = 1;
      } else if ((l - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(m + 1, k, l - 1, j,
                                                           x, z, s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= 0) && (l <= m - 1) &&
               currMode_eval <
                   mkljModeFinder_Evaluator(m - 1, k, l, j, x, z, s, t, o)) {
      m_ud = -1;
      if ((l + 1 <= m - 1) &&
          (currMode_eval <
           mkljModeFinder_Evaluator(m - 1, k, l + 1, j, x, z, s, t, o))) {
        l_ud = 1;
      } else if ((l - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(m - 1, k, l - 1, j,
                                                           x, z, s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= 0) && (l - 1 >= 0) &&
               (currMode_eval <
                mkljModeFinder_Evaluator(m - 1, k, l - 1, j, x, z, s, t, o))) {
      m_ud = -1;
      l_ud = -1;
    } else {
      m_ud = 0;
      if (m == 0) {
        l_ud = 0;
      } else {
        if ((l + 1 <= m) &&
            (currMode_eval <
             mkljModeFinder_Evaluator(m, k, l + 1, j, x, z, s, t, o))) {
          l_ud = 1;
        } else if ((l - 1 >= 0) &&
                   (currMode_eval <
                    mkljModeFinder_Evaluator(m, k, l - 1, j, x, z, s, t, o))) {
          l_ud = -1;
        } else {
          l_ud = 0;
        }
      }
    }

    m += m_ud;
    l += l_ud;
    currMode_eval = mkljModeFinder_Evaluator(m, k, l, j, x, z, s, t, o);

    if (currMode_eval <
        mkljModeFinder_Evaluator(m, k + 1, l, j, x, z, s, t, o)) {
      k_ud = 1;
      if (currMode_eval <
          mkljModeFinder_Evaluator(m, k + 1, l, j + 1, x, z, s, t, o)) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(m, k + 1, l, j - 1,
                                                           x, z, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j <= k - 1) &&
               currMode_eval <
                   mkljModeFinder_Evaluator(m, k - 1, l, j, x, z, s, t, o)) {
      k_ud = -1;
      if ((j + 1 <= k - 1) &&
          currMode_eval <
              mkljModeFinder_Evaluator(m, k - 1, l, j + 1, x, z, s, t, o)) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(m, k - 1, l, j - 1,
                                                           x, z, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j - 1 >= 0) &&
               (currMode_eval <
                mkljModeFinder_Evaluator(m, k - 1, l, j - 1, x, z, s, t, o))) {
      k_ud = -1;
      j_ud = -1;
    } else {
      k_ud = 0;
      if (k == 0) {
        j_ud = 0;
      } else {
        if ((j + 1 <= k) &&
            (currMode_eval <
             mkljModeFinder_Evaluator(m, k, l, j + 1, x, z, s, t, o))) {
          j_ud = 1;
        } else if ((j - 1 >= 0) &&
                   (currMode_eval <
                    mkljModeFinder_Evaluator(m, k, l, j - 1, x, z, s, t, o))) {
          j_ud = -1;
        } else {
          j_ud = 0;
        }
      }
    }

    k += k_ud;
    j += j_ud;
    currMode_eval = mkljModeFinder_Evaluator(m, k, l, j, x, z, s, t, o);

    /// Iterate until (m,k,l,j) told not to change any more
    if (!(m_ud > 0) && !(m_ud < 0) && !(k_ud > 0) && !(k_ud < 0) &&
        !(j_ud > 0) && !(j_ud < 0)) {
      stop = true;
      returnvec.push_back(m);
      returnvec.push_back(k);
      returnvec.push_back(l);
      returnvec.push_back(j);
    }
  }

  return returnvec;
}

vector<int> WrightFisher::mklModeFinder(
    double100 x, double100 z, double100 s, double100 t,
    const Options &o)  /// Routine for finding mode over (m,k,l)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(s).first)),
      k = static_cast<int>(floor(GriffithsParas(t - s).first)),
      bumper = (thetaP.empty() ? 1 : 0);
  int l = max(static_cast<int>(floor(static_cast<double100>(m) * x)), bumper);

  int m_ud = -1, k_ud = -1, l_ud = -1,
      mkLower = thetaP.empty()
                    ? 1
                    : 0;  /// Starting from the modes of qm, qk and Bin(m,x)
  double100 currMode_eval = mklModeFinder_Evaluator(m, k, l, x, z, s, t, o);

  bool rerun = false;

  /// Iteratively increment k and (m,l) depending on whether function increases
  /// or decreases with proposed move - (m,l) updated jointly to make sure 0 <=
  /// l <= m at all times
  while (!((m_ud == 0) && (k_ud == 0) && (l_ud == 0))) {
  Rerun:

    if (currMode_eval < mklModeFinder_Evaluator(m, k + 1, l, x, z, s, t, o)) {
      k_ud = 1;
    } else if ((k - 1 >= mkLower) &&
               currMode_eval <
                   mklModeFinder_Evaluator(m, k - 1, l, x, z, s, t, o)) {
      k_ud = -1;
    } else {
      k_ud = 0;
    }

    k += k_ud;
    currMode_eval = mklModeFinder_Evaluator(m, k, l, x, z, s, t, o);
    int lLower = (thetaP.empty() || (!(thetaP[0] > 0.0) && !(z > 0.0)) ? 1 : 0),
        lUpper =
            (thetaP.empty() || (!(thetaP[1] > 0.0) && !(z < 1.0)) ? m - 1 : m);

    if (currMode_eval < mklModeFinder_Evaluator(m + 1, k, l, x, z, s, t, o)) {
      m_ud = 1;
      if (currMode_eval <
          mklModeFinder_Evaluator(m + 1, k, l + 1, x, z, s, t, o)) {
        l_ud = 1;
      } else if ((l - 1 >= lLower) &&
                 (currMode_eval <
                  mklModeFinder_Evaluator(m + 1, k, l - 1, x, z, s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= mkLower) && (l <= lUpper - 1) &&
               currMode_eval <
                   mklModeFinder_Evaluator(m - 1, k, l, x, z, s, t, o)) {
      m_ud = -1;
      if ((l + 1 <= lUpper - 1) &&
          (currMode_eval <
           mklModeFinder_Evaluator(m - 1, k, l + 1, x, z, s, t, o))) {
        l_ud = 1;
      } else if ((l - 1 >= lLower) &&
                 (currMode_eval <
                  mklModeFinder_Evaluator(m - 1, k, l - 1, x, z, s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= mkLower) && (l - 1 >= lLower) &&
               (currMode_eval <
                mklModeFinder_Evaluator(m - 1, k, l - 1, x, z, s, t, o))) {
      m_ud = -1;
      l_ud = -1;
    } else {
      m_ud = 0;
      if (m == mkLower) {
        l_ud = 0;
      } else {
        if ((l + 1 <= lUpper) &&
            (currMode_eval <
             mklModeFinder_Evaluator(m, k, l + 1, x, z, s, t, o))) {
          l_ud = 1;
        } else if ((l - 1 >= lLower) &&
                   (currMode_eval <
                    mklModeFinder_Evaluator(m, k, l - 1, x, z, s, t, o))) {
          l_ud = -1;
        } else {
          l_ud = 0;
        }
      }
    }

    m += m_ud;
    l += l_ud;
    currMode_eval = mklModeFinder_Evaluator(m, k, l, x, z, s, t, o);

    if (rerun == true && !(k_ud == 0 && m_ud == 0 && l_ud == 0)) {
      rerun = false;
    }

    if (k_ud == 0 && m_ud == 0 && l_ud == 0 && rerun == false) {
      rerun = true;
      goto Rerun;
    }

    /// Iterate until told to not update (m,k,l) anymore
  }

  if ((thetaP.empty() || (!(thetaP[0] > 0.0) && (l == 0))) ||
      (thetaP.empty() || (!(thetaP[1] > 0.0) && (l == m)))) {
    l = ((thetaP.empty() || (!(thetaP[0] > 0.0) && (l == 0)))
             ? 1
             : m - 1);  /// Sometimes l == 0/m but cannot be for cases when
                        /// thetaP.empty()
  }

  returnvec.push_back(m);
  returnvec.push_back(k);
  returnvec.push_back(l);

  return returnvec;
}

double100 WrightFisher::CustomGammaRatio(double100 a, double100 b) {
  double100 answer;
  // If arguments are large, just use Stirling approximation for gamma function:
  // Gamma(x+1) = sqrt(2*pi*x)*(x/e)^x
  if (min(a, b) > 100.0) {
    answer = 0.5 * (log(a - 1.0) - log(b - 1.0)) + (a - 1.0) * log(a - 1.0) -
             (b - 1.0) * log(b - 1.0) + (b - a);
  } else if (a != b) {  // Otherwise compute sum of corresponding logs
    // Could use boost::math::tgamma_ratio here, but below implementation is
    // faster
    answer = 0.0;
    double100 high_val = max(a, b), low_val = min(a, b);
    double100 start_val = (low_val);
    int counter = 0;
    while (counter < static_cast<int>(floor(high_val) - floor(low_val))) {
      answer += log(start_val);
      start_val += 1.0;
      counter++;
    }
    if (high_val == b) {
      answer = -answer;
    }
  } else {
    answer = 0.0;
  }
  // NB: Returns log of gamma ratio!
  return answer;
}

pair<double100, int> WrightFisher::DrawBridgepoint(
    double100 x, double100 z, double100 t1, double100 t2, double100 s,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Routine to decide which bridge sampler to invoke for sampling
               /// at time s from neutral bridge diffusion started at x at time
               /// t1, ending at z in time t2, conditioned on non-absorption on
               /// (t1,t2)
{
  assert((x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > t1) &&
         (s < t2));
  double100 y, para1, para2;
  int m, k, j, l, coeffcount = -1;
  vector<int> mklj;

  if ((s - t1 <= o.bridgethreshold || t2 - s <= o.bridgethreshold) ||
      (((theta - 1.0) / (exp(0.5 * theta * (s - t1)))) +
           ((theta - 1.0) / (exp(0.5 * theta * (t2 - s)))) >
       260.0))  /// Use diffusion approximations
  {
    // Draw M, K from corresponding ancestral process
    m = DrawAncestralProcessG1984(s - t1, gen);
    k = DrawAncestralProcessG1984(t2 - s, gen);

    // Setup theta vector appropriately
    if (thetaP.empty()) {
      ThetaResetter();
    }
    double100 theta1 = thetaP.front(), theta2 = thetaP.back();

    boost::math::binomial_distribution<double100> BIN_m(m, x);

    // Compute probabilities for L, J using custom funciton for gamma ratios,
    // and utilising log-sum-exp trick to normalise probabilities
    vector<vector<double100>> probs;
    vector<double100> factors_l_i;
    for (int l_i = 0; l_i <= m; l_i++) {
      vector<double100> log_probs_l_i(k + 1);
      double100 bin_pmf = pdf(BIN_m, l_i);
      double100 maxProb =
          static_cast<double100>(-std::numeric_limits<double>::max());
      for (int j_i = 0; j_i <= k; j_i++) {
        boost::math::beta_distribution<double100> BETA_z(theta1 + j_i,
                                                         theta2 + k - j_i);
        double100 add_on =
            log(bin_pmf) + log(pdf(BETA_z, z)) -
            log(boost::math::factorial<double100>(
                static_cast<double100>(j_i))) +
            CustomGammaRatio(static_cast<double100>(k + 1),
                             static_cast<double100>(k + 1 - j_i)) +
            CustomGammaRatio(theta1 + static_cast<double100>(l_i + j_i),
                             theta1 + static_cast<double100>(l_i)) +
            CustomGammaRatio(theta2 + static_cast<double100>(m - l_i + k - j_i),
                             theta2 + static_cast<double100>(m - l_i)) +
            CustomGammaRatio(theta1 + theta2 + static_cast<double100>(m),
                             theta1 + theta2 + static_cast<double100>(m + k));
        maxProb = max(maxProb, add_on);
        log_probs_l_i[j_i] = add_on;
      }
      double100 factor_li = LogSumExp(log_probs_l_i, maxProb);
      factors_l_i.push_back(factor_li);
      probs.push_back(log_probs_l_i);
    }

    boost::random::discrete_distribution<> DISC_l(factors_l_i);
    l = DISC_l(gen);
    boost::random::discrete_distribution<> DISC_j(probs[l]);
    j = DISC_j(gen);

    para1 = static_cast<double100>(theta1 + l + j),
    para2 = static_cast<double100>(theta2 + m - l + k - j);

    boost::random::gamma_distribution<> GAMMA1(para1, 1.0), GAMMA2(para2, 1.0);

    y = -1.0;
    while (!(0.0 < y && y < 1.0)) {
      double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }
  } else  /// Else use bridge simulator
  {
    if (!(x > 0.0))  /// x = 0
    {
      if (!(z > 0.0))  /// x = 0 & z = 0
      {
        if (s - t1 <= o.g1984threshold &&
            t2 - s <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFSameMutationApprox(x, s - t1, t2 - t1, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj =
              DrawBridgePMFSameMutationOneQApprox(x, s - t1, t2 - t1, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFSameMutation(x, s - t1, t2 - t1, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = 0;
        j = 0;
      } else if (!(z < 1.0))  /// x = 0 & z = 1
      {
        if (s - t1 <= o.g1984threshold &&
            t2 - s <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj =
              DrawBridgePMFDifferentMutationApprox(s - t1, t2 - t1, x, o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFDifferentMutationOneQApprox(s - t1, t2 - t1, x, o,
                                                          gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFDifferentMutation(s - t1, t2 - t1, x, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = 0;
        j = k;
      } else  /// x = 0 & z in (0,1)
      {
        if (s - t1 <= o.g1984threshold &&
            t2 - s <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFInteriorMutationApprox(x, z, s - t1, t2 - t1, o,
                                                     gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(x, z, s - t1, t2 - t1,
                                                         o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFInteriorMutation(x, z, s - t1, t2 - t1, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = mklj[2];
        j = mklj[3];
      }
    } else if (!(x < 1.0))  /// x = 1
    {
      if (!(z > 0.0))  /// x = 1 & z = 0
      {
        if (s - t1 <= o.g1984threshold &&
            t2 - s <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj =
              DrawBridgePMFDifferentMutationApprox(s - t1, t2 - t1, x, o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFDifferentMutationOneQApprox(s - t1, t2 - t1, x, o,
                                                          gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFDifferentMutation(s - t1, t2 - t1, x, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = m;
        j = 0;
      } else if (!(z < 1.0))  /// x = 1 & z = 1
      {
        if (s - t1 <= o.g1984threshold &&
            t2 - s <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFSameMutationApprox(x, s - t1, t2 - t1, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj =
              DrawBridgePMFSameMutationOneQApprox(x, s - t1, t2 - t1, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFSameMutation(x, s - t1, t2 - t1, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = m;
        j = k;
      } else  /// x = 1 & z in (0,1)
      {
        if (s - t1 <= o.g1984threshold &&
            t2 - s <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFInteriorMutationApprox(x, z, s - t1, t2 - t1, o,
                                                     gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(x, z, s - t1, t2 - t1,
                                                         o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFInteriorMutation(x, z, s - t1, t2 - t1, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = mklj[2];
        j = mklj[3];
      }
    } else  /// x in (0,1)
    {
      if (!(z > 0.0))  /// x in (0,1) & z = 0
      {
        double100 newt2 = t2 - t1, news = t2 - s;

        if (news <= o.g1984threshold &&
            newt2 - news <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFInteriorMutationApprox(
              z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time bridge
        }            /// One time increment both below threshold
        else if ((news <= o.g1984threshold &&
                  newt2 - news > o.g1984threshold) ||
                 (news > o.g1984threshold &&
                  newt2 - news <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(
              z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time bridge
        } else       /// Time increments are large enough for alternating series
                     /// method
        {
          mklj = DrawBridgePMFInteriorMutation(
              z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time bridge
        }

        m = mklj[0];
        k = mklj[1];
        l = mklj[2];
        j = mklj[3];
      } else if (!(z < 1.0))  /// x in (0,1) & z = 1
      {
        double100 newt2 = t2 - t1, news = t2 - s;

        if (news <= o.g1984threshold &&
            newt2 - news <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFInteriorMutationApprox(
              z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time bridge
        }            /// One time increment both below threshold
        else if ((news <= o.g1984threshold &&
                  newt2 - news > o.g1984threshold) ||
                 (news > o.g1984threshold &&
                  newt2 - news <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(
              z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time bridge
        } else       /// Time increments are large enough for alternating series
                     /// method
        {
          mklj = DrawBridgePMFInteriorMutation(
              z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time bridge
        }

        m = mklj[0];
        k = mklj[1];
        l = mklj[2];
        j = mklj[3];
      } else  /// x in (0,1) & z in (0,1)
      {
        if (s - t1 <= o.g1984threshold &&
            t2 - s < o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFG1984(x, z, s - t1, t2 - t1, o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFOneQApprox(x, z, s - t1, t2 - t1, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMF(x, z, s - t1, t2 - t1, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = mklj[2];
        j = mklj[3];
      }
    }

    para1 = static_cast<double100>(thetaP[0] + l + j),
    para2 = static_cast<double100>(thetaP[1] + m - l + k - j);

    boost::random::gamma_distribution<> GAMMA1(para1, 1.0), GAMMA2(para2, 1.0);

    y = -1.0;
    while (!(0.0 < y && y < 1.0)) {
      double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }
  }

  return make_pair(y, coeffcount);
}

pair<double100, int> WrightFisher::DrawUnconditionedBridge(
    double100 x, double100 z, double100 t1, double100 t2, double100 s,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Routine to decide which bridge sampler to invoke for sampling
               /// at time s from neutral bridge diffusion started at x at time
               /// t1, ending at z in time t2, with potential absorption at any
               /// time in between [t1,t2]
{
  assert((x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > t1) &&
         (s < t2));
  double100 y, para1, para2;
  int m = -1, k = -1, j = -1, l = -1, coeffcount = -1;
  vector<int> mklj;

  if ((s - t1 <= o.bridgethreshold || t2 - s <= o.bridgethreshold) ||
      (((theta - 1.0) / (exp(0.5 * theta * (s - t1)))) +
           ((theta - 1.0) / (exp(0.5 * theta * (t2 - s)))) >
       260.0))  /// Use diffusion approximations
  {
    /// Last condition checks that corresponding m and k terms are not too
    /// large
    /// to create computational bottleneck
    double100 y1 =
        DrawEndpoint(x, t1, s, o, gen)
            .first;  /// Diffusion approximation for when times are too small
    /// and bridge takes too long to compute

    y = abs((y1 - x) + ((t2 - s) / (t2 - t1)) * x +
            ((s - t1) / (t2 - t1)) * z);  /// Ensures y remains positive

    if (y > 1.0)  /// Ensure y remains <= 1.0
    {
      y = 1.0 - abs(1.0 - y);
    }
  } else {
    if (thetaP.empty()) {
      if (!(x > 0.0)) {
        if (!(z > 0.0)) {
          return make_pair(0.0, coeffcount);
        } else {
          cerr << "Right endpoint for bridge is impossible for the given "
                  "configuration! Returned value corresponds to the left "
                  "endpoint!";
          return make_pair(0.0, coeffcount);
        }
      } else if (!(x < 1.0)) {
        if (!(z < 1.0)) {
          return make_pair(1.0, coeffcount);
        } else {
          cerr << "Right endpoint for bridge is impossible for the given "
                  "configuration! Returned value corresponds to the left "
                  "endpoint!";
          return make_pair(1.0, coeffcount);
        }
      } else {
        if (!(z > 0.0)) {
          vector<int> mklj;

          if (s - t1 <= o.g1984threshold && t2 - s <= o.g1984threshold) {
            mklj =
                DrawBridgePMFUnconditionalApprox(x, z, s - t1, t2 - t1, o, gen);
          } else if ((s - t1 <= o.g1984threshold &&
                      t2 - s > o.g1984threshold) ||
                     (s - t1 > o.g1984threshold &&
                      t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFUnconditionalOneQApprox(x, z, s - t1, t2 - t1,
                                                        o, gen);
          } else {
            mklj = DrawBridgePMFUnconditional(x, z, s - t1, t2 - t1, o, gen);
          }

          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = 0;

          if (l == 0) {
            return make_pair(0.0, coeffcount);
          }
        } else if (!(z < 1.0)) {
          vector<int> mklj;

          if (s - t1 <= o.g1984threshold && t2 - s <= o.g1984threshold) {
            mklj =
                DrawBridgePMFUnconditionalApprox(x, z, s - t1, t2 - t1, o, gen);
          } else if ((s - t1 <= o.g1984threshold &&
                      t2 - s > o.g1984threshold) ||
                     (s - t1 > o.g1984threshold &&
                      t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFUnconditionalOneQApprox(x, z, s - t1, t2 - t1,
                                                        o, gen);
          } else {
            mklj = DrawBridgePMFUnconditional(x, z, s - t1, t2 - t1, o, gen);
          }

          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = k;

          if (l == m) {
            return make_pair(1.0, coeffcount);
          }
        } else  /// Otherwise the conditions imposed imply the diffusion cannot
                /// be absorbed over [t1,t2], so we can use Drawbridgepoint
        {
          ThetaResetter();
          return DrawBridgepoint(x, z, t1, t2, s, o, gen);
        }
      }
    } else if (!(thetaP[0] > 0.0) || !(thetaP[1] > 0.0)) {
      if (!(thetaP[0] > 0.0)) {
        if (!(x > 0.0)) {
          if (!(z > 0.0)) {
            return make_pair(0.0, coeffcount);
          } else {
            cerr << "Right endpoint for bridge is impossible for the given "
                    "configuration! Returned value corresponds to the left "
                    "endpoint!";
            return make_pair(0.0, coeffcount);
          }
        } else {
          if (!(z > 0.0)) {
            vector<int> mklj;

            if (s - t1 <= o.g1984threshold && t2 - s <= o.g1984threshold) {
              mklj = DrawBridgePMFUnconditionalApprox(x, z, s - t1, t2 - t1, o,
                                                      gen);
            } else if ((s - t1 <= o.g1984threshold &&
                        t2 - s > o.g1984threshold) ||
                       (s - t1 > o.g1984threshold &&
                        t2 - s <= o.g1984threshold)) {
              mklj = DrawBridgePMFUnconditionalOneQApprox(x, z, s - t1, t2 - t1,
                                                          o, gen);
            } else {
              mklj = DrawBridgePMFUnconditional(x, z, s - t1, t2 - t1, o, gen);
            }

            m = mklj[0];
            k = mklj[1];
            l = mklj[2];
            j = 0;

            if (l == 0) {
              return make_pair(0.0, coeffcount);
            }
          } else {
            ThetaResetter();
            return DrawBridgepoint(x, z, t1, t2, s, o, gen);
          }
        }
      } else {
        if (!(x < 1.0)) {
          if (!(z < 1.0)) {
            return make_pair(1.0, coeffcount);
          } else {
            cerr << "Right endpoint for bridge is impossible for the given "
                    "configuration! Returned value corresponds to the left "
                    "endpoint!";
            return make_pair(1.0, coeffcount);
          }
        } else {
          if (!(z < 1.0)) {
            vector<int> mklj;

            if (s - t1 <= o.g1984threshold && t2 - s <= o.g1984threshold) {
              mklj = DrawBridgePMFUnconditionalApprox(x, z, s - t1, t2 - t1, o,
                                                      gen);
            } else if ((s - t1 <= o.g1984threshold &&
                        t2 - s > o.g1984threshold) ||
                       (s - t1 > o.g1984threshold &&
                        t2 - s <= o.g1984threshold)) {
              mklj = DrawBridgePMFUnconditionalOneQApprox(x, z, s - t1, t2 - t1,
                                                          o, gen);
            } else {
              mklj = DrawBridgePMFUnconditional(x, z, s - t1, t2 - t1, o, gen);
            }

            m = mklj[0];
            k = mklj[1];
            l = mklj[2];
            j = 0;

            if (l == m) {
              return make_pair(1.0, coeffcount);
            }
          } else {
            ThetaResetter();
            return DrawBridgepoint(x, z, t1, t2, s, o, gen);
          }
        }
      }
    } else  /// No absorption can happen, so we are in the same case as in
            /// DrawBridgepoint
    {
      ThetaResetter();
      return DrawBridgepoint(x, z, t1, t2, s, o, gen);
    }
  }

  para1 = (thetaP.empty() ? static_cast<double100>(l + j)
                          : (!(thetaP[0] > 0.0)
                                 ? static_cast<double100>(l + j)
                                 : static_cast<double100>(thetaP[0] + l + j)));
  para2 = (thetaP.empty()
               ? static_cast<double100>(m - l + k - j)
               : (!(thetaP[1] > 0.0)
                      ? static_cast<double100>(m - l + k - j)
                      : static_cast<double100>(thetaP[1] + m - l + k - j)));

  boost::random::gamma_distribution<> GAMMA1(para1, 1.0), GAMMA2(para2, 1.0);

  y = -1.0;
  while (!(0.0 < y && y < 1.0)) {
    double100 G1 = GAMMA1(gen), G2 = GAMMA2(gen);
    y = G1 / (G1 + G2);
  }

  return make_pair(y, coeffcount);
}

/// BRIDGE SIMULATION - NON-NEUTRAL PATHS

vector<vector<double100>> WrightFisher::NonNeutralDrawBridge(
    double100 x, double100 t1, double100 t2, double100 z, bool Absorption,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Draws of paths from non-neutral WF diffusion bridge started
               /// from x at time t1 and ending at z at time t2
{
  bool accept = false;
  vector<double100> paras{phiMin, phiMax, phiMax - phiMin};
  double100 kapmean = paras[2] * (t2 - t1);  /// Rate of Poisson point process
  boost::random::poisson_distribution<int> kap(static_cast<double>(kapmean));

  boost::random::uniform_01<double100> unift,
      unifm;  /// Set up uniform generators for points over [t1,t2] *
  /// [0,phiMax-phiMin] * [0,1]
  vector<vector<double100>> ptr;
  int rcount = 0;
  vector<double100> kappastore;

  while (!accept)  /// Until all skeleton points get accepted, keep going
  {
    int kappa = kap(gen);  /// Generate kappa ~ Pois
    double100 small_offset = 1.0e-14;
    vector<double100> path, times(kappa), marks(kappa), rejcount;
    auto gent = [&t1, &t2, &unift, &small_offset,
                 &gen]()  /// time stamps ~ Unif [t1,t2]
    { return (t1 + small_offset + ((t2 - t1) * unift(gen))); };
    auto genm = [&paras, &unifm, &gen]()  /// marks ~ Unif [0,phiMax-phiMin]
    { return (paras[2] * unifm(gen)); };
    std::generate(begin(times), end(times), gent);
    std::generate(begin(marks), end(marks), genm);
    sortVectorsAscending(times, times,
                         marks);  /// Sort vectors according to timestamps

    if (kappa == 0)  /// No skeleton points -> accept
    {
      kappastore.push_back(1.0 * kappa);
      ptr.push_back(path);
      ptr.push_back(times);
      ptr.push_back(marks);
      ptr.push_back(kappastore);
      accept = true;
    } else  /// kappa > 0 - generate skeleton points and check them
    {
      for (vector<double100>::iterator itt = times.begin(), itm = marks.begin();
           itt != times.end(); itt++, itm++) {
        if (itt == times.begin()) {
          if (Absorption) {
            path.push_back(
                DrawUnconditionedBridge(x, z, t1, t2, *itt, o, gen).first);
          } else {
            path.push_back(
                DrawBridgepoint(x, z, t1, t2, *itt, o, gen)
                    .first);  /// Generate skeleton points sequentially
          }

          if (Phitilde(path.back()) - paras[0] >
              *itm)  /// Test generated point is OK, otherwise can stop and
                     /// generate a new Poisson point process
          {
            rcount++;
            break;
          }

          if (kappa == 1)  /// We only needed to generate one skeleton point,
                           /// which we accepted, so we can exit
          {
            kappastore.push_back(1.0 * kappa);
            ptr.push_back(path);
            ptr.push_back(times);
            ptr.push_back(marks);
            ptr.push_back(kappastore);
            accept = true;
          }
        } else if (*itt !=
                   times.back())  /// There are more than 2 skeleton points, and
                                  /// we are not at the last one yet
        {
          if (Absorption) {
            path.push_back(DrawUnconditionedBridge(path.back(), z, *(itt - 1),
                                                   t2, *itt, o, gen)
                               .first);
          } else {
            path.push_back(
                DrawBridgepoint(path.back(), z, *(itt - 1), t2, *itt, o, gen)
                    .first);  /// Generate skeleton points sequentially
          }

          if (Phitilde(path.back()) - paras[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop and
                     /// generate a new Poisson point process
          {
            rcount++;
            break;
          }
        } else  /// We are at the last skeleton point
        {
          if (Absorption) {
            path.push_back(
                DrawUnconditionedBridge(path.back(), z, *(itt - 1), t2, *itt, o,
                                        gen)
                    .first);  /// Generate skeleton point sequentially
          } else {
            path.push_back(
                DrawBridgepoint(path.back(), z, *(itt - 1), t2, *itt, o, gen)
                    .first);  /// Generate skeleton point sequentially
          }

          if (Phitilde(path.back()) - paras[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop and
                     /// generate a new Poisson point process
          {
            rcount++;
            break;
          }
          kappastore.push_back(1.0 * kappa);
          ptr.push_back(path);
          ptr.push_back(times);
          ptr.push_back(marks);
          ptr.push_back(kappastore);
          accept = true;
        }
      }
    }
  }

  return ptr;
}

pair<double100, int> WrightFisher::NonNeutralDrawBridgepoint(
    double100 x, double100 t1, double100 t2, double100 z, double100 testT,
    bool Absorption, const Options &o,
    boost::random::mt19937
        &gen)  /// Invoke NonNeutralDrawBridge to generate a whole bridge path,
               /// conditionally on the generated path, generate a neutral draw
               /// at time testT
{
  vector<double100> bridgeSection, bridgeTimes;
  bridgeSection.push_back(x);
  bridgeTimes.push_back(t1);

  vector<vector<double100>> skeleton = NonNeutralDrawBridge(
      x, t1, t2, z, Absorption, o, gen);  /// Create skeleton points for bridge
  bridgeSection.insert(bridgeSection.end(), skeleton[0].begin(),
                       skeleton[0].end());
  bridgeTimes.insert(bridgeTimes.end(), skeleton[1].begin(), skeleton[1].end());

  bridgeSection.push_back(z);
  bridgeTimes.push_back(t2);

  vector<double100>::iterator timeIt = bridgeTimes.begin(),
                              xIt = bridgeSection.begin();
  while (
      *timeIt <
      testT)  /// Cycle through skeleton points to find appropriate end points
  {
    timeIt++;
    xIt++;
  }
  timeIt--;
  xIt--;
  if (Absorption) {
    return DrawUnconditionedBridge(*xIt, *(xIt + 1), *timeIt, *(timeIt + 1),
                                   testT, o, gen);
  } else {
    return DrawBridgepoint(*xIt, *(xIt + 1), *timeIt, *(timeIt + 1), testT, o,
                           gen);  /// Return neutral draw using corresponding
                                  /// endpoints and time increments
  }
}

/// SIMULATION RUNNER FUNCTIONS

void WrightFisher::DiffusionRunner(
    int nSim, double100 x, double100 startT, double100 endT, bool Absorption,
    string &Filename, double100 diffusion_threshold,
    double100 bridge_threshold)  /// Function to generate specified number of
                                 /// diffusion draws
{
  std::cout << "You've asked to generate " << nSim
            << " draws from the law of a Wright--Fisher diffusion with the "
               "following properties:"
            << std::endl;
  std::cout << "Start point: " << x << std::endl;
  std::cout << "Start time: " << startT << std::endl;
  std::cout << "Sampling time: " << endT << std::endl;
  if (Absorption) {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary, but the provided mutation rates "
                   "are strictly positive, so the diffusion cannot be absorbed"
                << std::endl;
    }
  } else {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion cannot be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout
          << "You have further specified that the diffusion cannot be "
             "absorbed at the boundary, but the provided mutation rates "
             "are strictly positive and thus already ensure this, so the "
             "resulting draws are coming from the *unconditioned* diffusion!"
          << std::endl;
    }
  }
  std::cout << "You've further specified the time threshold for Gaussian "
               "approximations at "
            << diffusion_threshold << std::endl;
  std::cout << "Output will be printed to file in " << Filename << std::endl;
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);

  boost::random::mt19937 gen;
  int nosamples = 1, loader = floor(0.1 * nSim), loader_count = 1;
  while (nosamples < nSim + 1) {
    if (non_neutral) {
      saveFile
          << NonNeutralDrawEndpoint(x, startT, endT, Absorption, o, gen).first
          << "\n";
    } else {
      if (Absorption) {
        saveFile << DrawUnconditionedDiffusion(x, endT - startT, o, gen).first
                 << "\n";
      } else {
        saveFile << DrawEndpoint(x, startT, endT, o, gen).first << "\n";
      }
    }

    nosamples++;
    if (nosamples % (loader * loader_count) == 0) {
      std::cout << "Simulated " << nosamples << " samples." << endl;
      loader_count++;
    }
  }

  std::cout << "Diffusion simulation complete." << endl;
}

void WrightFisher::BridgeDiffusionRunner(
    int nSim, double100 x, double100 z, double100 startT, double100 endT,
    double100 sampleT, bool Absorption, string &Filename,
    double100 diffusion_threshold,
    double100 bridge_threshold)  /// Function to generate specified number of
                                 /// bridge draws
{
  std::cout
      << "You've asked to generate " << nSim
      << " draws from the law of a Wright--Fisher diffusion bridge with the "
         "following properties:"
      << std::endl;
  std::cout << "Start point: " << x << std::endl;
  std::cout << "Start time: " << startT << std::endl;
  std::cout << "End point: " << z << std::endl;
  std::cout << "End time: " << endT << std::endl;
  std::cout << "Sampling time: " << sampleT << std::endl;
  if (Absorption) {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary, but the provided mutation rates "
                   "are strictly positive, so the diffusion cannot be absorbed"
                << std::endl;
    }
  } else {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion cannot be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout
          << "You have further specified that the diffusion cannot be "
             "absorbed at the boundary, but the provided mutation rates "
             "are strictly positive and thus already ensure this, so the "
             "resulting draws are coming from the *unconditioned* diffusion!"
          << std::endl;
    }
  }
  std::cout << "You've further specified the time threshold for Gaussian "
               "approximations at "
            << diffusion_threshold
            << ", whilst the bridge approximations threshold was set to "
            << bridge_threshold << std::endl;
  std::cout << "Output will be printed to file in " << Filename << std::endl;
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);

  boost::random::mt19937 gen;
  int nosamples = 1, loader = floor(0.1 * nSim), loader_count = 1;
  while (nosamples < nSim + 1) {
    if (non_neutral) {
      saveFile << NonNeutralDrawBridgepoint(x, startT, endT, z, sampleT,
                                            Absorption, o, gen)
                      .first
               << "\n";
    } else {
      if (Absorption) {
        saveFile << DrawUnconditionedBridge(x, z, startT, endT, sampleT, o, gen)
                        .first
                 << "\n";
      } else {
        ThetaResetter();
        saveFile << DrawBridgepoint(x, z, startT, endT, sampleT, o, gen).first
                 << "\n";
      }
    }

    nosamples++;
    if (nosamples % (loader * loader_count) == 0) {
      std::cout << "Simulated " << nosamples << " samples." << endl;
      loader_count++;
    }
  }

  std::cout << "Bridge simulation complete." << endl;
}

void WrightFisher::DiffusionDensityCalculator(
    int meshSize, double100 x, double100 startT, double100 endT,
    bool Absorption, string &Filename, double100 diffusion_threshold,
    double100 bridge_threshold)  /// Function to compute truncated diffusion
                                 /// transition density
{
  std::cout << "You've asked to compute the transition density of a "
               "Wright--Fisher diffusion with the "
               "following properties:"
            << std::endl;
  std::cout << "Start point: " << x << std::endl;
  std::cout << "Start time: " << startT << std::endl;
  std::cout << "Sampling time: " << endT << std::endl;
  if (Absorption) {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary, but the provided mutation rates "
                   "are strictly positive, so the diffusion cannot be absorbed"
                << std::endl;
    }
  } else {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion cannot be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout
          << "You have further specified that the diffusion cannot be "
             "absorbed at the boundary, but the provided mutation rates "
             "are strictly positive and thus already ensure this, so the "
             "resulting draws are coming from the *unconditioned* diffusion!"
          << std::endl;
    }
  }
  std::cout << "You've further specified the time threshold for Gaussian "
               "approximations at "
            << diffusion_threshold << std::endl;
  std::cout << "The pointwise computation will be performed over a mesh "
               "consisting of "
            << meshSize << " equally spaced intervals on [0,1]" << std::endl;
  std::cout << "Output will be printed to file in " << Filename << std::endl;
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);

  int counter = 0;
  double100 timeInc = endT - startT,
            yinc = 1.0 / static_cast<double100>(meshSize), y, ycount = 0.1;
  if (non_neutral) {
    cerr << "Truncated density cannot be computed for non-neutral case due "
            "to presence of intractable quantities!";
  } else {
    while (counter <= meshSize) {
      if (counter == 0) {
        y = 0.0;
      } else if (counter == meshSize) {
        y = 1.0;
      } else {
        y += yinc;
      }
      if (Absorption) {
        if ((x > 0.0) && (x < 1.0)) {
          saveFile << y << " "
                   << UnconditionedDiffusionDensity(x, y, timeInc, o) << "\n";
        } else {
          if (!(x > y) && !(x < y)) {
            saveFile << y << " " << 1.0 << "\n";
          } else {
            saveFile << y << " " << 0.0 << "\n";
          }
        }
      } else {
        saveFile << y << " " << DiffusionDensityApproximation(x, y, timeInc, o)
                 << "\n";
      }

      if (y >= ycount) {
        std::cout << "Calculated density up to y = " << ycount << endl;
        ycount += 0.1;
      }
      counter++;
    }

    std::cout << "Density calculation complete." << endl;

    saveFile.close();
  }
}

void WrightFisher::BridgeDiffusionDensityCalculator(
    int meshSize, double100 x, double100 z, double100 startT, double100 endT,
    double100 sampleT, bool Absorption, string &Filename,
    double100 diffusion_threshold,
    double100 bridge_threshold)  /// Function to compute truncated diffusion
                                 /// bridge transition density
{
  std::cout << "You've asked to compute the transition density of a "
               "Wright--Fisher diffusion bridge with the "
               "following properties:"
            << std::endl;
  std::cout << "Start point: " << x << std::endl;
  std::cout << "Start time: " << startT << std::endl;
  std::cout << "End point: " << z << std::endl;
  std::cout << "End time: " << endT << std::endl;
  std::cout << "Sampling time: " << sampleT << std::endl;
  if (Absorption) {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout << "You have further specified that the diffusion can be "
                   "absorbed at the boundary, but the provided mutation rates "
                   "are strictly positive, so the diffusion cannot be absorbed"
                << std::endl;
    }
  } else {
    if (thetaP.empty() || ((thetaP.front() == 0.0 && thetaP.back() != 0.0) ||
                           (thetaP.front() != 0.0 && thetaP.back() == 0.0))) {
      std::cout << "You have further specified that the diffusion cannot be "
                   "absorbed at the boundary"
                << std::endl;
    } else {
      std::cout
          << "You have further specified that the diffusion cannot be "
             "absorbed at the boundary, but the provided mutation rates "
             "are strictly positive and thus already ensure this, so the "
             "resulting draws are coming from the *unconditioned* diffusion!"
          << std::endl;
    }
  }
  std::cout << "You've further specified the time threshold for Gaussian "
               "approximations at "
            << diffusion_threshold
            << ", whilst the bridge approximations threshold was set to "
            << bridge_threshold << std::endl;
  std::cout << "The pointwise computation will be performed over a mesh "
               "consisting of "
            << meshSize << "equally spaced intervals on [0,1]" << std::endl;
  std::cout << "Output will be printed to file in " << Filename << std::endl;
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);

  int counter = 0;
  double100 timeInc = endT - startT, sampleInc = sampleT - startT,
            yinc = 1.0 / static_cast<double100>(meshSize), y, ycount = 0.1;
  while (counter <= meshSize) {
    if (counter == 0) {
      y = 0.0;
    } else if (counter == meshSize) {
      y = 1.0;
    } else {
      y += yinc;
    }
    if (non_neutral) {
      cerr << "Truncated density cannot be computed for non-neutral case due "
              "to presence of intractable quantities!";
    } else {
      if (Absorption) {
        if ((x > 0.0) && (x < 1.0)) {
          saveFile << y << " "
                   << UnconditionedBridgeDensity(x, z, y, sampleInc, timeInc, o)
                   << "\n";
        } else {
          if (!(x > y) && !(x < y)) {
            saveFile << y << " " << 1.0 << "\n";
          } else {
            saveFile << y << " " << 0.0 << "\n";
          }
        }
      } else {
        saveFile << y << " " << BridgeDensity(x, z, y, sampleInc, timeInc, o)
                 << "\n";
      }
    }

    if (y >= ycount) {
      std::cout << "Calculated density up to y = " << ycount << endl;
      ycount += 0.1;
    }
    counter++;
  }

  std::cout << "Density calculation complete." << endl;

  saveFile.close();
}
