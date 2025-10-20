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

void WrightFisher::ThetaSetter(size_t N_i)  /// Set theta depending on what
                                            /// thetaP is entered
{
  if (!thetaP[N_i].empty())  /// When thetaP is not empty, set theta to be the
                             /// L1 norm of thetaP
    theta[N_i] = accumulate(thetaP[N_i].begin(), thetaP[N_i].end(), 0.0);
  else
    theta[N_i] = 0.0;  /// Otherwise thetaP is empty and we set theta to be 0.0
}

void WrightFisher::ThetaResetter(size_t N_i)  /// For conditioned bridge ONLY -
                                              /// set thetaP empty to (2,2), or
                                              /// thetaP = (0,theta)/(theta,0)
                                              /// to (2,theta)/(theta,2)
{
  auto &tp = thetaP[N_i];
  if (tp.empty())  /// If no mutation present, we need to set both thetaP
                   /// entries to 2
  {
    tp.push_back(2.0);
    tp.push_back(2.0);

    theta[N_i] = 4.0;
  } else if (tp.size() >= 2 &&
             (tp[0] <= 0.0 ||
              tp[1] <= 0.0)) {  /// If one sided mutation, we need to set
                                /// corresponding parameter to 2
    if (tp[0] <= 0.0) tp[0] = 2.0;
    if (tp[1] <= 0.0) tp[1] = 2.0;
    theta[N_i] = std::accumulate(tp.begin(), tp.end(), 0.0);
  }
}

void WrightFisher::SelectionSetter(
    size_t N_i)  /// Setting up selection mechanisms
{
  if (!non_neutral) return;  /// If true, we don't need to run anything

  Polynomial GenicSel(-1.0, 1.0,
                      0.0);  /// Set up genic selection (which will
                             /// always be present) as Polynomial class
  double alpha1 = thetaP[N_i].empty() ? 0.0 : 0.5 * thetaP[N_i][0],
         alpha2 = thetaP[N_i].empty() ? 0.0 : -0.5 * theta[N_i];
  Polynomial Alpha(alpha2, alpha1);
  SelectionFunction.resize(changepts.size());
  PhiFunction.resize(changepts.size());
  AtildeFunction.resize(changepts.size());

  switch (SelectionSetup) {
    case 0: {  /// Genic Selection - sigma*x*(1-x)
      SelPolyDeg = 2;
      SelectionFunction[N_i].Copy(0.5 * sigma[N_i] * GenicSel);
      auto tmpPhi = 0.5 * (2.0 * 0.5 * sigma[N_i] * Alpha +
                           0.25 * sigma[N_i] * sigma[N_i] * GenicSel);
      Polynomial tmpA(0.5 * sigma[N_i], 0.0);
      PhiFunction[N_i].Copy(tmpPhi);
      AtildeFunction[N_i].Copy(tmpA);
      break;
    }
    case 1: {  /// Diploid Selection - sigma*x*(1-x)*(h+x*(1-2*h))
      Polynomial Eta(1.0 - 2.0 * dominanceParameter, dominanceParameter);
      SelectionFunction[N_i].Copy(0.5 * sigma[N_i] * GenicSel * Eta);
      auto tmpPhi =
          0.5 * (2.0 * 0.5 * sigma[N_i] * Alpha * Eta +
                 0.25 * sigma[N_i] * sigma[N_i] * GenicSel * Eta * Eta +
                 0.5 * sigma[N_i] * GenicSel * Eta.Derivative());
      Polynomial tmpA(0.25 * sigma[N_i] * (1.0 - 2.0 * dominanceParameter),
                      0.5 * sigma[N_i] * dominanceParameter, 0.0);
      PhiFunction[N_i].Copy(tmpPhi);
      AtildeFunction[N_i].Copy(tmpA);
      break;
    }
    case 2: {  /// General polynomial selection - sigma*x*(1-x)*eta(x)
      SelPolyDeg = static_cast<int>(selectionCoeffs.size()) - 1;
      Polynomial Eta(&selectionCoeffs[0], SelPolyDeg);
      auto tmpPhi =
          0.5 * (2.0 * 0.5 * sigma[N_i] * Alpha * Eta +
                 0.25 * sigma[N_i] * sigma[N_i] * GenicSel * Eta * Eta +
                 0.5 * sigma[N_i] * GenicSel * Eta.Derivative());
      vector<double> Acoefs(SelPolyDeg + 2);
      Acoefs[0] = 0.0;
      for (int k = 1; k <= SelPolyDeg + 1;
           ++k)  /// Scaling polynomials coefficients appropriately
        Acoefs[k] = (0.5 * sigma[N_i] * selectionCoeffs[k - 1]) / double(k);
      Polynomial tmpA(&Acoefs[0], SelPolyDeg + 1);
      SelectionFunction[N_i].Copy(0.5 * sigma[N_i] * GenicSel * Eta);
      PhiFunction[N_i].Copy(tmpPhi);
      AtildeFunction[N_i].Copy(tmpA);
      break;
    }
    default:  /// If values not within {0,1,2} complain!
      cout << "Please enter a valid selection setup!\n";
  }
  PhiSetter(N_i);
}

void WrightFisher::PhiSetter(size_t N_i)  /// Compute max and min of Phi &
                                          /// Atilde functions
{
  if (!non_neutral) return;  /// If true, we don't need to bother
  phiMin.resize(changepts.size());
  phiMax.resize(changepts.size());
  AtildeMin.resize(changepts.size());
  AtildeMax.resize(changepts.size());
  if (SelectionSetup == 0)  /// Genic selection case
  {
    auto mm = PhitildeMinMaxRange(N_i);
    phiMin[N_i] = mm[0];
    phiMax[N_i] = mm[1];
    AtildeMax[N_i] = 0.5 * sigma[N_i];
    AtildeMin[N_i] = 0.0;
  } else  /// Otherwise we set Phi & Atilde fns, then use the
          /// PolynomialRootFinder to compute their extrema
  {
    Polynomial temp1 = PhiFunction[N_i].Derivative(),
               temp2 = AtildeFunction[N_i].Derivative();
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
    double potMax = max(PhiFunction[N_i].EvaluateReal(0.0),
                        PhiFunction[N_i].EvaluateReal(1.0)),
           potMin = min(PhiFunction[N_i].EvaluateReal(0.0),
                        PhiFunction[N_i].EvaluateReal(1.0));
    double AMax = max(AtildeFunction[N_i].EvaluateReal(0.0),
                      AtildeFunction[N_i].EvaluateReal(1.0));
    double AMin = min(AtildeFunction[N_i].EvaluateReal(0.0),
                      AtildeFunction[N_i].EvaluateReal(1.0));
    for (vector<double>::iterator vRit = validRoots.begin();
         vRit != validRoots.end(); vRit++) {
      potMax = max(potMax, PhiFunction[N_i].EvaluateReal(*vRit));
      potMin = min(potMin, PhiFunction[N_i].EvaluateReal(*vRit));
    }
    phiMax[N_i] = potMax;
    phiMin[N_i] = potMin;
    for (vector<double>::iterator vRit = Aroots.begin(); vRit != Aroots.end();
         vRit++) {
      AMax = max(AMax, AtildeFunction[N_i].EvaluateReal(*vRit));
      AMin = min(AMin, AtildeFunction[N_i].EvaluateReal(*vRit));
    }
    AtildeMax[N_i] = AMax;
    AtildeMin[N_i] = AMin;
  }
}

vector<vector<double>> WrightFisher::get_Theta() { return thetaP; }

double WrightFisher::Phitilde(
    size_t N_i,
    double y)  /// Returns value of the quadratic function phitilde(y)
{
  assert((y >= 0.0) && (y <= 1.0));
  if (SelectionSetup == 0)  /// For genic selection
  {
    if (sigma[N_i] == 0.0) {
      return 0.0;
    } else {
      return ((0.25) * sigma[N_i] *
              (-(0.5) * sigma[N_i] * y * y +
               ((0.5 * sigma[N_i]) - theta[N_i]) * y +
               (thetaP[N_i].empty() ? 0.0 : thetaP[N_i][0])));
    }
  } else  /// Otherwise use Polynomial class to evaluate
  {
    return PhiFunction[N_i].EvaluateReal(y);
  }
}

vector<double> WrightFisher::PhitildeMinMaxRange(
    size_t N_i)  /// Returns the minimum, maximum and
                 /// range of phitilde respectively for
                 /// genic selection case
{
  double phiargmin,
      phiargmax = max(min(0.5 - (theta[N_i] / sigma[N_i]), 1.0),
                      0.0);  /// Find out max value of phitilde by plugging in
  /// derivative(phitilde)=0

  if (!(phiargmax > 0.0))  /// Figure out if max is at either boundary and set
                           /// min to be other
  {
    phiargmin = 1.0;
  } else if (!(phiargmax < 1.0)) {
    phiargmin = 0.0;
  } else  /// Otherwise we need to choose the smaller value at each endpoints
          /// as the min
  {
    double phizero = Phitilde(N_i, 0.0), phione = Phitilde(N_i, 1.0);

    if (min(phizero, phione) == phizero) {
      phiargmin = 0.0;
    } else {
      phiargmin = 1.0;
    }
  }

  vector<double> MinMaxRange;

  MinMaxRange.push_back(Phitilde(N_i, phiargmin));
  MinMaxRange.push_back(Phitilde(N_i, phiargmax));
  MinMaxRange.push_back(MinMaxRange[1] - MinMaxRange[0]);

  return MinMaxRange;
}

double WrightFisher::Atilde(size_t N_i,
                            double x)  /// Returns the value of Atilde(x)
{
  assert((x >= 0.0) && (x <= 1.0));
  if (SelectionSetup == 0)  /// Genic selection case
  {
    return 0.5 * (sigma[N_i] * x);
  } else  /// Otherwise use Polynomial class
  {
    return AtildeFunction[N_i].EvaluateReal(x);
  }
}

double WrightFisher::Atildeplus(size_t N_i)  /// Returns the max of Atilde
{
  return AtildeMax[N_i];
}

double WrightFisher::Atildeminus(size_t N_i)  /// Returns the min of Atilde
{
  return AtildeMin[N_i];
}

size_t WrightFisher::first_greater_index(double s) {
  // REQUIRE: times is sorted ascending
  auto it = upper_bound(changepts.begin(), changepts.end(), s);  // first > s
  if (it == changepts.end())
    return static_cast<size_t>((changepts.end() - 1) -
                               changepts.begin());  // none found
  return static_cast<size_t>(it - changepts.begin());
}

size_t WrightFisher::first_ge_index(double s) {
  auto it = lower_bound(changepts.begin(), changepts.end(), s);  // first >= s
  if (it == changepts.end())
    return static_cast<size_t>((changepts.end() - 1) - changepts.begin());
  return static_cast<size_t>(it - changepts.begin());
}

size_t WrightFisher::getIndex(double s) {
  auto it = lower_bound(changepts.begin(), changepts.end(), s);  // first >= x
  if (it == changepts.begin()) return static_cast<size_t>(0);    // none < x
  return static_cast<size_t>(distance(changepts.begin(), it) - 1);
}

pair<double, double> WrightFisher::GriffithsParas(
    size_t N_i,
    double t)  /// Compute parameters for Griffiths approximation
{
  assert(t > 0.0);
  double beta = (theta[N_i] - 1.0) * static_cast<double>(t) / 2.0;
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

int WrightFisher::radiate_from_mode(size_t N_i, int index, const double t)
    const  /// Moving from mode as opposed to starting from m = 0 and working
           /// upwards
{
  double beta = (theta[N_i] - 1.0) * static_cast<double>(t) / 2.0;
  double eta = (abs(beta) <= 2.5e-5 ? 1.0 : beta / (exp(beta) - 1.0));
  int mmode = static_cast<int>(round(2 * eta / static_cast<double>(t))),
      threshold = (thetaP[N_i].empty() ? 1 : 0), adjindex = index + threshold;
  if (adjindex > 2 * (mmode)-threshold) {
    return adjindex;
  } else {
    return mmode + (index % 2 == 0 ? 1 : -1) * ((index + 1) / 2);
  }
}

void WrightFisher::increment_on_mk(bool isDiffTheta, size_t N_i,
                                   vector<int> &mk, const double s,
                                   const double t)
    const  /// Incrementing (m,k) in bridge sampler using bijective fn
{
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  int &m_index = mk[2], &k_index = mk[3];
  --m_index;
  ++k_index;
  if (m_index < 0) {
    m_index = m_index + k_index + 1;
    k_index = 0;
  }
  mk[0] = radiate_from_mode(N_i_s, m_index, s);
  mk[1] = radiate_from_mode(N_i_t, k_index, t - s);
}

vector<int> WrightFisher::radiate2d_mode(bool isDiffTheta, size_t N_i,
                                         size_t idx, double s, double t) const {
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  double beta_m = (theta[N_i_s] - 1.0) * static_cast<double>(s) / 2.0;
  double eta_m = (abs(beta_m) <= 2.5e-5 ? 1.0 : beta_m / (exp(beta_m) - 1.0));
  int mmode = static_cast<int>(round(2 * eta_m / static_cast<double>(s)));
  double beta_k = (theta[N_i_t] - 1.0) * static_cast<double>(t - s) / 2.0;
  double eta_k = (abs(beta_k) <= 2.5e-5 ? 1.0 : beta_k / (exp(beta_k) - 1.0));
  int kmode = static_cast<int>(round(2 * eta_k / static_cast<double>(t - s)));
  if (idx == 0) return {mmode, kmode};

  // find ring r >= 1 such that (2r+1)^2 > idx
  int r = int(std::floor((std::sqrt(double(idx)) - 1.0) / 2.0)) + 1;
  int prev_count = (2 * (r - 1) + 1) * (2 * (r - 1) + 1);
  int offset = int(idx - prev_count - 1);  // position within this ring
  int side_len = 2 * r;
  int side = offset / side_len;  // 0=right,1=top,2=left,3=bottom
  int pos = offset % side_len;

  int dx, dy;
  switch (side) {
    case 0:  // right edge
      dx = r;
      dy = -r + 1 + pos;
      break;
    case 1:  // top edge
      dx = r - 1 - pos;
      dy = r;
      break;
    case 2:  // left edge
      dx = -r;
      dy = r - 1 - pos;
      break;
    default:  // bottom edge
      dx = -r + 1 + pos;
      dy = -r;
      break;
  }
  vector<int> ret_vec = {mmode + dx, kmode + dy};
  return ret_vec;
}

double WrightFisher::Getd(
    size_t N_i, vector<double> &d, int i, double x, double z,
    double t)  /// Compute contributions to denom for DrawBridgePMF cases
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset, c2 = m - j;
      boost::math::binomial_distribution<double> B(c2, x);
      double Expected_Dirichlet = 0.0;
      for (int l = 0; l <= c2; ++l) {
        double para1 = (thetaP[N_i].empty()
                            ? static_cast<double>(l)
                            : static_cast<double>(thetaP[N_i][0] + l)),
               para2 = (thetaP[N_i].empty()
                            ? static_cast<double>(c2 - l)
                            : static_cast<double>(thetaP[N_i][1] + c2 - l));

        boost::math::beta_distribution<double> D(para1, para2);
        Expected_Dirichlet += pdf(B, l) * pdf(D, z);
      }

      d[i] += exp(Getlogakm<double>(N_i, c1, c2) +
                  static_cast<double>(-c1 * (c1 + theta[N_i] - 1) * t / 2.0)) *
              Expected_Dirichlet;
      assert(exp(Getlogakm<double>(N_i, c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double WrightFisher::GetdDiffTheta(
    size_t N_i, vector<double> &d, int i, double x, double z, double s,
    double t,
    const Options
        &o)  /// Compute contributions to denom for DrawBridgePMFDiffTheta cases
{
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  if (s <= o.g1984threshold &&
      (t - s) <= o.g1984threshold) {  // Only used by density calculator!
    /// Compute denominator
    double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-50;
    int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
    int dmmode = static_cast<int>(ceil(GriffithsParas(N_i_s, s).first)),
        dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
    bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
    int dkmode = static_cast<int>(ceil(GriffithsParas(N_i_t, t - s).first));

    while (!dmSwitch) {
      emCOldInc = emCInc;
      double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
      int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
      bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
      while (!dkSwitch) {
        ekCOldInc = ekCInc;
        double addon = 0.0;
        precomputeA(N_i_s, N_i_t, dm, dk);
        double extra_factors =
            factorials[dm] + factorials[dk] + static_cast<double>(dm) * log1x +
            static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
            static_cast<double>(thetaP[N_i_t][1] + dk - 1.0) * log1z +
            lg_theta[N_i_s][dm] + lg_theta[N_i_t][dk] -
            lg_theta[N_i_s][dm + dk];
        for (int l = 0; l != dm; l++) {
          for (int j = 0; j != dk; j++) {
            addon +=
                exp(extra_factors - factorials[l] - factorials[dm - l] -
                    factorials[j] - factorials[dk - j] - lg_theta1[N_i_s][l] -
                    lg_theta2[N_i_s][dm - l] - lg_theta1[N_i_t][j] -
                    lg_theta2[N_i_t][dk - j] + lg_theta1[N_i_s][l + j] +
                    lg_theta2[N_i_s][dm - l + dk - j] +
                    static_cast<double>(l) * (logx - log1x) +
                    static_cast<double>(j) * (logz - log1z));
          }
        }
        ekCInc = QmApprox(N_i_t, dk, t - s, o) * addon;
        ekC += ekCInc;
        if (!(dkDownSwitch))  /// Switching mechanism for j
        {
          if (sgn(dk - dkmode) <= 0) {
            dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
          }
        }

        if (!(dkUpSwitch)) {
          if (sgn(dk - dkmode) >= 0) {
            dkUpSwitch = (ekCInc < eC_threshold);
          }
        }

        dkSwitch = (dkDownSwitch && dkUpSwitch);
        if (!dkSwitch) {
          if (dkFlip == 1) {
            dkU++;
            dk = dkmode + dkU;
            dkFlip *= (dkDownSwitch ? 1 : -1);
          } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
            dkD++;
            dk = dkmode - dkD;
            dkFlip *= (dkUpSwitch ? 1 : -1);
          }
        }
      }
      emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
      eC += emCInc;

      if (!(dmDownSwitch))  /// Switching mechanism for j
      {
        if (sgn(dm - dmmode) <= 0) {
          dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
        }
      }

      if (!(dmUpSwitch)) {
        if (sgn(dm - dmmode) >= 0) {
          dmUpSwitch = (emCInc < eC_threshold);
        }
      }

      dmSwitch = (dmDownSwitch && dmUpSwitch);
      if (!dmSwitch) {
        if (dmFlip == 1) {
          dmU++;
          dm = dmmode + dmU;
          dmFlip *= (dmDownSwitch ? 1 : -1);
        } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
          dmD++;
          dm = dmmode - dmD;
          dmFlip *= (dmUpSwitch ? 1 : -1);
        }
      }
    }
    return eC;
  } else {  // Used in all simulation routines
    if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
    if (d[i] < 0.0) {
      d[i] = 0.0;
      int m = i / 2, offset = i % 2;

      for (int j = 0; j <= m; ++j) {
        int c1 = m + j + offset, c2 = m - j;
        d[i] += computeAGammaLambda(N_i, c1, c2, x, z, s, t, o);
      }
    }
    assert(d[i] >= 0.0);
    return d[i];
  }
}

double WrightFisher::GetdDiffThetaBoundaries(
    size_t N_i, vector<double> &d, int i, double x, double z, double s,
    double t,
    const Options
        &o)  /// Compute contributions to denom for DrawBridgePMF cases
{
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  if (s <= o.g1984threshold &&
      (t - s) <= o.g1984threshold) {  // Only used by density calculator!
    /// Compute denominator
    double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-50;
    int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
    int dmmode = static_cast<int>(ceil(GriffithsParas(N_i_s, s).first)),
        dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
    bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
    int dkmode = static_cast<int>(ceil(GriffithsParas(N_i_t, t - s).first));

    while (!dmSwitch) {
      emCOldInc = emCInc;
      double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
      int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
      bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
      while (!dkSwitch) {
        ekCOldInc = ekCInc;
        double addon = 0.0;
        precomputeA(N_i_s, N_i_t, dm, dk);
        double beta1, beta2, beta3;
        if (!(x > 0.0)) {
          if (!(z > 0.0)) {  // x = 0 && z = 0
            beta1 = lg_theta1[N_i_s][0] + lg_theta2[N_i_s][dm + dk] -
                    lg_theta[N_i_s][dm + dk];
            beta2 = lg_theta[N_i_s][dm] - lg_theta1[N_i_s][0] -
                    lg_theta2[N_i_s][dm];
            beta3 = lg_theta[N_i_t][dk] - lg_theta1[N_i_t][0] -
                    lg_theta2[N_i_t][dk];
          } else {  // x = 0 && z = 1
            beta1 = lg_theta1[N_i_s][dk] + lg_theta2[N_i_s][dm] -
                    lg_theta[N_i_s][dm + dk];
            beta2 = lg_theta[N_i_s][dm] - lg_theta1[N_i_s][0] -
                    lg_theta2[N_i_s][dm];
            beta3 = lg_theta[N_i_t][dk] - lg_theta1[N_i_t][dk] -
                    lg_theta2[N_i_t][0];
          }
        } else {
          if (!(z > 0.0)) {  // x = 1 && z = 0
            beta1 = lg_theta1[N_i_s][dm] + lg_theta2[N_i_s][dk] -
                    lg_theta[N_i_s][dm + dk];
            beta2 = lg_theta[N_i_s][dm] - lg_theta1[N_i_s][dm] -
                    lg_theta2[N_i_s][0];
            beta3 = lg_theta[N_i_t][dk] - lg_theta1[N_i_t][0] -
                    lg_theta2[N_i_t][dk];
          } else {  // x = 1 && z = 1
            beta1 = lg_theta1[N_i_s][dk + dm] + lg_theta2[N_i_s][0] -
                    lg_theta[N_i_s][dm + dk];
            beta2 = lg_theta[N_i_s][dm] - lg_theta1[N_i_s][dm] -
                    lg_theta2[N_i_s][0];
            beta3 = lg_theta[N_i_t][dk] - lg_theta1[N_i_t][dk] -
                    lg_theta2[N_i_t][0];
          }
        }
        addon += exp(beta1 + beta2 + beta3);
        ekCInc = QmApprox(N_i_t, dk, t - s, o) * addon;
        ekC += ekCInc;
        if (!(dkDownSwitch))  /// Switching mechanism for k
        {
          if (sgn(dk - dkmode) <= 0) {
            dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
          }
        }

        if (!(dkUpSwitch)) {
          if (sgn(dk - dkmode) >= 0) {
            dkUpSwitch = (ekCInc < eC_threshold);
          }
        }

        dkSwitch = (dkDownSwitch && dkUpSwitch);
        if (!dkSwitch) {
          if (dkFlip == 1) {
            dkU++;
            dk = dkmode + dkU;
            dkFlip *= (dkDownSwitch ? 1 : -1);
          } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
            dkD++;
            dk = dkmode - dkD;
            dkFlip *= (dkUpSwitch ? 1 : -1);
          }
        }
      }
      emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
      eC += emCInc;

      if (!(dmDownSwitch))  /// Switching mechanism for m
      {
        if (sgn(dm - dmmode) <= 0) {
          dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
        }
      }

      if (!(dmUpSwitch)) {
        if (sgn(dm - dmmode) >= 0) {
          dmUpSwitch = (emCInc < eC_threshold);
        }
      }

      dmSwitch = (dmDownSwitch && dmUpSwitch);
      if (!dmSwitch) {
        if (dmFlip == 1) {
          dmU++;
          dm = dmmode + dmU;
          dmFlip *= (dmDownSwitch ? 1 : -1);
        } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
          dmD++;
          dm = dmmode - dmD;
          dmFlip *= (dmUpSwitch ? 1 : -1);
        }
      }
    }
    return eC;
  } else {  // Used in all simulation routines
    if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
    if (d[i] < 0.0) {
      d[i] = 0.0;
      int m = i / 2, offset = i % 2;

      for (int j = 0; j <= m; ++j) {
        int c1 = m + j + offset, c2 = m - j;
        d[i] += computeAGammaLambda(N_i, c1, c2, x, z, s, t, o);
      }
    }
    assert(d[i] >= 0.0);
    return d[i];
  }
}

double WrightFisher::GetdDiffThetaOneQApprox(
    size_t N_i, vector<double> &d, int i, double x, double z, double s,
    double t,
    const Options
        &o)  /// Compute denom when we have different theta but using OneQApprox
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  size_t N_i_approx = (s < (t - s)) ? N_i : N_i + 1;
  size_t N_i_noapprox = (N_i_approx == N_i) ? N_i + 1 : N_i;
  double tapprox = (s < t - s) ? s : t - s;
  double tnoapprox = (tapprox == s) ? t - s : s;
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset, c2 = m - j;
      double expectation_term =
          calculate_expectation_qApprox(N_i, c2, x, z, s, t, o);

      d[i] += exp(Getlogakm<double>(N_i_noapprox, c1, c2) +
                  static_cast<double>(-c1 * (c1 + theta[N_i_noapprox] - 1) *
                                      tnoapprox / 2.0) +
                  log(expectation_term));
      assert(exp(Getlogakm<double>(N_i_noapprox, c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double WrightFisher::calculate_expectation_qApprox(
    size_t N_i, int k, double x, double z, double s, double t,
    const Options &o) {  /// Calculate expectation by including discretised
                         /// normal as part of inner expectation
  size_t N_i_approx = (s < t - s) ? N_i : N_i + 1;
  double tapprox = (s < t - s) ? s : t - s;
  pair<double, double> paras_m = GriffithsParas(N_i_approx, tapprox);
  int lower_m_ind = max(int(floor(paras_m.first - 5.0 * paras_m.second)), 0);
  int upper_m_ind = ceil(paras_m.first + 5.0 * paras_m.second);
  double expectation = 0.0;
  for (int m_ind = lower_m_ind; m_ind <= upper_m_ind; m_ind++) {
    double exp_term = (s < t - s)
                          ? log(calculate_expectation(N_i, m_ind, k, x, z))
                          : log(calculate_expectation(N_i, k, m_ind, x, z));
    expectation += exp(log(QmApprox(N_i_approx, m_ind, tapprox, o)) + exp_term);
  }
  return expectation;
}

double WrightFisher::computeAGammaLambda(
    size_t N_i, int gamma, int lambda, double x, double z, double s, double t,
    const Options &o)  /// Function to compute A_{gamma, lambda}
{
  double A_gamma_lambda = 0.0;
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  if (s > o.g1984threshold && t - s > o.g1984threshold) {
    for (int m = 0; m <= lambda; m++) {
      double add_on = 0.0;
      for (int h = m; h <= m + gamma - lambda; h++) {
        add_on +=
            exp(Getlogakm<double>(N_i_s, h, m) +
                static_cast<double>(-h * (h + theta[N_i_s] - 1) * s / 2.0) +
                Getlogakm<double>(N_i_t, gamma - h, lambda - m) +
                static_cast<double>(-(gamma - h) *
                                    ((gamma - h) + theta[N_i_t] - 1) * (t - s) /
                                    2.0));
      }
      A_gamma_lambda +=
          (add_on * calculate_expectation(N_i, m, lambda - m, x, z));
    }
  } else {  // Only one of the intervals is small, otherwise we would have
            // invoked DrawBridgePMFDiffThetaApprox function, and this wouldn't
            // be called anyway
    if (s <= o.g1984threshold) {  // Then the small interval is over [0,s) and
                                  // we approximate q_m^{theta}(s) via
                                  // discretised normals
      for (int m = 0; m <= lambda; m++) {
        double add_on = 0.0;
        for (int h = m; h <= m + gamma - lambda; h++) {
          add_on += exp(Getlogakm<double>(N_i_t, gamma - h, lambda - m) +
                        static_cast<double>(-(gamma - h) *
                                            ((gamma - h) + theta[N_i_t] - 1) *
                                            (t - s) / 2.0));
        }
        A_gamma_lambda += QmApprox(N_i_s, m, s, o) * add_on *
                          calculate_expectation(N_i, m, lambda - m, x, z);
      }
    } else {
      for (int m = 0; m <= lambda; m++) {
        A_gamma_lambda +=
            (exp(Getlogakm<double>(N_i_s, gamma, m) +
                 static_cast<double>(-gamma * (gamma + theta[N_i_s] - 1) * s /
                                     2.0)) *
             QmApprox(N_i_t, lambda - m, t - s, o) *
             calculate_expectation(N_i, m, lambda - m, x, z));
      }
    }
  }
  return A_gamma_lambda;
}

double WrightFisher::calculate_expectation(
    size_t N_i, int m, int k, double x,
    double z) {  /// Function to compute expectation present in inner usms for
                 /// different thetas
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  precomputeA(N_i_s, N_i_t, m, k);
  double beta1, beta2, beta3, R = 0.0;
  if (!(x > 0.0)) {
    if (!(z > 0.0)) {  // x = 0 && z = 0
      beta1 = lg_theta1[N_i_s][0] + lg_theta2[N_i_s][m + k] -
              lg_theta[N_i_s][m + k];
      beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m];
      beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k];
      R = exp(beta1 + beta2 + beta3);
    } else if (!(z < 1.0)) {  // x = 0 && z = 1
      beta1 =
          lg_theta1[N_i_s][k] + lg_theta2[N_i_s][m] - lg_theta[N_i_s][m + k];
      beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m];
      beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0];
      R = exp(beta1 + beta2 + beta3);
    } else {  // x = 0 && z in (0,1)
      double constant_terms =
          factorials[k] - lg_theta[N_i_s][m + k] + lg_theta[N_i_s][m] -
          lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m] + lg_theta[N_i_t][k] +
          static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
          static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z;
      for (int j = 0; j <= k; j++) {
        double term = -factorials[j] - factorials[k - j] + lg_theta1[N_i_s][j] +
                      lg_theta2[N_i_s][m + k - j] - lg_theta1[N_i_t][j] -
                      lg_theta2[N_i_t][k - j] +
                      static_cast<double>(j) * (logz - log1z);
        R += exp(term);
      }
      R *= exp(constant_terms);
    }
  } else if (!(x < 1.0)) {
    if (!(z > 0.0)) {  // x = 1 && z = 0
      beta1 =
          lg_theta1[N_i_s][m] + lg_theta2[N_i_s][k] - lg_theta[N_i_s][m + k];
      beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0];
      beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k];
      R = exp(beta1 + beta2 + beta3);
    } else if (!(z < 1.0)) {  // x = 1 && z = 1
      beta1 = lg_theta1[N_i_s][k + m] + lg_theta2[N_i_s][0] -
              lg_theta[N_i_s][m + k];
      beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0];
      beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0];
      R = exp(beta1 + beta2 + beta3);
    } else {  // x = 1 && z in (0,1)
      double constant_terms =
          factorials[k] - lg_theta[N_i_s][m + k] + lg_theta[N_i_s][m] -
          lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0] + lg_theta[N_i_t][k] +
          static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
          static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z;
      for (int j = 0; j <= k; j++) {
        double term = -factorials[j] - factorials[k - j] +
                      lg_theta1[N_i_s][m + j] + lg_theta2[N_i_s][k - j] -
                      lg_theta1[N_i_t][j] - lg_theta2[N_i_t][k - j] +
                      static_cast<double>(j) * (logz - log1z);
        R += exp(term);
      }
      R *= exp(constant_terms);
    }
  } else {
    if (!(z > 0.0)) {  // x in (0,1), z = 0
      double constant_terms = factorials[m] - lg_theta[N_i_s][m + k] +
                              lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                              lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k] +
                              static_cast<double>(m) * log1x;
      for (int l = 0; l <= m; l++) {
        double term = -factorials[l] - factorials[m - l] + lg_theta1[N_i_s][l] +
                      lg_theta2[N_i_s][m - l + k] - lg_theta1[N_i_s][l] -
                      lg_theta2[N_i_s][m - l] +
                      static_cast<double>(l) * (logx - log1x);
        R += exp(term);
      }
      R *= exp(constant_terms);
    } else if (!(z < 1.0)) {  // x in (0,1), z = 1
      double constant_terms = factorials[m] - lg_theta[N_i_s][m + k] +
                              lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                              lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0] +
                              static_cast<double>(m) * log1x;
      for (int l = 0; l <= m; l++) {
        double term = -factorials[l] - factorials[m - l] +
                      lg_theta1[N_i_s][l + k] + lg_theta2[N_i_s][m - l] -
                      lg_theta1[N_i_s][l] - lg_theta2[N_i_s][m - l] +
                      static_cast<double>(l) * (logx - log1x);
        R += exp(term);
      }
      R *= exp(constant_terms);
    } else {  //  x in (0,1), z in (0,1)
      // outer loop over l
      for (int l = 0; l <= m; ++l) {
        // compute the inner sum in j
        for (int j = 0; j <= k; ++j) {
          double term = static_cast<double>(l) * (logx - log1x) +
                        static_cast<double>(j) * (logz - log1z) -
                        factorials[l] - factorials[m - l] - factorials[j] -
                        factorials[k - j] + lg_theta1[N_i_s][l + j] +
                        lg_theta2[N_i_s][m - l + k - j] - lg_theta1[N_i_s][l] -
                        lg_theta1[N_i_t][j] - lg_theta2[N_i_s][m - l] -
                        lg_theta2[N_i_t][k - j];
          R += exp(term);
        }
      }
      R *=
          exp(factorials[m] + factorials[k] + static_cast<double>(m) * log1x +
              (thetaP[N_i_t][0] - 1.0) * logz +
              (static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z) +
              lg_theta[N_i_s][m] + lg_theta[N_i_t][k] - lg_theta[N_i_s][m + k]);
    }
  }
  return R;
}

double WrightFisher::Getd2(size_t N_i, vector<double> &d, int i, double x,
                           double t)  /// Compute contributions to denom for
                                      /// AncestralProcessConditional function
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;
    int bumper = (thetaP[N_i].empty()) ? 2 : 1;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + bumper + j + offset, c2 = m + bumper - j;
      double binomialremainder =
          (bumper == 2) ? (1.0 - pow(1.0 - x, static_cast<double>(c2)) -
                           pow(x, static_cast<double>(c2)))
                        : (1.0 - pow(1.0 - x, static_cast<double>(c2)));
      boost::math::binomial_distribution<double> B(
          c2,
          x);  // Could make these distributions static for the duration of
      // this x and z
      d[i] += exp(Getlogakm<double>(N_i, c1, c2) +
                  static_cast<double>(-c1 * (c1 + theta[N_i] - 1) * t / 2.0)) *
              binomialremainder;
      assert(exp(Getlogakm<double>(N_i, c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double WrightFisher::GetdBridgeSame(
    size_t N_i, vector<double> &d, int i, double x,
    double t)  /// Compute contributions to denom for DrawBridgePMFSame function
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset, c2 = m - j;
      double addon =
          (thetaP[N_i][0] > 0.0 && thetaP[N_i][1] > 0.0
               ? (!(x > 0.0) ? (log(boost::math::tgamma_ratio(
                                    static_cast<double>(theta[N_i] + c2),
                                    static_cast<double>(thetaP[N_i][1] + c2))) -
                                log(boost::math::tgamma(thetaP[N_i][0])))
                             : (log(boost::math::tgamma_ratio(
                                    static_cast<double>(theta[N_i] + c2),
                                    static_cast<double>(thetaP[N_i][0] + c2))) -
                                log(boost::math::tgamma(thetaP[N_i][1]))))
               : (log(boost::math::tgamma_ratio(
                      static_cast<double>(theta[N_i] + c2),
                      static_cast<double>(c2))) -
                  log(boost::math::tgamma(theta[N_i]))));
      d[i] += exp(Getlogakm<double>(N_i, c1, c2) +
                  static_cast<double>(-c1 * (c1 + theta[N_i] - 1) * t / 2.0) +
                  addon);
      assert(exp(Getlogakm<double>(N_i, c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double WrightFisher::GetdBridgeInterior(
    size_t N_i, vector<double> &d, int i, double x, double z,
    double t)  /// Compute contributions to denom for DrawBridgePMFInterior
               /// function
{
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset, c2 = m - j;
      double para1, para2;

      if (!(x > 0.0)) {
        para1 = thetaP[N_i][0];
        para2 = thetaP[N_i][1] + c2;
      } else {
        para1 = thetaP[N_i][0] + c2;
        para2 = thetaP[N_i][1];
      }

      boost::math::beta_distribution<double> BETA(para1, para2);
      double zcontribution = log(pdf(BETA, z));
      d[i] += exp(Getlogakm<double>(N_i, c1, c2) +
                  static_cast<double>(-c1 * (c1 + theta[N_i] - 1) * t / 2.0) +
                  zcontribution);
      assert(exp(Getlogakm<double>(N_i, c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double WrightFisher::GetdBridgeUnconditional(
    size_t N_i, vector<double> &d, int i, double x, double z,
    double t)  /// Compute contributions to denom for
               /// DrawBridgePMFUnconditioned function
{
  int thetaDep = (thetaP[N_i].empty() ? 1 : 0);
  if (i > static_cast<int>(d.size()) - 1) d.resize(i + 1, -1.0);
  if (d[i] < 0.0) {
    d[i] = 0.0;
    int m = i / 2, offset = i % 2;

    for (int j = 0; j <= m; ++j) {
      int c1 = m + j + offset + thetaDep, c2 = m - j + thetaDep;
      double xcontribution =
          (!(z > 0.0) ? static_cast<double>(c2) * log(1.0 - x)
                      : static_cast<double>(c2) * log(x));
      d[i] += exp(Getlogakm<double>(N_i, c1, c2) +
                  static_cast<double>(-c1 * (c1 + theta[N_i] - 1) * t / 2.0) +
                  xcontribution);
      assert(exp(Getlogakm<double>(N_i, c1, c2)) > 0.0);
    }
  }
  assert(d[i] >= 0.0);
  return d[i];
}

double WrightFisher::computeA(
    size_t N_i, int m, int k, int l, int j, double x,
    double z)  /// Compute weights for bridge diffusion decomposition
{
  assert((m >= 0) && (k >= 0) && (l >= 0) && (j >= 0) && (m >= l) && (k >= j) &&
         (x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0));
  double theta1 = thetaP[N_i][0], theta2 = thetaP[N_i][1];
  boost::math::binomial_distribution<double> BIN(m, x);
  boost::math::beta_distribution<double> BETA(theta1 + j, theta2 + k - j);

  double betaBinPdfs;

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

  double log_bin = LogBinomialCoefficientCalculator(k, j);
  double log_beta_1 =
      log(boost::math::beta<double>(theta1 + l + j, theta2 + m - l + k - j));
  double log_beta_2 =
      log(boost::math::beta<double>(theta1 + l, theta2 + m - l));

  return betaBinPdfs * exp(log_bin + log_beta_1 - log_beta_2);
}

double WrightFisher::computeAUnconditional(
    size_t N_i, int m, int k, int l, double x,
    double z)  /// Compute weights for unconditioned bridge diffusion
               /// decomposition
{
  assert((m >= 0) && (k >= 0) && (l >= 0) && (m >= l) && (x > 0.0) &&
         (x < 1.0) && (!(z > 0.0) || !(z < 1.0)));
  boost::math::binomial_distribution<double> BIN(m, x);

  if ((l == 0) && (thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0))) {
    return pdf(BIN, l);
  } else if ((l == m) && (thetaP[N_i].empty() || !(thetaP[N_i][1] > 0.0))) {
    return pdf(BIN, l);
  }

  if (!(z > 0.0)) {
    return pdf(BIN, l) *
           exp(boost::math::lgamma(
                   static_cast<double>(theta[N_i] + m - l + k)) +
               boost::math::lgamma(static_cast<double>(theta[N_i] + m)) -
               boost::math::lgamma(static_cast<double>(theta[N_i] + m + k)) -
               boost::math::lgamma(static_cast<double>(theta[N_i] + m - l)));
  } else if (!(z < 1.0)) {
    return pdf(BIN, l) *
           exp(boost::math::lgamma(static_cast<double>(theta[N_i] + l + k)) +
               boost::math::lgamma(static_cast<double>(theta[N_i] + m)) -
               boost::math::lgamma(static_cast<double>(theta[N_i] + l)) -
               boost::math::lgamma(static_cast<double>(theta[N_i] + k + m)));
  } else {
    cerr << "z was not 0 or 1! Returning 0.";
    return 0.0;
  }
}

void WrightFisher::precomputeA(
    size_t N_i_s, size_t N_i_t, int m,
    int k) {  /// Function to precompute all relevant terms for A_{m,k,l,j}
  if ((int)factorials.size() <= max(m, k)) {
    int old = factorials.size();
    factorials.resize(max(m, k) + 1);
    if (old == 0) {
      factorials[0] = 0.0;
      old = 1;
    }
    for (int i = old; i <= max(m, k); ++i) {
      factorials[i] = factorials[i - 1] + std::log(double(i));
    }
  }

  if (N_i_s != N_i_t) {
    if ((int)lg_theta1.size() <= N_i_t) {
      lg_theta1.resize(N_i_t + 1);
      lg_theta2.resize(N_i_t + 1);
      lg_theta.resize(N_i_t + 1);
    }
  } else {
    if ((int)lg_theta1.size() <= N_i_s) {
      lg_theta1.resize(N_i_s + 1);
      lg_theta2.resize(N_i_s + 1);
      lg_theta.resize(N_i_s + 1);
    }
  }

  // 3) lg_theta1[N_i_s] / lg_theta2[N_i_s] up to m + k
  if ((int)lg_theta1[N_i_s].size() <= m + k) {
    int old = lg_theta1[N_i_s].size();
    lg_theta1[N_i_s].resize(m + k + 1);
    lg_theta2[N_i_s].resize(m + k + 1);
    // base
    if (old == 0) {
      lg_theta1[N_i_s][0] = boost::math::lgamma(thetaP[N_i_s][0]);
      lg_theta2[N_i_s][0] = boost::math::lgamma(thetaP[N_i_s][1]);
      old = 1;
    }
    for (int i = old; i <= m + k; ++i) {
      lg_theta1[N_i_s][i] =
          lg_theta1[N_i_s][i - 1] + std::log(thetaP[N_i_s][0] + (i - 1));
      lg_theta2[N_i_s][i] =
          lg_theta2[N_i_s][i - 1] + std::log(thetaP[N_i_s][1] + (i - 1));
    }
  }

  if (N_i_s != N_i_t) {
    // 4) Repeat for N_i_t
    if ((int)lg_theta1[N_i_t].size() <= m + k) {
      int old = lg_theta1[N_i_t].size();
      lg_theta1[N_i_t].resize(m + k + 1);
      lg_theta2[N_i_t].resize(m + k + 1);
      // base
      if (old == 0) {
        lg_theta1[N_i_t][0] = boost::math::lgamma(thetaP[N_i_t][0]);
        lg_theta2[N_i_t][0] = boost::math::lgamma(thetaP[N_i_t][1]);
        old = 1;
      }
      for (int i = old; i <= m + k; ++i) {
        lg_theta1[N_i_t][i] =
            lg_theta1[N_i_t][i - 1] + std::log(thetaP[N_i_t][0] + (i - 1));
        lg_theta2[N_i_t][i] =
            lg_theta2[N_i_t][i - 1] + std::log(thetaP[N_i_t][1] + (i - 1));
      }
    }
  }

  // 4) Repeat for lg_theta
  if ((int)lg_theta[N_i_s].size() <= m + k) {
    int old = lg_theta[N_i_s].size();
    lg_theta[N_i_s].resize(m + k + 1);
    // base
    if (old == 0) {
      lg_theta[N_i_s][0] = boost::math::lgamma(theta[N_i_s]);
      old = 1;
    }
    for (int i = old; i <= m + k; ++i) {
      lg_theta[N_i_s][i] =
          lg_theta[N_i_s][i - 1] + std::log(theta[N_i_s] + (i - 1));
    }
  }

  if (N_i_s != N_i_t) {
    if ((int)lg_theta[N_i_t].size() <= m + k) {
      int old = lg_theta[N_i_t].size();
      lg_theta[N_i_t].resize(m + k + 1);
      // base
      if (old == 0) {
        lg_theta[N_i_t][0] = boost::math::lgamma(theta[N_i_t]);
        old = 1;
      }
      for (int i = old; i <= m + k; ++i) {
        lg_theta[N_i_t][i] =
            lg_theta[N_i_t][i - 1] + std::log(theta[N_i_t] + (i - 1));
      }
    }
  }
}

int WrightFisher::computeC(
    size_t N_i, int m,
    pair<vector<int>, double> &C)  /// Compute the quantity C_m
{
  assert(m >= 0);
  if (m > static_cast<int>(C.first.size() - 1) || C.first[m] < 0) {
    C.first.resize(m + 1, -1);
    int i = 0;
    double bimnext =
               exp(Getlogakm<double>(N_i, i + m + 1, m) -
                   (i + m + 1) * (i + m + 1 + theta[N_i] - 1) * C.second / 2.0),
           bim = exp(Getlogakm<double>(N_i, i + m, m) -
                     (i + m) * (i + m + theta[N_i] - 1) * C.second / 2.0);
    while (bimnext > bim) {
      ++i;
      bim = bimnext;
      bimnext =
          exp(Getlogakm<double>(N_i, i + m + 1, m) -
              (i + m + 1) * (i + m + 1 + theta[N_i] - 1) * C.second / 2.0);
    }
    C.first[m] = i;
  }
  assert(C.first[m] >= 0);
  return C.first[m];
}

int WrightFisher::computeF(
    size_t N_i, int m,
    pair<vector<int>, double> &C)  /// Compute the quantity F_m
{
  assert(m >= 0);
  if (m > static_cast<int>(C.first.size() - 1) || C.first[m] < 0) {
    C.first.resize(m + 1, -1);
    int i = 0;
    double bimnext =
               exp(Getlogakm<double>(N_i, i + m + 1, m) -
                   (i + m + 1) * (i + m + 1 + theta[N_i] - 1) * C.second / 2.0),
           bim = exp(Getlogakm<double>(N_i, i + m, m) -
                     (i + m) * (i + m + theta[N_i] - 1) * C.second / 2.0);
    while ((2.0 * bimnext) > bim) {
      ++i;
      bim = bimnext;
      bimnext =
          exp(Getlogakm<double>(N_i, i + m + 1, m) -
              (i + m + 1) * (i + m + 1 + theta[N_i] - 1) * C.second / 2.0);
    }
    C.first[m] = i;
  }
  assert(C.first[m] >= 0);
  return C.first[m];
}

int WrightFisher::computeE(
    size_t N_i, pair<vector<int>, double> &C)  /// Compute the quantity E
{
  int next_constraint = computeC(N_i, 0, C), curr_row = 0;
  int diag_index = next_constraint + (next_constraint % 2);
  bool Efound = false;
  while (!Efound) {
    ++curr_row;
    next_constraint = computeC(N_i, curr_row, C) + curr_row;
    if (diag_index - curr_row < next_constraint) {
      diag_index = next_constraint + curr_row;
      diag_index += (diag_index % 2);
    }
    if (curr_row == diag_index / 2) Efound = true;
  }
  return curr_row;
}

int WrightFisher::computeG(
    size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Compute the quantity G for different theta
{
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  int Kx = static_cast<int>(
      max(theta[N_i_s] / thetaP[N_i_s][1], (1.0 + theta[N_i_s]) / (1.0 - x)) +
      ((1.0 + (1.0 / x)) *
       max((((x * theta[N_i_s]) + 1.0) / thetaP[N_i_s][0]), 1.0)));
  int Kz = static_cast<int>(
      max(theta[N_i_t] / thetaP[N_i_t][1], (1.0 + theta[N_i_t]) / (1.0 - z)) +
      ((1.0 + (1.0 / z)) *
       max((((z * theta[N_i_t]) + 1.0) / thetaP[N_i_t][0]), 1.0)));
  int Dx = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(s) - (theta[N_i_s] + 1.0) / 2.0)));
  while ((theta[N_i_s] + 2 * Dx + 1) *
             exp(-(2 * Dx + theta[N_i_s]) * s / 2.0) >=
         1 - o.eps)
    ++Dx;
  int Dz = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(t - s) - (theta[N_i_t] + 1.0) / 2.0)));
  while ((theta[N_i_t] + 2 * Dz + 1) *
             exp(-(2 * Dz + theta[N_i_t]) * (t - s) / 2.0) >=
         1 - o.eps)
    ++Dz;
  int start_g = max(max(max(Kx, Kz), Dx), Dz);
  pair<vector<int>, double> Fs, Ft;
  Fs.second = s;
  Ft.second = t - s;
  int Fgs = 0, Fgt = 0;
  int g_ind = 0;
  while (g_ind <= static_cast<int>(max(max(Fgs, start_g), Fgt))) {
    g_ind++;
    Fgs = max(Fgs, computeF(N_i_s, g_ind, Fs));
    Fgt = max(Fgt, computeF(N_i_t, g_ind, Ft));
  }

  return g_ind;
}

int WrightFisher::computeGBoundaries(
    size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Compute the quantity G for different theta when
                       /// starting from boundary
{
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  int Kx = 1.0;
  int Kz = (z == 1.0)
               ? static_cast<int>(ceil((thetaP[N_i_s][0] * theta[N_i_t] /
                                        (theta[N_i_s] * thetaP[N_i_t][0]))))
               : static_cast<int>(ceil((thetaP[N_i_s][1] * theta[N_i_t] /
                                        (theta[N_i_s] * thetaP[N_i_t][1]))));
  int Dx = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(s) - (theta[N_i_s] + 1.0) / 2.0)));
  while ((theta[N_i_s] + 2 * Dx + 1) *
             exp(-(2 * Dx + theta[N_i_s]) * s / 2.0) >=
         1 - o.eps)
    ++Dx;
  int Dz = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(t - s) - (theta[N_i_t] + 1.0) / 2.0)));
  while ((theta[N_i_t] + 2 * Dz + 1) *
             exp(-(2 * Dz + theta[N_i_t]) * (t - s) / 2.0) >=
         1 - o.eps)
    ++Dz;
  int start_g = max(max(max(Kx, Kz), Dx), Dz);
  pair<vector<int>, double> Fs, Ft;
  Fs.second = s;
  Ft.second = t - s;
  int Fgs = 0, Fgt = 0;
  int g_ind = 0;
  while (g_ind <= static_cast<int>(max(max(Fgs, start_g), Fgt))) {
    g_ind++;
    Fgs = max(Fgs, computeF(N_i_s, g_ind, Fs));
    Fgt = max(Fgt, computeF(N_i_t, g_ind, Ft));
  }

  return g_ind;
}

int WrightFisher::computeGZInterior(
    size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Compute the quantity G for different theta when x is
                       /// on the boundary
{
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  int Kx = 1.0;
  int Kz = static_cast<int>(
      max(theta[N_i_t] / thetaP[N_i_t][1], (1.0 + theta[N_i_t]) / (1.0 - z)) +
      ((1.0 + (1.0 / z)) *
       max((((z * theta[N_i_t]) + 1.0) / thetaP[N_i_t][0]), 1.0)));
  int Dx = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(s) - (theta[N_i_s] + 1.0) / 2.0)));
  while ((theta[N_i_s] + 2 * Dx + 1) *
             exp(-(2 * Dx + theta[N_i_s]) * s / 2.0) >=
         1 - o.eps)
    ++Dx;
  int Dz = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(t - s) - (theta[N_i_t] + 1.0) / 2.0)));
  while ((theta[N_i_t] + 2 * Dz + 1) *
             exp(-(2 * Dz + theta[N_i_t]) * (t - s) / 2.0) >=
         1 - o.eps)
    ++Dz;
  int start_g = max(max(max(Kx, Kz), Dx), Dz);
  pair<vector<int>, double> Fs, Ft;
  Fs.second = s;
  Ft.second = t - s;
  int Fgs = 0, Fgt = 0;
  int g_ind = 0;
  while (g_ind <= static_cast<int>(max(max(Fgs, start_g), Fgt))) {
    g_ind++;
    Fgs = max(Fgs, computeF(N_i_s, g_ind, Fs));
    Fgt = max(Fgt, computeF(N_i_t, g_ind, Ft));
  }

  return g_ind;
}

int WrightFisher::computeGXInterior(
    size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Compute the quantity G for different theta when z is
                       /// on the boundary
{
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  int Kx = static_cast<int>(
      max(theta[N_i_s] / thetaP[N_i_s][1], (1.0 + theta[N_i_s]) / (1.0 - x)) +
      ((1.0 + (1.0 / x)) *
       max((((x * theta[N_i_s]) + 1.0) / thetaP[N_i_s][0]), 1.0)));
  int Kz = (z == 1.0) ? ((thetaP[N_i_s][0] / theta[N_i_s]) *
                         (theta[N_i_t] / thetaP[N_i_t][0]))
                      : ((thetaP[N_i_s][1] / theta[N_i_s]) *
                         (theta[N_i_t] / thetaP[N_i_t][1]));
  int Dx = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(s) - (theta[N_i_s] + 1.0) / 2.0)));
  while ((theta[N_i_s] + 2 * Dx + 1) *
             exp(-(2 * Dx + theta[N_i_s]) * s / 2.0) >=
         1 - o.eps)
    ++Dx;
  int Dz = static_cast<int>(ceil(
      max(0.0, 1.0 / static_cast<double>(t - s) - (theta[N_i_t] + 1.0) / 2.0)));
  while ((theta[N_i_t] + 2 * Dz + 1) *
             exp(-(2 * Dz + theta[N_i_t]) * (t - s) / 2.0) >=
         1 - o.eps)
    ++Dz;
  int start_g = max(max(max(Kx, Kz), Dx), Dz);
  pair<vector<int>, double> Fs, Ft;
  Fs.second = s;
  Ft.second = t - s;
  int Fgs = 0, Fgt = 0;
  int g_ind = 0;
  while (g_ind <= static_cast<int>(max(max(Fgs, start_g), Fgt))) {
    g_ind++;
    Fgs = max(Fgs, computeF(N_i_s, g_ind, Fs));
    Fgt = max(Fgt, computeF(N_i_t, g_ind, Ft));
  }

  return g_ind;
}

double WrightFisher::NormCDF(double x, double m,
                             double v)  /// CDF for N(m,v) - v is the variance
{
  return 0.5 * erfc(-(x - m) / sqrt(2 * v));
}

double WrightFisher::DiscretisedNormCDF(
    size_t N_i, int m,
    double t)  /// CDF for discretised normal, binning mass into nearest integer
{
  double returnval;
  int threshold = (thetaP[N_i].empty()                                ? 2
                   : (thetaP[N_i][0] == 0 || !(thetaP[N_i][1] > 0.0)) ? 1
                                                                      : 0);
  double beta = (theta[N_i] - 1.0) * static_cast<double>(t) / 2.0;
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

double WrightFisher::LogBinomialCoefficientCalculator(
    int n, int k)  /// Calculate usual binomial, with approximations kicking
                   /// in for large n and k
{
  assert(n >= 0 && k >= 0 && k <= n);
  if (n <= 1000) {
    return log(boost::math::binomial_coefficient<double>(n, k));
  } else {  // Compute log approximation using Stirling's formula
    return static_cast<double>(n) * log(static_cast<double>(n)) -
           static_cast<double>(k) * log(static_cast<double>(k)) -
           static_cast<double>(n - k) * log(static_cast<double>(n - k)) +
           0.5 * (log(static_cast<double>(n)) - log(static_cast<double>(k)) -
                  log(static_cast<double>(n - k)) -
                  log(2 * boost::math::constants::pi<double>()));
  }
}

double WrightFisher::UnconditionedDiffusionDensity(
    size_t N_i, double x, double y, double t,
    const Options &o)  /// Compute truncation to diffusion transition density
                       /// where diffusion can be absorbed at any time
{
  assert((x > 0.0) && (x < 1.0) && (y >= 0.0) && (y <= 1.0) && (t > 0.0));

  int thetaDependent =
      (thetaP[N_i].empty() ? 1
                           : 0);  /// Check what thetaP configuration we have
  double density = 0.0, density_inc, threshold = 1.0e-12;
  int mMode = static_cast<int>(floor(GriffithsParas(N_i, t).first)), m = mMode,
      mFlip = 1, mU = 0, mD = 0;
  bool mSwitch = false, mUpSwitch = false, mDownSwitch = false;

  while (!mSwitch) {
    if (thetaP[N_i].empty()) {
      double addon =
          (!(y > 0.0) ? pow(1.0 - x, static_cast<double>(m))
                      : (!(y < 1.0) ? pow(x, static_cast<double>(m)) : 0.0));
      for (int l = 1; l <= m - 1; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<> BETA(static_cast<double>(l),
                                              static_cast<double>(m - l));
        addon += pdf(BIN, l) * pdf(BETA, y);
      }
      density_inc = QmApprox(N_i, m, t, o) * addon;
    } else if (!(thetaP[N_i][0] > 0.0)) {
      double addon = (!(y > 0.0) ? pow(1.0 - x, static_cast<double>(m)) : 0.0);
      for (int l = 1; l <= m; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<> BETA(
            static_cast<double>(l), static_cast<double>(theta[N_i] + m - l));
        addon += pdf(BIN, l) * pdf(BETA, y);
      }
      density_inc = QmApprox(N_i, m, t, o) * addon;
    } else {
      double addon = (!(y < 1.0) ? pow(x, static_cast<double>(m)) : 0.0);
      for (int l = 0; l <= m - 1; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<> BETA(
            static_cast<double>(theta[N_i] + l), static_cast<double>(m - l));
        addon += pdf(BIN, l) * pdf(BETA, y);
      }
      density_inc = QmApprox(N_i, m, t, o) * addon;
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

double WrightFisher::DiffusionDensityApproximationDenom(
    size_t N_i, double x, double t,
    const Options
        &o)  /// Compute denominator for truncation to diffusion transition
             /// density when diffusion conditioned on non-absorption
{
  assert((x >= 0.0) && (x <= 1.0) && (t > 0.0));
  int thetaDependent =
      (thetaP[N_i].empty()
           ? 2
           : ((!(thetaP[N_i][0] > 0.0) || !(thetaP[N_i][1] > 0.0)) ? 1 : 0));
  double denom;

  if ((thetaDependent == 1) ||
      (thetaDependent == 2))  /// Check whether we have a zero mutation entry
  {
    double denom_inc = 1.0;
    denom = 0.0;
    int dMode = static_cast<int>(ceil(GriffithsParas(N_i, t).first)), d = dMode,
        Dflip = 1, Djm = 0, Djp = 0;

    while (denom_inc >
           0.0)  /// As long as increments are positive, keep computing
    {
      if (!(x > 0.0) || !(x < 1.0)) {
        denom_inc = QmApprox(N_i, d, t, o) * static_cast<double>(d);
      } else {
        denom_inc =
            QmApprox(N_i, d, t, o) * (1.0 - pow(x, static_cast<double>(d)) -
                                      pow(1.0 - x, static_cast<double>(d)));
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

double WrightFisher::DiffusionDensityApproximation(
    size_t N_i, double x, double y, double t,
    const Options &o)  /// Compute truncation to diffusion transition density
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (t > 0.0));
  int thetaDependent =
      (thetaP[N_i].empty()
           ? 2
           : ((thetaP[N_i][0] > 0.0 && thetaP[N_i][1] > 0.0)
                  ? 0
                  : 1));  /// Check what thetaP configuration we have
  double density = 0.0, density_inc,
         denom = DiffusionDensityApproximationDenom(N_i, x, t, o),
         threshold = max(1.0e-12 * denom, 1.0e-50);
  int mMode = static_cast<int>(floor(GriffithsParas(N_i, t).first)), m = mMode,
      mFlip = 1, mU = 0, mD = 0;
  bool mSwitch = false, mUpSwitch = false, mDownSwitch = false;

  while (!mSwitch) {
    if (thetaP[N_i].empty()) {
      if (!(x > 0.0)) {
        density_inc = QmApprox(N_i, m, t, o) *
                      static_cast<double>(m * (m - 1)) *
                      pow(1.0 - y, static_cast<double>(m - 2));
      } else if (!(x < 1.0)) {
        density_inc = QmApprox(N_i, m, t, o) *
                      static_cast<double>(m * (m - 1)) *
                      pow(y, static_cast<double>(m - 2));
      } else {
        double addon = 0.0;
        for (int l = 1; l != m; l++) {
          boost::math::binomial_distribution<> BIN(m, x);
          boost::math::beta_distribution<double> BETA(
              static_cast<double>(l), static_cast<double>(m - l));
          if (l < 1 && y == 0.0) {
            addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
          } else if (m - l < 1 && y == 1.0) {
            addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
          } else {
            addon += pdf(BIN, l) * pdf(BETA, y);
          }
        }
        density_inc = QmApprox(N_i, m, t, o) * addon;
      }
    } else if (!(thetaP[N_i][0] > 0.0)) {
      if (!(x > 0.0)) {
        density_inc = QmApprox(N_i, m, t, o) *
                      static_cast<double>(m * (theta[N_i] + m - 1)) *
                      pow(1.0 - y, static_cast<double>(theta[N_i] + m - 2));
      } else {
        double addon = 0.0;
        for (int l = 1; l != m + 1; l++) {
          boost::math::binomial_distribution<> BIN(m, x);
          boost::math::beta_distribution<double> BETA(
              static_cast<double>(l), static_cast<double>(theta[N_i] + m - l));
          if (l < 1 && y == 0.0) {
            addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
          } else if (theta[N_i] + m - l < 1.0 && y == 1.0) {
            addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
          } else {
            addon += pdf(BIN, l) * pdf(BETA, y);
          }
        }
        density_inc = QmApprox(N_i, m, t, o) * addon;
      }
    } else if (!(thetaP[N_i][1] > 0.0)) {
      if (!(x < 1.0)) {
        density_inc = QmApprox(N_i, m, t, o) *
                      static_cast<double>(m * (theta[N_i] + m - 1)) *
                      pow(y, static_cast<double>(theta[N_i] + m - 2));
      } else {
        double addon = 0.0;
        for (int l = 0; l != m; l++) {
          boost::math::binomial_distribution<> BIN(m, x);
          boost::math::beta_distribution<double> BETA(
              static_cast<double>(theta[N_i] + l), static_cast<double>(m - l));
          if (theta[N_i] + l < 1.0 && y == 0.0) {
            addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
          } else if (m - l < 1 && y == 1.0) {
            addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
          } else {
            addon += pdf(BIN, l) * pdf(BETA, y);
          }
        }
        density_inc = QmApprox(N_i, m, t, o) * addon;
      }
    } else {
      double addon = 0.0;
      for (int l = 0; l != m + 1; l++) {
        boost::math::binomial_distribution<> BIN(m, x);
        boost::math::beta_distribution<double> BETA(
            static_cast<double>(thetaP[N_i][0] + l),
            static_cast<double>(thetaP[N_i][1] + m - l));
        if (thetaP[N_i][0] + l < 1.0 && y == 0.0) {
          addon += pdf(BIN, l) * pdf(BETA, y + 1.0e-12);
        } else if (thetaP[N_i][1] + m - l < 1.0 && y == 1.0) {
          addon += pdf(BIN, l) * pdf(BETA, y - 1.0e-12);
        } else {
          addon += pdf(BIN, l) * pdf(BETA, y);
        }
      }
      density_inc = QmApprox(N_i, m, t, o) * addon;
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

double WrightFisher::QmApprox(
    size_t N_i, int m, double t,
    const Options &o)  /// Compute an approximation to q_m(t)
{
  assert((m >= 0) && (t > 0.0));
  double qm = 0.0, qmold = -1.0;

  if (t <= o.g1984threshold)  /// If time increment is too small, use
                              /// discretised Gaussian
  {
    qm = DiscretisedNormCDF(N_i, m, t);
  } else {
    int mkIndex = m;
    while (abs(qm - qmold) > 1.0e-12 || qm < 0.0 ||
           qm > 1.0)  /// If increments are big enough, keep going
    {
      qmold = qm;
      qm += pow(-1.0, mkIndex - m) *
            exp(Getlogakm<double>(N_i, mkIndex, m) +
                static_cast<double>(-(mkIndex) * (mkIndex + theta[N_i] - 1) *
                                    (t) / 2.0));
      mkIndex++;

      if (!(qm > qmold) && !(qm < qmold) &&
          (qm < 0.0 || qm > 1.0))  /// We have lost precision, so use
                                   /// discretised normal approximation
      {
        return DiscretisedNormCDF(N_i, m, t);
      }
    }
  }
  assert(qm >= 0.0);

  return qm;
}

double WrightFisher::UnconditionedBridgeDensity(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options
        &o)  /// Compute an approximation to the transition density of the
             /// diffusion bridge when absorption can happen at any time
{
  assert((x > 0.0) && (x < 1.0) && (y >= 0.0) && (y <= 1.0) && (s > 0.0) &&
         (s < t));

  if (thetaP[N_i].empty() &&
      ((!(z > 0.0) && !(y < 1.0)) || (!(z < 1.0) && !(y > 0.0)))) {
    return 0.0;
  }

  double eC = 0.0, denom_inc = 1.0;
  int dmode = static_cast<int>(ceil(GriffithsParas(N_i, t).first)), d = dmode,
      Dflip = 1, Djm = 0, Djp = 0, mkdLower = (thetaP[N_i].empty() ? 1 : 0);

  while (denom_inc >
         0.0)  /// As long as increments are positive, keep computing
  {
    if (!(z > 0.0))  /// z = 0
    {
      denom_inc = QmApprox(N_i, d, t, o) * pow(1.0 - x, static_cast<double>(d));
    } else  /// z = 1
    {
      denom_inc = QmApprox(N_i, d, t, o) * pow(x, static_cast<double>(d));
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

  if ((!(y > 0.0) && ((thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0)))) ||
      (!(y < 1.0) && ((thetaP[N_i].empty() || !(thetaP[N_i][1] > 0.0))))) {
    int mMode = GriffithsParas(N_i, s).first,
        kMode =
            GriffithsParas(N_i, t - s).first;  /// Use these together eC to get
    /// estimate of suitable threshold
    double constcontr =
        ((!(y > 0.0) && ((thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0))))
             ? static_cast<double>(mMode) * log(1.0 - x)
             : static_cast<double>(mMode) * log(x));
    double density = 0.0,
           threshold =
               max(exp(log(max(1.0e-300, QmApprox(N_i, mMode, s, o))) +
                       log(max(1.0e-300, QmApprox(N_i, kMode, t - s, o))) +
                       constcontr - log(eC)) *
                       1.0e-20,
                   1.0e-50);

    int m = mMode, mFlip = 1, mD = 0, mU = 0;
    bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;

    while (!mSwitch) {
      int k = kMode, kFlip = 1, kD = 0, kU = 0;
      bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;

      while (!kSwitch) {
        double num_inc;
        if (!(y > 0.0) && ((thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0)))) {
          num_inc =
              exp(log(max(1.0e-300, QmApprox(N_i, m, s, o))) +
                  log(max(1.0e-300, QmApprox(N_i, k, t - s, o))) +
                  static_cast<double>(m) * log(1.0 - x) -
                  log(eC));  /// Putting all the separate contributions together
        } else {
          num_inc =
              exp(log(max(1.0e-300, QmApprox(N_i, m, s, o))) +
                  log(max(1.0e-300, QmApprox(N_i, k, t - s, o))) +
                  static_cast<double>(m) * log(x) -
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
  } else if (thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0)) {
    vector<int> modeGuess =
        mklModeFinder(false, N_i, x, z, s, t,
                      o);  /// Compute mode over (m,k,l) to use as start points
    int mMode = modeGuess[0], kMode = modeGuess[1],
        lMode = modeGuess[2];  /// Use these estimates together with eC as a
    /// gauge for suitable threshold
    double ycontr, xcontr;
    if (((lMode == 0) && (!(z > 0.0))) || ((lMode == mMode) && (!(z < 1.0)))) {
      xcontr =
          static_cast<double>(mMode) * (!(z > 0.0) ? log(1.0 - x) : log(x));
      ycontr = 0.0;
    } else {
      boost::math::binomial_distribution<> BIN(mMode, x);
      xcontr = log(pdf(BIN, lMode));
      double p1 = static_cast<double>(lMode);
      double p2 = static_cast<double>(theta[N_i] + mMode - lMode);
      boost::math::beta_distribution<double> BETA(p1, p2);
      ycontr = log(pdf(BETA, y)) + static_cast<double>(kMode) *
                                       (!(z > 0.0) ? log(1.0 - y) : log(y));
    }
    double density = 0.0,
           threshold =
               max(exp(log(max(1.0e-300, QmApprox(N_i, mMode, s, o))) +
                       log(max(1.0e-300, QmApprox(N_i, kMode, t - s, o))) +
                       xcontr + ycontr - log(eC)) *
                       1.0e-20,
                   1.0e-50);

    int m = mMode, mFlip = 1, mD = 0, mU = 0;
    bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
    double mContr_D = boost::math::lgamma(static_cast<double>(theta[N_i] + m)),
           mContr_U = mContr_D, mContr;
    double constContr = -log(y) - log(1.0 - y);
    /// Allows to avoid calling beta functions, and thus faster
    while (!mSwitch) {
      double qm = max(1.0e-300, QmApprox(N_i, m, s, o));
      if (m != mMode)  /// Increment m contributions accordingly
      {
        if (mU > mD) {
          mContr_U += log(static_cast<double>(theta[N_i] + (m - 1)));
          mContr = log(qm) + mContr_U;
        } else {
          mContr_D -= log(static_cast<double>(theta[N_i] + (m + 1) - 1));
          mContr = log(qm) + mContr_D;
        }
      } else {
        mContr = log(qm) + mContr_U;
      }

      int k = kMode, kFlip = 1, kD = 0, kU = 0;
      bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
      double kContr_D =
                 static_cast<double>(k) * (!(z > 0.0) ? log(1.0 - y) : log(y)),
             kContr_U = kContr_D, kContr;  /// Calculate k contributions

      while (!kSwitch) {
        double qk = max(1.0e-300, QmApprox(N_i, k, t - s, o));
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

        int lLower = (thetaP[N_i].empty()
                          ? 1
                          : ((!(thetaP[N_i][0] > 0.0) && !(z > 0.0)) ? 1 : 0)),
            lUpper =
                (thetaP[N_i].empty()
                     ? m - 1
                     : ((!(thetaP[N_i][1] > 0.0) && !(z < 1.0)) ? m - 1 : m));
        int lFlip = 1, newlMode = min(lMode, lUpper), l = newlMode, lU = 0,
            lD = 0;  /// Need to ensure l <= m as m changes!
        bool lSwitch = false, lDownSwitch = false,
             lUpSwitch = false;  /// Compute l contributions
        boost::math::binomial_distribution<double> BIN(m, x);
        double lContr_D =
            log(pdf(BIN, l)) + static_cast<double>(l) * log(y) +
            static_cast<double>(m - l) * log(1.0 - y) -
            boost::math::lgamma(static_cast<double>(l)) -
            boost::math::lgamma(static_cast<double>(theta[N_i] + m - l));
        double lContr_U = lContr_D, lContr;

        while (!lSwitch) {
          if (l != newlMode) {
            if (lU > lD) {
              lContr_U +=
                  log(static_cast<double>(m - (l - 1))) -
                  log(static_cast<double>((l - 1) + 1)) + log(x) -
                  log(1.0 - x) + log(y) - log(1.0 - y) +
                  log(static_cast<double>(theta[N_i] + m - (l - 1) - 1)) -
                  log(static_cast<double>((l - 1)));
              lContr = lContr_U;
            } else {
              lContr_D += log(static_cast<double>(l + 1)) -
                          log(static_cast<double>(m - (l + 1) + 1)) +
                          log(1.0 - x) - log(x) + log(1.0 - y) - log(y) +
                          log(static_cast<double>((l + 1) - 1)) -
                          log(static_cast<double>(theta[N_i] + m - (l + 1)));
              lContr = lContr_D;
            }
          } else {
            lContr = lContr_U;
          }
          double density_inc =
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
    vector<int> modeGuess =
        mklModeFinder(false, N_i, x, z, s, t,
                      o);  /// Compute mode over (m,k,l) to use as start points
    int mMode = modeGuess[0], kMode = modeGuess[1],
        lMode = modeGuess[2];  /// Use these estimates together with eC as a
    /// gauge for suitable threshold
    double xcontr, ycontr;
    if (lMode == 0) {
      xcontr = static_cast<double>(mMode) * log(x);
      ycontr = 0.0;
    } else {
      boost::math::binomial_distribution<> BIN(mMode, x);
      xcontr = log(pdf(BIN, lMode));
      double p1 = static_cast<double>(theta[N_i] + lMode),
             p2 = static_cast<double>(mMode - lMode);
      boost::math::beta_distribution<double> B1(p1, p2);
      ycontr = log(pdf(B1, y)) + static_cast<double>(kMode) * log(y);
    }
    double density = 0.0,
           threshold =
               max(exp(log(max(1.0e-300, QmApprox(N_i, mMode, s, o))) +
                       log(max(1.0e-300, QmApprox(N_i, kMode, t - s, o))) +
                       xcontr + ycontr - log(eC)) *
                       1.0e-20,
                   1.0e-50);

    int m = mMode, mFlip = 1, mD = 0, mU = 0;
    bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
    double mContr_D = boost::math::lgamma(static_cast<double>(theta[N_i] + m)),
           mContr_U = mContr_D, mContr;
    double constContr = -log(y) - log(1.0 - y);
    /// Allows to avoid calling beta functions, and thus faster
    while (!mSwitch) {
      double qm = max(1.0e-300, QmApprox(N_i, m, s, o));
      if (m != mMode)  /// Increment m contributions accordingly
      {
        if (mU > mD) {
          mContr_U += log(static_cast<double>(theta[N_i] + (m - 1)));
          mContr = log(qm) + mContr_U;
        } else {
          mContr_D -= log(static_cast<double>(theta[N_i] + (m + 1) - 1));
          mContr = log(qm) + mContr_D;
        }
      } else {
        mContr = log(qm) + mContr_U;
      }

      int k = kMode, kFlip = 1, kD = 0, kU = 0;
      bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
      double kContr_D = static_cast<double>(k) * log(y), kContr_U = kContr_D,
             kContr;  /// Calculate k contributions

      while (!kSwitch) {
        double qk = max(1.0e-300, QmApprox(N_i, k, t - s, o));
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
        int lLower = (thetaP[N_i].empty()
                          ? 1
                          : ((!(thetaP[N_i][0] > 0.0) && !(z > 0.0)) ? 1 : 0)),
            lUpper =
                (thetaP[N_i].empty()
                     ? m - 1
                     : ((!(thetaP[N_i][1] > 0.0) && !(z < 1.0)) ? m - 1 : m));
        bool lSwitch = false, lDownSwitch = false,
             lUpSwitch = false;  /// Compute l contributions
        boost::math::binomial_distribution<double> BIN(m, x);
        double lContr_D =
            log(pdf(BIN, l)) + static_cast<double>(l) * log(y) +
            static_cast<double>(m - l) * log(1.0 - y) -
            boost::math::lgamma(static_cast<double>(theta[N_i] + l)) -
            boost::math::lgamma(static_cast<double>(m - l));
        double lContr_U = lContr_D, lContr;

        while (!lSwitch) {
          if (l != newlMode) {
            if (lU > lD) {
              lContr_U += log(static_cast<double>(m - (l - 1))) -
                          log(static_cast<double>((l - 1) + 1)) + log(x) -
                          log(1.0 - x) + log(y) - log(1.0 - y) +
                          log(static_cast<double>(m - (l - 1) - 1)) -
                          log(static_cast<double>(theta[N_i] + (l - 1)));
              lContr = lContr_U;
            } else {
              lContr_D += log(static_cast<double>(l + 1)) -
                          log(static_cast<double>(m - (l + 1) + 1)) +
                          log(1.0 - x) - log(x) + log(1.0 - y) - log(y) +
                          log(static_cast<double>(theta[N_i] + (l + 1) - 1)) -
                          log(static_cast<double>(m - (l + 1)));
              lContr = lContr_D;
            }
          } else {
            lContr = lContr_U;
          }
          double density_inc =
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

double WrightFisher::BridgeDenom(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options &o)  /// Compute an approximation to the denominator for the
                       /// transition density of the diffusion bridge
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  double denom = 0.0;
  double inc = 1.0;
  int dmode = static_cast<int>(ceil(GriffithsParas(N_i, t).first)), d = dmode,
      Dflip = 1, Djm = 0, Djp = 0;

  while (inc > 0.0)  /// As long as increments are positive, keep computing
  {
    if ((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0))  /// x,z in (0,1)
    {
      double betabin = 0.0;
      for (int f = 0; f != d + 1; f++) {
        boost::math::binomial_distribution<double> BIN(d, x);
        double para1 = static_cast<double>(thetaP[N_i][0] + f),
               para2 = static_cast<double>(thetaP[N_i][1] + d - f);
        boost::math::beta_distribution<double> BETA(para1, para2);

        betabin += pdf(BIN, f) * pdf(BETA, z);
      }

      inc = QmApprox(N_i, d, t, o) * betabin;
    } else if ((z > 0.0) && (z < 1.0) &&
               (!(x > 0.0) || !(x < 1.0)))  /// x in {0,1}, z in (0,1)
    {
      double para1 = (!(x > 0.0) ? thetaP[N_i][0]
                                 : static_cast<double>(thetaP[N_i][0] + d)),
             para2 = (!(x > 0.0) ? static_cast<double>(thetaP[N_i][1] + d)
                                 : thetaP[N_i][1]);
      boost::math::beta_distribution<double> BETAZ(para1, para2);
      inc = QmApprox(N_i, d, t, o) * pdf(BETAZ, z);
    } else if ((x > 0.0) && (x < 1.0) &&
               (!(z > 0.0) || !(z < 1.0)))  /// x in (0,1), z in {0,1}
    {
      double para1 = (!(z > 0.0) ? thetaP[N_i][0]
                                 : static_cast<double>(thetaP[N_i][0] + d)),
             para2 = (!(z > 0.0) ? static_cast<double>(thetaP[N_i][1] + d)
                                 : thetaP[N_i][1]);
      double xcon = (!(z > 0.0) ? pow(1.0 - x, static_cast<double>(d))
                                : pow(x, static_cast<double>(d)));
      inc = QmApprox(N_i, d, t, o) * xcon *
            (1.0 / boost::math::beta<double>(para1, para2));
    } else if (!(x < z) && !(x > z))  /// x,z in {0,1} and x=z
    {
      double para1 = (!(x > 0.0) ? thetaP[N_i][0]
                                 : static_cast<double>(thetaP[N_i][0] + d)),
             para2 = (!(x > 0.0) ? static_cast<double>(thetaP[N_i][1] + d)
                                 : thetaP[N_i][1]);
      inc = QmApprox(N_i, d, t, o) / boost::math::beta<double>(para1, para2);
    } else  /// x,z in {0,1} and x!=z
    {
      return QmApprox(N_i, 0, t, o) /
             boost::math::beta<double>(thetaP[N_i][0], thetaP[N_i][1]);
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

double WrightFisher::ComputeDensity1(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options &o)  /// Compute bridge density when x,z in (0,1) (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }
  double logx = log(x), log1x = log(1.0 - x), logy = log(y),
         log1y = log(1.0 - y), logz = log(z), log1z = log(1.0 - z);
  vector<int> modeGuess =
      mkljModeFinder(false, N_i, x, z, s, t,
                     o);  /// Compute mode over (m,k,l,j) to use as start points
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2],
      jMode = modeGuess[3];  /// Use these estimates together with eC as a gauge
  /// for suitable threshold
  double qmmode = max(1.0e-300, QmApprox(N_i, mMode, s, o)),
         qkmode = max(1.0e-300, QmApprox(N_i, kMode, t - s, o));
  boost::math::binomial_distribution<> BINx(mMode, x), BINy(kMode, y);
  boost::math::beta_distribution<double> BETAz(
      static_cast<double>(thetaP[N_i][0] + jMode),
      static_cast<double>(thetaP[N_i][1] + kMode - jMode)),
      BETAy(static_cast<double>(thetaP[N_i][0] + lMode),
            static_cast<double>(thetaP[N_i][1] + mMode - lMode));
  double density = 0.0, eC = BridgeDenom(N_i, x, z, y, s, t, o),
         threshold = max(exp(log(qmmode) + log(qkmode) + log(pdf(BINx, lMode)) +
                             log(pdf(BINy, jMode)) + log(pdf(BETAz, z)) +
                             log(pdf(BETAy, y)) - log(eC)) *
                             1.0e-6,
                         1.0e-50);
  double constContr =
      static_cast<double>(thetaP[N_i][0] - 1.0) * (logz + logy) +
      static_cast<double>(thetaP[N_i][1] - 1.0) * (log1z + log1y);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double mContr;  /// Compute contributions depending only on m
  /// Allows to avoid calling beta functions, and thus faster
  while (!mSwitch) {
    double qm = QmApprox(N_i, m, s, o);
    precomputeA(N_i, N_i, m, kMode);
    mContr = factorials[m] + static_cast<double>(m) * (log1x + log1y) +
             lg_theta[N_i][m];
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double kContr;  /// Calculate k contributions

    while (!kSwitch) {
      double qk = QmApprox(N_i, k, t - s, o);
      precomputeA(N_i, N_i, m, k);
      kContr = factorials[k] + static_cast<double>(k) * (log1z + log1y) +
               lg_theta[N_i][k];
      kContr += log(qk);

      int lFlip = 1, lU = 0, lD = 0, newlMode = min(lMode, m),
          l = newlMode;  /// Need to ensure l <= m since m changes!
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;
      double lContr;

      while (!lSwitch) {
        assert((l >= 0) && (l <= m));
        lContr = -factorials[l] - factorials[m - l] +
                 static_cast<double>(l) * (logx + logy - log1x - log1y) -
                 lg_theta1[N_i][l] - lg_theta2[N_i][m - l];

        int jFlip = 1, jU = 0, jD = 0, newjMode = min(jMode, k),
            j = newjMode;  /// Need to ensure j <= k as k changes!
        bool jSwitch = false, jDownSwitch = false,
             jUpSwitch = false;  /// Compute j contributions

        double jContr;

        while (!jSwitch) {
          jContr = -factorials[j] - factorials[k - j] +
                   static_cast<double>(j) * (logz + logy - log1z - log1y) -
                   lg_theta1[N_i][j] - lg_theta2[N_i][k - j];

          double density_inc =
              exp(mContr + kContr + lContr + jContr + constContr -
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

          if (!jSwitch)  /// If we cannot, we need to move upwards or
                         /// downwards
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

        if (!(lDownSwitch))  /// Same procedure as for j contributions, but
                             /// now we don't need to worry about increments
                             /// being below threshold
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

double WrightFisher::ComputeDensity2(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options
        &o)  /// Compute bridge density when x in {0,1},z in (0,1) (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }
  double logy = log(y), log1y = log(1.0 - y), logz = log(z),
         log1z = log(1.0 - z);
  vector<int> modeGuess =
      mkjModeFinder(false, N_i, x, z, s, t,
                    o);  /// Compute mode over (m,k,j) to use as start points
  int mMode = modeGuess[0], kMode = modeGuess[1],
      jMode = modeGuess[2];  /// Use these estimates together with eC as a gauge
  /// for suitable threshold
  double qmmode = max(1.0e-300, QmApprox(N_i, mMode, s, o)),
         qkmode = max(1.0e-300, QmApprox(N_i, kMode, t - s, o));
  double p1 = static_cast<double>(!(x > 0.0) ? thetaP[N_i][0]
                                             : thetaP[N_i][0] + mMode),
         p2 = static_cast<double>(!(x < 1.0) ? thetaP[N_i][1] + mMode
                                             : thetaP[N_i][1]);
  boost::math::binomial_distribution<> BINy(kMode, y);
  boost::math::beta_distribution<double> BETAz(
      static_cast<double>(thetaP[N_i][0] + jMode),
      static_cast<double>(thetaP[N_i][1] + kMode - jMode)),
      BETAy(p1, p2);
  double density = 0.0, eC = BridgeDenom(N_i, x, z, y, s, t, o),
         threshold =
             max(exp(log(qmmode) + log(qkmode) + log(pdf(BINy, jMode)) +
                     log(pdf(BETAz, z)) + log(pdf(BETAy, y)) - log(eC)) *
                     1.0e-6,
                 1.0e-50);
  precomputeA(N_i, N_i, mMode, kMode);
  double constContr =
      static_cast<double>(thetaP[N_i][0] - 1.0) * (logz + logy) +
      static_cast<double>(thetaP[N_i][1] - 1.0) * (log1z + log1y) +
      ((!(x > 0.0)) ? -lg_theta1[N_i][0] : -lg_theta2[N_i][0]);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double mContr;
  /// Allows to avoid calling beta functions, and thus faster
  while (!mSwitch) {
    precomputeA(N_i, N_i, m, kMode);
    double qm = QmApprox(N_i, m, s, o);
    mContr = lg_theta[N_i][m] +
             ((!(x > 0.0)) ? static_cast<double>(m) * log1y - lg_theta2[N_i][m]
                           : static_cast<double>(m) * logy - lg_theta1[N_i][m]);
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double kContr;  /// Calculate k contributions

    while (!kSwitch) {
      precomputeA(N_i, N_i, m, k);
      double qk = QmApprox(N_i, k, t - s, o);
      kContr = factorials[k] + lg_theta[N_i][k] +
               static_cast<double>(k) * (log1z + log1y);
      kContr += log(qk);

      int jFlip = 1, newjMode = min(jMode, k), j = newjMode, jU = 0,
          jD = 0;  /// Need to ensure j <= k as k changes!
      bool jSwitch = false, jDownSwitch = false,
           jUpSwitch = false;  /// Compute j contributions

      double jContr;
      while (!jSwitch) {
        jContr = -factorials[j] - factorials[k - j] +
                 static_cast<double>(j) * (logz + logy - log1z - log1y) -
                 lg_theta1[N_i][j] - lg_theta2[N_i][k - j];
        double density_inc =
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

double WrightFisher::ComputeDensity3(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options
        &o)  /// Compute bridge density when x,z in {0,1}, x=z (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }
  double logy = log(y), log1y = log(1.0 - y);
  vector<int> modeGuess =
      mkModeFinder(false, N_i, x, z, s, t, o);  /// Find mode over (m,k)
  int mMode = modeGuess[0],
      kMode = modeGuess[1];  /// Use these together eC to get estimate of
  /// suitable threshold
  double qmmode = max(1.0e-300, QmApprox(N_i, mMode, s, o)),
         qkmode = max(1.0e-300, QmApprox(N_i, kMode, t - s, o));
  double th1 =
             static_cast<double>(!(x > 0.0) ? thetaP[N_i][0] : thetaP[N_i][1]),
         th2 = theta[N_i] - th1;
  double gammaratios =
      -boost::math::lgamma(th1) +
      boost::math::lgamma(static_cast<double>(theta[N_i] + mMode)) -
      boost::math::lgamma(static_cast<double>(th2 + mMode)) +
      boost::math::lgamma(static_cast<double>(theta[N_i] + kMode)) -
      boost::math::lgamma(static_cast<double>(th2 + kMode)) +
      boost::math::lgamma(static_cast<double>(th2 + mMode + kMode)) -
      boost::math::lgamma(static_cast<double>(theta[N_i] + mMode + kMode));
  boost::math::beta_distribution<double> BETA(
      static_cast<double>(!(x > 0.0) ? thetaP[N_i][0]
                                     : thetaP[N_i][0] + mMode + kMode),
      static_cast<double>(!(x > 0.0) ? thetaP[N_i][1] + mMode + kMode
                                     : thetaP[N_i][1]));
  double density = 0.0, eC = BridgeDenom(N_i, x, z, y, s, t, o),
         threshold = max(exp(log(qmmode) + log(qkmode) + gammaratios +
                             log(pdf(BETA, y)) - log(eC)) *
                             1.0e-6,
                         1.0e-50);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double constContr =
      static_cast<double>(thetaP[N_i][0] - 1.0) * logy +
      static_cast<double>(thetaP[N_i][1] - 1.0) * log1y +
      (!(x > 0.0) ? -2.0 * lg_theta1[N_i][0] : -2.0 * lg_theta2[N_i][0]);
  double mContr;  /// Contributions from m

  while (!mSwitch) {
    precomputeA(N_i, N_i, m, kMode);
    double qm = QmApprox(N_i, m, s, o);
    mContr = lg_theta[N_i][m] +
             ((!(x > 0.0)) ? static_cast<double>(m) * log1y - lg_theta2[N_i][m]
                           : static_cast<double>(m) * logy - lg_theta1[N_i][m]);
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double kContr;  /// Contributions from k

    while (!kSwitch) {
      precomputeA(N_i, N_i, m, k);
      double qk = QmApprox(N_i, k, t - s, o);
      kContr =
          lg_theta[N_i][k] +
          ((!(x > 0.0)) ? static_cast<double>(k) * log1y - lg_theta2[N_i][k]
                        : static_cast<double>(k) * logy - lg_theta1[N_i][k]);
      kContr += log(qk);

      double density_inc =
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

double WrightFisher::ComputeDensity4(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options
        &o)  /// Compute bridge density when x,z in {0,1}, x!=z (any theta!)
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }
  double logy = log(y), log1y = log(1.0 - y);
  vector<int> modeGuess =
      mkModeFinder(false, N_i, x, z, s, t, o);  /// Find mode over (m,k)
  int mMode = modeGuess[0],
      kMode = modeGuess[1];  /// Use together with eC to get suitable threshold
  precomputeA(N_i, N_i, mMode, kMode);
  double qmmode = max(1.0e-300, QmApprox(N_i, mMode, s, o)),
         qkmode = max(1.0e-300, QmApprox(N_i, kMode, t - s, o));
  double gammaratios =
      -lg_theta1[N_i][0] - lg_theta2[N_i][0] + lg_theta[N_i][mMode] +
      lg_theta[N_i][kMode] -
      (!(x > 0.0) ? lg_theta1[N_i][kMode] : lg_theta1[N_i][mMode]) -
      (!(x > 0.0) ? lg_theta2[N_i][mMode] : lg_theta2[N_i][kMode]);
  double density = 0.0,
         eC = (t <= o.g1984threshold ? DiscretisedNormCDF(N_i, 0, t)
                                     : QmApprox(N_i, 0, t, o)) /
              boost::math::beta<double>(static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])),
         threshold = 1.0e-100;
  //  max(exp(log(qmmode) + log(qkmode) + gammaratios +
  //          static_cast<double>(thetaP[N_i][0] +
  //                              (!(x > 0.0) ? kMode : mMode) - 1.0) *
  //              logy +
  //          static_cast<double>(thetaP[N_i][1] +
  //                              (!(x > 0.0) ? mMode : kMode) - 1.0) *
  //              log1y -
  //          log(eC)) *
  //          1.0e-6,
  //      1.0e-50);

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double constContr = static_cast<double>(thetaP[N_i][0] - 1.0) * logy +
                      static_cast<double>(thetaP[N_i][1] - 1.0) * log1y -
                      lg_theta1[N_i][0] - lg_theta2[N_i][0];
  double mContr;  /// Contributions from m

  while (!mSwitch) {
    precomputeA(N_i, N_i, m, kMode);
    double qm = QmApprox(N_i, m, s, o);
    mContr = lg_theta[N_i][m] +
             ((!(x > 0.0)) ? static_cast<double>(m) * log1y - lg_theta2[N_i][m]
                           : static_cast<double>(m) * logy - lg_theta1[N_i][m]);
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double kContr;  /// Contributions from k

    while (!kSwitch) {
      precomputeA(N_i, N_i, m, k);
      double qk = QmApprox(N_i, k, t - s, o);
      kContr =
          lg_theta[N_i][k] +
          ((!(x > 0.0)) ? static_cast<double>(k) * logy - lg_theta1[N_i][k]
                        : static_cast<double>(k) * log1y - lg_theta2[N_i][k]);
      kContr += log(qk);

      double density_inc =
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

double WrightFisher::BridgeDensity(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options &o)  /// Function to determine which case from the above
                       /// computeDensity to invoke
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));

  if ((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0))  /// x,z in (0,1)
  {
    return ComputeDensity1(N_i, x, z, y, s, t, o);
  } else if ((x > 0.0) && (x < 1.0))  /// x in (0,1), z in {0,1}
  {
    return ComputeDensity2(N_i, z, x, y, t - s, t,
                           o);  /// Equivalent to a reverse bridge going from z
                                /// to y in time t-s and ending at x at time t
  } else if ((z > 0.0) && (z < 1.0))  /// x in {0,1}, z in (0,1)
  {
    return ComputeDensity2(N_i, x, z, y, s, t, o);
  } else if (x == z)  /// x,z in {0,1}, x=z
  {
    return ComputeDensity3(N_i, x, z, y, s, t, o);
  } else  /// x,z in {0,1}, x!=z
  {
    return ComputeDensity4(N_i, x, z, y, s, t, o);
  }
}

double WrightFisher::BridgeDiffThetaDensity(size_t N_i, double x, double z,
                                            double y, double s, double t,
                                            const Options &o) {
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));

  if ((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0))  /// x,z in (0,1)
  {
    return ComputeDensityDiffTheta(N_i, x, z, y, s, t, o);
  } else if (((x > 0.0) && (x < 1.0)) ||
             ((z > 0.0) &&
              (z <
               1.0)))  /// Either x or z (but at most one) are on the boundary!
  {
    return ComputeDensityDiffThetaInterior(N_i, x, z, y, s, t, o);
  } else {  // Both x and z are on the boundary
    return ComputeDensityDiffThetaBoundaries(N_i, x, z, y, s, t, o);
  }
}

double WrightFisher::ComputeDensityDiffTheta(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options &o)  /// Compute bridge density for different thetas
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logy = log(y),
         log1y = log(1.0 - y), logz = log(z), log1z = log(1.0 - z);
  vector<int> modeGuess =
      mkljModeFinder(true, N_i, x, z, s, t,
                     o);  /// Compute mode over (m,k,l,j) to use as start points
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2],
      jMode = modeGuess[3];  /// Use these estimates together with eC as a gauge
  /// for suitable threshold
  double qmmode = max(1.0e-300, QmApprox(N_i_s, mMode, s, o)),
         qkmode = max(1.0e-300, QmApprox(N_i_t, kMode, t - s, o));
  boost::math::binomial_distribution<> BINx(mMode, x), BINy(kMode, y);
  boost::math::beta_distribution<double> BETAz(
      static_cast<double>(thetaP[N_i_t][0] + jMode),
      static_cast<double>(thetaP[N_i_t][1] + kMode - jMode)),
      BETAy(static_cast<double>(thetaP[N_i_s][0] + lMode),
            static_cast<double>(thetaP[N_i_s][1] + mMode - lMode));
  double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-100;
  int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
  int dmmode = static_cast<int>(ceil(GriffithsParas(N_i_s, s).first)),
      dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
  bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
  int dkmode = static_cast<int>(ceil(GriffithsParas(N_i_t, t - s).first));
  while (!dmSwitch) {
    emCOldInc = emCInc;
    double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
    int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
    bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
    while (!dkSwitch) {
      ekCOldInc = ekCInc;
      ekCInc = QmApprox(N_i_t, dk, t - s, o) *
               calculate_expectation(N_i, dm, dk, x, z);
      ekC += ekCInc;
      if (!(dkDownSwitch))  /// Switching mechanism for j
      {
        if (sgn(dk - dkmode) <= 0) {
          dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
        }
      }

      if (!(dkUpSwitch)) {
        if (sgn(dk - dkmode) >= 0) {
          dkUpSwitch = (ekCInc < eC_threshold);
        }
      }

      dkSwitch = (dkDownSwitch && dkUpSwitch);
      if (!dkSwitch) {
        if (dkFlip == 1) {
          dkU++;
          dk = dkmode + dkU;
          dkFlip *= (dkDownSwitch ? 1 : -1);
        } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
          dkD++;
          dk = dkmode - dkD;
          dkFlip *= (dkUpSwitch ? 1 : -1);
        }
      }
    }
    emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
    eC += emCInc;

    if (!(dmDownSwitch))  /// Switching mechanism for j
    {
      if (sgn(dm - dmmode) <= 0) {
        dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
      }
    }

    if (!(dmUpSwitch)) {
      if (sgn(dm - dmmode) >= 0) {
        dmUpSwitch = (emCInc < eC_threshold);
      }
    }

    dmSwitch = (dmDownSwitch && dmUpSwitch);
    if (!dmSwitch) {
      if (dmFlip == 1) {
        dmU++;
        dm = dmmode + dmU;
        dmFlip *= (dmDownSwitch ? 1 : -1);
      } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
        dmD++;
        dm = dmmode - dmD;
        dmFlip *= (dmUpSwitch ? 1 : -1);
      }
    }
  }
  double density = 0.0,
         threshold = max(exp(log(qmmode) + log(qkmode) + log(pdf(BINx, lMode)) +
                             log(pdf(BINy, jMode)) + log(pdf(BETAz, z)) +
                             log(pdf(BETAy, y)) - log(eC)) *
                             1.0e-8,
                         1.0e-100);

  double constContr =
      (thetaP[N_i_t][0] - 1.0) * logz + (thetaP[N_i_t][1] - 1.0) * log1z +
      (thetaP[N_i_s][0] - 1.0) * logy + (thetaP[N_i_s][1] - 1.0) * log1y;
  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  precomputeA(N_i_s, N_i_t, m, kMode);
  double mContr;  /// Compute contributions depending only on m
  /// Allows to avoid calling beta functions, and thus faster
  while (!mSwitch) {
    double qm = QmApprox(N_i_s, m, s, o);
    mContr = lg_theta[N_i_s][m] + factorials[m] +
             static_cast<double>(m) * (log1x + log1y);
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double kContr;  /// Calculate k contributions

    while (!kSwitch) {
      precomputeA(N_i_s, N_i_t, m, k);
      double qk = QmApprox(N_i_t, k, t - s, o);
      kContr = lg_theta[N_i_t][k] + factorials[k] +
               static_cast<double>(k) * (log1z + log1y);
      kContr += log(qk);

      int lFlip = 1, lU = 0, lD = 0, newlMode = min(lMode, m),
          l = newlMode;  /// Need to ensure l <= m since m changes!
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;
      boost::math::binomial_distribution<double> BINL(
          m, x);  /// Contributions from l
      double lContr;

      while (!lSwitch) {
        assert((l >= 0) && (l <= m));
        lContr = -factorials[l] - factorials[m - l] +
                 static_cast<double>(l) * (logx + logy - log1x - log1y) -
                 lg_theta1[N_i_s][l] - lg_theta2[N_i_s][m - l];

        int jFlip = 1, jU = 0, jD = 0, newjMode = min(jMode, k),
            j = newjMode;  /// Need to ensure j <= k as k changes!
        bool jSwitch = false, jDownSwitch = false,
             jUpSwitch = false;  /// Compute j contributions
        double jContr;

        while (!jSwitch) {
          jContr = -factorials[j] - factorials[k - j] +
                   static_cast<double>(j) * (logz + logy - log1z - log1y) -
                   lg_theta1[N_i_t][j] - lg_theta2[N_i_t][k - j];

          double density_inc =
              exp(mContr + kContr + lContr + jContr + constContr -
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

          if (!jSwitch)  /// If we cannot, we need to move upwards or
                         /// downwards
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

        if (!(lDownSwitch))  /// Same procedure as for j contributions, but
                             /// now we don't need to worry about increments
                             /// being below threshold
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

double WrightFisher::ComputeDensityDiffThetaBoundaries(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options &o)  /// Compute bridge density for different theta when x, z
                       /// on the boundary
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logy = log(y), log1y = log(1.0 - y);
  vector<int> modeGuess =
      mkModeFinder(true, N_i, x, z, s, t,
                   o);  /// Compute mode over (m,k,l,j) to use as start points

  int mMode = modeGuess[0],
      kMode = modeGuess[1];  /// Use these estimates together with eC as a gauge
  /// for suitable threshold
  double qmmode = max(1.0e-300, QmApprox(N_i_s, mMode, s, o)),
         qkmode = max(1.0e-300, QmApprox(N_i_t, kMode, t - s, o));
  boost::math::binomial_distribution<> BINx(mMode, x), BINy(kMode, y);
  double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-100;
  int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
  int dmmode = static_cast<int>(ceil(GriffithsParas(N_i_s, s).first)),
      dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
  bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
  int dkmode = static_cast<int>(ceil(GriffithsParas(N_i_t, t - s).first));

  while (!dmSwitch) {
    emCOldInc = emCInc;
    double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
    int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
    bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
    while (!dkSwitch) {
      ekCOldInc = ekCInc;
      ekCInc = QmApprox(N_i_t, dk, t - s, o) *
               calculate_expectation(N_i, dm, dk, x, z);
      ekC += ekCInc;
      if (!(dkDownSwitch))  /// Switching mechanism for j
      {
        if (sgn(dk - dkmode) <= 0) {
          dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
        }
      }

      if (!(dkUpSwitch)) {
        if (sgn(dk - dkmode) >= 0) {
          dkUpSwitch = (ekCInc < eC_threshold);
        }
      }

      dkSwitch = (dkDownSwitch && dkUpSwitch);
      if (!dkSwitch) {
        if (dkFlip == 1) {
          dkU++;
          dk = dkmode + dkU;
          dkFlip *= (dkDownSwitch ? 1 : -1);
        } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
          dkD++;
          dk = dkmode - dkD;
          dkFlip *= (dkUpSwitch ? 1 : -1);
        }
      }
    }
    emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
    eC += emCInc;

    if (!(dmDownSwitch))  /// Switching mechanism for j
    {
      if (sgn(dm - dmmode) <= 0) {
        dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
      }
    }

    if (!(dmUpSwitch)) {
      if (sgn(dm - dmmode) >= 0) {
        dmUpSwitch = (emCInc < eC_threshold);
      }
    }

    dmSwitch = (dmDownSwitch && dmUpSwitch);
    if (!dmSwitch) {
      if (dmFlip == 1) {
        dmU++;
        dm = dmmode + dmU;
        dmFlip *= (dmDownSwitch ? 1 : -1);
      } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
        dmD++;
        dm = dmmode - dmD;
        dmFlip *= (dmUpSwitch ? 1 : -1);
      }
    }
  }
  precomputeA(N_i_s, N_i_t, mMode, kMode);
  double innerterm;
  if (!(x > 0.0)) {
    if (!(z > 0.0)) {  // x = 0 && z = 0
      innerterm =
          -lg_theta1[N_i_s][0] - lg_theta2[N_i_s][mMode] +
          lg_theta[N_i_s][mMode] - lg_theta1[N_i_t][0] -
          lg_theta2[N_i_t][kMode] + lg_theta[N_i_t][kMode] +
          static_cast<double>(thetaP[N_i_s][0] - 1.0) * logy +
          static_cast<double>(thetaP[N_i_s][1] + mMode + kMode - 1.0) * log1y;
    } else {  // x = 0 && z = 1
      innerterm = -lg_theta1[N_i_s][0] - lg_theta2[N_i_s][mMode] +
                  lg_theta[N_i_s][mMode] - lg_theta1[N_i_t][kMode] -
                  lg_theta2[N_i_t][0] + lg_theta[N_i_t][kMode] +
                  static_cast<double>(thetaP[N_i_s][0] + kMode - 1.0) * logy +
                  static_cast<double>(thetaP[N_i_s][1] + mMode - 1.0) * log1y;
    }
  } else {
    if (!(z > 0.0)) {  // x = 1 && z = 0
      innerterm = -lg_theta1[N_i_s][mMode] - lg_theta2[N_i_s][0] +
                  lg_theta[N_i_s][mMode] - lg_theta1[N_i_t][0] -
                  lg_theta2[N_i_t][kMode] + lg_theta[N_i_t][kMode] +
                  static_cast<double>(thetaP[N_i_s][0] + mMode - 1.0) * logy +
                  static_cast<double>(thetaP[N_i_s][1] + mMode - 1.0) * log1y;
    } else {  // x = 1 && z = 1
      innerterm =
          -lg_theta1[N_i_s][mMode] - lg_theta2[N_i_s][0] +
          lg_theta[N_i_s][mMode] - lg_theta1[N_i_t][kMode] -
          lg_theta2[N_i_t][0] + lg_theta[N_i_t][kMode] +
          static_cast<double>(thetaP[N_i_s][0] + mMode + kMode - 1.0) * logy +
          static_cast<double>(thetaP[N_i_s][1] - 1.0) * log1y;
    }
  }
  double density = 0.0,
         threshold =
             max(exp(log(qmmode) + log(qkmode) + innerterm - log(eC)) * 1.0e-8,
                 1.0e-50);
  double constant =
      (thetaP[N_i_s][0] - 1.0) * logy + (thetaP[N_i_s][1] - 1.0) * log1y +
      ((x == z) ? ((!(x > 0.0)) ? -lg_theta1[N_i_s][0] - lg_theta1[N_i_t][0]
                                : -lg_theta2[N_i_s][0] - lg_theta2[N_i_t][0])
                : ((!(x > 0.0)) ? -lg_theta1[N_i_s][0] - lg_theta2[N_i_t][0]
                                : -lg_theta2[N_i_s][0] - lg_theta1[N_i_t][0]));
  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double mContr, kContr;  /// Compute contributions depending only on m
  /// Allows to avoid calling beta functions, and thus faster
  while (!mSwitch) {
    precomputeA(N_i_s, N_i_t, m, kMode);
    double qm = QmApprox(N_i_s, m, s, o);
    mContr =
        lg_theta[N_i_s][m] +
        ((!(x > 0.0)) ? static_cast<double>(m) * log1y - lg_theta2[N_i_s][m]
                      : static_cast<double>(m) * logy - lg_theta1[N_i_s][m]);
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;

    while (!kSwitch) {
      precomputeA(N_i_s, N_i_t, m, k);
      double qk = QmApprox(N_i_t, k, t - s, o);
      kContr =
          lg_theta[N_i_t][k] +
          ((!(z > 0.0)) ? static_cast<double>(k) * log1y - lg_theta2[N_i_t][k]
                        : static_cast<double>(k) * logy -
                              lg_theta1[N_i_t][k]);  /// Calculate k
                                                     /// contributions
      kContr += log(qk);

      double density_inc =
          exp(mContr + kContr + constant -
              log(eC));  /// Put all separate contributions together
      density += density_inc;

      if (!(kDownSwitch))  /// Same switching procedure as for l
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

double WrightFisher::ComputeDensityDiffThetaInterior(
    size_t N_i, double x, double z, double y, double s, double t,
    const Options &o)  /// Compute bridge density for different theta when x, z
                       /// are in the interior
{
  assert((x >= 0.0) && (x <= 1.0) && (y >= 0.0) && (y <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (t > 0.0) && (s > 0.0) && (s < t));
  if (!(y > 0.0) || !(y < 1.0)) {
    return 0.0;
  }
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z), logy = log(y), log1y = log(1.0 - y);
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  vector<int> modeGuess =
      ((!(x > 0.0)) || (!(x < 1.0)))
          ? mkjModeFinder(true, N_i, x, z, s, t, o)
          : mklModeFinder(
                true, N_i, x, z, s, t,
                o);  /// Compute mode over (m,k,l,j) to use as start points
  int mMode = modeGuess[0], kMode = modeGuess[1],
      indexMode = ((!(x > 0.0)) || (!(x < 1.0)))
                      ? modeGuess[3]
                      : modeGuess[2];  /// Use these estimates together with eC
                                       /// as a gauge
  /// for suitable threshold
  double qmmode = max(1.0e-300, QmApprox(N_i_s, mMode, s, o)),
         qkmode = max(1.0e-300, QmApprox(N_i_t, kMode, t - s, o));
  double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-100;
  int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
  int dmmode = static_cast<int>(ceil(GriffithsParas(N_i_s, s).first)),
      dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
  bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
  int dkmode = static_cast<int>(ceil(GriffithsParas(N_i_t, t - s).first));
  while (!dmSwitch) {
    emCOldInc = emCInc;
    double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
    int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
    bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
    while (!dkSwitch) {
      ekCOldInc = ekCInc;
      ekCInc = QmApprox(N_i_t, dk, t - s, o) *
               calculate_expectation(N_i, dm, dk, x, z);
      ekC += ekCInc;
      if (!(dkDownSwitch))  /// Switching mechanism for j
      {
        if (sgn(dk - dkmode) <= 0) {
          dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
        }
      }

      if (!(dkUpSwitch)) {
        if (sgn(dk - dkmode) >= 0) {
          dkUpSwitch = (ekCInc < eC_threshold);
        }
      }

      dkSwitch = (dkDownSwitch && dkUpSwitch);
      if (!dkSwitch) {
        if (dkFlip == 1) {
          dkU++;
          dk = dkmode + dkU;
          dkFlip *= (dkDownSwitch ? 1 : -1);
        } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
          dkD++;
          dk = dkmode - dkD;
          dkFlip *= (dkUpSwitch ? 1 : -1);
        }
      }
    }
    emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
    eC += emCInc;

    if (!(dmDownSwitch))  /// Switching mechanism for j
    {
      if (sgn(dm - dmmode) <= 0) {
        dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
      }
    }

    if (!(dmUpSwitch)) {
      if (sgn(dm - dmmode) >= 0) {
        dmUpSwitch = (emCInc < eC_threshold);
      }
    }

    dmSwitch = (dmDownSwitch && dmUpSwitch);
    if (!dmSwitch) {
      if (dmFlip == 1) {
        dmU++;
        dm = dmmode + dmU;
        dmFlip *= (dmDownSwitch ? 1 : -1);
      } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
        dmD++;
        dm = dmmode - dmD;
        dmFlip *= (dmUpSwitch ? 1 : -1);
      }
    }
  }
  double
      density = 0.0,
      threshold = max(
          exp(log(qmmode) + log(qkmode) +
              ((!(x > 0.0))
                   ? factorials[kMode] - factorials[indexMode] -
                         factorials[kMode - indexMode] +
                         lg_theta[N_i_s][mMode] - lg_theta1[N_i_s][0] -
                         lg_theta2[N_i_s][mMode] + lg_theta[N_i_t][kMode] -
                         lg_theta1[N_i_t][indexMode] -
                         lg_theta2[N_i_t][kMode - indexMode] +
                         static_cast<double>(thetaP[N_i_t][0] + indexMode -
                                             1.0) *
                             logz +
                         static_cast<double>(thetaP[N_i_t][1] + kMode -
                                             indexMode - 1.0) *
                             log1z +
                         static_cast<double>(thetaP[N_i_s][0] + indexMode -
                                             1.0) *
                             logy +
                         static_cast<double>(thetaP[N_i_s][1] + mMode + kMode -
                                             indexMode - 1.0) *
                             log1y
                   : (((!(x < 1.0))
                           ? factorials[kMode] - factorials[indexMode] -
                                 factorials[kMode - indexMode] +
                                 lg_theta[N_i_s][mMode] -
                                 lg_theta1[N_i_s][mMode] - lg_theta2[N_i_s][0] +
                                 lg_theta[N_i_t][kMode] -
                                 lg_theta1[N_i_t][indexMode] -
                                 lg_theta2[N_i_t][kMode - indexMode] +
                                 static_cast<double>(thetaP[N_i_t][0] +
                                                     indexMode - 1.0) *
                                     logz +
                                 static_cast<double>(thetaP[N_i_t][1] + kMode -
                                                     indexMode - 1.0) *
                                     log1z +
                                 static_cast<double>(thetaP[N_i_s][0] + mMode +
                                                     indexMode - 1.0) *
                                     logy +
                                 static_cast<double>(thetaP[N_i_s][1] + kMode -
                                                     indexMode - 1.0) *
                                     log1y
                           : ((!(z > 0.0))
                                  ? factorials[mMode] - factorials[indexMode] -
                                        factorials[mMode - indexMode] +
                                        lg_theta[N_i_s][mMode] -
                                        lg_theta1[N_i_s][indexMode] -
                                        lg_theta2[N_i_s][mMode - indexMode] +
                                        lg_theta[N_i_t][kMode] -
                                        lg_theta1[N_i_t][0] -
                                        lg_theta2[N_i_t][kMode] +
                                        static_cast<double>(indexMode) * logx +
                                        static_cast<double>(mMode - indexMode) *
                                            log1x +
                                        static_cast<double>(thetaP[N_i_s][0] +
                                                            indexMode - 1.0) *
                                            logy +
                                        static_cast<double>(thetaP[N_i_s][1] +
                                                            mMode + kMode -
                                                            indexMode - 1.0) *
                                            log1y
                                  : factorials[mMode] - factorials[indexMode] -
                                        factorials[mMode - indexMode] +
                                        lg_theta[N_i_s][mMode] -
                                        lg_theta1[N_i_s][indexMode] -
                                        lg_theta2[N_i_s][mMode - indexMode] +
                                        lg_theta[N_i_t][kMode] -
                                        lg_theta1[N_i_t][kMode] -
                                        lg_theta2[N_i_t][0] +
                                        static_cast<double>(indexMode) * logx +
                                        static_cast<double>(mMode - indexMode) *
                                            log1x +
                                        static_cast<double>(thetaP[N_i_s][0] +
                                                            kMode + indexMode -
                                                            1.0) *
                                            logy +
                                        static_cast<double>(thetaP[N_i_s][1] +
                                                            mMode - indexMode -
                                                            1.0) *
                                            log1y))))) *
              1.0e-8,
          1.0e-100);

  double constants =
      static_cast<double>(thetaP[N_i_s][0] - 1.0) * logy +
      static_cast<double>(thetaP[N_i_s][1] - 1.0) * log1y +
      ((!(x > 0.0))
           ? static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                 static_cast<double>(thetaP[N_i_t][1] - 1.0) * log1z -
                 lg_theta1[N_i_s][0]
           : ((!(x < 1.0))
                  ? static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                        static_cast<double>(thetaP[N_i_t][1] - 1.0) * log1z -
                        lg_theta2[N_i_s][0]
                  : (!(z > 0.0) ? -lg_theta1[N_i_t][0]
                                : -lg_theta2[N_i_t][0])));
  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double mContr;  /// Compute contributions depending only on m
  /// Allows to avoid calling beta functions, and thus faster
  while (!mSwitch) {
    double qm = QmApprox(N_i_s, m, s, o);
    mContr = lg_theta[N_i_s][m] +
             ((!(x > 0.0))
                  ? -lg_theta2[N_i_s][m] + static_cast<double>(m) * log1y
                  : ((!(x < 1.0))
                         ? -lg_theta1[N_i_s][m] + static_cast<double>(m) * logy
                         : factorials[m] + static_cast<double>(m) * log1x +
                               static_cast<double>(m) * log1y));
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    double kContr;  /// Calculate k contributions

    while (!kSwitch) {
      double qk = QmApprox(N_i_t, k, t - s, o);
      kContr =
          lg_theta[N_i_t][k] +
          ((!(z < 1.0))
               ? -lg_theta1[N_i_t][k] + static_cast<double>(k) * logy
               : (!(z > 0.0)
                      ? -lg_theta2[N_i_t][k] + static_cast<double>(k) * log1y
                      : factorials[k] + static_cast<double>(k) * log1z +
                            static_cast<double>(k) * log1y));
      kContr += log(qk);

      int index_upper = (!(x > 0.0) || !(x < 1.0)) ? k : m;
      int other_upper = (index_upper == m) ? k : m;
      int indexFlip = 1, indexU = 0, indexD = 0,
          newindexMode = min(indexMode, index_upper),
          index = newindexMode;  /// Need to ensure l <= m since m changes!
      bool indexSwitch = false, indexDownSwitch = false, indexUpSwitch = false;
      double indexContr;

      while (!indexSwitch) {
        assert((index >= 0) && (index <= index_upper));
        indexContr = -factorials[index] - factorials[index_upper - index] +
                     ((!(x > 0.0) || !(x < 1.0))
                          ? -lg_theta1[N_i_t][index] -
                                lg_theta2[N_i_t][index_upper - index] +
                                static_cast<double>(index) *
                                    (logz - log1z + logy - log1y)
                          : -lg_theta1[N_i_s][index] -
                                lg_theta2[N_i_s][index_upper - index] +
                                static_cast<double>(index) *
                                    (logx - log1x + logy - log1y));

        double density_inc =
            exp(mContr + kContr + indexContr + constants -
                log(eC));  /// Put all separate contributions together
        density += density_inc;

        if (!(indexDownSwitch))  /// Check whether we can still move downwards
        {
          if (sgn(index - newindexMode) <= 0) {
            indexDownSwitch =
                ((density_inc < threshold) || (newindexMode - indexD - 1) < 0);
          }
        }

        if (!(indexUpSwitch))  /// Check whether we can still move downwards
        {
          if (sgn(index - newindexMode) >= 0) {
            indexUpSwitch = ((density_inc < threshold) ||
                             (newindexMode + indexU + 1) > index_upper);
          }
        }

        indexSwitch =
            (indexDownSwitch &&
             indexUpSwitch);  /// Decide if we can move out and change l

        if (!indexSwitch)  /// If we cannot, we need to move upwards or
                           /// downwards
        {
          if ((indexFlip == 1 && (newindexMode + indexU + 1 <= index_upper)) ||
              (indexDownSwitch && !(indexUpSwitch))) {
            indexU++;
            index = newindexMode + indexU;
            indexFlip *= (indexDownSwitch ? 1 : -1);
          } else if ((indexFlip == -1 && (newindexMode - indexD - 1 >= 0)) ||
                     (indexUpSwitch && !(indexDownSwitch))) {
            indexD++;
            index = newindexMode - indexD;
            indexFlip *= (indexUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(kDownSwitch))  /// Same switching procedure as for l
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch =
              (((indexU == 0) && (indexD == 0)) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = ((indexU == 0) && (indexD == 0));
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

double WrightFisher::LogSumExp(
    vector<double> &vecProb,
    double maxProb)  /// Routine to perform log-sum-exp trick
{
  double sumexp = 0.0;
  for (vector<double>::iterator vPi = vecProb.begin(); vPi != vecProb.end();
       vPi++) {
    sumexp += exp(*vPi - maxProb);
  }
  for (vector<double>::iterator vPi = vecProb.begin(); vPi != vecProb.end();
       vPi++) {
    *vPi -= maxProb + log(sumexp);
    *vPi = exp(*vPi);
  }

  return exp(maxProb + log(sumexp));
}

/// DIFFUSION SIMULATION - NEUTRAL PATHS

pair<int, int> WrightFisher::DrawAncestralProcess(
    size_t N_i, double t, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of the Ancestral
                                  /// Process using upper and lower sums
{
  assert(t > 0.0);
  /// Pre-computing all necessary quantities
  int mup = -1, mdown = 0, mmode = GriffithsParas(N_i, t).first, m = -1;
  bool m_found = false;

  /// Setting up all the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen), currsumU = 0.0, currsumL = 0.0;
  map<int, double> dL, dU;
  map<int, int> n_used;

  int threshold = ((thetaP[N_i].empty()) ? 1 : 0);

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
      double newcoefficientU =
          (((n - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double>(N_i, n, m) +
              static_cast<double>(-n * (n + theta[N_i] - 1) * t / 2.0));
      double newcoefficientL =
          (((n + 1 - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double>(N_i, n + 1, m) +
              static_cast<double>(-(n + 1) * ((n + 1) + theta[N_i] - 1) * t /
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
      return make_pair(DrawAncestralProcessG1984(N_i, t, gen), -1);
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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      for (int k = mmode - mdown; k <= mmode + mup; ++k) {
        ++n_used[k];
        dU[k] =
            dL[k] +
            (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                exp(Getlogakm<double>(N_i, n_used[k], k) +
                    static_cast<double>(
                        -n_used[k] * (n_used[k] + theta[N_i] - 1) * t / 2.0));

        ++n_used[k];
        dL[k] =
            dU[k] +
            (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                exp(Getlogakm<double>(N_i, n_used[k], k) +
                    static_cast<double>(
                        -n_used[k] * (n_used[k] + theta[N_i] - 1) * t / 2.0));

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
    size_t N_i, double t, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of the Ancestral
                                  /// Process when conditioning the diffusion on
                                  /// non-absorption but started from boundary
{
  assert(t > 0.0);
  /// Pre-computing all necessary quantities
  int mup = -1, mdown = 0,
      mmode = static_cast<int>(round(GriffithsParas(N_i, t).first)), m = -1,
      eCindex_computed = 0, q1index = 0;
  bool m_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen), currsumU = 0.0, currsumL = 0.0, eCL = 0.0, eCU,
         q1L = 0.0, q1U = 0.0;
  double eCnew = (thetaP[N_i].empty()
                      ? 1.0
                      : static_cast<double>((theta[N_i] + 1.0)) *
                            exp(-t * static_cast<double>(theta[N_i]) * 0.5)),
         consForD = (thetaP[N_i].empty()) ? 3.0 : 5.0;
  map<int, double> eAL, eAU, dL, dU;
  map<int, int> n_used, v_used;

  int threshold = ((thetaP[N_i].empty()) ? 2 : 1);

  while (!m_found) {
    if (mup > mdown &&
        mmode - mdown >
            threshold)  /// Checking whether down moves still allow - makes
                        /// sure m does not go below allowed values
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
    pair<vector<int>, double> Ct;
    Ct.second = t;
    int F1 = computeC(N_i, m, Ct),
        F2 = static_cast<int>(ceil(
            max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
        F3 = static_cast<int>(ceil((log(consForD) / t) - (theta[N_i] / 2.0)));
    while ((theta[N_i] + static_cast<double>(2 * F2 + 1)) *
               exp(-(static_cast<double>(2 * F2) + theta[N_i]) *
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
      double newcoefficientU =
          (((n - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double>(N_i, n, m) +
              static_cast<double>(-n * (n + theta[N_i] - 1) * t / 2.0));
      double newcoefficientL =
          (((n + 1 - m) % 2 != 0) ? -1.0 : 1.0) *
          exp(Getlogakm<double>(N_i, n + 1, m) +
              static_cast<double>(-(n + 1) * ((n + 1) + theta[N_i] - 1) * t /
                                  2.0));

      if ((thetaP[N_i].empty()) &&
          (2 * (n1 + 1) > q1index))  /// If mutation vector is empty, we need
                                     /// to subtract q_1 from normalising
                                     /// constant (because m can only be >= 2)
      {
        double newcoefficientUq1 =
                   (((n1 - 1) % 2 != 0) ? -1.0 : 1.0) *
                   exp(Getlogakm<double>(N_i, n1, 1) +
                       static_cast<double>(-n1 * (n1 + theta[N_i] - 1) * t /
                                           2.0)),
               newcoefficientLq1 =
                   ((n1 % 2 != 0) ? -1.0 : 1.0) *
                   exp(Getlogakm<double>(N_i, n1 + 1, 1) +
                       static_cast<double>(
                           -(n1 + 1) * ((n1 + 1) + theta[N_i] - 1) * t / 2.0));
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
        eCnew = exp(static_cast<double>(-(v + 1) * (v + theta[N_i]) * t / 2.0) +
                    log(static_cast<double>(2 * (v + 1) + theta[N_i] - 1)));
        eCU = eCL +
              (eCnew / (1.0 - consForD * exp(-static_cast<double>(
                                                 v + 1 + (theta[N_i] / 2.0)) *
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
      return make_pair(DrawAncestralProcessG1984(N_i, t, gen), -1);
    }

    dU[m] =
        ((!(eAU[m] > 0.0) || eCL < 0.0 || eCL < q1U)
             ? static_cast<double>(nan(""))
             : exp(log(static_cast<double>(m)) + log(eAU[m])) / (eCL - q1U));
    dL[m] =
        ((!(eAL[m] > 0.0) || eCU < q1L)
             ? static_cast<double>(0.0)
             : exp(log(static_cast<double>(m)) + log(eAL[m])) / (eCU - q1L));

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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      for (int k = max(mmode, threshold) - mdown;
           k <= max(mmode, threshold) + mup; ++k) {
        ++n_used[k];
        ++v_used[k];  /// Adding on more contributions to sharpen bounds

        if (thetaP[N_i].empty() && ((n_used[k] + 1) > q1index)) {
          double newcoefficientUq1 =
                     (((n_used[k] - 1) % 2 != 0) ? -1.0 : 1.0) *
                     exp(Getlogakm<double>(N_i, n_used[k], 1) +
                         static_cast<double>(-n_used[k] *
                                             (n_used[k] + theta[N_i] - 1) * t /
                                             2.0)),
                 newcoefficientLq1 =
                     (((n_used[k]) % 2 != 0) ? -1.0 : 1.0) *
                     exp(Getlogakm<double>(N_i, n_used[k] + 1, 1) +
                         static_cast<double>(
                             -(n_used[k] + 1) *
                             ((n_used[k] + 1) + theta[N_i] - 1) * t / 2.0));
          q1U = q1L + newcoefficientUq1;
          q1L = q1U + newcoefficientLq1;
        }

        if ((v_used[k] + 1) > eCindex_computed) {
          eCL += eCnew;
          eCnew = exp(
              static_cast<double>(-(v_used[k] + 2) *
                                  (v_used[k] + theta[N_i] + 1) * t / 2.0) +
              log(static_cast<double>(2 * (v_used[k] + 2) + theta[N_i] - 1)));
          eCU = eCL +
                (eCnew /
                 (1.0 - consForD * exp(-static_cast<double>(v_used[k] + 2 +
                                                            theta[N_i] / 2.0) *
                                       t)));
          eCindex_computed += 1;
        }

        eAU[k] += (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                  exp(Getlogakm<double>(N_i, n_used[k], k) +
                      static_cast<double>(
                          -n_used[k] * (n_used[k] + theta[N_i] - 1) * t / 2.0));
        ++n_used[k];
        eAL[k] += (((n_used[k] - k) % 2 != 0) ? -1.0 : 1.0) *
                  exp(Getlogakm<double>(N_i, n_used[k], k) +
                      static_cast<double>(
                          -n_used[k] * (n_used[k] + theta[N_i] - 1) * t / 2.0));

        dU[k] = ((eCL < 0.0 || !(eAU[k] > 0.0) || eCL < q1U)
                     ? static_cast<double>(nan(""))
                     : exp(log(static_cast<double>(k)) + log(eAU[k])) /
                           (eCL - q1U));
        dL[k] = ((eAL[k] < 0.0 || eCU < q1L)
                     ? static_cast<double>(0.0)
                     : exp(log(static_cast<double>(k)) + log(eAL[k])) /
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
    size_t N_i, double t, double x, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of the Ancestral Process when
               /// conditioning the diffusion on non-absorption but started
               /// from inside (0,1)
{
  assert((x > 0.0) && (x <= 1.0) && (t > 0.0));
  /// Pre-computing all necessary quantities
  int mup = -1, mdown = 0,
      mmode = static_cast<int>(round(GriffithsParas(N_i, t).first)), m = -1,
      eCindex_computed = 0, lstore;
  bool ml_found = false;

  /// Setting up all the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU = Getd2(N_i, eCvec, 0, x, t);
  map<int, double> eAL, eAU, dL, dU;

  map<int, int> n_used, v_used;

  int threshold = ((thetaP[N_i].empty()) ? 2 : 1);

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
    pair<vector<int>, double> Ct;
    Ct.second = t;
    int F1 = computeC(N_i, m, Ct),
        F2 = static_cast<int>(ceil(
            max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
        F3 = static_cast<int>(ceil((log(3.0) / t) - (theta[N_i] / 2.0)));
    while ((theta[N_i] + static_cast<double>(2 * F2 + 1)) *
               exp(-(static_cast<double>(2 * F2) + theta[N_i]) *
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
      double newcoefficientU =
                 (((n - m) % 2 != 0) ? -1.0 : 1.0) *
                 exp(Getlogakm<double>(N_i, n, m) +
                     static_cast<double>(-n * (n + theta[N_i] - 1) * t / 2.0)),
             newcoefficientL =
                 (((n + 1 - m) % 2 != 0) ? -1.0 : 1.0) *
                 exp(Getlogakm<double>(N_i, n + 1, m) +
                     static_cast<double>(-(n + 1) * ((n + 1) + theta[N_i] - 1) *
                                         t / 2.0));

      if (2 * (v + 1) > eCindex_computed)  /// Computing denominator
      {
        assert(2 * v == eCindex_computed);
        eCL = eCU - Getd2(N_i, eCvec, 2 * v + 1, x, t);
        eCU = eCL + Getd2(N_i, eCvec, 2 * v + 2, x, t);
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
      m = DrawAncestralProcessG1984(N_i, t, gen);
      vector<double> binomialBins;
      boost::math::binomial_distribution<double> BINOMIAL(m, x);
      int mlimit = (threshold == 2) ? m - 1 : m;
      for (int k = 1; k <= mlimit; ++k) {
        binomialBins.push_back(pdf(BINOMIAL, k));
      }
      boost::random::discrete_distribution<> CATEGORICAL(binomialBins.begin(),
                                                         binomialBins.end());
      int l = CATEGORICAL(gen) + 1;

      return make_pair(make_pair(m, l), -1);
    }

    dU[m] = (eCL < 0.0 ? static_cast<double>(nan("")) : eAU[m] / eCL);
    dL[m] = (eAL[m] < 0.0 ? static_cast<double>(0.0) : eAL[m] / eCU);

    int llimit = (threshold == 2) ? m - 1 : m;
    boost::math::binomial_distribution<double> BIN(m, x);
    for (int l = 1; l <= llimit && !ml_found;
         ++l)  /// Adding on the corresponding binomial terms (minus the edge
               /// cases i.e. l=0 and potentially l=m (depending on theta))
    {
      double currsummaddon =
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
        double currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const map<int, double> dUold(dU), dLold(dL);
        for (int k = max(mmode, threshold) - mdown;
             k <= max(mmode, threshold) + mup; ++k) {
          ++v_used[k];
          ++n_used[k];  /// Adding on new contributions to sharpen bounds
          double newcoefficientU =
                     exp(Getlogakm<double>(N_i, k + 2 * n_used[k], k) +
                         static_cast<double>(
                             -(k + 2 * v_used[k]) *
                             (k + 2 * n_used[k] + theta[N_i] - 1) * t / 2.0)),
                 newcoefficientL = exp(
                     Getlogakm<double>(N_i, k + 2 * n_used[k] + 1, k) +
                     static_cast<double>(
                         -(k + 2 * n_used[k] + 1) *
                         (k + 2 * n_used[k] + 1 + theta[N_i] - 1) * t / 2.0));

          eAU[k] = eAL[k] + newcoefficientU;
          eAL[k] = eAU[k] - newcoefficientL;

          if (2 * v_used[k] + 2 > eCindex_computed) {
            assert(2 * v_used[k] == eCindex_computed);
            eCL = eCU - Getd2(N_i, eCvec, 2 * v_used[k] + 1, x, t);
            eCU = eCL + Getd2(N_i, eCvec, 2 * v_used[k] + 2, x, t);
            eCindex_computed += 2;
          }

          dU[k] = (eCL < 0.0 ? static_cast<double>(nan("")) : eAU[k] / eCL);
          dL[k] = (eAL[k] < 0.0 ? static_cast<double>(0.0) : eAL[k] / eCU);

          int l1limit = (threshold == 2) ? k - 1 : k;
          boost::math::binomial_distribution<double> BIN2(k, x);
          for (int l1 = 1; l1 <= l1limit; ++l1) {
            double addon = ((!(x < 1.0) || !(x > 0.0)) ? 1.0 : pdf(BIN2, l1));
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
          // cerr << "Error: currsumLold = " << currsumLold << " > " <<
          // currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,l) = (" << n << "," << m
                    << "," << l << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumLold = " << currsumLold << " > " <<
          // currsumL
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
    size_t N_i, double t,
    boost::random::mt19937
        &gen)  /// Draws from the law of the Ancestral Process
               /// when time increment falls below threshold
{
  assert(t > 0.0);
  int threshold = (thetaP[N_i].empty()                                    ? 2
                   : (!(thetaP[N_i][0] > 0.0) || !(thetaP[N_i][1] > 0.0)) ? 1
                                                                          : 0);
  /// Pre-compute all necessary quantities
  double mean = GriffithsParas(N_i, t).first, v = GriffithsParas(N_i, t).second;
  assert(v > 0.0);
  boost::random::normal_distribution<double> NORMAL(mean, sqrt(v));
  return max(static_cast<int>(round(NORMAL(gen))),
             threshold);  // Return discretised Gaussian approximation
}

int WrightFisher::DrawSizebiasedAncestralProcess(
    size_t N_i, int d, double t,
    boost::random::mt19937
        &gen)  /// Draws from the size-biased distribution of the Ancestral
               /// process when the time increment falls below threshold
{
  assert((d > 0) && (t > 0.0));
  /// Pre-compute all necessary quantities
  double mean = GriffithsParas(N_i, t).first, v = GriffithsParas(N_i, t).second;
  assert(v > 0.0);
  /// Normalisation obtained from integral_{0}^{inf} of x*Gaussian pdf(x) dx
  double normalisation =
      0.5 * sqrt(v) *
      (exp(-pow(mean, 2.0) * 0.5 * (1.0 / v)) *
           sqrt(2.0 / boost::math::constants::pi<double>()) +
       (mean / sqrt(v)) * (1.0 + boost::math::erf(mean / (sqrt(2.0)))));

  bool decision = false;
  int flip = 1, mlimit = thetaP[N_i].empty() ? 2 : 1,
      mode = static_cast<int>(round(mean)), m = mode, jp = 0,
      jm = 0;  /// Starting m from mode to speed up routine
  double mlim = static_cast<double>(mlimit) + 0.5, summqm = 0.0;

  boost::random::uniform_01<double> U01;
  double u = U01(gen);

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

pair<double, int> WrightFisher::DrawEndpoint(
    size_t N_i, double x, double t1, double t2, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of a Wright-Fisher
                                  /// diffusion conditioned on non-absorption
{
  // std::cout << "x is " << x << ", t1 is " << t1 << ", t2 is " << t2
  //           << std::endl;
  assert((x >= 0.0) && (x <= 1.0) && (t1 < t2));
  double y;
  int m, coeffcount;
  pair<int, int> m_and_coeffcount;

  if (thetaP[N_i]
          .empty())  /// No mutation - condition on avoiding either boundary
  {
    if (!(x > 0.0))  /// Diffusion conditioned on non-absorption but started
                     /// from 0
    {
      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(N_i, 1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount =
            DrawAncestralProcessConditionalZero(N_i, t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(1.0, 1.0),
          GAMMA2(static_cast<double>(m) - 1.0, 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    } else if (!(x < 1.0))  /// Diffusion conditioned on non-absorption but
                            /// started from 1
    {
      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(N_i, 1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount =
            DrawAncestralProcessConditionalZero(N_i, t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(static_cast<double>(m) - 1.0,
                                                 1.0),
          GAMMA2(1.0, 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
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
        m = DrawAncestralProcessG1984(N_i, t2 - t1, gen);
        coeffcount = -1;
        vector<double> binomialBins;
        boost::math::binomial_distribution<double> BINOMIAL(
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
            DrawAncestralProcessConditionalInterior(N_i, t2 - t1, x, o, gen);
        m = (storer.first).first;
        l = (storer.first).second;
        coeffcount = storer.second;
      }

      boost::random::gamma_distribution<> GAMMA1(static_cast<double>(l), 1.0),
          GAMMA2(static_cast<double>(m - l), 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    }
  } else if ((!(thetaP[N_i][0] > 0.0)) ||
             (!(thetaP[N_i][1] > 0.0)))  /// One sided mutation - condition on
                                         /// avoiding one boundary only
  {
    if ((!(thetaP[N_i][0] > 0.0)) &&
        (!(x >
           0.0)))  /// Diffusion conditioned on avoiding 0 but started from 0
    {
      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(N_i, 1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount =
            DrawAncestralProcessConditionalZero(N_i, t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(1.0, 1.0),
          GAMMA2(thetaP[N_i][1] + static_cast<double>(m) - 1.0, 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    } else if ((!(thetaP[N_i][1] > 0.0)) &&
               (!(x < 1.0)))  /// Diffusion conditioned on avoiding 1 but
                              /// started from 1 - can translate to diffusion
                              /// with mutation entries swapped avoiding 0
                              /// started from 0 by symmetry
    {
      iter_swap(thetaP[N_i].begin(),
                thetaP[N_i].begin() + 1);  /// Swap mutation entries

      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawSizebiasedAncestralProcess(N_i, 1, t2 - t1, gen);
        coeffcount = -1;
      } else  /// Otherwise use alternating series technique
      {
        m_and_coeffcount =
            DrawAncestralProcessConditionalZero(N_i, t2 - t1, o, gen);
        m = m_and_coeffcount.first;
        coeffcount = m_and_coeffcount.second;
      }

      boost::random::gamma_distribution<> GAMMA1(1.0, 1.0),
          GAMMA2(thetaP[N_i][1] + static_cast<double>(m) - 1.0, 1.0);

      double y1 = -1.0;
      while (!(0.0 < y1 && y1 < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y1 = G1 / (G1 + G2);
      }

      y = 1.0 - y1;

      iter_swap(thetaP[N_i].begin(),
                thetaP[N_i].begin() + 1);  /// Swap back mutation entries
    } else if (!(thetaP[N_i][1] >
                 0.0))  /// Diffusion conditioned on avoiding 1 but started
                        /// from x in [0,1) - can translate to diffusion with
                        /// mutation entries swapped avoiding 0 started from
                        /// 1-x in interior (0,1] by symmetry
    {
      iter_swap(thetaP[N_i].begin(),
                thetaP[N_i].begin() + 1);  /// Swap mutation entries

      int l;

      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawAncestralProcessG1984(N_i, t2 - t1, gen);
        coeffcount = -1;
        vector<double> binomialBins;
        boost::math::binomial_distribution<double> BINOMIAL(m, 1.0 - x);
        for (int k = 1; k <= m; ++k) {
          binomialBins.push_back(pdf(BINOMIAL, k));
        }
        boost::random::discrete_distribution<> CATEGORICAL(binomialBins.begin(),
                                                           binomialBins.end());
        l = CATEGORICAL(gen) + 1;
      } else  /// Otherwise use alternating series technique
      {
        pair<pair<int, int>, int> storer =
            DrawAncestralProcessConditionalInterior(N_i, t2 - t1, 1.0 - x, o,
                                                    gen);
        m = (storer.first).first;
        l = (storer.first).second;
        coeffcount = storer.second;
      }

      boost::random::gamma_distribution<> GAMMA1(
          thetaP[N_i][0] + static_cast<double>(l), 1.0),
          GAMMA2(thetaP[N_i][1] + static_cast<double>(m - l), 1.0);

      double y1 = -1.0;
      while (!(0.0 < y1 && y1 < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y1 = G1 / (G1 + G2);
      }
      y = 1.0 - y1;

      iter_swap(thetaP[N_i].begin(),
                thetaP[N_i].begin() + 1);  /// Swap mutation entries back
    } else  /// Diffusion conditioned on avoiding 0 but started from x in
            /// (0,1]
    {
      int l;

      if (t2 - t1 <=
          o.g1984threshold)  /// Time increment smaller than threshold, use
                             /// Gaussian approximations
      {
        m = DrawAncestralProcessG1984(N_i, t2 - t1, gen);
        coeffcount = -1;
        vector<double> binomialBins;
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
            DrawAncestralProcessConditionalInterior(N_i, t2 - t1, x, o, gen);
        m = (storer.first).first;
        l = (storer.first).second;
        coeffcount = storer.second;
      }

      boost::random::gamma_distribution<> GAMMA1(
          thetaP[N_i][0] + static_cast<double>(l), 1.0),
          GAMMA2(thetaP[N_i][1] + static_cast<double>(m - l), 1.0);

      y = -1.0;
      while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
      {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    }
  } else  /// Strictly positive entries in mutation vector, so unconditioned
          /// diffusion
  {
    if (t2 - t1 <= o.g1984threshold)  /// Time increment smaller than threshold,
                                      /// use Gaussian approximations
    {
      m = DrawAncestralProcessG1984(N_i, t2 - t1, gen);
      coeffcount = -1;
    } else  /// Otherwise use alternating series technique
    {
      m_and_coeffcount = DrawAncestralProcess(N_i, t2 - t1, o, gen);
      m = m_and_coeffcount.first;
      coeffcount = m_and_coeffcount.second;
    }

    boost::random::binomial_distribution<int> BINOMIAL(m,
                                                       static_cast<double>(x));
    int l = BINOMIAL(gen);

    boost::random::gamma_distribution<> GAMMA1(
        thetaP[N_i][0] + static_cast<double>(l), 1.0),
        GAMMA2(thetaP[N_i][1] + static_cast<double>(m - l), 1.0);

    y = -1.0;
    while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
    {
      double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }
  }

  return make_pair(y, coeffcount);
}

pair<double, int> WrightFisher::DrawUnconditionedDiffusion(
    size_t N_i, double x, double t, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of an unconditioned
                                  /// Wright-Fisher diffusion
{
  assert((x >= 0.0) && (x <= 1.0) && (t > 0.0));

  if (!(x > 0.0) &&
      (thetaP[N_i].empty() ||
       (!(thetaP[N_i][0] > 0.0))))  /// If we start at 0 and there is no
                                    /// mutation, diffusion stays there
  {
    return make_pair(0.0, -1);
  } else if (!(x < 1.0) &&
             (thetaP[N_i].empty() ||
              (!(thetaP[N_i][1] >
                 0.0))))  /// Similarly if started from 1 with no mutation
  {
    return make_pair(1.0, -1);
  } else  /// Otherwise we first draw M ~ {q_m}
  {
    int m, coefcount;
    if (t <= o.g1984threshold) {
      m = DrawAncestralProcessG1984(N_i, t, gen);
      coefcount = -1;
    } else {
      pair<int, int> temp = DrawAncestralProcess(N_i, t, o, gen);
      m = temp.first;
      coefcount = temp.second;
    }

    double para1, para2;
    boost::random::binomial_distribution<> BIN(m, x);  /// Draw L ~ Bin(m,x)
    int l = BIN(gen);

    if (thetaP[N_i]
            .empty())  /// Sort out cases based on what mutation parameters is
    {
      if (l == 0) {
        return make_pair(0.0, coefcount);
      } else if (l == m) {
        return make_pair(1.0, coefcount);
      } else {
        para1 = static_cast<double>(l);
        para2 = static_cast<double>(m - l);
      }
    } else if (!(thetaP[N_i][0] > 0.0)) {
      if (l == 0) {
        return make_pair(0.0, coefcount);
      } else {
        para1 = static_cast<double>(l);
        para2 = static_cast<double>(thetaP[N_i][1] + m - l);
      }
    } else if (!(thetaP[N_i][1] > 0.0)) {
      if (l == m) {
        return make_pair(1.0, coefcount);
      } else {
        para1 = static_cast<double>(thetaP[N_i][0] + l);
        para2 = static_cast<double>(m - l);
      }
    } else {
      para1 = static_cast<double>(thetaP[N_i][0] + l);
      para2 = static_cast<double>(thetaP[N_i][0] + m - l);
    }

    boost::random::gamma_distribution<double> GAMMA1(para1, 1.0),
        GAMMA2(para2, 1.0);
    double y = -1.0;
    while (!(0.0 < y && y < 1.0))  /// Occasionally get y == 0.0 or y == 1.0
    {
      double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }

    return make_pair(y, coefcount);
  }
}

/// DIFFUSION SIMULATION - NON-NEUTRAL PATHS

vector<vector<double>> WrightFisher::NonNeutralDraw(
    size_t N_i, double x, double t1, double t2, bool Absorption,
    const Options &o,
    boost::random::mt19937 &gen)  /// Draws of paths from non-neutral WF
                                  /// diffusion started from x at time t1
{
  assert((x >= 0.0) && (x <= 1.0) && (t1 < t2));
  bool accept = false;
  vector<double> paras{phiMin[N_i], phiMax[N_i], phiMax[N_i] - phiMin[N_i]};
  vector<vector<double>> ptr;
  double kapmean = paras[2] * (t2 - t1);  /// Rate of Poisson point process

  boost::random::poisson_distribution<int> kap(static_cast<double>(kapmean));

  boost::random::uniform_01<double> unift, unifm,
      unifU;  /// Set up uniform generators for points over [t1,t2] *
  /// [0,phiMax-phiMin] * [0,1]
  int rcount = 0;
  while (!accept)  /// Until all skeleton points get accepted, keep going
  {
    int kappa = kap(gen);  /// Generate kappa ~ Pois
    double u = unifU(gen);
    double small_offset = 1.0e-14;

    vector<double> path, times(kappa), marks(kappa), rejcount;
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
    marks.push_back(u);  /// Add on end point - to be checked differently to
                         /// skeleton points

    for (vector<double>::iterator itt = times.begin(), itm = marks.begin();
         itt != times.end(); itt++, itm++) {
      if (kappa == 0)  /// No skeleton points -> check end point directly
      {
        if (Absorption) {
          path.push_back(DrawUnconditionedDiffusion(N_i, x, t2 - t1, o, gen)
                             .first);  /// Generate endpoint
        } else {
          path.push_back(DrawEndpoint(N_i, x, t1, t2, o, gen)
                             .first);  /// Generate endpoint
        }

        if (exp(Atilde(N_i, path.back()) - Atildeplus(N_i)) <
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
                DrawUnconditionedDiffusion(N_i, x, (*itt) - t1, o, gen)
                    .first);  /// Generate skeleton points sequentially
          } else {
            path.push_back(
                DrawEndpoint(N_i, x, t1, *itt, o, gen)
                    .first);  /// Generate skeleton points sequentially
          }

          if (Phitilde(N_i, path.back()) - paras[0] >
              *itm)  /// Test generated point is OK, otherwise can stop and
                     /// generate a new Poisson point process
          {
            rcount++;
            break;
          }

        } else if (*itt != t2)  /// Separate first time stamp and rest to
                                /// ensure correct sequential sampling
        {
          if (Absorption) {
            path.push_back(DrawUnconditionedDiffusion(
                               N_i, path.back(), (*itt) - *(itt - 1), o, gen)
                               .first);
          } else {
            path.push_back(
                DrawEndpoint(N_i, path.back(), *(itt - 1), *itt, o, gen).first);
          }

          if (Phitilde(N_i, path.back()) - paras[0] > *itm) {
            rcount++;
            break;
          }

        } else  /// Endpoint draw
        {
          if (Absorption) {
            path.push_back(DrawUnconditionedDiffusion(
                               N_i, path.back(), (*itt) - *(itt - 1), o, gen)
                               .first);
          } else {
            path.push_back(
                DrawEndpoint(N_i, path.back(), *(itt - 1), *itt, o, gen).first);
          }

          if (exp(Atilde(N_i, path.back()) - Atildeplus(N_i)) <
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

pair<double, int> WrightFisher::NonNeutralDrawEndpoint(
    size_t N_i, double x, double t1, double t2, bool Absorption,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Invoke NonNeutralDraw to generate a whole path, but retains
               /// only the endpoint at time t2
{
  vector<vector<double>> ptr =
      NonNeutralDraw(N_i, x, t1, t2, Absorption, o, gen);
  return make_pair(ptr[0].back(), -1);
}

/// BRIDGE SIMULATION - NEUTRAL PATHS

vector<int> WrightFisher::DrawBridgePMFUnconditional(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting within (0,1),
               /// ending at some boundary point (absorption can happen at any
               /// time) with time increments being large enough
{
  assert((x > 0.0) && (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (s > 0.0) &&
         (t > 0.0));
  double logx = log(x), log1x = log(1.0 - x);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkl_found = false;

  /// Setting up necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdBridgeUnconditional(N_i, eCvec, 0, x, z, t));
  vector<int> v_used;
  double Amkl;
  int n = -1, Fmkl = 0, eCindex_computed = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkl_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, m, Cs), F2 = computeC(N_i, k, Cts);
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
      double newcoefficientU =
          exp(Getlogakm<double>(N_i, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) * (m + 2 * v + theta[N_i] - 1) *
                                  s / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, m + 2 * v + 1, m) +
              static_cast<double>(-(m + 2 * v + 1) *
                                  (m + 2 * v + 1 + theta[N_i] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double>(N_i, k + 2 * v, k) +
              static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i] - 1) *
                                  (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeUnconditional(N_i, eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeUnconditional(N_i, eCvec, 2 * v + 2, x, z, t);
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
          return DrawBridgePMFUnconditionalApprox(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
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
            (((thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0)) && !(z > 0.0))
                 ? 0
                 : 1),
        lupper =
            (((thetaP[N_i].empty() || !(thetaP[N_i][1] > 0.0)) && !(z < 1.0))
                 ? m
                 : m - 1);
    precomputeA(N_i, N_i, m, k);
    for (int l = llower; l <= lupper && !mkl_found; ++l) {
      double addon =
          (!(z > 0.0))
              ? ((l == 0) ? static_cast<double>(m) * log1x
                          : factorials[m] - factorials[l] - factorials[m - l] +
                                lg_theta[N_i][m - l + k] + lg_theta[N_i][m] -
                                lg_theta[N_i][m + k] - lg_theta[N_i][m - l])
              : ((l == lupper)
                     ? static_cast<double>(m) * logx
                     : factorials[m] - factorials[l] - factorials[m - l] +
                           lg_theta[N_i][l + k] + lg_theta[N_i][m] -
                           lg_theta[N_i][m + k] - lg_theta[N_i][l]);
      Amkl = exp(addon);
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
        return DrawBridgePMFUnconditionalApprox(N_i, x, z, s, t, o, gen);
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
        double currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
          ++v_used[i];
          double newcoefficientU =
              exp(Getlogakm<double>(N_i, mi + 2 * v_used[i], mi) +
                  static_cast<double>(-(mi + 2 * v_used[i]) *
                                      (mi + 2 * v_used[i] + theta[N_i] - 1) *
                                      s / 2.0));
          double newcoefficientL =
              exp(Getlogakm<double>(N_i, mi + 2 * v_used[i] + 1, mi) +
                  static_cast<double>(
                      -(mi + 2 * v_used[i] + 1) *
                      (mi + 2 * v_used[i] + 1 + theta[N_i] - 1) * s / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          newcoefficientU =
              exp(Getlogakm<double>(N_i, ki + 2 * v_used[i], ki) +
                  static_cast<double>(-(ki + 2 * v_used[i]) *
                                      (ki + 2 * v_used[i] + theta[N_i] - 1) *
                                      (t - s) / 2.0));
          newcoefficientL = exp(
              Getlogakm<double>(N_i, ki + 2 * v_used[i] + 1, ki) +
              static_cast<double>(-(ki + 2 * v_used[i] + 1) *
                                  (ki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                  (t - s) / 2.0));

          eBU[i] = eBL[i] + newcoefficientU;
          eBL[i] = eBU[i] - newcoefficientL;

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU - GetdBridgeUnconditional(N_i, eCvec, 2 * v_used[i] + 1,
                                                x, z, t);
            eCU = eCL + GetdBridgeUnconditional(N_i, eCvec, 2 * v_used[i] + 2,
                                                x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                             : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
          dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                       ? static_cast<double>(0.0)
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
            return DrawBridgePMFUnconditionalApprox(N_i, x, z, s, t, o, gen);
          }

          int lilower = (((thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0)) &&
                          !(z > 0.0))
                             ? 0
                             : 1),
              liupper = (((thetaP[N_i].empty() || !(thetaP[N_i][1] > 0.0)) &&
                          !(z < 1.0))
                             ? mi
                             : mi - 1);
          precomputeA(N_i, N_i, mi, ki);
          for (int l2 = lilower; l2 <= liupper; ++l2) {
            addon =
                (!(z > 0.0))
                    ? ((l2 == 0)
                           ? static_cast<double>(mi) * log1x
                           : factorials[mi] - factorials[l2] -
                                 factorials[mi - l2] +
                                 lg_theta[N_i][mi - l2 + ki] +
                                 lg_theta[N_i][mi] - lg_theta[N_i][mi + ki] -
                                 lg_theta[N_i][mi - l2])
                    : ((l2 == liupper)
                           ? static_cast<double>(mi) * logx
                           : factorials[mi] - factorials[l2] -
                                 factorials[mi - l2] + lg_theta[N_i][l2 + ki] +
                                 lg_theta[N_i][mi] - lg_theta[N_i][mi + ki] -
                                 lg_theta[N_i][l2]);
            Amkl = exp(addon);
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
          // cerr << "Error: currsumLold = " << currsumLold << " > " <<
          // currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,l) = (" << n << "," << m
                    << "," << k << "," << l << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " <<
          // currsumU
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting within (0,1),
               /// ending at some boundary point (absorption can happen at any
               /// time) with one time increment falling below the threshold
{
  assert((x > 0.0) && (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (t > s) &&
         (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);
  double logx = log(x), log1x = log(1.0 - x);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkl_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdBridgeUnconditional(N_i, eCvec, 0, x, z, t)), qApprox = 0.0,
         runningMax = 1.0e-300;
  vector<int> v_used;
  double Amkl;
  int n = -1, Fmkl = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i, s).first,
         vm = GriffithsParas(N_i, s).second,
         kmode = GriffithsParas(N_i, t - s).first,
         vk = GriffithsParas(N_i, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkl_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, mork, Csorts);
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
      double newcoefficientU = exp(
          Getlogakm<double>(N_i, mork + 2 * v, mork) +
          static_cast<double>(-(mork + 2 * v) *
                              (mork + 2 * v + theta[N_i] - 1) * sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i] - 1) * sorts /
                                  2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i, morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeUnconditional(N_i, eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeUnconditional(N_i, eCvec, 2 * v + 2, x, z, t);
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
          return DrawBridgePMFUnconditionalApprox(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double>(0.0)
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
            (((thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0)) && !(z > 0.0))
                 ? 0
                 : 1),
        lupper =
            (((thetaP[N_i].empty() || !(thetaP[N_i][1] > 0.0)) && !(z < 1.0))
                 ? m
                 : m - 1);
    precomputeA(N_i, N_i, m, k);
    for (int l = llower; l <= lupper && !mkl_found; ++l) {
      double addon =
          (!(z > 0.0))
              ? ((l == 0) ? static_cast<double>(m) * log1x
                          : factorials[m] - factorials[l] - factorials[m - l] +
                                lg_theta[N_i][m - l + k] + lg_theta[N_i][m] -
                                lg_theta[N_i][m + k] - lg_theta[N_i][m - l])
              : ((l == lupper)
                     ? static_cast<double>(m) * logx
                     : factorials[m] - factorials[l] - factorials[m - l] +
                           lg_theta[N_i][l + k] + lg_theta[N_i][m] -
                           lg_theta[N_i][m + k] - lg_theta[N_i][l]);
      Amkl = exp(addon);
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
        return DrawBridgePMFUnconditionalApprox(N_i, x, z, s, t, o, gen);
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
        double currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
              morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
          ++v_used[i];
          double newcoefficientU =
              exp(Getlogakm<double>(N_i, morki + 2 * v_used[i], morki) +
                  static_cast<double>(-(morki + 2 * v_used[i]) *
                                      (morki + 2 * v_used[i] + theta[N_i] - 1) *
                                      sorts / 2.0));
          double newcoefficientL = exp(
              Getlogakm<double>(N_i, morki + 2 * v_used[i] + 1, morki) +
              static_cast<double>(-(morki + 2 * v_used[i] + 1) *
                                  (morki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                  sorts / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          qApprox = DiscretisedNormCDF(N_i, morkapproxi, sortsapprox);

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU - GetdBridgeUnconditional(N_i, eCvec, 2 * v_used[i] + 1,
                                                x, z, t);
            eCU = eCL + GetdBridgeUnconditional(N_i, eCvec, 2 * v_used[i] + 2,
                                                x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                             : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
          dL[i] = (eAL[i] < 0.0 ? static_cast<double>(0.0)
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
            return DrawBridgePMFUnconditionalApprox(N_i, x, z, s, t, o, gen);
          }

          int lilower = (((thetaP[N_i].empty() || !(thetaP[N_i][0] > 0.0)) &&
                          !(z > 0.0))
                             ? 0
                             : 1),
              liupper = (((thetaP[N_i].empty() || !(thetaP[N_i][1] > 0.0)) &&
                          !(z < 1.0))
                             ? mi
                             : mi - 1);
          precomputeA(N_i, N_i, mi, ki);
          for (int l2 = lilower; l2 <= liupper; ++l2) {
            addon =
                (!(z > 0.0))
                    ? ((l2 == 0)
                           ? static_cast<double>(mi) * log1x
                           : factorials[mi] - factorials[l2] -
                                 factorials[mi - l2] +
                                 lg_theta[N_i][mi - l2 + ki] +
                                 lg_theta[N_i][mi] - lg_theta[N_i][mi + ki] -
                                 lg_theta[N_i][mi - l2])
                    : ((l == liupper)
                           ? static_cast<double>(mi) * logx
                           : factorials[mi] - factorials[l2] -
                                 factorials[mi - l2] + lg_theta[N_i][l2 + ki] +
                                 lg_theta[N_i][mi] - lg_theta[N_i][mi + ki] -
                                 lg_theta[N_i][l2]);
            Amkl = exp(addon);
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
          // cerr << "Error: currsumLold = " << currsumLold << " > " <<
          // currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,l) = (" << n << "," << m
                    << "," << k << "," << l << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " <<
          // currsumU
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

    if (counter == totalpts)  /// Gaussian approximation leads to currsum
                              /// summing to < 1.0, so we renormalise and
                              /// sample from the computed quantities
    {
      LogSumExp(currsumStore, runningMax);
      double sum = 0.0;
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting within (0,1),
               /// ending at some boundary point (absorption can happen at any
               /// time) with both time increments falling below the threshold
{
  assert((x > 0.0) && (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (t > s) &&
         (s > 0.0));
  double logx = log(x), log1x = log(1.0 - x);
  vector<int> returnvec;
  vector<double> currsumStore;
  vector<int> mkljStore;

  /// Compute denominator
  double eC = 0.0, eCInc = 1.0, eCOldInc = 1.0;
  int Dflip = 1, Djm = 0, Djp = 0, mkLower = (thetaP[N_i].empty() ? 1 : 0);
  int dmode = static_cast<int>(ceil(GriffithsParas(N_i, t).first)), d = dmode;

  while (max(eCOldInc, eCInc) > 0.0 || (!(eC > 0.0))) {
    eCOldInc = eCInc;
    double xcont = (!(z > 0.0) ? static_cast<double>(d) * log1x
                               : static_cast<double>(d) * logx);

    eCInc = exp(log(QmApprox(N_i, d, t, o)) + xcont);
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

  vector<int> modeGuess =
      mklModeFinder(false, N_i, x, z, s, t,
                    o);  /// Get a guess to location of mode over (m,k,l)
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2];

  boost::random::uniform_01<double>
      U01;  /// Use these guesses & eC to set a suitable threshold for
  /// subsequent computations
  double currsum = 0.0, u = U01(gen),
         threshold = exp(mklModeFinder_Evaluator(false, N_i, mMode, kMode,
                                                 lMode, x, z, s, t, o) -
                         log(eC)) *
                     1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0, kFlip = 1, kD = 0, kU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double mContr,
      runningMax = -1.0e100;  /// Computing m contributions
  precomputeA(N_i, N_i, mMode, kMode);
  while (!mSwitch) {
    precomputeA(N_i, N_i, m, kMode);
    double qm = QmApprox(N_i, m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to
                      /// very small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    mContr = factorials[m] + static_cast<double>(m) * log1x + lg_theta[N_i][m];
    mContr += log(qm);

    int k = kMode, j = (!(z > 0.0) ? 0 : k);
    kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Computing k contributions
    double kContr;

    while (!kSwitch) {
      precomputeA(N_i, N_i, m, k);
      double qk = QmApprox(N_i, k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      kContr = -lg_theta[N_i][m + k];
      kContr += log(qk);

      int lFlip = 1, newlMode = min(m, lMode), l = newlMode, lU = 0,
          lD = 0;  /// Need to redefine lMode as m might be too small!
      int lLower = ((thetaP[N_i].empty() && !(z < 1.0)) ? 1 : 0),
          lUpper = ((thetaP[N_i].empty() && !(z > 0.0)) ? m : m - 1);
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;

      boost::math::binomial_distribution<> BIN(m, x);
      double lContr =
          (!(z > 0.0))
              ? ((l == 0) ? static_cast<double>(m) * log1x
                          : -factorials[l] - factorials[m - l] +
                                static_cast<double>(l) * (logx - log1x) +
                                lg_theta[N_i][m - l + k] - lg_theta[N_i][m - l])
              : ((l == lUpper) ? static_cast<double>(m) * logx
                               : -factorials[l] - factorials[m - l] +
                                     static_cast<double>(l) * (logx - log1x) +
                                     factorials[l + k - 1] - factorials[l - 1]);

      while (!lSwitch) {
        double currsum_inc =
            (!(z > 0.0)) ? ((l == 0) ? lContr - log(eC)
                                     : mContr + kContr + lContr - log(eC))
                         : ((l == lUpper) ? lContr - log(eC)
                                          : mContr + kContr + lContr - log(eC));
        runningMax = max(currsum_inc,
                         runningMax);  /// Running max of log probabilities
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
  double sum = 0.0;
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
    size_t N_i, double x, double s, double t, const Options &o,
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
  boost::random::uniform_01<double> U01;
  vector<double> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double u = U01(gen), currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdBridgeSame(N_i, eCvec, 0, x, t));
  vector<int> v_used;
  map<vector<int>, double> Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0;

  /// Compute F ensuring theoretical convergence of lower and upper sums - F1
  /// and F2 depend on m and k but F3-5 do not. None depend on j.
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, m, Cs), F2 = computeC(N_i, k, Cts);
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
      double newcoefficientU =
          exp(Getlogakm<double>(N_i, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) * (m + 2 * v + theta[N_i] - 1) *
                                  s / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, m + 2 * v + 1, m) +
              static_cast<double>(-(m + 2 * v + 1) *
                                  (m + 2 * v + 1 + theta[N_i] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double>(N_i, k + 2 * v, k) +
              static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i] - 1) *
                                  (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeSame(N_i, eCvec, 2 * v + 1, x, t);
        eCU = eCL + GetdBridgeSame(N_i, eCvec, 2 * v + 2, x, t);
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
          return DrawBridgePMFSameMutationApprox(N_i, x, s, t, gen);
        }

        break;
      }
    }

    double addon;  /// Compute the appropriate additional terms
    precomputeA(N_i, N_i, m, k);
    if (!(x > 0.0)) {
      addon = lg_theta1[N_i][0] + lg_theta2[N_i][m + k] - lg_theta[N_i][m + k] +
              lg_theta[N_i][m] - lg_theta1[N_i][0] - lg_theta2[N_i][m] +
              lg_theta[N_i][k] - lg_theta1[N_i][0] - lg_theta2[N_i][k];
    } else {
      addon = lg_theta1[N_i][m + k] + lg_theta2[N_i][0] - lg_theta[N_i][m + k] +
              lg_theta[N_i][m] - lg_theta1[N_i][m] - lg_theta2[N_i][0] +
              lg_theta[N_i][k] - lg_theta1[N_i][k] - lg_theta2[N_i][0];
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) + addon - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        ++v_used[i];
        double newcoefficientU =
                   exp(Getlogakm<double>(N_i, mi + 2 * v_used[i], mi) +
                       static_cast<double>(
                           -(mi + 2 * v_used[i]) *
                           (mi + 2 * v_used[i] + theta[N_i] - 1) * s / 2.0)),
               newcoefficientL = exp(
                   Getlogakm<double>(N_i, mi + 2 * v_used[i] + 1, mi) +
                   static_cast<double>(
                       -(mi + 2 * v_used[i] + 1) *
                       (mi + 2 * v_used[i] + 1 + theta[N_i] - 1) * s / 2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        newcoefficientU =
            exp(Getlogakm<double>(N_i, ki + 2 * v_used[i], ki) +
                static_cast<double>(-(ki + 2 * v_used[i]) *
                                    (ki + 2 * v_used[i] + theta[N_i] - 1) *
                                    (t - s) / 2.0));
        newcoefficientL =
            exp(Getlogakm<double>(N_i, ki + 2 * v_used[i] + 1, ki) +
                static_cast<double>(-(ki + 2 * v_used[i] + 1) *
                                    (ki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                    (t - s) / 2.0));

        eBU[i] = eBL[i] + newcoefficientU;
        eBL[i] = eBU[i] - newcoefficientL;

        if (2 * (v_used[i] + 1) > eCindex_computed) {
          assert(2 * v_used[i] == eCindex_computed);
          eCL = eCU - GetdBridgeSame(N_i, eCvec, 2 * v_used[i] + 1, x, t);
          eCU = eCL + GetdBridgeSame(N_i, eCvec, 2 * v_used[i] + 2, x, t);
          eCindex_computed += 2;
        }

        precomputeA(N_i, N_i, mi, ki);
        if (!(x > 0.0)) {
          addon = lg_theta1[N_i][0] + lg_theta2[N_i][mi + ki] -
                  lg_theta[N_i][mi + ki] + lg_theta[N_i][mi] -
                  lg_theta1[N_i][0] - lg_theta2[N_i][mi] + lg_theta[N_i][ki] -
                  lg_theta1[N_i][0] - lg_theta2[N_i][ki];
        } else {
          addon = lg_theta1[N_i][mi + ki] + lg_theta2[N_i][0] -
                  lg_theta[N_i][mi + ki] + lg_theta[N_i][mi] -
                  lg_theta1[N_i][mi] - lg_theta2[N_i][0] + lg_theta[N_i][ki] -
                  lg_theta1[N_i][ki] - lg_theta2[N_i][0];
        }

        dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                           : exp(log(eAU[i]) + log(eBU[i]) + addon - log(eCL)));
        dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                     ? static_cast<double>(0.0)
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
          return DrawBridgePMFSameMutationApprox(N_i, x, s, t, gen);
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
    size_t N_i, double x, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending at the
               /// same boundary point, but one of the time increments falls
               /// below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdBridgeSame(N_i, eCvec, 0, x, t)), qApprox = 0.0,
         runningMax = 1.0e-300;
  vector<int> v_used;
  map<vector<int>, double> Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0,
      threshold =
          (thetaP[N_i].empty()
               ? 2
               : ((!(thetaP[N_i][0] > 0.0) || !(thetaP[N_i][1] > 0.0)) ? 1
                                                                       : 0));

  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i, s).first,
         vm = GriffithsParas(N_i, s).second;
  double kmode = GriffithsParas(N_i, t - s).first,
         vk = GriffithsParas(N_i, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(threshold, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(threshold, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, mork, Csorts);
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
      double newcoefficientU = exp(
          Getlogakm<double>(N_i, mork + 2 * v, mork) +
          static_cast<double>(-(mork + 2 * v) *
                              (mork + 2 * v + theta[N_i] - 1) * sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i] - 1) * sorts /
                                  2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i, morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeSame(N_i, eCvec, 2 * v + 1, x, t);
        eCU = eCL + GetdBridgeSame(N_i, eCvec, 2 * v + 2, x, t);
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
          return DrawBridgePMFSameMutationApprox(N_i, x, s, t, gen);
        }

        break;
      }
    }

    double addon;  /// Computing the appropriate additional contributions

    precomputeA(N_i, N_i, m, k);
    if (!(x > 0.0)) {
      addon = lg_theta1[N_i][0] + lg_theta2[N_i][m + k] - lg_theta[N_i][m + k] +
              lg_theta[N_i][m] - lg_theta1[N_i][0] - lg_theta2[N_i][m] +
              lg_theta[N_i][k] - lg_theta1[N_i][0] - lg_theta2[N_i][k];
    } else {
      addon = lg_theta1[N_i][m + k] + lg_theta2[N_i][0] - lg_theta[N_i][m + k] +
              lg_theta[N_i][m] - lg_theta1[N_i][m] - lg_theta2[N_i][0] +
              lg_theta[N_i][k] - lg_theta1[N_i][k] - lg_theta2[N_i][0];
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) + addon - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double>(0.0)
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
      double sum = 0.0;
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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        int morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
        ++v_used[i];
        double newcoefficientU =
            exp(Getlogakm<double>(N_i, morki + 2 * v_used[i], morki) +
                static_cast<double>(-(morki + 2 * v_used[i]) *
                                    (morki + 2 * v_used[i] + theta[N_i] - 1) *
                                    sorts / 2.0));
        double newcoefficientL = exp(
            Getlogakm<double>(N_i, morki + 2 * v_used[i] + 1, morki) +
            static_cast<double>(-(morki + 2 * v_used[i] + 1) *
                                (morki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                sorts / 2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        qApprox = DiscretisedNormCDF(N_i, morkapproxi, sortsapprox);

        if (2 * (v_used[i] + 1) > eCindex_computed) {
          assert(2 * v_used[i] == eCindex_computed);
          eCL = eCU - GetdBridgeSame(N_i, eCvec, 2 * v_used[i] + 1, x, t);
          eCU = eCL + GetdBridgeSame(N_i, eCvec, 2 * v_used[i] + 2, x, t);
          eCindex_computed += 2;
        }

        precomputeA(N_i, N_i, mi, ki);
        if (!(x > 0.0)) {
          addon = lg_theta1[N_i][0] + lg_theta2[N_i][mi + ki] -
                  lg_theta[N_i][mi + ki] + lg_theta[N_i][mi] -
                  lg_theta1[N_i][0] - lg_theta2[N_i][mi] + lg_theta[N_i][ki] -
                  lg_theta1[N_i][0] - lg_theta2[N_i][ki];
        } else {
          addon = lg_theta1[N_i][mi + ki] + lg_theta2[N_i][0] -
                  lg_theta[N_i][mi + ki] + lg_theta[N_i][mi] -
                  lg_theta1[N_i][mi] - lg_theta2[N_i][0] + lg_theta[N_i][ki] -
                  lg_theta1[N_i][ki] - lg_theta2[N_i][0];
        }

        dU[i] =
            (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[i]) + log(qApprox) + addon - log(eCL)));
        dL[i] =
            (eAL[i] < 0.0 ? static_cast<double>(0.0)
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
          return DrawBridgePMFSameMutationApprox(N_i, x, s, t, gen);
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
    size_t N_i, double x, double s, double t,
    boost::random::mt19937 &gen)  /// Draws from the law of a bridge starting
                                  /// and ending at the same boundary point, but
                                  /// both time increments fall below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  vector<int> returnvec;
  bool accept = false;
  int m, k;

  /// Return (m,k) by running a rejection sampler
  double mean = GriffithsParas(N_i, t).first, v = GriffithsParas(N_i, t).second;
  int l = max(static_cast<int>(round(mean - 4.0 * sqrt(v))), 1);
  /// Precompute an upper bound for the rejection sampler
  double norm =
      boost::math::lgamma(static_cast<double>(thetaP[N_i][1] + (2 * l))) +
      boost::math::lgamma(
          static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + l)) +
      boost::math::lgamma(
          static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + l)) -
      boost::math::lgamma(
          static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + (2 * l))) -
      boost::math::lgamma(static_cast<double>(thetaP[N_i][1] + l)) -
      boost::math::lgamma(static_cast<double>(thetaP[N_i][1] + l));

  while (!accept) {
    m = DrawSizebiasedAncestralProcess(N_i, 2, s, gen),
    k = DrawSizebiasedAncestralProcess(
        N_i, 2, t - s,
        gen);  /// Draw (m,k) from a size-biased Ancestral Process
    boost::random::uniform_01<double>
        U01;  /// Compute alpha and run an accept/reject step
    double u = U01(gen),
           alpha =
               boost::math::lgamma(
                   static_cast<double>(thetaP[N_i][1] + m + k)) +
               boost::math::lgamma(
                   static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + m)) +
               boost::math::lgamma(
                   static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + k)) -
               boost::math::lgamma(static_cast<double>(
                   thetaP[N_i][0] + thetaP[N_i][1] + m + k)) -
               boost::math::lgamma(static_cast<double>(thetaP[N_i][1] + m)) -
               boost::math::lgamma(static_cast<double>(thetaP[N_i][1] + k));

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
    size_t N_i, double s, double t, double x, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of a bridge starting
                                  /// and ending at different boundary points
                                  /// and time increments are large enough
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eAL, eAU, eBL, eBU, eCU, eCL, dL, dU;
  double currsumU = 0.0, currsumL = 0.0;
  vector<int> v_used;
  int n = -1, Fmk = 0,
      denomqindex = ((thetaP[N_i][0] > 0.0 && thetaP[N_i][1] > 0.0) ? 0 : 1);

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F4 = static_cast<int>(
      ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0)));
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, m, Cs), F2 = computeC(N_i, k, Cts),
        F3 = computeC(N_i, denomqindex, Ct);
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
      double newcoefficientU =
          exp(Getlogakm<double>(N_i, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) * (m + 2 * v + theta[N_i] - 1) *
                                  s / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, m + 2 * v + 1, m) +
              static_cast<double>(-(m + 2 * v + 1) *
                                  (m + 2 * v + 1 + theta[N_i] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double>(N_i, k + 2 * v, k) +
              static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i] - 1) *
                                  (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double>(N_i, denomqindex + 2 * v, denomqindex) +
              static_cast<double>(-(denomqindex + 2 * v) *
                                  (denomqindex + 2 * v + theta[N_i] - 1) * (t) /
                                  2.0));  // Computing q_2
      newcoefficientL =
          exp(Getlogakm<double>(N_i, denomqindex + 2 * v + 1, denomqindex) +
              static_cast<double>(-(denomqindex + 2 * v + 1) *
                                  (denomqindex + 2 * v + 1 + theta[N_i] - 1) *
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
          return DrawBridgePMFDifferentMutationApprox(N_i, s, t, x, o, gen);
        }

        break;
      }
    }

    /// Compute the corresponding additional contributions
    precomputeA(N_i, N_i, m, k);
    double addon = lg_theta[N_i][m] + lg_theta[N_i][k] - lg_theta[N_i][m + k] -
                   lg_theta1[N_i][0] - lg_theta2[N_i][0];

    dU[n] = (eCL[n] < 0.0
                 ? static_cast<double>(nan(""))
                 : exp(log(eAU[n]) + log(eBU[n]) + addon -
                       log(eCL[n] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
                 : exp(log(eAL[n]) + log(eBL[n]) + addon -
                       log(eCU[n] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));
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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        ++v_used[i];
        double newcoefficientU =
                   exp(Getlogakm<double>(N_i, mi + 2 * v_used[i], mi) +
                       static_cast<double>(
                           -(mi + 2 * v_used[i]) *
                           (mi + 2 * v_used[i] + theta[N_i] - 1) * s / 2.0)),
               newcoefficientL = exp(
                   Getlogakm<double>(N_i, mi + 2 * v_used[i] + 1, mi) +
                   static_cast<double>(
                       -(mi + 2 * v_used[i] + 1) *
                       (mi + 2 * v_used[i] + 1 + theta[N_i] - 1) * s / 2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        newcoefficientU =
            exp(Getlogakm<double>(N_i, ki + 2 * v_used[i], ki) +
                static_cast<double>(-(ki + 2 * v_used[i]) *
                                    (ki + 2 * v_used[i] + theta[N_i] - 1) *
                                    (t - s) / 2.0));
        newcoefficientL =
            exp(Getlogakm<double>(N_i, ki + 2 * v_used[i] + 1, ki) +
                static_cast<double>(-(ki + 2 * v_used[i] + 1) *
                                    (ki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                    (t - s) / 2.0));

        eBU[i] = eBL[i] + newcoefficientU;
        eBL[i] = eBU[i] - newcoefficientL;

        newcoefficientU = exp(
            Getlogakm<double>(N_i, denomqindex + 2 * v_used[i], denomqindex) +
            static_cast<double>(-(denomqindex + 2 * v_used[i]) *
                                (denomqindex + 2 * v_used[i] + theta[N_i] - 1) *
                                (t) / 2.0));
        newcoefficientL =
            exp(Getlogakm<double>(N_i, denomqindex + 2 * v_used[i] + 1,
                                  denomqindex) +
                static_cast<double>(
                    -(denomqindex + 2 * v_used[i] + 1) *
                    (denomqindex + 2 * v_used[i] + 1 + theta[N_i] - 1) * (t) /
                    2.0));

        eCU[i] = eCL[i] + newcoefficientU;
        eCL[i] = eCU[i] - newcoefficientL;

        precomputeA(N_i, N_i, mi, ki);
        addon = lg_theta[N_i][mi] + lg_theta[N_i][ki] - lg_theta[N_i][mi + ki] -
                lg_theta1[N_i][0] - lg_theta2[N_i][0];

        dU[i] =
            (eCL[i] < 0.0
                 ? static_cast<double>(nan(""))
                 : exp(log(eAU[i]) + log(eBU[i]) + addon -
                       log(eCL[i] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));
        dL[i] =
            ((eAL[i] < 0.0 || eBL[i] < 0.0)
                 ? static_cast<double>(0.0)
                 : exp(log(eAL[i]) + log(eBL[i]) + addon -
                       log(eCU[i] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));

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
          return DrawBridgePMFDifferentMutationApprox(N_i, s, t, x, o, gen);
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
    size_t N_i, double s, double t, double x, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending at
               /// different boundary points, but one of the time increments
               /// falls below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eAL, eAU, eCU, eCL, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, qApprox = 0.0, runningMax = 1.0e-300;
  vector<int> v_used;
  int n = -1, Fmk = 0,
      denomqindex = ((thetaP[N_i][0] > 0.0 && thetaP[N_i][1] > 0.0) ? 0 : 1);
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i, s).first,
         vm = GriffithsParas(N_i, s).second,
         kmode = GriffithsParas(N_i, t - s).first,
         vk = GriffithsParas(N_i, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F4 = static_cast<int>(
      ceil(max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0)));
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, mork, Csorts), F3 = computeC(N_i, denomqindex, Ct);
    if (o.debug > 2) cerr << "(F1,F3) = (" << F1 << "," << F3 << ")" << endl;
    Fmk = max(max(F1, F3), F4);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << s << "," << t << ") = " << Fmk << endl;
    }

    while (2 * v < Fmk) {
      ++v;
      double newcoefficientU = exp(
          Getlogakm<double>(N_i, mork + 2 * v, mork) +
          static_cast<double>(-(mork + 2 * v) *
                              (mork + 2 * v + theta[N_i] - 1) * sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i] - 1) * sorts /
                                  2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i, morkapprox, sortsapprox);

      newcoefficientU =
          exp(Getlogakm<double>(N_i, denomqindex + 2 * v, denomqindex) +
              static_cast<double>(-(denomqindex + 2 * v) *
                                  (denomqindex + 2 * v + theta[N_i] - 1) * (t) /
                                  2.0));  // Computing q_2
      newcoefficientL =
          exp(Getlogakm<double>(N_i, denomqindex + 2 * v + 1, denomqindex) +
              static_cast<double>(-(denomqindex + 2 * v + 1) *
                                  (denomqindex + 2 * v + 1 + theta[N_i] - 1) *
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
          return DrawBridgePMFDifferentMutationApprox(N_i, s, t, x, o, gen);
        }

        break;
      }
    }

    /// Compute the appropriate additional contributions
    precomputeA(N_i, N_i, m, k);
    double addon = lg_theta[N_i][m] + lg_theta[N_i][k] - lg_theta[N_i][m + k] -
                   lg_theta1[N_i][0] - lg_theta2[N_i][0];

    dU[n] = (eCL[n] < 0.0
                 ? static_cast<double>(nan(""))
                 : exp(log(eAU[n]) + log(qApprox) + addon -
                       log(eCL[n] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));
    dL[n] = (eAL[n] < 0.0
                 ? static_cast<double>(0.0)
                 : exp(log(eAL[n]) + log(qApprox) + addon -
                       log(eCU[n] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));
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
      double sum = 0.0;
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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1], morki = (ind1 ? mi : ki),
            morkapproxi = (morki == mi ? ki : mi);
        ++v_used[i];
        double newcoefficientU = exp(
                   Getlogakm<double>(N_i, morki + 2 * v_used[i], morki) +
                   static_cast<double>(
                       -(morki + 2 * v_used[i]) *
                       (morki + 2 * v_used[i] + theta[N_i] - 1) * sorts / 2.0)),
               newcoefficientL = exp(
                   Getlogakm<double>(N_i, morki + 2 * v_used[i] + 1, morki) +
                   static_cast<double>(
                       -(morki + 2 * v_used[i] + 1) *
                       (morki + 2 * v_used[i] + 1 + theta[N_i] - 1) * sorts /
                       2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        qApprox = DiscretisedNormCDF(N_i, morkapproxi, sortsapprox);

        newcoefficientU = exp(
            Getlogakm<double>(N_i, denomqindex + 2 * v_used[i], denomqindex) +
            static_cast<double>(-(denomqindex + 2 * v_used[i]) *
                                (denomqindex + 2 * v_used[i] + theta[N_i] - 1) *
                                (t) / 2.0));
        newcoefficientL =
            exp(Getlogakm<double>(N_i, denomqindex + 2 * v_used[i] + 1,
                                  denomqindex) +
                static_cast<double>(
                    -(denomqindex + 2 * v_used[i] + 1) *
                    (denomqindex + 2 * v_used[i] + 1 + theta[N_i] - 1) * (t) /
                    2.0));

        eCU[i] = eCL[i] + newcoefficientU;
        eCL[i] = eCU[i] - newcoefficientL;

        precomputeA(N_i, N_i, mi, ki);
        addon = lg_theta[N_i][mi] + lg_theta[N_i][ki] - lg_theta[N_i][mi + ki] -
                lg_theta1[N_i][0] - lg_theta2[N_i][0];

        dU[i] =
            (eCL[i] < 0.0
                 ? static_cast<double>(nan(""))
                 : exp(log(eAU[i]) + log(qApprox) + addon -
                       log(eCL[i] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));
        dL[i] =
            (eAL[i] < 0.0
                 ? static_cast<double>(0.0)
                 : exp(log(eAL[i]) + log(qApprox) + addon -
                       log(eCU[i] / boost::math::beta<double>(
                                        static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1])))));

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
          return DrawBridgePMFDifferentMutationApprox(N_i, s, t, x, o, gen);
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
    size_t N_i, double s, double t, double x, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending at
               /// different boundary points, but both time increments fall
               /// below threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  vector<int> returnvec;
  vector<double> currsumStore;
  vector<int> mkStore;

  /// Compute denominator depending on value of time increment t
  double eC = (t <= o.g1984threshold ? DiscretisedNormCDF(N_i, 0, t)
                                     : QmApprox(N_i, 0, t, o)) /
              boost::math::beta<double>(static_cast<double>(thetaP[N_i][0]),
                                        static_cast<double>(thetaP[N_i][1]));

  /// Compute a guess to the mode over (m,k)
  vector<int> modeGuess = mkModeFinder(false, N_i, x, 1.0 - x, s, t, o);
  int mMode = modeGuess[0],
      kMode = modeGuess[1];  /// Use this guess & eC to obtain a suitable
  /// threshold for subsequent calculations
  precomputeA(N_i, N_i, mMode, kMode);
  double constContr = -lg_theta1[N_i][0] - lg_theta2[N_i][0];
  double currsum = 0.0,
         threshold = exp(mkModeFinder_Evaluator(false, N_i, mMode, kMode, x,
                                                1.0 - x, s, t, o)) *
                     1.0e-8,
         runningMax = -1.0e100;

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  while (!mSwitch) {
    int kFlip = 1, k = kMode, kU = 0, kD = 0;
    bool kSwitch = false, kDownSwitch = false, kUpSwitch = false;
    precomputeA(N_i, N_i, m, kMode);
    while (!kSwitch) {
      precomputeA(N_i, N_i, m, k);
      double currsum_inc = log(QmApprox(N_i, m, s, o)) +
                           log(QmApprox(N_i, k, t - s, o)) + lg_theta[N_i][k] +
                           lg_theta[N_i][m] - lg_theta[N_i][m + k] +
                           constContr - log(eC);

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
  double sum = 0.0;
  int index, ind = 0;
  boost::random::uniform_01<double> U01;
  double u = U01(gen);

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
      index = indexing[ind];  /// Figure out what the correct index for
                              /// mkljStore is
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting at a boundary point
               /// but ending in the interior of (0,1) with both time
               /// increments large enough
{
  assert((!(x > 0.0) || !(x < 1.0)) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  double logz = log(z), log1z = log(1.0 - z);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkj_found = false;

  /// Setting up necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdBridgeInterior(N_i, eCvec, 0, x, z, t));
  vector<int> v_used;
  double Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
      F5 = (!(x > 0.0)
                ? static_cast<int>(floor((theta[N_i] + 2.0) /
                                             (o.eps * (thetaP[N_i][1] + 1.0)) -
                                         1)) +
                      1
                : static_cast<int>(floor((theta[N_i] + 2.0) /
                                             (o.eps * (thetaP[N_i][0] + 1.0)) -
                                         1)) +
                      1);
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, m, Cs), F2 = computeC(N_i, k, Cts);
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
      double newcoefficientU =
          exp(Getlogakm<double>(N_i, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) * (m + 2 * v + theta[N_i] - 1) *
                                  s / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, m + 2 * v + 1, m) +
              static_cast<double>(-(m + 2 * v + 1) *
                                  (m + 2 * v + 1 + theta[N_i] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double>(N_i, k + 2 * v, k) +
              static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i] - 1) *
                                  (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeInterior(N_i, eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeInterior(N_i, eCvec, 2 * v + 2, x, z, t);
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
          return DrawBridgePMFInteriorMutationApprox(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
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

    precomputeA(N_i, N_i, m, k);
    double constant_factors =
        (!(x > 0.0))
            ? (factorials[k] - lg_theta[N_i][m + k] + lg_theta[N_i][m] +
               lg_theta[N_i][k] - lg_theta1[N_i][0] - lg_theta2[N_i][m] +
               static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
               static_cast<double>(thetaP[N_i][1] + k - 1.0) * log1z)
            : (factorials[k] - lg_theta[N_i][m + k] + lg_theta[N_i][m] +
               lg_theta[N_i][k] - lg_theta1[N_i][m] - lg_theta2[N_i][0] +
               static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
               static_cast<double>(thetaP[N_i][1] + k - 1.0) * log1z);
    /// Add on j contributions
    int lindex = (!(x < 1.0) ? m : 0);
    for (int j = 0; j <= k && !mkj_found; ++j) {
      double var_stuff = (!(x > 0.0))
                             ? (-factorials[j] - factorials[k - j] +
                                lg_theta1[N_i][j] + lg_theta2[N_i][m + k - j] -
                                lg_theta1[N_i][j] - lg_theta2[N_i][k - j] +
                                static_cast<double>(j) * (logz - log1z))
                             : (-factorials[j] - factorials[k - j] +
                                lg_theta1[N_i][m + j] + lg_theta2[N_i][k - j] -
                                lg_theta1[N_i][j] - lg_theta2[N_i][k - j] +
                                static_cast<double>(j) * (logz - log1z));
      Amkj = exp(constant_factors + var_stuff);
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
        return DrawBridgePMFInteriorMutationApprox(N_i, x, z, s, t, o, gen);
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
        double currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
          ++v_used[i];
          double newcoefficientU =
              exp(Getlogakm<double>(N_i, mi + 2 * v_used[i], mi) +
                  static_cast<double>(-(mi + 2 * v_used[i]) *
                                      (mi + 2 * v_used[i] + theta[N_i] - 1) *
                                      s / 2.0));
          double newcoefficientL =
              exp(Getlogakm<double>(N_i, mi + 2 * v_used[i] + 1, mi) +
                  static_cast<double>(
                      -(mi + 2 * v_used[i] + 1) *
                      (mi + 2 * v_used[i] + 1 + theta[N_i] - 1) * s / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          newcoefficientU =
              exp(Getlogakm<double>(N_i, ki + 2 * v_used[i], ki) +
                  static_cast<double>(-(ki + 2 * v_used[i]) *
                                      (ki + 2 * v_used[i] + theta[N_i] - 1) *
                                      (t - s) / 2.0));
          newcoefficientL = exp(
              Getlogakm<double>(N_i, ki + 2 * v_used[i] + 1, ki) +
              static_cast<double>(-(ki + 2 * v_used[i] + 1) *
                                  (ki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                  (t - s) / 2.0));

          eBU[i] = eBL[i] + newcoefficientU;
          eBL[i] = eBU[i] - newcoefficientL;

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU -
                  GetdBridgeInterior(N_i, eCvec, 2 * v_used[i] + 1, x, z, t);
            eCU = eCL +
                  GetdBridgeInterior(N_i, eCvec, 2 * v_used[i] + 2, x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                             : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
          dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                       ? static_cast<double>(0.0)
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
            return DrawBridgePMFInteriorMutationApprox(N_i, x, z, s, t, o, gen);
          }

          precomputeA(N_i, N_i, mi, ki);
          double constant_factors_i =
              (!(x > 0.0))
                  ? (factorials[ki] - lg_theta[N_i][mi + ki] +
                     lg_theta[N_i][mi] + lg_theta[N_i][ki] - lg_theta1[N_i][0] -
                     lg_theta2[N_i][mi] +
                     static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
                     static_cast<double>(thetaP[N_i][1] + ki - 1.0) * log1z)
                  : (factorials[ki] - lg_theta[N_i][mi + ki] +
                     lg_theta[N_i][mi] + lg_theta[N_i][ki] -
                     lg_theta1[N_i][mi] - lg_theta2[N_i][0] +
                     static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
                     static_cast<double>(thetaP[N_i][1] + ki - 1.0) * log1z);
          int liindex = (!(x < 1.0) ? mi : 0);
          for (int j2 = 0; j2 <= ki; ++j2) {
            double var_stuff_i =
                (!(x > 0.0))
                    ? (-factorials[j2] - factorials[ki - j2] +
                       lg_theta1[N_i][j2] + lg_theta2[N_i][mi + ki - j2] -
                       lg_theta1[N_i][j2] - lg_theta2[N_i][ki - j2] +
                       static_cast<double>(j2) * (logz - log1z))
                    : (-factorials[j2] - factorials[ki - j2] +
                       lg_theta1[N_i][mi + j2] + lg_theta2[N_i][ki - j2] -
                       lg_theta1[N_i][j2] - lg_theta2[N_i][ki - j2] +
                       static_cast<double>(j2) * (logz - log1z));
            Amkj = exp(constant_factors_i + var_stuff_i);
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
          // cerr << "Error: currsumLold = " << currsumLold << " > " <<
          // currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,j) = (" << n << "," << m
                    << "," << k << "," << j << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " <<
          // currsumU
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting at a boundary point
               /// but ending in the interior of (0,1) but one time increment
               /// is below the threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);
  double logz = log(z), log1z = log(1.0 - z);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mkj_found = false;

  /// Setting up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdBridgeInterior(N_i, eCvec, 0, x, z, t)), qApprox = 0.0;
  vector<int> v_used;
  double Amkj, runningMax = 1.0e-300;
  int n = -1, Fmkj = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i, s).first,
         vm = GriffithsParas(N_i, s).second,
         kmode = GriffithsParas(N_i, t - s).first,
         vk = GriffithsParas(N_i, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(2.0 / o.eps) + 1);
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = 2 * max(max(F3, F4), F5);

  while (!mkj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, mork, Csorts);
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
      double newcoefficientU = exp(
          Getlogakm<double>(N_i, mork + 2 * v, mork) +
          static_cast<double>(-(mork + 2 * v) *
                              (mork + 2 * v + theta[N_i] - 1) * sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i] - 1) * sorts /
                                  2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i, morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdBridgeInterior(N_i, eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + GetdBridgeInterior(N_i, eCvec, 2 * v + 2, x, z, t);
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
          return DrawBridgePMFInteriorMutationApprox(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double>(0.0)
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

    precomputeA(N_i, N_i, m, k);
    double constant_factors =
        (!(x > 0.0))
            ? (factorials[k] - lg_theta[N_i][m + k] + lg_theta[N_i][m] +
               lg_theta[N_i][k] - lg_theta1[N_i][0] - lg_theta2[N_i][m] +
               static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
               static_cast<double>(thetaP[N_i][1] + k - 1.0) * log1z)
            : (factorials[k] - lg_theta[N_i][m + k] + lg_theta[N_i][m] +
               lg_theta[N_i][k] - lg_theta1[N_i][m] - lg_theta2[N_i][0] +
               static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
               static_cast<double>(thetaP[N_i][1] + k - 1.0) * log1z);
    /// Add on j contributions
    int lindex = (!(x < 1.0) ? m : 0);
    for (int j = 0; j <= k && !mkj_found; ++j) {
      double var_stuff = (!(x > 0.0))
                             ? (-factorials[j] - factorials[k - j] +
                                lg_theta1[N_i][j] + lg_theta2[N_i][m + k - j] -
                                lg_theta1[N_i][j] - lg_theta2[N_i][k - j] +
                                static_cast<double>(j) * (logz - log1z))
                             : (-factorials[j] - factorials[k - j] +
                                lg_theta1[N_i][m + j] + lg_theta2[N_i][k - j] -
                                lg_theta1[N_i][j] - lg_theta2[N_i][k - j] +
                                static_cast<double>(j) * (logz - log1z));
      Amkj = exp(constant_factors + var_stuff);
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
        return DrawBridgePMFInteriorMutationApprox(N_i, x, z, s, t, o, gen);
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
        double currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
              morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
          ++v_used[i];
          double newcoefficientU =
              exp(Getlogakm<double>(N_i, morki + 2 * v_used[i], morki) +
                  static_cast<double>(-(morki + 2 * v_used[i]) *
                                      (morki + 2 * v_used[i] + theta[N_i] - 1) *
                                      sorts / 2.0));
          double newcoefficientL = exp(
              Getlogakm<double>(N_i, morki + 2 * v_used[i] + 1, morki) +
              static_cast<double>(-(morki + 2 * v_used[i] + 1) *
                                  (morki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                  sorts / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          qApprox = DiscretisedNormCDF(N_i, morkapproxi, sortsapprox);

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU -
                  GetdBridgeInterior(N_i, eCvec, 2 * v_used[i] + 1, x, z, t);
            eCU = eCL +
                  GetdBridgeInterior(N_i, eCvec, 2 * v_used[i] + 2, x, z, t);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                             : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
          dL[i] = (eAL[i] < 0.0 ? static_cast<double>(0.0)
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
            return DrawBridgePMFInteriorMutationApprox(N_i, x, z, s, t, o, gen);
          }

          precomputeA(N_i, N_i, mi, ki);
          double constant_factors_i =
              (!(x > 0.0))
                  ? (factorials[ki] - lg_theta[N_i][mi + ki] +
                     lg_theta[N_i][mi] + lg_theta[N_i][ki] - lg_theta1[N_i][0] -
                     lg_theta2[N_i][mi] +
                     static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
                     static_cast<double>(thetaP[N_i][1] + ki - 1.0) * log1z)
                  : (factorials[ki] - lg_theta[N_i][mi + ki] +
                     lg_theta[N_i][mi] + lg_theta[N_i][ki] -
                     lg_theta1[N_i][mi] - lg_theta2[N_i][0] +
                     static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
                     static_cast<double>(thetaP[N_i][1] + ki - 1.0) * log1z);
          int liindex = (!(x < 1.0) ? mi : 0);
          for (int j2 = 0; j2 <= ki; ++j2) {
            double var_stuff_i =
                (!(x > 0.0))
                    ? (-factorials[j2] - factorials[ki - j2] +
                       lg_theta1[N_i][j2] + lg_theta2[N_i][mi + ki - j2] -
                       lg_theta1[N_i][j2] - lg_theta2[N_i][ki - j2] +
                       static_cast<double>(j2) * (logz - log1z))
                    : (-factorials[j2] - factorials[ki - j2] +
                       lg_theta1[N_i][mi + j2] + lg_theta2[N_i][ki - j2] -
                       lg_theta1[N_i][j2] - lg_theta2[N_i][ki - j2] +
                       static_cast<double>(j2) * (logz - log1z));
            Amkj = exp(constant_factors_i + var_stuff_i);
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
          // cerr << "Error: currsumLold = " << currsumLold << " > " <<
          // currsumL
          std::cout << "Error: currsumLold = " << currsumLold << " > "
                    << currsumL << " = currsumL (n,m,k,l,j) = (" << n << ","
                    << m << "," << k << "," << j << ")." << endl;
          // exit(1);
        }
        if (currsumUold < currsumU) {
          // cerr << "Error: currsumUold = " << currsumUold << " < " <<
          // currsumU
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
      double sum = 0.0;
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting at a boundary point
               /// but ending in the interior of (0,1) with both time
               /// increments below the threshold
{
  assert((!(x > 0.0) || !(x < 1.0)) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  vector<int> returnvec;
  vector<double> currsumStore;
  vector<int> mkljStore;

  /// Compute denominator
  double eC = 0.0, eCInc = 1.0, eCOldInc = 1.0;
  int Dflip = 1, Djm = 0, Djp = 0;
  int dmode = static_cast<int>(ceil(GriffithsParas(N_i, t).first)), d = dmode;

  while (max(eCOldInc, eCInc) > 0.0 || (!(eC > 0.0))) {
    eCOldInc = eCInc;
    double para1 = (!(x > 0.0) ? static_cast<double>(thetaP[N_i][0])
                               : static_cast<double>(thetaP[N_i][0] + d));
    double para2 = (!(x > 0.0) ? static_cast<double>(thetaP[N_i][1] + d)
                               : static_cast<double>(thetaP[N_i][1]));
    double zcont = (!(x > 0.0) ? static_cast<double>(d) * log(1.0 - z)
                               : static_cast<double>(d) * log(z));

    eCInc = exp(log(QmApprox(N_i, d, t, o)) + zcont -
                boost::math::lgamma(static_cast<double>(para1)) -
                boost::math::lgamma(static_cast<double>(para2)) +
                boost::math::lgamma(static_cast<double>(para1 + para2)));
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

  vector<int> modeGuess =
      mkjModeFinder(false, N_i, x, z, s, t,
                    o);  /// Get a guess to location of mode over (m,k,j)
  int mMode = modeGuess[0], kMode = modeGuess[1], jMode = modeGuess[2];

  boost::random::uniform_01<double>
      U01;  /// Use these guesses & eC to set a suitable threshold for
  /// subsequent computations
  double currsum = 0.0, u = U01(gen),
         threshold = exp(mkjModeFinder_Evaluator(false, N_i, mMode, kMode,
                                                 jMode, x, z, s, t, o) -
                         log(eC)) *
                     1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0, kFlip = 1, kD = 0, kU = 0, l;
  bool mSwitch = false, mDownSwitch = false, mUpSwitch = false;
  double constContr =
      (static_cast<double>(thetaP[N_i][0] - 1.0)) * log(z) +
      (static_cast<double>(thetaP[N_i][1] - 1.0)) * log(1.0 - z) +
      (!(x > 0.0) ? -boost::math::lgamma(static_cast<double>(thetaP[N_i][0]))
                  : -boost::math::lgamma(static_cast<double>(thetaP[N_i][1])));
  double mContr_D =
      boost::math::lgamma(
          static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + m)) +
      (!(x > 0.0)
           ? -boost::math::lgamma(static_cast<double>(thetaP[N_i][1] + m))
           : -boost::math::lgamma(static_cast<double>(thetaP[N_i][0] + m)));
  double mContr_U = mContr_D, mContr,
         runningMax = -1.0e100;  /// Computing m contributions

  while (!mSwitch) {
    double qm = QmApprox(N_i, m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to
                      /// very small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    l = (!(x > 0.0) ? 0 : m);
    if (m != mMode) {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                    (m - 1))) +
            (!(x > 0.0) ? -log(static_cast<double>(thetaP[N_i][1] + (m - 1)))
                        : -log(static_cast<double>(thetaP[N_i][0] + (m - 1))));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D +=
            -log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + (m + 1) -
                                     1)) +
            (!(x > 0.0)
                 ? log(static_cast<double>(thetaP[N_i][1] + (m + 1) - 1))
                 : log(static_cast<double>(thetaP[N_i][0] + (m + 1) - 1)));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode;
    kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Computing k contributions
    double kContr_D = (boost::math::lgamma(static_cast<double>(
                           thetaP[N_i][0] + thetaP[N_i][1] + k)) -
                       boost::math::lgamma(static_cast<double>(
                           thetaP[N_i][0] + thetaP[N_i][1] + m + k))) +
                      static_cast<double>(k) * log(1.0 - z),
           kContr_U = kContr_D, kContr;

    while (!kSwitch) {
      double qk = QmApprox(N_i, k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      if (k != kMode) {
        if (kU > kD) {
          kContr_U += log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              (k - 1))) -
                      log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              m + (k - 1))) +
                      log(1.0 - z);
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D += log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              m + (k + 1) - 1)) -
                      log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              (k + 1) - 1)) -
                      log(1.0 - z);
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      int jFlip = 1, newjMode = min(k, jMode), j = newjMode, jU = 0,
          jD = 0;  /// Need to redefine jMode as k might be too small!
      bool jSwitch = false, jDownSwitch = false, jUpSwitch = false;

      double jContr_D =
          LogBinomialCoefficientCalculator(k, j) +
          (!(x > 0.0) ? boost::math::lgamma(
                            static_cast<double>(thetaP[N_i][1] + m + k - j)) -
                            boost::math::lgamma(
                                static_cast<double>(thetaP[N_i][1] + k - j))
                      : boost::math::lgamma(
                            static_cast<double>(thetaP[N_i][0] + m + j)) -
                            boost::math::lgamma(
                                static_cast<double>(thetaP[N_i][0] + j))) +
          static_cast<double>(j) * log(z) +
          static_cast<double>(-j) * log(1.0 - z);
      double jContr_U = jContr_D, jContr;

      while (!jSwitch) {
        if (j != newjMode) {
          if (jU > jD) {
            jContr_U +=
                log(static_cast<double>(k - (j - 1))) -
                log(static_cast<double>((j - 1) + 1)) + log(z) - log(1.0 - z) +
                (!(x > 0.0)
                     ? log(static_cast<double>(thetaP[N_i][1] + k - (j - 1) -
                                               1)) -
                           log(static_cast<double>(thetaP[N_i][1] + m + k -
                                                   (j - 1) - 1))
                     : log(static_cast<double>(thetaP[N_i][0] + m + (j - 1))) -
                           log(static_cast<double>(thetaP[N_i][0] + (j - 1))));
            jContr = jContr_U;
          } else {
            jContr_D +=
                log(static_cast<double>(j + 1)) -
                log(static_cast<double>(k - (j + 1) + 1)) + log(1.0 - z) -
                log(z) +
                (!(x > 0.0)
                     ? log(static_cast<double>(thetaP[N_i][1] + m + k -
                                               (j + 1))) -
                           log(static_cast<double>(thetaP[N_i][1] + k -
                                                   (j + 1)))
                     : log(static_cast<double>(thetaP[N_i][0] + (j + 1) - 1)) -
                           log(static_cast<double>(thetaP[N_i][0] + m +
                                                   (j + 1) - 1)));
            jContr = jContr_D;
          }
        } else {
          jContr = jContr_U;
        }
        double currsum_inc = constContr + mContr + kContr + jContr - log(eC);
        runningMax = max(currsum_inc,
                         runningMax);  /// Running max of log probabilities
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
  double sum = 0.0;
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending in the
               /// interior of (0,1) with both time increments large enough
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(Getd(N_i, eCvec, 0, x, z, t));
  vector<int> v_used;
  double Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F5;
  if (thetaP[N_i].empty()) {
    F5 = static_cast<int>(
        ceil(((1.0 - x) / (1.0 - z)) + (x / z) * (1.0 + z) * o.eps));
  } else {
    F5 = static_cast<int>(ceil(
        2.0 *
        (max((theta[N_i] / thetaP[N_i][1]) * (1.0 - static_cast<double>(z)),
             (1.0 + theta[N_i]) / (1.0 - static_cast<double>(z))) *
             (1.0 - static_cast<double>(x)) +
         (1.0 + (1.0 / static_cast<double>(z))) *
             (max(
                 ((static_cast<double>(z) * theta[N_i]) + 1.0) / thetaP[N_i][0],
                 1.0)) *
             static_cast<double>(x)) /
        o.eps));
  }
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0)));
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, m, Cs), F2 = computeC(N_i, k, Cts);
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
      double newcoefficientU =
          exp(Getlogakm<double>(N_i, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) * (m + 2 * v + theta[N_i] - 1) *
                                  s / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, m + 2 * v + 1, m) +
              static_cast<double>(-(m + 2 * v + 1) *
                                  (m + 2 * v + 1 + theta[N_i] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU =
          exp(Getlogakm<double>(N_i, k + 2 * v, k) +
              static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i] - 1) *
                                  (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - Getd(N_i, eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + Getd(N_i, eCvec, 2 * v + 2, x, z, t);
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
          return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
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

    precomputeA(N_i, N_i, m, k);
    double extra_factors =
        factorials[m] + factorials[k] + static_cast<double>(m) * log1x +
        static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
        static_cast<double>(thetaP[N_i][1] + k - 1.0) * log1z +
        lg_theta[N_i][m] + lg_theta[N_i][k] - lg_theta[N_i][m + k];
    for (int l = 0; l <= m && !mklj_found; ++l) {
      for (int j = 0; j <= k && !mklj_found; ++j) {
        Amklj = exp(extra_factors - factorials[l] - factorials[m - l] -
                    factorials[j] - factorials[k - j] - lg_theta1[N_i][l] -
                    lg_theta2[N_i][m - l] - lg_theta1[N_i][j] -
                    lg_theta2[N_i][k - j] + lg_theta1[N_i][l + j] +
                    lg_theta2[N_i][m - l + k - j] +
                    static_cast<double>(l) * (logx - log1x) +
                    static_cast<double>(j) * (logz - log1z));
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
          return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
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
          double currsumLold = currsumL, currsumUold = currsumU;
          currsumL = 0.0;
          currsumU = 0.0;

          const vector<double> dUold(dU), dLold(dL);
          for (int i = 0; i <= n; ++i) {
            int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
            ++v_used[i];
            double newcoefficientU = exp(
                       Getlogakm<double>(N_i, mi + 2 * v_used[i], mi) +
                       static_cast<double>(
                           -(mi + 2 * v_used[i]) *
                           (mi + 2 * v_used[i] + theta[N_i] - 1) * s / 2.0)),
                   newcoefficientL =
                       exp(Getlogakm<double>(N_i, mi + 2 * v_used[i] + 1, mi) +
                           static_cast<double>(
                               -(mi + 2 * v_used[i] + 1) *
                               (mi + 2 * v_used[i] + 1 + theta[N_i] - 1) * s /
                               2.0));

            eAU[i] = eAL[i] + newcoefficientU;
            eAL[i] = eAU[i] - newcoefficientL;

            newcoefficientU =
                exp(Getlogakm<double>(N_i, ki + 2 * v_used[i], ki) +
                    static_cast<double>(-(ki + 2 * v_used[i]) *
                                        (ki + 2 * v_used[i] + theta[N_i] - 1) *
                                        (t - s) / 2.0));
            newcoefficientL = exp(
                Getlogakm<double>(N_i, ki + 2 * v_used[i] + 1, ki) +
                static_cast<double>(-(ki + 2 * v_used[i] + 1) *
                                    (ki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                                    (t - s) / 2.0));

            eBU[i] = eBL[i] + newcoefficientU;
            eBL[i] = eBU[i] - newcoefficientL;

            if (2 * v_used[i] + 2 > eCindex_computed) {
              assert(2 * v_used[i] == eCindex_computed);
              eCL = eCU - Getd(N_i, eCvec, 2 * v_used[i] + 1, x, z, t);
              eCU = eCL + Getd(N_i, eCvec, 2 * v_used[i] + 2, x, z, t);
              eCindex_computed += 2;
            }

            dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                               : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
            dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                         ? static_cast<double>(0.0)
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
              return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
            }
            precomputeA(N_i, N_i, mi, ki);
            double extra_factors_i =
                factorials[mi] + factorials[ki] +
                static_cast<double>(mi) * log1x +
                static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
                static_cast<double>(thetaP[N_i][1] + ki - 1.0) * log1z +
                lg_theta[N_i][mi] + lg_theta[N_i][ki] - lg_theta[N_i][mi + ki];
            int l_upper = (i == n ? l : mi);
            for (int l2 = 0; l2 <= l_upper; ++l2) {
              int j_upper = ((i == n && l2 == l_upper) ? j : ki);
              for (int j2 = 0; j2 <= j_upper; ++j2) {
                Amklj = exp(extra_factors_i - factorials[l2] -
                            factorials[mi - l2] - factorials[j2] -
                            factorials[ki - j2] - lg_theta1[N_i][l2] -
                            lg_theta2[N_i][mi - l2] - lg_theta1[N_i][j2] -
                            lg_theta2[N_i][ki - j2] + lg_theta1[N_i][l2 + j2] +
                            lg_theta2[N_i][mi - l2 + ki - j2] +
                            static_cast<double>(l2) * (logx - log1x) +
                            static_cast<double>(j2) * (logz - log1z));
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen)  /// Draws from the law of a bridge starting
                                  /// and ending in the interior of (0,1), but
                                  /// one time increment is below the threshold
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));

  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(Getd(N_i, eCvec, 0, x, z, t)), qApprox = 0.0,
         runningMax = 1.0e-300;
  vector<int> v_used;
  double Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i, s).first,
         vm = GriffithsParas(N_i, s).second,
         kmode = GriffithsParas(N_i, t - s).first,
         vk = GriffithsParas(N_i, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F5;
  if (thetaP[N_i].empty()) {
    F5 = static_cast<int>(
        ceil(((1.0 - x) / (1.0 - z)) + (x / z) * (1.0 + z) * o.eps));
  } else {
    F5 = static_cast<int>(ceil(
        2.0 *
        (max((theta[N_i] / thetaP[N_i][1]) * (1.0 - static_cast<double>(z)),
             (1.0 + theta[N_i]) / (1.0 - static_cast<double>(z))) *
             (1.0 - static_cast<double>(x)) +
         (1.0 + (1.0 / static_cast<double>(z))) *
             (max(
                 ((static_cast<double>(z) * theta[N_i]) + 1.0) / thetaP[N_i][0],
                 1.0)) *
             static_cast<double>(x)) /
        o.eps));
  }
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0)));
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(false, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, mork, Csorts);
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
      double newcoefficientU = exp(
          Getlogakm<double>(N_i, mork + 2 * v, mork) +
          static_cast<double>(-(mork + 2 * v) *
                              (mork + 2 * v + theta[N_i] - 1) * sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i] - 1) * sorts /
                                  2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i, morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - Getd(N_i, eCvec, 2 * v + 1, x, z, t);
        eCU = eCL + Getd(N_i, eCvec, 2 * v + 2, x, z, t);
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
          return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double>(0.0)
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
    precomputeA(N_i, N_i, m, k);
    double extra_factors =
        factorials[m] + factorials[k] + static_cast<double>(m) * log1x +
        static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
        static_cast<double>(thetaP[N_i][1] + k - 1.0) * log1z +
        lg_theta[N_i][m] + lg_theta[N_i][k] - lg_theta[N_i][m + k];
    for (int l = 0; l <= m && !mklj_found; ++l) {
      for (int j = 0; j <= k && !mklj_found; ++j) {
        Amklj = exp(extra_factors - factorials[l] - factorials[m - l] -
                    factorials[j] - factorials[k - j] - lg_theta1[N_i][l] -
                    lg_theta2[N_i][m - l] - lg_theta1[N_i][j] -
                    lg_theta2[N_i][k - j] + lg_theta1[N_i][l + j] +
                    lg_theta2[N_i][m - l + k - j] +
                    static_cast<double>(l) * (logx - log1x) +
                    static_cast<double>(j) * (logz - log1z));
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
          return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
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
          double currsumLold = currsumL, currsumUold = currsumU;
          currsumL = 0.0;
          currsumU = 0.0;

          const vector<double> dUold(dU), dLold(dL);
          for (int i = 0; i <= n; ++i) {
            int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
                morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
            ;
            ++v_used[i];
            double newcoefficientU = exp(
                       Getlogakm<double>(N_i, morki + 2 * v_used[i], morki) +
                       static_cast<double>(
                           -(morki + 2 * v_used[i]) *
                           (morki + 2 * v_used[i] + theta[N_i] - 1) * sorts /
                           2.0)),
                   newcoefficientL =
                       exp(Getlogakm<double>(N_i, morki + 2 * v_used[i] + 1,
                                             morki) +
                           static_cast<double>(
                               -(morki + 2 * v_used[i] + 1) *
                               (morki + 2 * v_used[i] + 1 + theta[N_i] - 1) *
                               sorts / 2.0));

            eAU[i] = eAL[i] + newcoefficientU;
            eAL[i] = eAU[i] - newcoefficientL;

            qApprox = DiscretisedNormCDF(N_i, morkapproxi, sortsapprox);

            if (2 * v_used[i] + 2 > eCindex_computed) {
              assert(2 * v_used[i] == eCindex_computed);
              eCL = eCU - Getd(N_i, eCvec, 2 * v_used[i] + 1, x, z, t);
              eCU = eCL + Getd(N_i, eCvec, 2 * v_used[i] + 2, x, z, t);
              eCindex_computed += 2;
            }

            dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                               : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
            dL[i] = (eAL[i] < 0.0 ? static_cast<double>(0.0)
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
              return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
            }

            precomputeA(N_i, N_i, mi, ki);
            double extra_factors_i =
                factorials[mi] + factorials[ki] +
                static_cast<double>(mi) * log1x +
                static_cast<double>(thetaP[N_i][0] - 1.0) * logz +
                static_cast<double>(thetaP[N_i][1] + ki - 1.0) * log1z +
                lg_theta[N_i][mi] + lg_theta[N_i][ki] - lg_theta[N_i][mi + ki];
            int l_upper = (i == n ? l : mi);
            for (int l2 = 0; l2 <= l_upper; ++l2) {
              int j_upper = ((i == n && l2 == l_upper) ? j : ki);
              for (int j2 = 0; j2 <= j_upper; ++j2) {
                Amklj = exp(extra_factors_i - factorials[l2] -
                            factorials[mi - l2] - factorials[j2] -
                            factorials[ki - j2] - lg_theta1[N_i][l2] -
                            lg_theta2[N_i][mi - l2] - lg_theta1[N_i][j2] -
                            lg_theta2[N_i][ki - j2] + lg_theta1[N_i][l2 + j2] +
                            lg_theta2[N_i][mi - l2 + ki - j2] +
                            static_cast<double>(l2) * (logx - log1x) +
                            static_cast<double>(j2) * (logz - log1z));
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
      double sum = 0.0;
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
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending in the
               /// interior of (0,1), but both time increments are below the
               /// threshold
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  vector<int> returnvec;
  vector<double> currsumStore;
  vector<int> mkljStore;
  bool mklj_found = false,
       earlyStop =
           (abs(x - z) <= 0.6);  /// earlyStop gauges whether we can use currsum
  /// as is, or whether we should compute all
  /// probabilities and then sample

  /// Compute denominator
  double eC = 0.0, eCInc = 1.0, eCOldInc = 1.0;
  int Dflip = 1, Djm = 0, Djp = 0;
  int dmode = static_cast<int>(ceil(GriffithsParas(N_i, t).first)), d = dmode;

  while (max(eCOldInc, eCInc) > 0.0 || (!(eC > 0.0))) {
    eCOldInc = eCInc;
    double addon = 0.0;
    for (int f = 0; f != d; f++) {
      double para1 = static_cast<double>(thetaP[N_i][0] + f),
             para2 = static_cast<double>(thetaP[N_i][1] + d - f);
      boost::math::binomial_distribution<double> BIN(d, x);
      boost::math::beta_distribution<double> BETA(para1, para2);

      addon += pdf(BIN, f) * pdf(BETA, z);
    }
    eCInc = QmApprox(N_i, d, t, o) * addon;
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

  vector<int> modeGuess = mkljModeFinder(
      false, N_i, x, z, s, t, o);  /// Get a guess on mode over (m,k,l,j)
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2],
      jMode = modeGuess[3];

  boost::random::uniform_01<double>
      U01;  /// Use this guess & eC to compute a suitable threshold for
  /// subsequent computations
  double currsum = 0.0, u = U01(gen),
         threshold = exp(mkljModeFinder_Evaluator(false, N_i, mMode, kMode,
                                                  lMode, jMode, x, z, s, t, o) -
                         log(eC)) *
                     1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false,
       mUpSwitch = false;  /// Compute m contributions
  double mContr_D = boost::math::lgamma(
             static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + m)),
         mContr_U = mContr_D, mContr, runningMax = -1.0e100;

  while (!mSwitch) {
    double qm = QmApprox(N_i, m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to
                      /// very small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    if (m != mMode) {
      if (mU > mD) {
        mContr_U +=
            log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + (m - 1)));
        mContr = log(qm) + mContr_U;
      } else {
        mContr_D -= log(
            static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] + (m + 1) - 1));
        mContr = log(qm) + mContr_D;
      }
    } else {
      mContr = log(qm) + mContr_U;
    }

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Compute k contributions
    double kContr_D = (boost::math::lgamma(static_cast<double>(
                           thetaP[N_i][0] + thetaP[N_i][1] + k)) -
                       boost::math::lgamma(static_cast<double>(
                           thetaP[N_i][0] + thetaP[N_i][1] + m + k))),
           kContr_U = kContr_D, kContr;

    while (!kSwitch) {
      double qk = QmApprox(N_i, k, t - s, o);
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
          kContr_U += log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              (k - 1))) -
                      log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              m + (k - 1)));
          kContr = log(qk) + kContr_U;
        } else {
          kContr_D += log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              m + (k + 1) - 1)) -
                      log(static_cast<double>(thetaP[N_i][0] + thetaP[N_i][1] +
                                              (k + 1) - 1));
          kContr = log(qk) + kContr_D;
        }
      } else {
        kContr = log(qk) + kContr_U;
      }

      int lFlip = 1, lU = 0, lD = 0, newlMode = min(lMode, m),
          l = newlMode;  /// Redefine lMode in case m is too small!
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;
      boost::math::binomial_distribution<double> BINL(
          m, x);  /// Compute l contributions
      double lContr_D =
                 (log(pdf(BINL, l)) -
                  boost::math::lgamma(static_cast<double>(thetaP[N_i][0] + l)) -
                  boost::math::lgamma(
                      static_cast<double>(thetaP[N_i][1] + m - l))),
             lContr_U = lContr_D, lContr;

      while (!lSwitch) {
        assert((l >= 0) && (l <= m));
        if (l != newlMode) {
          if (lU > lD) {
            lContr_U +=
                log(static_cast<double>(m - (l - 1))) -
                log(static_cast<double>((l - 1) + 1)) + log(x) - log(1.0 - x) +
                log(static_cast<double>(thetaP[N_i][1] + m - (l - 1) - 1)) -
                log(static_cast<double>(thetaP[N_i][0] + (l - 1)));
            lContr = lContr_U;
          } else {
            lContr_D += log(static_cast<double>(l + 1)) -
                        log(static_cast<double>(m - (l + 1) + 1)) +
                        log(1.0 - x) - log(x) +
                        log(static_cast<double>(thetaP[N_i][0] + (l + 1) - 1)) -
                        log(static_cast<double>(thetaP[N_i][1] + m - (l + 1)));
            lContr = lContr_D;
          }
        } else {
          lContr = lContr_U;
        }

        int jFlip = 1, jU = 0, jD = 0, newjMode = min(jMode, k),
            j = newjMode;  /// Redefine jMode in case k is too small!
        bool jSwitch = false, jDownSwitch = false,
             jUpSwitch = false;  /// Compute j contributions

        double jContr_D =
            LogBinomialCoefficientCalculator(k, j) +
            boost::math::lgamma(static_cast<double>(thetaP[N_i][0] + l + j)) -
            boost::math::lgamma(static_cast<double>(thetaP[N_i][0] + j)) +
            boost::math::lgamma(
                static_cast<double>(thetaP[N_i][1] + m - l + k - j)) -
            boost::math::lgamma(static_cast<double>(thetaP[N_i][1] + k - j)) +
            static_cast<double>(thetaP[N_i][0] + j - 1) * log(z) +
            static_cast<double>(thetaP[N_i][1] + k - j - 1) * log(1.0 - z);
        double jContr_U = jContr_D, jContr;

        while (!jSwitch) {
          if (j != newjMode) {
            if (jU > jD) {
              jContr_U +=
                  log(static_cast<double>(k - (j - 1))) -
                  log(static_cast<double>((j - 1) + 1)) + log(z) -
                  log(1.0 - z) +
                  log(static_cast<double>(thetaP[N_i][0] + l + (j - 1))) -
                  log(static_cast<double>(thetaP[N_i][0] + (j - 1))) +
                  log(static_cast<double>(thetaP[N_i][1] + k - (j - 1) - 1)) -
                  log(static_cast<double>(thetaP[N_i][1] + m - l + k - (j - 1) -
                                          1));
              jContr = jContr_U;
            } else {
              jContr_D +=
                  log(static_cast<double>(j + 1)) -
                  log(static_cast<double>(k - (j + 1) + 1)) + log(1.0 - z) -
                  log(z) +
                  log(static_cast<double>(thetaP[N_i][0] + (j + 1) - 1)) -
                  log(static_cast<double>(thetaP[N_i][0] + l + (j + 1) - 1)) +
                  log(static_cast<double>(thetaP[N_i][1] + m - l + k -
                                          (j + 1))) -
                  log(static_cast<double>(thetaP[N_i][1] + k - (j + 1)));
              jContr = jContr_D;
            }
          } else {
            jContr = jContr_U;
          }
          double currsum_inc = exp(mContr + kContr + lContr + jContr - log(eC));
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
    double sum = 0.0;
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

vector<int> WrightFisher::DrawBridgePMFDiffTheta(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending in the
               /// interior of (0,1) with both time increments large enough
               /// but with different thetas either side of the sampling time
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdDiffTheta(N_i, eCvec, 0, x, z, s, t, o));
  vector<int> v_used;
  double Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int G = computeG(N_i, x, z, s, t, o);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(true, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, m, Cs), F2 = computeC(N_i, k, Cts);
    if (o.debug > 2)
      cerr << "(F1,F2,G) = (" << F1 << "," << F2 << "," << G << ")" << endl;
    Fmklj = max(max(F1, F2), G);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmklj << endl;
    }

    while (2 * v < Fmklj || eAU[n] > 1.0 || eBU[n] > 1.0) {
      ++v;
      double newcoefficientU =
          exp(Getlogakm<double>(N_i_s, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) *
                                  (m + 2 * v + theta[N_i_s] - 1) * s / 2.0));
      double newcoefficientL = exp(
          Getlogakm<double>(N_i_s, m + 2 * v + 1, m) +
          static_cast<double>(-(m + 2 * v + 1) *
                              (m + 2 * v + 1 + theta[N_i_s] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU = exp(
          Getlogakm<double>(N_i_t, k + 2 * v, k) +
          static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i_t] - 1) *
                              (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i_t, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i_t] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdDiffTheta(N_i, eCvec, 2 * v + 1, x, z, s, t, o);
        eCU = eCL + GetdDiffTheta(N_i, eCvec, 2 * v + 2, x, z, s, t, o);
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
          return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
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
    precomputeA(N_i_s, N_i_t, m, k);
    double extra_factors =
        factorials[m] + factorials[k] + static_cast<double>(m) * log1x +
        static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
        static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z +
        lg_theta[N_i_s][m] + lg_theta[N_i_t][k] - lg_theta[N_i_s][m + k];
    for (int l = 0; l <= m && !mklj_found; ++l) {
      for (int j = 0; j <= k && !mklj_found; ++j) {
        Amklj = exp(extra_factors - factorials[l] - factorials[m - l] -
                    factorials[j] - factorials[k - j] - lg_theta1[N_i_s][l] -
                    lg_theta2[N_i_s][m - l] - lg_theta1[N_i_t][j] -
                    lg_theta2[N_i_t][k - j] + lg_theta1[N_i_s][l + j] +
                    lg_theta2[N_i_s][m - l + k - j] +
                    static_cast<double>(l) * (logx - log1x) +
                    static_cast<double>(j) * (logz - log1z));
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
          return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
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
          double currsumLold = currsumL, currsumUold = currsumU;
          currsumL = 0.0;
          currsumU = 0.0;

          const vector<double> dUold(dU), dLold(dL);
          for (int i = 0; i <= n; ++i) {
            int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
            ++v_used[i];
            double newcoefficientU = exp(
                       Getlogakm<double>(N_i_s, mi + 2 * v_used[i], mi) +
                       static_cast<double>(
                           -(mi + 2 * v_used[i]) *
                           (mi + 2 * v_used[i] + theta[N_i_s] - 1) * s / 2.0)),
                   newcoefficientL = exp(
                       Getlogakm<double>(N_i_s, mi + 2 * v_used[i] + 1, mi) +
                       static_cast<double>(
                           -(mi + 2 * v_used[i] + 1) *
                           (mi + 2 * v_used[i] + 1 + theta[N_i_s] - 1) * s /
                           2.0));

            eAU[i] = eAL[i] + newcoefficientU;
            eAL[i] = eAU[i] - newcoefficientL;

            newcoefficientU = exp(
                Getlogakm<double>(N_i_t, ki + 2 * v_used[i], ki) +
                static_cast<double>(-(ki + 2 * v_used[i]) *
                                    (ki + 2 * v_used[i] + theta[N_i_t] - 1) *
                                    (t - s) / 2.0));
            newcoefficientL =
                exp(Getlogakm<double>(N_i_t, ki + 2 * v_used[i] + 1, ki) +
                    static_cast<double>(
                        -(ki + 2 * v_used[i] + 1) *
                        (ki + 2 * v_used[i] + 1 + theta[N_i_t] - 1) * (t - s) /
                        2.0));

            eBU[i] = eBL[i] + newcoefficientU;
            eBL[i] = eBU[i] - newcoefficientL;

            if (2 * v_used[i] + 2 > eCindex_computed) {
              assert(2 * v_used[i] == eCindex_computed);
              eCL = eCU -
                    GetdDiffTheta(N_i, eCvec, 2 * v_used[i] + 1, x, z, s, t, o);
              eCU = eCL +
                    GetdDiffTheta(N_i, eCvec, 2 * v_used[i] + 2, x, z, s, t, o);
              eCindex_computed += 2;
            }

            dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                               : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
            dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                         ? static_cast<double>(0.0)
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
              return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
            }
            precomputeA(N_i_s, N_i_t, mi, ki);
            double extra_factors_i =
                factorials[mi] + factorials[ki] +
                static_cast<double>(mi) * log1x +
                static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                static_cast<double>(thetaP[N_i_t][1] + ki - 1.0) * log1z +
                lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                lg_theta[N_i_s][mi + ki];
            int l_upper = (i == n ? l : mi);
            for (int l2 = 0; l2 <= l_upper; ++l2) {
              int j_upper = ((i == n && l2 == l_upper) ? j : ki);
              for (int j2 = 0; j2 <= j_upper; ++j2) {
                Amklj =
                    exp(extra_factors_i - factorials[l2] - factorials[mi - l2] -
                        factorials[j2] - factorials[ki - j2] -
                        lg_theta1[N_i_s][l2] - lg_theta2[N_i_s][mi - l2] -
                        lg_theta1[N_i_t][j2] - lg_theta2[N_i_t][ki - j2] +
                        lg_theta1[N_i_s][l2 + j2] +
                        lg_theta2[N_i_s][mi - l2 + ki - j2] +
                        static_cast<double>(l2) * (logx - log1x) +
                        static_cast<double>(j2) * (logz - log1z));
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
            cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
                 << " = currsumL (n,m,k,l,j) = (" << n << "," << m << "," << k
                 << "," << l << "," << j << ")." << endl;
            exit(1);
          }
          if (currsumUold < currsumU) {
            cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
                 << " = currsumU (n,m,k,l,j) = (" << n << "," << m << "," << k
                 << "," << l << "," << j << ")." << endl;
            exit(1);
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaOneQApprox(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen) {
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  // ofstream outFile;
  // outFile.open("debugging.txt");
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);
  size_t N_i_sorts = (sorts == s) ? N_i_s : N_i_t,
         N_i_sortsapprox = (N_i_sorts == N_i_s) ? N_i_t : N_i_s;

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdDiffThetaOneQApprox(N_i, eCvec, 0, x, z, s, t, o)),
         qApprox = 0.0, runningMax = 1.0e-300;
  vector<int> v_used;
  double Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i_s, s).first,
         vm = GriffithsParas(N_i_s, s).second,
         kmode = GriffithsParas(N_i_t, t - s).first,
         vk = GriffithsParas(N_i_t, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F5;
  if (thetaP[N_i].empty()) {
    F5 = static_cast<int>(
        ceil(((1.0 - x) / (1.0 - z)) + (x / z) * (1.0 + z) * o.eps));
  } else {
    F5 = static_cast<int>(ceil(
        2.0 *
        (max((theta[N_i] / thetaP[N_i][1]) * (1.0 - static_cast<double>(z)),
             (1.0 + theta[N_i]) / (1.0 - static_cast<double>(z))) *
             (1.0 - static_cast<double>(x)) +
         (1.0 + (1.0 / static_cast<double>(z))) *
             (max(
                 ((static_cast<double>(z) * theta[N_i]) + 1.0) / thetaP[N_i][0],
                 1.0)) *
             static_cast<double>(x)) /
        o.eps));
  }
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0)));
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int G = max(max(F3, F4), F5);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(true, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i_sorts, mork, Csorts);

    if (o.debug > 2) cerr << "(F1,G) = (" << F1 << "," << G << ")" << endl;
    Fmklj = max(F1, G);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmklj << endl;
    }

    while (2 * v < Fmklj) {
      ++v;
      double newcoefficientU =
          exp(Getlogakm<double>(N_i_sorts, mork + 2 * v, mork) +
              static_cast<double>(-(mork + 2 * v) *
                                  (mork + 2 * v + theta[N_i_sorts] - 1) *
                                  sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i_sorts, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i_sorts] - 1) *
                                  sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i_sortsapprox, morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL =
            eCU - GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v + 1, x, z, s, t, o);
        eCU =
            eCL + GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v + 2, x, z, s, t, o);
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
          return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double>(0.0)
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
    precomputeA(N_i_s, N_i_t, m, k);
    double extra_factors =
        factorials[m] + factorials[k] + static_cast<double>(m) * log1x +
        static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
        static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z +
        lg_theta[N_i_s][m] + lg_theta[N_i_t][k] - lg_theta[N_i_s][m + k];
    for (int l = 0; l <= m && !mklj_found; ++l) {
      for (int j = 0; j <= k && !mklj_found; ++j) {
        Amklj = exp(extra_factors - factorials[l] - factorials[m - l] -
                    factorials[j] - factorials[k - j] - lg_theta1[N_i_s][l] -
                    lg_theta2[N_i_s][m - l] - lg_theta1[N_i_t][j] -
                    lg_theta2[N_i_t][k - j] + lg_theta1[N_i_s][l + j] +
                    lg_theta2[N_i_s][m - l + k - j] +
                    static_cast<double>(l) * (logx - log1x) +
                    static_cast<double>(j) * (logz - log1z));
        if (o.debug > 2)
          cerr << "Adding to currsums with A(n,m,k,l,j) = A(" << n << "," << m
               << "," << k << "," << l << "," << j << ") = " << endl;

        if (Amklj * dL[n] > 1.0 || Amklj * dU[n] < 0.0 || eAU[n] < 0.0 ||
            eAL[n] > 1.0 || eCU < 0.0) {
          cerr << "Numerical error detected: (m,k,l,j) = " << m << "," << k
               << "," << l << "," << j << "), Amklj = " << Amklj << ", dL[" << n
               << "] = " << dL[n] << ", dU[" << n << "] = " << dU[n] << ", eAU["
               << n << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
               << ", eCU = " << eCU << ", eCL = " << eCL
               << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
               << "," << z << "," << s << "," << t << ") ..." << endl;
          return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
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
          double currsumLold = currsumL, currsumUold = currsumU;
          currsumL = 0.0;
          currsumU = 0.0;

          const vector<double> dUold(dU), dLold(dL);
          for (int i = 0; i <= n; ++i) {
            int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
                morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
            ;
            ++v_used[i];
            double newcoefficientU =
                       exp(Getlogakm<double>(N_i_sorts, morki + 2 * v_used[i],
                                             morki) +
                           static_cast<double>(
                               -(morki + 2 * v_used[i]) *
                               (morki + 2 * v_used[i] + theta[N_i_sorts] - 1) *
                               sorts / 2.0)),
                   newcoefficientL =
                       exp(Getlogakm<double>(N_i_sorts,
                                             morki + 2 * v_used[i] + 1, morki) +
                           static_cast<double>(-(morki + 2 * v_used[i] + 1) *
                                               (morki + 2 * v_used[i] + 1 +
                                                theta[N_i_sorts] - 1) *
                                               sorts / 2.0));

            eAU[i] = eAL[i] + newcoefficientU;
            eAL[i] = eAU[i] - newcoefficientL;

            qApprox =
                DiscretisedNormCDF(N_i_sortsapprox, morkapproxi, sortsapprox);

            if (2 * v_used[i] + 2 > eCindex_computed) {
              assert(2 * v_used[i] == eCindex_computed);
              eCL = eCU - GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v_used[i] + 1,
                                                  x, z, s, t, o);
              eCU = eCL + GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v_used[i] + 2,
                                                  x, z, s, t, o);
              eCindex_computed += 2;
            }

            dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                               : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
            dL[i] = (eAL[i] < 0.0 ? static_cast<double>(0.0)
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
              return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
            }
            double extra_factors_i =
                factorials[mi] + factorials[ki] +
                static_cast<double>(mi) * log1x +
                static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                static_cast<double>(thetaP[N_i_t][1] + ki - 1.0) * log1z +
                lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                lg_theta[N_i_s][mi + ki];
            int l_upper = (i == n ? l : mi);
            for (int l2 = 0; l2 <= l_upper; ++l2) {
              int j_upper = ((i == n && l2 == l_upper) ? j : ki);
              for (int j2 = 0; j2 <= j_upper; ++j2) {
                Amklj =
                    exp(extra_factors_i - factorials[l2] - factorials[mi - l2] -
                        factorials[j2] - factorials[ki - j2] -
                        lg_theta1[N_i_s][l2] - lg_theta2[N_i_s][mi - l2] -
                        lg_theta1[N_i_t][j2] - lg_theta2[N_i_t][ki - j2] +
                        lg_theta1[N_i_s][l2 + j2] +
                        lg_theta2[N_i_s][mi - l2 + ki - j2] +
                        static_cast<double>(l2) * (logx - log1x) +
                        static_cast<double>(j2) * (logz - log1z));
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
            cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
                 << " = currsumL (n,m,k,l,j) = (" << n << "," << m << "," << k
                 << "," << l << "," << j << ")." << endl;
            exit(1);
          }
          if (currsumUold < currsumU) {
            cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
                 << " = currsumU (n,m,k,l,j) = (" << n << "," << m << "," << k
                 << "," << l << "," << j << ")." << endl;
            exit(1);
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
      double sum = 0.0;
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaApprox(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937
        &gen)  /// Draws from the law of a bridge starting and ending in the
               /// interior of (0,1), but both time increments are below the
               /// threshold and the mutation parameters differ either side of
               /// s
{
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  vector<int> returnvec;
  vector<double> currsumStore;
  vector<int> mkljStore;
  bool mklj_found = false,
       earlyStop =
           (abs(x - z) <= 0.6);  /// earlyStop gauges whether we can use currsum
  /// as is, or whether we should compute all
  /// probabilities and then sample

  /// Compute denominator
  double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-50;
  int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
  int dmmode = static_cast<int>(ceil(GriffithsParas(N_i_s, s).first)),
      dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
  bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
  int dkmode = static_cast<int>(ceil(GriffithsParas(N_i_t, t - s).first));

  while (!dmSwitch) {
    emCOldInc = emCInc;
    double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
    int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
    bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
    while (!dkSwitch) {
      ekCOldInc = ekCInc;
      ekCInc = QmApprox(N_i_t, dk, t - s, o) *
               calculate_expectation(N_i, dm, dk, x, z);
      ekC += ekCInc;
      if (!(dkDownSwitch))  /// Switching mechanism for j
      {
        if (sgn(dk - dkmode) <= 0) {
          dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
        }
      }

      if (!(dkUpSwitch)) {
        if (sgn(dk - dkmode) >= 0) {
          dkUpSwitch = (ekCInc < eC_threshold);
        }
      }

      dkSwitch = (dkDownSwitch && dkUpSwitch);
      if (!dkSwitch) {
        if (dkFlip == 1) {
          dkU++;
          dk = dkmode + dkU;
          dkFlip *= (dkDownSwitch ? 1 : -1);
        } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
          dkD++;
          dk = dkmode - dkD;
          dkFlip *= (dkUpSwitch ? 1 : -1);
        }
      }
    }
    emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
    eC += emCInc;

    if (!(dmDownSwitch))  /// Switching mechanism for j
    {
      if (sgn(dm - dmmode) <= 0) {
        dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
      }
    }

    if (!(dmUpSwitch)) {
      if (sgn(dm - dmmode) >= 0) {
        dmUpSwitch = (emCInc < eC_threshold);
      }
    }

    dmSwitch = (dmDownSwitch && dmUpSwitch);
    if (!dmSwitch) {
      if (dmFlip == 1) {
        dmU++;
        dm = dmmode + dmU;
        dmFlip *= (dmDownSwitch ? 1 : -1);
      } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
        dmD++;
        dm = dmmode - dmD;
        dmFlip *= (dmUpSwitch ? 1 : -1);
      }
    }
  }

  vector<int> modeGuess = mkljModeFinder(
      true, N_i, x, z, s, t, o);  /// Get a guess on mode over (m,k,l,j)
  int mMode = modeGuess[0], kMode = modeGuess[1], lMode = modeGuess[2],
      jMode = modeGuess[3];

  boost::random::uniform_01<double>
      U01;  /// Use this guess & eC to compute a suitable threshold for
  /// subsequent computations
  double currsum = 0.0, u = U01(gen),
         threshold = exp(mkljModeFinder_Evaluator(true, N_i, mMode, kMode,
                                                  lMode, jMode, x, z, s, t, o) -
                         log(eC)) *
                     1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false,
       mUpSwitch = false;  /// Compute m contributions
  double mContr, runningMax = -1.0e100;
  double constant = static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                    static_cast<double>(thetaP[N_i_t][1] - 1.0) * log1z;

  while (!mSwitch) {
    precomputeA(N_i_s, N_i_t, m, kMode);
    double qm = QmApprox(N_i_s, m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to
                      /// very small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    mContr =
        factorials[m] + lg_theta[N_i_s][m] + static_cast<double>(m) * log1x;
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Compute k contributions
    double kContr;

    while (!kSwitch) {
      precomputeA(N_i_s, N_i_t, m, k);
      double qk = QmApprox(N_i_t, k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      if (!(qk > 0.0)) {
        break;
      }
      kContr = factorials[k] + lg_theta[N_i_t][k] - lg_theta[N_i_s][m + k] +
               static_cast<double>(k) * log1z;
      kContr += log(qk);

      int lFlip = 1, lU = 0, lD = 0, newlMode = min(lMode, m),
          l = newlMode;  /// Redefine lMode in case m is too small!
      bool lSwitch = false, lDownSwitch = false, lUpSwitch = false;
      /// Compute l contributions
      double lContr;

      while (!lSwitch) {
        assert((l >= 0) && (l <= m));
        lContr = -factorials[l] - factorials[m - l] +
                 static_cast<double>(l) * (logx - log1x) - lg_theta1[N_i_s][l] -
                 lg_theta2[N_i_s][m - l];

        int jFlip = 1, jU = 0, jD = 0, newjMode = min(jMode, k),
            j = newjMode;  /// Redefine jMode in case k is too small!
        bool jSwitch = false, jDownSwitch = false,
             jUpSwitch = false;  /// Compute j contributions

        double jContr;

        while (!jSwitch) {
          jContr = -factorials[j] - factorials[k - j] +
                   lg_theta1[N_i_s][l + j] - lg_theta1[N_i_t][j] +
                   lg_theta2[N_i_s][m - l + k - j] - lg_theta2[N_i_t][k - j] +
                   static_cast<double>(j) * (logz - log1z);
          double currsum_inc =
              exp(mContr + kContr + lContr + jContr + constant - log(eC));
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
    double sum = 0.0;
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaBoundaries(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen) {
  assert((!(x > 0.0) || !(x < 1.0)) && (t > s) && (s > 0.0));
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double> U01;
  vector<double> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double u = U01(gen), currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdDiffThetaBoundaries(N_i, eCvec, 0, x, z, s, t, o));
  vector<int> v_used;
  map<vector<int>, double> Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0;

  /// Compute F ensuring theoretical convergence of lower and upper sums - F1
  /// and F2 depend on m and k but F3-5 do not. None depend on j.
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int F3 = max(computeE(N_i_s, Cs), computeE(N_i_t, Cts)),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(s) - (theta[N_i_s] + 1.0) / 2.0))),
      F5 = static_cast<int>(ceil(max(
          0.0, 1.0 / static_cast<double>(t - s) - (theta[N_i_t] + 1.0) / 2.0))),
      F6 = computeGBoundaries(N_i, x, z, s, t, o);
  while ((theta[N_i_s] + 2 * F4 + 1) *
             exp(-(2 * F4 + theta[N_i_s]) * s / 2.0) >=
         1 - o.eps)
    ++F4;
  while ((theta[N_i_t] + 2 * F5 + 1) *
             exp(-(2 * F5 + theta[N_i_t]) * (t - s) / 2.0) >=
         1 - o.eps)
    ++F5;
  int F7 = max(max(max(F3, F4), F5), F6);

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(true, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i_s, m, Cs), F2 = computeC(N_i_t, k, Cts);
    if (o.debug > 2)
      cerr << "(F1,F2,F3,F4,F5) = (" << F1 << "," << F2 << "," << F3 << ","
           << F4 << "," << F5 << ")" << endl;
    Fmkj = max(max(F1, F2), F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << s << "," << t << ") = " << Fmkj << endl;
    }

    while (2 * v < Fmkj || eAU[n] > 1.0 || eBU[n] > 1.0 || eCL < 0.0) {
      ++v;
      double newcoefficientU =
          exp(Getlogakm<double>(N_i_s, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) *
                                  (m + 2 * v + theta[N_i_s] - 1) * s / 2.0));
      double newcoefficientL = exp(
          Getlogakm<double>(N_i_s, m + 2 * v + 1, m) +
          static_cast<double>(-(m + 2 * v + 1) *
                              (m + 2 * v + 1 + theta[N_i_s] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU = exp(
          Getlogakm<double>(N_i_t, k + 2 * v, k) +
          static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i_t] - 1) *
                              (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i_t, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i_t] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL =
            eCU - GetdDiffThetaBoundaries(N_i, eCvec, 2 * v + 1, x, z, s, t, o);
        eCU =
            eCL + GetdDiffThetaBoundaries(N_i, eCvec, 2 * v + 2, x, z, s, t, o);
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
          return DrawBridgePMFDiffThetaBoundariesApprox(N_i, x, z, s, t, o,
                                                        gen);
        }

        break;
      }
    }

    /// Compute the appropriate additional terms

    precomputeA(N_i_s, N_i_t, m, k);
    double beta1, beta2, beta3;
    if (!(x > 0.0)) {
      if (!(z > 0.0)) {  // x = 0 && z = 0
        beta1 = lg_theta1[N_i_s][0] + lg_theta2[N_i_s][m + k] -
                lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k];
      } else {  // x = 0 && z = 1
        beta1 =
            lg_theta1[N_i_s][k] + lg_theta2[N_i_s][m] - lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0];
      }
    } else {
      if (!(z > 0.0)) {  // x = 1 && z = 0
        beta1 =
            lg_theta1[N_i_s][m] + lg_theta2[N_i_s][k] - lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k];
      } else {  // x = 1 && z = 1
        beta1 = lg_theta1[N_i_s][k + m] + lg_theta2[N_i_s][0] -
                lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0];
      }
    }
    double addon = beta1 + beta2 + beta3;

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) + addon - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        ++v_used[i];
        double newcoefficientU =
                   exp(Getlogakm<double>(N_i_s, mi + 2 * v_used[i], mi) +
                       static_cast<double>(
                           -(mi + 2 * v_used[i]) *
                           (mi + 2 * v_used[i] + theta[N_i_s] - 1) * s / 2.0)),
               newcoefficientL = exp(
                   Getlogakm<double>(N_i_s, mi + 2 * v_used[i] + 1, mi) +
                   static_cast<double>(
                       -(mi + 2 * v_used[i] + 1) *
                       (mi + 2 * v_used[i] + 1 + theta[N_i_s] - 1) * s / 2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        newcoefficientU =
            exp(Getlogakm<double>(N_i_t, ki + 2 * v_used[i], ki) +
                static_cast<double>(-(ki + 2 * v_used[i]) *
                                    (ki + 2 * v_used[i] + theta[N_i_t] - 1) *
                                    (t - s) / 2.0));
        newcoefficientL = exp(
            Getlogakm<double>(N_i_t, ki + 2 * v_used[i] + 1, ki) +
            static_cast<double>(-(ki + 2 * v_used[i] + 1) *
                                (ki + 2 * v_used[i] + 1 + theta[N_i_t] - 1) *
                                (t - s) / 2.0));

        eBU[i] = eBL[i] + newcoefficientU;
        eBL[i] = eBU[i] - newcoefficientL;

        if (2 * (v_used[i] + 1) > eCindex_computed) {
          assert(2 * v_used[i] == eCindex_computed);
          eCL = eCU - GetdDiffThetaBoundaries(N_i, eCvec, 2 * v_used[i] + 1, x,
                                              z, s, t, o);
          eCU = eCL + GetdDiffThetaBoundaries(N_i, eCvec, 2 * v_used[i] + 2, x,
                                              z, s, t, o);
          eCindex_computed += 2;
        }

        if (!(x > 0.0)) {
          if (!(z > 0.0)) {  // x = 0 && z = 0
            beta1 = lg_theta1[N_i_s][0] + lg_theta2[N_i_s][mi + ki] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][0] -
                    lg_theta2[N_i_s][mi];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][0] -
                    lg_theta2[N_i_t][ki];
          } else {  // x = 0 && z = 1
            beta1 = lg_theta1[N_i_s][ki] + lg_theta2[N_i_s][mi] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][0] -
                    lg_theta2[N_i_s][mi];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][ki] -
                    lg_theta2[N_i_t][0];
          }
        } else {
          if (!(z > 0.0)) {  // x = 1 && z = 0
            beta1 = lg_theta1[N_i_s][mi] + lg_theta2[N_i_s][ki] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][mi] -
                    lg_theta2[N_i_s][0];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][0] -
                    lg_theta2[N_i_t][ki];
          } else {  // x = 1 && z = 1
            beta1 = lg_theta1[N_i_s][ki + mi] + lg_theta2[N_i_s][0] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][mi] -
                    lg_theta2[N_i_s][0];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][ki] -
                    lg_theta2[N_i_t][0];
          }
        }
        addon = beta1 + beta2 + beta3;

        dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                           : exp(log(eAU[i]) + log(eBU[i]) + addon - log(eCL)));
        dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                     ? static_cast<double>(0.0)
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
          return DrawBridgePMFDiffThetaBoundariesApprox(N_i, x, z, s, t, o,
                                                        gen);
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
    if (!(z > 0.0)) {
      curr_mk[n][2] = 0;  /// Setting l & j to be 0
      curr_mk[n][3] = 0;
    } else {
      curr_mk[n][2] = 0;  /// Setting l to be 0, j to be k
      curr_mk[n][3] = curr_mk[n][1];
    }
  } else {
    if (!(z > 0.0)) {
      curr_mk[n][2] = curr_mk[n][0];  /// Setting l to be m, j to be 0
      curr_mk[n][3] = 0;
    } else {
      curr_mk[n][2] = curr_mk[n][0];  /// Setting l & j to be m & k respectively
      curr_mk[n][3] = curr_mk[n][1];
    }
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaBoundariesOneQApprox(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen) {
  assert((!(x > 0.0) || !(x < 1.0)) && (!(z > 0.0) || !(z < 1.0)) && (t > s) &&
         (s > 0.0));
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);
  size_t N_i_sorts = (sorts == s) ? N_i : N_i + 1,
         N_i_sortsapprox = (N_i_sorts == N_i) ? N_i + 1 : N_i, N_i_s = N_i,
         N_i_t = N_i + 1;

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mk_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdDiffThetaOneQApprox(N_i, eCvec, 0, x, z, s, t, o)),
         qApprox = 0.0, runningMax = 1.0e-300;
  vector<int> v_used;
  map<vector<int>, double> Amkj;
  int n = -1, Fmkj = 0, eCindex_computed = 0;

  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i_s, s).first,
         vm = GriffithsParas(N_i_s, s).second;
  double kmode = GriffithsParas(N_i_t, t - s).first,
         vk = GriffithsParas(N_i_t, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring theoretical convergence of upper and lower bounds
  pair<vector<int>, double> Csorts;
  Csorts.second = sorts;
  int F3 = computeE(N_i_s, Csorts),
      F4 = static_cast<int>(ceil(max(0.0, 1.0 / static_cast<double>(sorts) -
                                              (theta[N_i_sorts] + 1.0) / 2.0))),
      F5 = computeGBoundaries(N_i, x, z, s, t, o);
  while ((theta[N_i_sorts] + 2 * F4 + 1) *
             exp(-(2 * F4 + theta[N_i_sorts]) * sorts / 2.0) >=
         1 - o.eps)
    ++F4;
  int F6 = max(max(F3, F4), F5);

  while (!mk_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(true, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, mork, Csorts);
    if (o.debug > 2)
      cerr << "(F1,F3,F4,F5) = (" << F1 << "," << F3 << "," << F4 << "," << F5
           << ")" << endl;
    Fmkj = max(F1, F6);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << s << "," << t << ") = " << Fmkj << endl;
    }

    while (2 * v < Fmkj || eAU[n] > 1.0 || eCL < 0.0) {
      ++v;
      double newcoefficientU =
          exp(Getlogakm<double>(N_i_sorts, mork + 2 * v, mork) +
              static_cast<double>(-(mork + 2 * v) *
                                  (mork + 2 * v + theta[N_i_sorts] - 1) *
                                  sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i_sorts, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i_sorts] - 1) *
                                  sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i_sortsapprox, morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL =
            eCU - GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v + 1, x, z, s, t, o);
        eCU =
            eCL + GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v + 2, x, z, s, t, o);
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
          return DrawBridgePMFDiffThetaBoundariesApprox(N_i, x, z, s, t, o,
                                                        gen);
        }

        break;
      }
    }

    /// Computing the appropriate additional contributions

    precomputeA(N_i_s, N_i_t, m, k);
    double beta1, beta2, beta3;
    if (!(x > 0.0)) {
      if (!(z > 0.0)) {  // x = 0 && z = 0
        beta1 = lg_theta1[N_i_s][0] + lg_theta2[N_i_s][m + k] -
                lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k];
      } else {  // x = 0 && z = 1
        beta1 =
            lg_theta1[N_i_s][k] + lg_theta2[N_i_s][m] - lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0];
      }
    } else {
      if (!(z > 0.0)) {  // x = 1 && z = 0
        beta1 =
            lg_theta1[N_i_s][m] + lg_theta2[N_i_s][k] - lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k];
      } else {  // x = 1 && z = 1
        beta1 = lg_theta1[N_i_s][k + m] + lg_theta2[N_i_s][0] -
                lg_theta[N_i_s][m + k];
        beta2 = lg_theta[N_i_s][m] - lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0];
        beta3 = lg_theta[N_i_t][k] - lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0];
      }
    }
    double addon = beta1 + beta2 + beta3;

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) + addon - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double>(0.0)
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
    mkljStore[3 + (4 * (currsumStore.size() - 1))] = (!(z > 0.0) ? 0 : k);

    if (counter == totalpts)  /// Gaussian approximation leads to currsum
                              /// summing to < 1.0, so we renormalise and sample
    {
      LogSumExp(currsumStore, runningMax);
      double sum = 0.0;
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
      double currsumLold = currsumL, currsumUold = currsumU;
      currsumL = 0.0;
      currsumU = 0.0;

      const vector<double> dUold(dU), dLold(dL);
      for (int i = 0; i <= n; ++i) {
        int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
        int morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
        ++v_used[i];
        double newcoefficientU = exp(
            Getlogakm<double>(N_i_sorts, morki + 2 * v_used[i], morki) +
            static_cast<double>(-(morki + 2 * v_used[i]) *
                                (morki + 2 * v_used[i] + theta[N_i_sorts] - 1) *
                                sorts / 2.0));
        double newcoefficientL =
            exp(Getlogakm<double>(N_i_sorts, morki + 2 * v_used[i] + 1, morki) +
                static_cast<double>(
                    -(morki + 2 * v_used[i] + 1) *
                    (morki + 2 * v_used[i] + 1 + theta[N_i_sorts] - 1) * sorts /
                    2.0));

        eAU[i] = eAL[i] + newcoefficientU;
        eAL[i] = eAU[i] - newcoefficientL;

        qApprox = DiscretisedNormCDF(N_i_sortsapprox, morkapproxi, sortsapprox);

        if (2 * (v_used[i] + 1) > eCindex_computed) {
          assert(2 * v_used[i] == eCindex_computed);
          eCL = eCU - GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v_used[i] + 1, x,
                                              z, s, t, o);
          eCU = eCL + GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v_used[i] + 2, x,
                                              z, s, t, o);
          eCindex_computed += 2;
        }

        if (!(x > 0.0)) {
          if (!(z > 0.0)) {  // x = 0 && z = 0
            beta1 = lg_theta1[N_i_s][0] + lg_theta2[N_i_s][mi + ki] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][0] -
                    lg_theta2[N_i_s][mi];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][0] -
                    lg_theta2[N_i_t][ki];
          } else {  // x = 0 && z = 1
            beta1 = lg_theta1[N_i_s][ki] + lg_theta2[N_i_s][mi] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][0] -
                    lg_theta2[N_i_s][mi];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][ki] -
                    lg_theta2[N_i_t][0];
          }
        } else {
          if (!(z > 0.0)) {  // x = 1 && z = 0
            beta1 = lg_theta1[N_i_s][mi] + lg_theta2[N_i_s][ki] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][mi] -
                    lg_theta2[N_i_s][0];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][0] -
                    lg_theta2[N_i_t][ki];
          } else {  // x = 1 && z = 1
            beta1 = lg_theta1[N_i_s][ki + mi] + lg_theta2[N_i_s][0] -
                    lg_theta[N_i_s][mi + ki];
            beta2 = lg_theta[N_i_s][mi] - lg_theta1[N_i_s][mi] -
                    lg_theta2[N_i_s][0];
            beta3 = lg_theta[N_i_t][ki] - lg_theta1[N_i_t][ki] -
                    lg_theta2[N_i_t][0];
          }
        }
        addon = beta1 + beta2 + beta3;

        dU[i] =
            (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[i]) + log(qApprox) + addon - log(eCL)));
        dL[i] =
            (eAL[i] < 0.0 ? static_cast<double>(0.0)
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
          return DrawBridgePMFDiffThetaBoundariesApprox(N_i, x, z, s, t, o,
                                                        gen);
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
    if (!(z > 0.0)) {
      curr_mk[n][2] = 0;  /// Setting l & j to be 0
      curr_mk[n][3] = 0;
    } else {
      curr_mk[n][2] = 0;  /// Setting l to be 0, j to be k
      curr_mk[n][3] = curr_mk[n][1];
    }
  } else {
    if (!(z > 0.0)) {
      curr_mk[n][2] = curr_mk[n][0];  /// Setting l to be m, j to be 0
      curr_mk[n][3] = 0;
    } else {
      curr_mk[n][2] = curr_mk[n][0];  /// Setting l & j to be m & k respectively
      curr_mk[n][3] = curr_mk[n][1];
    }
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaBoundariesApprox(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen) {
  assert((!(x > 0.0) || !(x < 1.0)) && (!(z > 0.0) || !(z < 1.0)) && (t > s) &&
         (s > 0.0));
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  bool mklj_found = false;
  vector<int> returnvec;
  vector<double> currsumStore;
  vector<int> mkljStore;
  /// Compute denominator
  double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-50;
  int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
  int dmmode = static_cast<int>(ceil(GriffithsParas(N_i_s, s).first)),
      dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
  bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
  int dkmode = static_cast<int>(ceil(GriffithsParas(N_i_t, t - s).first));

  while (!dmSwitch) {
    emCOldInc = emCInc;
    double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
    int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
    bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
    while (!dkSwitch) {
      ekCOldInc = ekCInc;
      ekCInc = QmApprox(N_i_t, dk, t - s, o) *
               calculate_expectation(N_i, dm, dk, x, z);
      ekC += ekCInc;
      if (!(dkDownSwitch))  /// Switching mechanism for j
      {
        if (sgn(dk - dkmode) <= 0) {
          dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
        }
      }

      if (!(dkUpSwitch)) {
        if (sgn(dk - dkmode) >= 0) {
          dkUpSwitch = (ekCInc < eC_threshold);
        }
      }

      dkSwitch = (dkDownSwitch && dkUpSwitch);
      if (!dkSwitch) {
        if (dkFlip == 1) {
          dkU++;
          dk = dkmode + dkU;
          dkFlip *= (dkDownSwitch ? 1 : -1);
        } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
          dkD++;
          dk = dkmode - dkD;
          dkFlip *= (dkUpSwitch ? 1 : -1);
        }
      }
    }
    emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
    eC += emCInc;

    if (!(dmDownSwitch))  /// Switching mechanism for j
    {
      if (sgn(dm - dmmode) <= 0) {
        dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
      }
    }

    if (!(dmUpSwitch)) {
      if (sgn(dm - dmmode) >= 0) {
        dmUpSwitch = (emCInc < eC_threshold);
      }
    }

    dmSwitch = (dmDownSwitch && dmUpSwitch);
    if (!dmSwitch) {
      if (dmFlip == 1) {
        dmU++;
        dm = dmmode + dmU;
        dmFlip *= (dmDownSwitch ? 1 : -1);
      } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
        dmD++;
        dm = dmmode - dmD;
        dmFlip *= (dmUpSwitch ? 1 : -1);
      }
    }
  }

  vector<int> modeGuess = mkModeFinder(
      true, N_i, x, z, s, t, o);  /// Get a guess on mode over (m,k,l,j)
  int mMode = modeGuess[0], kMode = modeGuess[1];

  boost::random::uniform_01<double>
      U01;  /// Use this guess & eC to compute a suitable threshold for
  /// subsequent computations
  double currsum = 0.0, u = U01(gen),
         threshold = exp(mkModeFinder_Evaluator(true, N_i, mMode, kMode, x, z,
                                                s, t, o) -
                         log(eC)) *
                     1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false,
       mUpSwitch = false;  /// Compute m contributions
  double mContr, runningMax = -1.0e100;
  double constant =
      (!(x > 0.0) ? ((!(z > 0.0)) ? -lg_theta1[N_i_t][0]
                                  : -lg_theta1[N_i_s][0] - lg_theta2[N_i_t][0])
                  : ((!(z > 0.0)) ? -lg_theta1[N_i_t][0] - lg_theta2[N_i_s][0]
                                  : -lg_theta1[N_i_t][0]));

  while (!mSwitch) {
    double qm = QmApprox(N_i_s, m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to
                      /// very small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    mContr = lg_theta[N_i_s][m] + ((x == z) ? -lg_theta2[N_i_s][m] : 0.0);
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Compute k contributions
    double kContr;

    while (!kSwitch) {
      double qk = QmApprox(N_i_t, k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      if (!(qk > 0.0)) {
        break;
      }
      kContr = lg_theta[N_i_t][k] - lg_theta[N_i_s][m + k] +
               ((x == z) ? lg_theta2[N_i_s][m + k] - lg_theta2[N_i_t][k]
                         : ((!(x > 0.0))
                                ? lg_theta1[N_i_s][k] - lg_theta1[N_i_t][k]
                                : lg_theta2[N_i_s][k] - lg_theta2[N_i_t][k]));
      kContr += log(qk);

      double currsum_inc = exp(mContr + kContr + constant - log(eC));
      runningMax =
          max(currsum_inc,
              runningMax);  /// Running max needed for log-sum-exp trick later
      currsum += currsum_inc;
      currsumStore.push_back(currsum);

      mkljStore.resize(mkljStore.size() + 4, -1);
      mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
      mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
      mkljStore[2 + (4 * (currsumStore.size() - 1))] = ((!(x > 0.0)) ? 0 : m);
      mkljStore[3 + (4 * (currsumStore.size() - 1))] = ((!(z > 0.0)) ? 0 : k);

      if ((currsum > u))  /// if earlyStop is allowed, we can
                          /// stop once currsum exceeds u
      {
        returnvec.push_back(m);
        returnvec.push_back(k);
        returnvec.push_back(((!(x > 0.0)) ? 0 : m));
        returnvec.push_back(((!(z > 0.0)) ? 0 : k));

        mklj_found = true;
        goto End;
      }

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch = ((currsum_inc < threshold) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = (currsum_inc < threshold);
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
    double sum = 0.0;
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaInterior(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen) {
  assert(((!(x > 0.0)) || (!(x < 1.0)) || (!(z > 0.0)) || (!(z < 1.0))) &&
         (t > s) && (s > 0.0));
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0);
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, eBL, eBU, dL, dU;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdDiffTheta(N_i, eCvec, 0, x, z, s, t, o));
  vector<int> v_used;
  double Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double> Cs, Cts, Ct;
  Cs.second = s;
  Cts.second = t - s;
  Ct.second = t;
  int G = computeG(N_i, x, z, s, t, o);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(true, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i, m, Cs), F2 = computeC(N_i, k, Cts);
    if (o.debug > 2)
      cerr << "(F1,F2,G) = (" << F1 << "," << F2 << "," << G << ")" << endl;
    Fmklj = max(max(F1, F2), G);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmklj << endl;
    }

    while (2 * v < Fmklj || eAU[n] > 1.0 || eBU[n] > 1.0) {
      ++v;
      double newcoefficientU =
          exp(Getlogakm<double>(N_i_s, m + 2 * v, m) +
              static_cast<double>(-(m + 2 * v) *
                                  (m + 2 * v + theta[N_i_s] - 1) * s / 2.0));
      double newcoefficientL = exp(
          Getlogakm<double>(N_i_s, m + 2 * v + 1, m) +
          static_cast<double>(-(m + 2 * v + 1) *
                              (m + 2 * v + 1 + theta[N_i_s] - 1) * s / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      newcoefficientU = exp(
          Getlogakm<double>(N_i_t, k + 2 * v, k) +
          static_cast<double>(-(k + 2 * v) * (k + 2 * v + theta[N_i_t] - 1) *
                              (t - s) / 2.0));
      newcoefficientL =
          exp(Getlogakm<double>(N_i_t, k + 2 * v + 1, k) +
              static_cast<double>(-(k + 2 * v + 1) *
                                  (k + 2 * v + 1 + theta[N_i_t] - 1) * (t - s) /
                                  2.0));

      eBU[n] = eBL[n] + newcoefficientU;
      eBL[n] = eBU[n] - newcoefficientL;

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL = eCU - GetdDiffTheta(N_i, eCvec, 2 * v + 1, x, z, s, t, o);
        eCU = eCL + GetdDiffTheta(N_i, eCvec, 2 * v + 2, x, z, s, t, o);
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
          return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(eBU[n]) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 || eBL[n] < 0.0
                 ? static_cast<double>(0.0)
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
    precomputeA(N_i_s, N_i_t, m, k);
    double constant_factors =
        (!(x > 0.0))
            ? (factorials[k] - lg_theta[N_i_s][m + k] + lg_theta[N_i_s][m] +
               lg_theta[N_i_t][k] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m] +
               static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
               static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z)
            : ((!(x < 1.0))
                   ? (factorials[k] - lg_theta[N_i_s][m + k] +
                      lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                      lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0] +
                      static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                      static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z)
                   : ((!(z > 0.0))
                          ? (factorials[m] - lg_theta[N_i_s][m + k] +
                             lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                             lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k] +
                             static_cast<double>(m) * log1x)
                          : (factorials[m] - lg_theta[N_i_s][m + k] +
                             lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                             lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0] +
                             static_cast<double>(m) * log1x)));
    int upper_bound = ((!(x > 0.0)) || (!(x < 1.0))) ? k : m;
    for (int index = 0; index <= upper_bound && !mklj_found; ++index) {
      double var_stuff =
          (!(x > 0.0))
              ? (-factorials[index] - factorials[k - index] +
                 lg_theta1[N_i_s][index] + lg_theta2[N_i_s][m + k - index] -
                 lg_theta1[N_i_t][index] - lg_theta2[N_i_t][k - index] +
                 static_cast<double>(index) * (logz - log1z))
              : ((!(x < 1.0))
                     ? (-factorials[index] - factorials[k - index] +
                        lg_theta1[N_i_s][m + index] +
                        lg_theta2[N_i_s][k - index] - lg_theta1[N_i_t][index] -
                        lg_theta2[N_i_t][k - index] +
                        static_cast<double>(index) * (logz - log1z))
                     : ((!(z > 0.0))
                            ? (-factorials[index] - factorials[m - index] +
                               lg_theta1[N_i_s][index] +
                               lg_theta2[N_i_s][m - index + k] -
                               lg_theta1[N_i_s][index] -
                               lg_theta2[N_i_s][m - index] +
                               static_cast<double>(index) * (logx - log1x))
                            : (-factorials[index] - factorials[m - index] +
                               lg_theta1[N_i_s][index + k] +
                               lg_theta2[N_i_s][m - index] -
                               lg_theta1[N_i_s][index] -
                               lg_theta2[N_i_s][m - index] +
                               static_cast<double>(index) * (logx - log1x))));
      Amklj = exp(constant_factors + var_stuff);
      if (o.debug > 2)
        cerr << "Adding to currsums with A(n,m,k,index) = A(" << n << "," << m
             << "," << k << "," << index << ") = " << Amklj << endl;

      if (Amklj * dL[n] > 1.0 || Amklj * dU[n] < 0.0 || eAU[n] < 0.0 ||
          eBU[n] < 0.0 || eAL[n] > 1.0 || eBL[n] > 1.0 || eCU < 0.0) {
        cerr << "Numerical error detected: (m,k,index) = " << m << "," << k
             << "," << index << "), Amklj = " << Amklj << ", dL[" << n
             << "] = " << dL[n] << ", dU[" << n << "] = " << dU[n] << ", eAU["
             << n << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
             << ",  eBU[" << n << "] = " << eBU[n] << ", eBL[" << n
             << "] = " << eBL[n] << ", eCU = " << eCU << ", eCL = " << eCL
             << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
             << "," << z << "," << s << "," << t << ") ..." << endl;
        return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
      }

      currsumL += Amklj * dL[n];
      currsumU += Amklj * dU[n];
      bool decision_on_mklj_made = (currsumL > u || currsumU < u);

      if (currsumL > currsumU) {
        cerr << "Error: currsumL = " << currsumL << " > " << currsumU
             << " = currsumU (n,m,k,index) = (" << n << "," << m << "," << k
             << "," << index << ")." << endl;
        exit(1);
      }

      while (!decision_on_mklj_made)  /// Refine upper and lower bounds
      {
        double currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1];
          ++v_used[i];
          double newcoefficientU = exp(
                     Getlogakm<double>(N_i_s, mi + 2 * v_used[i], mi) +
                     static_cast<double>(
                         -(mi + 2 * v_used[i]) *
                         (mi + 2 * v_used[i] + theta[N_i_s] - 1) * s / 2.0)),
                 newcoefficientL =
                     exp(Getlogakm<double>(N_i_s, mi + 2 * v_used[i] + 1, mi) +
                         static_cast<double>(
                             -(mi + 2 * v_used[i] + 1) *
                             (mi + 2 * v_used[i] + 1 + theta[N_i_s] - 1) * s /
                             2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          newcoefficientU =
              exp(Getlogakm<double>(N_i_t, ki + 2 * v_used[i], ki) +
                  static_cast<double>(-(ki + 2 * v_used[i]) *
                                      (ki + 2 * v_used[i] + theta[N_i_t] - 1) *
                                      (t - s) / 2.0));
          newcoefficientL = exp(
              Getlogakm<double>(N_i_t, ki + 2 * v_used[i] + 1, ki) +
              static_cast<double>(-(ki + 2 * v_used[i] + 1) *
                                  (ki + 2 * v_used[i] + 1 + theta[N_i_t] - 1) *
                                  (t - s) / 2.0));

          eBU[i] = eBL[i] + newcoefficientU;
          eBL[i] = eBU[i] - newcoefficientL;

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU -
                  GetdDiffTheta(N_i, eCvec, 2 * v_used[i] + 1, x, z, s, t, o);
            eCU = eCL +
                  GetdDiffTheta(N_i, eCvec, 2 * v_used[i] + 2, x, z, s, t, o);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                             : exp(log(eAU[i]) + log(eBU[i]) - log(eCL)));
          dL[i] = ((eAL[i] < 0.0 || eBL[i] < 0.0)
                       ? static_cast<double>(0.0)
                       : exp(log(eAL[i]) + log(eBL[i]) - log(eCU)));
          if (dLold[i] > dL[i] || dUold[i] < dU[i] || dL[i] > dU[i] ||
              static_cast<double>(eAU[i]) < 0.0 ||
              static_cast<double>(eBU[i]) < 0.0 ||
              static_cast<double>(eAL[i]) > 1.0 ||
              static_cast<double>(eBL[i]) > 1.0 ||
              static_cast<double>(eCU) < 0.0) {
            cerr << "Numerical error detected: (m,k,l,j) = " << mi << "," << ki
                 << ", *, *), dL[" << i << "] = " << dL[i] << ", dU[" << i
                 << "] = " << dU[i] << ", eAU[" << i << "] = " << eAU[i]
                 << ", eAL[" << i << "] = " << eAL[i] << ",  eBU[" << i
                 << "] = " << eBU[i] << ", eBL[" << i << "] = " << eBL[i]
                 << ", eCU = " << eCU << ", eCL = " << eCL
                 << ". Resorting to G1984-style approximation (x,z,s,t) = ("
                 << x << "," << z << "," << s << "," << t << ") ..." << endl;
            return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
          }
          precomputeA(N_i_s, N_i_t, mi, ki);
          double constant_factors_i =
              (!(x > 0.0))
                  ? (factorials[ki] - lg_theta[N_i_s][mi + ki] +
                     lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                     lg_theta1[N_i_s][0] - lg_theta2[N_i_s][mi] +
                     static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                     static_cast<double>(thetaP[N_i_t][1] + ki - 1.0) * log1z)
                  : ((!(x < 1.0))
                         ? (factorials[ki] - lg_theta[N_i_s][mi + ki] +
                            lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                            lg_theta1[N_i_s][mi] - lg_theta2[N_i_s][0] +
                            static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                            static_cast<double>(thetaP[N_i_t][1] + ki - 1.0) *
                                log1z)
                         : ((!(z > 0.0))
                                ? (factorials[mi] - lg_theta[N_i_s][mi + ki] +
                                   lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                                   lg_theta1[N_i_t][0] - lg_theta2[N_i_t][ki] +
                                   static_cast<double>(mi) * log1x)
                                : (factorials[mi] - lg_theta[N_i_s][mi + ki] +
                                   lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                                   lg_theta1[N_i_t][ki] - lg_theta2[N_i_t][0] +
                                   static_cast<double>(mi) * log1x)));
          int index_upper =
              (i == n ? index : ((!(x > 0.0) || !(x < 1.0)) ? ki : mi));
          for (int index2 = 0; index2 <= index_upper; ++index2) {
            double var_stuff_i =
                (!(x > 0.0))
                    ? (-factorials[index2] - factorials[ki - index2] +
                       lg_theta1[N_i_s][index2] +
                       lg_theta2[N_i_s][mi + ki - index2] -
                       lg_theta1[N_i_t][index2] -
                       lg_theta2[N_i_t][ki - index2] +
                       static_cast<double>(index2) * (logz - log1z))
                    : ((!(x < 1.0))
                           ? (-factorials[index2] - factorials[ki - index2] +
                              lg_theta1[N_i_s][mi + index2] +
                              lg_theta2[N_i_s][ki - index2] -
                              lg_theta1[N_i_t][index2] -
                              lg_theta2[N_i_t][ki - index2] +
                              static_cast<double>(index2) * (logz - log1z))
                           : ((!(z > 0.0))
                                  ? (-factorials[index2] -
                                     factorials[mi - index2] +
                                     lg_theta1[N_i_s][index2] +
                                     lg_theta2[N_i_s][mi - index2 + ki] -
                                     lg_theta1[N_i_s][index2] -
                                     lg_theta2[N_i_s][mi - index2] +
                                     static_cast<double>(index2) *
                                         (logx - log1x))
                                  : (-factorials[index2] -
                                     factorials[mi - index2] +
                                     lg_theta1[N_i_s][index2 + ki] +
                                     lg_theta2[N_i_s][mi - index2] -
                                     lg_theta1[N_i_s][index2] -
                                     lg_theta2[N_i_s][mi - index2] +
                                     static_cast<double>(index2) *
                                         (logx - log1x))));
            Amklj = exp(constant_factors_i + var_stuff_i);
            currsumL += Amklj * dL[i];
            currsumU += Amklj * dU[i];
            if (o.debug > 3)
              cerr << "Recomputing currsums with A(n,m,k,index) = A(" << i
                   << "," << mi << "," << ki << "," << index2 << ") = " << Amklj
                   << endl;
          }

          if (currsumL > currsumU) {
            cerr << "Error: currsumL = " << currsumL << " > " << currsumU
                 << " = currsumU (n,m,k,index) = (" << n << "," << m << "," << k
                 << "," << index << ")." << endl;
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
          cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
               << " = currsumL (n,m,k,index) = (" << n << "," << m << "," << k
               << "," << index << ")." << endl;
          exit(1);
        }
        if (currsumUold < currsumU) {
          cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
               << " = currsumU (n,m,k,index) = (" << n << "," << m << "," << k
               << "," << index << ")." << endl;
          exit(1);
        }

        decision_on_mklj_made = (currsumL > u || currsumU < u);
      }

      mklj_found = (currsumL > u);
      if (mklj_found) {
        curr_mk[n][2] = (!(x > 0.0)) ? 0 : ((!(x < 1.0)) ? m : index);
        curr_mk[n][3] = (!(z > 0.0)) ? 0 : ((!(z < 1.0)) ? k : index);
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaInteriorOneQApprox(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen) {
  assert((x > 0.0) && (x < 1.0) && (z > 0.0) && (z < 1.0) && (t > s) &&
         (s > 0.0));
  // ofstream outFile;
  // outFile.open("debugging.txt");
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  bool ind1 =
      (s > o.g1984threshold);  /// Figure out whether to approximate q_m or q_k
  double sorts = (ind1 ? s : t - s), sortsapprox = (sorts == s ? t - s : s);
  size_t N_i_sorts = (sorts == s) ? N_i_s : N_i_t,
         N_i_sortsapprox = (N_i_sorts == N_i_s) ? N_i_t : N_i_s;

  vector<vector<int>> curr_mk;
  vector<int> first_mk(4, 0), mkljStore;
  first_mk[2] = 1;
  first_mk[3] = -1;
  curr_mk.push_back(first_mk);
  bool mklj_found = false;

  /// Set up the necessary quantities
  boost::random::uniform_01<double> U01;
  double u = U01(gen);
  vector<double> eCvec, eAL, eAU, dL, dU, currsumStore;
  double currsumU = 0.0, currsumL = 0.0, eCL = 0.0,
         eCU(GetdDiffThetaOneQApprox(N_i, eCvec, 0, x, z, s, t, o)),
         qApprox = 0.0, runningMax = 1.0e-300;
  vector<int> v_used;
  double Amklj;
  int n = -1, Fmklj = 0, eCindex_computed = 0;
  /// Mechanism to check whether we have run out of precision due to Gaussian
  /// approximations used for one side of the bridge interval - very rare to
  /// be needed
  double mmode = GriffithsParas(N_i_s, s).first,
         vm = GriffithsParas(N_i_s, s).second,
         kmode = GriffithsParas(N_i_t, t - s).first,
         vk = GriffithsParas(N_i_t, t - s).second;
  int mlimU = static_cast<int>(ceil(mmode + 5.0 * sqrt(vm))),
      klimU = static_cast<int>(ceil(kmode + 5.0 * sqrt(vk))),
      mlimL = max(0, static_cast<int>(floor(mmode - 5.0 * sqrt(vm)))),
      klimL = max(0, static_cast<int>(floor(kmode - 5.0 * sqrt(vk))));
  int totalpts = (mlimU - mlimL) * (klimU - klimL), counter = 0;

  /// Compute F ensuring convergence of upper and lower bounds
  pair<vector<int>, double> Csorts, Ct;
  Csorts.second = sorts;
  Ct.second = t;
  int F5;
  if (thetaP[N_i].empty()) {
    F5 = static_cast<int>(
        ceil(((1.0 - x) / (1.0 - z)) + (x / z) * (1.0 + z) * o.eps));
  } else {
    F5 = static_cast<int>(ceil(
        2.0 *
        (max((theta[N_i] / thetaP[N_i][1]) * (1.0 - static_cast<double>(z)),
             (1.0 + theta[N_i]) / (1.0 - static_cast<double>(z))) *
             (1.0 - static_cast<double>(x)) +
         (1.0 + (1.0 / static_cast<double>(z))) *
             (max(
                 ((static_cast<double>(z) * theta[N_i]) + 1.0) / thetaP[N_i][0],
                 1.0)) *
             static_cast<double>(x)) /
        o.eps));
  }
  int F3 = computeE(N_i, Ct),
      F4 = static_cast<int>(ceil(
          max(0.0, 1.0 / static_cast<double>(t) - (theta[N_i] + 1.0) / 2.0)));
  while ((theta[N_i] + 2 * F4 + 1) * exp(-(2 * F4 + theta[N_i]) * t / 2.0) >=
         1 - o.eps)
    ++F4;
  int G = max(max(F3, F4), F5);

  while (!mklj_found) {
    ++n;
    if (n > 0) curr_mk.push_back(curr_mk.back());
    increment_on_mk(true, N_i, curr_mk.back(), s, t);
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
    int F1 = computeC(N_i_sorts, mork, Csorts);

    if (o.debug > 2) cerr << "(F1,G) = (" << F1 << "," << G << ")" << endl;
    Fmklj = max(F1, G);

    int v = -1;
    if (o.debug > 2) {
      cerr << "F(" << setprecision(0) << m << "," << k << "," << setprecision(4)
           << x << "," << z << "," << s << "," << t << ") = " << Fmklj << endl;
    }

    while (2 * v < Fmklj) {
      ++v;
      double newcoefficientU =
          exp(Getlogakm<double>(N_i_sorts, mork + 2 * v, mork) +
              static_cast<double>(-(mork + 2 * v) *
                                  (mork + 2 * v + theta[N_i_sorts] - 1) *
                                  sorts / 2.0));
      double newcoefficientL =
          exp(Getlogakm<double>(N_i_sorts, mork + 2 * v + 1, mork) +
              static_cast<double>(-(mork + 2 * v + 1) *
                                  (mork + 2 * v + 1 + theta[N_i_sorts] - 1) *
                                  sorts / 2.0));

      eAU[n] = eAL[n] + newcoefficientU;
      eAL[n] = eAU[n] - newcoefficientL;

      qApprox = DiscretisedNormCDF(N_i_sortsapprox, morkapprox, sortsapprox);

      if (2 * v + 2 > eCindex_computed) {
        assert(2 * v == eCindex_computed);
        eCL =
            eCU - GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v + 1, x, z, s, t, o);
        eCU =
            eCL + GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v + 2, x, z, s, t, o);
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
          return DrawBridgePMFG1984(N_i, x, z, s, t, o, gen);
        }

        break;
      }
    }

    dU[n] = (eCL < 0.0 ? static_cast<double>(nan(""))
                       : exp(log(eAU[n]) + log(qApprox) - log(eCL)));
    dL[n] = (eAL[n] < 0.0 ? static_cast<double>(0.0)
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
    precomputeA(N_i_s, N_i_t, m, k);
    double constant_factors =
        (!(x > 0.0))
            ? (factorials[k] - lg_theta[N_i_s][m + k] + lg_theta[N_i_s][m] +
               lg_theta[N_i_t][k] - lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m] +
               static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
               static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z)
            : ((!(x < 1.0))
                   ? (factorials[k] - lg_theta[N_i_s][m + k] +
                      lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                      lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0] +
                      static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                      static_cast<double>(thetaP[N_i_t][1] + k - 1.0) * log1z)
                   : ((!(z > 0.0))
                          ? (factorials[m] - lg_theta[N_i_s][m + k] +
                             lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                             lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k] +
                             static_cast<double>(m) * log1x)
                          : (factorials[m] - lg_theta[N_i_s][m + k] +
                             lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                             lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0] +
                             static_cast<double>(m) * log1x)));
    int upper_bound = ((!(x > 0.0)) || (!(x < 1.0))) ? k : m;
    for (int index = 0; index <= upper_bound && !mklj_found; ++index) {
      double var_stuff =
          (!(x > 0.0))
              ? (-factorials[index] - factorials[k - index] +
                 lg_theta1[N_i_s][index] + lg_theta2[N_i_s][m + k - index] -
                 lg_theta1[N_i_t][index] - lg_theta2[N_i_t][k - index] +
                 static_cast<double>(index) * (logz - log1z))
              : ((!(x < 1.0))
                     ? (-factorials[index] - factorials[k - index] +
                        lg_theta1[N_i_s][m + index] +
                        lg_theta2[N_i_s][k - index] - lg_theta1[N_i_t][index] -
                        lg_theta2[N_i_t][k - index] +
                        static_cast<double>(index) * (logz - log1z))
                     : ((!(z > 0.0))
                            ? (-factorials[index] - factorials[m - index] +
                               lg_theta1[N_i_s][index] +
                               lg_theta2[N_i_s][m - index + k] -
                               lg_theta1[N_i_s][index] -
                               lg_theta2[N_i_s][m - index] +
                               static_cast<double>(index) * (logx - log1x))
                            : (-factorials[index] - factorials[m - index] +
                               lg_theta1[N_i_s][index + k] +
                               lg_theta2[N_i_s][m - index] -
                               lg_theta1[N_i_s][index] -
                               lg_theta2[N_i_s][m - index] +
                               static_cast<double>(index) * (logx - log1x))));
      Amklj = exp(constant_factors + var_stuff);
      if (o.debug > 2)
        cerr << "Adding to currsums with A(n,m,k,index) = A(" << n << "," << m
             << "," << k << "," << index << ") = " << endl;

      if (Amklj * dL[n] > 1.0 || Amklj * dU[n] < 0.0 || eAU[n] < 0.0 ||
          eAL[n] > 1.0 || eCU < 0.0) {
        cerr << "Numerical error detected: (m,k,index) = " << m << "," << k
             << "," << index << "), Amklj = " << Amklj << ", dL[" << n
             << "] = " << dL[n] << ", dU[" << n << "] = " << dU[n] << ", eAU["
             << n << "] = " << eAU[n] << ", eAL[" << n << "] = " << eAL[n]
             << ", eCU = " << eCU << ", eCL = " << eCL
             << ". Resorting to G1984-style approximation (x,z,s,t) = (" << x
             << "," << z << "," << s << "," << t << ") ..." << endl;
        return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
      }

      currsumL += Amklj * dL[n];
      currsumU += Amklj * dU[n];

      currsumStore.push_back(log(Amklj) + log(dL[n]));
      runningMax = max(runningMax, currsumStore.back());
      mkljStore.resize(mkljStore.size() + 4, -1);
      mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
      mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
      mkljStore[2 + (4 * (currsumStore.size() - 1))] =
          (!(x > 0.0)) ? 0 : ((!(x < 1.0)) ? m : index);
      mkljStore[3 + (4 * (currsumStore.size() - 1))] =
          (!(z > 0.0)) ? 0 : ((!(z < 1.0)) ? k : index);
      ;

      bool decision_on_mklj_made = (currsumL > u || currsumU < u);

      if (currsumL > currsumU) {
        cerr << "Error: currsumL = " << currsumL << " > " << currsumU
             << " = currsumU (n,m,k,index) = (" << n << "," << m << "," << k
             << "," << index << ")." << endl;
        exit(1);
      }

      while (!decision_on_mklj_made)  /// Refine upper and lower bounds
      {
        double currsumLold = currsumL, currsumUold = currsumU;
        currsumL = 0.0;
        currsumU = 0.0;

        const vector<double> dUold(dU), dLold(dL);
        for (int i = 0; i <= n; ++i) {
          int &mi = curr_mk[i][0], &ki = curr_mk[i][1],
              morki = (ind1 ? mi : ki), morkapproxi = (morki == mi ? ki : mi);
          ;
          ++v_used[i];
          double newcoefficientU =
                     exp(Getlogakm<double>(N_i_sorts, morki + 2 * v_used[i],
                                           morki) +
                         static_cast<double>(
                             -(morki + 2 * v_used[i]) *
                             (morki + 2 * v_used[i] + theta[N_i_sorts] - 1) *
                             sorts / 2.0)),
                 newcoefficientL =
                     exp(Getlogakm<double>(N_i_sorts, morki + 2 * v_used[i] + 1,
                                           morki) +
                         static_cast<double>(-(morki + 2 * v_used[i] + 1) *
                                             (morki + 2 * v_used[i] + 1 +
                                              theta[N_i_sorts] - 1) *
                                             sorts / 2.0));

          eAU[i] = eAL[i] + newcoefficientU;
          eAL[i] = eAU[i] - newcoefficientL;

          qApprox =
              DiscretisedNormCDF(N_i_sortsapprox, morkapproxi, sortsapprox);

          if (2 * v_used[i] + 2 > eCindex_computed) {
            assert(2 * v_used[i] == eCindex_computed);
            eCL = eCU - GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v_used[i] + 1,
                                                x, z, s, t, o);
            eCU = eCL + GetdDiffThetaOneQApprox(N_i, eCvec, 2 * v_used[i] + 2,
                                                x, z, s, t, o);
            eCindex_computed += 2;
          }

          dU[i] = (eCL < 0.0 ? static_cast<double>(nan(""))
                             : exp(log(eAU[i]) + log(qApprox) - log(eCL)));
          dL[i] = (eAL[i] < 0.0 ? static_cast<double>(0.0)
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
            return DrawBridgePMFDiffThetaApprox(N_i, x, z, s, t, o, gen);
          }
          precomputeA(N_i_s, N_i_t, mi, ki);
          double constant_factors_i =
              (!(x > 0.0))
                  ? (factorials[ki] - lg_theta[N_i_s][mi + ki] +
                     lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                     lg_theta1[N_i_s][0] - lg_theta2[N_i_s][mi] +
                     static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                     static_cast<double>(thetaP[N_i_t][1] + ki - 1.0) * log1z)
                  : ((!(x < 1.0))
                         ? (factorials[ki] - lg_theta[N_i_s][mi + ki] +
                            lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                            lg_theta1[N_i_s][mi] - lg_theta2[N_i_s][0] +
                            static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                            static_cast<double>(thetaP[N_i_t][1] + ki - 1.0) *
                                log1z)
                         : ((!(z > 0.0))
                                ? (factorials[mi] - lg_theta[N_i_s][mi + ki] +
                                   lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                                   lg_theta1[N_i_t][0] - lg_theta2[N_i_t][ki] +
                                   static_cast<double>(mi) * log1x)
                                : (factorials[mi] - lg_theta[N_i_s][mi + ki] +
                                   lg_theta[N_i_s][mi] + lg_theta[N_i_t][ki] -
                                   lg_theta1[N_i_t][ki] - lg_theta2[N_i_t][0] +
                                   static_cast<double>(mi) * log1x)));
          int index_upper =
              (i == n ? index : ((!(x > 0.0) || !(x < 1.0)) ? ki : mi));
          for (int index2 = 0; index2 <= index_upper; ++index2) {
            double var_stuff_i =
                (!(x > 0.0))
                    ? (-factorials[index2] - factorials[ki - index2] +
                       lg_theta1[N_i_s][index2] +
                       lg_theta2[N_i_s][mi + ki - index2] -
                       lg_theta1[N_i_t][index2] -
                       lg_theta2[N_i_t][ki - index2] +
                       static_cast<double>(index2) * (logz - log1z))
                    : ((!(x < 1.0))
                           ? (-factorials[index2] - factorials[ki - index2] +
                              lg_theta1[N_i_s][mi + index2] +
                              lg_theta2[N_i_s][ki - index2] -
                              lg_theta1[N_i_t][index2] -
                              lg_theta2[N_i_t][ki - index2] +
                              static_cast<double>(index2) * (logz - log1z))
                           : ((!(z > 0.0))
                                  ? (-factorials[index2] -
                                     factorials[mi - index2] +
                                     lg_theta1[N_i_s][index2] +
                                     lg_theta2[N_i_s][mi - index2 + ki] -
                                     lg_theta1[N_i_s][index2] -
                                     lg_theta2[N_i_s][mi - index2] +
                                     static_cast<double>(index2) *
                                         (logx - log1x))
                                  : (-factorials[index2] -
                                     factorials[mi - index2] +
                                     lg_theta1[N_i_s][index2 + ki] +
                                     lg_theta2[N_i_s][mi - index2] -
                                     lg_theta1[N_i_s][index2] -
                                     lg_theta2[N_i_s][mi - index2] +
                                     static_cast<double>(index2) *
                                         (logx - log1x))));
            Amklj = exp(constant_factors_i + var_stuff_i);
            currsumL += Amklj * dL[i];
            currsumU += Amklj * dU[i];
            currsumStore[i] = log(Amklj) + log(dL[i]);
            runningMax = max(runningMax, currsumStore[i]);
            if (o.debug > 3)
              cerr << "Recomputing currsums with A(n,m,k,index) = A(" << i
                   << "," << mi << "," << ki << "," << index2 << ") = " << Amklj
                   << endl;
          }
        }

        if (currsumL > currsumU) {
          cerr << "Error: currsumL = " << currsumL << " > " << currsumU
               << " = currsumU (n,m,k,index) = (" << n << "," << m << "," << k
               << "," << index << ")." << endl;
          exit(1);
        }

        if (o.debug > 2) {
          cerr << "\ndL ";
          printVec(dL, cerr);
          cerr << "\ndU ";
          printVec(dU, cerr);
        }

        if (currsumLold > currsumL) {
          cerr << "Error: currsumLold = " << currsumLold << " > " << currsumL
               << " = currsumL (n,m,k,index) = (" << n << "," << m << "," << k
               << "," << index << ")." << endl;
          exit(1);
        }
        if (currsumUold < currsumU) {
          cerr << "Error: currsumUold = " << currsumUold << " < " << currsumU
               << " = currsumU (n,m,k,index) = (" << n << "," << m << "," << k
               << "," << index << ")." << endl;
          exit(1);
        }

        decision_on_mklj_made = (currsumL > u || currsumU < u);
      }

      mklj_found = (currsumL > u);
      if (mklj_found) {
        curr_mk[n][2] = (!(x > 0.0)) ? 0 : ((!(x < 1.0)) ? m : index);
        curr_mk[n][3] = (!(z > 0.0)) ? 0 : ((!(z < 1.0)) ? k : index);
        ;
      }
    }
  }

  if (counter == totalpts)  /// Gaussian approximation leads to currsum
                            /// summing to < 1.0, so we renormalise and sample
  {
    LogSumExp(currsumStore, runningMax);
    double sum = 0.0;
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

vector<int> WrightFisher::DrawBridgePMFDiffThetaInteriorApprox(
    size_t N_i, double x, double z, double s, double t, const Options &o,
    boost::random::mt19937 &gen) {
  assert(((!(x > 0.0)) || (!(x < 1.0)) || (!(z > 0.0)) || (!(z < 1.0))) &&
         (t > s) && (s > 0.0));
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  double logx = log(x), log1x = log(1.0 - x), logz = log(z),
         log1z = log(1.0 - z);
  vector<int> returnvec;
  vector<double> currsumStore;
  vector<int> mkljStore;
  bool mklj_found = false;  /// earlyStop gauges whether we can use currsum
  /// as is, or whether we should compute all
  /// probabilities and then sample
  vector<int> modeGuess =
      (!(x > 0.0) || !(x < 1.0))
          ? mkjModeFinder(true, N_i, x, z, s, t, o)
          : mklModeFinder(true, N_i, x, z, s, t,
                          o);  /// Get a guess on mode over (m,k,l/j)

  /// Compute denominator
  double eC = 0.0, emCInc = 1.0, emCOldInc = 1.0, eC_threshold = 1.0e-50;
  int Dmflip = 1, Dkflip = 1, Dmm = 0, Dmp = 0, Dkm = 0, Dkp = 0;
  int dmmode = modeGuess[0],  // static_cast<int>(ceil(GriffithsParas(N_i_s,
                              // s).first)),
      dm = dmmode, dmFlip = 1, dmD = 0, dmU = 0;
  bool dmSwitch = false, dmUpSwitch = false, dmDownSwitch = false;
  int dkmode = modeGuess[1];  // static_cast<int>(ceil(GriffithsParas(N_i_t, t -
                              // s).first));

  while (!dmSwitch) {
    emCOldInc = emCInc;
    double ekC = 0.0, ekCInc = 1.0, ekCOldInc = 1.0;
    int dk = dkmode, dkFlip = 1, dkD = 0, dkU = 0;
    bool dkSwitch = false, dkUpSwitch = false, dkDownSwitch = false;
    while (!dkSwitch) {
      ekCOldInc = ekCInc;
      ekCInc = QmApprox(N_i_t, dk, t - s, o) *
               calculate_expectation(N_i, dm, dk, x, z);
      ekC += ekCInc;
      if (!(dkDownSwitch))  /// Switching mechanism for j
      {
        if (sgn(dk - dkmode) <= 0) {
          dkDownSwitch = ((ekCInc < eC_threshold) || (dkmode - dkD - 1) < 0);
        }
      }

      if (!(dkUpSwitch)) {
        if (sgn(dk - dkmode) >= 0) {
          dkUpSwitch = (ekCInc < eC_threshold);
        }
      }

      dkSwitch = (dkDownSwitch && dkUpSwitch);
      if (!dkSwitch) {
        if (dkFlip == 1) {
          dkU++;
          dk = dkmode + dkU;
          dkFlip *= (dkDownSwitch ? 1 : -1);
        } else if ((dkFlip == -1) && (dkmode - dkD - 1 >= 0)) {
          dkD++;
          dk = dkmode - dkD;
          dkFlip *= (dkUpSwitch ? 1 : -1);
        }
      }
    }
    emCInc = QmApprox(N_i_s, dm, s, o) * ekC;
    eC += emCInc;

    if (!(dmDownSwitch))  /// Switching mechanism for j
    {
      if (sgn(dm - dmmode) <= 0) {
        dmDownSwitch = ((emCInc < eC_threshold) || (dmmode - dmD - 1) < 0);
      }
    }

    if (!(dmUpSwitch)) {
      if (sgn(dm - dmmode) >= 0) {
        dmUpSwitch = (emCInc < eC_threshold);
      }
    }

    dmSwitch = (dmDownSwitch && dmUpSwitch);
    if (!dmSwitch) {
      if (dmFlip == 1) {
        dmU++;
        dm = dmmode + dmU;
        dmFlip *= (dmDownSwitch ? 1 : -1);
      } else if ((dmFlip == -1) && (dmmode - dmD - 1 >= 0)) {
        dmD++;
        dm = dmmode - dmD;
        dmFlip *= (dmUpSwitch ? 1 : -1);
      }
    }
  }

  int mMode = modeGuess[0], kMode = modeGuess[1], indexMode = modeGuess[2];

  boost::random::uniform_01<double>
      U01;  /// Use this guess & eC to compute a suitable threshold for
  /// subsequent computations
  double currsum = 0.0, u = U01(gen),
         thr_term = (!(x > 0.0) || !(x < 1.0))
                        ? mkjModeFinder_Evaluator(true, N_i, mMode, kMode,
                                                  indexMode, x, z, s, t, o)
                        : mklModeFinder_Evaluator(true, N_i, mMode, kMode,
                                                  indexMode, x, z, s, t, o),
         threshold = exp(thr_term - log(eC)) * 1.0e-4;

  int m = mMode, mFlip = 1, mD = 0, mU = 0;
  bool mSwitch = false, mDownSwitch = false,
       mUpSwitch = false;  /// Compute m contributions
  double mContr, runningMax = -1.0e100;
  double constant =
      (!(x > 0.0))
          ? -lg_theta1[N_i_s][0] +
                static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                static_cast<double>(thetaP[N_i_t][1] - 1.0) * log1z
          : ((!(x < 1.0))
                 ? -lg_theta2[N_i_s][0] +
                       static_cast<double>(thetaP[N_i_t][0] - 1.0) * logz +
                       static_cast<double>(thetaP[N_i_t][1] - 1.0) * log1z
                 : ((!(z > 0.0)) ? -lg_theta1[N_i_t][0]
                                 : -lg_theta2[N_i_t][0]));

  while (!mSwitch) {
    precomputeA(N_i_s, N_i_t, m, kMode);
    double qm = QmApprox(N_i_s, m, s, o);
    if (!(qm > 0.0))  /// This should not trigger, but if it does, sets to
                      /// very small value (taking logs later so cannot be 0!)
    {
      qm = 1.0e-300;
    }
    mContr =
        lg_theta[N_i_s][m] +
        ((!(x > 0.0))
             ? -lg_theta2[N_i_s][m]
             : ((!(x < 1.0)) ? -lg_theta1[N_i_s][m]
                             : factorials[m] + static_cast<double>(m) * log1x));
    mContr += log(qm);

    int k = kMode, kFlip = 1, kD = 0, kU = 0;
    bool kSwitch = false, kDownSwitch = false,
         kUpSwitch = false;  /// Compute k contributions
    double kContr;

    while (!kSwitch) {
      precomputeA(N_i_s, N_i_t, m, k);
      double qk = QmApprox(N_i_t, k, t - s, o);
      if (!(qk > 0.0))  /// This should not trigger, but if it does, sets to
                        /// very small value (taking logs later so cannot be 0!)
      {
        qk = 1.0e-300;
      }
      if (!(qk > 0.0)) {
        break;
      }
      kContr = lg_theta[N_i_t][k] - lg_theta[N_i_s][m + k] +
               ((!(z > 0.0))
                    ? -lg_theta2[N_i_t][k]
                    : ((!(z < 1.0))
                           ? -lg_theta1[N_i_t][k]
                           : factorials[k] + static_cast<double>(k) * log1z));
      kContr += log(qk);

      int index_upper = (!(x > 0.0) || !(x < 1.0)) ? k : m;
      int other_upper = (index_upper == m) ? k : m;
      int indexFlip = 1, indexU = 0, indexD = 0,
          newindexMode = min(indexMode, index_upper),
          index = newindexMode;  /// Redefine lMode in case m is too small!
      bool indexSwitch = false, indexDownSwitch = false, indexUpSwitch = false;
      /// Compute l contributions
      double indexContr;

      while (!indexSwitch) {
        indexContr =
            -factorials[index] - factorials[index_upper - index] +
            ((!(x > 0.0))
                 ? lg_theta1[N_i_s][index] +
                       lg_theta2[N_i_s][other_upper + index_upper - index] -
                       lg_theta1[N_i_t][index] -
                       lg_theta2[N_i_t][index_upper - index] +
                       static_cast<double>(index) * (logz - log1z)
                 : ((!(x < 1.0))
                        ? lg_theta1[N_i_s][other_upper + index] +
                              lg_theta2[N_i_s][index_upper - index] -
                              lg_theta1[N_i_t][index] -
                              lg_theta2[N_i_t][index_upper - index] +
                              static_cast<double>(index) * (logz - log1z)
                        : ((!(z > 0.0))
                               ? lg_theta1[N_i_s][index] +
                                     lg_theta2[N_i_s][other_upper +
                                                      index_upper - index] -
                                     lg_theta1[N_i_s][index] -
                                     lg_theta2[N_i_s][index_upper - index] +
                                     static_cast<double>(index) * (logx - log1x)
                               : lg_theta1[N_i_s][other_upper + index] +
                                     lg_theta2[N_i_s][index_upper - index] -
                                     lg_theta1[N_i_s][index] -
                                     lg_theta2[N_i_s][index_upper - index] +
                                     static_cast<double>(index) *
                                         (logx - log1x))));
        double currsum_inc =
            exp(mContr + kContr + indexContr + constant - log(eC));
        runningMax =
            max(currsum_inc,
                runningMax);  /// Running max needed for log-sum-exp trick later
        currsum += currsum_inc;
        currsumStore.push_back(currsum);

        mkljStore.resize(mkljStore.size() + 4, -1);
        mkljStore[0 + (4 * (currsumStore.size() - 1))] = m;
        mkljStore[1 + (4 * (currsumStore.size() - 1))] = k;
        mkljStore[2 + (4 * (currsumStore.size() - 1))] =
            (!(x > 0.0)) ? 0 : ((!(x < 1.0)) ? m : index);
        mkljStore[3 + (4 * (currsumStore.size() - 1))] =
            (!(z > 0.0)) ? 0 : ((!(z < 1.0)) ? k : index);

        if (currsum > u)  /// if earlyStop is allowed, we can
                          /// stop once currsum exceeds u
        {
          returnvec.push_back(m);
          returnvec.push_back(k);
          returnvec.push_back((!(x > 0.0)) ? 0 : ((!(x < 1.0)) ? m : index));
          returnvec.push_back((!(z > 0.0)) ? 0 : ((!(z < 1.0)) ? k : index));

          mklj_found = true;
          goto End;
        }

        if (!(indexDownSwitch))  /// Switching mechanism for j
        {
          if (sgn(index - newindexMode) <= 0) {
            indexDownSwitch =
                ((currsum_inc < threshold) || (newindexMode - indexD - 1) < 0);
          }
        }

        if (!(indexUpSwitch)) {
          if (sgn(index - newindexMode) >= 0) {
            indexUpSwitch = ((currsum_inc < threshold) ||
                             (newindexMode + indexU + 1) > index_upper);
          }
        }

        indexSwitch = (indexDownSwitch && indexUpSwitch);

        if (!indexSwitch) {
          if ((indexFlip == 1 && (newindexMode + indexU + 1 <= index_upper)) ||
              (indexDownSwitch && !(indexUpSwitch))) {
            indexU++;
            index = newindexMode + indexU;
            indexFlip *= (indexDownSwitch ? 1 : -1);
          } else if ((indexFlip == -1 && (newindexMode - indexD - 1 >= 0)) ||
                     (indexUpSwitch && !(indexDownSwitch))) {
            indexD++;
            index = newindexMode - indexD;
            indexFlip *= (indexUpSwitch ? 1 : -1);
          }
        }
      }

      if (!(kDownSwitch))  /// Switching mechanism for k
      {
        if (sgn(k - kMode) <= 0) {
          kDownSwitch =
              (((indexU == 0) && (indexD == 0)) || (kMode - kD - 1 < 0));
        }
      }

      if (!(kUpSwitch)) {
        if (sgn(k - kMode) >= 0) {
          kUpSwitch = ((indexU == 0) && (indexD == 0));
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
    double sum = 0.0;
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

double WrightFisher::mkModeFinder_Evaluator(
    bool isDiffTheta, size_t N_i, int m, int k, double x, double z, double s,
    double t,
    const Options &o)  /// Evaluation function for finding mode over (m,k)
{
  assert((m >= 0) && (k >= 0) && (x >= 0.0) && (x <= 1.0) && (z >= 0.0) &&
         (z <= 1.0) && (s > 0.0) && (s < t));
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  precomputeA(N_i_s, N_i_t, m, k);
  double qm = QmApprox(N_i_s, m, s, o), qk = QmApprox(N_i_t, k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If
                    /// they are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }

  double ret_val = log(qm) + log(qk);
  if (x != z) {
    if (isDiffTheta) {
      if (!(x > 0.0)) {  // x = 0 && z = 1
        ret_val +=
            lg_theta1[N_i_s][k] + lg_theta2[N_i_s][m] - lg_theta[N_i_s][m + k] -
            lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m] + lg_theta[N_i_s][m] -
            lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0] + lg_theta[N_i_t][k];
      } else {  // x = 1 && z =0
        ret_val +=
            lg_theta1[N_i_s][m] + lg_theta2[N_i_s][k] - lg_theta[N_i_s][m + k] -
            lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0] + lg_theta[N_i_s][m] -
            lg_theta1[N_i_t][0] - lg_theta2[N_i_t][k] + lg_theta[N_i_t][k];
      }
    } else {
      ret_val += lg_theta[N_i][m] + lg_theta[N_i_s][k] - lg_theta[N_i_s][m + k];
    }
  } else if (!(x > 0.0)) {  // x = 0 && z = 0
    if (isDiffTheta) {
      ret_val += lg_theta1[N_i_s][0] + lg_theta2[N_i_s][m + k] -
                 lg_theta[N_i_s][m + k] - lg_theta1[N_i_s][0] -
                 lg_theta2[N_i_s][m] + lg_theta[N_i_s][m] -
                 lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0] + lg_theta[N_i_t][k];
    } else {
      ret_val += lg_theta2[N_i_s][m + k] + lg_theta[N_i_s][m] +
                 lg_theta[N_i][k] - lg_theta[N_i_s][m] - lg_theta2[N_i_s][k] -
                 lg_theta[N_i_s][m + k];
    }
  } else {  // x = 1 && z = 1
    if (isDiffTheta) {
      ret_val += lg_theta1[N_i_s][m + k] + lg_theta2[N_i_s][0] -
                 lg_theta[N_i_s][m + k] - lg_theta1[N_i_s][m] -
                 lg_theta2[N_i_s][0] + lg_theta[N_i_s][m] -
                 lg_theta1[N_i_t][k] - lg_theta2[N_i_t][0] + lg_theta[N_i_t][k];
    } else {
      ret_val += lg_theta1[N_i_s][m + k] + lg_theta[N_i_s][m] +
                 lg_theta[N_i_s][k] - lg_theta1[N_i_s][m] -
                 lg_theta1[N_i_s][k] - lg_theta[N_i_s][m + k];
    }
  }
  return ret_val;
}

double WrightFisher::mkjModeFinder_Evaluator(
    bool isDiffTheta, size_t N_i, int m, int k, int j, double x, double z,
    double s, double t,
    const Options &o)  /// Evaluation function for finding mode over (m,k,j)
{
  assert((m >= 0) && (k >= 0) && (j >= 0) && (j <= k) && (x >= 0.0) &&
         (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > 0.0) && (s < t));
  double logz = log(z), log1z = log(1.0 - z);
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  precomputeA(N_i_s, N_i_t, m, k);
  boost::math::binomial_distribution<> BIN(k, z);
  double qm = QmApprox(N_i_s, m, s, o), qk = QmApprox(N_i_t, k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If
                    /// they are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }
  double ret_val = log(qm) + log(qk) + lg_theta[N_i_s][m] + lg_theta[N_i_s][k] -
                   lg_theta[N_i_s][m + k] + factorials[k] - factorials[j] -
                   factorials[k - j];
  if (!(x > 0.0)) {
    ret_val += lg_theta1[N_i_s][j] + lg_theta2[N_i_s][m + k - j] -
               lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m] - lg_theta1[N_i_t][j] -
               lg_theta2[N_i_t][k - j] +
               static_cast<double>(thetaP[N_i_t][0] + j - 1.0) * logz +
               static_cast<double>(thetaP[N_i_t][1] + k - j - 1.0) * log1z;
  } else {
    ret_val += lg_theta1[N_i_s][m + j] + lg_theta2[N_i_s][k - j] -
               lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0] - lg_theta1[N_i_t][j] -
               lg_theta2[N_i_t][k - j] +
               static_cast<double>(thetaP[N_i_t][0] + j - 1.0) * logz +
               static_cast<double>(thetaP[N_i_t][1] + k - j - 1.0) * log1z;
  }
  return ret_val;
}

double WrightFisher::mkjDensityModeFinder_Evaluator(
    bool isDiffTheta, size_t N_i, int m, int k, int j, double x, double z,
    double y, double s, double t,
    const Options &o)  /// Evaluation function for finding mode over (m,k,j)
{
  assert((m >= 0) && (k >= 0) && (j >= 0) && (j <= k) && (x >= 0.0) &&
         (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > 0.0) && (s < t));
  double logz = log(z), log1z = log(1.0 - z), logy = log(y),
         log1y = log(1.0 - y);
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  precomputeA(N_i_s, N_i_t, m, k);
  boost::math::binomial_distribution<> BIN(k, z);
  double qm = QmApprox(N_i_s, m, s, o), qk = QmApprox(N_i_t, k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If
                    /// they are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }
  double ret_val = log(qm) + log(qk) + factorials[k] - factorials[j] -
                   factorials[k - j] + lg_theta[N_i_s][m] + lg_theta[N_i_t][k];
  if (!(x > 0.0)) {
    ret_val += -lg_theta1[N_i_s][0] - lg_theta2[N_i_s][m] -
               lg_theta1[N_i_t][j] - lg_theta2[N_i_t][k - j] +
               static_cast<double>(thetaP[N_i_t][0] + j - 1.0) * logz +
               static_cast<double>(thetaP[N_i_t][1] + k - j - 1.0) * log1z +
               static_cast<double>(thetaP[N_i_s][0] + j - 1.0) * logy +
               static_cast<double>(thetaP[N_i_s][1] + m + k - j - 1.0) * log1y;
  } else {
    ret_val += -lg_theta1[N_i_s][m] - lg_theta2[N_i_s][0] -
               lg_theta1[N_i_t][j] - lg_theta2[N_i_t][k - j] +
               static_cast<double>(thetaP[N_i_t][0] + j - 1.0) * logz +
               static_cast<double>(thetaP[N_i_t][1] + k - j - 1.0) * log1z +
               static_cast<double>(thetaP[N_i_s][0] + m + j - 1.0) * logy +
               static_cast<double>(thetaP[N_i_s][1] + k - j - 1.0) * log1y;
  }
  return ret_val;
}

double WrightFisher::mkljModeFinder_Evaluator(
    bool isDiffTheta, size_t N_i, int m, int k, int l, int j, double x,
    double z, double s, double t,
    const Options &o)  /// Evaluation function for finding mode over (m,k,l,j)
{
  assert((m >= 0) && (k >= 0) && (j >= 0) && (j <= k) && (l >= 0) && (l <= m) &&
         (x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > 0.0) &&
         (s < t));
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  boost::math::binomial_distribution<> BIN(m, x), BINZ(k, z);
  double qm = QmApprox(N_i_s, m, s, o), qk = QmApprox(N_i_t, k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If
                    /// they are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }
  return log(qm) + log(qk) + log(pdf(BIN, l)) + log(pdf(BINZ, j)) +
         boost::math::lgamma(static_cast<double>(thetaP[N_i_s][0] + l + j)) +
         boost::math::lgamma(
             static_cast<double>(thetaP[N_i_s][1] + m - l + k - j)) +
         boost::math::lgamma(static_cast<double>(theta[N_i_s] + m)) -
         boost::math::lgamma(static_cast<double>(theta[N_i_s] + m + k)) -
         boost::math::lgamma(static_cast<double>(thetaP[N_i_s][0] + l)) -
         boost::math::lgamma(static_cast<double>(thetaP[N_i_s][1] + m - l)) +
         boost::math::lgamma(static_cast<double>(theta[N_i_t] + k)) -
         boost::math::lgamma(static_cast<double>(thetaP[N_i_t][0] + j)) -
         boost::math::lgamma(static_cast<double>(thetaP[N_i_t][1] + k - j));
}

double WrightFisher::mklModeFinder_Evaluator(
    bool isDiffTheta, size_t N_i, int m, int k, int l, double x, double z,
    double s, double t,
    const Options &o)  /// Evaluation function for finding mode over (m,k,l)
{
  assert((m >= 0) && (k >= 0) && (l >= 0) && (l <= m) && (x > 0.0) &&
         (x < 1.0) && (!(z > 0.0) || !(z < 1.0)) && (s > 0.0) && (s < t));
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  precomputeA(N_i_s, N_i_t, m, k);
  double logx = log(x), log1x = log(1.0 - x);
  double qm = QmApprox(N_i, m, s, o), qk = QmApprox(N_i, k, t - s, o);
  if (!(qm > 0.0))  /// Ensure qm and qk are not zero when taking logs! If
                    /// they are, set to a very small positive value
  {
    qm = 1.0e-300;
  }
  if (!(qk > 0.0)) {
    qk = 1.0e-300;
  }
  double ret_val = log(qm) + log(qk) + lg_theta[N_i_s][m] + lg_theta[N_i_t][k] -
                   lg_theta[N_i_s][m + k] + factorials[m] - factorials[l] -
                   factorials[m - l] + static_cast<double>(l) * logx +
                   static_cast<double>(m - l) * log1x;
  if (!(z > 0.0)) {
    ret_val += lg_theta1[N_i_s][l] + lg_theta2[N_i_s][m - l + k] -
               lg_theta1[N_i_s][l] - lg_theta2[N_i_s][m - l] -
               lg_theta1[N_i_t][0] - lg_theta[N_i_t][k];
  } else {
    ret_val += lg_theta1[N_i_s][l + k] + lg_theta2[N_i_s][m - l] -
               lg_theta1[N_i_s][l] - lg_theta2[N_i_s][m - l] -
               lg_theta1[N_i_t][k] - lg_theta[N_i_t][0];
  }
  return ret_val;
}

vector<int> WrightFisher::mkModeFinder(
    bool isDiffTheta, size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Routine for finding mode over (m,k)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(N_i, s).first)),
      k = static_cast<int>(floor(GriffithsParas(N_i, t - s).first));

  int m_ud, k_ud;  /// Start at mode from qm and qk
  double currMode_eval =
      mkModeFinder_Evaluator(isDiffTheta, N_i, m, k, x, z, s, t, o);

  bool stop = false;

  /// Iteratively increment m and k depending on whether function increases or
  /// decreases with proposed move
  while (!stop) {
    if (currMode_eval <
        mkModeFinder_Evaluator(isDiffTheta, N_i, m + 1, k, x, z, s, t, o)) {
      m_ud = 1;
    } else if ((m - 1 >= 0) &&
               currMode_eval < mkModeFinder_Evaluator(isDiffTheta, N_i, m - 1,
                                                      k, x, z, s, t, o)) {
      m_ud = -1;
    } else {
      m_ud = 0;
    }

    m += m_ud;
    currMode_eval =
        mkModeFinder_Evaluator(isDiffTheta, N_i, m, k, x, z, s, t, o);
    if (currMode_eval <
        mkModeFinder_Evaluator(isDiffTheta, N_i, m, k + 1, x, z, s, t, o)) {
      k_ud = 1;
    } else if ((k - 1 >= 0) &&
               currMode_eval < mkModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                      k - 1, x, z, s, t, o)) {
      k_ud = -1;
    } else {
      k_ud = 0;
    }

    k += k_ud;
    currMode_eval =
        mkModeFinder_Evaluator(isDiffTheta, N_i, m, k, x, z, s, t, o);

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
    bool isDiffTheta, size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Routine for finding mode over (m,k,j)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(N_i, s).first)),
      k = static_cast<int>(floor(GriffithsParas(N_i, t - s).first));
  int j = static_cast<int>(floor(static_cast<double>(k) * z));

  int m_ud, k_ud, j_ud;  /// Starting from the modes of qm, qk and Bin(k,z)
  double currMode_eval =
      mkjModeFinder_Evaluator(isDiffTheta, N_i, m, k, j, x, z, s, t, o);

  bool stop = false;

  /// Iteratively increment m and (k,j) depending on whether function
  /// increases or decreases with proposed move - (k,j) updated jointly to
  /// make sure 0 <= j <= k at all times
  while (!stop) {
    if (currMode_eval <
        mkjModeFinder_Evaluator(isDiffTheta, N_i, m + 1, k, j, x, z, s, t, o)) {
      m_ud = 1;
    } else if ((m - 1 >= 0) &&
               currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m - 1,
                                                       k, j, x, z, s, t, o)) {
      m_ud = -1;
    } else {
      m_ud = 0;
    }

    m += m_ud;
    currMode_eval =
        mkjModeFinder_Evaluator(isDiffTheta, N_i, m, k, j, x, z, s, t, o);

    if (currMode_eval <
        mkjModeFinder_Evaluator(isDiffTheta, N_i, m, k + 1, j, x, z, s, t, o)) {
      k_ud = 1;
      if (currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m, k + 1,
                                                  j + 1, x, z, s, t, o)) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                          k + 1, j - 1, x, z, s,
                                                          t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j <= k - 1) &&
               currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                       k - 1, j, x, z, s, t,
                                                       o)) {
      k_ud = -1;
      if ((j + 1 <= k - 1) &&
          (currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m, k - 1,
                                                   j + 1, x, z, s, t, o))) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                          k - 1, j - 1, x, z, s,
                                                          t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j - 1 >= 0) &&
               (currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                        k - 1, j - 1, x, z, s,
                                                        t, o))) {
      k_ud = -1;
      j_ud = -1;
    } else {
      k_ud = 0;
      if (k == 0) {
        j_ud = 0;
      } else {
        if ((j + 1 <= k) &&
            (currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m, k,
                                                     j + 1, x, z, s, t, o))) {
          j_ud = 1;
        } else if ((j - 1 >= 0) &&
                   (currMode_eval < mkjModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                            k, j - 1, x, z, s,
                                                            t, o))) {
          j_ud = -1;
        } else {
          j_ud = 0;
        }
      }
    }

    k += k_ud;
    j += j_ud;
    currMode_eval =
        mkjModeFinder_Evaluator(isDiffTheta, N_i, m, k, j, x, z, s, t, o);

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

vector<int> WrightFisher::mkjDensityModeFinder(
    bool isDiffTheta, size_t N_i, double x, double z, double y, double s,
    double t,
    const Options &o)  /// Routine for finding mode over (m,k,j)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(N_i, s).first)),
      k = static_cast<int>(floor(GriffithsParas(N_i, t - s).first));
  int j = static_cast<int>(floor(static_cast<double>(k) * z));

  int m_ud, k_ud, j_ud;  /// Starting from the modes of qm, qk and Bin(k,z)
  double currMode_eval = mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m, k,
                                                        j, x, z, y, s, t, o);

  bool stop = false;

  /// Iteratively increment m and (k,j) depending on whether function
  /// increases or decreases with proposed move - (k,j) updated jointly to
  /// make sure 0 <= j <= k at all times
  while (!stop) {
    if (currMode_eval < mkjDensityModeFinder_Evaluator(
                            isDiffTheta, N_i, m + 1, k, j, x, z, y, s, t, o)) {
      m_ud = 1;
    } else if ((m - 1 >= 0) && currMode_eval < mkjDensityModeFinder_Evaluator(
                                                   isDiffTheta, N_i, m - 1, k,
                                                   j, x, z, y, s, t, o)) {
      m_ud = -1;
    } else {
      m_ud = 0;
    }

    m += m_ud;
    currMode_eval = mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m, k, j, x,
                                                   z, y, s, t, o);

    if (currMode_eval < mkjDensityModeFinder_Evaluator(
                            isDiffTheta, N_i, m, k + 1, j, x, z, y, s, t, o)) {
      k_ud = 1;
      if (currMode_eval < mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                         k + 1, j + 1, x, z, y,
                                                         s, t, o)) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval <
                  mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m, k + 1,
                                                 j - 1, x, z, y, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j <= k - 1) &&
               currMode_eval < mkjDensityModeFinder_Evaluator(isDiffTheta, N_i,
                                                              m, k - 1, j, x, z,
                                                              y, s, t, o)) {
      k_ud = -1;
      if ((j + 1 <= k - 1) && (currMode_eval < mkjDensityModeFinder_Evaluator(
                                                   isDiffTheta, N_i, m, k - 1,
                                                   j + 1, x, z, y, s, t, o))) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval <
                  mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m, k - 1,
                                                 j - 1, x, z, y, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j - 1 >= 0) &&
               (currMode_eval <
                mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m, k - 1,
                                               j - 1, x, z, y, s, t, o))) {
      k_ud = -1;
      j_ud = -1;
    } else {
      k_ud = 0;
      if (k == 0) {
        j_ud = 0;
      } else {
        if ((j + 1 <= k) && (currMode_eval < mkjDensityModeFinder_Evaluator(
                                                 isDiffTheta, N_i, m, k, j + 1,
                                                 x, z, y, s, t, o))) {
          j_ud = 1;
        } else if ((j - 1 >= 0) &&
                   (currMode_eval <
                    mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m, k,
                                                   j - 1, x, z, y, s, t, o))) {
          j_ud = -1;
        } else {
          j_ud = 0;
        }
      }
    }

    k += k_ud;
    j += j_ud;
    currMode_eval = mkjDensityModeFinder_Evaluator(isDiffTheta, N_i, m, k, j, x,
                                                   z, y, s, t, o);

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
    bool isDiffTheta, size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Routine for finding mode over (m,k,l,j)
{
  vector<int> returnvec;
  size_t N_i_s = N_i, N_i_t = (isDiffTheta) ? N_i + 1 : N_i;
  int m = static_cast<int>(floor(GriffithsParas(N_i_s, s).first)),
      k = static_cast<int>(floor(GriffithsParas(N_i_t, t - s).first));
  int l = static_cast<int>(floor(static_cast<double>(m) * x)),
      j = static_cast<int>(floor(static_cast<double>(k) * z));

  int m_ud, k_ud, l_ud,
      j_ud;  /// Initialise at modes from qm, qk, Bin(m,x), Bin(k,z)
  double currMode_eval =
      mkljModeFinder_Evaluator(isDiffTheta, N_i, m, k, l, j, x, z, s, t, o);

  bool stop = false;

  /// Iteratively increment (m,l) and (k,j) depending on whether function
  /// increases or decreases with proposed move - (m,l) & (k,j) updated
  /// jointly to make sure 0 <= l <= m & 0 <= j <= k at all times
  while (!stop) {
    if (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m + 1, k, l,
                                                 j, x, z, s, t, o)) {
      m_ud = 1;
      if (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m + 1, k,
                                                   l + 1, j, x, z, s, t, o)) {
        l_ud = 1;
      } else if ((l - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i,
                                                           m + 1, k, l - 1, j,
                                                           x, z, s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= 0) && (l <= m - 1) &&
               currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m - 1,
                                                        k, l, j, x, z, s, t,
                                                        o)) {
      m_ud = -1;
      if ((l + 1 <= m - 1) &&
          (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m - 1, k,
                                                    l + 1, j, x, z, s, t, o))) {
        l_ud = 1;
      } else if ((l - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i,
                                                           m - 1, k, l - 1, j,
                                                           x, z, s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= 0) && (l - 1 >= 0) &&
               (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i,
                                                         m - 1, k, l - 1, j, x,
                                                         z, s, t, o))) {
      m_ud = -1;
      l_ud = -1;
    } else {
      m_ud = 0;
      if (m == 0) {
        l_ud = 0;
      } else {
        if ((l + 1 <= m) && (currMode_eval < mkljModeFinder_Evaluator(
                                                 isDiffTheta, N_i, m, k, l + 1,
                                                 j, x, z, s, t, o))) {
          l_ud = 1;
        } else if ((l - 1 >= 0) &&
                   (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i,
                                                             m, k, l - 1, j, x,
                                                             z, s, t, o))) {
          l_ud = -1;
        } else {
          l_ud = 0;
        }
      }
    }

    m += m_ud;
    l += l_ud;
    currMode_eval =
        mkljModeFinder_Evaluator(isDiffTheta, N_i, m, k, l, j, x, z, s, t, o);

    if (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m, k + 1, l,
                                                 j, x, z, s, t, o)) {
      k_ud = 1;
      if (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m, k + 1,
                                                   l, j + 1, x, z, s, t, o)) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                           k + 1, l, j - 1, x,
                                                           z, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j <= k - 1) &&
               currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                        k - 1, l, j, x, z, s, t,
                                                        o)) {
      k_ud = -1;
      if ((j + 1 <= k - 1) &&
          currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m, k - 1,
                                                   l, j + 1, x, z, s, t, o)) {
        j_ud = 1;
      } else if ((j - 1 >= 0) &&
                 (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                           k - 1, l, j - 1, x,
                                                           z, s, t, o))) {
        j_ud = -1;
      } else {
        j_ud = 0;
      }
    } else if ((k - 1 >= 0) && (j - 1 >= 0) &&
               (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                         k - 1, l, j - 1, x, z,
                                                         s, t, o))) {
      k_ud = -1;
      j_ud = -1;
    } else {
      k_ud = 0;
      if (k == 0) {
        j_ud = 0;
      } else {
        if ((j + 1 <= k) &&
            (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i, m, k, l,
                                                      j + 1, x, z, s, t, o))) {
          j_ud = 1;
        } else if ((j - 1 >= 0) &&
                   (currMode_eval < mkljModeFinder_Evaluator(isDiffTheta, N_i,
                                                             m, k, l, j - 1, x,
                                                             z, s, t, o))) {
          j_ud = -1;
        } else {
          j_ud = 0;
        }
      }
    }

    k += k_ud;
    j += j_ud;
    currMode_eval =
        mkljModeFinder_Evaluator(isDiffTheta, N_i, m, k, l, j, x, z, s, t, o);

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
    bool isDiffTheta, size_t N_i, double x, double z, double s, double t,
    const Options &o)  /// Routine for finding mode over (m,k,l)
{
  vector<int> returnvec;
  int m = static_cast<int>(floor(GriffithsParas(N_i, s).first)),
      k = static_cast<int>(floor(GriffithsParas(N_i, t - s).first)),
      bumper = (thetaP[N_i].empty() ? 1 : 0);
  int l = max(static_cast<int>(floor(static_cast<double>(m) * x)), bumper);

  int m_ud = -1, k_ud = -1, l_ud = -1,
      mkLower = thetaP[N_i].empty()
                    ? 1
                    : 0;  /// Starting from the modes of qm, qk and Bin(m,x)
  double currMode_eval =
      mklModeFinder_Evaluator(isDiffTheta, N_i, m, k, l, x, z, s, t, o);

  bool rerun = false;

  /// Iteratively increment k and (m,l) depending on whether function
  /// increases or decreases with proposed move - (m,l) updated jointly to
  /// make sure 0 <= l <= m at all times
  while (!((m_ud == 0) && (k_ud == 0) && (l_ud == 0))) {
  Rerun:

    if (currMode_eval <
        mklModeFinder_Evaluator(isDiffTheta, N_i, m, k + 1, l, x, z, s, t, o)) {
      k_ud = 1;
    } else if ((k - 1 >= mkLower) &&
               currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                       k - 1, l, x, z, s, t,
                                                       o)) {
      k_ud = -1;
    } else {
      k_ud = 0;
    }

    k += k_ud;
    currMode_eval =
        mklModeFinder_Evaluator(isDiffTheta, N_i, m, k, l, x, z, s, t, o);
    int lLower = (thetaP[N_i].empty() || (!(thetaP[N_i][0] > 0.0) && !(z > 0.0))
                      ? 1
                      : 0),
        lUpper = (thetaP[N_i].empty() || (!(thetaP[N_i][1] > 0.0) && !(z < 1.0))
                      ? m - 1
                      : m);

    if (currMode_eval <
        mklModeFinder_Evaluator(isDiffTheta, N_i, m + 1, k, l, x, z, s, t, o)) {
      m_ud = 1;
      if (currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i, m + 1, k,
                                                  l + 1, x, z, s, t, o)) {
        l_ud = 1;
      } else if ((l - 1 >= lLower) &&
                 (currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i,
                                                          m + 1, k, l - 1, x, z,
                                                          s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= mkLower) && (l <= lUpper - 1) &&
               currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i, m - 1,
                                                       k, l, x, z, s, t, o)) {
      m_ud = -1;
      if ((l + 1 <= lUpper - 1) &&
          (currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i, m - 1, k,
                                                   l + 1, x, z, s, t, o))) {
        l_ud = 1;
      } else if ((l - 1 >= lLower) &&
                 (currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i,
                                                          m - 1, k, l - 1, x, z,
                                                          s, t, o))) {
        l_ud = -1;
      } else {
        l_ud = 0;
      }
    } else if ((m - 1 >= mkLower) && (l - 1 >= lLower) &&
               (currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i, m - 1,
                                                        k, l - 1, x, z, s, t,
                                                        o))) {
      m_ud = -1;
      l_ud = -1;
    } else {
      m_ud = 0;
      if (m == mkLower) {
        l_ud = 0;
      } else {
        if ((l + 1 <= lUpper) &&
            (currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i, m, k,
                                                     l + 1, x, z, s, t, o))) {
          l_ud = 1;
        } else if ((l - 1 >= lLower) &&
                   (currMode_eval < mklModeFinder_Evaluator(isDiffTheta, N_i, m,
                                                            k, l - 1, x, z, s,
                                                            t, o))) {
          l_ud = -1;
        } else {
          l_ud = 0;
        }
      }
    }

    m += m_ud;
    l += l_ud;
    currMode_eval =
        mklModeFinder_Evaluator(isDiffTheta, N_i, m, k, l, x, z, s, t, o);

    if (rerun == true && !(k_ud == 0 && m_ud == 0 && l_ud == 0)) {
      rerun = false;
    }

    if (k_ud == 0 && m_ud == 0 && l_ud == 0 && rerun == false) {
      rerun = true;
      goto Rerun;
    }

    /// Iterate until told to not update (m,k,l) anymore
  }

  if ((thetaP[N_i].empty() || (!(thetaP[N_i][0] > 0.0) && (l == 0))) ||
      (thetaP[N_i].empty() || (!(thetaP[N_i][1] > 0.0) && (l == m)))) {
    l = ((thetaP[N_i].empty() || (!(thetaP[N_i][0] > 0.0) && (l == 0)))
             ? 1
             : m - 1);  /// Sometimes l == 0/m but cannot be for cases when
                        /// thetaP[N_i].empty()
  }

  returnvec.push_back(m);
  returnvec.push_back(k);
  returnvec.push_back(l);

  return returnvec;
}

pair<double, int> WrightFisher::DrawBridgepoint(
    size_t N_i, double x, double z, double t1, double t2, double s,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Routine to decide which bridge sampler to invoke for
               /// sampling at time s from neutral bridge diffusion started at
               /// x at time t1, ending at z in time t2, conditioned on
               /// non-absorption on (t1,t2)
{
  assert((x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > t1) &&
         (s < t2));
  double y, para1, para2;
  int m, k, j, l, coeffcount = -1;
  vector<int> mklj;

  if ((s - t1 <= o.bridgethreshold || t2 - s <= o.bridgethreshold) ||
      (((theta[N_i] - 1.0) / (exp(0.5 * (theta[N_i] - 1.0) * (s - t1)) - 1.0)) +
           ((theta[N_i] - 1.0) /
            (exp(0.5 * (theta[N_i] - 1.0) * (t2 - s)) - 1.0)) >
       260.0))  /// Use diffusion approximations
  {
    /// Last condition checks that corresponding m and k terms are not too
    /// large
    /// to create computational bottleneck
    double y1 =
        DrawEndpoint(N_i, x, t1, s, o, gen)
            .first;  /// Diffusion approximation for when times are too small
    /// and bridge takes too long to compute

    y = abs((y1 - x) + ((t2 - s) / (t2 - t1)) * x +
            ((s - t1) / (t2 - t1)) * z);  /// Ensures y remains positive

    if (y > 1.0)  /// Ensure y remains <= 1.0
    {
      y = 1.0 - abs(1.0 - y);
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
          mklj = DrawBridgePMFSameMutationApprox(N_i, x, s - t1, t2 - t1, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFSameMutationOneQApprox(N_i, x, s - t1, t2 - t1, o,
                                                     gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFSameMutation(N_i, x, s - t1, t2 - t1, o, gen);
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
          mklj = DrawBridgePMFDifferentMutationApprox(N_i, s - t1, t2 - t1, x,
                                                      o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFDifferentMutationOneQApprox(N_i, s - t1, t2 - t1,
                                                          x, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj =
              DrawBridgePMFDifferentMutation(N_i, s - t1, t2 - t1, x, o, gen);
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
          mklj = DrawBridgePMFInteriorMutationApprox(N_i, x, z, s - t1, t2 - t1,
                                                     o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(N_i, x, z, s - t1,
                                                         t2 - t1, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj =
              DrawBridgePMFInteriorMutation(N_i, x, z, s - t1, t2 - t1, o, gen);
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
          mklj = DrawBridgePMFDifferentMutationApprox(N_i, s - t1, t2 - t1, x,
                                                      o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFDifferentMutationOneQApprox(N_i, s - t1, t2 - t1,
                                                          x, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj =
              DrawBridgePMFDifferentMutation(N_i, s - t1, t2 - t1, x, o, gen);
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
          mklj = DrawBridgePMFSameMutationApprox(N_i, x, s - t1, t2 - t1, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFSameMutationOneQApprox(N_i, x, s - t1, t2 - t1, o,
                                                     gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMFSameMutation(N_i, x, s - t1, t2 - t1, o, gen);
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
          mklj = DrawBridgePMFInteriorMutationApprox(N_i, x, z, s - t1, t2 - t1,
                                                     o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(N_i, x, z, s - t1,
                                                         t2 - t1, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj =
              DrawBridgePMFInteriorMutation(N_i, x, z, s - t1, t2 - t1, o, gen);
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
        double newt2 = t2 - t1, news = t2 - s;

        if (news <= o.g1984threshold &&
            newt2 - news <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFInteriorMutationApprox(
              N_i, z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time
                     /// bridge
        }            /// One time increment both below threshold
        else if ((news <= o.g1984threshold &&
                  newt2 - news > o.g1984threshold) ||
                 (news > o.g1984threshold &&
                  newt2 - news <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(
              N_i, z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time
                     /// bridge
        } else       /// Time increments are large enough for alternating series
                     /// method
        {
          mklj = DrawBridgePMFInteriorMutation(
              N_i, z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time
                     /// bridge
        }

        m = mklj[0];
        k = mklj[1];
        l = mklj[2];
        j = mklj[3];
      } else if (!(z < 1.0))  /// x in (0,1) & z = 1
      {
        double newt2 = t2 - t1, news = t2 - s;

        if (news <= o.g1984threshold &&
            newt2 - news <=
                o.g1984threshold)  /// Time increments both below threshold
        {
          mklj = DrawBridgePMFInteriorMutationApprox(
              N_i, z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time
                     /// bridge
        }            /// One time increment both below threshold
        else if ((news <= o.g1984threshold &&
                  newt2 - news > o.g1984threshold) ||
                 (news > o.g1984threshold &&
                  newt2 - news <= o.g1984threshold)) {
          mklj = DrawBridgePMFInteriorMutationOneQApprox(
              N_i, z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time
                     /// bridge
        } else       /// Time increments are large enough for alternating series
                     /// method
        {
          mklj = DrawBridgePMFInteriorMutation(
              N_i, z, x, news, newt2, o,
              gen);  /// Time reversal! Flip x and z because reverse time
                     /// bridge
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
          mklj = DrawBridgePMFG1984(N_i, x, z, s - t1, t2 - t1, o, gen);
        }  /// One time increment both below threshold
        else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                 (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
          mklj = DrawBridgePMFOneQApprox(N_i, x, z, s - t1, t2 - t1, o, gen);
        } else  /// Time increments are large enough for alternating series
                /// method
        {
          mklj = DrawBridgePMF(N_i, x, z, s - t1, t2 - t1, o, gen);
        }

        m = mklj[0];
        k = mklj[1];
        l = mklj[2];
        j = mklj[3];
      }
    }

    para1 = static_cast<double>(thetaP[N_i][0] + l + j),
    para2 = static_cast<double>(thetaP[N_i][1] + m - l + k - j);

    boost::random::gamma_distribution<> GAMMA1(para1, 1.0), GAMMA2(para2, 1.0);

    y = -1.0;
    while (!(0.0 < y && y < 1.0)) {
      double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }
  }

  return make_pair(y, coeffcount);
}

pair<double, int> WrightFisher::DrawBridgepointDiffTheta(
    size_t N_i, double x, double z, double t1, double t2, double s,
    const Options &o, boost::random::mt19937 &gen) {
  assert((x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > t1) &&
         (s < t2));
  double y, para1, para2;
  int m, k, j, l, coeffcount = -1;
  vector<int> mklj;
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  if (s <
      changepts[N_i_t]) {  // We need to sample the diffusion at the changepoint
                           // and then simulate the point at time s from the
                           // bridge having mutation parameter thetaP[N_i_s]
    pair<double, int> intermediate =
        DrawBridgepointDiffTheta(N_i, x, z, t1, t2, changepts[N_i_t], o, gen);
    return DrawBridgepoint(N_i, x, intermediate.first, t1, changepts[N_i_t], s,
                           o, gen);
  } else if (s > changepts[N_i_t]) {
    // We need to sample the diffusion at the changepoint
    // and then simulate the point at time s from the
    // bridge having mutation parameter thetaP[N_i_t]
    pair<double, int> intermediate =
        DrawBridgepointDiffTheta(N_i, x, z, t1, t2, changepts[N_i_t], o, gen);
    return DrawBridgepoint(N_i_t, intermediate.first, z, changepts[N_i_t], t2,
                           s, o, gen);
  } else {
    if ((s - t1 <= o.bridgethreshold || t2 - s <= o.bridgethreshold) ||
        (((theta[N_i_s] - 1.0) /
          (exp(0.5 * (theta[N_i_s] - 1.0) * (s - t1)) - 1.0)) +
             ((theta[N_i_t] - 1.0) /
              (exp(0.5 * (theta[N_i_t] - 1.0) * (t2 - s)) - 1.0)) >
         260.0))  /// Use diffusion approximations
    {
      /// Last condition checks that corresponding m and k terms are not too
      /// large
      /// to create computational bottleneck
      double y1 =
          DrawEndpoint(N_i, x, t1, s, o, gen)
              .first;  /// Diffusion approximation for when times are too small
      /// and bridge takes too long to compute

      y = abs((y1 - x) + ((t2 - s) / (t2 - t1)) * x +
              ((s - t1) / (t2 - t1)) * z);  /// Ensures y remains positive

      if (y > 1.0)  /// Ensure y remains <= 1.0
      {
        y = 1.0 - abs(1.0 - y);
      }
    } else  /// Else use bridge simulator
    {
      if (!(x > 0.0)) {
        if (!(z > 0.0) || (!(z < 1.0))) {  // x = 0 && (z = 0 || z = 1)
          if (s - t1 <= o.g1984threshold &&
              t2 - s <=
                  o.g1984threshold)  /// Time increments both below threshold
          {
            mklj = DrawBridgePMFDiffThetaBoundariesApprox(N_i, x, z, s - t1,
                                                          t2 - t1, o, gen);
          }  /// One time increment both below threshold
          else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                   (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFDiffThetaBoundariesOneQApprox(N_i, x, z, s - t1,
                                                              t2 - t1, o, gen);
          } else  /// Time increments are large enough for alternating series
                  /// method
          {
            mklj = DrawBridgePMFDiffThetaBoundaries(N_i, x, z, s - t1, t2 - t1,
                                                    o, gen);
          }
          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = mklj[3];
        } else {  // x = 0 && z in (0,1)
          if (s - t1 <= o.g1984threshold &&
              t2 - s <=
                  o.g1984threshold)  /// Time increments both below threshold
          {
            mklj = DrawBridgePMFDiffThetaInteriorApprox(N_i, x, z, s - t1,
                                                        t2 - t1, o, gen);
          }  /// One time increment both below threshold
          else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                   (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFDiffThetaInteriorOneQApprox(N_i, x, z, s - t1,
                                                            t2 - t1, o, gen);
          } else  /// Time increments are large enough for alternating series
                  /// method
          {
            mklj = DrawBridgePMFDiffThetaInterior(N_i, x, z, s - t1, t2 - t1, o,
                                                  gen);
          }
          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = mklj[3];
        }
      } else if (!(x < 1.0)) {
        if (!(z > 0.0) || (!(z < 1.0))) {  // x = 1 && (z = 0 || z = 1)
          if (s - t1 <= o.g1984threshold &&
              t2 - s <=
                  o.g1984threshold)  /// Time increments both below threshold
          {
            mklj = DrawBridgePMFDiffThetaBoundariesApprox(N_i, x, z, s - t1,
                                                          t2 - t1, o, gen);
          }  /// One time increment both below threshold
          else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                   (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFDiffThetaBoundariesOneQApprox(N_i, x, z, s - t1,
                                                              t2 - t1, o, gen);
          } else  /// Time increments are large enough for alternating series
                  /// method
          {
            mklj = DrawBridgePMFDiffThetaBoundaries(N_i, x, z, s - t1, t2 - t1,
                                                    o, gen);
          }
          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = mklj[3];
        } else {  // x = 1 && z in (0,1)
          if (s - t1 <= o.g1984threshold &&
              t2 - s <=
                  o.g1984threshold)  /// Time increments both below threshold
          {
            mklj = DrawBridgePMFDiffThetaInteriorApprox(N_i, x, z, s - t1,
                                                        t2 - t1, o, gen);
          }  /// One time increment both below threshold
          else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                   (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFDiffThetaInteriorOneQApprox(N_i, x, z, s - t1,
                                                            t2 - t1, o, gen);
          } else  /// Time increments are large enough for alternating series
                  /// method
          {
            mklj = DrawBridgePMFDiffThetaInterior(N_i, x, z, s - t1, t2 - t1, o,
                                                  gen);
          }
          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = mklj[3];
        }
      } else {
        if ((!(z > 0.0)) || (!(z < 1.0))) {  // x in (0,1) && (z = 0 || z = 1)
          if (s - t1 <= o.g1984threshold &&
              t2 - s <=
                  o.g1984threshold)  /// Time increments both below threshold
          {
            mklj = DrawBridgePMFDiffThetaInteriorApprox(N_i, x, z, s - t1,
                                                        t2 - t1, o, gen);
          }  /// One time increment both below threshold
          else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                   (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFDiffThetaInteriorOneQApprox(N_i, x, z, s - t1,
                                                            t2 - t1, o, gen);
          } else  /// Time increments are large enough for alternating series
                  /// method
          {
            mklj = DrawBridgePMFDiffThetaInterior(N_i, x, z, s - t1, t2 - t1, o,
                                                  gen);
          }
          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = mklj[3];
        } else {  // x in (0,1) && z in (0,1)
          if (s - t1 <= o.g1984threshold &&
              t2 - s <=
                  o.g1984threshold)  /// Time increments both below threshold
          {
            mklj = DrawBridgePMFDiffThetaApprox(N_i, x, z, s - t1, t2 - t1, o,
                                                gen);
          }  /// One time increment both below threshold
          else if ((s - t1 <= o.g1984threshold && t2 - s > o.g1984threshold) ||
                   (s - t1 > o.g1984threshold && t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFDiffThetaOneQApprox(N_i, x, z, s - t1, t2 - t1,
                                                    o, gen);
          } else  /// Time increments are large enough for alternating series
                  /// method
          {
            mklj = DrawBridgePMFDiffTheta(N_i, x, z, s - t1, t2 - t1, o, gen);
          }
          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = mklj[3];
        }
      }

      para1 = static_cast<double>(thetaP[N_i][0] + l + j),
      para2 = static_cast<double>(thetaP[N_i][1] + m - l + k - j);
      boost::random::gamma_distribution<> GAMMA1(para1, 1.0),
          GAMMA2(para2, 1.0);
      y = -1.0;
      while (!(0.0 < y && y < 1.0)) {
        double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
        y = G1 / (G1 + G2);
      }
    }

    return make_pair(y, coeffcount);
  }
}

pair<double, int> WrightFisher::DrawUnconditionedBridge(
    size_t N_i, double x, double z, double t1, double t2, double s,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Routine to decide which bridge sampler to invoke for
               /// sampling at time s from neutral bridge diffusion started at
               /// x at time t1, ending at z in time t2, with potential
               /// absorption at any time in between [t1,t2]
{
  assert((x >= 0.0) && (x <= 1.0) && (z >= 0.0) && (z <= 1.0) && (s > t1) &&
         (s < t2));
  double y, para1, para2;
  int m = -1, k = -1, j = -1, l = -1, coeffcount = -1;
  vector<int> mklj;

  if ((s - t1 <= o.bridgethreshold || t2 - s <= o.bridgethreshold) ||
      (((theta[N_i] - 1.0) / (exp(0.5 * theta[N_i] * (s - t1)))) +
           ((theta[N_i] - 1.0) / (exp(0.5 * theta[N_i] * (t2 - s)))) >
       260.0))  /// Use diffusion approximations
  {
    /// Last condition checks that corresponding m and k terms are not too
    /// large
    /// to create computational bottleneck
    double y1 =
        DrawEndpoint(N_i, x, t1, s, o, gen)
            .first;  /// Diffusion approximation for when times are too small
    /// and bridge takes too long to compute

    y = abs((y1 - x) + ((t2 - s) / (t2 - t1)) * x +
            ((s - t1) / (t2 - t1)) * z);  /// Ensures y remains positive

    if (y > 1.0)  /// Ensure y remains <= 1.0
    {
      y = 1.0 - abs(1.0 - y);
    }
  } else {
    if (thetaP[N_i].empty()) {
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
            mklj = DrawBridgePMFUnconditionalApprox(N_i, x, z, s - t1, t2 - t1,
                                                    o, gen);
          } else if ((s - t1 <= o.g1984threshold &&
                      t2 - s > o.g1984threshold) ||
                     (s - t1 > o.g1984threshold &&
                      t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFUnconditionalOneQApprox(N_i, x, z, s - t1,
                                                        t2 - t1, o, gen);
          } else {
            mklj =
                DrawBridgePMFUnconditional(N_i, x, z, s - t1, t2 - t1, o, gen);
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
            mklj = DrawBridgePMFUnconditionalApprox(N_i, x, z, s - t1, t2 - t1,
                                                    o, gen);
          } else if ((s - t1 <= o.g1984threshold &&
                      t2 - s > o.g1984threshold) ||
                     (s - t1 > o.g1984threshold &&
                      t2 - s <= o.g1984threshold)) {
            mklj = DrawBridgePMFUnconditionalOneQApprox(N_i, x, z, s - t1,
                                                        t2 - t1, o, gen);
          } else {
            mklj =
                DrawBridgePMFUnconditional(N_i, x, z, s - t1, t2 - t1, o, gen);
          }

          m = mklj[0];
          k = mklj[1];
          l = mklj[2];
          j = k;

          if (l == m) {
            return make_pair(1.0, coeffcount);
          }
        } else  /// Otherwise the conditions imposed imply the diffusion
                /// cannot be absorbed over [t1,t2], so we can use
                /// Drawbridgepoint
        {
          ThetaResetter(N_i);
          return DrawBridgepoint(N_i, x, z, t1, t2, s, o, gen);
        }
      }
    } else if (!(thetaP[N_i][0] > 0.0) || !(thetaP[N_i][1] > 0.0)) {
      if (!(thetaP[N_i][0] > 0.0)) {
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
              mklj = DrawBridgePMFUnconditionalApprox(N_i, x, z, s - t1,
                                                      t2 - t1, o, gen);
            } else if ((s - t1 <= o.g1984threshold &&
                        t2 - s > o.g1984threshold) ||
                       (s - t1 > o.g1984threshold &&
                        t2 - s <= o.g1984threshold)) {
              mklj = DrawBridgePMFUnconditionalOneQApprox(N_i, x, z, s - t1,
                                                          t2 - t1, o, gen);
            } else {
              mklj = DrawBridgePMFUnconditional(N_i, x, z, s - t1, t2 - t1, o,
                                                gen);
            }

            m = mklj[0];
            k = mklj[1];
            l = mklj[2];
            j = 0;

            if (l == 0) {
              return make_pair(0.0, coeffcount);
            }
          } else {
            ThetaResetter(N_i);
            return DrawBridgepoint(N_i, x, z, t1, t2, s, o, gen);
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
              mklj = DrawBridgePMFUnconditionalApprox(N_i, x, z, s - t1,
                                                      t2 - t1, o, gen);
            } else if ((s - t1 <= o.g1984threshold &&
                        t2 - s > o.g1984threshold) ||
                       (s - t1 > o.g1984threshold &&
                        t2 - s <= o.g1984threshold)) {
              mklj = DrawBridgePMFUnconditionalOneQApprox(N_i, x, z, s - t1,
                                                          t2 - t1, o, gen);
            } else {
              mklj = DrawBridgePMFUnconditional(N_i, x, z, s - t1, t2 - t1, o,
                                                gen);
            }

            m = mklj[0];
            k = mklj[1];
            l = mklj[2];
            j = 0;

            if (l == m) {
              return make_pair(1.0, coeffcount);
            }
          } else {
            ThetaResetter(N_i);
            return DrawBridgepoint(N_i, x, z, t1, t2, s, o, gen);
          }
        }
      }
    } else  /// No absorption can happen, so we are in the same case as in
            /// DrawBridgepoint
    {
      ThetaResetter(N_i);
      return DrawBridgepoint(N_i, x, z, t1, t2, s, o, gen);
    }

    para1 = (thetaP[N_i].empty()
                 ? static_cast<double>(l + j)
                 : (!(thetaP[N_i][0] > 0.0)
                        ? static_cast<double>(l + j)
                        : static_cast<double>(thetaP[N_i][0] + l + j)));
    para2 = (thetaP[N_i].empty()
                 ? static_cast<double>(m - l + k - j)
                 : (!(thetaP[N_i][1] > 0.0)
                        ? static_cast<double>(m - l + k - j)
                        : static_cast<double>(thetaP[N_i][1] + m - l + k - j)));

    boost::random::gamma_distribution<> GAMMA1(para1, 1.0), GAMMA2(para2, 1.0);

    y = -1.0;
    while (!(0.0 < y && y < 1.0)) {
      double G1 = GAMMA1(gen), G2 = GAMMA2(gen);
      y = G1 / (G1 + G2);
    }
  }
  return make_pair(y, coeffcount);
}

/// BRIDGE SIMULATION - NON-NEUTRAL PATHS

vector<vector<double>> WrightFisher::NonNeutralDrawBridge(
    size_t N_i, double x, double t1, double t2, double z, bool Absorption,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Draws of paths from non-neutral WF diffusion bridge started
               /// from x at time t1 and ending at z at time t2
{
  bool accept = false;
  vector<double> paras{phiMin[N_i], phiMax[N_i], phiMax[N_i] - phiMin[N_i]};
  double kapmean = paras[2] * (t2 - t1);  /// Rate of Poisson point process
  boost::random::poisson_distribution<int> kap(static_cast<double>(kapmean));

  boost::random::uniform_01<double> unift,
      unifm;  /// Set up uniform generators for points over [t1,t2] *
  /// [0,phiMax-phiMin] * [0,1]
  vector<vector<double>> ptr;
  int rcount = 0;
  vector<double> kappastore;

  while (!accept)  /// Until all skeleton points get accepted, keep going
  {
    int kappa = kap(gen);  /// Generate kappa ~ Pois
    double small_offset = 1.0e-14;
    vector<double> path, times(kappa), marks(kappa), rejcount;
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
      for (vector<double>::iterator itt = times.begin(), itm = marks.begin();
           itt != times.end(); itt++, itm++) {
        if (itt == times.begin()) {
          if (Absorption) {
            path.push_back(
                DrawUnconditionedBridge(N_i, x, z, t1, t2, *itt, o, gen).first);
          } else {
            path.push_back(
                DrawBridgepoint(N_i, x, z, t1, t2, *itt, o, gen)
                    .first);  /// Generate skeleton points sequentially
          }

          if (Phitilde(N_i, path.back()) - paras[0] >
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
                   times.back())  /// There are more than 2 skeleton points,
                                  /// and we are not at the last one yet
        {
          if (Absorption) {
            path.push_back(DrawUnconditionedBridge(N_i, path.back(), z,
                                                   *(itt - 1), t2, *itt, o, gen)
                               .first);
          } else {
            path.push_back(
                DrawBridgepoint(N_i, path.back(), z, *(itt - 1), t2, *itt, o,
                                gen)
                    .first);  /// Generate skeleton points sequentially
          }

          if (Phitilde(N_i, path.back()) - paras[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop
                     /// and generate a new Poisson point process
          {
            rcount++;
            break;
          }
        } else  /// We are at the last skeleton point
        {
          if (Absorption) {
            path.push_back(
                DrawUnconditionedBridge(N_i, path.back(), z, *(itt - 1), t2,
                                        *itt, o,
                                        gen)
                    .first);  /// Generate skeleton point sequentially
          } else {
            path.push_back(
                DrawBridgepoint(N_i, path.back(), z, *(itt - 1), t2, *itt, o,
                                gen)
                    .first);  /// Generate skeleton point sequentially
          }

          if (Phitilde(N_i, path.back()) - paras[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop
                     /// and generate a new Poisson point process
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

pair<double, vector<vector<double>>>
WrightFisher::NonNeutralDrawBridgeDiffTheta(
    size_t N_i, double x, double t1, double t2, double z, double s,
    const Options &o,
    boost::random::mt19937
        &gen)  /// Draws of paths from non-neutral WF diffusion bridge started
               /// from x at time t1 and ending at z at time t2
{
  bool accept = false;
  size_t N_i_s = N_i, N_i_t = N_i + 1;
  vector<double> paras1{phiMin[N_i_s], phiMax[N_i_s],
                        phiMax[N_i_s] - phiMin[N_i_s]},
      paras2{phiMin[N_i_t], phiMax[N_i_t], phiMax[N_i_t] - phiMin[N_i_t]};
  double kapmean1 = paras1[2] * (s - t1),
         kapmean2 = paras2[2] * (t2 - s);  /// Rate of Poisson point process
  boost::random::poisson_distribution<int> kap1(static_cast<double>(kapmean1)),
      kap2(static_cast<double>(kapmean2));

  boost::random::uniform_01<double> unift1, unift2, unifm1,
      unifm2;  /// Set up uniform generators for points over [t1,t2] *
  /// [0,phiMax-phiMin] * [0,1]

  double Xs;
  int rcount = 0;
  vector<double> kappastore;
  vector<vector<double>> ptr;

  while (!accept)  /// Until all skeleton points get accepted, keep going
  {
    vector<vector<double>> tmp;
    Xs = DrawBridgepointDiffTheta(N_i, x, z, t1, t2, s, o, gen).first;
    bool accept1 = false, accept2 = false;
    int kappa1 = kap1(gen);  /// Generate kappa ~ Pois
    double small_offset = 1.0e-14;
    vector<double> path, times(kappa1), marks(kappa1), rejcount, kapstore;
    auto gent1 = [&t1, &s, &unift1, &small_offset,
                  &gen]()  /// time stamps ~ Unif [t1,s]
    { return (t1 + small_offset + ((s - t1) * unift1(gen))); };
    auto genm1 = [&paras1, &unifm1, &gen]()  /// marks ~ Unif [0,phiMax-phiMin]
    { return (paras1[2] * unifm1(gen)); };
    std::generate(begin(times), end(times), gent1);
    std::generate(begin(marks), end(marks), genm1);
    sortVectorsAscending(times, times,
                         marks);  /// Sort vectors according to timestamps

    if (kappa1 == 0)  /// No skeleton points -> accept
    {
      kapstore.push_back(1.0 * kappa1);
      accept1 = true;
    } else  /// kappa > 0 - generate skeleton points and check them
    {
      for (vector<double>::iterator itt = times.begin(), itm = marks.begin();
           itt != times.end(); itt++, itm++) {
        if (itt == times.begin()) {
          path.push_back(DrawBridgepoint(N_i, x, Xs, t1, s, *itt, o, gen)
                             .first);  /// Generate skeleton points sequentially

          if (Phitilde(N_i, path.back()) - paras1[0] >
              *itm)  /// Test generated point is OK, otherwise can stop and
                     /// generate a new Poisson point process
          {
            rcount++;
            break;
          }

          if (kappa1 == 1)  /// We only needed to generate one skeleton point,
                            /// which we accepted, so we can exit
          {
            kapstore.push_back(1.0 * kappa1);
            accept1 = true;
          }
        } else if (*itt !=
                   times.back())  /// There are more than 2 skeleton points,
                                  /// and we are not at the last one yet
        {
          path.push_back(DrawBridgepoint(N_i, path.back(), Xs, *(itt - 1), s,
                                         *itt, o,
                                         gen)
                             .first);  /// Generate skeleton points sequentially

          if (Phitilde(N_i, path.back()) - paras1[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop
                     /// and generate a new Poisson point process
          {
            rcount++;
            break;
          }
          kapstore.push_back(1.0 * kappa1);
        } else  /// We are at the last skeleton point
        {
          path.push_back(DrawBridgepoint(N_i, path.back(), Xs, *(itt - 1), s,
                                         *itt, o,
                                         gen)
                             .first);  /// Generate skeleton point sequentially

          if (Phitilde(N_i, path.back()) - paras1[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop
                     /// and generate a new Poisson point process
          {
            rcount++;
            break;
          }
          kapstore.push_back(1.0 * kappa1);
          accept1 = true;
        }
      }
    }

    int kappa2 = kap2(gen);  /// Generate kappa ~ Pois
    vector<double> path2, times2(kappa2), marks2(kappa2), rejcount2;
    auto gent2 = [&s, &t2, &unift2, &small_offset,
                  &gen]()  /// time stamps ~ Unif [t1,s]
    { return (s + small_offset + ((t2 - s) * unift2(gen))); };
    auto genm2 = [&paras2, &unifm2, &gen]()  /// marks ~ Unif [0,phiMax-phiMin]
    { return (paras2[2] * unifm2(gen)); };
    std::generate(begin(times2), end(times2), gent2);
    std::generate(begin(marks2), end(marks2), genm2);
    sortVectorsAscending(times2, times2,
                         marks2);  /// Sort vectors according to timestamps

    if (kappa2 == 0)  /// No skeleton points -> accept
    {
      kapstore.push_back(1.0 * kappa2);
      tmp.push_back(path);
      tmp.push_back(times);
      tmp.push_back(marks);
      tmp.push_back(kapstore);
      accept2 = true;
    } else  /// kappa > 0 - generate skeleton points and check them
    {
      for (vector<double>::iterator itt = times2.begin(), itm = marks2.begin();
           itt != times2.end(); itt++, itm++) {
        if (itt == times2.begin()) {
          path.push_back(DrawBridgepoint(N_i + 1, Xs, z, s, t2, *itt, o, gen)
                             .first);  /// Generate skeleton points sequentially

          if (Phitilde(N_i + 1, path.back()) - paras2[0] >
              *itm)  /// Test generated point is OK, otherwise can stop and
                     /// generate a new Poisson point process
          {
            rcount++;
            break;
          }

          if (kappa2 == 1)  /// We only needed to generate one skeleton point,
                            /// which we accepted, so we can exit
          {
            kapstore.push_back(1.0 * kappa2);
            tmp.push_back(path);
            times.insert(times.end(), times2.begin(), times2.end());
            tmp.push_back(times);
            marks.insert(marks.end(), marks2.begin(), marks2.end());
            tmp.push_back(marks2);
            tmp.push_back(kapstore);
            accept2 = true;
          }
        } else if (*itt !=
                   times2.back())  /// There are more than 2 skeleton points,
                                   /// and we are not at the last one yet
        {
          path.push_back(DrawBridgepoint(N_i + 1, path.back(), z, *(itt - 1),
                                         t2, *itt, o,
                                         gen)
                             .first);  /// Generate skeleton points sequentially

          if (Phitilde(N_i + 1, path.back()) - paras2[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop
                     /// and generate a new Poisson point process
          {
            rcount++;
            break;
          }
          kapstore.push_back(1.0 * kappa2);
        } else  /// We are at the last skeleton point
        {
          path.push_back(DrawBridgepoint(N_i + 1, path.back(), z, *(itt - 1),
                                         t2, *itt, o,
                                         gen)
                             .first);  /// Generate skeleton point sequentially

          if (Phitilde(N_i + 1, path.back()) - paras2[0] >
              *itm)  /// Check the generated point is OK, otherwise can stop
                     /// and generate a new Poisson point process
          {
            rcount++;
            break;
          }
          kapstore.push_back(1.0 * kappa2);
          tmp.push_back(path);
          times.insert(times.end(), times2.begin(), times2.end());
          tmp.push_back(times);
          marks.insert(marks.end(), marks2.begin(), marks2.end());
          tmp.push_back(marks);
          tmp.push_back(kapstore);
          accept2 = true;
        }
      }
    }
    accept = (accept1 && accept2);
    if (accept) {
      ptr = tmp;
    }
  }

  return make_pair(Xs, ptr);
}

pair<double, int> WrightFisher::NonNeutralDrawBridgepoint(
    size_t N_i, double x, double t1, double t2, double z, double testT,
    bool Absorption, const Options &o,
    boost::random::mt19937
        &gen)  /// Invoke NonNeutralDrawBridge to generate a whole bridge
               /// path, conditionally on the generated path, generate a
               /// neutral draw at time testT
{
  vector<double> bridgeSection, bridgeTimes;
  bridgeSection.push_back(x);
  bridgeTimes.push_back(t1);

  vector<vector<double>> skeleton =
      NonNeutralDrawBridge(N_i, x, t1, t2, z, Absorption, o,
                           gen);  /// Create skeleton points for bridge
  bridgeSection.insert(bridgeSection.end(), skeleton[0].begin(),
                       skeleton[0].end());
  bridgeTimes.insert(bridgeTimes.end(), skeleton[1].begin(), skeleton[1].end());

  bridgeSection.push_back(z);
  bridgeTimes.push_back(t2);

  vector<double>::iterator timeIt = bridgeTimes.begin(),
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
    return DrawUnconditionedBridge(N_i, *xIt, *(xIt + 1), *timeIt,
                                   *(timeIt + 1), testT, o, gen);
  } else {
    return DrawBridgepoint(N_i, *xIt, *(xIt + 1), *timeIt, *(timeIt + 1), testT,
                           o,
                           gen);  /// Return neutral draw using corresponding
                                  /// endpoints and time increments
  }
}

pair<double, int> WrightFisher::NonNeutralDrawBridgepointDiffTheta(
    size_t N_i, double x, double t1, double t2, double z, double testT,
    const Options &o,
    boost::random::mt19937 &gen)  ///
{
  vector<double> bridgeSection, bridgeTimes;
  bridgeSection.push_back(x);
  bridgeTimes.push_back(t1);

  double s = changepts[getIndex(t1) + 1];
  pair<double, vector<vector<double>>> out =
      NonNeutralDrawBridgeDiffTheta(N_i, x, t1, t2, z, s, o, gen);
  vector<vector<double>> skeleton = out.second;

  bridgeSection.insert(bridgeSection.end(), skeleton[0].begin(),
                       skeleton[0].end());
  bridgeTimes.insert(bridgeTimes.end(), skeleton[1].begin(), skeleton[1].end());

  bridgeSection.push_back(out.first);
  bridgeTimes.push_back(s);

  bridgeSection.push_back(z);
  bridgeTimes.push_back(t2);
  sortVectorsAscending(bridgeTimes, bridgeTimes, bridgeSection);

  vector<double>::iterator timeIt = bridgeTimes.begin(),
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
  if ((*timeIt < s) && (*(timeIt + 1) <= s)) {
    return DrawBridgepoint(N_i, *xIt, *(xIt + 1), *timeIt, *(timeIt + 1), testT,
                           o, gen);
  } else if ((*timeIt >= s) && (*(timeIt + 1) > s)) {
    return DrawBridgepoint(N_i + 1, *xIt, *(xIt + 1), *timeIt, *(timeIt + 1),
                           testT, o, gen);
  } else {
    return DrawBridgepointDiffTheta(
        N_i, *xIt, *(xIt + 1), *timeIt, *(timeIt + 1), testT, o,
        gen);  /// Return neutral draw using corresponding
               /// endpoints and time increments
  }
}

/// SIMULATION RUNNER FUNCTIONS

void WrightFisher::DiffusionRunner(
    int nSim, double x, double startT, double endT, bool Absorption,
    string &Filename, double diffusion_threshold,
    double bridge_threshold)  /// Function to generate specified number of
                              /// diffusion draws
{
  size_t N_i = getIndex(startT);
  std::cout << "You've asked to generate " << nSim
            << " draws from the law of a Wright--Fisher diffusion with the "
               "following properties:"
            << std::endl;
  std::cout << "Start point: " << x << std::endl;
  std::cout << "Start time: " << startT << std::endl;
  std::cout << "Sampling time: " << endT << std::endl;
  if (Absorption) {
    if (thetaP[N_i].empty() ||
        ((thetaP[N_i].front() == 0.0 && thetaP[N_i].back() != 0.0) ||
         (thetaP[N_i].front() != 0.0 && thetaP[N_i].back() == 0.0))) {
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
    if (thetaP[N_i].empty() ||
        ((thetaP[N_i].front() == 0.0 && thetaP[N_i].back() != 0.0) ||
         (thetaP[N_i].front() != 0.0 && thetaP[N_i].back() == 0.0))) {
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

  int nosamples = 1, loader = max(static_cast<int>(floor(0.1 * nSim)), 1),
      loader_count = 1;
  while (nosamples < nSim + 1) {
    if (non_neutral) {
      saveFile << NonNeutralDrawEndpoint(N_i, x, startT, endT, Absorption, o,
                                         WF_gen)
                      .first
               << "\n";
    } else {
      if (Absorption) {
        saveFile << DrawUnconditionedDiffusion(N_i, x, endT - startT, o, WF_gen)
                        .first
                 << "\n";
      } else {
        saveFile << DrawEndpoint(N_i, x, startT, endT, o, WF_gen).first << "\n";
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

void WrightFisher::DiffusionRunnerVector(
    int nSim, vector<double> x, double startT, double endT, bool Absorption,
    string &Filename, double diffusion_threshold,
    double bridge_threshold)  /// Function to generate specified number of
                              /// diffusion draws
{
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);
  size_t N_i = first_ge_index(startT);

  for (vector<double>::iterator x_val = x.begin(); x_val != x.end(); x_val++) {
    int nosamples = 1;
    while (nosamples < nSim + 1) {
      if (non_neutral) {
        saveFile << NonNeutralDrawEndpoint(N_i, *x_val, startT, endT,
                                           Absorption, o, WF_gen)
                        .first
                 << "\t";
      } else {
        if (Absorption) {
          saveFile << DrawUnconditionedDiffusion(N_i, *x_val, endT - startT, o,
                                                 WF_gen)
                          .first
                   << "\t";
        } else {
          saveFile << DrawEndpoint(N_i, *x_val, startT, endT, o, WF_gen).first
                   << "\t";
        }
      }
      nosamples++;
    }
    saveFile << "\n";
  }
}

void WrightFisher::DiffusionTrajectoryVector(
    int nSim, double x, vector<double> times, bool Absorption, string &Filename,
    double diffusion_threshold,
    double bridge_threshold)  /// Function to generate specified number of
                              /// diffusion draws
{
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);

  int nosamples = 1;
  while (nosamples < nSim + 1) {
    double start_x = x, next_val;
    for (vector<double>::iterator t = times.begin(); t != times.end() - 1;
         t++) {
      size_t N_i = getIndex(*t), N_ip1 = getIndex(*(t + 1));
      if (N_i == N_ip1) {
        if (non_neutral) {
          next_val = NonNeutralDrawEndpoint(N_i, start_x, *t, *(t + 1),
                                            Absorption, o, WF_gen)
                         .first;
        } else {
          if (Absorption) {
            next_val = DrawUnconditionedDiffusion(N_i, start_x, *(t + 1) - *t,
                                                  o, WF_gen)
                           .first;
          } else {
            next_val =
                DrawEndpoint(N_i, start_x, *t, *(t + 1), o, WF_gen).first;
          }
        }
        saveFile << next_val << "\t";
        start_x = next_val;
      } else if (N_ip1 - N_i == 1) {
        double intermediate_t = changepts[N_ip1], intermediate_val;
        if (non_neutral) {
          intermediate_val =
              NonNeutralDrawEndpoint(N_i, start_x, *t, intermediate_t,
                                     Absorption, o, WF_gen)
                  .first;
          next_val =
              NonNeutralDrawEndpoint(N_ip1, intermediate_val, intermediate_t,
                                     *(t + 1), Absorption, o, WF_gen)
                  .first;
        } else {
          if (Absorption) {
            intermediate_val = DrawUnconditionedDiffusion(
                                   N_i, start_x, intermediate_t - *t, o, WF_gen)
                                   .first;
            next_val =
                DrawUnconditionedDiffusion(N_ip1, intermediate_val,
                                           *(t + 1) - intermediate_t, o, WF_gen)
                    .first;
          } else {
            intermediate_val =
                DrawEndpoint(N_i, start_x, *t, intermediate_t, o, WF_gen).first;
            next_val = DrawEndpoint(N_ip1, intermediate_val, intermediate_t,
                                    *(t + 1), o, WF_gen)
                           .first;
          }
        }
      } else {
        cerr << "Non-adjacent epoch between sampling times " << *t << " and "
             << *(t + 1) << "!" << endl;
        exit(1);
      }
    }
    nosamples++;
    saveFile << "\n";
  }
}

void WrightFisher::BridgeDiffusionRunner(
    int nSim, double x, double z, double startT, double endT, double sampleT,
    bool Absorption, string &Filename, bool verbose, double diffusion_threshold,
    double bridge_threshold)  /// Function to generate specified number of
                              /// bridge draws
{
  if (verbose) {
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
  }
  size_t N_i_start = getIndex(startT), N_i_end = getIndex(endT),
         N_i_sample = getIndex(sampleT);
  if ((N_i_end < N_i_start) || (N_i_sample < N_i_start) ||
      (N_i_end < N_i_sample)) {
    cerr << "Incorrect start, end and sampling times inserted!" << endl;
    exit(1);
  } else if (N_i_start == N_i_end) {
    if (verbose) {
      cout << "Start and end points both in the same epoch!" << endl;
    }
    if (Absorption) {
      if (thetaP[N_i_start].empty() || ((thetaP[N_i_start].front() == 0.0 &&
                                         thetaP[N_i_start].back() != 0.0) ||
                                        (thetaP[N_i_start].front() != 0.0 &&
                                         thetaP[N_i_start].back() == 0.0))) {
        if (verbose) {
          std::cout << "You have further specified that the diffusion can be "
                       "absorbed at the boundary"
                    << std::endl;
        }
      } else {
        if (verbose) {
          std::cout
              << "You have further specified that the diffusion can be "
                 "absorbed at the boundary, but the provided mutation rates "
                 "are strictly positive, so the diffusion cannot be absorbed"
              << std::endl;
        }
      }
    } else {
      if (thetaP[N_i_start].empty() || ((thetaP[N_i_start].front() == 0.0 &&
                                         thetaP[N_i_start].back() != 0.0) ||
                                        (thetaP[N_i_start].front() != 0.0 &&
                                         thetaP[N_i_start].back() == 0.0))) {
        if (verbose) {
          std::cout
              << "You have further specified that the diffusion cannot be "
                 "absorbed at the boundary"
              << std::endl;
        }
      } else {
        if (verbose) {
          std::cout
              << "You have further specified that the diffusion cannot be "
                 "absorbed at the boundary, but the provided mutation rates "
                 "are strictly positive and thus already ensure this, so the "
                 "resulting draws are coming from the *unconditioned* "
                 "diffusion!"
              << std::endl;
        }
      }
    }
    if (verbose) {
      std::cout << "You've further specified the time threshold for Gaussian "
                   "approximations at "
                << diffusion_threshold
                << ", whilst the bridge approximations threshold was set to "
                << bridge_threshold << std::endl;
      std::cout << "Output will be printed to file in " << Filename
                << std::endl;
    }
    const Options o(diffusion_threshold, bridge_threshold);
    ofstream saveFile;
    saveFile.open(Filename);

    int nosamples = 1, loader = max(static_cast<int>(floor(0.1 * nSim)), 1),
        loader_count = 1;
    while (nosamples < nSim + 1) {
      if (non_neutral) {
        saveFile << NonNeutralDrawBridgepoint(N_i_start, x, startT, endT, z,
                                              sampleT, Absorption, o, WF_gen)
                        .first
                 << "\n";
      } else {
        if (Absorption) {
          saveFile << DrawUnconditionedBridge(N_i_start, x, z, startT, endT,
                                              sampleT, o, WF_gen)
                          .first
                   << "\n";
        } else {
          ThetaResetter(N_i_start);
          saveFile << DrawBridgepoint(N_i_start, x, z, startT, endT, sampleT, o,
                                      WF_gen)
                          .first
                   << "\n";
        }
      }

      nosamples++;
      if ((nosamples % static_cast<int>((loader * loader_count)) == 0) &&
          (verbose)) {
        std::cout << "Simulated " << nosamples << " samples." << endl;
        loader_count++;
      }
    }
    if (verbose) {
      std::cout << "Bridge simulation complete." << endl;
    }
  } else if (N_i_end - N_i_start == 1) {
    if (verbose) {
      std::cout << "The start and end point lie in adjacent epochs!"
                << std::endl;
      std::cout << "You've specified the time threshold for Gaussian "
                   "approximations at "
                << diffusion_threshold
                << ", whilst the bridge approximations threshold was set to "
                << bridge_threshold << std::endl;
      std::cout << "Output will be printed to file in " << Filename
                << std::endl;
    }
    const Options o(diffusion_threshold, bridge_threshold);
    ofstream saveFile;
    saveFile.open(Filename);

    int nosamples = 1, loader = max(static_cast<int>(floor(0.1 * nSim)), 1),
        loader_count = 1;
    while (nosamples < nSim + 1) {
      if (non_neutral) {
        saveFile << NonNeutralDrawBridgepointDiffTheta(
                        N_i_start, x, startT, endT, z, sampleT, o, WF_gen)
                        .first
                 << "\n";
      } else {
        saveFile << DrawBridgepointDiffTheta(N_i_start, x, z, startT, endT,
                                             sampleT, o, WF_gen)
                        .first
                 << "\n";
      }

      nosamples++;
      if ((nosamples % static_cast<int>((loader * loader_count)) == 0) &&
          (verbose)) {
        std::cout << "Simulated " << nosamples << " samples." << endl;
        loader_count++;
      }
    }
    if (verbose) {
      std::cout << "Bridge simulation complete." << endl;
    }
  } else {
    cerr << "The chosen start and end point are in non-adjacent epochs!"
         << endl;
    exit(1);
  }
}

void WrightFisher::DiffusionDensityCalculator(
    int meshSize, double x, double startT, double endT, bool Absorption,
    string &Filename, bool verbose, double diffusion_threshold,
    double bridge_threshold)  /// Function to compute truncated diffusion
                              /// transition density
{
  size_t N_i = getIndex(startT);
  if (verbose) {
    std::cout << "You've asked to compute the transition density of a "
                 "Wright--Fisher diffusion with the "
                 "following properties:"
              << std::endl;
    std::cout << "Start point: " << x << std::endl;
    std::cout << "Start time: " << startT << std::endl;
    std::cout << "Sampling time: " << endT << std::endl;
  }
  if (Absorption) {
    if (thetaP[N_i].empty() ||
        ((thetaP[N_i].front() == 0.0 && thetaP[N_i].back() != 0.0) ||
         (thetaP[N_i].front() != 0.0 && thetaP[N_i].back() == 0.0))) {
      if (verbose) {
        std::cout << "You have further specified that the diffusion can be "
                     "absorbed at the boundary"
                  << std::endl;
      }
    } else {
      if (verbose) {
        std::cout
            << "You have further specified that the diffusion can be "
               "absorbed at the boundary, but the provided mutation rates "
               "are strictly positive, so the diffusion cannot be absorbed"
            << std::endl;
      }
    }
  } else {
    if (thetaP[N_i].empty() ||
        ((thetaP[N_i].front() == 0.0 && thetaP[N_i].back() != 0.0) ||
         (thetaP[N_i].front() != 0.0 && thetaP[N_i].back() == 0.0))) {
      if (verbose) {
        std::cout << "You have further specified that the diffusion cannot be "
                     "absorbed at the boundary"
                  << std::endl;
      }
    } else {
      if (verbose) {
        std::cout
            << "You have further specified that the diffusion cannot be "
               "absorbed at the boundary, but the provided mutation rates "
               "are strictly positive and thus already ensure this, so the "
               "resulting draws are coming from the *unconditioned* diffusion!"
            << std::endl;
      }
    }
  }
  if (verbose) {
    std::cout << "You've further specified the time threshold for Gaussian "
                 "approximations at "
              << diffusion_threshold << std::endl;
    std::cout << "The pointwise computation will be performed over a mesh "
                 "consisting of "
              << meshSize << " equally spaced intervals on [0,1]" << std::endl;
    std::cout << "Output will be printed to file in " << Filename << std::endl;
  }
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);

  int counter = 0;
  double timeInc = endT - startT, yinc = 1.0 / static_cast<double>(meshSize), y,
         ycount = 0.1;
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
                   << UnconditionedDiffusionDensity(N_i, x, y, timeInc, o)
                   << "\n";
        } else {
          if (!(x > y) && !(x < y)) {
            saveFile << y << " " << 1.0 << "\n";
          } else {
            saveFile << y << " " << 0.0 << "\n";
          }
        }
      } else {
        saveFile << y << " "
                 << DiffusionDensityApproximation(N_i, x, y, timeInc, o)
                 << "\n";
      }

      if ((y >= ycount) && (verbose)) {
        std::cout << "Calculated density up to y = " << ycount << endl;
        ycount += 0.1;
      }
      counter++;
    }
    if (verbose) {
      std::cout << "Density calculation complete." << endl;
    }

    saveFile.close();
  }
}

void WrightFisher::BridgeDiffusionDensityCalculator(
    int meshSize, double x, double z, double startT, double endT,
    double sampleT, bool Absorption, string &Filename, bool verbose,
    double diffusion_threshold,
    double bridge_threshold)  /// Function to compute truncated diffusion
                              /// bridge transition density
{
  size_t N_i = getIndex(startT);
  size_t N_i_start = getIndex(startT), N_i_end = getIndex(endT),
         N_i_sample = getIndex(sampleT);
  ofstream saveFile;
  saveFile.open(Filename);
  if (verbose) {
    std::cout << "You've asked to compute the transition density of a "
                 "Wright--Fisher diffusion bridge with the "
                 "following properties:"
              << std::endl;
    std::cout << "Start point: " << x << std::endl;
    std::cout << "Start time: " << startT << std::endl;
    std::cout << "End point: " << z << std::endl;
    std::cout << "End time: " << endT << std::endl;
    std::cout << "Sampling time: " << sampleT << std::endl;
  }
  if ((N_i_end < N_i_start) || (N_i_sample < N_i_start) ||
      (N_i_end < N_i_sample)) {
    cerr << "Incorrect start, end and sampling times inserted!" << endl;
    exit(1);
  } else if (N_i_start == N_i_end) {
    if (Absorption) {
      if (thetaP[N_i].empty() ||
          ((thetaP[N_i].front() == 0.0 && thetaP[N_i].back() != 0.0) ||
           (thetaP[N_i].front() != 0.0 && thetaP[N_i].back() == 0.0))) {
        if (verbose) {
          std::cout << "You have further specified that the diffusion can be "
                       "absorbed at the boundary"
                    << std::endl;
        }
      } else {
        if (verbose) {
          std::cout
              << "You have further specified that the diffusion can be "
                 "absorbed at the boundary, but the provided mutation rates "
                 "are strictly positive, so the diffusion cannot be absorbed"
              << std::endl;
        }
      }
    } else {
      if (thetaP[N_i].empty() ||
          ((thetaP[N_i].front() == 0.0 && thetaP[N_i].back() != 0.0) ||
           (thetaP[N_i].front() != 0.0 && thetaP[N_i].back() == 0.0))) {
        if (verbose) {
          std::cout
              << "You have further specified that the diffusion cannot be "
                 "absorbed at the boundary"
              << std::endl;
        }
      } else {
        if (verbose) {
          std::cout
              << "You have further specified that the diffusion cannot be "
                 "absorbed at the boundary, but the provided mutation rates "
                 "are strictly positive and thus already ensure this, so the "
                 "resulting draws are coming from the *unconditioned* "
                 "diffusion!"
              << std::endl;
        }
      }
    }
    if (verbose) {
      std::cout << "You've further specified the time threshold for Gaussian "
                   "approximations at "
                << diffusion_threshold
                << ", whilst the bridge approximations threshold was set to "
                << bridge_threshold << std::endl;
      std::cout << "The pointwise computation will be performed over a mesh "
                   "consisting of "
                << meshSize << "equally spaced intervals on [0,1]" << std::endl;
      std::cout << "Output will be printed to file in " << Filename
                << std::endl;
    }
    const Options o(diffusion_threshold, bridge_threshold);

    int counter = 0;
    double timeInc = endT - startT, sampleInc = sampleT - startT,
           yinc = 1.0 / static_cast<double>(meshSize), y, ycount = 0.1;
    while (counter <= meshSize) {
      if (counter == 0) {
        y = (thetaP[N_i][0] >= 1.0) ? 0.0 : 1.0e-3;
      } else if (counter == meshSize) {
        y = (thetaP[N_i][1] >= 1.0) ? 1.0 : 1.0 - 1.0e-3;
        ;
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
                     << UnconditionedBridgeDensity(N_i, x, z, y, sampleInc,
                                                   timeInc, o)
                     << "\n";
          } else {
            if (!(x > y) && !(x < y)) {
              saveFile << y << " " << 1.0 << "\n";
            } else {
              saveFile << y << " " << 0.0 << "\n";
            }
          }
        } else {
          saveFile << y << " "
                   << BridgeDensity(N_i, x, z, y, sampleInc, timeInc, o)
                   << "\n";
        }
      }

      if ((y >= ycount) && (verbose)) {
        std::cout << "Calculated density up to y = " << ycount << endl;
        ycount += 0.1;
      }
      counter++;
    }
  } else if (N_i_end - N_i_start == 1) {
    if (verbose) {
      std::cout << "You've further specified the time threshold for Gaussian "
                   "approximations at "
                << diffusion_threshold
                << ", whilst the bridge approximations threshold was set to "
                << bridge_threshold << std::endl;
      std::cout << "The pointwise computation will be performed over a mesh "
                   "consisting of "
                << meshSize << "equally spaced intervals on [0,1]" << std::endl;
      std::cout << "Output will be printed to file in " << Filename
                << std::endl;
    }

    const Options o(diffusion_threshold, bridge_threshold);
    ofstream saveFile;
    saveFile.open(Filename);

    int counter = 0;
    double timeInc = endT - startT, sampleInc = sampleT - startT,
           yinc = 1.0 / static_cast<double>(meshSize), y, ycount = 0.1;
    while (counter <= meshSize) {
      if (counter == 0) {
        y = (thetaP[N_i][0] >= 1.0) ? 0.0 : 1.0e-3;
        ;
      } else if (counter == meshSize) {
        y = (thetaP[N_i][1] >= 1.0) ? 1.0 : 1.0 - 1.0e-3;
        ;
      } else {
        y += yinc;
      }
      if (non_neutral) {
        cerr << "Truncated density cannot be computed for non-neutral case due "
                "to presence of intractable quantities!";
      } else {
        saveFile << y << " "
                 << BridgeDiffThetaDensity(N_i, x, z, y, sampleInc, timeInc, o)
                 << "\n";
      }

      if ((y >= ycount) && (verbose)) {
        std::cout << "Calculated density up to y = " << ycount << endl;
        ycount += 0.1;
      }
      counter++;
    }
  } else {
    cerr << "The chosen start and end point are in non-adjacent epochs!"
         << endl;
    exit(1);
  }

  if (verbose) {
    std::cout << "Density calculation complete." << endl;
  }

  saveFile.close();
}

void WrightFisher::DrawBridgeDiffTheta(string Filename, int nSim, size_t N_i,
                                       double x, double z, double startT,
                                       double endT, double sampleT,
                                       double diffusion_threshold,
                                       double bridge_threshold) {
  const Options o(diffusion_threshold, bridge_threshold);
  ofstream saveFile;
  saveFile.open(Filename);
  int nosamples = 1, loader = max(static_cast<int>(floor(0.1 * nSim)), 1),
      loader_count = 1;
  while (nosamples < nSim + 1) {
    if (non_neutral) {
      saveFile << NonNeutralDrawBridgepointDiffTheta(N_i, x, z, startT, endT,
                                                     sampleT, o, WF_gen)
                      .first
               << "\n";
    } else {
      saveFile << DrawBridgepointDiffTheta(N_i, x, z, startT, endT, sampleT, o,
                                           WF_gen)
                      .first
               << "\n";
    }
    nosamples++;
    if (nosamples % (loader * loader_count) == 0) {
      std::cout << "Simulated " << nosamples << " samples." << endl;
      loader_count++;
    }
  }

  saveFile.close();
}