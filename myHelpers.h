#ifndef EA_myHelpers_h
#define EA_myHelpers_h

#include <boost/math/constants/constants.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <cmath>
#include <numeric>
#include <string>

using namespace std;
const double pi = boost::math::constants::pi<double>();

/// typedef boost::multiprecision::cpp_dec_float_100 double100; /// Can switch
/// to this for increased precision, but makes computations much slower!!
typedef double double100;

template <class T1, class T2>
bool EqualFirstCoord(const pair<T1, T2> &P1, const pair<T1, T2> &P2) {
  return P1.first == P2.first;
}

struct Options {
  Options()
      : w(20), g1984threshold(0.08), bridgethreshold(0.008), eps(0.99),
        debug(0){};
  int w; /// w is the field width for outputting
  double100 g1984threshold, bridgethreshold, eps;
  vector<double> thetaP;
  int debug;
};

template <typename T> void printVec(const vector<T> &vec, ostream &o = cout) {
  for (unsigned int i = 0; i < vec.size(); i++) {
    o << vec[i] << " ";
  }
  cout << endl;
}
template <typename T, typename Compare>
void getSortPermutation(std::vector<unsigned> &out, const std::vector<T> &v,
                        Compare compare = std::less<T>()) {
  out.resize(v.size());
  std::iota(out.begin(), out.end(), 0);

  std::sort(out.begin(), out.end(),
            [&](unsigned i, unsigned j) { return compare(v[i], v[j]); });
}

template <typename T>
void applyPermutation(const std::vector<unsigned> &order, std::vector<T> &t) {
  assert(order.size() == t.size());
  std::vector<T> st(t.size());
  for (unsigned i = 0; i < t.size(); i++) {
    st[i] = t[order[i]];
  }
  t = st;
}

template <typename T, typename... S>
void applyPermutation(const std::vector<unsigned> &order, std::vector<T> &t,
                      std::vector<S> &... s) {
  applyPermutation(order, t);
  applyPermutation(order, s...);
}

template <typename T, typename Compare, typename... SS>
void sortVectors(const std::vector<T> &t, Compare comp,
                 std::vector<SS> &... ss) {
  std::vector<unsigned> order;
  getSortPermutation(order, t, comp);
  applyPermutation(order, ss...);
}

template <typename T, typename... SS>
void sortVectorsAscending(const std::vector<T> &t, std::vector<SS> &... ss) {
  sortVectors(t, std::less<T>(), ss...);
}

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

#endif
