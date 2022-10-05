#include "myHelpers.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>

ostream &operator<<(ostream &s, const Options &o) {
  s << "g1984\t" << o.g1984threshold << endl;
  s << "eps\t" << o.eps << endl;
  s << "Debug\t" << o.debug << "\n" << endl;
  return s;
}
