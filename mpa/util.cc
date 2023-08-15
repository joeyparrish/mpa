#include "util.h"

#include <cmath>

namespace mpa {

double rad2deg(double radians) { return (180.0 / M_PI) * radians; }

double deg2rad(double degrees) { return (M_PI / 180.0) * degrees; }

double limit_degrees(double degrees) {
  double limited;

  degrees /= 360.0;
  limited = 360.0 * (degrees - floor(degrees));
  if (limited < 0) limited += 360.0;

  return limited;
}

double third_order_polynomial(double a, double b, double c, double d,
                              double x) {
  return ((a * x + b) * x + c) * x + d;
}

}  // namespace mpa
