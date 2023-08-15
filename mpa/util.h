#pragma once

namespace mpa {

// Utility functions

double deg2rad(double degrees);
double rad2deg(double radians);
double limit_degrees(double degrees);
double third_order_polynomial(double a, double b, double c, double d, double x);

}  // namespace mpa
