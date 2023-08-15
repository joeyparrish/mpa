#include "util.h"

#include <cmath>

#define PI 3.1415926535897932384626433832795028841971
#define SUN_RADIUS 0.26667

namespace mpa {

double rad2deg(double radians) { return (180.0 / PI) * radians; }

double deg2rad(double degrees) { return (PI / 180.0) * degrees; }

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

double geocentric_right_ascension(double lamda, double epsilon, double beta) {
  double lamda_rad = deg2rad(lamda);
  double epsilon_rad = deg2rad(epsilon);

  return limit_degrees(rad2deg(atan2(
      sin(lamda_rad) * cos(epsilon_rad) - tan(deg2rad(beta)) * sin(epsilon_rad),
      cos(lamda_rad))));
}

double geocentric_declination(double beta, double epsilon, double lamda) {
  double beta_rad = deg2rad(beta);
  double epsilon_rad = deg2rad(epsilon);

  return rad2deg(asin(sin(beta_rad) * cos(epsilon_rad) +
                      cos(beta_rad) * sin(epsilon_rad) * sin(deg2rad(lamda))));
}

double observer_hour_angle(double nu, double longitude, double alpha_deg) {
  return limit_degrees(nu + longitude - alpha_deg);
}

void right_ascension_parallax_and_topocentric_dec(double latitude,
                                                  double elevation, double xi,
                                                  double h, double delta,
                                                  double *delta_alpha,
                                                  double *delta_prime) {
  double delta_alpha_rad;
  double lat_rad = deg2rad(latitude);
  double xi_rad = deg2rad(xi);
  double h_rad = deg2rad(h);
  double delta_rad = deg2rad(delta);
  double u = atan(0.99664719 * tan(lat_rad));
  double y = 0.99664719 * sin(u) + elevation * sin(lat_rad) / 6378140.0;
  double x = cos(u) + elevation * cos(lat_rad) / 6378140.0;

  delta_alpha_rad = atan2(-x * sin(xi_rad) * sin(h_rad),
                          cos(delta_rad) - x * sin(xi_rad) * cos(h_rad));

  *delta_prime =
      rad2deg(atan2((sin(delta_rad) - y * sin(xi_rad)) * cos(delta_alpha_rad),
                    cos(delta_rad) - x * sin(xi_rad) * cos(h_rad)));

  *delta_alpha = rad2deg(delta_alpha_rad);
}

double topocentric_local_hour_angle(double h, double delta_alpha) {
  return h - delta_alpha;
}

double topocentric_elevation_angle(double latitude, double delta_prime,
                                   double h_prime) {
  double lat_rad = deg2rad(latitude);
  double delta_prime_rad = deg2rad(delta_prime);

  return rad2deg(
      asin(sin(lat_rad) * sin(delta_prime_rad) +
           cos(lat_rad) * cos(delta_prime_rad) * cos(deg2rad(h_prime))));
}

double atmospheric_refraction_correction(double pressure, double temperature,
                                         double atmos_refract, double e0) {
  double del_e = 0;

  if (e0 >= -1 * (SUN_RADIUS + atmos_refract))
    del_e = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 1.02 /
            (60.0 * tan(deg2rad(e0 + 10.3 / (e0 + 5.11))));

  return del_e;
}

double topocentric_elevation_angle_corrected(double e0, double delta_e) {
  return e0 + delta_e;
}

double topocentric_zenith_angle(double e) { return 90.0 - e; }

double topocentric_azimuth_angle_astro(double h_prime, double latitude,
                                       double delta_prime) {
  double h_prime_rad = deg2rad(h_prime);
  double lat_rad = deg2rad(latitude);

  return limit_degrees(rad2deg(
      atan2(sin(h_prime_rad), cos(h_prime_rad) * sin(lat_rad) -
                                  tan(deg2rad(delta_prime)) * cos(lat_rad))));
}

double topocentric_azimuth_angle(double azimuth_astro) {
  return limit_degrees(azimuth_astro + 180.0);
}

}  // namespace mpa
