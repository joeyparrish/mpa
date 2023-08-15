#pragma once

namespace mpa {

// Utility functions

double deg2rad(double degrees);
double rad2deg(double radians);
double limit_degrees(double degrees);
double third_order_polynomial(double a, double b, double c, double d, double x);

double geocentric_right_ascension(double lamda, double epsilon, double beta);
double geocentric_declination(double beta, double epsilon, double lamda);
double observer_hour_angle(double nu, double longitude, double alpha_deg);
void right_ascension_parallax_and_topocentric_dec(double latitude,
                                                  double elevation, double xi,
                                                  double h, double delta,
                                                  double *delta_alpha,
                                                  double *delta_prime);
double topocentric_local_hour_angle(double h, double delta_alpha);
double topocentric_elevation_angle(double latitude, double delta_prime,
                                   double h_prime);
double atmospheric_refraction_correction(double pressure, double temperature,
                                         double atmos_refract, double e0);
double topocentric_elevation_angle_corrected(double e0, double delta_e);
double topocentric_zenith_angle(double e);
double topocentric_azimuth_angle_astro(double h_prime, double latitude,
                                       double delta_prime);
double topocentric_azimuth_angle(double azimuth_astro);

}  // namespace mpa
