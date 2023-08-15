/////////////////////////////////////////////
//          HEADER FILE for SPA.C          //
//                                         //
//      Solar Position Algorithm (SPA)     //
//                   for                   //
//        Solar Radiation Application      //
//                                         //
//               May 12, 2003              //
//                                         //
//   Filename: SPA.H                       //
//                                         //
//   Afshin Michael Andreas                //
//   afshin_andreas@nrel.gov (303)384-6383 //
//                                         //
//   Measurement & Instrumentation Team    //
//   Solar Radiation Research Laboratory   //
//   National Renewable Energy Laboratory  //
//   1617 Cole Blvd, Golden, CO 80401      //
/////////////////////////////////////////////

#pragma once

namespace mpa {

typedef struct {
  //----------------------INPUT VALUES------------------------

  int year;    // 4-digit year,      valid range: -2000 to 6000, error code: 1
  int month;   // 2-digit month,         valid range: 1 to  12,  error code: 2
  int day;     // 2-digit day,           valid range: 1 to  31,  error code: 3
  int hour;    // Observer local hour,   valid range: 0 to  24,  error code: 4
  int minute;  // Observer local minute, valid range: 0 to  59,  error code: 5
  double
      second;  // Observer local second, valid range: 0 to <60,  error code: 6

  double longitude;  // Observer longitude (negative west of Greenwich)
                     // valid range: -180  to  180 degrees, error code: 9

  double latitude;  // Observer latitude (negative south of equator)
                    // valid range: -90   to   90 degrees, error code: 10

  double
      elevation;  // Observer elevation [meters]
                  // valid range: -6500000 or higher meters,    error code: 11

  double pressure;  // Annual average local pressure [millibars]
                    // valid range:    0 to 5000 millibars,       error code: 12

  double
      temperature;  // Annual average local temperature [degrees Celsius]
                    // valid range: -273 to 6000 degrees Celsius, error code; 13

  double atmos_refract;  // Atmospheric refraction at sunrise and sunset (0.5667
                         // deg is typical) valid range: -5   to   5 degrees,
                         // error code: 16

  //-----------------Intermediate OUTPUT VALUES--------------------
  double jd;  // Julian day
  double jc;  // Julian century
  double jm;  // Julian millennium

  double l;  // earth heliocentric longitude [degrees]
  double b;  // earth heliocentric latitude [degrees]
  double r;  // earth radius vector [Astronomical Units, AU]

  double theta;  // geocentric longitude [degrees]
  double beta;   // geocentric latitude [degrees]

  double del_psi;      // nutation longitude [degrees]
  double del_epsilon;  // nutation obliquity [degrees]
  double epsilon0;     // ecliptic mean obliquity [arc seconds]
  double epsilon;      // ecliptic true obliquity  [degrees]

  double del_tau;  // aberration correction [degrees]
  double lamda;    // apparent sun longitude [degrees]
  double nu0;      // Greenwich mean sidereal time [degrees]
  double nu;       // Greenwich sidereal time [degrees]

  double alpha;  // geocentric sun right ascension [degrees]
  double delta;  // geocentric sun declination [degrees]

  double h;            // observer hour angle [degrees]
  double xi;           // sun equatorial horizontal parallax [degrees]
  double del_alpha;    // sun right ascension parallax [degrees]
  double delta_prime;  // topocentric sun declination [degrees]
  double h_prime;      // topocentric local hour angle [degrees]

  double e0;     // topocentric elevation angle (uncorrected) [degrees]
  double del_e;  // atmospheric refraction correction [degrees]
  double e;      // topocentric elevation angle (corrected) [degrees]
  double azimuth_astro;  // topocentric azimuth angle (westward from south) [for
                         // astronomers]

  //---------------------Final OUTPUT VALUES------------------------
  double zenith;         // topocentric zenith angle [degrees]
  double azimuth;        // topocentric azimuth angle (eastward from north) [for
                         // navigators and solar radiation]
} spa_data;

// Utility functions for other applications (such as NREL's SAMPA)
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

// Calculate SPA output values (in structure) based on input values passed in
// structure
int spa_calculate(spa_data *spa);

}  // namespace mpa
