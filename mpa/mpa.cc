/**
 * Moon Position Algorithm (MPA)
 *
 * Based on SAMPA (Solar And Moon Position Algorithm) and SPA (Solar Position
 * Algorithm), downloaded from the US Department of Energy, stripped to the
 * barest essentials for the moon, and ported to C++.
 *
 * SPA and SAMPA were written in 2003 and 2012, respectively, by Afshin Michael
 * Andreas (Afshin.Andreas@NREL.gov), Measurement & Instrumentation Team, Solar
 * Radiation Research Laboratory, National Renewable Energy Laboratory (NREL).
 * They were based on the NREL technical report "Solar Position Algorithm for
 * Solar Radiation Application" by Ibrahim Reda & Afshin Andreas, and the NREL
 * technical report "Solar Eclipse Monitoring for Solar Energy Applications
 * using the Solar and Moon Position Algorithms" by Ibrahim Reda.
 *
 * Ported by Joey Parrish <joey.parrish@gmail.com>.
 *
 * Original copyright notices:
 *
 * Copyright (C) 2008-2012 Alliance for Sustainable Energy, LLC, All Rights
 * Reserved.
 *
 * The Solar Position Algorithm ("Software") is code in development prepared by
 * employees of the Alliance for Sustainable Energy, LLC, (hereinafter the
 * "Contractor"), under Contract No. DE-AC36-08GO28308 ("Contract") with the
 * U.S. Department of Energy (the "DOE"). The United States Government has been
 * granted for itself and others acting on its behalf a paid-up, non-
 * exclusive, irrevocable, worldwide license in the Software to reproduce,
 * prepare derivative works, and perform publicly and display publicly.
 * Beginning five (5) years after the date permission to assert copyright is
 * obtained from the DOE, and subject to any subsequent five (5) year renewals,
 * the United States Government is granted for itself and others acting on its
 * behalf a paid-up, non-exclusive, irrevocable, worldwide license in the
 * Software to reproduce, prepare derivative works, distribute copies to the
 * public, perform publicly and display publicly, and to permit others to do
 * so. If the Contractor ceases to make this computer software available, it
 * may be obtained from DOE's Office of Scientific and Technical Information's
 * Energy Science and Technology Software Center (ESTSC) at P.O. Box 1020, Oak
 * Ridge, TN 37831-1020.
 *
 * THIS SOFTWARE IS PROVIDED BY THE CONTRACTOR "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE CONTRACTOR OR THE U.S. GOVERNMENT BE LIABLE FOR ANY SPECIAL,
 * INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER, INCLUDING BUT
 * NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF DATA OR PROFITS, WHICH MAY
 * RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR OTHER TORTIOUS CLAIM THAT
 * ARISES OUT OF OR IN CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS
 * SOFTWARE.
 *
 * The Software is being provided for internal, noncommercial purposes only and
 * shall not be re-distributed. Please contact Jennifer Ramsey
 * (Jennifer.Ramsey@nrel.gov) in the NREL Commercialization and Technology
 * Transfer Office for information concerning a commercial license to use the
 * Software.
 *
 * As a condition of using the Software in an application, the developer of the
 * application agrees to reference the use of the Software and make this Notice
 * readily accessible to any end-user in a Help|About screen or equivalent
 * manner.
 */

#include "mpa.h"

#include <cmath>
#include <cstring>

#define COUNT 60
#define Y_COUNT 63
#define SUN_RADIUS 0.26667

namespace mpa {

namespace {

enum { TERM_X0, TERM_X1, TERM_X2, TERM_X3, TERM_X4, SPA_X_TERM_COUNT };
enum { TERM_PSI_A, TERM_PSI_B, TERM_EPS_C, TERM_EPS_D, SPA_PE_TERM_COUNT };
enum { TERM_D, TERM_M, TERM_MPR, TERM_F, TERM_LB, TERM_R, MPA_TERM_COUNT };

//  Periodic Terms for the nutation in longitude and obliquity
const int Y_TERMS[Y_COUNT][SPA_X_TERM_COUNT] = {
    {0, 0, 0, 0, 1},   {-2, 0, 0, 2, 2},  {0, 0, 0, 2, 2},   {0, 0, 0, 0, 2},
    {0, 1, 0, 0, 0},   {0, 0, 1, 0, 0},   {-2, 1, 0, 2, 2},  {0, 0, 0, 2, 1},
    {0, 0, 1, 2, 2},   {-2, -1, 0, 2, 2}, {-2, 0, 1, 0, 0},  {-2, 0, 0, 2, 1},
    {0, 0, -1, 2, 2},  {2, 0, 0, 0, 0},   {0, 0, 1, 0, 1},   {2, 0, -1, 2, 2},
    {0, 0, -1, 0, 1},  {0, 0, 1, 2, 1},   {-2, 0, 2, 0, 0},  {0, 0, -2, 2, 1},
    {2, 0, 0, 2, 2},   {0, 0, 2, 2, 2},   {0, 0, 2, 0, 0},   {-2, 0, 1, 2, 2},
    {0, 0, 0, 2, 0},   {-2, 0, 0, 2, 0},  {0, 0, -1, 2, 1},  {0, 2, 0, 0, 0},
    {2, 0, -1, 0, 1},  {-2, 2, 0, 2, 2},  {0, 1, 0, 0, 1},   {-2, 0, 1, 0, 1},
    {0, -1, 0, 0, 1},  {0, 0, 2, -2, 0},  {2, 0, -1, 2, 1},  {2, 0, 1, 2, 2},
    {0, 1, 0, 2, 2},   {-2, 1, 1, 0, 0},  {0, -1, 0, 2, 2},  {2, 0, 0, 2, 1},
    {2, 0, 1, 0, 0},   {-2, 0, 2, 2, 2},  {-2, 0, 1, 2, 1},  {2, 0, -2, 0, 1},
    {2, 0, 0, 0, 1},   {0, -1, 1, 0, 0},  {-2, -1, 0, 2, 1}, {-2, 0, 0, 0, 1},
    {0, 0, 2, 2, 1},   {-2, 0, 2, 0, 1},  {-2, 1, 0, 2, 1},  {0, 0, 1, -2, 0},
    {-1, 0, 1, 0, 0},  {-2, 1, 0, 0, 0},  {1, 0, 0, 0, 0},   {0, 0, 1, 2, 0},
    {0, 0, -2, 2, 2},  {-1, -1, 1, 0, 0}, {0, 1, 1, 0, 0},   {0, -1, 1, 2, 2},
    {2, -1, -1, 2, 2}, {0, 0, 3, 2, 2},   {2, -1, 0, 2, 2},
};

const double PE_TERMS[Y_COUNT][SPA_PE_TERM_COUNT] = {
    {-171996, -174.2, 92025, 8.9},
    {-13187, -1.6, 5736, -3.1},
    {-2274, -0.2, 977, -0.5},
    {2062, 0.2, -895, 0.5},
    {1426, -3.4, 54, -0.1},
    {712, 0.1, -7, 0},
    {-517, 1.2, 224, -0.6},
    {-386, -0.4, 200, 0},
    {-301, 0, 129, -0.1},
    {217, -0.5, -95, 0.3},
    {-158, 0, 0, 0},
    {129, 0.1, -70, 0},
    {123, 0, -53, 0},
    {63, 0, 0, 0},
    {63, 0.1, -33, 0},
    {-59, 0, 26, 0},
    {-58, -0.1, 32, 0},
    {-51, 0, 27, 0},
    {48, 0, 0, 0},
    {46, 0, -24, 0},
    {-38, 0, 16, 0},
    {-31, 0, 13, 0},
    {29, 0, 0, 0},
    {29, 0, -12, 0},
    {26, 0, 0, 0},
    {-22, 0, 0, 0},
    {21, 0, -10, 0},
    {17, -0.1, 0, 0},
    {16, 0, -8, 0},
    {-16, 0.1, 7, 0},
    {-15, 0, 9, 0},
    {-13, 0, 7, 0},
    {-12, 0, 6, 0},
    {11, 0, 0, 0},
    {-10, 0, 5, 0},
    {-8, 0, 3, 0},
    {7, 0, -3, 0},
    {-7, 0, 0, 0},
    {-7, 0, 3, 0},
    {-7, 0, 3, 0},
    {6, 0, 0, 0},
    {6, 0, -3, 0},
    {6, 0, -3, 0},
    {-6, 0, 3, 0},
    {-6, 0, 3, 0},
    {5, 0, 0, 0},
    {-5, 0, 3, 0},
    {-5, 0, 3, 0},
    {-5, 0, 3, 0},
    {4, 0, 0, 0},
    {4, 0, 0, 0},
    {4, 0, 0, 0},
    {-4, 0, 0, 0},
    {-4, 0, 0, 0},
    {-4, 0, 0, 0},
    {3, 0, 0, 0},
    {-3, 0, 0, 0},
    {-3, 0, 0, 0},
    {-3, 0, 0, 0},
    {-3, 0, 0, 0},
    {-3, 0, 0, 0},
    {-3, 0, 0, 0},
    {-3, 0, 0, 0},
};

// Moon's Periodic Terms for Longitude and Distance
const double ML_TERMS[COUNT][MPA_TERM_COUNT] = {
    {0, 0, 1, 0, 6288774, -20905355},
    {2, 0, -1, 0, 1274027, -3699111},
    {2, 0, 0, 0, 658314, -2955968},
    {0, 0, 2, 0, 213618, -569925},
    {0, 1, 0, 0, -185116, 48888},
    {0, 0, 0, 2, -114332, -3149},
    {2, 0, -2, 0, 58793, 246158},
    {2, -1, -1, 0, 57066, -152138},
    {2, 0, 1, 0, 53322, -170733},
    {2, -1, 0, 0, 45758, -204586},
    {0, 1, -1, 0, -40923, -129620},
    {1, 0, 0, 0, -34720, 108743},
    {0, 1, 1, 0, -30383, 104755},
    {2, 0, 0, -2, 15327, 10321},
    {0, 0, 1, 2, -12528, 0},
    {0, 0, 1, -2, 10980, 79661},
    {4, 0, -1, 0, 10675, -34782},
    {0, 0, 3, 0, 10034, -23210},
    {4, 0, -2, 0, 8548, -21636},
    {2, 1, -1, 0, -7888, 24208},
    {2, 1, 0, 0, -6766, 30824},
    {1, 0, -1, 0, -5163, -8379},
    {1, 1, 0, 0, 4987, -16675},
    {2, -1, 1, 0, 4036, -12831},
    {2, 0, 2, 0, 3994, -10445},
    {4, 0, 0, 0, 3861, -11650},
    {2, 0, -3, 0, 3665, 14403},
    {0, 1, -2, 0, -2689, -7003},
    {2, 0, -1, 2, -2602, 0},
    {2, -1, -2, 0, 2390, 10056},
    {1, 0, 1, 0, -2348, 6322},
    {2, -2, 0, 0, 2236, -9884},
    {0, 1, 2, 0, -2120, 5751},
    {0, 2, 0, 0, -2069, 0},
    {2, -2, -1, 0, 2048, -4950},
    {2, 0, 1, -2, -1773, 4130},
    {2, 0, 0, 2, -1595, 0},
    {4, -1, -1, 0, 1215, -3958},
    {0, 0, 2, 2, -1110, 0},
    {3, 0, -1, 0, -892, 3258},
    {2, 1, 1, 0, -810, 2616},
    {4, -1, -2, 0, 759, -1897},
    {0, 2, -1, 0, -713, -2117},
    {2, 2, -1, 0, -700, 2354},
    {2, 1, -2, 0, 691, 0},
    {2, -1, 0, -2, 596, 0},
    {4, 0, 1, 0, 549, -1423},
    {0, 0, 4, 0, 537, -1117},
    {4, -1, 0, 0, 520, -1571},
    {1, 0, -2, 0, -487, -1739},
    {2, 1, 0, -2, -399, 0},
    {0, 0, 2, -2, -381, -4421},
    {1, 1, 1, 0, 351, 0},
    {3, 0, -2, 0, -340, 0},
    {4, 0, -3, 0, 330, 0},
    {2, -1, 2, 0, 327, 0},
    {0, 2, 1, 0, -323, 1165},
    {1, 1, -1, 0, 299, 0},
    {2, 0, 3, 0, 294, 0},
    {2, 0, -1, -2, 0, 8752}};

// Moon's Periodic Terms for Latitude
const double MB_TERMS[COUNT][MPA_TERM_COUNT] = {
    {0, 0, 0, 1, 5128122, 0}, {0, 0, 1, 1, 280602, 0},
    {0, 0, 1, -1, 277693, 0}, {2, 0, 0, -1, 173237, 0},
    {2, 0, -1, 1, 55413, 0},  {2, 0, -1, -1, 46271, 0},
    {2, 0, 0, 1, 32573, 0},   {0, 0, 2, 1, 17198, 0},
    {2, 0, 1, -1, 9266, 0},   {0, 0, 2, -1, 8822, 0},
    {2, -1, 0, -1, 8216, 0},  {2, 0, -2, -1, 4324, 0},
    {2, 0, 1, 1, 4200, 0},    {2, 1, 0, -1, -3359, 0},
    {2, -1, -1, 1, 2463, 0},  {2, -1, 0, 1, 2211, 0},
    {2, -1, -1, -1, 2065, 0}, {0, 1, -1, -1, -1870, 0},
    {4, 0, -1, -1, 1828, 0},  {0, 1, 0, 1, -1794, 0},
    {0, 0, 0, 3, -1749, 0},   {0, 1, -1, 1, -1565, 0},
    {1, 0, 0, 1, -1491, 0},   {0, 1, 1, 1, -1475, 0},
    {0, 1, 1, -1, -1410, 0},  {0, 1, 0, -1, -1344, 0},
    {1, 0, 0, -1, -1335, 0},  {0, 0, 3, 1, 1107, 0},
    {4, 0, 0, -1, 1021, 0},   {4, 0, -1, 1, 833, 0},
    {0, 0, 1, -3, 777, 0},    {4, 0, -2, 1, 671, 0},
    {2, 0, 0, -3, 607, 0},    {2, 0, 2, -1, 596, 0},
    {2, -1, 1, -1, 491, 0},   {2, 0, -2, 1, -451, 0},
    {0, 0, 3, -1, 439, 0},    {2, 0, 2, 1, 422, 0},
    {2, 0, -3, -1, 421, 0},   {2, 1, -1, 1, -366, 0},
    {2, 1, 0, 1, -351, 0},    {4, 0, 0, 1, 331, 0},
    {2, -1, 1, 1, 315, 0},    {2, -2, 0, -1, 302, 0},
    {0, 0, 1, 3, -283, 0},    {2, 1, 1, -1, -229, 0},
    {1, 1, 0, -1, 223, 0},    {1, 1, 0, 1, 223, 0},
    {0, 1, -2, -1, -220, 0},  {2, 1, -1, -1, -220, 0},
    {1, 0, 1, 1, -185, 0},    {2, -1, -2, -1, 181, 0},
    {0, 1, 2, 1, -177, 0},    {4, 0, -2, -1, 176, 0},
    {4, -1, -1, -1, 166, 0},  {1, 0, 1, -1, -164, 0},
    {4, 0, 1, -1, 132, 0},    {1, 0, -1, -1, -119, 0},
    {4, -1, 0, -1, 115, 0},   {2, -2, 0, 1, 107, 0}};

int integer(double value) { return value; }

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

double julian_day(int year, int month, int day, int hour, int minute,
                  double second) {
  double day_decimal, julian_day, a;

  day_decimal = day + (hour + (minute + second / 60.0) / 60.0) / 24.0;

  if (month < 3) {
    month += 12;
    year--;
  }

  julian_day = integer(365.25 * (year + 4716.0)) +
               integer(30.6001 * (month + 1)) + day_decimal - 1524.5;

  if (julian_day > 2299160.0) {
    a = integer(year / 100);
    julian_day += (2 - a + integer(a / 4));
  }

  return julian_day;
}

double julian_century(double jd) { return (jd - 2451545.0) / 36525.0; }

double julian_millennium(double jc) { return (jc / 10.0); }

double mean_elongation_moon_sun(double jc) {
  return third_order_polynomial(1.0 / 189474.0, -0.0019142, 445267.11148,
                                297.85036, jc);
}

double mean_anomaly_sun(double jc) {
  return third_order_polynomial(-1.0 / 300000.0, -0.0001603, 35999.05034,
                                357.52772, jc);
}

double mean_anomaly_moon(double jc) {
  return third_order_polynomial(1.0 / 56250.0, 0.0086972, 477198.867398,
                                134.96298, jc);
}

double argument_latitude_moon(double jc) {
  return third_order_polynomial(1.0 / 327270.0, -0.0036825, 483202.017538,
                                93.27191, jc);
}

double ascending_longitude_moon(double jc) {
  return third_order_polynomial(1.0 / 450000.0, 0.0020708, -1934.136261,
                                125.04452, jc);
}

double xy_term_summation(int i, double x[SPA_X_TERM_COUNT]) {
  int j;
  double sum = 0;

  for (j = 0; j < SPA_X_TERM_COUNT; j++) sum += x[j] * Y_TERMS[i][j];

  return sum;
}

void nutation_longitude_and_obliquity(double jc, double x[SPA_X_TERM_COUNT],
                                      double *del_psi, double *del_epsilon) {
  int i;
  double xy_term_sum, sum_psi = 0, sum_epsilon = 0;

  for (i = 0; i < Y_COUNT; i++) {
    xy_term_sum = deg2rad(xy_term_summation(i, x));
    sum_psi += (PE_TERMS[i][TERM_PSI_A] + jc * PE_TERMS[i][TERM_PSI_B]) *
               sin(xy_term_sum);
    sum_epsilon += (PE_TERMS[i][TERM_EPS_C] + jc * PE_TERMS[i][TERM_EPS_D]) *
                   cos(xy_term_sum);
  }

  *del_psi = sum_psi / 36000000.0;
  *del_epsilon = sum_epsilon / 36000000.0;
}

double ecliptic_mean_obliquity(double jm) {
  double u = jm / 10.0;

  return 84381.448 +
         u * (-4680.93 +
              u * (-1.55 +
                   u * (1999.25 +
                        u * (-51.38 +
                             u * (-249.67 +
                                  u * (-39.05 +
                                       u * (7.12 +
                                            u * (27.87 +
                                                 u * (5.79 + u * 2.45)))))))));
}

double ecliptic_true_obliquity(double delta_epsilon, double epsilon0) {
  return delta_epsilon + epsilon0 / 3600.0;
}

double greenwich_mean_sidereal_time(double jd, double jc) {
  return limit_degrees(280.46061837 + 360.98564736629 * (jd - 2451545.0) +
                       jc * jc * (0.000387933 - jc / 38710000.0));
}

double greenwich_sidereal_time(double nu0, double delta_psi, double epsilon) {
  return nu0 + delta_psi * cos(deg2rad(epsilon));
}
double fourth_order_polynomial(double a, double b, double c, double d, double e,
                               double x) {
  return (((a * x + b) * x + c) * x + d) * x + e;
}

double moon_mean_longitude(double jc) {
  return limit_degrees(fourth_order_polynomial(-1.0 / 65194000, 1.0 / 538841,
                                               -0.0015786, 481267.88123421,
                                               218.3164477, jc));
}

double moon_mean_elongation(double jc) {
  return limit_degrees(fourth_order_polynomial(-1.0 / 113065000, 1.0 / 545868,
                                               -0.0018819, 445267.1114034,
                                               297.8501921, jc));
}

double sun_mean_anomaly(double jc) {
  return limit_degrees(third_order_polynomial(1.0 / 24490000, -0.0001536,
                                              35999.0502909, 357.5291092, jc));
}

double moon_mean_anomaly(double jc) {
  return limit_degrees(fourth_order_polynomial(-1.0 / 14712000, 1.0 / 69699,
                                               0.0087414, 477198.8675055,
                                               134.9633964, jc));
}

double moon_latitude_argument(double jc) {
  return limit_degrees(fourth_order_polynomial(1.0 / 863310000, -1.0 / 3526000,
                                               -0.0036539, 483202.0175233,
                                               93.2720950, jc));
}

void moon_periodic_term_summation(double d, double m, double m_prime, double f,
                                  double jc,
                                  const double terms[COUNT][MPA_TERM_COUNT],
                                  double *sin_sum, double *cos_sum) {
  int i;
  double e_mult, trig_arg;
  double e = 1.0 - jc * (0.002516 + jc * 0.0000074);

  *sin_sum = 0;
  if (cos_sum != NULL) *cos_sum = 0;
  for (i = 0; i < COUNT; i++) {
    e_mult = pow(e, fabs(terms[i][TERM_M]));
    trig_arg = deg2rad(terms[i][TERM_D] * d + terms[i][TERM_M] * m +
                       terms[i][TERM_F] * f + terms[i][TERM_MPR] * m_prime);
    *sin_sum += e_mult * terms[i][TERM_LB] * sin(trig_arg);
    if (cos_sum != NULL) *cos_sum += e_mult * terms[i][TERM_R] * cos(trig_arg);
  }
}

void moon_longitude_and_latitude(double jc, double l_prime, double f,
                                 double m_prime, double l, double b,
                                 double *lamda_prime, double *beta) {
  double a1 = 119.75 + 131.849 * jc;
  double a2 = 53.09 + 479264.290 * jc;
  double a3 = 313.45 + 481266.484 * jc;
  double delta_l = 3958 * sin(deg2rad(a1)) + 318 * sin(deg2rad(a2)) +
                   1962 * sin(deg2rad(l_prime - f));
  double delta_b = -2235 * sin(deg2rad(l_prime)) + 175 * sin(deg2rad(a1 - f)) +
                   127 * sin(deg2rad(l_prime - m_prime)) +
                   382 * sin(deg2rad(a3)) + 175 * sin(deg2rad(a1 + f)) -
                   115 * sin(deg2rad(l_prime + m_prime));

  *lamda_prime = limit_degrees(l_prime + (l + delta_l) / 1000000);
  *beta = limit_degrees((b + delta_b) / 1000000);
}

double moon_earth_distance(double r) { return 385000.56 + r / 1000; }

double moon_equatorial_horiz_parallax(double delta) {
  return rad2deg(asin(6378.14 / delta));
}

double apparent_moon_longitude(double lamda_prime, double del_psi) {
  return lamda_prime + del_psi;
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

}  // namespace

void compute_mpa(const mpa_input &input, mpa_output *output) {
  // Julian day
  double jd = julian_day(input.year, input.month, input.day, input.hour,
                         input.minute, input.second);
  // Julian century
  double jc = julian_century(jd);
  // Julian millenium
  double jm = julian_millennium(jc);

  double x[SPA_X_TERM_COUNT];
  x[TERM_X0] = mean_elongation_moon_sun(jc);
  x[TERM_X1] = mean_anomaly_sun(jc);
  x[TERM_X2] = mean_anomaly_moon(jc);
  x[TERM_X3] = argument_latitude_moon(jc);
  x[TERM_X4] = ascending_longitude_moon(jc);

  // nutation longitude [degrees]
  double del_psi;
  // nutation obliquity [degrees]
  double del_epsilon;
  nutation_longitude_and_obliquity(jc, x, &del_psi, &del_epsilon);

  // ecliptic mean obliquity [arc seconds]
  double epsilon0 = ecliptic_mean_obliquity(jm);
  double epsilon = ecliptic_true_obliquity(del_epsilon, epsilon0);

  // Greenwich mean sidereal time [degrees]
  double nu0 = greenwich_mean_sidereal_time(jd, jc);
  // Greenwich sidereal time [degrees]
  double nu = greenwich_sidereal_time(nu0, del_psi, epsilon);
  // FIXME: check

  // moon mean longitude [degrees]
  double l_prime = moon_mean_longitude(jc);
  // moon mean elongation [degrees]
  double d = moon_mean_elongation(jc);
  // sun mean anomaly [degrees]
  double m = sun_mean_anomaly(jc);
  // moon mean anomaly [degrees]
  double m_prime = moon_mean_anomaly(jc);
  // moon argument of latitude [degrees]
  double f = moon_latitude_argument(jc);

  double l, r;
  moon_periodic_term_summation(d, m, m_prime, f, jc, ML_TERMS, &l, &r);
  double b;
  moon_periodic_term_summation(d, m, m_prime, f, jc, MB_TERMS, &b, NULL);

  // moon longitude [degrees]
  double lamda_prime;
  // moon latitude [degrees]
  double beta;
  moon_longitude_and_latitude(jc, l_prime, f, m_prime, l, b, &lamda_prime,
                              &beta);

  // distance from earth to moon [kilometers]
  double cap_delta = moon_earth_distance(r);
  // moon equatorial horizontal parallax [degrees]
  double pi = moon_equatorial_horiz_parallax(cap_delta);

  // apparent moon longitude [degrees]
  double lamda = apparent_moon_longitude(lamda_prime, del_psi);

  // geocentric moon right ascension [degrees]
  double alpha = geocentric_right_ascension(lamda, epsilon, beta);
  // geocentric moon declination [degrees]
  double delta = geocentric_declination(beta, epsilon, lamda);

  // observer hour angle [degrees]
  double h = observer_hour_angle(nu, input.longitude, alpha);

  // moon right ascension parallax [degrees]
  double del_alpha;
  // topocentric moon declination [degrees]
  double delta_prime;
  right_ascension_parallax_and_topocentric_dec(
      input.latitude, input.elevation, pi, h, delta, &del_alpha, &delta_prime);
  // topocentric local hour angle [degrees]
  double h_prime = topocentric_local_hour_angle(h, del_alpha);

  // topocentric elevation angle (uncorrected) [degrees]
  double e0 = topocentric_elevation_angle(input.latitude, delta_prime, h_prime);
  // atmospheric refraction correction [degrees]
  double del_e = atmospheric_refraction_correction(
      input.pressure, input.temperature, input.atmos_refract, e0);
  // topocentric elevation angle (corrected) [degrees]
  double e = topocentric_elevation_angle_corrected(e0, del_e);

  output->zenith = topocentric_zenith_angle(e);
  // topocentric azimuth angle (westward from south) [for astronomers]
  double azimuth_astro =
      topocentric_azimuth_angle_astro(h_prime, input.latitude, delta_prime);
  output->azimuth = topocentric_azimuth_angle(azimuth_astro);
}

}  // namespace mpa
