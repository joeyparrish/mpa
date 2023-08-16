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

#include "spa.h"

#include <cmath>

#include "util.h"

#define Y_COUNT 63

#define TERM_Y_COUNT TERM_X_COUNT

namespace mpa {

namespace {

enum { TERM_X0, TERM_X1, TERM_X2, TERM_X3, TERM_X4, TERM_X_COUNT };
enum { TERM_PSI_A, TERM_PSI_B, TERM_EPS_C, TERM_EPS_D, TERM_PE_COUNT };

//  Periodic Terms for the nutation in longitude and obliquity
const int Y_TERMS[Y_COUNT][TERM_Y_COUNT] = {
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

const double PE_TERMS[Y_COUNT][TERM_PE_COUNT] = {
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

int integer(double value) { return value; }

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

double xy_term_summation(int i, double x[TERM_X_COUNT]) {
  int j;
  double sum = 0;

  for (j = 0; j < TERM_Y_COUNT; j++) sum += x[j] * Y_TERMS[i][j];

  return sum;
}

void nutation_longitude_and_obliquity(double jc, double x[TERM_X_COUNT],
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

}  // namespace

void spa_calculate(const mpa_input& input, spa_data *spa) {
  // Julian day
  double jd = julian_day(input.year, input.month, input.day, input.hour,
                         input.minute, input.second);
  // Julian century
  spa->jc = julian_century(jd);
  // Julian millenium
  double jm = julian_millennium(spa->jc);

  double x[TERM_X_COUNT];
  x[TERM_X0] = mean_elongation_moon_sun(spa->jc);
  x[TERM_X1] = mean_anomaly_sun(spa->jc);
  x[TERM_X2] = mean_anomaly_moon(spa->jc);
  x[TERM_X3] = argument_latitude_moon(spa->jc);
  x[TERM_X4] = ascending_longitude_moon(spa->jc);

  // nutation obliquity [degrees]
  double del_epsilon;
  nutation_longitude_and_obliquity(spa->jc, x, &spa->del_psi,
                                   &del_epsilon);

  // ecliptic mean obliquity [arc seconds]
  double epsilon0 = ecliptic_mean_obliquity(jm);
  // // ecliptic true obliquity  [degrees]
  spa->epsilon = ecliptic_true_obliquity(del_epsilon, epsilon0);

  // Greenwich mean sidereal time [degrees]
  double nu0 = greenwich_mean_sidereal_time(jd, spa->jc);
  // Greenwich sidereal time [degrees]
  spa->nu = greenwich_sidereal_time(nu0, spa->del_psi, spa->epsilon);
}

}  // namespace mpa
