///////////////////////////////////////////////
// Solar and Moon Position Algorithm (SAMPA) //
//                   for                     //
//        Solar Radiation Application        //
//                                           //
//              August 1, 2012               //
//                                           //
//   Filename: SAMPA.C                       //
//                                           //
//   Afshin Michael Andreas                  //
//   Afshin.Andreas@NREL.gov (303)384-6383   //
//                                           //
//   Solar Resource and Forecasting Group    //
//   Solar Radiation Research Laboratory     //
//   National Renewable Energy Laboratory    //
//   15013 Denver W Pkwy, Golden, CO 80401   //
///////////////////////////////////////////////

///////////////////////////////////////////////
//  See the SAMPA.H header file for usage    //
//                                           //
//  This code is based on the NREL           //
//  technical report "Solar Eclipse          //
//  Monitoring for Solar Energy Applications //
//  using the Solar and Moon Position        //
//  Algorithms" by Ibrahim Reda              //
///////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////
//
//   NOTICE
//   Copyright (C) 2012 the Alliance for Sustainable Energy, LLC, All Rights
//   Reserved
//
// This computer software is prepared by the Alliance for Sustainable Energy,
// LLC, (hereinafter the "Contractor"), under Contract DE-AC36-08GO28308
// ("Contract") with the Department of Energy ("DOE"). The United States
// Government has been granted for itself and others acting on its behalf a
// paid-up, non-exclusive, irrevocable, worldwide license in the Software to
// reproduce, prepare derivative works, and perform publicly and display
// publicly. Beginning five (5) years after the date permission to assert
// copyright is obtained from DOE, and subject to any subsequent five (5) year
// renewals, the United States Government is granted for itself and others
// acting on its behalf a paid-up, non-exclusive, irrevocable, worldwide license
// in the Software to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do
// so. If the Contractor ceases to make this computer software available, it may
// be obtained from DOE's Office of Scientific and Technical Information's
// Energy Science and Technology Software Center (ESTSC) at PO Box 1020, Oak
// Ridge, TN 37831-1020. THIS SOFTWARE IS PROVIDED BY THE CONTRACTOR "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE CONTRACTOR, DOE, OR THE U.S GOVERNMENT BE
// LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER, INCLUDING BUT NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF
// DATA OR PROFITS, WHICH MAY RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR
// OTHER TORTIOUS CLAIM THAT ARISES OUT OF OR IN CONNECTION WITH THE ACCESS, USE
// OR PERFORMANCE OF THIS SOFTWARE.
//
// The software is being provided for internal, noncommercial purposes only and
// shall not be re-distributed. Please contact Jennifer Ramsey
// (Jennifer.Ramsey@nrel.gov) in the NREL Commercialization and Technology
// Transfer Office for information concerning a commercial license to use the
// Software.
//
// As a condition of using the software in an application, the developer of the
// application agrees to reference the use of the software and make this notice
// readily accessible to any end-user in a Help|About screen or equivalent
// manner.
//
///////////////////////////////////////////////////////////////////////////////////////////////

#include "mpa.h"

#include <cmath>
#include <cstring>

#include "spa.h"
#include "util.h"

#define COUNT 60

namespace mpa {

namespace {

enum { TERM_D, TERM_M, TERM_MPR, TERM_F, TERM_LB, TERM_R, TERM_COUNT };

struct mpa_intermediate {
  double l_prime;      // moon mean longitude [degrees]
  double d;            // moon mean elongation [degrees]
  double m;            // sun mean anomaly [degrees]
  double m_prime;      // moon mean anomaly [degrees]
  double f;            // moon argument of latitude [degrees]
  double l;            // term l
  double r;            // term r
  double b;            // term b
  double lamda_prime;  // moon longitude [degrees]
  double beta;         // moon latitude [degrees]
  double cap_delta;    // distance from earth to moon [kilometers]
  double pi;           // moon equatorial horizontal parallax [degrees]
  double lamda;        // apparent moon longitude [degrees]

  double alpha;  // geocentric moon right ascension [degrees]
  double delta;  // geocentric moon declination [degrees]

  double h;            // observer hour angle [degrees]
  double del_alpha;    // moon right ascension parallax [degrees]
  double delta_prime;  // topocentric moon declination [degrees]
  double h_prime;      // topocentric local hour angle [degrees]

  double e0;             // topocentric elevation angle (uncorrected) [degrees]
  double del_e;          // atmospheric refraction correction [degrees]
  double e;              // topocentric elevation angle (corrected) [degrees]
  double azimuth_astro;  // topocentric azimuth angle (westward from south) [for
                         // astronomers]
};

///////////////////////////////////////////////////////
///  Moon's Periodic Terms for Longitude and Distance
///////////////////////////////////////////////////////
const double ML_TERMS[COUNT][TERM_COUNT] = {{0, 0, 1, 0, 6288774, -20905355},
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
///////////////////////////////////////////////////////
///  Moon's Periodic Terms for Latitude
///////////////////////////////////////////////////////
const double MB_TERMS[COUNT][TERM_COUNT] = {
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
///////////////////////////////////////////////////////////////////////////////////////////////

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
                                  const double terms[COUNT][TERM_COUNT],
                                  double *sin_sum, double *cos_sum) {
  int i;
  double e_mult, trig_arg;
  double e = 1.0 - jc * (0.002516 + jc * 0.0000074);

  *sin_sum = 0;
  if (cos_sum != 0) *cos_sum = 0;
  for (i = 0; i < COUNT; i++) {
    e_mult = pow(e, fabs(terms[i][TERM_M]));
    trig_arg = deg2rad(terms[i][TERM_D] * d + terms[i][TERM_M] * m +
                       terms[i][TERM_F] * f + terms[i][TERM_MPR] * m_prime);
    *sin_sum += e_mult * terms[i][TERM_LB] * sin(trig_arg);
    if (cos_sum != 0) *cos_sum += e_mult * terms[i][TERM_R] * cos(trig_arg);
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

}  // namespace

void compute_mpa(const mpa_input& input, mpa_output* output) {
  // FIXME: share input structure
  spa_data spa;
  memset(&spa, 0, sizeof(spa));
  spa.year = input.year;
  spa.month = input.month;
  spa.day = input.day;
  spa.hour = input.hour;
  spa.minute = input.minute;
  spa.second = input.second;
  spa.latitude = input.latitude;
  spa.longitude = input.longitude;
  spa.elevation = input.elevation;

  spa_calculate(&spa);

  mpa_intermediate mpa;
  memset(&mpa, 0, sizeof(mpa));

  mpa.l_prime = moon_mean_longitude(spa.jc);
  mpa.d = moon_mean_elongation(spa.jc);
  mpa.m = sun_mean_anomaly(spa.jc);
  mpa.m_prime = moon_mean_anomaly(spa.jc);
  mpa.f = moon_latitude_argument(spa.jc);

  moon_periodic_term_summation(mpa.d, mpa.m, mpa.m_prime, mpa.f, spa.jc,
                               ML_TERMS, &mpa.l, &mpa.r);
  moon_periodic_term_summation(mpa.d, mpa.m, mpa.m_prime, mpa.f, spa.jc,
                               MB_TERMS, &mpa.b, 0);

  moon_longitude_and_latitude(spa.jc, mpa.l_prime, mpa.f, mpa.m_prime,
                              mpa.l, mpa.b, &mpa.lamda_prime, &mpa.beta);

  mpa.cap_delta = moon_earth_distance(mpa.r);
  mpa.pi = moon_equatorial_horiz_parallax(mpa.cap_delta);

  mpa.lamda = apparent_moon_longitude(mpa.lamda_prime, spa.del_psi);

  mpa.alpha = geocentric_right_ascension(mpa.lamda, spa.epsilon, mpa.beta);
  mpa.delta = geocentric_declination(mpa.beta, spa.epsilon, mpa.lamda);

  mpa.h = observer_hour_angle(spa.nu, spa.longitude, mpa.alpha);

  right_ascension_parallax_and_topocentric_dec(
      spa.latitude, spa.elevation, mpa.pi, mpa.h, mpa.delta,
      &(mpa.del_alpha), &(mpa.delta_prime));
  mpa.h_prime = topocentric_local_hour_angle(mpa.h, mpa.del_alpha);

  mpa.e0 = topocentric_elevation_angle(spa.latitude, mpa.delta_prime,
                                        mpa.h_prime);
  mpa.del_e = atmospheric_refraction_correction(
      spa.pressure, spa.temperature, spa.atmos_refract, mpa.e0);
  mpa.e = topocentric_elevation_angle_corrected(mpa.e0, mpa.del_e);

  output->zenith = topocentric_zenith_angle(mpa.e);
  mpa.azimuth_astro = topocentric_azimuth_angle_astro(
      mpa.h_prime, spa.latitude, mpa.delta_prime);
  output->azimuth = topocentric_azimuth_angle(mpa.azimuth_astro);
}

}  // namespace mpa
