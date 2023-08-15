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

#include "mpa.h"

namespace mpa {

typedef struct {
  // Values reused by mpa
  double jc;  // Julian century
  double del_psi;      // nutation longitude [degrees]
  double epsilon;      // ecliptic true obliquity  [degrees]
  double nu;       // Greenwich sidereal time [degrees]

  // Local values only
  double jd;  // Julian day
  double jm;  // Julian millennium

  double l;  // earth heliocentric longitude [degrees]
  double b;  // earth heliocentric latitude [degrees]
  double r;  // earth radius vector [Astronomical Units, AU]

  double theta;  // geocentric longitude [degrees]
  double beta;   // geocentric latitude [degrees]

  double del_epsilon;  // nutation obliquity [degrees]
  double epsilon0;     // ecliptic mean obliquity [arc seconds]

  double del_tau;  // aberration correction [degrees]
  double lamda;    // apparent sun longitude [degrees]
  double nu0;      // Greenwich mean sidereal time [degrees]

  double alpha;  // geocentric sun right ascension [degrees]
  double delta;  // geocentric sun declination [degrees]

  double h;            // observer hour angle [degrees]
  double xi;           // sun equatorial horizontal parallax [degrees]
  double del_alpha;    // sun right ascension parallax [degrees]
  double delta_prime;  // topocentric sun declination [degrees]
  double h_prime;      // topocentric local hour angle [degrees]

  double e0;             // topocentric elevation angle (uncorrected) [degrees]
  double del_e;          // atmospheric refraction correction [degrees]
  double e;              // topocentric elevation angle (corrected) [degrees]
  double azimuth_astro;  // topocentric azimuth angle (westward from south) [for
                         // astronomers]
} spa_data;

// Calculate SPA output values (in structure) based on input values passed in
// structure
void spa_calculate(const mpa_input& input, spa_data *spa);

}  // namespace mpa
