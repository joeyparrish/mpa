///////////////////////////////////////////////
//         HEADER FILE for SAMPA.C           //
//                                           //
// Solar and Moon Position Algorithm (SAMPA) //
//                   for                     //
//        Solar Radiation Application        //
//                                           //
//              August 1, 2012               //
//                                           //
//   Filename: SAMPA.H                       //
//                                           //
//   Afshin Michael Andreas                  //
//   Afshin.Andreas@NREL.gov (303)384-6383   //
//                                           //
//   Solar Resource and Forecasting Group    //
//   Solar Radiation Research Laboratory     //
//   National Renewable Energy Laboratory    //
//   15013 Denver W Pkwy, Golden, CO 80401   //
//                                           //
//  This code is based on the NREL           //
//  technical report "Solar Eclipse          //
//  Monitoring for Solar Energy Applications //
//  using the Solar and Moon Position 		 //
//  Algorithms" by Ibrahim Reda              //
///////////////////////////////////////////////

#pragma once

#include "spa.h"

namespace mpa {

typedef struct {
  //-----------------Intermediate MPA OUTPUT VALUES--------------------

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
  double alpha_prime;  // topocentric moon right ascension [degrees]
  double h_prime;      // topocentric local hour angle [degrees]

  double e0;     // topocentric elevation angle (uncorrected) [degrees]
  double del_e;  // atmospheric refraction correction [degrees]
  double e;      // topocentric elevation angle (corrected) [degrees]
  double azimuth_astro;  // topocentric azimuth angle (westward from south) [for
                         // astronomers]

  //---------------------Final MPA OUTPUT VALUES------------------------
  double zenith;         // topocentric zenith angle [degrees]
  double azimuth;        // topocentric azimuth angle (eastward from north) [for
                         // navigators and solar radiation]
} mpa_data;  // Moon Position Algorithm (MPA) structure

typedef struct {
  spa_data spa;  // Enter required INPUT VALUES into SPA structure (see SPA.H)
                 // spa.function will be forced to SPA_ZA, therefore slope &
                 // azm_rotation not required)

  mpa_data mpa;  // Moon Position Algorithm structure (defined above)
} sampa_data;  // Solar and Moon Position Algorithm (SAMPA) structure

// Calculate SAMPA output values (in structure) based on input values passed in
// structure
int sampa_calculate(sampa_data *sampa);

}  // namespace mpa
