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
} spa_data;

// Calculate SPA output values (in structure) based on input values passed in
// structure
void spa_calculate(const mpa_input& input, spa_data *spa);

}  // namespace mpa
