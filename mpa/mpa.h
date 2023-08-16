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

#pragma once

#include <ctime>

namespace mpa {

struct mpa_input {
  int year;       // 4-digit year,   valid range: -2000 to 6000
  int month;      // 2-digit month,  valid range: 1 to  12
  int day;        // 2-digit day,    valid range: 1 to  31
  int hour;       // UTC hour,       valid range: 0 to  24
  int minute;     // UTC minute,     valid range: 0 to  59
  double second;  // UTC second,     valid range: 0 to <60

  double longitude;  // Observer longitude (degrees east of prime meridian)
                     // valid range: -180  to  180 degrees
  double latitude;   // Observer latitude (degrees north of equator)
                     // valid range: -90   to   90 degrees
  double elevation;  // Observer elevation (meters)
                     // valid range: -6500000 or higher meters

  double pressure;       // Annual average local pressure (millibars)
                         // valid range: 0 to 5000 millibars
  double temperature;    // Annual average local temperature (degrees Celsius)
                         // valid range: -273 to 6000 degrees Celsius
  double atmos_refract;  // Atmospheric refraction at sunrise and sunset
                         // (Typically 0.5667)
                         // valid range: -5   to   5 degrees
};

struct mpa_output {
  double zenith;   // topocentric zenith angle (degrees above horizon)
  double azimuth;  // topocentric azimuth angle (degrees east of true north)
};

void set_mpa_time(const time_t time, mpa_input* input);
void set_mpa_time(const struct tm& tm, mpa_input* input);
void compute_mpa(const mpa_input& input, mpa_output* output);

}  // namespace mpa
