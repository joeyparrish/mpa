#pragma once

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

void compute_mpa(const mpa_input& input, mpa_output* output);

}  // namespace mpa
