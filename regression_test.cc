/**
 * Moon Position Algorithm (MPA) Regression Tests
 *
 * Copyright 2023 Joey Parrish <joey.parrish@gmail.com>
 *
 * Ensures that this ported, stripped-down version of MPA is compatible with
 * the original SAMPA implementation from the DoE and NREL.
 */

#include <cinttypes>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "gtest/gtest.h"

#include "sampa/sampa.h"
#include "mpa/mpa.h"

namespace {

int64_t roll_int(int64_t min, int64_t max) {  // inclusive
  // NOTE: This is not uniformly distributed, and that's fine for the
  // randomized testing we're doing here.
  return min + (random() % (max - min + 1));
}

double roll_double(double min, double max) {  // inclusive
  // NOTE: This is not uniformly distributed, and that's fine for the
  // randomized testing we're doing here.
  double zero_to_one =
      static_cast<double>(random()) / static_cast<double>(RAND_MAX);
  return min + ((max - min) * zero_to_one);
}

void compare_algorithms(time_t now, double latitude, double longitude, double elevation) {
  sampa::sampa_data old_data;
  mpa::mpa_input new_input;
  mpa::mpa_output new_output;

  memset(&old_data, 0, sizeof(old_data));
  memset(&new_input, 0, sizeof(new_input));
  memset(&new_output, 0, sizeof(new_output));

  // Break out the unix timestamp into fields.
  struct tm utc;
  gmtime_r(&now, &utc);

  // Fill in the data to drive the old algorithm.
  old_data.spa.year = 1900 + utc.tm_year;
  old_data.spa.month = utc.tm_mon + 1;
  old_data.spa.day = utc.tm_mday;

  old_data.spa.hour = utc.tm_hour;
  old_data.spa.minute = utc.tm_min;
  old_data.spa.second = utc.tm_sec;

  old_data.spa.latitude = latitude;
  old_data.spa.longitude = longitude;
  old_data.spa.elevation = elevation;

  // Fill in the data to drive the new algorithm.
  set_mpa_time(utc, &new_input);
  new_input.latitude = latitude;
  new_input.longitude = longitude;
  new_input.elevation = elevation;

  // If the old one fails, the test data was invalid.
  ASSERT_EQ(0, sampa::sampa_calculate(&old_data));
  // The new one returns void.
  mpa::compute_mpa(new_input, &new_output);

  EXPECT_DOUBLE_EQ(old_data.mpa.azimuth, new_output.azimuth) << "Azimuth mismatch";
  EXPECT_DOUBLE_EQ(old_data.mpa.zenith, new_output.zenith) << "Zenith mismatch";
}

}  // namespace

TEST(PositionAlgorithm, RandomInputs) {
  for (int i = 0; i < 1000; ++i) {
    // 2023 - 2083.
    time_t now = roll_int(1672552800LL, 3597544799LL);

    // Anywhere on earth.
    double latitude = roll_double(-90.0, 90.0);
    double longitude = roll_double(-180.0, 180.0);

    // From the Dead Sea to Everest.
    double elevation = roll_double(-420.0, 8849.0);

    compare_algorithms(now, latitude, longitude, elevation);
  }
}
