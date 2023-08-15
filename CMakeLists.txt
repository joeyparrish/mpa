cmake_minimum_required(VERSION 3.22)
project(libmpa)

set (CMAKE_CXX_STANDARD 11)

add_library(sampa
  sampa/spa.cc
  sampa/spa.h
  sampa/sampa.cc
  sampa/sampa.h
  sampa/bird.cc
  sampa/bird.h
)

add_library(mpa
  mpa/spa.cc
  mpa/spa.h
  mpa/sampa.cc
  mpa/sampa.h
)

# Defines the libraries "gtest" and "gtest_main"
add_subdirectory(googletest)

add_executable(regression_test
  regression_test.cc
)

target_link_libraries(regression_test
  sampa
  mpa
  gtest
  gtest_main
)
