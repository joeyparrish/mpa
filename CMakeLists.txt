cmake_minimum_required(VERSION 3.22)
project(libmpa)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Wall -Wunused-function -Werror")

add_library(sampa
  sampa/spa.cc
  sampa/spa.h
  sampa/sampa.cc
  sampa/sampa.h
  sampa/bird.cc
  sampa/bird.h
)

add_library(mpa
  mpa/util.cc
  mpa/util.h
  mpa/spa.cc
  mpa/spa.h
  mpa/mpa.cc
  mpa/mpa.h
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
