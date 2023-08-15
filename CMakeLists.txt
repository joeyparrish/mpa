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
  mpa/bird.cc
  mpa/bird.h
)
