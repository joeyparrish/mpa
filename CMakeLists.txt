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
