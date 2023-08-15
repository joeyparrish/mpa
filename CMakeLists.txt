cmake_minimum_required(VERSION 3.22)
project(libmpa)

set (CMAKE_CXX_STANDARD 11)

add_library(sampa
  sampa/spa.c
  sampa/spa.h
  sampa/sampa.c
  sampa/sampa.h
  sampa/bird.c
  sampa/bird.h
)
