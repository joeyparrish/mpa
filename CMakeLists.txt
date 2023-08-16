# Moon Position Algorithm (MPA)
#
# Based on SAMPA (Solar And Moon Position Algorithm) and SPA (Solar Position
# Algorithm), downloaded from the US Department of Energy, stripped to the
# barest essentials for the moon, and ported to C++.
#
# SPA and SAMPA were written in 2003 and 2012, respectively, by Afshin Michael
# Andreas (Afshin.Andreas@NREL.gov), Measurement & Instrumentation Team, Solar
# Radiation Research Laboratory, National Renewable Energy Laboratory (NREL).
# They were based on the NREL technical report "Solar Position Algorithm for
# Solar Radiation Application" by Ibrahim Reda & Afshin Andreas, and the NREL
# technical report "Solar Eclipse Monitoring for Solar Energy Applications
# using the Solar and Moon Position Algorithms" by Ibrahim Reda.
#
# Ported by Joey Parrish <joey.parrish@gmail.com>.
#
# Original copyright notices:
#
# Copyright (C) 2008-2012 Alliance for Sustainable Energy, LLC, All Rights
# Reserved.
#
# The Solar Position Algorithm ("Software") is code in development prepared by
# employees of the Alliance for Sustainable Energy, LLC, (hereinafter the
# "Contractor"), under Contract No. DE-AC36-08GO28308 ("Contract") with the
# U.S. Department of Energy (the "DOE"). The United States Government has been
# granted for itself and others acting on its behalf a paid-up, non-
# exclusive, irrevocable, worldwide license in the Software to reproduce,
# prepare derivative works, and perform publicly and display publicly.
# Beginning five (5) years after the date permission to assert copyright is
# obtained from the DOE, and subject to any subsequent five (5) year renewals,
# the United States Government is granted for itself and others acting on its
# behalf a paid-up, non-exclusive, irrevocable, worldwide license in the
# Software to reproduce, prepare derivative works, distribute copies to the
# public, perform publicly and display publicly, and to permit others to do
# so. If the Contractor ceases to make this computer software available, it
# may be obtained from DOE's Office of Scientific and Technical Information's
# Energy Science and Technology Software Center (ESTSC) at P.O. Box 1020, Oak
# Ridge, TN 37831-1020.
#
# THIS SOFTWARE IS PROVIDED BY THE CONTRACTOR "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL THE CONTRACTOR OR THE U.S. GOVERNMENT BE LIABLE FOR ANY SPECIAL,
# INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER, INCLUDING BUT
# NOT LIMITED TO CLAIMS ASSOCIATED WITH THE LOSS OF DATA OR PROFITS, WHICH MAY
# RESULT FROM AN ACTION IN CONTRACT, NEGLIGENCE OR OTHER TORTIOUS CLAIM THAT
# ARISES OUT OF OR IN CONNECTION WITH THE ACCESS, USE OR PERFORMANCE OF THIS
# SOFTWARE.
#
# The Software is being provided for internal, noncommercial purposes only and
# shall not be re-distributed. Please contact Jennifer Ramsey
# (Jennifer.Ramsey@nrel.gov) in the NREL Commercialization and Technology
# Transfer Office for information concerning a commercial license to use the
# Software.
#
# As a condition of using the Software in an application, the developer of the
# application agrees to reference the use of the Software and make this Notice
# readily accessible to any end-user in a Help|About screen or equivalent
# manner.

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
