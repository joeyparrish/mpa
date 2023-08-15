#!/bin/bash

set -e

dos2unix -q "$@"
clang-format --style=Google -i "$@"
