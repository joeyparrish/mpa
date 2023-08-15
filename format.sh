#!/bin/bash

set -e

dos2unix "$@"
clang-format --style=Google -i "$@"
