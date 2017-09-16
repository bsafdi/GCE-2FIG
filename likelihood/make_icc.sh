#!/bin/bash

# Script to compile cython

rm -r build/
rm *.so *.c
python setup_icc.py build_ext --inplace
