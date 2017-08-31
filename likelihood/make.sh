#!/bin/bash

# Script to compile cython

rm -r build/
rm *.so *.c
python setup.py build_ext --inplace
