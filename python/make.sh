#!/bin/bash

# Script to compile cython

rm *.so *.c
python setup.py build_ext --inplace
