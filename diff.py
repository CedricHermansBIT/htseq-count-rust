#!/usr/bin/python

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

with open(file1, 'r') as f1, open(file2, 'r') as f2:
    for line1, line2 in zip(f1, f2):
        spline1 = line1.split('\t')
        spline2 = line2.split('\t')
        if spline1[-1] != spline2[-1]:
            # split the lines on tab, print the first, fourth, sixth and last fields side by side for each file
            print('\t'.join([spline1[0], spline1[2], spline1[3], spline1[5], spline1[-1].strip()]), '\t', line2.split('\t')[-1].strip())