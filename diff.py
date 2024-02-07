#!/usr/bin/python

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

with open(file1, 'r') as f1, open(file2, 'r') as f2:
    for line1, line2 in zip(f1, f2):
        if line1 != line2:
            # split the lines on tab, print the first, fourth, sixth and last fields side by side for each file
            print('\t'.join([line1.split('\t')[0], line1.split('\t')[3], line1.split('\t')[5], line1.split('\t')[-1].strip()]), '\t', line2.split('\t')[-1].strip())