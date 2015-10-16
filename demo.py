#!/usr/bin/python2

import mlstools

POLY = 0x11d   # 8
#POLY = 0x409   # 10
#POLY = 0x1053  # 12

mls = mlstools.generate_mls(POLY)
p1, p2 = mlstools.generate_permutations(POLY)

print "Polynomial: %s\n" % bin(POLY)
print "Degree: %d" % mlstools.poly_degree(POLY)
print "Period: %d" % len(mls)
print "Sequence: %s\n" % str(mls)
print "C-arrays:\n"
print mlstools.generate_c(POLY)
