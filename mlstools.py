#!/usr/bin/python2

"""functions to generate maximum length sequences and to
calculate the permutations needed for cross-correlating
them with the the Fast Hadamard Transform. It has been
written to assist in preparing data for use in embedded 
applications. This module can generate C-code (only the 
arrays, not the code) for the permutation tables. It also 
shows how to perform the permutation on the fly during the 
butterfly to reduce memory usage on the target platform.

(c) 2015 Bernd Kreuss <prof7bit@gmail.com>
"""

import sys
import unittest

def poly_degree(poly):
    """Return the degree of the polynominal"""
    l = -1
    while poly:
        poly >>= 1
        l += 1
    return l


def mls_length_from_poly(poly):
    """Return the period length of the generated mls.
    This function does NOT check whether it really is
    a primitive polynominal, it will just calculate
    how long such an MLS would be."""
    return (2 ** poly_degree(poly)) - 1


def generate_mls(poly):
    """Implement the LFSR and return a list containing 
    one period of the generated MLS with elements consisting 
    of 1 or 0. This function will not check whether the 
    polynominal is really primitive, it will just blindly 
    assume it is."""
    deg = poly_degree(poly)
    length = mls_length_from_poly(poly)
    register = [1] * deg
    mls = []

    # make a list of taps
    taps = []
    n = 1
    poly >>= 1
    while poly > 1:
        if poly & 1:
            taps.append(n)
        poly >>= 1
        n += 1
                
    # iterate over the expected output length while 
    # performing the feedback and the register rotation
    for i in range(length):
        feedback = register[0]
        for tap in taps:
            feedback ^= register[tap]
        register.append(feedback)
        mls.append(register.pop(0))

    return mls


def generate_input_permutation(poly):
    """Generate the input permutation for the Fast Hadamard.
    This function will return a permutation vector of length 
    Period + 1 so that its length is a power of 2. This is 
    because it is meant to be applied to the input vector AFTER 
    the input vector has already been prepended with 0 to be 
    compatible with the implementation of the in-place FHT in 
    this module. The first element of the vector will be 0 
    (meaning it will not permute the first element of the 
    input vector).
    """
    mls = generate_mls(poly)
    P = len(mls)
    N = poly_degree(poly)
    
    # M = LS where
    # L is a matrix of order P x N  
    # S is the N x P matrix formed by the first N rows of M
    #
    # now lets create the matrix S
    # there are P columns with N elements,
    S = [[mls[(i + j) % P] for i in range(N)] for j in range(P)]
    
    #print S
    
    # interpret each column as a binary number and
    # make a list of these numbers, along with their
    # associated column indices. Note that the first
    # row shall contain the most significant bit.
    list = []
    for i in range(P):
        num = 0
        for j in range(N):
            num += S[i][N - j - 1] * (1 << j)
        list.append((num, i))
    
    #print list
    
    # now we sort this list by the binary number and this 
    # will then also reorder the attached column index which 
    # is then the desired permutation for the input samples.
    # Because our inplace butterfly will need the permutation
    # and the samples vector to have a power of 2 elements we
    # prepend a 0 and make the permutation indices begin at 1
    list.sort()
    perm = [0]
    for num,idx in list:
        perm.append(idx + 1)
    
    return perm


def generate_output_permutation(poly):
    """not yet implemented"""
    pass


def inplace_permuted_butterfly(x, perm):
    """Apply a 'permuted' butterfly to the input vector x.
    This function will operate directly on the elements of
    the input vector x. It will NOT move the positions of its 
    elements, instead it will apply the permutation to the 
    butterfly algorithm itself on the fly when accessing the 
    list. If you  need this function to calculate an impulse 
    response you must perform  a combined permutation of both 
    (input and output) to the vector x AFTER this function has
    been applied. If you only need the highest amplitude you 
    don't need to apply any permutations at all."""
    assert(len(x) & (len(x)-1) == 0), "length must be power of 2"
    assert(len(x) == len(perm)), "x must have same length as perm"
    assert(perm[0] == 0), "first element of perm must be 0"
    assert(min(perm) == 0), "no negative elements allowed in perm"
    assert(max(perm) == len(perm) - 1), "greatest element must be len-1"
    
    samples = len(x)
    span = samples
    while span > 1:
        next = span
        span = span >> 1
        for j in range(0, span):
            i = j
            while i < samples:
                temp                = x[perm[i + span]]
                x[perm[i + span]]   = x[perm[i]] - temp
                x[perm[i]]          = x[perm[i]] + temp
                i += next    
                

def generate_c(poly):
    """generate C code for an array that holds the permutation
    vector. It does not generate any other code, you will have 
    to come up with a port of the inplace_permuted_butterfly() 
    function suitable and optimized for your specific application
    yourself."""
    perm = generate_input_permutation(poly)
    count = len(perm)
    
    if count > 256:
        type = "uint16_t"
        format = "%5i"
    else:
        type = "uint8_t"
        format = "%3i"
    
    cols = 16
    lines = []
    line = []
    for x in range(count):
        line.append(format % perm[x])
        if len(line) == cols:
            lines.append(list(line))
            line = []
    if len(line):
        lines.append(list(line))
    
    lines = [", ".join(line) for line in lines]
    datablock = "    " + ",\n    ".join(lines)
    return ("const %s permutation[%g] = {\n" % (type, count)) + datablock + "\n};\n"



# ***********************************
# * unit tests (and usage examples) *
# ***********************************

class TestCase(unittest.TestCase):
    def test_poly_degree(self):
        self.assertEqual(poly_degree(0x0b), 3)
        self.assertEqual(poly_degree(0x1053), 12)

    def test_mls_length_from_poly(self):
        self.assertEqual(mls_length_from_poly(0x0b), 7)
        self.assertEqual(mls_length_from_poly(0x1053), 4095)

    def test_generate_mls(self):
        POLY3 = 0x0B
        POLY4 = 0x19
        POLY5 = 0x25
        m3 = generate_mls(POLY3)
        m4 = generate_mls(POLY4)
        m5 = generate_mls(POLY5)
        self.assertEqual(m3, [1, 1, 1, 0, 0, 1, 0])
        self.assertEqual(m4, [1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0])
        self.assertEqual(m5, [1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0])

    def test_generate_input_permutation(self):
        POLY3 = 0x0b
        perm = generate_input_permutation(POLY3)
        self.assertEqual(perm, [0, 4, 5, 7, 3, 6, 2, 1])

    def test_inplace_permuted_butterfly(self):
        # test vectors manually created and found in other examples
        # x and perm prepended with 0 to make them 2^n elements
        x1    = [0,  1,  1,  1, -1,  1, -1, -1]
        x2    = [0, -1,  1,  1,  1, -1,  1, -1]
        x3    = [0, -1, -1,  1,  1,  1, -1,  1]
        x4    = [0,  1, -1, -1,  1,  1,  1, -1]
        x5    = [0, -1,  1, -1, -1,  1,  1,  1]
        x6    = [0,  1, -1,  1, -1, -1,  1,  1]
        x7    = [0,  1,  1, -1,  1, -1, -1,  1]
        perm = [0,  3,  2,  7,  1,  4,  6,  5]
        inplace_permuted_butterfly(x1, perm)
        inplace_permuted_butterfly(x2, perm)
        inplace_permuted_butterfly(x3, perm)
        inplace_permuted_butterfly(x4, perm)
        inplace_permuted_butterfly(x5, perm)
        inplace_permuted_butterfly(x6, perm)
        inplace_permuted_butterfly(x7, perm)
        self.assertEqual(x1, [1, 1, 1, 1, 1, -7, 1, 1])
        self.assertEqual(x2, [1, 1, 1, 1, 1, 1, 1, -7])
        self.assertEqual(x3, [1, 1, 1, -7, 1, 1, 1, 1])
        self.assertEqual(x4, [1, -7, 1, 1, 1, 1, 1, 1])
        self.assertEqual(x5, [1, 1, -7, 1, 1, 1, 1, 1])
        self.assertEqual(x6, [1, 1, 1, 1, -7, 1, 1, 1])
        self.assertEqual(x7, [1, 1, 1, 1, 1, 1, -7, 1])
        
    def test_inplace_permuted_butterfly_with_own_generated_perms(self):
        poly = 0x25
        mls = generate_mls(poly)
        self.assertEqual(mls, [1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0])
        
        # We need a power of 2 number of measurement samples for the
        # inplace-butterfly. We create them directly from the mls and
        # prepend 0 as the first element.
        samples = [0] + [2 * x - 1 for x in mls]
        self.assertEqual(samples, [0, 1, 1, 1, 1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1])
        
        # the generated permutation already has the 0 prepended for the correct size
        perm = generate_input_permutation(poly)
        self.assertEqual(perm, [0, 19, 20, 6, 21, 24, 7, 30, 17, 22, 15, 25, 27, 8, 11, 31, 18, 5, 23, 29, 16, 14, 26, 10, 4, 28, 13, 9, 3, 12, 2, 1])
        
        # all information about the mls is encoded in the permutation,
        # the butterfly will use it to know how to fold the samples
        inplace_permuted_butterfly(samples, perm)
        self.assertEqual(samples, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -31, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        
    def test_inplace_permuted_butterfly_with_own_generated_perms_large(self):
        poly = 0x1107 # this will be quite large
        mls = generate_mls(poly)
        samples = [0] + [2 * x - 1 for x in mls]
        perm = generate_input_permutation(poly)
        inplace_permuted_butterfly(samples, perm)
        
        # all of them will be 1 except one single element:
        self.assertEqual(samples, [1]*2288 + [-4095] + [1]*1807)
        

if __name__ == '__main__':    
    if len(sys.argv) == 1:
        unittest.main(exit=False)
        print "\nusage: ./mlstools.py <function> [<argument> [<argument>]]"
        print "example: ./mlstools.py generate_c 0x211"
    else:
        print eval(sys.argv[1] + "(" + ",".join(sys.argv[2:]) + ")")
