#!/usr/bin/python2

"""functions to generate maximum length sequences and to
calculate the permutations needed for cross-correlating
them with the the Fast Hadamard Transform. It has been
written to assist in preparing data for use in embedded 
applications. This module can generate C-code (only the 
arrays, not the code) for the permutation tables.

(c) 2015 Bernd Kreuss <prof7bit@gmail.com>
"""

import math
import unittest

def debug_print(vect, n=2):
    """print a list or a list of lists"""
    for x in vect:
        if type(x) is list:
            debug_print(x, n)
        else:
            print "%*g" % (n,x),
    print


def poly_degree(poly):
    """Return the degree of the polynomial"""
    l = -1
    while poly:
        poly >>= 1
        l += 1
    return l


def rotate(vect, x):
    """return a rotated list containing the elements of vect 
    rotated by x elements to the left. Positive values of x
    rotate to the left, negative values rotate to the right
    """
    P = len(vect)
    out = [0] * P
    for i in range(P):
        out[i] = vect[(i + x) % P]
    return out


def mls_length_from_poly(poly):
    """Return the period length of the generated mls.
    This function does NOT check whether it really is
    a primitive polynomial, it will just calculate
    how long such an MLS would be."""
    return (2 ** poly_degree(poly)) - 1


def generate_mls(poly):
    """Implement the LFSR and return a list containing 
    one period of the generated MLS with elements consisting 
    of 1 or 0. This function will not check whether the 
    polynomial is really primitive, it will just blindly 
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


def invert_permutation(perm):
    """Return the permutation vector for the inversion of
    perm, that is the permutation that when applied after
    perm would would undo perm and bring all elements back 
    to where they have been.
    """
    inv = [0] * len(perm)
    for i in range(len(perm)):
        inv[perm[i]] = i
    return inv


def permute(vect, perm):
    """apply the permutation perm to vect and return a 
    new list with the permuted elements of vect
    """
    out = [0] * len(vect)
    for i in range(len(vect)):
        out[i] = vect[perm[i]]
    return out


def generate_permutations(poly):
    """Generate the permutations for the Fast Hadamard.
    This function will return two permutation vectors of length 
    Period + 1 so that their length is a power of 2. This is 
    because it is meant to be applied after prepending the 0
    and before removing it again to be compatible with the 
    implementation of the butterfly in this module. The 
    first element of the vectors will be 0 (meaning it will 
    not permute the first element of the vector).
    
    This function implements the algorithm described in
    J. Borish & J. B. Angell, "An Efficient Algorithm for 
    Measuring the Impulse Response Using Pseudrandom Noise," 
    J. Audio Eng. Soc., vol. 31, pp. 478-489 (1983) 
    """
    mls = generate_mls(poly)
    P = len(mls)
    N = poly_degree(poly)

    # As illustrated in Borish's paper, fig. 5 we
    # construct the PxP matrix of shifted mls and
    # determine the 'tags' from the first N rows.
    tags = [0]
    for col in range(P):
        tag = 0
        bitmask = 1
        for row in range(N):
            x = mls[(col - row) % P]
            
            # interpret each column as a binary number, 
            # call them 'tags' like in the Borish paper
            if x > 0:
                tag |= bitmask
            bitmask <<= 1
        
        # and make a list of them
        tags.append(tag)

    # these tags are the inverted permutation vector for 
    # the input samples, so all we need to do is invert it
    # and we are done. Note that we have prepended a 0 in
    # front of it and all other elements are > 0. This is
    # very convenient because later we want to apply it to
    # sampled data that has a 0 prepended too because of
    # the 2**N requirement of the butterfly algorithm.
    perm_in = invert_permutation(tags)
        

    # now we also need the output permutation
    #
    # This is illustrated in Borish's paper in fig. 7.
    # Again we start with the matrix M like above but this time
    # we rearrange the columns with the above column permutation
    # and from that matrix we pick only the columns
    # 1, 2, 4, ..., 2**(N-1) to supply the bits for the tags.
    # The following nested loops will do this in one go:
    perm_out = [0]
    for row in range(P):
        tag = 0
        bitmask = 1
        for n in range(N):
            perm_col = 2**n                     # the neeed column in the permuted matrix
            real_col = perm_in[perm_col] - 1    # corresponding column in original matrix
            x = mls[(real_col - row) % P]
            
            # interpret each row as a binary number,
            # (calling it 'tag' like in the Borish paper)
            if x > 0:
                tag |= bitmask
            bitmask <<= 1
        
        # and here these tags ARE the permutation vector!
        perm_out.append(tag)

    return (perm_in, perm_out)


def butterfly(x, shave_bits=0):
    """Apply the butterfly algorithm to the input vector x.
    This function will operate directly on the elements of
    the input vector x. The samples must first be permuted
    with the input permutation, then the butterfly can be
    applied and afterwards they need to be permuted with 
    the output permutation to bring the result into the 
    correct chronologiocal order again.
    
    When shave_bits > 0 then it will shave off the least 
    significant bit with a division by two during each 
    of the last n rounds, where n=shave_bits. This trick
    is not needed when running this algorithm in python on
    a pc but it is of great usefulness when implementing this 
    algorithm on a small microcontroller where only short 
    integers can be used to hold the huge sample arrays.
    """
    assert(len(x) & (len(x)-1) == 0), "length must be power of 2"
    samples = len(x)
    span = samples
    while span > 1:
        next = span
        span = span >> 1
        shave_last_bit = (span < (1 << shave_bits))
        if shave_last_bit:
            for j in range(0, span):
                i = j
                while i < samples:
                    temp            = x[i + span]
                    x[i + span]     = (x[i] - temp) / 2
                    x[i]            = (x[i] + temp) / 2
                    i += next
        else:
            for j in range(0, span):
                i = j
                while i < samples:
                    temp            = x[i + span]
                    x[i + span]     = x[i] - temp
                    x[i]            = x[i] + temp
                    i += next
                
                
def find_primitive_poly(poly, degree):
    """Find primitive polynomial.
    Original Pascal code written by Hagen Reddmann, 
    published in https://www.mikrocontroller.net/topic/279499#2950484
    ported to Python more or lesss without modifications.
    
    you should not call this directly but instead
    call the function find_all_primitive_polys()
    which will make use of this code.
    
    Please don't ask me any questions about this 
    code, ask Hagen instead because he wrote it.'
    """
    assert(degree > 0)
    assert(degree < 32)
    
    def even_parity(poly):
        """returns TRUE if count of bits set to 1 in Poly is even"""
        p = True
        while poly:
            if poly & 1:
                p = not p
            poly >>= 1
        return p
    
    def poly_init(M, poly, degree):
        l = 0x80000000
        M[degree - 1] = l
        for i in range(degree - 1):
            l >>= 1
            M[i] = l
        for i in range(degree - 1):
            if poly & 1:
                M[i] |= 0x80000000
            poly >>= 1
            
    def poly_copy(D, S):
        for i in range(len(S)):
            D[i] = S[i]
    
    def poly_mul(R, M, degree):
        T = [0] * 32 # TPolyMatrix
        for i in range(degree):
            n = M[i]
            d = 0
            for j in range(degree):
                if (n & 0x80000000) != 0:
                    d = d ^ R[j]
                n <<= 1
            T[i] = d
        poly_copy(R, T)
        
    def poly_pow_mod(R, M, n, degree):
        poly_copy(R, M)
        l = 0x80000000
        while (n & l) == 0:
            l >>= 1
        while l > 1:
            l >>= 1
            poly_mul(R, R, degree)
            if (l & n) != 0:
                poly_mul(R, M, degree)
    
    def poly_is_primitive(poly, degree, factors, factor_count):
        P = [0] * 32   # TPolyMatrix in original pascal code
        M = [0] * 32
        poly_init(M, poly, degree)
        poly_copy(P, M)
        state = P[0]
        for i in range(1, degree+1):
            poly_mul(P, P, degree)
            if P[0] == state:
                if i == degree:
                    for j in range(factor_count):
                        poly_pow_mod(P, M, factors[j], degree)
                        if P[0] == 0x80000000:
                            return False
                    return True
                else:
                    return False
        return False
    
    def factor_order(factors, degree):
        """find factors of 2^Degree-1 = possible Order of Polynom
        can be surrely more optimized, but the runtime here is not important yet
        as example:
            For all orders from 2^0 -1 upto 2^32-1 exists only 33 possible primefactors.
            Instead to looping trough all odd numbers as factors we could reduce the
            loopcount if we use a lookuptable of all the 33 possible primefactors.
        """
        result = 0
        last_factor = 0
        prime = 3
        order = 0xffffffff >> (32 - degree)
        rest = order
        bound = int(round(math.sqrt(rest)))
        while (rest != 1) and (prime < bound):
            if rest % prime == 0:
                rest = rest // prime
                factor = order // prime
                if factor != last_factor:
                    last_factor = factor
                    factors[result] = factor
                    result += 1
                bound = int(round(math.sqrt(rest)))
            else:
                prime += 2
        if result > 0: 
            # only if 2^Degree-1 itself isn't prime
            factors[result] = order // rest
            result += 1
        return result
    
    factors = [0] * 6
    mask = 0xffffffff >> (32 - degree)
    poly = (poly & mask) | (1 << (degree - 1))
    factor_count = factor_order(factors, degree)
    while poly <= mask:
        if poly_is_primitive(poly, degree, factors, factor_count):
            return poly
        else:
            while True:
                poly += 1
                if (poly > mask) or even_parity(poly):
                    break


def find_all_primitive_polys(degree):
    """Find all primitive polynomials of certain degree.
    Each polynomial is represented by a binary number 
    where the nth bit represents the term x^n. for 
    example consider the polynomial:
    
        x^4 + x + 1
        
    this would be represented by the binary number
    
        10011
        
    all functions in this module are using the same
    binary representation for the polynomials.
    """
    poly = 0
    results = []
    while True:
        poly = find_primitive_poly(poly, degree)
        if poly > 0:
            poly_with_1 = (poly << 1) | 1
            results.append(poly_with_1)
            poly += 1
        else:
            return results


def generate_c_array(name, arr):
    count = len(arr)
    
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
        line.append(format % arr[x])
        if len(line) == cols:
            lines.append(list(line))
            line = []
    if len(line):
        lines.append(list(line))
    
    lines = [", ".join(line) for line in lines]
    datablock = "    " + ",\n    ".join(lines)
    return ("const %s %s[%g] = {\n" % (type, name, count)) + datablock + "\n};\n\n"


def generate_c(poly):
    """generate C code for an arrays that holds the permutation
    vectors. It does not generate any other code, you will have 
    to come up with a port of the butterfly() function suitable 
    and optimized for your specific application yourself."""

    pi, po = generate_permutations(poly)
    code =  generate_c_array("input_permutation", pi)
    code += generate_c_array("output_permutation", po)
    return code


# ***********************************
# * unit tests (and usage examples) *
# ***********************************

class TestCase(unittest.TestCase):
    def test_invert_permutation(self):
        p = [0,19,20,6,21,24,7,30,17,22,15,25,27,8,11,31,18,5,23,29,16,14,26,10,4,28,13,9,3,12,2,1]
        pi = invert_permutation(p)
        pn = permute(p, pi)
        self.assertEqual(pi, [0,31,30,28,24,17,3,6,13,27,23,14,29,26,21,10,20,8,16,1,2,4,9,18,5,11,22,12,25,19,7,15])
        self.assertEqual(pn, [i for i in range(32)])
    
    def test_rotate(self):
        x = [1,2,3,4,5,6,7,8]
        self.assertEqual(rotate(x,-2), [7,8,1,2,3,4,5,6])
        self.assertEqual(rotate(x, 2), [3,4,5,6,7,8,1,2])
    
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
    
    def test_generate_permutations(self):
        POLY3 = 0x0b
        pi, po = generate_permutations(POLY3)
        self.assertEqual(pi, [0, 6, 7, 2, 5, 1, 4, 3])
        self.assertEqual(po, [0, 1, 2, 4, 5, 7, 3, 6])
        
    def test_butterfly(self):
        poly = 0x25
        #poly = 0x0b
        mls = generate_mls(poly)
        P = mls_length_from_poly(poly)
        
        # generate some sample data by rotating the mls
        # P times, prepending the 0 and also generate the
        # expected correlation results for each of them
        smplist = []
        reslist = []
        P = len(mls)
        smp0 = [2 * x - 1 for x in mls]
        res0 = [-P] + [1] * (P - 1)
        for i in range(P):
            smp = [0] + rotate(smp0, -i)
            res = [1] + rotate(res0, -i)
            smplist.append(smp)
            reslist.append(res)
            
        #debug_print(smplist, 3)
        #debug_print(reslist, 3)
            
        p1, p2 = generate_permutations(poly)
                
        for i in range(len(smplist)):
            samples = permute(smplist[i], p1)
            butterfly(samples)
            result = permute(samples, p2)
            self.assertEqual(result, reslist[i])
        
    def test_butterfly_large_with_shave(self):
        poly = 0x1107 # this will be quite large
        mls = generate_mls(poly)
        
        # simulate samples that were made with a 12bit ADC
        # and use the full range. This would cause the maximum
        # possible correlation result to overflow an int16_t
        # by orders of magnitude if we would not already begin
        # shaving off bits  during the butterfly.
        samples = [0] + [4095 * x - 2048 for x in mls]
        
        pi, po = generate_permutations(poly)
        samples = permute(samples, pi)
        
        # we must shave off at least 8 bits to make it not 
        # overflow a 16 bit signed integer.
        butterfly(samples, 8)
        
        #debug_print(samples)
        self.assertEqual(min(samples), -32752)
        
    def test_hagens_code(self):
        degree = 8
        results = find_all_primitive_polys(degree)
        self.assertEqual(results, [285, 299, 301, 333, 351, 355, 357, 361, 369, 391, 397, 425, 451, 463, 487, 501])
        
if __name__ == '__main__':    
    unittest.main(exit=False)
