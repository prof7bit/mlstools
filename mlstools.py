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
import math
import unittest

def poly_degree(poly):
    """Return the degree of the polynomial"""
    l = -1
    while poly:
        poly >>= 1
        l += 1
    return l


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
    inv = [0] * len(perm)
    for i in range(len(perm)):
        inv[perm[i]] = i
    return inv

def permute(vect, perm):
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
    implementation of the in-place FHT in this module. The 
    first element of the vectors will be 0 (meaning it will 
    not permute the first element of the vector).
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
    perm_in = [0]
    for num,idx in list:
        perm_in.append(idx + 1)
    
    # Now that we have the input permutation
    # we also need the output permutation. We get
    # it by constructing L and reordering its rows in a similar way
    #
    # now lets create the matrix L
    # there are P rows with N elements
    # 
    # we construct it by starting again with M
    # (for convenience this time we would express M as a list of rows)
    # such a matrix would look like this:
    #
    # M = [[mls[(j + i) % P] for i in range(P)] for j in range(P)]
    #
    # now reorder the columns (index i) according to perm_in,
    # we could do the reordering on the fly, so we would replace 
    # the above list comprehension with this one;
    # (note that the indices in perm_in[] begin at 1 and not at 0)
    #
    # M1 = [[mls[(j + perm_in[i + 1] - 1) % P] for i in range(P)] for j in range(P)]
    # 
    # but we are still not done yet, we only need the 1st, 2nd, 4th, ...,  x**(N-1)th 
    # column of M1 to form the matrix L so that the upper N*N elements of L would
    # form an identity matrix. We can put all of the above in one monster
    # list comprehension and our matrix L now looks like this:
    # (a list of P rows with L elements each):
    L = [[mls[(j + perm_in[(2**i)] - 1) % P] for i in range(N)] for j in range(P)]
 
    # for j in range(P): print L[j];
    
    # now we determine the correct row order like we did it above by 
    # interpreting each row as a binary number, left column (j=0) shall 
    # be the least significant bit. These values are the permutation vector
    perm_out = [0]
    for i in range(P):
        num = 0
        for j in range(N):
            num += L[i][j] * (1 << j)
        perm_out.insert(1, num)

    return (perm_in, perm_out)

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


def generate_c(poly):
    """generate C code for an array that holds the permutation
    vector. It does not generate any other code, you will have 
    to come up with a port of the inplace_permuted_butterfly() 
    function suitable and optimized for your specific application
    yourself."""
    pi, po = generate_permutations(poly)
    count = len(pi)
    
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
        line.append(format % pi[x])
        if len(line) == cols:
            lines.append(list(line))
            line = []
    if len(line):
        lines.append(list(line))
    
    lines = [", ".join(line) for line in lines]
    datablock = "    " + ",\n    ".join(lines)
    return ("const %s input_permutation[%g] = {\n" % (type, count)) + datablock + "\n};\n"



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
    
    def test_generate_permutations(self):
        POLY3 = 0x0b
        pi, po = generate_permutations(POLY3)
        self.assertEqual(pi, [0, 4, 5, 7, 3, 6, 2, 1])
        self.assertEqual(po, [0, 5, 7, 3, 6, 1, 2, 4])
        
    def test_inplace_permuted_butterfly(self):
        poly = 0x25
        #poly = 0x0b
        mls = generate_mls(poly)
        P = mls_length_from_poly(poly)
        
        # generate some sample data by rotating the mls
        # P times, prepending the 0 and also generate the
        # expected correlation results for each of them
        siglist = []
        reslist = []
        P = len(mls)
        for i in range(P):
            sig = [0]
            for j in range(P):
                sig.append(2 * mls[(i+j) % P] - 1)

            res = [1]
            res += [1]*(P-i-1)
            res += [-P] 
            res += [1]*i 
                
            siglist.append(sig)
            reslist.append(res)
            
            # example siglist
            # [0, 1, 1, 1, -1, -1, 1, -1]
            # [0, 1, 1, -1, -1, 1, -1, 1]
            # [0, 1, -1, -1, 1, -1, 1, 1]
            # [0, -1, -1, 1, -1, 1, 1, 1]
            # [0, -1, 1, -1, 1, 1, 1, -1]
            # [0, 1, -1, 1, 1, 1, -1, -1]
            # [0, -1, 1, 1, 1, -1, -1, 1]
            #
            # example reslist:
            # [1, 1, 1, 1, 1, 1, 1, -7]
            # [1, 1, 1, 1, 1, 1, -7, 1]
            # [1, 1, 1, 1, 1, -7, 1, 1]
            # [1, 1, 1, 1, -7, 1, 1, 1]
            # [1, 1, 1, -7, 1, 1, 1, 1]
            # [1, 1, -7, 1, 1, 1, 1, 1]
            # [1, -7, 1, 1, 1, 1, 1, 1]
        
        p1, p2 = generate_permutations(poly)
        
        # combined permutation: first p1 then p2
        p3 = permute(p1, p2)
        
        for i in range(len(siglist)):
            sig = siglist[i]
            inplace_permuted_butterfly(sig, p1)
            res = permute(sig, p3)
            self.assertEqual(res, reslist[i])

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
        pi, po = generate_permutations(poly)
        self.assertEqual(pi, [0, 19, 20, 6, 21, 24, 7, 30, 17, 22, 15, 25, 27, 8, 11, 31, 18, 5, 23, 29, 16, 14, 26, 10, 4, 28, 13, 9, 3, 12, 2, 1])
        
        # all information about the mls is encoded in the permutation,
        # the butterfly will use it to know how to fold the samples
        inplace_permuted_butterfly(samples, pi)
        self.assertEqual(samples, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -31, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        
    def test_inplace_permuted_butterfly_with_own_generated_perms_large(self):
        poly = 0x1107 # this will be quite large
        mls = generate_mls(poly)
        samples = [0] + [2 * x - 1 for x in mls]
        pi, po = generate_permutations(poly)
        inplace_permuted_butterfly(samples, pi)
        
        # all of them will be 1 except one single element:
        self.assertEqual(samples, [1]*2288 + [-4095] + [1]*1807)
        
    def test_hagens_code(self):
        degree = 8
        results = find_all_primitive_polys(degree)
        self.assertEqual(results, [285, 299, 301, 333, 351, 355, 357, 361, 369, 391, 397, 425, 451, 463, 487, 501])
        

if __name__ == '__main__':    
    if len(sys.argv) == 1:
        unittest.main(exit=False)
        print "\nusage: ./mlstools.py <function> [<argument> [<argument>]]"
        print "example: ./mlstools.py generate_c 0x211"
    else:
        print eval(sys.argv[1] + "(" + ",".join(sys.argv[2:]) + ")")
