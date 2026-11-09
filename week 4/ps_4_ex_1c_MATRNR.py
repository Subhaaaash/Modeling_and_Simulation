"""Problem sheet 4, exercise 1(c).

Find the largest x = 2**(-n) for which 1.0 + x == 1.0 in double precision,
compare it with NumPy's machine epsilon and derive the number of mantissa bits.
"""
import numpy as np
import matplotlib.pyplot as plt  # required header of the sheet (not used in 1c)

# Halve x until adding it to 1.0 no longer changes 1.0.  The powers of two are
# scanned from large to small, so the first hit is the largest such x.
n = 0
x = 1.0
while 1.0 + x != 1.0:
    n += 1
    x = 2.0 ** (-n)

eps = np.finfo(float).eps
mantissa_bits_stored = -int(np.log2(eps))  # eps = 2**(-(p-1)) = 2**(-N_m)
precision = mantissa_bits_stored + 1       # + implicit leading bit

print(f"largest x = 2^-n with 1.0 + x == 1.0: n = {n}, x = {x:.6e}")
print(f"1.0 + 2^-{n}     == 1.0 : {1.0 + 2.0 ** (-n) == 1.0}")
print(f"1.0 + 2^-{n - 1} == 1.0 : {1.0 + 2.0 ** (-(n - 1)) == 1.0}")
print(f"numpy.finfo(float).eps  = {eps:.6e}")
print(f"eps / x                 = {eps / x:.1f}")
print(f"stored mantissa bits  N_m = {mantissa_bits_stored} "
      f"(numpy.finfo(float).nmant = {np.finfo(float).nmant})")
print(f"precision p = N_m + 1     = {precision} (including the implicit bit)")
