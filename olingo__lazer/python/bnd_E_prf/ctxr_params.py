from math import sqrt, exp
from sage.all import var, solve

vname = "param"

# sec=128
q = 2101003516513793
qhat = 633544653872304603170707504554316801

# sec=192
# q = 2435264670919169
# qhat = 1665928971017885857758033605809930753

# sec=256
#q = 3040015283376641
#qhat = 4015297092641499797387712543978426881

# sec=128 
l = 22
k = 25
n = 1
# sec=192
# l = 30
# k = 33
# n = 1
# sec=256
#l = 38
#k = 41
#n = 1

deg = 256  # ring Rp degree d
mod = qhat  # ring Rp modulus p
dim = (n + l, k + l)

sigma_tdec = 1208

tail_bound_2 = 1.14

# Partition of s
wpart = [list(range(k)), list(range(k, k+l))]


# l2-norm bounds
wl2 = [0, 6 * tail_bound_2 * sigma_tdec * sqrt(deg)]


# binary coefficients (0 = no restriction)
wbin = [1, 0]

# rejection sampling flags (0 = not used)
wrej = [0, 0]
