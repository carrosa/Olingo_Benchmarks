from math import sqrt, exp
from sage.all import var, solve

vname = "param"

secpar = 256

if secpar == 128:
    ringdim = 256
    q = 81623040078337
    q_hat = 23473657271235664258052971185590273
    moddim = 8
    secdim = 8
    sigma_y = 2^10
    sigma_w = 2^29
    kappa_y = 29
    kappa_w = 32
    nu = 23
    moddim_hat = secdim_hat = m = 21
    sigma_tdec = 123298420299561636
elif secpar == 192:
    ringdim = 256
    q = 174243429534721
    q_hat = 115233010602376640884019849787898369
    moddim = 11
    secdim = 11
    sigma_y = 2^10
    sigma_w = 2^30
    kappa_y = 29
    kappa_w = 33
    nu = 39
    moddim_hat = secdim_hat = m = 29
    sigma_tdec = 283537181095073949
elif secpar == 256:
    ringdim = 256
    q = 254716596705793
    q_hat = 327595548451764648065623801250983937
    moddim = 14
    secdim = 14
    sigma_y = 2^10
    sigma_w = 2^30
    kappa_y = 29
    kappa_w = 33
    nu = 60
    moddim_hat = secdim_hat = m = 37
    sigma_tdec = 551404451102749325


cmt_l = secdim_hat
cmt_n = 1
cmt_k = cmt_l + cmt_n + 1


deg = 256  # ring Rp degree d
mod = q_hat  # ring Rp modulus p
dim = (cmt_l + cmt_n, cmt_k + cmt_l)  # dimensions of A in Rp^(m,n)

tail_bound_2 = 1.14

# Partition of s
wpart = [list(range(cmt_k)), list(range(cmt_k, cmt_k + cmt_l))]


# l2-norm bounds
wl2 = [
    0,
    0,
    #    6*tail_bound_2 * sigma_tdec * sqrt(deg)
]


# binary coefficients (0 = no restriction)
wbin = [1, 1]

# rejection sampling flags (0 = not used)
wrej = [0 for _ in range(m + 1)]
