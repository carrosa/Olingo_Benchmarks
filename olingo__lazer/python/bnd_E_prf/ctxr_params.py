from math import sqrt, exp
from sage.all import var, solve

vname = "param"

secpar=192

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


# sec=128

its = 1

l = secdim_hat*(1 + its)
n = 1
k = l + n + 1

deg = 256  # ring Rp degree d
mod = q_hat  # ring Rp modulus p
dim = (n + l, k + l)

# sigma_tdec = 1208
tail_bound_2 = 1.14
# B = sigma_tdec*tail_bound_2*sqrt(deg)
B = 0


for i in range(its):
    B = tail_bound_2 * sigma_tdec * sqrt(deg)
    beta = var('x')
    solution = solve(sqrt(deg) * beta**2 + sqrt(deg)*beta - B == 0, beta, solution_dict=True)
    beta_2 = solution[0][beta] if solution[0][beta] > 0 else solution[1][beta]
    beta_2 = int(beta_2)
    B = beta_2 * sqrt(deg)


# Partition of s
wpart = [list(range(k)), list(range(k, k + l))]


# l2-norm bounds
wl2 = [0, 5 * B]


# binary coefficients (0 = no restriction)
wbin = [1, 0]

# rejection sampling flags (0 = not used)
wrej = [0, 0]
