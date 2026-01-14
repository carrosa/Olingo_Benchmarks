import os
#os.chdir("/Users/phmbp2025/Desktop/lwe-estimator")

from estimator import *


#security parameter

secpar = 256

print("Security level: ", secpar)


# threshold sizes

n = 1024
t = 1023

print("Threshold (n, t)=", "(", n, ",", t,")")

# ring dimension, module dimension, sec dimension, tail bound

ringdim = 256

moddim = 0
secdim = 0

if secpar == 128: 
	moddim = 7
	secdim = 10

if secpar == 192: 
	moddim = 10
	secdim = 13

if secpar == 256: 
	moddim = 12
	secdim = 16

# Setting cut-off point for truncated gaussian sampling

rej_rate = 0.01

tail_bound = 0.01
tail_bound_2 = 1.01

while (2 * exp((-tail_bound**2) / 2)) > rej_rate:
    tail_bound = tail_bound + 0.01

while ((tail_bound_2**ringdim)*exp((ringdim/2)*(1-tail_bound_2^2))) > rej_rate:
    tail_bound_2 = tail_bound_2 + 0.01


# Total number of signatures

Q = 2**60

# Bit-dropping parameters 

if secpar==128:
	kappa_y = 38
	kappa_w = 41
if secpar==192:
	kappa_y = 36
	kappa_w = 40
if secpar==256:
	kappa_y = 35
	kappa_w = 41

# Gaussian widths for s and r, set according to Raccoon analysis


if secpar==128:
	sigma_r = 2**38
	sigma_s = 2**20
if secpar==192:
	sigma_r = 2**38
	sigma_s = 2**20
if secpar==256:
	sigma_r = 2**38
	sigma_s = 2**20


# Defining challenge space 1-norm

nu = 0
while ((2**nu) * binomial(ringdim,nu) <2**secpar):
    nu = nu + 1


# Hint-MLWE bound

B_HMLWE = Q * nu * (1+ringdim*(secpar+1+2 * log(ringdim,2))/sqrt(Q))

assert(nu*2**kappa_y+2**(kappa_w+1) <= sigma_r * sqrt(t * (1+moddim/secdim)))

assert(sigma_s <= sigma_r/nu)


# Gaussian width for MLWE security

sigma = sqrt(1/(2/(n*sigma_s**2) + 2*B_HMLWE/(t*sigma_r**2)))

B = e**(1/4) * (nu * sigma_s * n  + t * sigma_r)*sqrt(ringdim * (moddim + secdim)) + (nu * (2**kappa_y) + 2**(kappa_w+1)) * sqrt(ringdim * moddim)

#print("B_2z =", log(B,2).n())

B_STMSIS = B + sqrt(nu) + (nu * (2**(kappa_y)) + 2**(kappa_w + 1)) * sqrt(ringdim * moddim)
#print("B_STMSIS = ", log(B_STMSIS,2).n(), "bits")

B_MSIS = B_STMSIS - nu #Following Raccoon parameter selection recommendation

print("B_MSIS bits =", log(B_MSIS,2).n())

B_inf =  sqrt(ringdim) * tail_bound * sigma_s * sqrt(n) + tail_bound * sigma_r * sqrt(t)

# Setting signature Modulus

if secpar == 128:
	logq = 56

if secpar == 192:
	logq = 56

if secpar == 256:
	logq = 56


q = round(2^logq)
while  (q % (2*ringdim) != 1) or (gcd(q-1, 2*ringdim) == 1):   #setting for NTT friendliness
	q = next_prime(q)

print("Signature parameters: ", "( N =", ringdim, ", q bits=", round(log(q,2)), ", #sigs = 2^",log(Q,2), ", chalspace 1-norm =", nu, ", k_y =", kappa_y, ", k_w =", kappa_w, ", mod dim =", moddim, ", sec dim =", secdim, ", sigma_y =", log(sigma_s,2), ", sigma_w =", log(sigma_r,2), ")")

# Security check

#LWE for partial signing key
y_lwe_param = LWE.Parameters(n=ringdim*moddim, q=q, Xs=ND.DiscreteGaussian(sigma_s), Xe=ND.DiscreteGaussian(sigma_s), tag='Partial Signature Key LWE')

#LWE for partial commitment
w_lwe_param = LWE.Parameters(n=ringdim*moddim, q=q, Xs=ND.DiscreteGaussian(sigma_r), Xe=ND.DiscreteGaussian(sigma_r), tag='Partial Commitment LWE')

#LWE for combined signature 
lwe_param = LWE.Parameters(n=ringdim*moddim, q=q, Xs=ND.DiscreteGaussian(sigma), Xe=ND.DiscreteGaussian(sigma), tag='Combined Signature LWE')
print()

#print(log(B_MSIS,2).n())
assert(B_MSIS < (q-1)/2)

#SIS for security
sis_param = SIS.Parameters(n=(secdim)*ringdim, q=q, length_bound=B_MSIS, m=(moddim+secdim+1)*ringdim, norm=2, tag='Signature SIS')

print("----------")
print("Start LWE estimation")
print("----------")
print()

print("LWE Sec. of partial signing key")
y_lwe = LWE.estimate.rough(y_lwe_param)
print()
print("LWE Sec. of partial commitment")
w_lwe = LWE.estimate.rough(w_lwe_param)
print()

print("LWE Sec. of combined signature")
r_lwe = LWE.estimate.rough(lwe_param)
print()

print("Combined LWE for Hint LWE:")
print(r_lwe)
print()

print("----------")
print("Start SIS estimation")
print("----------")
print()

r_sis = SIS.estimate.rough(sis_param)

print("SIS for SelfTargetSIS")
print(r_sis)
print()


# Set parameters for encryption scheme
print("----------")
print("Start HE Parameters")
print("----------")
print()
p = q

# params for granular adjustment of security level
expansionfactor = 4 
adjust_factor = 0

if secpar == 128:
	adjust_factor = 17
if secpar == 192:
	adjust_factor = 21
if secpar == 256:
	adjust_factor = 25

#secdimhat = secdim * expansionfactor
#moddimhat = moddim * expansionfactor

m = round(secdim * expansionfactor - adjust_factor)

secdimhat = m
moddimhat = m

sigma_ctx = sqrt(2/3)
sigma_B = sqrt(2/3)

B_2 = nu * n * moddimhat * tail_bound_2 * sigma_B * ringdim * tail_bound_2 * sigma_ctx + n * nu * p * tail_bound_2 * sigma_ctx*sqrt(ringdim) + p * nu * n * secdimhat * tail_bound_2 *sigma_B * ringdim * tail_bound_2 * sigma_ctx + t * tail_bound_2 * sigma_B * ringdim * tail_bound_2 * sigma_ctx * moddimhat + p * t * tail_bound_2 * sigma_ctx * sqrt(ringdim) + t * p * secdimhat * tail_bound_2 * sigma_B * ringdim *tail_bound_2 * sigma_ctx

sigma_tdec = sqrt(Q * B_2**2 / p**2)

B_dec = nu * n * tail_bound_2 * sigma_s * sqrt(ringdim * secdim) + t * tail_bound_2 * sigma_r * sqrt(ringdim * secdim) +B_2 + p * t * tail_bound_2 * sigma_tdec 
q_enc = ((B_dec) * 2).n()
while  (q_enc % (2*ringdim) != 1) or (gcd(q_enc-1, 2*ringdim) == 1):
	q_enc = next_prime(q_enc)
q_enc_actual = q_enc
log_q_enc = log(q_enc_actual,2)
#print("q_enc:",N(q_enc_actual))
#print("q_enc:",str(q_enc_actual))
#print("log q_enc:", N(log_q_enc))

print("HE parameters: ", "( N =", ringdim, ", qhat bits=", N(log_q_enc,4), ", secdimhat =", secdimhat, ", moddimhat = ", moddimhat, ", m = ", m, ", secdim = ", secdim, ", moddim = ", moddim, ")")


v_enc_lwe_param = LWE.Parameters(n=ceil(ringdim*m), q=q_enc_actual, Xs=ND.DiscreteGaussian(sigma_ctx), Xe=ND.DiscreteGaussian(sigma_ctx),m=2**60, tag='HE LWE v')
u_enc_lwe_param = LWE.Parameters(n=ceil(ringdim*secdimhat), q=q_enc_actual, Xs=ND.DiscreteGaussian(sigma_ctx), Xe=ND.DiscreteGaussian(sigma_ctx),m=2**60, tag='HE LWE u')
pk_enc_lwe_param = LWE.Parameters(n=ceil(ringdim*moddimhat), q=q_enc_actual, Xs=ND.DiscreteGaussian(sigma_B), Xe=ND.DiscreteGaussian(sigma_B), tag='HE LWE Public Key')

print()
print("----------")
print("Start HE LWE estimation")
print("----------") 
print()

print("LWE Sec. of cipher text element 'v'")
v_lwe = LWE.estimate.rough(v_enc_lwe_param)
print()

print("LWE Sec. of cipher text element 'u'")
u_lwe = LWE.estimate.rough(u_enc_lwe_param)
print()

print("LWE Sec. encryption pk")
pk_lwe = LWE.estimate.rough(pk_enc_lwe_param)
print()


print("----------")
print("BDLOP Commitment Params")
print("----------")
print()

print("Computing hiding hardness:")
sig_com_rand = 0.667
A1_height = 1
A2_height = m
A1_A2_width = A1_height+A2_height+1
com_hiding_param = LWE.Parameters(n=ceil(ringdim*(A1_height+A2_height)), q=q_enc_actual, Xs=ND.DiscreteGaussian(sig_com_rand), Xe=ND.DiscreteGaussian(sig_com_rand), tag='com hiding')
com_hiding = LWE.estimate.rough(com_hiding_param)
print(com_hiding)


print()
print("Computing binding hardness:")
B_com_rand = 16*sig_com_rand*sqrt(nu*ringdim)
print(f"{B_com_rand=}")
sis_param = SIS.Parameters(n=(A1_height)*ringdim, q=q_enc_actual, length_bound=B_com_rand, m=ringdim*A1_A2_width, norm=2, tag='Signature SIS')
com_sis = SIS.estimate.rough(sis_param)
print(com_sis)



print()
print("----------")
print("Computing sizes")
print("----------")
print()

# Compute sizes in bits

challenge = ringdim * 1
response = ringdim * (secdim) * ceil(log(B_inf,2))
hint = ringdim * moddim * ceil(log(B/(2**kappa_w * sqrt(ringdim * (moddim))),2))

signature = challenge + response + hint
pk = ringdim * moddim * (logq - kappa_y) + secpar

print("SIG size in KB:", N(signature/8000,digits=4))
print("PK size in KB:", N(pk/8000,digits=4))
print()

##### DKG Passive sizing sizes #######

h_i = 256
B_i = moddimhat * m * ringdim * log(q,2)
B_ij = moddimhat * m * ringdim * log(q,2)
S_ij = secdimhat * m * ringdim * log(q,2)
E_ij = moddimhat * m * ringdim * log(q,2)

DKG_comcostbits = h_i + B_i + n*B_ij + n*S_ij + n*E_ij 

DKG_comcost = (DKG_comcostbits/8000000).n()

print("DKG passive com cost=", DKG_comcost, "MB")
print()

##### Passive signing sizing sizes #######

 
hash_w_size = ringdim
w_size = moddim*ringdim*logq
ctx_r_size = (secdimhat+m)*ringdim*log_q_enc
decshare_size = m*ringdim*log_q_enc

decshare_size_bits=(decshare_size/8000).n()

print("DecShare_Cost =" f"{decshare_size_bits:.2f}", "KB")

ComCost_Passive = ((round(hash_w_size + w_size + ctx_r_size + decshare_size))/8000).n()
print()
print("Passive signing communication cost = ", f"{ComCost_Passive:.2f}", "KB")
print()


###### pi_dsi proof size ######
sig_pi_dsi = 2659
com_E = N((m+1)*ringdim*log(q_enc,2)/8000,digits=4)
#print("Com(E) in KB:", com_E)
z_E = (m+3)*ringdim*log(tail_bound*sig_pi_dsi,2)
#print("z_E cost in KB: ", N(z_E/8000,digits=4))
z_S = (m+3)* m *ringdim * log(tail_bound*sig_pi_dsi,2)
#print("z_S cost in KB: ", N(z_S/8000,digits=4))
pi_LIN_dsi = N((challenge + z_E + z_S)/8000,digits=4)
#print("pi_LIN_dsi cost in KB: ", pi_LIN_dsi)


# LNP proofs calculated externally using LaZer library

LNP_dsi=0
if secpar == 128: 
	LNP_dsi=89
if secpar == 192: 
	LNP_dsi=105.22
if secpar == 256: 
	LNP_dsi=120.22

LNP_ri=0
if secpar == 128: 
	LNP_ri=193
if secpar == 192: 
	LNP_ri=250.75
if secpar == 256: 
	LNP_ri=303

pi_dsi = com_E + pi_LIN_dsi + LNP_dsi
print("Final round proof cost in KB: ", pi_dsi)

print("Active signing cost =", f"{ComCost_Passive + pi_dsi + LNP_ri:.2f}", "KB") 
print("Online signing cost =", f"{decshare_size_bits + pi_dsi :.2f}", "KB") 
print("Online signing cost with LaBRADOR =", f"{decshare_size_bits + com_E + pi_LIN_dsi + 66 :.2f}", "KB") 
print("Active signing cost using LaBRADOR =", f"{ComCost_Passive + 66 + LNP_ri:.2f}", "KB") 
print("Optimistic signing cost =", f"{decshare_size_bits:.2f}", "KB")




print("Done")

