import galois
import numpy as np

# Part 2
print("----------------------------------------- Part 2 -----------------------------------------")

# p = 71
p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
FP = galois.GF(p)

# input arguments
x = FP(2)
y = FP(3)

v1 = x * x
v2 = y * y
v3 = 5 * x * v1
v4 = 4 * v1 * v2
out = 5*x**3 - 4*x**2*y**2 + 13*x*y**2 + x**2 - 10*y

w = FP([1, out, x, y, v1, v2, v3, v4])

print("w =", w)

R = FP([[0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 5, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 4, 0, 0, 0],
         [0, 0, 13, 0, 0, 0, 0, 0]])

L = FP([[0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0]])

O = FP([[0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 1],
         [0, 1, 0, 10, FP(p - 1), 0, FP(p - 1), 1]])

Lw = np.dot(L, w)
Rw = np.dot(R, w)
Ow = np.dot(O, w)

print("Lw =", Lw)
print("Rw =", Rw)

LwRw = np.multiply(Lw, Rw)

print("Lw * Rw =", LwRw)

print("Ow =     ", Ow)

assert np.all(LwRw == Ow)

# Part 3
print("----------------------------------------- Part 3 -----------------------------------------")

mtxs = [L, R, O]
poly_m = []

for m in mtxs:
    poly_list = []
    for i in range(0, m.shape[1]):
        points_x = FP(np.zeros(m.shape[0], dtype=int)) # [0, 0, 0, 0, 0]
        points_y = FP(np.zeros(m.shape[0], dtype=int)) #[0, 0, 0, 0, 0]
        for j in range(0, m.shape[0]):
            points_x[j] = FP(j+1)
            points_y[j] = m[j][i]

        # points_x = [1, 2, 3, 4, 5]
        # points_y = [1, 0, 5, 0, 13]
        poly = galois.lagrange_poly(points_x, points_y)
        coef = poly.coefficients()[::-1] # reverse
        if len(coef) < m.shape[0]:
            coef = np.append(coef, np.zeros(m.shape[0] - len(coef), dtype=int))
        poly_list.append(coef)

    poly_m.append(FP(poly_list))

Lp = poly_m[0]
Rp = poly_m[1]
Op = poly_m[2]

print(f'''Lp
{Lp}
''')

print(f'''Rp
{Rp}
''')

print(f'''Op
{Op}
''')

U = galois.Poly((w @ Lp)[::-1])
V = galois.Poly((w @ Rp)[::-1])
W = galois.Poly((w @ Op)[::-1])

print("U = ", U)
print("V = ", V)
print("W = ", W)

tau = FP(20)
# assert U(tau)*V(tau) == W(tau) # will fail
print(f"τ = {tau}")

T = galois.Poly([1, p-1], field=FP)  #T(x) = x - 1
for i in range(2, L.shape[0] + 1):
    T *= galois.Poly([1, p-i], field=FP)

# T(x) = (x-1)(x-2)(x-3)(x-4)(x-5)
print("\nT = ", T)

# T(1)=0, T(2)=0, T(3)=0, T(4)=0, T(5)=0, T(6)=49
for i in range(1, L.shape[0] + 2):
    print(f"T({i}) = ", T(i))
    if i == L.shape[0]:
        print("-"*10)
        
T_tau = T(tau)
print(f"\nT(τ) = {T_tau}")

H = (U * V - W) // T
rem = (U * V - W) % T

print("H = ", H)
print("rem = ", rem)

assert rem == 0

u = U(tau)
v = V(tau)
_w = W(tau) # sorry, w is taken by witness vector
ht = H(tau)*T_tau

assert u * v - _w == ht, f"{u} * {v} - {_w} != {ht}" # this equation should hold

# # # Part 4
print("----------------------------------------- Part 4 -----------------------------------------")
print("Setup phase")
print("-"*10)
print("Toxic waste:")
tau = FP(20)

print(f"τ = {tau}")

T_tau = T(tau)
print(f"\nT(τ) = {T_tau}")

# ...
from py_ecc.optimized_bn128 import multiply, G1, G2, add, pairing, neg, normalize

# G1[τ^0], G1[τ^1], ..., G1[τ^d]
tau_G1 = [multiply(G1, int(tau**i)) for i in range(0, T.degree)]
# G1[τ^0 * T(τ)], G1[τ^1 * T(τ)], ..., G1[τ^d-1 * T(τ)]
target_G1 = [multiply(G1, int(tau**i * T_tau)) for i in range(0, T.degree - 1)]

# G2[τ^0], G2[τ^1], ..., G2[τ^d]
tau_G2 = [multiply(G2, int(tau**i)) for i in range(0, T.degree)]

print("Trusted setup:")
print("-"*10)
print(f"[τ]G1 = {[normalize(point) for point in tau_G1]}")
print(f"[T(τ)]G1 = {[normalize(point) for point in target_G1]}")

print(f"\n[τ]G2 = {[normalize(point) for point in tau_G2]}")

def evaluate_poly(poly, trusted_points, verbose=False):
    coeff = poly.coefficients()[::-1]

    assert len(coeff) == len(trusted_points), "Polynomial degree mismatch!"

    if verbose:
        [print(normalize(point)) for point in trusted_points]

    # zip: pairs each element in trusted_points with the corresponding coefficient in coeff.
    # For each pair (point, coeff), it computes multiply(point, int(coeff))
    # This creates a list terms where each element is point * coeff.
    terms = [multiply(point, int(coeff)) for point, coeff in zip(trusted_points, coeff)]  
    evaluation = terms[0]
    for i in range(1, len(terms)):
        # Accumulate the Sum
        evaluation = add(evaluation, terms[i])

    if verbose:
        print("-"*10)
        print(normalize(evaluation))
    return evaluation

print("\nProof generation:")
print("-"*10)
# G1[u0 * τ^0] + G1[u1 * τ^1] + ... + G1[ud-1 * τ^d-1]
A_G1 = evaluate_poly(U, tau_G1)
# G2[v0 * τ^0] + G2[v1 * τ^1] + ... + G2[vd-1 * τ^d-1]
B_G2 = evaluate_poly(V, tau_G2)
# # G1[w0 * τ^0] + G1[w1 * τ^1] + ... + G1[wd-1 * τ^d-1]
# B_G1 = evaluate_poly(V, tau_G1)
# G1[w0 * τ^0] + G1[w1 * τ^1] + ... + G1[wd-1 * τ^d-1]
Cw_G1 = evaluate_poly(W, tau_G1)
# G1[h0 * τ^0 * T(τ)] + G1[h1 * τ^1 * T(τ)] + ... + G1[hd-2 * τ^d-2 * T(τ)]
HT_G1 = evaluate_poly(H, target_G1)

C_G1 = add(Cw_G1, HT_G1)

print(f"[A]G1 = {normalize(A_G1)}")
print(f"[B]G2 = {normalize(B_G2)}")
print(f"[C]G1 = {normalize(C_G1)}")

print("\nProof verification:")
print("-"*10)
# e(A, B) == e(C, G2[1])
assert pairing(B_G2, A_G1) == pairing(G2, C_G1)
print("Pairing check passed!")

# Part 5
print("----------------------------------------- Part 5 -----------------------------------------")

# considering: (A + α)(B + β) = αβ + (βA + αB + C)  instead of A * B = C
print("Adding alpha and beta to the equation to determine if A,B,C are genuine points in results of an evaluation. Prover cannot cheat!")
def evaluate_poly_list(poly_list, x):
    results = []
    for poly in poly_list:
        results.append(poly(x))
    return results

def print_evaluation(name, results):
    print(f'\n{name} polynomial evaluations:')
    for i in range(0, len(results)):
        print(f'{name}_{i} = {results[i]}')

# converts each row of a matrix into a polynomial object
def to_poly(mtx):
    poly_list = []
    for i in range(0, mtx.shape[0]):
        poly_list.append( galois.Poly(mtx[i][::-1]) )
    return poly_list

def print_poly(name, poly_list):
    print(f'\n{name} polynomials:')
    for i in range(0, len(poly_list)):
        print(f'{name}_{i} = {poly_list[i]}')

print("Setup phase")
print("-"*10)
print("Toxic waste:")
alpha = FP(2)
beta = FP(3)
tau = FP(20)

print(f"α = {alpha}")
print(f"β = {beta}")
print(f"τ = {tau}")

beta_L = beta * Lp
alpha_R = alpha * Rp
K = beta_L + alpha_R + Op # preimage of [βA + αB + C]

Kp = to_poly(K)
print_poly("K", Kp)

print("K evaluations:")
K_eval = evaluate_poly_list(Kp, tau)
print([int(k) for k in K_eval])

from py_ecc.optimized_bn128 import multiply, G1, G2, add, pairing, neg, normalize, eq, is_on_curve

# G1[α]
alpha_G1 = multiply(G1, int(alpha))
# G1[β]
beta_G1 = multiply(G1, int(beta))
# G1[τ^0], G1[τ^1], ..., G1[τ^d]
tau_G1 = [multiply(G1, int(tau**i)) for i in range(0, T.degree)]
# G1[βU0(τ) + αV0(τ) + W0(τ)], G1[βU1(τ) + αV1(τ) + W1(τ)], ..., G1[βUd(τ) + αVd(τ) + Wd(τ)]
k_G1 = [multiply(G1, int(k)) for k in K_eval]
# G1[τ^0 * T(τ)], G1[τ^1 * T(τ)], ..., G1[τ^d-1 * T(τ)]
target_G1 = [multiply(G1, int(tau**i * T_tau)) for i in range(0, T.degree - 1)]

# G2[α]
beta_G2 = multiply(G2, int(beta))
# G2[τ^0], G2[τ^1], ..., G2[τ^d-1]
tau_G2 = [multiply(G2, int(tau**i)) for i in range(0, T.degree)]

print("Trusted setup:")
print("-"*10)
print(f"[α]G1 = {normalize(alpha_G1)}")
print(f"[β]G1 = {normalize(beta_G1)}")
print(f"[τ]G1 = {[normalize(point) for point in tau_G1]}")
print(f"[k]G1 = {[normalize(point) for point in k_G1]}")
print(f"[τT(τ)]G1 = {[normalize(point) for point in target_G1]}")

print(f"\n[β]G2 = {normalize(beta_G2)}")
print(f"[τ]G2 = {[normalize(point) for point in tau_G2]}")

print("Prover's Part")
print("\nProof generation:")
print("-"*10)

U = galois.Poly((w @ Lp)[::-1])
V = galois.Poly((w @ Rp)[::-1])

# G1[u0 * τ^0] + G1[u1 * τ^1] + ... + G1[ud-1 * τ^d-1]
A_G1 = evaluate_poly(U, tau_G1)
# G1[A] = G1[A] + G1[α]
A_G1 = add(A_G1, alpha_G1)
# G2[v0 * τ^0] + G2[v1 * τ^1] + ... + G2[vd-1 * τ^d-1]
B_G2 = evaluate_poly(V, tau_G2)
# G2[B] = G2[B] + G2[β]
B_G2 = add(B_G2, beta_G2)
# G1[h0 * τ^0 * T(τ)] + G1[h1 * τ^1 * T(τ)] + ... + G1[hd-2 * τ^d-2 * T(τ)]
HT_G1 = evaluate_poly(H, target_G1)
assert len(w) == len(k_G1), "Polynomial degree mismatch!"
# w0 * G1[k0] + w1 * G1[k1] + ... + wd-1 * G1[kd-1]
K_G1_terms = [multiply(point, int(scaler)) for point, scaler in zip(k_G1, w)]
K_G1 = K_G1_terms[0]
for i in range(1, len(K_G1_terms)):
    K_G1 = add(K_G1, K_G1_terms[i])

C_G1 = add(HT_G1, K_G1)

print(f"[A]G1 = {normalize(A_G1)}")
print(f"[B]G2 = {normalize(B_G2)}")
print(f"[C]G1 = {normalize(C_G1)}")
print("-" * 10)
print("Verifier uses:")
print(f"[α]G1 = {normalize(alpha_G1)}")
print(f"[β]G1 = {normalize(beta_G1)}")

print("Verifier's Part")
# A = A + α
# B = B + β
# C = βA + αB + C
# AB == αβ + [βA + αB + C]
# assert pairing(B_G2, A_G1) == pairing(beta_G1, alpha_G1) + pairing(G2, C_G1)

print("Public and private inputs separation")
def split_poly(poly):
    coef = [int(c) for c in poly.coefficients()]
    p1 = coef[-2:] # [1, out]
    p2 = coef[:-2] + [0] * 2 

    return galois.Poly(p1, field=FP), galois.Poly(p2, field=FP)

u = U(tau)
v = V(tau)
_w = W(tau) # w taken by witness vector
ht = H(tau)*T_tau

U1, U2 = split_poly(U)
V1, V2 = split_poly(V)
W1, W2 = split_poly(W)

w1 = W1(tau)
w2 = W2(tau)

u1 = U1(tau)
u2 = U2(tau)

v1 = V1(tau)
v2 = V2(tau)

c = (beta * u2 + alpha * v2 + w2) + ht  # c = [βU+αV+W+HT] 
k = (beta * u1 + alpha * v1 + w1) # k = [βU+αV+W] 

assert (u + alpha) * (v + beta) == alpha * beta + k + c # should be equal

k_pub_G1, k_priv_G1 = k_G1[:2], k_G1[2:]
pub_input, priv_input = w[:2], w[2:]

print(f"[k_pub]G1 = {[normalize(point) for point in k_pub_G1]}")
print(f"[k_priv]G1 = {[normalize(point) for point in k_priv_G1]}")

print(f"pub_input = {pub_input}")
print(f"priv_input = {priv_input}")

print("Prover's part")
K_priv_G1_terms = [multiply(point, int(scaler)) for point, scaler in zip(k_priv_G1, priv_input)]
K_priv_G1 = K_priv_G1_terms[0]
for i in range(1, len(K_priv_G1_terms)):
    K_priv_G1 = add(K_priv_G1, K_priv_G1_terms[i])

C_G1 = add(HT_G1, K_priv_G1)
print(f"C_G1 = {C_G1}")

# # Part 6
print("----------------------------------------- Part 6 -----------------------------------------")
# Adding gamma and delta
print("Setup phase")
print("-"*10)
print("Toxic waste:")
alpha = FP(2)
beta = FP(3)
gamma = FP(4)
delta = FP(5)
tau = FP(20)

print(f"α = {alpha}")
print(f"β = {beta}")
print(f"γ = {gamma}")
print(f"δ = {delta}")
print(f"τ = {tau}")

def split_poly(poly):
    coef = [int(c) for c in poly.coefficients()]
    p1 = coef[-2:]
    p2 = coef[:-2] + [0] * 2

    return galois.Poly(p1, field=FP), galois.Poly(p2, field=FP)

u = U(tau)
v = V(tau)
ht = H(tau)*T_tau

U1, U2 = split_poly(U)
V1, V2 = split_poly(V)
W1, W2 = split_poly(W)

w1 = W1(tau)
w2 = W2(tau)

u1 = U1(tau)
u2 = U2(tau)

v1 = V1(tau)
v2 = V2(tau)

c = (beta * u2 + alpha * v2 + w2) * delta**-1 + ht * delta**-1
k = (beta * u1 + alpha * v1 + w1) * gamma**-1

a = u + alpha
b = v + beta

assert a * b == alpha * beta + k * gamma + c * delta # should be equal.

tau_G1 = [multiply(G1, int(tau**i)) for i in range(0, T.degree)]
tau_G2 = [multiply(G2, int(tau**i)) for i in range(0, T.degree)]

powers_tauTtau_div_delta = [(tau**i * T_tau) / delta for i in range(0, T.degree - 1)]
target_G1 = [multiply(G1, int(pTd)) for pTd in powers_tauTtau_div_delta]

assert len(target_G1) == len(H.coefficients()), f"target_G1 length mismatch! {len(target_G1)} != {len(H.coefficients())}"

K_gamma, K_delta = [k/gamma for k in K_eval[:2]], [k/delta for k in K_eval[2:]]
K_gamma_G1 = [multiply(G1, int(k)) for k in K_gamma]
K_delta_G1 = [multiply(G1, int(k)) for k in K_delta]

# G2[γ]
gamma_G2 = multiply(G2, int(gamma))
# G2[δ]
delta_G1 = multiply(G1, int(delta))
# G2[δ]
delta_G2 = multiply(G2, int(delta))

print("Trusted setup:")
print("-"*10)
print(f"[α]G1 = {normalize(alpha_G1)}")
print(f"[β]G2 = {normalize(beta_G2)}")
print(f"[γ]G2 = {normalize(gamma_G2)}")
print(f"[δ]G2 = {normalize(delta_G2)}")
print(f"[τ]G1 = {[normalize(point) for point in tau_G1]}")
print(f"[τ]G2 = {[normalize(point) for point in tau_G2]}")
print(f"[τT(τ)/δ]G1 = {[normalize(point) for point in target_G1]}")
print(f"[K/γ]G1 = {[normalize(point) for point in K_gamma_G1]}")
print(f"[K/δ]G1 = {[normalize(point) for point in K_delta_G1]}")

# Prover
A_G1 = evaluate_poly(U, tau_G1)
A_G1 = add(A_G1, alpha_G1)
B_G2 = evaluate_poly(V, tau_G2)
B_G2 = add(B_G2, beta_G2)
HT_G1 = evaluate_poly(H, target_G1)

w_pub, w_priv = w[:2], w[2:]

# [K/δ*w_priv]G1
Kw_delta_G1_terms = [multiply(point, int(scaler)) for point, scaler in zip(K_delta_G1, w_priv)]
Kw_delta_G1 = Kw_delta_G1_terms[0]
for i in range(1, len(Kw_delta_G1_terms)):
    Kw_delta_G1 = add(Kw_delta_G1, Kw_delta_G1_terms[i])

C_G1 = add(Kw_delta_G1, HT_G1)

# Adding r and s
print("Prover picks random r, s")
r = FP(12) # should be random
s = FP(13) # should be random

print(f"r = {r}")
print(f"s = {s}")

def split_poly(poly):
    coef = [int(c) for c in poly.coefficients()]
    p1 = coef[-2:]
    p2 = coef[:-2] + [0] * 2

    return galois.Poly(p1, field=FP), galois.Poly(p2, field=FP)

u = U(tau)
v = V(tau)
ht = H(tau)*T_tau

U1, U2 = split_poly(U)
V1, V2 = split_poly(V)
W1, W2 = split_poly(W)

w1 = W1(tau)
w2 = W2(tau)

u1 = U1(tau)
u2 = U2(tau)

v1 = V1(tau)
v2 = V2(tau)

a = u + alpha + r * delta
b = v + beta + s * delta

c = ((beta * u2 + alpha * v2 + w2) * delta**-1 + ht * delta**-1) + s * a + r * b - r * s * delta
k = (beta * u1 + alpha * v1 + w1) * gamma**-1

assert a * b == alpha * beta + k * gamma + c * delta # should be equal

r_delta_G1 = multiply(delta_G1, int(r))
s_delta_G1 = multiply(delta_G1, int(s))
s_delta_G2 = multiply(delta_G2, int(s))

A_G1 = evaluate_poly(U, tau_G1)
A_G1 = add(A_G1, alpha_G1)
A_G1 = add(A_G1, r_delta_G1)

B_G2 = evaluate_poly(V, tau_G2)
B_G2 = add(B_G2, beta_G2)
B_G2 = add(B_G2, s_delta_G2)

B_G1 = evaluate_poly(V, tau_G1)
B_G1 = add(B_G1, beta_G1)
B_G1 = add(B_G1, s_delta_G1)

As_G1 = multiply(A_G1, int(s))
Br_G1 = multiply(B_G1, int(r))
rs_delta_G1 = multiply(delta_G1, int(-r*s))

C_G1 = add(Kw_delta_G1, HT_G1)
C_G1 = add(C_G1, As_G1)
C_G1 = add(C_G1, Br_G1)
C_G1 = add(C_G1, rs_delta_G1)

print(f"[A]G1 = {normalize(A_G1)}")
print(f"[B]G2 = {normalize(B_G2)}")
print(f"[C]G1 = {normalize(C_G1)}")