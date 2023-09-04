import numpy as np

M = 8
F = 1024
G = 2562

def real(F = 1):
    return  F*0

def real_addition(F = 1):
    return  F*1

def real_multiplication(F = 1):
    return  F*1

def real_inverse(F = 1):
    return  F*1

def real_negation(F = 1):
    return  F*1

def assignment(F = 1):
    return  F*0

def comparison(F = 1):
    return  F*0

# complex addition
    # 2 real additions
def complex_addition(F):
    return F*2*real_addition()

# complex multiplication
    # 4 real multiplications
    # 2 real additions
    # 2 assignments
def complex_multiplication(F):
    return  F*(
            4*real_multiplication() \
        +   2*real_addition() \
        +   2*assignment()
    )

# complex magnitude
    # 4 real multiplications
    # 1 real addition
    # 1 square-root
    # 1 assignment
def complex_magnitude(F):
    return  F*(
            4*real_multiplication() \
        +   1*real_addition() \
        +   1*assignment()
    )

# complex conjugate
    # 1 negation
    # 2 assignments
def complex_conjugate(F):
    return  F*(
            1*real_negation()
        +   2*assignment() 
    )

# F/2*log_2(F) complex additions
# F*log_2(F) complex multiplications
def FFT(F):
    return  F/2*np.log2(F) * complex_addition(1) \
        +   F*np.log2(F) * complex_multiplication(1)
def inverse_FFT(F):
    return FFT(F)

# set ψ = 1 / (mag(X_0) * mag(X_1))
# set χ = multiply(
#     ψ, 
#     multiply(
#         X_0, 
#         conj(X_1) 
#     ) 
# ) 
# set R_01 = real(invFFT(χ))
def compute_GCC (F):
    return  2*complex_magnitude(F) \
        +   1*real_multiplication(F) \
        +   1*real_inverse(F) \
        +   1*assignment(F) \
        +   1*complex_conjugate(F) \
        +   2*complex_multiplication(F) \
        +   1*assignment(F) \
        +   1*inverse_FFT(F) \
        +   1*real(F) \
        +   1*assignment(F)

# E = 0
# for (set i = 0; i < M; i++) do
#     for (set j = 0; j < M; j++) do
#         if i == j do
#             continue
#         end
#         set E = E + R[i,j,τ[g,i,j]]
#     end
# end
def compute_beam_energy(F):
    count = 0
    for i in range(M):
        for j in range(M):
            count = count + comparison()
            if i == j:
                continue
            count = count \
                + real_addition() \
                + assignment()
    return count

# foreach (x_i in x) do
#     set X = FFT(x)
# end

# for (set i = 0; i < M; i++) do
#     for (set j = 0; j < M; j++) do
#         if i == j do
#             continue
#         end
#         set R[i,j,:] = compute_GCC(X[:,i], X[:,j])
#     end
# end

# set g_max = 0
# set E_max = 0
# for (set g = 0; g < G; g++) do
#     set E = compute_beam_energy(x, g, R)
#     if E >= E_max do
#         set g_max = g
#         set E_max = E
#     end
# end
def main():
    count = 0
    old_count = count
    for i in range(M):
        count = count + FFT(F)
    print("FFT did {} FLOs.".format(count - old_count))

    old_count = count
    for i in range(M):
        for j in range(M):
            if i == j:
                continue
            count = count + compute_GCC(F) + assignment(F)
    print("GCC did {} FLOs.".format(count - old_count))
    
    old_count = count
    for g in range(G):
        count = count + compute_beam_energy(F) + assignment()
        count = count + comparison()
        count = count + 2*assignment()
    print("Beams did {} FLOs.".format(count - old_count))
    print("All did {} FLOs".format(count))
    return count

f_s = 48*10**3
T_F = F/f_s
print("Frame-period: {} s".format(T_F))
f_I = 900*10**6
n_I = f_I*T_F
n_FLO = n_I/10
print("Number of FLO: {}".format(n_FLO))

count = main()
f_FLO = count/T_F
print("{} FLOPS".format(f_FLO))