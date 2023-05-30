import math
from mpmath import mpf, besseli, besselk
import numpy as np
import matplotlib.pyplot as plt

cv = 1
time = 1
R = 1
steps = 10
Betta = 1
t = 1
G = 1
K=10
q = 10
M=10
alpha=1
r = np.linspace(0.01, R, steps)
def pbar(r,s,f_const):
    sigma = (s * R ** 2) / cv
    ro = r / R
      # Crear un arreglo para almacenar los resultados este se inicializa en cero, es como crear una lista pero array
    I_cero_ro = besseli(0, ro * (sigma) ** (1 / 2))
    K_cero_ro = besselk(0, ro * (sigma) ** (1 / 2))
    I_cero = besseli(0, (sigma) ** (1 / 2))
    K_cero = besselk(0, (sigma) ** (1 / 2))
    term_1 = (Betta * G * f_const / q) * (1 - I_cero_ro / I_cero)
    term_2 = (1 / sigma) * K_cero_ro
    term_3 = (K_cero / I_cero) * I_cero_ro
    v_ray = (term_1 - term_2 - term_3)
    return v_ray
def f_coef(r,s,f_const):
    sigma = (s * R ** 2) / cv
    part_1 = (2*alpha/sigma**2)*(1-1/(besseli(0,np.sqrt(sigma))))
    part_2 = 3*K/G + alpha*Betta*(1-(2*besseli(1,np.sqrt(sigma)))/(np.sqrt(sigma)*besseli(0,np.sqrt(sigma))))
    return (part_1/part_2)*q/G

def talbot(M, t,Fun,r=0,f_const=0):
    v=0
    for k in range(M - 1):
        if k == 0:
            delta =2 * M / 5
            gamma = np.exp(delta) / 2
        else:
            cot = 1 / np.tan(k * np.pi / M)
            delta = (2 * k * np.pi / 5) * (cot + 1j)
            gamma = (1 + 1j * (k * np.pi / M) *(1+ (cot) ** 2) - 1j* cot) * np.exp(delta)
        term_fun = delta / t
        termino = (gamma * Fun(r,term_fun,f_const))
        v += termino.real
    v_final=((2) / (5 * t)) * v
    return v_final
def analitical(r, R, G, q, cv,M,t):
    v_final = []
    f_const = talbot(M, t, f_coef)
    print(f_const)
    for l in range(len(r)):   
        v = 0
        v_final.append([talbot(M, t, pbar,r[l],f_const)])
    return np.array(v_final)
t=0.1

for i in range(3):
    resultado = analitical(r, R, G, q, cv,M,t)
    plt.plot(r, resultado)
    t +=1000000
plt.plot(r,np.log(r/R))


