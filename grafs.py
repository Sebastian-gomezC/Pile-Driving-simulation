import numpy as np
import matplotlib.pyplot as plt 

types=['BDF1','BDF2','BDF2op','BDF3']
dtdot=0.00826
L=[]
for scheme in types:
    L.append(np.loadtxt('resultados/L2norm_dt%s_fine_%s.out'%(dtdot,scheme)))
    
ig, ax = plt.subplots()
for i in range(len(L)):
    line1, =ax.plot(L[i][:,1],L[i][:,0], "-",label=types[i])

ax.set_yscale('log') 
ax.set_xscale('log') 
plt.xlabel("$t*$ ")
plt.ylabel("$L^2 error$")
ax.legend(loc= 'upper right')
plt.savefig('resultados/L2norm_dt%s_fine_%s.png'%(dtdot,'compare'),dpi=300)