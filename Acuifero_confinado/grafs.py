import numpy as np
import matplotlib.pyplot as plt 

types=['BDF1','BDF2','BDF3']
dtdot=0.00101
L=[]
for scheme in types:
    L.append(np.loadtxt('resultados_Terzaghi_multi/L2norm_dt%s_grosse_%s.out'%(dtdot,scheme)))
    
ig, ax = plt.subplots()
for i in range(len(L)):
    line1, =ax.plot(L[i][:,1],np.sqrt(L[i][:,0]), "-",label=types[i])

ax.set_yscale('log') 
ax.set_xscale('log') 
plt.xlabel("$t*$ ",fontsize=20)
plt.ylabel("$L^2 error$",fontsize=20)
ax.legend(loc= 'upper right')
plt.grid(True,which="both",alpha=0.3, ls='-')
plt.savefig('resultados/L2norm_dt%s_grosse_%s.png'%(dtdot,'compare'),dpi=300)
