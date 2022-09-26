import numpy as np 
from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
R =0.1 #radio de la tuberia
N_s= [4,20,30,50,100,200,300,400,500,1000]
E=np.zeros(10)
j=0
for d in N_s:
    
    N=d#numero de volumenes usados
    dr=(R)/N # usamos celdas uniformes 
    mu=0.001
    Px=1
    r_i=np.zeros(N) # posición del centro de cada celda 
    
    def r_f(n,posicion): #con esta función hallamos la posicion de las caras este u oeste para cada celda n
        if posicion == 'w':
            return(r_i[n]+dr/2)
        elif posicion == 'e':
            return(r_i[n]-dr/2)
    for k in range(N):
            r_i[k]=dr*(k+1)-dr/2
    #definimos los coeficeintes de la forma discreta p
    def aw(n):
        return(r_f(n,'w'))
    def ap(n):
        if n==0:
            return(-(r_f(n,'w')+r_f(n,'e')))
        elif n==N-1 :
            return(-(2*r_f(n,'w')+r_f(n,'e')))
        else :
            return(-(r_f(n,'w')+r_f(n,'e')))
    def ae(n):
        if n== N:
            return(0.0)
        else:
            return(r_f(n, 'e'))
    #definimos los coeficeintes del lado derecho
    def S(n):
        return((-Px/mu)*r_i[n]*dr**2)
    #creamos la matriz tridiagonal
    Tdm=np.zeros((N,N))
    for n in range(N):
        if n!=0:
            Tdm[n][n-1]=ae(n)
        Tdm[n][n]=ap(n)
        if n!=N-1:
            Tdm[n][n+1]=aw(n)
    #creamos el vector de terminos fuentes
    So=np.zeros(N)
    
    for i in range(N):
        So[i]=S(i)
    #calculamos la matriz inversa y multiplicamos por el vector de terminos fuente
    
    U=np.dot(np.linalg.inv(Tdm),So)
    
    #solucion analitica 
    y=np.linspace(0, 0.1,100)
    
    def u_a(r):
        return((Px/(4*mu))*(R**2-r**2))
    e=0
    for i in range(N):
        e=e+abs((U[i]-u_a(r_i[i]))/u_a(r_i[i]))
    E[j]=e
    j=j+1
    
    print(e)
    fig, ax = plt.subplots()
    line1, = ax.plot(y,u_a(y),'-', label='Analitica')
    line2, = ax.plot(r_i,U,'--', label='FVM',color='red')
    ax.legend(handler_map={line1: HandlerLine2D(numpoints=4)})
    plt.title("Perfiles de velocidad ")
    plt.xlabel("posición radial r")
    plt.ylabel("velocidad axial m/s")
    plt.savefig('flow_pipe_n %s.png'%(d))
    plt.close()
plt.loglog( N_s,E)

plt.title("")
plt.xlabel("Numero de celdas")
plt.ylabel("Error relativo")