import tkinter
from tkinter import ttk
from tkinter import filedialog as fd

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import os
from dolfin import*
from fenics import*
from mshr import *
from ufl import nabla_div, nabla_grad
import numpy as np
from matplotlib import pyplot
def test_mesh(h,l):
    geo_test="""
    SetFactory("OpenCASCADE");
    h=%s;
    l=%s;
    H=60;
    Point(1) = {0,0,0};
    Point(2) = {0,20,0};
    Point(3) = {H/2-l/2,20,0};
    Point(4) = {H/2-l/2,20-h,0};
    Point(5) = {H/2+l/2,20-h,0};
    Point(6) = {H/2+l/2,20,0};
    Point(7) = {H,20,0};
    Point(8) = {H,0,0};
    
    Line(1) = {1, 2};
    
    Line(2) = {2, 3};
    
    Line(3) = {3, 4};
    //+
    Line(4) = {4, 5};
    //+
    Line(5) = {5, 6};
    //+
    Line(6) = {6, 7};
    //+
    Line(7) = {7, 8};
    //+
    Line(8) = {8, 1};
    //+
    Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
    //+
    Plane Surface(1) = {1};
    Physical Curve("a",1)={2};
    Physical Curve("bordes",3)={1,7};
    Physical Curve("neuman0",4)={3,4,5,8};
    Physical Curve("c",6)={6};
    Physical Surface("superficie",5)={1};
    
    Mesh 2;
    Coherence Mesh;
    
    RefineMesh;
    RefineMesh;
    RefineMesh;
    Mesh.MshFileVersion = 2.2;
    Save StrCat(StrPrefix(General.FileName), ".msh");
    """%(h,l)

    nombre = "proyecto_de_simulacion"
    if os.path.exists("%s.mesh_test"%(nombre)):
            a=os.path.exists("%s.mesh"%(nombre))
    else:
            os.mkdir("%s.mesh_test"%(nombre))
    with open("%s.mesh_test/%s.geo"%(nombre,nombre), 'w') as filed:
        filed.write(geo_test)
    #creacion de malla de prueba 
    os.system('gmsh -2  {}.mesh_test/{}.geo -format msh2'.format(nombre,nombre))
    os.system('dolfin-convert -i gmsh {}.mesh_test/{}.msh {}.mesh_test/{}.xml'.format(nombre, nombre,nombre,nombre))
    mesh=Mesh("%s.mesh_test/%s.xml"%(nombre,nombre))
    con = MeshFunction("size_t", mesh, "%s.mesh_test/"%(nombre)+nombre+"_facet_region.xml")
    sub= MeshFunction("size_t",mesh,"%s.mesh_test/"%(nombre)+nombre+"_physical_region.xml")
    return mesh,con,sub
def mallado(geo):

    nombre = "proyecto_de_simulacion"
    if os.path.exists("%s.mesh_extern"%(nombre)):
            a=os.path.exists("%s.mesh_extern"%(nombre))
    else:
            os.mkdir("%s.mesh_extern"%(nombre))
    with open("%s.mesh_extern/%s.geo"%(nombre,nombre), 'w') as filed:
        filed.write(geo)
    #creacion de malla de prueba 
    os.system('gmsh -2  {}.mesh_extern/{}.geo -format msh2'.format(nombre,nombre))
    os.system('dolfin-convert -i gmsh {}.mesh_extern/{}.msh {}.mesh_extern/{}.xml'.format(nombre, nombre,nombre,nombre))
    mesh=Mesh("{}.mesh_extern/{}.xml".format(nombre,nombre))
    con = MeshFunction("size_t", mesh, "%s.mesh_extern/"%(nombre)+nombre+"_facet_region.xml")
    sub= MeshFunction("size_t",mesh,"%s.mesh_extern/"%(nombre)+nombre+"_physical_region.xml")
    return mesh,con,sub
    
    return 
def solucionador(M,C,su,bs_in,bs_out,kx,ky):
    
    mesh=M
    contorno=C
    subdominio=su
    
    V = FunctionSpace(mesh, "CG", 2)
    u = TrialFunction(V)
    v = TestFunction(V)
    n=FacetNormal(mesh)
    #condiciones de contorno y constantes del problema 
    q= Constant(0)#flujo en las fronteras impermehables
    k= Constant(((kx,0),(0,ky)))#permeabilidad del suelo
    bc = DirichletBC(V, Constant((bs_in)),contorno,1)#potencial en la presa
    bc2 = DirichletBC(V, Constant((bs_out)),contorno,6)#potencial a la salida de la presa
    
    #-----------------------------------------------------------------------
    I = inner(k*grad(u),grad(v))*dx #problema variacional
    
    D=v*q*ds(subdomain_id=4, domain=mesh, subdomain_data=contorno) # condicion de flujo del problema variacional
    
    u = Function(V)
    solve( I==D,u,[bc,bc2]) #solucionador
    #______________________________________________________________________________
    #calculo vectores de flujo
    f=-k*grad(u)
    Vv=VectorFunctionSpace(mesh,"P",1)
    f=project(f,Vv)
    return u,f


Ventana = tkinter.Tk()
Ventana.title("Simulador de flujo para suelos saturados")
Ventana.geometry('1000x500')



Etiqueta0=tkinter.Label(Ventana,text="Geometria de la presa de muestra")
Etiqueta0.place(relx=0.1,rely=0.05)
Etiqueta1=tkinter.Label(Ventana,text="h(m)")
Etiqueta1.place(relx=0.07,rely=0.1)
caja1 = ttk.Entry(Ventana)
caja1.place(relx=0.13, rely=0.09, relwidth=0.2, relheight=0.05)

Etiqueta2=tkinter.Label(Ventana,text="l(m)")
Etiqueta2.place(relx=0.07,rely=0.2)
caja2 = ttk.Entry(Ventana)
caja2.place(relx=0.13, rely=0.19, relwidth=0.2, relheight=0.05)

#__________________________________________________________________________________

Etiqueta0=tkinter.Label(Ventana,text="Permeabilidad del suelo")
Etiqueta0.place(relx=0.1,rely=0.27)
Etiqueta3=tkinter.Label(Ventana,text="kx(m/s)")
Etiqueta3.place(relx=0.07,rely=0.35)
caja3 = ttk.Entry()
caja3.place(relx=0.13, rely=0.32, relwidth=0.2, relheight=0.05)
Etiqueta4=tkinter.Label(Ventana,text="ky(m/s)")
Etiqueta4.place(relx=0.07,rely=0.45)
caja4 = ttk.Entry()
caja4.place(relx=0.13, rely=0.42, relwidth=0.2, relheight=0.05)
#___________________________________________________________________________________________________________

Etiqueta0=tkinter.Label(Ventana,text="columna de agua antes de la presa")
Etiqueta0.place(relx=0.1,rely=0.55)


Etiqueta5=tkinter.Label(Ventana,text="H1(m)")
Etiqueta5.place(relx=0.07,rely=0.60)
caja5 = ttk.Entry()
caja5.place(relx=0.13, rely=0.58, relwidth=0.2, relheight=0.05)


Etiqueta0=tkinter.Label(Ventana,text="columna de agua despues de la presa")
Etiqueta0.place(relx=0.1,rely=0.65)


Etiqueta6=tkinter.Label(Ventana,text="H2(m)")
Etiqueta6.place(relx=0.07,rely=0.69)
caja6 = ttk.Entry()
caja6.place(relx=0.13, rely=0.68, relwidth=0.2, relheight=0.05)





#--------funciones---------


 
#recibir datos de cajas

def EjecutarGrafico():

    h = caja1.get()
    l = caja2.get()
    kx = caja3.get()
    ky = caja4.get()
    H1 = caja5.get()
    H2 = caja6.get()
    malla, contorno, subdominio = test_mesh(h, l)
    u, f = solucionador(malla, contorno, subdominio, H1, H2, kx, ky)
    p=plot(u)
    flux=plot(f)
    p.set_cmap ('RdBu')
    pyplot . colorbar ( p )
    pyplot.title('Campo de presión de poros y de flujo',loc='center')
    pyplot.show()
#generar grafico

#recibir datos de archivo
def AbrirArchivo():
    data=fd.askopenfilename(initialdir = "/",title = "Seleccione archivo",filetypes = (("geo files","*.geo"),("todos los archivos","*.*")))
    with open ("/home/sebastian/Desktop/file1.geo", "r") as myfile:
        geo=myfile.read()
    kx = caja3.get()
    ky = caja4.get()
    H1 = caja5.get()
    H2 = caja6.get()
    malla, contorno, subdominio = mallado(geo)
    u, f = solucionador(malla, contorno, subdominio, H1, H2, kx, ky)
    p=plot(u)
    flux=plot(f)
    p.set_cmap ('RdBu')
    pyplot . colorbar ( p )
    pyplot.title('Campo de presión de poros y de flujo',loc='center')
    pyplot.show()    

Etiqueta0=tkinter.Label(Ventana,text=""""Si desea ejecutar un caso externo puede dibujar 
la geometria con gmsh respetando los indices 
para las condiciones de contorno. """)
Etiqueta0.place(relx=0.05,rely=0.75)

boton1_ejecutar = ttk.Button(Ventana,text="abrir archivo" , command = AbrirArchivo)
boton1_ejecutar.place(relx=0.13, rely=0.85)

boton2_ejecutar = ttk.Button(Ventana,text="Ejecutar" , command = EjecutarGrafico)
boton2_ejecutar.place(relx=0.7, rely=0.60)
Etiqueta0=tkinter.Label(Ventana,text=""""El presente software nos permite calcular la presión de poros y el 
campo de flujo para un dominó determinado de una presa como se ilustra  
en el de la figura presentada u otro dominio cargado desde una geometría de gmsh""")
Etiqueta0.place(relx=0.45,rely=0.75)



Ventana.mainloop()