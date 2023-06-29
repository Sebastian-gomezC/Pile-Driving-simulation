from dolfin import *
from numpy.linalg import norm as np_norm
import matplotlib.pyplot as plt
def solve_linear_elasticity(mesh, boundaries, d):
    c = Constant(d)

    V = VectorFunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    E, nu = 10.0, 0.3
    mu = E/(2.0*(1.0 + nu))
    lmbda = E*nu/((1.0 + nu)*(1.0 -2.0*nu))
    def sigma(u):
        return 2*mu*sym(grad(u)) + lmbda*tr(grad(u))*Identity(2)
    F = inner(sigma(u), grad(v))*dx 
    a, L = lhs(F), rhs(F)

    bcs = [DirichletBC(V, Constant((0.0, 0.0)), boundaries, 1),
           DirichletBC(V.sub(0), c, boundaries, 2)]

    displacement = Function(V)
    solve(a==L, displacement, bcs)
    s = sigma(displacement) - (1./3)*tr(sigma(displacement))*Identity(2)
    von_Mises = sqrt(3./2*inner(s, s))
    V = FunctionSpace(mesh, 'P', 1)
    von_Mises = project(von_Mises, V)

    return displacement , von_Mises

def update_mesh(mesh, displacement, boundaries):

    new_mesh = Mesh(mesh)
    new_boundaries = MeshFunction("size_t", new_mesh, 1)
    new_boundaries.set_values(boundaries.array())
    ALE.move(new_mesh, displacement)
    print("1 good")
    return new_mesh, new_boundaries


# Original mesh
mesh = RectangleMesh(Point(0,0), Point(1,1),8,8)
subdomain1 = CompiledSubDomain("near(x[1], 0)")
subdomain2 = CompiledSubDomain("near(x[1], 1)")
boundaries = MeshFunction("size_t", mesh, 1)
boundaries.set_all(0)
subdomain1.mark(boundaries, 1)
subdomain2.mark(boundaries, 2)
plot(mesh)

# First iteration (accepted)
d=0
for i in range(100):
    d=0.01
    displacement,von = solve_linear_elasticity(mesh, boundaries, d)
    mesh, boundaries = update_mesh(mesh, displacement, boundaries)
    plot(von, title = "First iteration (accepted)")
    plt.savefig('ale1%s.png'%(i),dpi=300,)
    plt.close()
