from pyoomph import *
from pyoomph.expressions import *
import numpy as np
import sys

# Switch ON/OFF the driving forces
switch_bulk = 0
switch_inhomo = 1
# Parameters (additional parameters are in the initial condition)
ha = 0.5
eps = 0.001

if len(sys.argv) > 1:   # Set epsilon from cmd line
      eps = float(sys.argv[1])

def wetting_pot(h):
        return (ha**3/5*h**(-5) - 1/2*h**(-2))
def theta(x,x0,sigma=0.1):
      return 0.5*(1+np.tanh((x-x0)/sigma))

class ThinfilmEquation(Equations):
     def __init__(self,mu=1):
             super(ThinfilmEquation, self).__init__()
             self.mu=mu # Viscosity (in the mobility Q)
             t = var("time")
             # Activate bulk force and inhomogeneous Hawmaker constant for t>t0
             t0=1e3
             activation_fun = theta(t,t0)
             # Select wether activate c, A or both drivings
             amplitude_C = 1
             amplitude_A = 1
             if switch_bulk > 0:
                self.c = activation_fun*amplitude_C
             else:
                self.c = 0
             if switch_inhomo > 0:
                self.A0 = activation_fun*amplitude_A
             else:
                self.A0 = 0  


     def wetting_pot(self,h):
             return (ha**3/5*h**(-5) - 1/2*h**(-2))
     def dwetting_pot(self,h):
             return (- ha**3*h**(-6) + h**(-3))
     # Bulk force on x axis
     def bulk_force(self,x,h):
             return -eps*self.c*x*h
     def dbulk_force(self,x,h):
             return -eps*self.c*x

     def define_fields(self):
             self.define_scalar_field("h","C2")
             self.define_scalar_field("p","C2")

     def define_residuals(self):
             h,v=var_and_test("h")
             p,eta=var_and_test("p")
             x = var("coordinate_x")
             Q=(h**3)/(3*self.mu) # Mobility
             A = 1+self.A0*eps*x  # Hawmaker constant, linearly along y-direction
             dW = A*self.dwetting_pot(h) + self.dbulk_force(x,h)  # W'(h)

             self.add_residual(weak(partial_t(h),v)+weak(Q*grad(p),grad(v)))
             self.add_residual(weak(p,eta)-weak(grad(h),grad(eta))-weak(dW,eta))

class ThinfilmProblem(Problem):
    def __init__(self, L=20, N=1000):
        super().__init__()
        self.L = L
        self.N = N
    def define_problem(self):
        mesh = LineMesh(minimum=-self.L/2, size=self.L, N=self.N)
        # Add the mesh (default name is "domain" with boundaries "left" and "right")
        self.add_mesh(mesh)
        x = var("coordinate_x")
        h = var('h')

        # Assemble the system
        equations = ThinfilmEquation()  # create a Poisson equation with source g=1
        equations += TextFileOutput()  # Add a simple text file output
        equations += InitialCondition(h=1+0.5*np.cos(np.pi*2*x/self.L))
        equations += SpatialErrorEstimator(h=1.0)
        #equations += PeriodicBC("right", offset=self.L) @ "left"

        # Measure gradh close to domain boundaries (and average)
        toll = 0.95
        boundary_function_left = theta(x,self.L/2*toll)
        boundary_function_right = theta(-x, self.L/2*toll)
        equations += IntegralObservables(gradh = (grad(h)*boundary_function_left+grad(h)*boundary_function_right)/(self.L*(1-toll)))
        equations += IntegralObservableOutput(filename='obs')      # Stationary state
        self.add_equations(equations @ "domain")  # Add the equation system on the domain named "domain"
        # The maximum height and its x-coordinate (droplet position) are measured from domain*.txt


if __name__ == "__main__":
    with ThinfilmProblem(L=25,N=1000) as problem:

        # Maximum refinement level
        problem.max_refinement_level = 1
        problem.max_permitted_error = 0.0005
        problem.min_permitted_error = 0.00005
        problem.run(endtime=10000, startstep=1, outstep=True, temporal_error=1, spatial_adapt=problem.max_refinement_level)

        #problem.solve()  # Solve the problem
        print(problem.get_mesh("domain").evaluate_all_observables())
        problem.output()  # Write output
