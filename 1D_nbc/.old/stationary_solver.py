from pyoomph import *
from pyoomph.expressions import *
import numpy as np
import sys

ha = 0.5

if len(sys.argv) > 2:   # From cmd line
      h0out = float(sys.argv[1])
      Axmax = float(sys.argv[2])
else:
      sys.exit('Not enough input arguments: h0out, A(xmax)')

def wetting_pot(h):
        return (ha**3/5*h**(-5) - 1/2*h**(-2))
def theta(x,x0,sigma=0.1):
      return 0.5*(1+np.tanh((x-x0)/sigma))

class StationaryThinfilmEquation(Equations):
     def __init__(self,mu=1):
             super(StationaryThinfilmEquation, self).__init__()
             self.mu=mu # Viscosity (in the mobility Q)

     def dwetting_pot(self,h):
             return (- ha**3*h**(-6) + h**(-3))

     def define_fields(self):
             self.define_scalar_field("h","C2")

     def define_residuals(self):
             h,v=var_and_test("h")
             dWin = Axmax*self.dwetting_pot(h)
             dWout = Axmax*self.dwetting_pot(h0out)
             self.add_residual(-weak(grad(h),grad(v))-weak(dWin, v)+0.4486452)

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
        equations = StationaryThinfilmEquation()  # create a Poisson equation with source g=1
        equations += TextFileOutput()  # Add a simple text file output
        equations += SpatialErrorEstimator(h=1.0)
        #equations+=DirichletBC(h=h0out)@"left"
        #equations+=DirichletBC(h=h0out)@"right"
        equations += InitialCondition(h=h0out-0.3+0.3*np.cos(np.pi*2*x/(self.L)))
        sigma=self.L/4
        #equations += InitialCondition(h=h0out+0.1*exp(-x**2/(2*sigma**2)))
        #equations += PeriodicBC("right", offset=self.L) @ "left"

        # Measure gradh close to domain boundaries (and average)
        hhat = h-h0out
        whatin = wetting_pot(h)-wetting_pot(h0out)
        equations += IntegralObservables(Omega = whatin)
        equations += IntegralObservables(I = h0out**3*hhat/(h**3))
        equations += IntegralObservables(K = hhat**2/h**3)
        equations += IntegralObservableOutput(filename='obs')      # Stationary state
        self.add_equations(equations @ "domain")  # Add the equation system on the domain named "domain"
        # The maximum height and its x-coordinate (droplet position) are measured from domain*.txt


if __name__ == "__main__":
    with ThinfilmProblem(L=25,N=1000) as problem:

        # Maximum refinement level
        problem.max_refinement_level = 1
        problem.max_permitted_error = 0.0005
        problem.min_permitted_error = 0.00005
        problem.solve(spatial_adapt=problem.max_refinement_level)  # Solve the problem
        #problem.run(endtime=3000, startstep=1, outstep=True, temporal_error=1, spatial_adapt=problem.max_refinement_level)

        print(problem.get_mesh("domain").evaluate_all_observables())
        problem.output()  # Write output
