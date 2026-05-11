from pyoomph import *
from pyoomph.expressions import *
import numpy as np
import sys

# Switch ON/OFF the driving forces
switch_bulk = 0
switch_inhomo = 1
# Parameters (additional parameters are in the initial condition)
ha = 0.5
eps = 0.005

if len(sys.argv) > 1:   # Set epsilon from cmd line
      eps = float(sys.argv[1])

      
def wetting_pot(h):
        return (ha**3/5*h**(-5) - 1/2*h**(-2))

class ThinfilmEquation(Equations):
     def __init__(self,mu=1):
             super(ThinfilmEquation, self).__init__()
             self.mu=mu # Viscosity (in the mobility Q)
             t = var("time")
             # Activate bulk force and inhomogeneous Hawmaker constant for t>t0
             t = var("time")
             # Activate bulk force and inhomogeneous Hawmaker constant for t>t0
             t0=1e3
             activation_fun = 0.5*(1+np.tanh(t-t0))
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
             y = var("coordinate_y")
             Q=(h**3)/(3*self.mu) # Mobility
             A = 1+self.A0*eps*x  # Hawmaker constant, linearly along y-direction
             dW = A*self.dwetting_pot(h)  # W'(h)
             dW = dW + self.dbulk_force(x,h)
             self.add_residual(weak(partial_t(h),v)+weak(Q*grad(p),grad(v)))
             self.add_residual(weak(p,eta)-weak(grad(h),grad(eta))-weak(dW,eta))

class ThinfilmProblem(Problem):
    def __init__(self, L=20, N=1000):
        super().__init__()
        self.L = L
        self.N = N
    def define_problem(self):
        mesh=RectangularQuadMesh(N=[self.N,self.N], size=[self.L,self.L],lower_left=[-self.L/2,-self.L/2])
        # Add the mesh (default name is "domain" with boundaries "left" and "right")
        self.add_mesh(mesh)
        x = var("coordinate")

        # Assemble the system
        equations = ThinfilmEquation()  # create a Poisson equation with source g=1
        equations += TextFileOutput()  # Add a simple text file output
        alpha = 0.5    # Fraction of system size occupied by the droplet
        r = np.sqrt(dot(x,x))
        initial_state = ha+np.exp(-r/(alpha*self.L/2))*ha*(2+np.cos(np.pi*2*(dot(x,x)**0.5)/(np.sqrt(2)*alpha*self.L)))
        #equations += InitialCondition(h=1+0.5*np.cos(np.pi*2*(dot(x,x)**0.5)/(np.sqrt(2)*self.L)))
        equations += InitialCondition(h=initial_state)
        equations += SpatialErrorEstimator(h=1.0)

        # Observables
        h = var('h')
        #hout = weak()
        hout = ha # TODO: Measure the outer!
        #hout = var('hout')

        hhat = h-hout
        whatin = wetting_pot(h)-wetting_pot(hout)
        equations += IntegralObservables(excess_mass = hhat)
        equations += IntegralObservables(Omega = whatin)
        equations += IntegralObservables(I = hout**3*hhat/(h**3))
        equations += IntegralObservables(K = hhat**2/h**3)
        x = var('coordinate_x')
        y = var('coordinate_y')
        equations += IntegralObservables(xCMm = x*hhat)        # Need to divide by the excess mass to find vCM
        equations += IntegralObservables(yCMm = y*hhat)
        if switch_bulk > 0 or switch_inhomo > 0:
                equations += IntegralObservableOutput(filename='obs_eps='+str(eps)+'bulk='+str(switch_bulk)+'inhomo='+str(switch_inhomo))
        else:
                equations += IntegralObservableOutput(filename='obs_stat')      # Stationary state
        self.add_equations(equations @ "domain")  # Add the equation system on the domain named "domain"
        
        self.add_equations(equations @ "domain")  # Add the equation system on the domain named "domain"
        


if __name__ == "__main__":
    with ThinfilmProblem(L=25,N=100) as problem:

        # Maximum refinement level
        problem.max_refinement_level = 1
        problem.max_permitted_error = 0.0005
        problem.min_permitted_error = 0.00005
        problem.run(endtime=1e4, startstep=1, outstep=True, temporal_error=1, spatial_adapt=problem.max_refinement_level)

        #problem.solve()  # Solve the problem
        print(problem.get_mesh("domain").evaluate_all_observables())
        problem.output()  # Write output
