from pyoomph import *
from pyoomph.expressions import *
import numpy as np
import sys

# Switch ON/OFF the driving forces
amplitude_C = 1
amplitude_A = 1
switch_bulk = 1
switch_inhomo = 0

# Parameters (additional parameters are in the initial condition)
ha = 0.5
eps = 0.001

try:
    h_center
except NameError:
    h_center = 1.0

if len(sys.argv) > 1:  # Set epsilon from cmd line
    eps = float(sys.argv[1])


def wetting_pot(h):
    return ha**3 / 5 * h**(-5) - 1 / 2 * h**(-2)


def theta(x, x0, sigma=0.1):
    return 0.5 * (1 + np.tanh((x - x0) / sigma))


class ThinfilmEquation(Equations):
    def __init__(self, mu=1):
        super(ThinfilmEquation, self).__init__()
        self.mu = mu  # Viscosity (in the mobility Q)


        t = var("time")

        # Activate bulk force and inhomogeneous Hawmaker constant for t > t0
        t0 = 6
        activation_fun = theta(t, t0)

        # Select whether activate c, A or both drivings
        amplitude_C = 1
        amplitude_A = 1

        if switch_bulk > 0:
            self.c = activation_fun * amplitude_C
        else:
            self.c = 0

        if switch_inhomo > 0:
            self.A0 = activation_fun * amplitude_A
        else:
            self.A0 = 0

    def wetting_pot(self, h):
        return ha**3 / 5 * h**(-5) - 1 / 2 * h**(-2)

    def dwetting_pot(self, h):
        return -ha**3 * h**(-6) + h**(-3)

    # Bulk force on x axis
    def bulk_force(self, x, h):
        return -eps * self.c * x * h

    #def dbulk_force(self, x, h):
    #    return -eps * self.c * x
    
    def dbulk_force(self,x,h):
            return self.c*np.sin(2*np.pi*x*eps)

    def define_fields(self):
        self.define_scalar_field("h", "C2")
        self.define_scalar_field("p", "C2")

    def define_residuals(self):
        h, v = var_and_test("h")
        p, eta = var_and_test("p")
        x = var("coordinate_x")

        Q = (h**3) / (3 * self.mu)  # Mobility
        A = 1 + self.A0 * eps * x  # Hawmaker constant, linearly along y-direction
        dW = A * self.dwetting_pot(h) + self.dbulk_force(x, h)  # W'(h)

        self.add_residual(
            weak(partial_t(h), v) + weak(Q * grad(p), grad(v))
        )
        self.add_residual(
            weak(p, eta) - weak(grad(h), grad(eta)) - weak(dW, eta)
        )


class ThinfilmProblem(Problem):
    def __init__(self, L=20, N=1000):
        super().__init__()
        self.L = L
        self.N = N
        self.h_center = h_center
        self.R = 1

    def define_problem(self):        
        # Add the mesh (default name is "domain" with boundaries "left" and "right")
        self.add_mesh(RectangularQuadMesh(N=[self.N,self.N],size=[self.L,self.L],lower_left=[-self.L*0.5,-self.L*0.5])) 
		


        h = var("h")  # Building disjoining pressure
        x = var("coordinate_x")
        r = var("coordinate")

        # Assemble the system
        equations = ThinfilmEquation()
        equations += TextFileOutput()

        # equations += InitialCondition(h=1+0.5*np.cos(np.pi*2*x/self.L))

        h_init = self.h_center * (1 - x**2 / self.R**2)  # Height functions of the droplets
        h_init = maximum(h_init, ha)-(amplitude_C/amplitude_A)*(ha**4/3)*np.sin(2*np.pi*x*eps)      # Parabolic cup with relaxed outer film
        
        equations += InitialCondition(h=h_init)

        equations += SpatialErrorEstimator(h=1.0)
        equations += PeriodicBC("right", offset=[self.L, 0]) @ "left"
        equations += PeriodicBC("top", offset=[0, self.L]) @ "bottom"

        self.add_equations(equations @ "domain")

        # The maximum height and its x-coordinate (droplet position)
        # are measured from domain*.txt
        #equations += IntegralObservables(gradh = (eps**(-1)*grad(h)*boundary_function_left+eps**(-1)*grad(h)*boundary_function_right)/(self.L*(1-toll)))  # Measuring the chi derivative
        equations += IntegralObservableOutput(filename=f'obs')      # Stationary state


if __name__ == "__main__":
    with ThinfilmProblem(L=120, N=80) as problem:

        # Maximum refinement level
        problem.max_refinement_level = 6
        #problem.max_permitted_error = 0.0005
        #problem.min_permitted_error = 0.00005
        problem.set_output_directory(f"pbc_hcenter_{h_center:g}_c={amplitude_C:g}_eps={eps:g}_relaxed")
        #problem.set_output_directory(f"single_droplet_hcenter_{h_center:g}_c=1")
        problem.run(
            endtime=10,
            startstep=0.01,
            maxstep=0.1,
            outstep=True,
            temporal_error=1,
            spatial_adapt=problem.max_refinement_level,
        )

        # problem.solve()  # Solve the problem
        print(problem.get_mesh("domain").evaluate_all_observables())
        problem.output()  # Write output