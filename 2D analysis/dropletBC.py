from pyoomph import *
from pyoomph.equations.additional import RefineAccordingToElement
from pyoomph.expressions import *
import numpy as np
import sys

# Switch ON/OFF the driving forces
amplitude_C = 0
amplitude_A = 0
switch_bulk = 0
switch_inhomo = 0
t0 = 1e6+1

# Parameters (additional parameters are in the initial condition or among the cmd inputs)
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
        self.R = (2*h_center/(np.sqrt(3/5)*ha**(-1)))*(1-ha/h_center)**(1/2)
        

    def define_problem(self):        
        # Add the mesh (default name is "domain" with boundaries "left" and "right")
        self.add_mesh(RectangularQuadMesh(N=[self.N,self.N],size=[self.L,self.L],lower_left=[-self.L*0.5,-self.L*0.5])) 
		
        # Assemble the system
        equations = ThinfilmEquation()
        equations += TextFileOutput()

        # Global spatial error estimator
        equations += SpatialErrorEstimator(
            h=1,
            group="global"
        )
        # Additional absolute-error estimator for the tail
        equations += SpatialErrorEstimator(
            h=1,
            group="tail",
            normalize_relative=0,
            weight=1
        )
        # Spatial error estimator inside the droplet and in the region close to the contact line
        lam = ha**2 / np.sqrt(3)
        R = self.R
        def refine_tail(e):
            x, y = e.get_Eulerian_midpoint()
            r = np.sqrt(x**2 + y**2)

            if r < R + 6*lam:
                return 3
            else:
                return 0

        equations += RefineAccordingToElement(
            refine_tail,
            prevent_unrefinement=True
        )
        equations += TemporalErrorEstimator(h=1, p=1)

        # Initial state
        r = var("coordinate")
        h = var("h")
        h_init = self.h_center * (1 - dot(r, r) / self.R**2)  # Height functions of the droplets
        deltaHout = h_center**(-1.19562246)*np.exp(-1.11011093)
        deltaHout=0
        hout0 = ha*(1+deltaHout)
        h_init = maximum(h_init, hout0)  # Initial height: maximum of h1, h2 and precursor
        equations += InitialCondition(h=h_init)

        equations += PeriodicBC("right", offset=[self.L, 0]) @ "left"
        equations += PeriodicBC("top", offset=[0, self.L]) @ "bottom"

        equations += IntegralObservables(V=var("h"))  # Total volume
        equations += IntegralObservableOutput(filename=f'obs')      # Stationary state
        self.add_equations(equations @ "domain")


if __name__ == "__main__":
    with ThinfilmProblem(L=eps**(-1), N=500) as problem:

        # Maximum refinement level
        #problem.max_refinement_level = 6
        #problem.max_permitted_error = 0.0005
        #problem.min_permitted_error = 0.00005

        problem.set_output_directory(f"pbc_hcenter_{h_center:g}_c={amplitude_C:g}_eps={eps:g}")
        # Maximum refinement level
        #problem.max_refinement_level = 1
        #problem.max_permitted_error = 0.0005
        #problem.min_permitted_error = 0.00005
        #problem.run(endtime=5e5, startstep=0.01, maxstep=1e4, outstep=True, temporal_error=1, spatial_adapt=problem.max_refinement_level)
               
        problem.max_refinement_level = 4
        problem.max_permitted_error = 5e-4
        problem.min_permitted_error = 5e-5

        problem.run(
            endtime=t0-1,
            startstep=0.01,
            maxstep=1e4,
            temporal_error=1e-1,
            outstep=True,
            spatial_adapt=3
        ) 
        """
        # Activate driving force
        problem.run(
            endtime=t0+10,
            startstep=0.01,
            maxstep=0.1,
            temporal_error=1,
            outstep=True,
            spatial_adapt=5
        ) 
        """
        """         
problem.max_refinement_level = 0

        problem.run(
            endtime=1,
            startstep=0.01,
            maxstep=0.01,
            outstep=True,
            spatial_adapt=0
        )
        """
        #problem.run(endtime=5e5, timestep=1e-2, outstep=True, maxstep=1e4, temporal_error=1e-2,spatial_adapt=15,newton_solver_tolerance=1e-10)
        #problem.quiet()
        # Activate bulk force
        #problem.run(endtime=1e5+11, timestep=1e-2, outstep=True, maxstep=0.1, temporal_error=1e-2,spatial_adapt=15,newton_solver_tolerance=1e-10)
        #problem.quiet()

        # problem.solve()  # Solve the problem
        print(problem.get_mesh("domain").evaluate_all_observables())
        problem.output()  # Write output