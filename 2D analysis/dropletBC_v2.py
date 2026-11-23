from pyoomph import *
from pyoomph.expressions import *
from pyoomph.equations.additional import RefineAccordingToElement
import numpy as np

# ------------------------------------------------------------------ user settings
ha = 0.5            # precursor height
FINE = 3            # refinement level at the contact line (try 7 for a convergence check)
BASE_N = 100        # base mesh N x N -> element size L/BASE_N = 10 for EPS=1e-3
T_RELAX = 1e5       # relaxation time (c = 0)
T_DRIVE = T_RELAX+10      # e.g. switch on the drive afterwards, None = no drive
RAMP_TIME = 0.1 

c=0.1
A=1

try:
    h_center
except NameError:
    h_center = 6.0

try:
    eps
except NameError:
    eps = 0.01

def theta(x,x0,sigma=0.1):
    return 0.5*(1+np.tanh((x-x0)/sigma))
def wetting_pot(h):
        return (ha**3/5*h**(-5) - 1/2*h**(-2))

class ThinFilmEquation(Equations):
    def __init__(self, mu=1):
        super().__init__()
        self.mu = mu

    def define_fields(self):
        self.define_scalar_field("h", "C2")
        self.define_scalar_field("p", "C2")
    def define_residuals(self):
        h, v = var_and_test("h")
        p, eta = var_and_test("p")
        x = var("coordinate_x")
        t = var("time")
        ramp = theta(t,T_RELAX,RAMP_TIME)

        Q = h**3 / (3 * self.mu)                                        # mobility
        dW = A*(-ha**3 * h**(-6) + h**(-3)) + c * ramp * sin(2 * pi * eps * x)

        self.add_residual(weak(partial_t(h), v) + weak(Q * grad(p), grad(v)))
        self.add_residual(weak(p, eta) - weak(grad(h), grad(eta)) - weak(dW, eta))

def equilibrium_hout(ha, hmax, R):
     """The value is computed requiring the pressure to be equal at the maximum height of the drop and far from the drop;
        For the pressure at the maximum height, the disjoining pressure contribution is not neglected, while the disjoining pressure
        outside the drop is expanded to leading order in (Hout-1)^2"""

     return ha*(1+2*ha**3*hmax/(3*R**2))


class ThinFilmProblem(Problem):
    def __init__(self, base_N=BASE_N):
        super().__init__()
        self.h_center = h_center
        self.fine = FINE
        self.base_N = base_N
        self.L = 1 / eps
        # centre of the drop
        self.xc, self.yc = 0.0, 0.0

    def define_problem(self):
        L, hc = self.L, self.h_center
        R = 2*hc*ha/np.sqrt(3 / 5)*(1-ha/hc)**0.5
        F = self.fine

        self.add_mesh(RectangularQuadMesh(N=[self.base_N] * 2, size=[L, L],
                                          lower_left=[-L / 2, -L / 2]))

        eqs = ThinFilmEquation()
        eqs += TextFileOutput()

        # ---- geometric refinement: element size = (L/base_N) / 2**level
        graded = [(1.0, F), (2.5, F - 1), (5.0, F - 2), (10.0, F - 3), (20.0, F - 4)]

        # Recompute R. In the relaxation run it is useless as it would give the same value, but after the relaxation hc changes so it would give a better estimate during driving
        #R = 
        def refine_target(e):
            x, y = e.get_Eulerian_midpoint()
            d = np.sqrt((x - self.xc)**2 + (y - self.yc)**2) - R
            for dmax, lev in graded:
                if abs(d) < dmax:
                    return max(lev, 0)
            if d < 0:
                return 2          # smooth interior of the drop
            return 1 if d < 60 else 0

        eqs += RefineAccordingToElement(refine_target)
        eqs += TemporalErrorEstimator(h=1, p=1)    # needed, otherwise temporal_error= is ignored

        # ---- initial condition: parabolic cap on a film at the equilibrium height
        h_out = equilibrium_hout(ha,hc,R)
        r = var("coordinate")
        cap = hc * (1 - dot(r, r) / R**2)
        eqs += InitialCondition(h=maximum(cap, h_out))

        eqs += PeriodicBC("right", offset=[L, 0]) @ "left"
        eqs += PeriodicBC("top", offset=[0, L]) @ "bottom"

        eqs += IntegralObservables(V=var("h"))
        eqs += IntegralObservableOutput(filename="obs")

        self.add_equations(eqs @ "domain")


if __name__ == "__main__":
    with ThinFilmProblem(base_N=BASE_N) as problem:
        problem.set_output_directory(f"pbc_hc{problem.h_center:g}_eps{eps:g}_F{problem.fine}")
        problem.max_refinement_level = problem.fine

        # Stage 1: relaxation without drive
        problem.run(endtime=T_RELAX-1, startstep=0.01, maxstep=1e4,
                    temporal_error=1, outstep=1e4, spatial_adapt=problem.fine)
        
        #problem.initialise()
        #problem.load_state("pbc_hc5_eps0.01_F3/_states/state_000006.dump", ignore_outstep=True)

        problem.fine = 6
        problem.max_refinement_level = problem.fine
        # Stage 2 (optional): switch the bulk force on and continue
        if T_DRIVE is not None:
            problem.run(endtime=T_DRIVE, startstep=0.01, maxstep=0.1,
                        temporal_error=1e-1, outstep=0.1, spatial_adapt=problem.fine)


        print(problem.get_mesh("domain").evaluate_all_observables())
        problem.output()
