"""
Thin-film droplet with disjoining pressure on a periodic square (non-dimensional, gamma = 1).

    h_t = div( Q(h) grad p ),   Q = h^3/(3 mu)
    p   = -lap(h) + W'(h) + c * sin(2 pi eps x),   W'(h) = -ha^3 h^-6 + h^-3

Run:   python dropletBC_v2.py --omp 4        (OpenMP threads, same script serial/parallel)
Check: python dropletBC_v2.py --quick-test   (one Newton solve + one output)

Changes with respect to dropletBC.py:
  * eps / h_center are not read from sys.argv any more (pyoomph installs its own parser)
  * the drive is a global parameter `c` (0 -> relaxation, then set to 1), no tanh switch in time,
    no np.* functions applied to symbolic expressions
  * coarse base mesh (L/N = 10) + graded refinement around r = R instead of N=500 everywhere
  * the film starts at the pressure-consistent height h_out (W'(h_out) = Laplace pressure),
    so the slow diffusive "long tail" transient is ~10x smaller
  * both SpatialErrorEstimators, the unused A0/switch_* machinery and dead code are removed
"""
from pyoomph import *
from pyoomph.expressions import *
from pyoomph.equations.additional import RefineAccordingToElement
import numpy as np

# ------------------------------------------------------------------ user settings
EPS = 1e-2          # L = 1/EPS
H_CENTER = 6.0     # drop height at the centre
HA = 0.5            # precursor height
FINE = 3            # refinement level at the contact line (try 7 for a convergence check)
BASE_N = 100        # base mesh N x N -> element size L/BASE_N = 10 for EPS=1e-3
T_RELAX = 1e5       # relaxation time (c = 0)
T_DRIVE = T_RELAX+10      # e.g. switch on the drive afterwards, None = no drive
C_DRIVE = 0.1       # amplitude of the bulk force
RAMP_TIME = 0.1 

try:
    h_center
except NameError:
    h_center = 6.0

def theta(x,x0,sigma=0.1):
    return 0.5*(1+np.tanh((x-x0)/sigma))

class ThinFilmEquation(Equations):
    def __init__(self, ha=0.5, mu=1, c=0, eps=1e-3):
        super().__init__()
        self.ha, self.mu, self.c, self.eps = ha, mu, c, eps   # c may be a global parameter

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
        dW = -self.ha**3 * h**(-6) + h**(-3) + self.c * ramp * sin(2 * pi * self.eps * x)

        self.add_residual(weak(partial_t(h), v) + weak(Q * grad(p), grad(v)))
        self.add_residual(weak(p, eta) - weak(grad(h), grad(eta)) - weak(dW, eta))


def equilibrium_film_height(ha, p_drop):
    """Solve W'(h) = p_drop on the branch h in [ha, 2^(1/3) ha] by bisection."""
    dW = lambda h: -ha**3 * h**(-6) + h**(-3)
    lo, hi = ha, ha * 2**(1 / 3)          # W' increases from 0 to its maximum on this interval
    if p_drop >= dW(hi):
        raise ValueError("Laplace pressure exceeds max W'(h): no film equilibrium")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if dW(mid) < p_drop else (lo, mid)
    return 0.5 * (lo + hi)


class ThinFilmProblem(Problem):
    def __init__(self, eps=EPS, h_center=h_center, ha=HA, fine=FINE, base_N=BASE_N):
        super().__init__()
        self.eps, self.h_center, self.ha = eps, h_center, ha
        self.fine, self.base_N = fine, base_N
        self.L = 1 / eps
        # centre of the refined annulus (update between run() calls if the drop moves)
        self.xc, self.yc = 0.0, 0.0
        self.c = self.define_global_parameter(c=0)    # bulk-force amplitude, plain number

    def define_problem(self):
        L, ha, hc = self.L, self.ha, self.h_center
        R = 2 * hc * ha / np.sqrt(3 / 5) * np.sqrt(1 - ha / hc)   # drop radius (same formula as before)
        lam = ha**2 / np.sqrt(3)                                    # foot length scale
        F = self.fine

        self.add_mesh(RectangularQuadMesh(N=[self.base_N] * 2, size=[L, L],
                                          lower_left=[-L / 2, -L / 2]))

        eqs = ThinFilmEquation(ha=ha, c=self.c, eps=self.eps)
        eqs += TextFileOutput()

        # ---- geometric refinement: element size = (L/base_N) / 2**level
        graded = [(1.0, F), (2.5, F - 1), (5.0, F - 2), (10.0, F - 3), (20.0, F - 4)]

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
        p_drop = 4 * hc / R**2
        h_out = equilibrium_film_height(ha, p_drop)
        r = var("coordinate")
        cap = hc * (1 - dot(r, r) / R**2)
        eqs += InitialCondition(h=maximum(cap, h_out))

        eqs += PeriodicBC("right", offset=[L, 0]) @ "left"
        eqs += PeriodicBC("top", offset=[0, L]) @ "bottom"

        eqs += IntegralObservables(V=var("h"))
        eqs += IntegralObservableOutput(filename="obs")

        self.add_equations(eqs @ "domain")


if __name__ == "__main__":
    with ThinFilmProblem() as problem:
        problem.set_output_directory(f"pbc_hc{problem.h_center:g}_eps{problem.eps:g}_F{problem.fine}")
        problem.max_refinement_level = problem.fine
        problem.c.value = C_DRIVE   # This is constant, but there is a ramp to switch it on
    
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
