from pyoomph import *
from pyoomph.expressions import *
import numpy as np
import sys
from scipy.optimize import root_scalar

ha = 0.5            # precursor height
T_RELAX = 3e6               # relaxation time (c = 0)
T_DRIVE = T_RELAX+10        # e.g. switch on the drive afterwards, None = no drive
RAMP_TIME = 0.1 

A=1
c=1

try:
    h_center
except NameError:
    h_center = 3.0

try:
    eps
except NameError:
    eps = 0.001

def theta(x,x0,sigma=0.1):
      return 0.5*(1+np.tanh((x-x0)/sigma))
def wetting_pot(h):
        return (ha**3/5*h**(-5) - 1/2*h**(-2))
def dwetting_pot(h):
        return -ha**3*h**(-6)+h**(-3)
def domega(H):
        return (-H**(-6) + H**(-3))
def omega(H):
        return (1/5*H**(-5) - 1/2*H**(-2))
class ThinfilmEquation(Equations):
     def __init__(self,mu=1):
        super(ThinfilmEquation, self).__init__()
        self.mu=mu # Viscosity (in the mobility Q)
 
     def define_fields(self):
             self.define_scalar_field("h","C2")
             self.define_scalar_field("p","C2")

     def define_residuals(self):
             h,v=var_and_test("h")
             p,eta=var_and_test("p")
             x = var("coordinate_x")
             t = var("time")
             ramp = theta(t,T_RELAX,RAMP_TIME)

             Q=(h**3)/(3*self.mu) # Mobility
             dW = A*dwetting_pot(h) + c*ramp*sin(2*pi*eps*x)  # W'(h)

             self.add_residual(weak(partial_t(h),v)+weak(Q*grad(p),grad(v)))
             self.add_residual(weak(p,eta)-weak(grad(h),grad(eta))-weak(dW,eta))


def equilibrium_hout(ha, hmax, R):
    """The value is computed requiring the pressure to be equal at the maximum height of the drop and far from the drop;
        We consider the smaller root of dhWout=Pin, since we know that Hout decreases with the size of the drop in the continuation branch"""

    #return ha*(1+2*ha**3*hmax/(3*R**2))

    f = lambda x: domega(x)-2*hmax*ha**3/R**2
    result = root_scalar(
        f,
        bracket=[1, 2**(1/3)],
        method='brentq'
    )
    return result.root*ha

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

class ThinfilmProblem(Problem):
    def __init__(self, L=20, N=1000):
        super().__init__()
        self.h_center = h_center
        self.L = L
        self.N = N

    def define_problem(self):
        L, hc = self.L, self.h_center
        h_out=ha    # First guess to compute R, then it will be computed later from the R and hc
        tanTheta=np.sqrt(-2*omega(h_out/ha))
        #R = 2*hc*ha/tanTheta*(1-ha/hc)**0.5
        R = 2*hc*ha/np.sqrt(3 / 5)*(1-ha/hc)**0.5
        mesh = LineMesh(minimum=-self.L/2, size=self.L, N=self.N)
        # Add the mesh (default name is "domain" with boundaries "left" and "right")
        self.add_mesh(mesh)
        x = var("coordinate_x")
        
        # Assemble the system
        equations = ThinfilmEquation()  # create a Poisson equation with source g=1
        equations += TextFileOutput()  # Add a simple text file output
        # Since hout depends on R and viceversa, we repeat some times up to convergence
        """
        for i in np.arange(0,100):
            h_out = equilibrium_hout(ha,hc,R)
            tanTheta=np.sqrt(-2*omega(h_out/ha))
            R=2*hc*ha/tanTheta*(1-ha/hc)**0.5
        """
        h_out = equilibrium_hout(ha,hc,R)
        #h_out = equilibrium_film_height(ha, 2*hc/R**2)
        cap = hc * (1 - x**2 / R**2)  # Height functions of the droplets
        h_init = maximum(cap, h_out)  # Initial height: maximum of h1, h2 and precursor
        equations += InitialCondition(h=h_init)
        equations += SpatialErrorEstimator(h=1.0)
        equations += PeriodicBC("right", offset=self.L) @ "left"
        equations += IntegralObservables(V=var("h"))  # Total volume
        equations += IntegralObservableOutput(filename=f'obs')      # Stationary state
        self.add_equations(equations @ "domain")  # Add the equation system on the domain named "domain"
        # The maximum height and its x-coordinate (droplet position) are measured from domain*.txt


if __name__ == "__main__":
    with ThinfilmProblem(L=eps**(-1),N=int(0.5*1e3/15*eps**(-1))) as problem:
                #problem.set_output_directory(f"pbc_hcenter_{h_center:g}_c={amplitude_C:g}_eps={eps:g}")
                problem.set_output_directory(f"pbc_hc{problem.h_center:g}_eps{eps:g}_c={c:g}")
        
                # Maximum refinement level
                """
                problem.max_refinement_level = 1
                problem.max_permitted_error = 0.0005
                problem.min_permitted_error = 0.00005
                problem.run(endtime=5e5, startstep=0.01, maxstep=1e4, outstep=True, temporal_error=1, spatial_adapt=problem.max_refinement_level)
                """

                problem.max_refinement_level = 1
                problem.max_permitted_error = 5e-4
                problem.min_permitted_error = 5e-5
                problem.run(
                    endtime=T_RELAX-1,
                    startstep=0.01,
                    maxstep=1e5,
                    temporal_error=1,
                    outstep=True,
                    spatial_adapt=1
                ) 
                problem.quiet()
                problem.run(
                    endtime=T_DRIVE,
                    startstep=0.01,
                    maxstep=0.1,
                    temporal_error=1e-1,
                    outstep=True,
                    spatial_adapt=5
                ) 
                #problem.run(endtime=20, timestep=1e-2, outstep=True, maxstep=1e4, temporal_error=1e-2,spatial_adapt=15,newton_solver_tolerance=1e-10)
                problem.quiet()
                #problem.solve()  # Solve the problem
                print(problem.get_mesh("domain").evaluate_all_observables())
                problem.output()  # Write output
