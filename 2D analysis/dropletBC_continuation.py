from pyoomph import *
from pyoomph.expressions import *
from pyoomph.equations.additional import RefineAccordingToElement
from scipy.optimize import root_scalar
import numpy as np
import os
import time

# ------------------------------------------------------------------ user settings
ha = 0.5            # precursor height
FINE = 5            # refinement level at the contact line (try 7 for a convergence check)
BASE_N = 100        # base mesh N x N -> element size L/BASE_N = 10 for EPS=1e-3
T_RELAX = 1e5       # relaxation time (c = 0)
#T_DRIVE = T_RELAX+10      # e.g. switch on the drive afterwards, None = no drive
#RAMP_TIME = 0.1 

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
def dwetting_pot(h,ha):
        return (-ha**3 * h**(-6) + h**(-3))
def domega(H):
        return (-H**(-6) + H**(-3))
def omega(H):
        return (1/5*H**(-5) - 1/2*H**(-2))

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
        dW = A*dwetting_pot(h,ha)

        self.add_residual(weak(partial_t(h), v) + weak(Q * grad(p), grad(v)))
        self.add_residual(weak(p, eta) - weak(grad(h), grad(eta)) - weak(dW, eta))


def equilibrium_hout(ha, hmax, R):
    """The value is computed requiring the pressure to be equal at the maximum height of the drop and far from the drop;"""

    #return ha*(1+2*ha**3*hmax/(3*R**2))

    f = lambda x: domega(x)-4*hmax*ha**3/R**2
    result = root_scalar(
        f,
        bracket=[1, 2**(1/3)],
        method='brentq'
    )
    return result.root*ha


class ThinFilmProblem(Problem):
    def __init__(self, base_N=BASE_N):
        super().__init__()
        self.h_center = h_center
        self.fine = FINE
        self.base_N = base_N
        self.L = 1 / eps
        # centre of the drop
        self.xc, self.yc = 0.0, 0.0
        ####Continuation
        self.xCOM = self.define_global_parameter(xCOM=0) 
        self.yCOM = self.define_global_parameter(yCOM=0) 
        self.V1 = self.define_global_parameter(V1=1e4)
        #self.vx   = self.get_global_parameter("vx")

    def define_problem(self):
        L, hc = self.L, self.h_center
        h_out=ha    # First guess to compute R, then it will be computed later from the R and hc
        tanTheta=np.sqrt(-2*omega(h_out/ha))
        R = 2*hc*ha/tanTheta*(1-ha/hc)**0.5
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


        # ---- initial condition: parabolic cap on a film at the equilibrium height
        # Since hout depends on R and viceversa, we repeat some times up to convergence
        for i in np.arange(0,100):
            h_out = equilibrium_hout(ha,hc,R)
            tanTheta=np.sqrt(-2*omega(h_out/ha))
            R=2*hc*ha/tanTheta*(1-ha/hc)**0.5
        
        r = var("coordinate")
        cap = hc * (1 - dot(r, r) / R**2)
        eqs += InitialCondition(h=maximum(cap, h_out))

        eqs += PeriodicBC("right", offset=[L, 0]) @ "left"
        eqs += PeriodicBC("top", offset=[0, L]) @ "bottom"

        #eqs += IntegralObservables(V=var("h"))
        eqs += IntegralObservableOutput(filename="obs")

        #eqs += RefineAccordingToElement(refine_target)
        eqs += SpatialErrorEstimator(h=10,p=1)
        eqs += TemporalErrorEstimator(h=10, p=1)    # needed, otherwise temporal_error= is ignored



        # Continuation
        eqs += IntegralObservables(V1_integral=var("h")-ha)
        eqs += IntegralObservables(liq1_comx=var("coordinate_x")*(var("h")-ha)/self.get_global_parameter("V1"))
        eqs += IntegralObservables(liq1_comy=var("coordinate_y")*(var("h")-ha)/self.get_global_parameter("V1"))
        #eqs += IntegralObservables(h_edge=var("h"), _len=1) @ "outer_boundary"
        h = var("h")

        eqs += IntegralObservables(
            # Far-field height

            # Helper integrals
            h_int=h,
            _area=1,
            invh=1/h,
            invh2=1/h**2,
            invh3=1/h**3,

            _extra_area = theta(dot(r,r)**0.5,L/2,sigma=0.00001),
            _hout=var("h") * theta(dot(r,r)**0.5,L/2,sigma=0.00001),
            hout = lambda _hout, _extra_area: _hout/_extra_area,
            #V=lambda _h_int,hout,_area: _h_int-hout*_area,
            #V=lambda _h_int,_area,hout:
            #    _h_int-hout*_area,

            #IB0=lambda _inv_h2,_inv_h3,hout:
            #    hout**3*(_inv_h2-hout*_inv_h3),

            #KB0=lambda _inv_h,_inv_h2,_inv_h3,hout:
            #    _inv_h-2*hout*_inv_h2+hout**2*_inv_h3
        )
        

        eqs += IntegralObservables(_area=self.L**2,_h_integral=var("h")-ha,_h_sqr_integral=(var("h")-ha)**2)
        eqs += IntegralObservables(l2norm=lambda _area,_h_sqr_integral : 1/_area*square_root(_h_sqr_integral))
        #V1_lagr,V2_lagr haave to be adaptive on runtime: global variables
        V1_lagr, V1_lagr_test = self.add_global_dof("V1_lagr", equation_contribution=-self.V1, domain="vol1_constraint", only_for_stationary_solve=True)
        vx_lagr, vx_lagr_test = self.add_global_dof("vx_lagr", equation_contribution=-self.xCOM, domain="vx_constraint", only_for_stationary_solve=True)
        vy_lagr, vy_lagr_test = self.add_global_dof("vy_lagr", equation_contribution=-self.yCOM, domain="vy_constraint", only_for_stationary_solve=True)
        #we add the integral of the heights to the test space of the global Lagrange multipliers
        eqs += WeakContribution(var("h")-ha, V1_lagr_test,dimensional_dx=False)
        eqs += WeakContribution(var("coordinate_x")*(var("h")-ha)/self.get_global_parameter("V1"),vx_lagr_test, dimensional_dx=False)
        eqs += WeakContribution(var("coordinate_y")*(var("h")-ha)/self.get_global_parameter("V1"),vy_lagr_test, dimensional_dx=False)
        #link back to actual equations as additional pressures within the two drops
        eqs += WeakContribution(-V1_lagr+vx_lagr*dot(grad(var("coordinate_x")),grad(var("h")))+vy_lagr*dot(grad(var("coordinate_y")),grad(var("h"))), testfunction("h"))


        self.add_equations(eqs @ "domain")


if __name__ == "__main__":
    with ThinFilmProblem(base_N=BASE_N) as problem:
        problem.set_output_directory(f"pbc_hc{problem.h_center:g}_cont")
        problem.max_refinement_level = problem.fine

        # Stage 1: relaxation without drive
        problem.run(endtime=T_RELAX-1, startstep=0.01, maxstep=1e4, outstep=1e4, temporal_error=1)
        
        #problem.initialise()
        #problem.load_state("pbc_hc5_eps0.01_F3/_states/state_000006.dump", ignore_outstep=True)

        # Continuation
        problem.xCOM.value = problem.get_mesh("domain").evaluate_all_observables()["liq1_comx"] 
        problem.yCOM.value = problem.get_mesh("domain").evaluate_all_observables()["liq1_comy"] 
        problem.V1.value = problem.get_mesh("domain").evaluate_all_observables()["V1_integral"] 
        #print("\n\nWe have set the center-of-mass for liquid 1 to ",problem.xCOM.value)

        time.sleep(1)

        n_eig_vct = 10
        cont_file=open(os.path.join(problem.get_output_directory(),"cont.txt"),"w")
        def write_to_cont_file():
            # this one somehow does not work with slepc
            #eigenvals,eigenvects=problem.solve_eigenproblem(n=n_eig_vct)#number of eigenvalues and igenvectors to compute.
            vx_num=problem.get_ode("vx_constraint").get_value("vx_lagr")
            vy_num=problem.get_ode("vy_constraint").get_value("vy_lagr")
            l2norm = problem.get_mesh("domain").evaluate_all_observables()["l2norm"]
            print("V1,vx,l2norm = ", problem.V1.value,vx_num,vy_num,l2norm)
            #h_rms=problem.get_mesh("domain").evaluate_observable("h_rms") # get the root mean square
            line=[problem.V1.value,vx_num,vy_num,l2norm]#,h_rms,problem.V1.value,problem.V2.value]
            cont_file.write("\t".join(map(str,line))+"\n")
            cont_file.flush()


        write_to_cont_file() # write the first state
        ds=1e-4
        stop_counter=0
        while (stop_counter<5):
            stop_counter +=1
            ds=problem.arclength_continuation(problem.V1,
                                            step=ds,
                                            spatial_adapt=1,
                                            max_ds=1e+4,
                                            newton_solver_tolerance=1e-10,
                                            min_ds=1e-13
                                            )
            problem.output()
            write_to_cont_file() # write the updated state
        ds=-1e-2
        stop_counter=0
        while (stop_counter<500):
            stop_counter +=1
            ds=problem.arclength_continuation(problem.V1,
                                            step=ds,
                                            spatial_adapt=1,
                                            max_ds=1e+4,
                                            newton_solver_tolerance=1e-10,
                                            min_ds=1e-13
                                            )
            problem.output()
            write_to_cont_file() # write the updated state
        """
        problem.fine = 6
        problem.max_refinement_level = problem.fine
        # Stage 2 (optional): switch the bulk force on and continue
        if T_DRIVE is not None:
            problem.run(endtime=T_DRIVE, startstep=0.01, maxstep=0.1,
                        temporal_error=1e-1, outstep=0.1, spatial_adapt=problem.fine)

        """
        print(problem.get_mesh("domain").evaluate_all_observables())
        problem.output()
