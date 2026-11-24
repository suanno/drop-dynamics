from pyoomph import *
from pyoomph.expressions import *
import numpy as np
import sys
from scipy.optimize import root_scalar
import os
import time


ha = 0.5            # precursor height
T_RELAX = 3e6               # relaxation time (c = 0)

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
def dwetting_pot(h,ha):
        return (ha**3/5*h**(-5) - 1/2*h**(-2))
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

             Q=(h**3)/(3*self.mu) # Mobility
             dW = A*dwetting_pot(h,ha)  # No driving

             self.add_residual(weak(partial_t(h),v)+weak(Q*grad(p),grad(v)))
             self.add_residual(weak(p,eta)-weak(grad(h),grad(eta))-weak(dW,eta))

def equilibrium_hout(ha, hmax, R):
    """The value is computed requiring the pressure to be equal at the maximum height of the drop and far from the drop;
        We consider the smaller root of dhWout=Pin, since we know that Hout decreases with the size of the drop in the continuation branch"""

    return ha*(1+1.1*2*ha**3*hmax/(3*R**2))

    f = lambda x: domega(x)-2*hmax*ha**3/R**2
    result = root_scalar(
        f,
        bracket=[1, 2**(1/3)],
        method='brentq'
    )
    return result.root*ha

class ThinfilmProblem(Problem):
    def __init__(self, L=20, N=1000):
        super().__init__()
        self.h_center = h_center
        self.L = L
        self.N = N
        ####Continuation
        self.xCOM = self.define_global_parameter(xCOM=0) 
        self.V1 = self.define_global_parameter(V1=1e4)
        #self.vx   = self.get_global_parameter("vx")

    def define_problem(self):
        L, hc = self.L, self.h_center
        h_out=ha    # First guess to compute R, then it will be computed later from the R and hc
        tanTheta=np.sqrt(-2*omega(h_out/ha))
        R = 2*hc*ha/tanTheta*(1-ha/hc)**0.5
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
            
            tanTheta=np.sqrt(-2*omega(h_out/ha))
            R=2*hc*ha/tanTheta*(1-ha/hc)**0.5
        """
        h_out = equilibrium_hout(ha,hc,R)

        cap = hc * (1 - x**2 / R**2)  # Height functions of the droplets
        h_init = maximum(cap, h_out)  # Initial height: maximum of h1, h2 and precursor
        equations += InitialCondition(h=h_init)
        equations += SpatialErrorEstimator(h=1.0)
        equations += PeriodicBC("right", offset=self.L) @ "left"
        equations += IntegralObservables(V=var("h"))  # Total volume


        # Continuation
        equations += IntegralObservables(V1_integral=var("h")-ha)
        equations += IntegralObservables(liq1_com=var("coordinate_x")*(var("h")-ha)/self.get_global_parameter("V1"))

        equations += IntegralObservables(_area=self.L,_h_integral=var("h")-ha,_h_sqr_integral=(var("h")-ha)**2)
        equations += IntegralObservables(l2norm=lambda _area,_h_sqr_integral : 1/_area*square_root(_h_sqr_integral))
        #V1_lagr,V2_lagr haave to be adaptive on runtime: global variables
        V1_lagr, V1_lagr_test = self.add_global_dof("V1_lagr", equation_contribution=-self.V1, domain="vol1_constraint", only_for_stationary_solve=True)
        vx_lagr, vx_lagr_test = self.add_global_dof("vx_lagr", equation_contribution=-self.xCOM, domain="vx_constraint", only_for_stationary_solve=True)
        #we add the integral of the heights to the test space of the global Lagrange multipliers
        equations += WeakContribution(var("h")-ha, V1_lagr_test,dimensional_dx=False)
        equations += WeakContribution(var("coordinate_x")*(var("h")-ha)/self.get_global_parameter("V1"),vx_lagr_test, dimensional_dx=False)
        #link back to actual equations as additional pressures within the two drops
        equations += WeakContribution(-V1_lagr+vx_lagr*dot(grad(var("coordinate_x")),grad(var("h"))), testfunction("h"))

        # Save
        equations += IntegralObservableOutput(filename=f'obs')      # Stationary state
        self.add_equations(equations @ "domain")  # Add the equation system on the domain named "domain"

if __name__ == "__main__":
    with ThinfilmProblem(L=eps**(-1),N=int(0.5*1e3/15*eps**(-1))) as problem:
                #problem.set_output_directory(f"pbc_hcenter_{h_center:g}_c={amplitude_C:g}_eps={eps:g}")
                problem.set_output_directory(f"pbc_hc{problem.h_center:g}_cont")
        
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
                    endtime=T_RELAX ,
                    startstep=0.01,
                    maxstep=1e5,
                    temporal_error=1,
                    outstep=True,
                    spatial_adapt=1
                ) 
                problem.quiet()
                print(problem.get_mesh("domain").evaluate_all_observables())
                problem.output()  # Write output


                # Continuation
                problem.xCOM.value = problem.get_mesh("domain").evaluate_all_observables()["liq1_com"] 
                problem.V1.value = problem.get_mesh("domain").evaluate_all_observables()["V1_integral"] 
                print("\n\nWe have set the center-of-mass for liquid 1 to ",problem.xCOM.value)

                time.sleep(1)

                n_eig_vct = 10
                cont_file=open(os.path.join(problem.get_output_directory(),"cont.txt"),"w")
                def write_to_cont_file():
                    # this one somehow does not work with slepc
                    #eigenvals,eigenvects=problem.solve_eigenproblem(n=n_eig_vct)#number of eigenvalues and igenvectors to compute.
                    vx_num=problem.get_ode("vx_constraint").get_value("vx_lagr")
                    l2norm = problem.get_mesh("domain").evaluate_all_observables()["l2norm"]
                    print("V1,vx,l2norm = ", problem.V1.value,vx_num,l2norm)
                    #h_rms=problem.get_mesh("domain").evaluate_observable("h_rms") # get the root mean square
                    line=[problem.V1.value,vx_num,l2norm]#,h_rms,problem.V1.value,problem.V2.value]
                    cont_file.write("\t".join(map(str,line))+"\n")
                    cont_file.flush()


                write_to_cont_file() # write the first state
                ds=1e-4
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


