from pyoomph import *
from pyoomph.expressions import *
import numpy as np
import sys
import os
import time

# Switch ON/OFF the driving forces
amplitude_C = 0
amplitude_A = 0
switch_bulk = 0
switch_inhomo = 0
# Parameters (additional parameters are in the initial condition)
ha = 0.5
try:
    h_center
except NameError:
    h_center = 100.0

try:
    eps
except NameError:
    eps = 1e-3


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
             t0=5
             activation_fun = theta(t,t0)
             # Select wether activate c, A or both drivings
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
             return self.c*np.sin(2*np.pi*x*eps)*h
     def dbulk_force(self,x,h):
             return self.c*np.sin(2*np.pi*x*eps)

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
    def __init__(self, L=20, N=8):
        super().__init__()
        self.L = L
        self.N = N
        self.R = 50
        self.h_center = h_center
        ####Continuation
        self.xCOM = self.define_global_parameter(xCOM=0) 
        self.V1 = self.define_global_parameter(V1=1e4)
        #self.vx   = self.get_global_parameter("vx")
    
    def omega(self,H):
        return (H**(-5)/5 - H**(-2)/2)
    def define_problem(self):
        mesh = LineMesh(minimum=-self.L/2, size=self.L, N=self.N)
        # Add the mesh (default name is "domain" with boundaries "left" and "right")
        self.add_mesh(mesh)
        x = var("coordinate_x")
        h = var('h')

        # Assemble the system
        equations = ThinfilmEquation()  # create a Poisson equation with source g=1
        equations += TextFileOutput()  # Add a simple text file output
        #equations += InitialCondition(h=1+0.5*np.cos(np.pi*2*x/self.L))
        h_init = self.h_center * (1 - x**2 / self.R**2)  # Height functions of the droplets
        hout_teo = ha*(1+(1/3)*(self.omega(self.h_center/ha)-self.omega(1))/(self.h_center/ha-1)) 
        h_init = maximum(h_init, hout_teo)  # Initial height: maximum of h1, h2 and precursor
        equations += InitialCondition(h=h_init)
        equations += SpatialErrorEstimator(h=1.0)
        equations += PeriodicBC("right", offset=self.L) @ "left"
        equations += IntegralObservableOutput(filename=f'obs')      # Stationary state

        ####Continuation
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


        ####end of Continuation


        equations += IntegralObservableOutput(file_ext=["pvd"])
        # Also add an output. This VTU output can be viewed with Paraview
        equations += MeshFileOutput()

        self.add_equations(equations @ "domain")  # Add the equation system on the domain named "domain"
        # The maximum height and its x-coordinate (droplet position) are measured from domain*.txt


if __name__ == "__main__":
    with ThinfilmProblem(L=200,N=512) as problem:#N=int(0.5*1e3/15*eps**(-1))) as problem:
                #problem.set_output_directory(f"pbc_hcenter_{h_center:g}_c={amplitude_C:g}_eps={eps:g}")
                problem.set_output_directory(f"./continuation_testing/")
                # Maximum refinement level
                #problem.max_refinement_level = 15
                #problem.max_permitted_error = 0.0005
                #problem.min_permitted_error = 0.00005
                #problem.run(endtime=20, startstep=0.01, maxstep=0.1, outstep=True, temporal_error=1, spatial_adapt=problem.max_refinement_level)
                problem.run(endtime=5e5, timestep=1e-2, outstep=True, maxstep=1e4, temporal_error=1e-2,spatial_adapt=15,newton_solver_tolerance=1e-10)
                problem.quiet()

                #problem.solve()  # Solve the problem
                print(problem.get_mesh("domain").evaluate_all_observables())
                problem.output()  # Write output

                #end of transient

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
                                                    spatial_adapt=10,
                                                    max_ds=1e+4,
                                                    newton_solver_tolerance=1e-10,
                                                    min_ds=1e-13
                                                    )
                    problem.output()
                    write_to_cont_file() # write the updated state
