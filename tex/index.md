
# Summary {#sec:ch7-sec9}

The chapter presents space-time finite element formulations for the
problems involving dynamic response of solids and structures with
nonlinear stress-strain relationships. The problem of dynamic
interaction between the concrete gravity dam and reservoir is taken as a
model problem, in which a generalized nonlinear stress-strain
relationship is used to describe the material behavior of concrete in
the dam. The foundation underneath the dam-reservoir (DR) system is
assumed to be perfectly rigid. The governing equations describing the
dynamic interaction between dam and reservoir constitute a system of
linear-nonlinear coupled equations, in which linear equations govern the
reservoir domain and nonlinear equations govern the solid domain.
Subsequently, v-ST/FEM weak form is derived and then discretized by
using the space-time finite elements. Accordingly, an unsymmetrical
system of linear-nonlinear algebraic equations describes the discrete
form of the v-ST/FEM weak form. A block-iterative scheme is devised to
enforce the coupling between the solid and fluid domain. In each
iteration of v-ST/FEM with the block-iterative scheme, the linearized
equations of the solid domain are first solved to compute the increments
in the velocity field. Subsequently, the total velocities are corrected
and then used for computing the trial values of hydrodynamic pressures
in the reservoir by solving the linear equation for the reservoir
domain. Iterations are performed until the convergence in the solutions
is achieved. In each iteration of the proposed scheme, therefore, linear
equations for the solid and fluid domain are solved, separately, which
significantly decreases the computation cost.

In the present problem, nonlinearity is caused by only the presence of
stress term in the v-ST/FEM weak form. Accordingly, two v-ST/FEM schemes
are proposed for the time integration of the space-time nodal vectors
and matrices comprising the stress term. The first scheme uses the
three-point Gauss-Lobatto quadrature rule for time-integration, whereas
the second scheme uses the two-point Gauss-Legendre quadrature rule.
Further, the Gauss-Lobatto quadrature rules include the endpoints of the
time interval, while the latter approach includes the points only
located inside the time interval.

Subsequently, to evaluate the numerical performance of the v-ST/FEM
schemes dynamic fracture analyses of the concrete dam are performed.
Numerical simulations are performed for following two cases: the first
case ignores the coupling between the dam and reservoir, and the second
case includes the coupling between dam and reservoir. In the former
case, hydrodynamic pressures are set to zero in the discretized
equations of the solid domain, consequently, block-iterative scheme
reduces to the Newton method.

A co-axially rotating crack model (CRCM) with exponential strain
softening rule is employed to model the fracture of the concrete. Finite
element implementation of the CRCM requires the average of integration
point strains to determine the cracking behavior of the element as a
whole. In v-ST/FEM, one can opt the average of either spatial
integration point strains or the space-time integration point strains.
Accordingly, three v-ST/FEM schemes are devised for the dynamic fracture
analysis of the concrete dam: *v-ST/FEM-1* corresponds to the
three-point Gauss-Lobatto rule and the spatially averaged strains,
*v-ST/FEM-2* corresponds to the two-point Gauss-Legendre rule and the
spatially averaged strains, and *v-ST/FEM-3* corresponds to the
two-point Gauss-Legendre rule and the space-time averaged strains.

Numerical simulations presented in this chapter confirms that all
v-ST/FEM schemes are consistent with each other. All schemes can
successfully simulate the crack propagation in the concrete dam during
the earthquake loading. Dynamic fracture analysis of concrete dam
involves the rapid change of the stiffness of the dam due to the crack
opening-closing-reopening cycles. The results indicates that
*v-ST/FEM-1* is the most robust algorithm among the three v-ST/FEM
schemes discussed here. Further, total number of iteration increases at
the instant of crack closing and smaller time steps may be required to
achieve the convergence. During the crack opening and reopening,
however, number of iterations are relatively low which can be attributed
to the high order accuracy of the v-ST/FEM schemes.


[^2]: In a typical numerical integration method based on quadrature
    rules, the integrand is evaluated a finite set of points called the
    *integration points* and a weighted sum of these values is used to
    approximate the integral. The choice of integration points and
    weights depends upon the specific method used and the accuracy
    required from the approximation.

[^3]: The CRCM model was first proposed by [@Bhattacharjee1993] and
    recently modified by [@Calayir2005]

[^4]: The fracture energy is defined as the energy per unit area
    required to form a fracture surface

