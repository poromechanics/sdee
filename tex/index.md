# Summary {#sec:ch5_sec8}

This chapter extends the v-ST/FEM presented in Chapter
[\[ch:ch3\]](#ch:ch3){reference-type="ref" reference="ch:ch3"} for the
dynamic soil-structure interaction problems. In dynamic SSI problem,
unbounded soil domain is truncated by placing artificial boundaries at
some distance from the area of interest. Viscous boundary condition of
Lysmer-Kuhlemeyer [@Lysmer1969] is modified by introducing additional
boundary terms related to the free-field response of unbounded soil
domain to facilitate the energy flow from far field to computation
domain. In the computer program, modified viscous boundary condition is
treated as a combination of various traction boundary conditions;
traction boundary condition due to dashpots, free-field motion, and
input seismic motion. It is found that the traction boundary condition
due to dashpots introduces space-time dashpot matrices
$\left[ \mathbf{C}_{\infty}\right]$ which contribute to the space-time
tangent matrix. The traction boundary condition due to input seismic
motion and free-field motion introduce corresponding space-time nodal
vectors. Furthermore, it is shown how multiple soil-column problems can
be solved (with various boundary conditions) in order to compute the
free-field response. Moreover, the computation of free-field response
does not depend upon the total response of soil and structure. Thanks to
this weak coupling, the soil-column problem can be solved first, and
then the total response of soil and structure can be computed by using
this free-field response. In this way, the proposed space-time finite
element formulation is applicable to a wide class of soil-structure
interaction problem.

Afterwards, a dynamic dam-soil interaction problem is considered to
validate the formulation and computer implementation of v-ST/FEM. In
this problem, a dam (without the reservoir) resting on an elastic
half-space is subjected to the horizontal component of the earthquake
motion. The material damping in both dam and soil domain is modeled
using the Rayleigh damping. The results obtained by proposed scheme are
validated by solving the same problem using the semi-discrete FEM with
classical Newmark-$\beta$ method. Results obtained by two methods are
compared and found to be nearly identical.



[^3]: Recall that Newmark-$\beta$ method has zero algorithmic damping
