# Numerical examples {#sec:ch3-sec7}

## Primary wave propagation in homogeneous linear elastic medium

A theoretical model of homogenous isotropic linear elastic media
occupying a square domain,
$\Omega=\left[0,L\right]\times\left[0,L\right]$, is considered in this
section for a numerical analysis of v-ST/FEM. The governing equations of
the problem are given by Eqs.
([\[eq:ch3-eq-123\]](#eq:ch3-eq-123){reference-type="ref"
reference="eq:ch3-eq-123"}--[\[eq:ch3-eq-127\]](#eq:ch3-eq-127){reference-type="ref"
reference="eq:ch3-eq-127"}). Further, it is assumed that the elastic
media is subjected to the periodic boundary conditions with no external
body force (i.e., $b_{i}=0$ in Eq.
[\[eq:ch3-eq-123\]](#eq:ch3-eq-123){reference-type="ref"
reference="eq:ch3-eq-123"}). The initial conditions for displacement and
velocity field corresponding to the Eq.
[\[eq:ch3-eq-126\]](#eq:ch3-eq-126){reference-type="eqref"
reference="eq:ch3-eq-126"} and Eq.
[\[eq:ch3-eq-127\]](#eq:ch3-eq-127){reference-type="eqref"
reference="eq:ch3-eq-127"} are given as: $$\begin{aligned}
{{\mathbf{u}}_{0}}\left({\mathbf{x}}\right) & ={{\mathbf{d}}_{0}}\cos\left({k{\mathbf{x}}\cdot{\mathbf{\hat{r}}}}\right), & {{\mathbf{v}}_{0}}\left({\mathbf{x}}\right) & =ck{{\mathbf{d}}_{0}}\sin\left({k{\mathbf{x}}\cdot{\mathbf{\hat{r}}}}\right)\label{eq:ch3-eq-162}
\end{aligned}$$ where $k$ is the wave-number, $\hat{\mathbf{r}}$ is the
direction vector of the wave propagation,
$\mathbf{d}_{0}\in\mathbb{R}^{2}$ denotes the direction of the motion of
particles of the medium, and $c$ is related to the speed of the wave in
the medium.

The above-mentioned initial conditions create a plane P-wave if
$\hat{\mathbf{r}}$ vector is parallel to $\mathbf{d}_{0}$. The
analytical solutions for displacement and velocity can be given by
[@Achenbach1973a]: $$\begin{aligned}
{\mathbf{u}}\left({{\mathbf{x}},t}\right) & ={{\mathbf{d}}_{0}}\cos\left[{k\left({{\mathbf{x}}\cdot{\mathbf{\hat{r}}}-{c_{p}}t}\right)}\right] & {\mathbf{v}}\left({{\mathbf{x}},t}\right)={c_{p}}k{{\mathbf{d}}_{0}}\sin\left[{k\left({{\mathbf{x}}\cdot{\mathbf{\hat{r}}}-{c_{p}}t}\right)}\right]\label{eq:ch3-eq-163}
\end{aligned}$$ where $${c_{p}}=\sqrt{\frac{{\lambda+2\mu}}{\rho}}$$ is
the speed of the P-wave.

::: {#tab:ch3-tab-2}
  $L$               $\lambda$   $\mu$   $\rho$         $k$            $\mathbf{d}_{0}$            $\hat{\mathbf{r}}$           $c_{p}$
  ---------------- ----------- ------- -------- ----------------- ------------------------ --------------------------------- ---------
  \[6pt\] $10.0$      $2.0$     $1.0$   $1.0$    $\sqrt{2}\pi/5$   $\left[1,1\right]^{T}$   $\left[1,1\right]^{T}/\sqrt{2}$      $2.0$
  \[6pt\]

  :  List of constants and parameter values used for the P-wave
  propagation problem. All variables are dimensionless.
:::

The values of constants and the parameters used for solving the problem
are given in Table [2](#tab:ch3-tab-2){reference-type="ref"
reference="tab:ch3-tab-2"}, and Fig.
[19](#fig:ch3-fig-19){reference-type="ref" reference="fig:ch3-fig-19"}
depicts the physical dimensions of the problem. All variables have been
made dimensionless. Letting $\mathbf{u}^{h}(\mathbf{x},t)$ and
$\mathbf{v}^{h}(\mathbf{x},t)$ denote the numerical solutions computed
by employing the v-ST/FEM, the error in displacement field
${E_{u}}\left({{t_{n}}}\right)$ and the error in velocity field
${E_{v}}\left({{t_{n}}}\right)$ are then defined as: $$\begin{aligned}
{E_{u}}\left({{t_{n}}}\right) & :={\left\Vert {{{\mathbf{u}}^{h}}\left({{\mathbf{x}},{t_{n}}}\right)-{\mathbf{u}}\left({{\mathbf{x}},{t_{n}}}\right)}\right\Vert _{2}}, & {E_{v}}\left({{t_{n}}}\right) & :={\left\Vert {{{\mathbf{v}}^{h}}\left({{\mathbf{x}},t_{n}^{-}}\right)-{\mathbf{v}}\left({{\mathbf{x}},{t_{n}}}\right)}\right\Vert _{2}}\label{eq:ch3-eq-164}
\end{aligned}$$ where $\left\Vert {\,\cdot\,}\right\Vert _{2}$ denotes
the $L_{2}$ norm given by
$${\left\Vert {\mathbf{u}}\right\Vert _{2}}={\left[{\int_{\Omega}^ {}{{\mathbf{u}}\left({{\mathbf{x}},t}\right)\cdot{\mathbf{u}}\left({{\mathbf{x}},t}\right)d\Omega}}\right]^{1/2}}\label{eq:ch3-eq-165}$$

![(a) Schematic diagram of the spatial domain for P-wave propagation
problem, (b) bilinear quadrilateral element (Quad4), (c) linear
triangular element (Tria3)](figures/ch3-fig-19){#fig:ch3-fig-19}

![ Rate of convergence of the solutions in space computed at time
$t=1.0$ sec: (a) displacement, and (b) velocity
](figures/ch3-fig-20){#fig:ch3-fig-20 width="100%"}

![ Rate of convergence of the solutions in time domain at time $t=35.0$
sec: (a) displacement, and (b) velocity.
](figures/ch3-fig-21){#fig:ch3-fig-21 width="100%"}

The numerical experiments for determining the convergence rate of the
solution in the space domain are performed on two sequences of regular
linear triangular (Tria3) and bilinear quadrilateral (Quad4) meshes (see
Fig. [19](#fig:ch3-fig-19){reference-type="ref"
reference="fig:ch3-fig-19"}), while keeping the time-step fixed
$\Delta t=0.1$ sec. Each sequence consists of four meshes with a
decreasing mesh size. In Fig. [20](#fig:ch3-fig-20){reference-type="ref"
reference="fig:ch3-fig-20"}, the L2 norm of the errors in the
displacement and velocity fields at time $t=1.0$ sec are given in
relation to mesh spacing parameter $h$. Based on the convergence
results, it can be stated that the v-ST/FEM formulation is nearly
second-order accurate in the space for both Quad4 and Tria3 elements. In
Fig. [20](#fig:ch3-fig-20){reference-type="ref"
reference="fig:ch3-fig-20"}, it is observed that the error for the
triangular spatial mesh (Tria3) is less than that of the quadrilateral
spatial element (Quad4) for the same mesh spacing. This can be
attributed to the perfect alignment of the diagonal of the triangular
elements with the characteristic lines of the wave propagation.

![ (a) Exact displacement (x-component) waveforms, (b) displacement
(x-component) waveforms obtained by using the v-ST/FEM, (c) exact
velocity (x-component) waveforms, (d) velocity (x-component) waveforms
obtained by using the v-ST/FEM. ](figures/ch3-fig-22){#fig:ch3-fig-22
width="90%"}

Furthermore, the convergence of the solutions (displacement and
velocity) in the time domain is illustrated in Fig.
[21](#fig:ch3-fig-21){reference-type="ref" reference="fig:ch3-fig-21"}.
It can be seen that both displacement and velocity fields computed by
the present method are third order accurate in time. Note that the
results presented in Fig. [21](#fig:ch3-fig-21){reference-type="ref"
reference="fig:ch3-fig-21"} are consistent with Eq.
[\[eq:ch3-eq-114\]](#eq:ch3-eq-114){reference-type="eqref"
reference="eq:ch3-eq-114"}.

Fig. [22](#fig:ch3-fig-22){reference-type="ref"
reference="fig:ch3-fig-22"} illustrates the spatial variation of the
displacement and velocity fields obtained by v-ST/FEM. The results are
obtained at time $t=35$ seconds with linear triangular mesh of size
$h=0.25$ and uniform time-step of size $\Delta t=0.1$ sec. The results
advocate the ability of v-ST/FEM to maintain the high-order accuracy of
the long-term solutions, especially the displacement and velocity
waveforms. This can be attributed to the very low numerical dissipation
and dispersion characteristics of the present method (see also Fig.
[15](#fig:ch3-fig-15){reference-type="ref" reference="fig:ch3-fig-15"}
and Fig. [16](#fig:ch3-fig-16){reference-type="ref"
reference="fig:ch3-fig-16"}). To further emphasize the accuracy of the
long-term solutions, the time histories of the computed displacement and
velocity at the midpoint P1 (see Fig.
[19](#fig:ch3-fig-19){reference-type="ref" reference="fig:ch3-fig-19"})
and the corresponding relative errors are presented in Fig.
[23](#fig:ch3-fig-23){reference-type="ref" reference="fig:ch3-fig-23"}.

![ Temporal variation of the (a) displacement, (b) relative error in
displacement, (c) velocity, and (d) relative error in velocity, at the
midpoint P1 computed by using v-ST/FEM.
](figures/ch3-fig-23){#fig:ch3-fig-23 width="90%"}

## Impulsive response of a fixed-free pile

In this section, we consider a pile of length $L=50$ m with a unit
cross-section area having one fixed end and one end loaded by an axial
impulsive force given by a step function, as shown in Fig.
[24](#fig:ch3-fig-24){reference-type="ref" reference="fig:ch3-fig-24"}.
The mass density $\rho$ is  $2500.0$ $kg/m^{3}$  and the Young's modulus
$E$ is  $1.0\times10^{10}$ $N/m^{2}$ . The magnitude of the implusive
force is $1.0\times10^{6}$ $N$. Under these circumstances, the
analytical solutions for the stress and velocity fields are
discontinuous in the spatial domain and given by the step functions.
Furthermore, the displacements are given by the piecewise continuous
linear functions [@Cormeau1991; @Li1998; @Verruijt2009].

![ Geometry, boundary conditions, and impulse loading for fixed-free
pile problem. ](figures/ch3-fig-24){#fig:ch3-fig-24}

To solve the problem by using v-ST/FEM, uniform linear spatial elements
of size $h=0.1m$ and a uniform time-step size $\Delta t=10^{-4}$ sec
have been adopted for discretizing the space and time domain,
respectively. To assess the performance of the present method, the
problem is also solved with semi-discretized FEM techniques using the
same mesh parameters. In the latter case, the trapezoidal rule and the
HHT-$\alpha$ method (with $\alpha=-1/3$) have been used as the
time-stepping algorithms. The results of the stress field, velocity
field, and displacement field obtained by different schemes are compared
in Fig. [25](#fig:ch3-fig-25){reference-type="ref"
reference="fig:ch3-fig-25"}, Fig.
[26](#fig:ch3-fig-26){reference-type="ref" reference="fig:ch3-fig-26"},
and Fig. [27](#fig:ch3-fig-27){reference-type="ref"
reference="fig:ch3-fig-27"}, respectively.

It is observed that the solutions for the stress and velocity fields
obtained by the semi-discrete algorithms contain severe oscillations.
Even worse, no improvement whatsoever is obtained when refining the mesh
spacing or the time-step size. Further, in the case of the Newmark-beta
method, these oscillations are present in the whole spatial domain; and
thus, the accuracy of the stress and velocity fields deteriorate over
time (see Fig. [25](#fig:ch3-fig-25){reference-type="ref"
reference="fig:ch3-fig-25"} and Fig.
[26](#fig:ch3-fig-26){reference-type="ref" reference="fig:ch3-fig-26"}).
The poor performance with the Newmark-beta method is due to the absence
of algorithmic damping, as discussed in previous sections. The
HHT-$\alpha$ method improves the solutions by attenuating higher
frequencies; however, the results are not satisfactory for the selected
time-step size as the oscillations are still present. It is remarkable
that v-ST/FEM completely localizes the oscillations in the stress and
velocity fields near the point of discontinuity and yields very accurate
solutions. The localization phenomenon can be attributed to the presence
of jump discontinuity in the velocity field. The jump in the velocity
field adds artificial viscous damping to the system, subsequently,
increasing the accuracy and stabilizing the solutions
[@Hulbert1990; @Johnson1993]. Furthermore, the presence of overshooting
and undershooting around the point of discontinuity is related to the
famous Gibbs phenomenon in the Fourier analysis [@Olver2016].

![ Stress field computed by employing the Newmark-beta method (First
column), HHT-$\alpha$ (Second column), and v-ST/FEM (Third column) at
various timesteps with linear spatial elements. Direction of wave
propagation is denoted by the arrow, and dotted lines represent the
analytical solutions. ](figures/ch3-fig-25){#fig:ch3-fig-25
width="100%"}

![ Velocity field computed by employing the Newmark-beta method (First
column), HHT-$\alpha$ (Second column), and v-ST/FEM (Third column) at
various time-steps with linear spatial elements. Direction of wave
propagation is denoted by the arrow, and dotted lines represent the
analytical solutions. ](figures/ch3-fig-26){#fig:ch3-fig-26
width="100%"}

![ Displacement field computed by employing the Newmark-beta method
(First column), HHT-$\alpha$ (Second column), and v-ST/FEM (Third
column) at various timesteps with linear spatial elements. Direction of
wave propagation is denoted by the arrow.
](figures/ch3-fig-27){#fig:ch3-fig-27 width="100%"}

## Dynamic plate load test (DPLT)

In this section, an attempt is made to validate v-ST/FEM by simulating
the dynamic plate loading test (DPLT) using a light falling weight
deflectometer (LFWD). DPLT using LFWD is a non-destructive technique for
a quick assessment of the field compaction quality. In DPLT, a rigid
circular loading plate of radius $0.15$ m is placed on the soil surface
and subjected to an impulse load generated by the falling weight from a
specified height onto the plate. Subsequently, the induced soil
movements are recorded and the dynamic resilience modulus of the tested
material is computed by some empirical relations. A complete description
of the test can be found elsewhere [@Adam2009; @Tawfik2017].

![ Geometry, boundary conditions and spatial mesh adopted for simulating
DPLT using v-ST/FEM. ](figures/ch3-fig-28){#fig:ch3-fig-28 width="100%"}

In the present study, DPLT is simulated by a two-dimensional
axisymmetric v-ST/FEM model, where the center of the loading plate is
positioned along the axis of symmetry. The geometry, boundary
conditions, and spatial mesh of the finite element model are depicted in
Fig. [28](#fig:ch3-fig-28){reference-type="ref"
reference="fig:ch3-fig-28"}. There are $5201$ linear triangular elements
(Tria3) and $2688$ nodes present in the spatial mesh. The impulse load
due to the falling weight is modeled by using an equivalent uniform
vertical stress pulse of amplitude $100$ kPa and a time duration of $20$
ms [@Adam2009]. The stress pulse $h(t)$ acting on the loading plate is
defined by a half sine wave as:
$$h(t)=-10^{5}\sin(50\pi t)\,\text{N/m}{}^{2}$$

The total simulation time $T$ is set to $30$ ms, and the linear time
elements of size $\Delta t=1$ ms have been adopted for discretizing the
time domain. Subsequently, the results computed by the proposed method
are compared with the two DPLT studies available in the literature. The
first study is denoted as S1 where the in-situ LFWD test was conducted
by @Tawfik2017, and the second study is denoted by S2 in which the
numerical investigation was conducted by @Adam2009 using the Boundary
Element Method (BEM). In these studies (i.e., S1 and S2), the soils have
different values for Young's modulus $E$ and common values for Poisson's
ratio $\nu$ and mass density $\rho$, as shown in Table
[3](#tab:ch3-tab-3){reference-type="ref" reference="tab:ch3-tab-3"}.

::: {#tab:ch3-tab-3}
  Elastic parameters              Method         $d_{max}$     $p_{max}$
  ------------------------ -------------------- ----------- ------------
  \[6pt\] $E=65$ Mpa             v-ST/FEM        $0.33$ mm    $90.0$ kPa
  $\nu=0.212$               In-situ study (S1)   $0.35$ mm     $100$ kPa
  $\rho=2000$ $kg/m^{3}$
  \[6pt\] $E=32$ Mpa             v-ST/FEM        $0.67$ mm    $90.0$ kPa
  $\nu=0.212$                 BEM study (S2)     $0.65$ mm    $93.1$ kPa
  $\rho=2000$ $kg/m^{3}$
  \[6pt\]

  :  List of elastic material parameters, and the results of maximum
  plate deflection and maximum soil-plate normal contact stress obtained
  by different schemes.
:::

The results of the time histories of the vertical displacement and the
velocity of the loading plate, and the soil-plate contact stress are
plotted in Fig. [29](#fig:ch3-fig-29){reference-type="ref"
reference="fig:ch3-fig-29"}. The maximum deflection of the plate
$d_{max}$ obtained by v-ST/FEM is about $0.33$ mm and $0.67$ mm
corresponding to the material parameters used in the studies S1 and S2,
respectively. From Table [3](#tab:ch3-tab-3){reference-type="ref"
reference="tab:ch3-tab-3"}, it is evident that these computed values for
$d_{max}$ are in good agreement with those reported in studies S1 and
S2. In the present study, the maximum soil-plate normal contact stress
$p_{max}$ is about $90$ kPa; which lower than the values reported in
studies S1 and S2 (see Table [3](#tab:ch3-tab-3){reference-type="ref"
reference="tab:ch3-tab-3"}). This may be due to the use of linear
triangular elements and the absence of the soil-plate interface elements
in the present formulation. Furthermore, the maximum deflection of the
plate occurs at time $t_{d}=11$ ms, and the maximum contact stress
$p_{max}$ occurs at time $t_{p}=10$ ms. The displacement contours at
these time steps are depicted in Fig.
[30](#fig:ch3-fig-30){reference-type="ref" reference="fig:ch3-fig-30"},
where displacements are normalized with respect to the maximum
deflection of the plate.

![ Temporal variation of the vertical displacement of plate (top-left),
the vertical velocity of plate (top-right), and soil-plate contact
stress (bottom-left), and variation of normal contact stress with the
plate displacement (bottom-right) obtained by v-ST/FEM.
](figures/ch3-fig-29){#fig:ch3-fig-29 width="100%"}

![ Normalized contours of the vertical displacements computed by
v-ST/FEM at various time-steps. ](figures/ch3-fig-30){#fig:ch3-fig-30
width="100%"}

# Summary {#sec:ch3-sec8}

In this chapter, the concept of time-discontinuous Galerkin (TDG) method
is presented. The details of single-field and two-field TDG schemes,
such as stability, the order of accuracy, and algorithmic damping, are
discussed. Subsequently, a velocity-based single-field TDG space-time
finite element method has been been developed for elastodynamics
problems. The main characteristics of this method is summarized below.

1.  In v-ST/FEM velocity is the primary unknown which is discontinuous
    in time and continuous in space. The time-continuity of the velocity
    field is satisfied in a weak sense.

2.  The displacement field is obtained by time-integration of the
    velocity field in a post-processing step. These displacement field
    are then used to compute stress-field. Both stress and displacement
    field remain continuous in time.

3.  The advantage of v-ST/FEM is that it involves less number of
    unknowns--unlike other ST/FEM-- which makes v-ST/FEM applicable to
    the large-scale practical problems at relatively low computational
    cost.

4.  It is demonstrated that the present method is unconditionally stable
    and third-order accurate in time for linear interpolation of
    velocity in time.

5.  The numerical dissipation (amplitude decay) and dispersion
    (phase-delay) of v-ST/FEM is smaller than the two-field ST/FEM and
    semi-discrete algorithms.

6.  In the case of a shock problem, the performance of v-ST/FEM was
    found to be superior to HHT-$\alpha$ and Newmark-beta method.
    v-ST/FEM is able to remove the spurious oscillations unlike
    HHT-$\alpha$ and Newmark-beta algorithms.

7.  The numerical characteristics of the v-ST/FEM scheme, therefore,
    make it highly suitable for computing the response of bodies
    subjected to dynamic loading conditions, such as fast-moving loads,
    impulsive loading, long-duration seismic loading, among others.

8.  It is shown that the proposed scheme is non-dissipative in
    high-frequency regime, and can attenuate only the middle band of
    frequencies. Therefore, the only drawback to the present methodology
    is that it does not include a parameter to control the numerical
    dissipation of the high-frequencies.

[^1]: Eq. [\[eq:ch3-eq-1\]](#eq:ch3-eq-1){reference-type="eqref"
    reference="eq:ch3-eq-1"} characterizes the modal equations of a
    linear parabolic partial differential equation (e.g., heat diffusion
    equation)

[^2]: The boundary of a one-dimensional finite element comprises
    end-points only.

[^3]: Informally speaking, the purpose of introducing the representative
    value $u^{*}$ in TDG/FEM is to facilitate the weak enforcement of
    the initial condition in the time domain. Further, $u^{*}$ connects
    the adjacent time-slabs by transferring the information about the
    solution from one time-slab to another time-slab. In TDG/FEM, this
    so-called connection between the time-slabs is weakly enforced which
    in turn allows the computation to be performed locally in an
    individual time-slab.

[^4]: The internal nodes may be located at specific locations in $I_{n}$
    which are determined from the zeros of the orthonormal polynomials
    like the Jacobi polynomials. This approach can be used for
    developing the high-order time accurate spectral finite element
    schemes [@Hesthaven2007 see].

[^5]: Homogeneous form of Eq.
    [\[eq:ch3-eq-1\]](#eq:ch3-eq-1){reference-type="eqref"
    reference="eq:ch3-eq-1"} is obtained by setting $f=0$, i.e.,
    $$\frac{{du}}{{dt}}+\lambda u=0\quad t\in\left[{0,T}\right],u(0)=u_{0}$$

[^6]: Lax equivalence theorem may be stated as *consistency and
    stability of an algorithm are the necessary and sufficient
    conditions for the convergence of an algorithm.*

[^7]: Henceforth uv-TDG/FEM will be used to denote the two-field TDG/FEM

[^8]: To prove that total energy of the spring-mass system is constant,
    multiply Eq.
    [\[eq:ch3-eq-93\]](#eq:ch3-eq-93){reference-type="eqref"
    reference="eq:ch3-eq-93"} with velocity $v=\frac{du}{dt}$,
    $$v\frac{{dv}}{{dt}}+\omega_{n}^{2}u\frac{{du}}{{dt}}=0$$ then
    $$\frac{d}{{dt}}\left({\frac{1}{2}{v^{2}}}\right)+\frac{d}{{dt}}\left({\frac{1}{2}\omega_{n}^{2}{u^{2}}}\right)=0$$
    or
    $$\frac{d}{{dt}}\left({\frac{1}{2}{v^{2}}+\frac{1}{2}\omega_{n}^{2}{u^{2}}}\right)=0$$
    therefore,
    $$\frac{1}{2}{v^{2}}+\frac{1}{2}\omega_{n}^{2}{u^{2}}={\text{constant}}$$

[^9]: Normalized total energy at any time $t$ is defined by total energy
    at time $t$ divided by initial total energy, i.e., $$$$ = $$$$

[^10]: The term *semi-discrete algorithms* is used for collectively
    referencing the second order accurate implicit time-stepping
    algorithms viz., the trapezoidal rule (Newmark method with
    $\gamma=0.5$ and $\beta=0.25$), the Wilson-$\theta$ method with
    $\theta=1.4$ @Bathe1976, the HHT-$\alpha$ method with $\alpha=-0.3$
    @Hilber1977, and Houbolt's method @Houbolt1950. For details about
    these semi-discrete algorithms see @Hughes1983 [@Hughes2012].

[^11]: Eq. [\[eq:ch3-eq-93\]](#eq:ch3-eq-93){reference-type="eqref"
    reference="eq:ch3-eq-93"} is given by following
    $$\frac{d^{2}u}{dt^{2}}+\omega_{n}^{2}=0$$
