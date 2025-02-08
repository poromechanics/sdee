# Wave propagation in elastic solids {#sec-elastic-waves-in-solid}

Linear and nonlinear dynamic response of the deformable media can be
described by a wave propagation phenomenon. A deformable media subjected
to a transient loading condition produces mechanical waves; as elements
of the medium are deformed the disturbance is transmitted from one point
in space to the next. The local mechanical disturbances of a medium is
not instantaneously detected at positions that are at a distance from
the region of excitation. It takes time for a disturbance to propagate
from its source to other positions. In this way as the disturbance
propagates through the medium it carries along amounts of energy due to
motions of particles of medium about their equilibrium position.

Further, the transmission of mechanical wave depends mainly upon the
deformation characteristics and inertia of the medium. As it will be
seen later both stiffness and inertia of a medium tend decrease the wave
speed. A rigid medium, (e.g. rocks with very high stiffness) deforms
insignificantly and the mechanical disturbance travels almost
instantaneously in it. Similarly, a hypothetical massless medium
(i.e. no inertia) allows mechanical wave to travel without any delay.
These concepts have been used widely to model the interaction between
the soil/rock foundation and the earth structures such as concrete
gravity dams. It is important to note that the inertia of a system first
offers resistance to motion, but once the medium is in motion inertia
with the resilience of the medium tends to sustain the motion
[@Achenbach1973a].

## Governing equations of motion

Consider a body $\mathfrak{B}$ occupying a regular domain
$\Omega \subset \mathbb{R}^{n_{sd}}$ in the space which may be bounded
or unbounded. $n_{sd}$ is the number of spatial dimensions. Let $\Gamma$
denotes the boundary of the domain, and
$\bar \Omega = \Omega \cup \Gamma$ be the closure of $\Omega$. Let
indices ${i,j,k,l}$ take values from $1,\cdots, n_{sd}$. The system of
equations describing the motion of a homogeneous, isotropic, linearly
elastic body consists of the Cauchy's equation of the motion,
generalized Hooke's law and linearized strain-displacement relationship:

$$\label{eq:Momentum}
\frac{{\partial {\sigma _{ij}}}}{{\partial {x_j}}} + \rho {b_i} = \rho \frac{{{\partial ^2}{u_i}}}{{\partial {t^2}}}$$

$$\label{eq:Hooke}
\sigma_{ij} = \lambda \epsilon_{kk} \delta_{ij} + 2 \mu \epsilon_{ij}$$

$$\label{eq:Strain-Disp}
\epsilon_{ij} = \frac{1}{2}\left( {\frac{{\partial {u_i}}}{{\partial {x_j}}} + \frac{{\partial {u_j}}}{{\partial {x_i}}}} \right)$$

Here, $\lambda$ and $\mu$ are the Lame's parameters. Subsequently, the
displacement equation of motion Eq.
[\[eq:Momentum-Disp-Form\]](#eq:Momentum-Disp-Form){reference-type="eqref"
reference="eq:Momentum-Disp-Form"} can be obtained by using Eq.
[\[eq:Hooke\]](#eq:Hooke){reference-type="eqref" reference="eq:Hooke"}
and [\[eq:Strain-Disp\]](#eq:Strain-Disp){reference-type="eqref"
reference="eq:Strain-Disp"} in Eq.
[\[eq:Momentum\]](#eq:Momentum){reference-type="eqref"
reference="eq:Momentum"}

$$\label{eq:Momentum-Disp-Form}
\mu \frac{{{\partial ^2}{u_i}}}{{\partial x_j^2}} + \left( {\lambda  + \mu } \right)\frac{{{\partial ^2}{u_j}}}{{\partial {x_j}\partial {x_i}}} + \rho {b_i} = \rho \frac{{{\partial ^2}{u_i}}}{{\partial {t^2}}}$$

To obtain a unique solution of the problem boundary conditions on the
boundary $\Gamma$ and initial state of the body must be prescribed. Some
of the commonly used boundary conditions are mentioned below.

1.  *Displacement boundary conditions*: Components of displacement field
    $u_i(\mathbf{x},t)$ are prescribed on the boundary

    $$u_i(\mathbf{x},t) = g_i(\mathbf{x},t) \quad \textit{on} \quad \Gamma$$

2.  *Traction boundary conditions*: This boundary condition relates the
    stress $\sigma_{ij}$ to the externally applied surface force $t_i$
    using the Cauchy's formula

    $$\sigma_{ij}{n_j} = t_{i} \quad \textit{on} \quad \Gamma$$

3.  *Mixed boundary conditions*: On a part of the boundary
    $\Gamma_{i}^{g}\subset \Gamma$ displacement $u_i$ is prescribed and
    on the remaining part of the boundary
    $\Gamma_{i}^{h} \subset \Gamma$ traction boundary condition is
    imposed.

    $$u_i(\mathbf{x},t) = g_i(\mathbf{x},t) \quad \textit{on} \quad \Gamma_{i}^{g}$$

    $$\sigma_{ij}{n_j} = t_{i} \quad \textit{on} \quad \Gamma_{i}^{h}$$

    with

    $$\Gamma_{i}^{g} \cup \Gamma_{i}^{h} = \Gamma \quad \textit{and} \quad \Gamma_{i}^{g} \cap \Gamma_{i}^{h} = \phi$$

The initial state of the body, the displacement and velocity field at
time $t=0$, must be well defined to complete the problem.

$$u_i(\mathbf{x},0) = u_{i}^{0}(\mathbf{x}) \quad \forall \mathbf{x} \in \bar{\Omega}$$

$$\dot{u_i}(\mathbf{x},0) = v_{i}^{0}(\mathbf{x}) \quad \forall \mathbf{x} \in \bar{\Omega}$$

The ordered pair $\left( {{\mathbf{u}},\sigma } \right)$ defines the
elastodynamic state on $\left( \bar{\Omega}  \times T \right)$ with the
displacement field $\mathbf{u}$ and the stress field $\sigma$,
corresponding to given external body force density $\mathbf{b}$, the
mass density $\rho \in \mathbb{R}^{+}$ and the Lame parameters
$\lambda, \mu \in \mathbb{R}^{+}$, if
${\mathbf{u}} \in {C^{2}}\left( {\Omega  \times T} \right) \cap {C^1}\left( {\bar{\Omega}  \times T} \right)$,
$\sigma  \in {C^1}\left( {\bar \Omega  \times T} \right)$ and
${\mathbf{b}} \in {C^1}\left( {\bar \Omega  \times T} \right)$,and Eq.
([\[eq:Momentum\]](#eq:Momentum){reference-type="ref"
reference="eq:Momentum"} --
[\[eq:Strain-Disp\]](#eq:Strain-Disp){reference-type="ref"
reference="eq:Strain-Disp"} ) is satisfied for all
$(\mathbf{x},t) \in \Omega \times T$ along with prescribed initial and
boundary conditions.

## Displacement potentials

The displacement equation of motion Eq.
[\[eq:Momentum-Disp-Form\]](#eq:Momentum-Disp-Form){reference-type="eqref"
reference="eq:Momentum-Disp-Form"} couples all the displacement
components $u_i(\mathbf{x},t)$. This equation, however, can be decoupled
by introducing the concept of displacement potentials [@Achenbach1973a].
For brevity, Eq.
[\[eq:Momentum-Disp-Form\]](#eq:Momentum-Disp-Form){reference-type="eqref"
reference="eq:Momentum-Disp-Form"} is first written in its vector form.

$$\label{eq: Momentum-Disp-Vec-Form}
\mu {\nabla ^2}{\mathbf{u}} + \left( {\lambda  + \mu } \right)\nabla  \otimes \left( {\nabla  \cdot {\mathbf{u}}} \right) + \rho {\mathbf{b}} = \rho \frac{{{\partial ^2}{\mathbf{u}}}}{{\partial {t^2}}}$$

Let ${\mathbf{u}} \in {C^2}\left( {\Omega  \times T} \right)$ and
${\mathbf{b}} \in {C^{1}}\left( {\Omega  \times T} \right)$ satisfies
Eq.
[\[eq: Momentum-Disp-Vec-Form\]](#eq: Momentum-Disp-Vec-Form){reference-type="eqref"
reference="eq: Momentum-Disp-Vec-Form"} for all
$\left( {{\mathbf{x}},t} \right) \in \Omega  \times T$. Let the
Helmholtz decomposition of the external body force density
$\mathbf{b}$[^1] is given by

$$\label{eq: Helmholtz-b}
\mathbf{b} = c_{L}^{2}\nabla{F} + c_{T}^{2}\nabla \times \mathbf{G}$$

Then there exists a scalar function
$\phi :\Omega  \times T \to \mathbb{R}$ and a vector-valued function
$\bm{\psi}:\Omega  \times T \to \mathbb{R}^{n_{sd}}$ such that the
displacement field can be described by

$$\label{eq: Helmholz-Disp-Vec-Form}
{\mathbf{u}} = \nabla \phi  + \nabla  \times \bm{\psi}$$

Further, the displacement potential $\bm{\psi}$ usually satisfies the
following constraint equation.

$$\nabla \cdot \bm{\psi} = 0$$

After substituting Eq.
([\[eq: Helmholtz-b\]](#eq: Helmholtz-b){reference-type="ref"
reference="eq: Helmholtz-b"} --
[\[eq: Helmholz-Disp-Vec-Form\]](#eq: Helmholz-Disp-Vec-Form){reference-type="ref"
reference="eq: Helmholz-Disp-Vec-Form"}) into Eq.
[\[eq: Momentum-Disp-Vec-Form\]](#eq: Momentum-Disp-Vec-Form){reference-type="eqref"
reference="eq: Momentum-Disp-Vec-Form"} one can obtain the following two
uncoupled wave equations.

$$\label{eq: wave-eq-phi}
{\nabla ^2}\phi + F = \frac{1}{{c_L^2}}\frac{{{\partial ^2}\phi }}{{\partial {t^2}}}$$

and

$$\label{eq: wave-eq-psi}
{\nabla ^2}\bm{\psi} + \mathbf{G}  = \frac{1}{{c_T^2}}\frac{{{\partial ^2}\bm{\psi} }}{{\partial {t^2}}}$$

Eq. [\[eq: wave-eq-phi\]](#eq: wave-eq-phi){reference-type="eqref"
reference="eq: wave-eq-phi"} is a scalar wave equation and Eq.
[\[eq: wave-eq-psi\]](#eq: wave-eq-psi){reference-type="eqref"
reference="eq: wave-eq-psi"} represents wave equation in $n_{sd}$
components of $\bm{\psi}$. The characteristic-speed $c_L$ and $c_T$ is
given by Eq.
[\[eq: longi-wave-speed\]](#eq: longi-wave-speed){reference-type="eqref"
reference="eq: longi-wave-speed"} and
[\[eq: trans-wave-speed\]](#eq: trans-wave-speed){reference-type="eqref"
reference="eq: trans-wave-speed"}, respectively.

$$\label{eq:longi-wave-speed}
{c_L} = \sqrt {\frac{{\lambda  + 2\mu }}{\rho }}$$

$$\label{eq: trans-wave-speed}
{c_T} = \sqrt {\frac{{\mu }}{\rho }}$$

## Longitudinal and transverse plane waves

Consider a plane displacement wave traveling with phase velocity $c$ in
the direction $\mathbf{p}$ is given by

$$\label{eq:plane-wave-1}
{\mathbf{u}} = {\mathbf{d}}f\left( {\mathbf{x}} \cdot {\mathbf{p}} - ct  \right)$$

where $\mathbf{d}$ is the direction of motion of the particles,
$\mathbf{x}$ is the position vector of the particle in the space and the
argument of $f(\cdot)$ is called the phase of the wave and given by

$$\label{eq:phase-eta}
\eta  = {\mathbf{x}} \cdot {\mathbf{p}} - ct$$

It is evident from Eq.
[\[eq:plane-wave-1\]](#eq:plane-wave-1){reference-type="eqref"
reference="eq:plane-wave-1"} that at any given time $t$ the plane of
constant phase, i.e. $\mathbf{x}\cdot \mathbf{p} =$ constant, are normal
to the direction of wave propagation. Further, the planes moves in the
direction $\mathbf{p}$ with phase-velocity $c$.

Furthermore, the displacement given by Eq.
[\[eq:plane-wave-1\]](#eq:plane-wave-1){reference-type="eqref"
reference="eq:plane-wave-1"} can satisfy the wave equation Eq.
[\[eq:Momentum-Disp-Form\]](#eq:Momentum-Disp-Form){reference-type="eqref"
reference="eq:Momentum-Disp-Form"} only in two cases:

1.  *Longitudinal wave*[^2]: In this case the motion of particles is in
    the direction of the wave propagation (see Fig.
    [2](#fig:Plane-wave-motion){reference-type="ref"
    reference="fig:Plane-wave-motion"}) and phase-velocity is given by
    Eq.
    [\[eq: longi-wave-speed\]](#eq: longi-wave-speed){reference-type="eqref"
    reference="eq: longi-wave-speed"}. Mathematically speaking,

    $$\label{eq:longi-wave-cond}
    {\mathbf{d}} =  \pm {\mathbf{p}}; \quad c = c_L$$

    and

    $$\label{eq:plane-wave-longi}
    {\mathbf{u}} = \pm {\mathbf{p}}f\left( {\mathbf{x}} \cdot {\mathbf{p}} - c_{L}t  \right)$$

2.  *Transverse wave*[^3]: In this case the motion of particles of media
    is restricted in the plane that is normal to the direction of wave
    propagation (see Fig.
    [2](#fig:Plane-wave-motion){reference-type="ref"
    reference="fig:Plane-wave-motion"}) and the phase velocity $c$ is
    given by Eq.
    [\[eq: trans-wave-speed\]](#eq: trans-wave-speed){reference-type="eqref"
    reference="eq: trans-wave-speed"}.

    $${\mathbf{p}} \cdot {\mathbf{d}} = 0; \quad c = {c_T}$$

The transverse waves further categorized into *SV-waves* and *SH-waves*
depending upon the plane in which particle motion resides. Let's say the
wave is traveling in $\left( x_1, x_2 \right)$ plane then in case of the
*SV-waves* particle motion resides in the same $\left( x_1, x_2 \right)$
plane but normal to the direction of wave propagation. Therefore,
$\mathbf{d}$ can be represented by

$$\label{eq:SV-wave-cond}
{\mathbf{d}} = {{\mathbf{e}}_3} \times {\mathbf{p}}$$

However, in the case of *SH-waves* particles will move in the $x_3$
direction, i.e.,

$$\label{eq:SH-wave-cond}
\mathbf{d} = \mathbf{e}_{3}$$

where ${{\mathbf{e}}_3} = {\left[ {0,0,1} \right]^T}$.

Finally, a complete classification of plane waves can be found in Fig
.[1](#fig:Plane-wave-types){reference-type="ref"
reference="fig:Plane-wave-types"}. For clarity the plane wave
propagation in unbounded elastic domain and the corresponding particle
motions are depicted in Fig.
[2](#fig:Plane-wave-motion){reference-type="ref"
reference="fig:Plane-wave-motion"}.

<figure id="fig:Plane-wave-types">
<img src="./figures/Plane-wave-1" />
<figcaption>Classification of plane waves in unbounded elastic
media</figcaption>
</figure>

<figure id="fig:Plane-wave-motion">
<embed src="./figures/Plane-wave-2.eps" />
<figcaption>Schematic diagram of longitudinal and transverse wave motion
in unbounded elastic media</figcaption>
</figure>

A special case of plane waves is plane harmonic waves; material points
on the plane of constant phase performs harmonic motion. The studies on
plane harmonic waves in a linearly elastic medium are of interest by the
virtue of the applicability of linear superimposition. By the use of
Fourier series, harmonic waves can be used to describe the propagation
of periodic disturbances. For plane harmonic waves Eq.
[\[eq:plane-wave-1\]](#eq:plane-wave-1){reference-type="eqref"
reference="eq:plane-wave-1"} becomes

$$\label{eq:plane-harmonic-1}
{\mathbf{u}} = A{\mathbf{d}}\exp \left[ {ik\left( {{\mathbf{x}} \cdot {\mathbf{p}} - ct} \right)} \right]$$

where $A \in \mathbb{R}$ is the amplitude of the particle motion,
$i = \sqrt{-1}$, $k=\omega c \in \mathbb{R}$ is called the wave number,
$c \in \mathbb{R}$ is the phase velocity, $\omega \in \mathbb{R}$ is the
circular frequency of the harmonic motion of the material points,
$\mathbf{d} \in \mathbb{R}^{n_{sd}}$ is the direction of motion of
particle and $\mathbf{p}$ is direction of propagation of wavefront.
Since the phase velocity does not depend upon the wave number (or the
wavelength), the plane harmonic waves in unbounded homogeneous,
isotropic, linearly elastic media are not dispersive.

Noting Eq.
[\[eq:longi-wave-cond\]](#eq:longi-wave-cond){reference-type="eqref"
reference="eq:longi-wave-cond"} and Eq.
[\[eq:plane-harmonic-1\]](#eq:plane-harmonic-1){reference-type="eqref"
reference="eq:plane-harmonic-1"} the longitudinal plane harmonic wave
traveling in $\left( x_{1}, x_{2}\right)$ plane in the direction
$\mathbf{p}=\left[ \sin\theta, \cos\theta,0 \right]$ can be described as
follows

$$\label{eq:longi-plane-harmonic}
{\mathbf{u}} = A\left[ {\begin{array}{cc}
  {\sin \theta } \\
  {\cos \theta } \\
  0
\end{array}} \right]\exp \left[ {ik\left( {{x_1}\sin {\theta _1} + {x_2}\cos {\theta _2} - {c_L}t} \right)} \right]$$

Similarly, using Eq.
[\[eq:SV-wave-cond\]](#eq:SV-wave-cond){reference-type="eqref"
reference="eq:SV-wave-cond"} in Eq.
[\[eq:plane-harmonic-1\]](#eq:plane-harmonic-1){reference-type="eqref"
reference="eq:plane-harmonic-1"}, the *SV-wave* traveling in
$\left( x_{1}, x_{2}\right)$ plane and in the direction
$\mathbf{p}=\left[ \sin\theta, \cos\theta,0 \right]$ is given by

$$\label{eq:SV-plane-harmonic}
{\mathbf{u}} = A\left[ {\begin{array}{cc}
  { - \cos \theta } \\
  {\sin \theta } \\
  0
\end{array}} \right]\exp \left[ {ik\left( {{x_1}\sin \theta  + {x_2}\cos \theta  - {c_T}t} \right)} \right]$$

The equation of SH-wave propagating in the direction
$\mathbf{p}=\left[ \sin\theta, \cos\theta,0 \right]$ can be obtained by
using Eq. [\[eq:SH-wave-cond\]](#eq:SH-wave-cond){reference-type="eqref"
reference="eq:SH-wave-cond"}

$$\label{eq:SH-plane-harmonic}
{u_3} = A\exp \left[ {ik\left( {{x_1}\sin \theta  + {x_2}\cos \theta  - {c_T}t} \right)} \right]$$

For the displacement given by Eq.
[\[eq:plane-harmonic-1\]](#eq:plane-harmonic-1){reference-type="eqref"
reference="eq:plane-harmonic-1"} the stresses generated in the elastic
body are computed using the Hooke's law (see Eq.
([\[eq:Hooke\]](#eq:Hooke){reference-type="ref" reference="eq:Hooke"} --
[\[eq:Strain-Disp\]](#eq:Strain-Disp){reference-type="ref"
reference="eq:Strain-Disp"}))

$$\label{eq:plane-harmonic-stress}
\sigma  = ikA\left[ {\lambda {\mathbf{d}} \cdot {\mathbf{p}} + \mu \left( {{\mathbf{d}} \otimes {\mathbf{p}} + {\mathbf{p}} \otimes {\mathbf{d}}} \right)} \right]\exp \left[ {ik\left( {{\mathbf{x}} \cdot {\mathbf{p}} - ct} \right)} \right]$$

## Reflection and refraction of plane waves

The presence of different materials significantly affects the systems of
waves propagating through that medium. When waves reach the interface
between two medium with different material properties, part of the wave
is reflected and the part is transmitted through the interface
[@Achenbach1973a]. The ratio of the mechanical impedances of two media
completely determines the nature of the reflection and the transmission
at the interface.

In such system of material discontinuity, system of plane waves ---
reflected and refracted waves at the interface --- can be superposed to
represent an incident wave. For a given incident wave, the amplitude
$A$, the unit propagation vectors $\mathbf{p}$ and the wave numbers $k$
of the reflected and refracted waves are computed by satisfying the
continuity conditions on the displacements and stress at the interface
between two media.

Consider two joined elastic half-spaces in $\left( x_1, x_2 \right)$
plane. Let $x_2 = 0$ be the interface plane between these two media. The
material properties of the medium carrying the incident and reflected
waves (i.e. $x_2<0$) are the Lame elastic constants $\lambda$ and $\mu$,
velocity of longitudinal wave $c_L$, velocity of transverse wave $c_T$
and the mass density $\rho$. Similarly, the material constants of the
medium into which refraction takes place are
$\lambda^{B}, \mu^{B}, c^{B}_{L}, c^{B}_{T}, \rho^{B}$. The superscript
with a number $n=1,2,3,4$ enclosed in parenthesis is used for denoting
the reflected and refracted waves. The number 0 is used for denoting the
incidental wave.

#### *Reflection and refraction of SH-wave*:

Consider the following incidental *SH-wave* traveling in the direction
$\mathbf{p}=\left[ \sin\theta, \cos\theta,0 \right]$ (see Fig.
[3](#fig:SH-wave-reflection){reference-type="ref"
reference="fig:SH-wave-reflection"}).

$$\label{eq:SH-wave-in}
u_3^{\left( 0 \right)} = {A_0}\exp \left[ {i{k_0}\left( {{x_1}\sin {\theta _0} + {x_2}\cos {\theta _0} - {c_T}t} \right)} \right]$$

<figure id="fig:SH-wave-reflection">
<img src="./figures/SH-Wave-Reflection" />
<figcaption>Reflection of SH-wave at the boundary</figcaption>
</figure>

The reflection and refraction[^4] of *SH-wave* at the interface
generates the *SH-wave*. The system of waves should satisfy the
continuity of the displacement and stress at the interface; displacement
(stress) due to the incidental and reflected wave should be equal to the
displacement(stress) due to the refracted wave. Accordingly, the
equations of reflected and refracted wave are given by

$$\label{eq:SH-wave-reflected}
u_3^{\left( 2 \right)} =  {A_2}\exp \left[ {i{k_0}\left( {{x_1}\sin {\theta_0} + {x_2}\cos {\theta_0} - {c_T}t} \right)} \right]$$

$$\label{eq:SH-wave-refracted}
u_3^{\left( 4 \right)} =  {A_4}\exp \left[ {i{k_4}\left( {{x_1}\sin {\theta_4} + {x_2}\cos {\theta_4} - {c_{T}^{B}}t} \right)} \right]$$

where $A_2, A_4, \theta_4, k_4$ are given by following relations

$$\label{eq:ch4-eq-28}
\frac{{{A_2}}}{{{A_0}}} = \frac{{\mu \cos {\theta _0} - {\mu ^B}\left( {{c_T}/c_T^B} \right)\cos {\theta _4}}}{{\mu \cos {\theta _0} + {\mu ^B}\left( {{c_T}/c_T^B} \right)\cos {\theta _4}}}$$

$$\frac{{{A_4}}}{{{A_0}}} = \frac{{2\mu \cos {\theta _0}}}{{\mu \cos {\theta _0} + {\mu ^B}\left( {{c_T}/c_T^B} \right)\cos {\theta _4}}}$$

$${k_4} = \frac{{{c_T}}}{{c_T^B}}{k_0}$$

$$\label{eq:ch4-eq-31}
\sin {\theta _4} = \frac{{c_T^B}}{{{c_T}}}\sin {\theta _0}$$

The inspection of Eq.
([\[eq:ch3-eq-25\]](#eq:ch3-eq-25){reference-type="ref"
reference="eq:ch3-eq-25"}--[\[eq:ch3-eq-31\]](#eq:ch3-eq-31){reference-type="ref"
reference="eq:ch3-eq-31"}) leads to several observations:

1.  The reflected wave is in phase with the incident wave and the
    wave-number and phase-velocity of incident and reflected wave are
    the same.

2.  Refracted wave separates from the incident wave while moving away
    from the normal if $c_{T}^{B}>c_{T}$, and it moves towards the
    normal in case $c_{T}^{B}<c_{T}$.

3.  The wave is completely transmitted (i.e. $A_2=0$) if

    $$\mu \cos {\theta _0} - {\mu ^B}\left( {{c_T}/c_T^B} \right)\cos {\theta _4} = 0$$

    which leads to

    $$\cos {\theta _0} = \sqrt {\frac{{1 - {{\left( {{c_T}/c_T^B} \right)}^2}}}{{1 - {{\left( {\mu /{\mu ^B}} \right)}^2}}}}$$

    Therefore, a combination of angle of incidence and material
    properties is possible for which no *SH-wave* is reflected.

4.  If the half-space $x_{2}>0$ is vacuum (i.e., the interface is
    boundary of elastic half-space) then there will be no reflected
    waves. Further, if the boundary condition condition at $x_2=0$ is
    such that total displacement vanishes then the reflected wave will
    have $180\,^{\circ}$ phase difference with the incident wave.
    Furthermore, if the total stress at the boundary vanishes then the
    reflected wave is in phase with the incident wave.

#### *Reflection and refraction of longitudinal wave*:

As mentioned earlier, in case of longitudinal wave (or *P-wave*)
material points move in the direction of wave propagation (see Eq.
[\[eq:plane-wave-longi\]](#eq:plane-wave-longi){reference-type="ref"
reference="eq:plane-wave-longi"}). The equation of incidental
longitudinal wave traveling in the direction
$\mathbf{p}=\left[ \sin\theta, \cos\theta,0 \right]$ is given by

$$\label{eq:P-wave-in}
{{\mathbf{u}}^{\left( 0 \right)}} = \left[ {\begin{array}{cc}
  {\sin {\theta _0}} \\
  {\cos {\theta _0}} \\
  0
\end{array}} \right]{A_0}\exp \left[ {i{k_0}\left( {{x_1}\sin {\theta _0} + {x_2}\cos {\theta _0} - {c_L}t} \right)} \right]$$

<figure id="fig:P-wave-reflection">
<img src="./figures/P-Wave-Reflection" />
<figcaption>Reflection of longitudinal plane wave</figcaption>
</figure>

When a *P-wave* encounters a material interface two reflected (a
*P-wave* and a *SV-wave*) and two refracted waves (a *P-wave* and a
*SV-wave*) are possible. The motion of system of waves is depicted in
Fig. [4](#fig:P-wave-reflection){reference-type="ref"
reference="fig:P-wave-reflection"}. The equation of reflected *P-wave*
and *SV-wave* is given by Eq.
[\[eq:P-wave-reflected\]](#eq:P-wave-reflected){reference-type="eqref"
reference="eq:P-wave-reflected"} and Eq.
[\[eq:SV-wave-reflected\]](#eq:SV-wave-reflected){reference-type="eqref"
reference="eq:SV-wave-reflected"}, respectively.

$$\label{eq:P-wave-reflected}
{{\mathbf{u}}^{\left( 1 \right)}} = \left[ {\begin{array}{cc}
  {\sin {\theta _1}} \\
  { - \cos {\theta _1}} \\
  0
\end{array}} \right]{A_1}\exp \left[ {i{k_1}\left( {{x_1}\sin {\theta _1} - {x_2}\cos {\theta _1} - {c_L}t} \right)} \right]$$

$$\label{eq:SV-wave-reflected}
{{\mathbf{u}}^{\left( 2 \right)}} = \left[ {\begin{array}{cc}
  {\cos {\theta _2}} \\
  {\sin {\theta _2}} \\
  {0}
\end{array}} \right]{A_2}\exp \left[ {i{k_2}\left( {{x_1}\sin {\theta _2} - {x_2}\cos {\theta _2} - {c_T}t} \right)} \right]$$

The equation of refracted *P-wave* and *SV-wave* is given by Eq.
[\[eq:P-wave-refracted\]](#eq:P-wave-refracted){reference-type="eqref"
reference="eq:P-wave-refracted"} and Eq.
[\[eq:SV-wave-refracted\]](#eq:SV-wave-refracted){reference-type="eqref"
reference="eq:SV-wave-refracted"} respectively.

$$\label{eq:P-wave-refracted}
{{\mathbf{u}}^{\left( 3 \right)}} = \left[ {\begin{array}{cc}
  {\sin {\theta _3}} \\
  {\cos {\theta _3}} \\
  0
\end{array}} \right]{A_3}\exp \left[ {i{k_3}\left( {{x_1}\sin {\theta _3} + {x_2}\cos {\theta _3} - c_L^Bt} \right)} \right]$$

$$\label{eq:SV-wave-refracted}
{{\mathbf{u}}^{\left( 4 \right)}} = \left[ {\begin{array}{cc}
  { - \cos {\theta _4}} \\
  {\sin {\theta _4}} \\
  0
\end{array}} \right]{A_4}\exp \left[ {i{k_4}\left( {{x_1}\sin {\theta _4} + {x_2}\cos {\theta _4} - c_T^Bt} \right)} \right]$$

This system of waves should satisfy four continuity equations at the
interface; one for each $u_{1}, u_{2}, \sigma_{12}, \sigma_{22}$. In
other words, total displacement (total stress) due to system of incident
and reflected waves must be equal to the total displacement (total
stress) due to the system of refracted waves. Accordingly, one can
obtain following results for the wave-numbers and direction of wave
propagation.

$${k_1} = {k_0};\quad \,{\theta _1} = {\theta _0}$$

$${k_2} = {k_0}\frac{{{c_L}}}{{{c_T}}};\quad \,\sin {\theta _2} = \frac{{{c_T}}}{{{c_L}}}\sin {\theta _0}$$

$${k_3} = {k_0}\frac{{{c_L}}}{{c_L^B}};\quad \,\sin {\theta _3} = \frac{{c_L^B}}{{{c_L}}}\sin {\theta _0}$$

$${k_4} = {k_0}\frac{{{c_L}}}{{c_T^B}};\quad \,\sin {\theta _4} = \frac{{c_T^B}}{{{c_L}}}\sin {\theta _0}$$

Subsequently, the aforementioned four continuity equations results in
four linear equations for the amplitudes $A_{1}, A_{2}, A_{3}$ and
$A_{4}$. In matrix form the system can be represented by

$$\mathbf{T}\mathbf{A}=\mathbf{A}_{0}$$

where

$${\mathbf{T}} = \left[ {\begin{array}{cc}
  { - \sin {\theta _1}}&{ - \cos {\theta _2}}&{\sin {\theta _3}}&{ - \cos {\theta _4}} \\
  {\cos {\theta _1}}&{ - \sin {\theta _2}}&{\cos {\theta _3}}&{\sin {\theta _4}} \\
  {\sin 2{\theta _1}}&{\frac{{{c_L}}}{{{c_T}}}\cos 2{\theta _2}}&{\frac{{{\mu ^B}}}{\mu }\frac{{{c_L}}}{{c_L^B}}\sin 2{\theta _3}}&{ - \frac{{{\mu ^B}}}{\mu }\frac{{{c_L}}}{{c_T^B}}\cos 2{\theta _4}} \\
  { - {{\left( {\frac{{{c_L}}}{{{c_T}}}} \right)}^2}\cos 2{\theta _2}}&{\frac{{{c_L}}}{{{c_T}}}\sin 2{\theta _2}}&{\frac{{{\mu ^B}}}{\mu }\frac{{{c_L}{c^B}_L}}{{{{\left( {c_T^B} \right)}^2}}}\cos 2{\theta _4}}&{\frac{{{\mu ^B}}}{\mu }\frac{{{c_L}}}{{c_T^B}}\sin 2{\theta _4}}
\end{array}} \right]$$

$${\mathbf{A}} = {\left[ {\begin{array}{cc}
  {{A_1}}&{{A_2}}&{{A_3}}&{{A_4}}
\end{array}} \right]^T}$$

$${{\mathbf{A}}_0} = {A_0}\left[ {\begin{array}{cc}
  {\sin {\theta _0}}&{\cos {\theta _0}}&{\sin 2{\theta _0}}&{{{\left( {\frac{{{c_L}}}{{{c_T}}}} \right)}^2}\cos 2{\theta _2}}
\end{array}} \right]^{T}$$

In case of normal incident of *P-wave* (i.e., $\theta_{0}=0$) the
reflected and refracted waves contains only a *P-wave*, i.e.

$$A_{2}=A_{4}=0$$

and

$$\label{eq:ch4-eq-44}
\frac{{{A_1}}}{{{A_0}}} = \frac{{{\rho ^B}c_L^B - \rho {c_L}}}{{{\rho ^B}c_L^B + \rho {c_L}}}$$

$$\label{eq:ch4-eq-45}
\frac{{{A_3}}}{{{A_0}}} = \frac{{2\rho {c_L}}}{{{\rho ^B}c_L^B + \rho {c_L}}}$$

From Eq. ([\[eq:ch4-eq-44\]](#eq:ch4-eq-44){reference-type="ref"
reference="eq:ch4-eq-44"}--[\[eq:ch4-eq-45\]](#eq:ch4-eq-45){reference-type="ref"
reference="eq:ch4-eq-45"}) it is clear that in case of normal incidence
of *P-wave* the amplitude of reflected and refracted *P-wave* depends
entirely on the mechanical impedance, the product of mass density and
phase speed, of the two media. If both media have the same mechanical
impedance then there will be no reflected wave. Further, if upper
half-space has more impedance then lower half-space then the reflected
*P-wave* will be out of phase with the incident wave.

In case upper half-space is vacuum -- the interface becomes the boundary
of elastic half-space -- then there will be no refracted waves and
system of reflected wave will consist a P-wave and a SV-wave. The nature
of reflected waves then depends upon the boundary conditions. If the
total displacement at the boundary vanishes then the reflected P-wave
and SV-wave are given by Eq.
[\[eq:P-wave-reflected\]](#eq:P-wave-reflected){reference-type="eqref"
reference="eq:P-wave-reflected"} and Eq.
[\[eq:SV-wave-reflected\]](#eq:SV-wave-reflected){reference-type="eqref"
reference="eq:SV-wave-reflected"}, respectively, with amplitudes,

$$\label{eq:P-wave-clamp-A1}
\frac{{{A_1}}}{{{A_0}}} = \frac{{\cos \left( {{\theta _0} + {\theta _2}} \right)}}{{\cos \left( {{\theta _0} - {\theta _2}} \right)}}$$

$$\label{eq:P-wave-clamp-A2}
\frac{{{A_2}}}{{{A_0}}} = \frac{{\sin 2{\theta _0}}}{{\cos \left( {{\theta _0} - {\theta _2}} \right)}}$$

If the total stress at the boundary vanishes (i.e., stress free surface)
then the amplitude ratio becomes

$$\label{eq:P-wave-free-A1}
\frac{{{A_1}}}{{{A_0}}} = \frac{{\sin 2{\theta _0}\sin 2{\theta _2} - {\kappa ^2}{{\cos }^2}{\theta _2}}}{{\sin 2{\theta _0}\sin 2{\theta _2} + {\kappa ^2}{{\cos }^2}2{\theta _2}}}$$

$$\label{eq:P-wave-free-A2}
\frac{{{A_2}}}{{{A_0}}} = \frac{{2\kappa \sin 2{\theta _0}\cos 2{\theta _2}}}{{\sin 2{\theta _0}\sin 2{\theta _2} + {\kappa ^2}{{\cos }^2}2{\theta _2}}}$$

From Eq.
([\[eq:P-wave-free-A1\]](#eq:P-wave-free-A1){reference-type="ref"
reference="eq:P-wave-free-A1"}--[\[eq:P-wave-free-A2\]](#eq:P-wave-free-A2){reference-type="ref"
reference="eq:P-wave-free-A2"}) it is observed that normal incident,
$\theta_{0}=0$, and grazing incident, $\theta_{0}=\pi/2$, of *P-wave*
generates only a reflected *P-wave* as $A_2=0$. In former case the
reflected *P-wave* is in the phase, and in later case reflected *P-wave*
is out of phase with incident wave.

## Surface waves {#sec:ch4_sec2-5}

So far we have discussed the homogeneous plane harmonic waves where the
amplitude of wave motion remains constant in the plane of constant
phase. There exists another type of wave motion for which the amplitude
changes in the plane of constant phase. Consequently, one can define the
plane of constant amplitude. It turns out that the two planes are normal
to each other. Moreover, the plane of constant phase moves in the
direction of wave propagation, therefore, the amplitude remains constant
in the wave propagation direction.

Surface waves are inhomogeneous plane waves for which amplitude of
disturbance exponentially decays as one moves away from the surface.
However, the amplitude of motion remains constant in the wave direction
of wave propagation. From an earthquake engineering viewpoint two type
of surface waves are of primary importance; *Rayleigh wave* and *Love
wave*.

Consider the in-plane motion of plane waves traveling in $x_1$-direction
in a homogeneous elastic half-space ($x_{2} \le 0$) with free surface at
$x_2 = 0$. The Motion of particles as the wave passes by can be
described as,

$$\label{eq:ch4-eq-50}
{u_1} = \left[ {{A_1}\exp \left( {{b_1}{x_2}} \right) + {A_2}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]$$

$$\label{eq:ch4-eq-51}
{u_2} = \left[ { - {A_1}\frac{{{b_1}}}{{i{k_R}}}\exp \left( {{b_1}{x_2}} \right) + {A_2}\frac{{i{k_R}}}{{{b_2}}}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]$$

where $A_{1},A_{2}$ are constants to be determined, $k_{R}$ and $c_R$
are the wave number and phase-velocity of the *Rayleigh wave*,
respectively. $b_{1}$ and $b_{2}$ are given by

$${b_1} = {k_R}{\left( {1 - \frac{{c_R^2}}{{c_L^2}}} \right)^{1/2}}$$

$${b_2} = {k_R}{\left( {1 - \frac{{c_R^2}}{{c_T^2}}} \right)^{1/2}}$$

Note that for $b_{1}, b_{2}$ are real valued if $c_{R}<c_{T}<c_{L}$. The
mathematical expressions for $A_{1}, A_{2}$ and $c_{R}$ are obtained by
satisfying the stress free conditions (i.e. $\sigma_{22}=\sigma_{12}=0$)
at $x_2=0$. Subsequently, one can obtain the following relationship
between the constants $A_{1}$ and $A_{2}$

$$\label{eq:ch4-eq-54}
{A_2} =  - \frac{{2{b_1}{b_2}}}{{\left( {k_R^2 + b_2^2} \right)}}{A_1}$$

The approximate value for $c_R$ can be given by [@Achenbach1973a]

$$\label{eq:ch4-eq-55}
{c_R} = \frac{{0.862 + 1.14\nu }}{{1 + \nu }}{c_T}$$

where $\nu$ denotes the Poission's ratio of the elastic half-space.

Using Eq. ([\[eq:ch4-eq-54\]](#eq:ch4-eq-54){reference-type="ref"
reference="eq:ch4-eq-54"}--[\[eq:ch4-eq-55\]](#eq:ch4-eq-55){reference-type="ref"
reference="eq:ch4-eq-55"}) in Eq.
([\[eq:ch4-eq-50\]](#eq:ch4-eq-50){reference-type="ref"
reference="eq:ch4-eq-50"}--[\[eq:ch4-eq-51\]](#eq:ch4-eq-51){reference-type="ref"
reference="eq:ch4-eq-51"}) following expression for displacement
components can be obtained.

$$\label{eq:ch4-eq-56}
{u_1} = {A_1}\left[ {\exp \left( {{b_1}{x_2}} \right) - \frac{{2{b_1}{b_2}}}{{k_R^2 + b_2^2}}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]$$

$$\label{eq:ch4-eq-57}
{u_2} = i{A_1}\left[ {\frac{{{b_1}}}{{{k_R}}}\exp \left( {{b_1}{x_2}} \right) - \frac{{2{b_1}{k_R}}}{{k_R^2 + b_2^2}}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]$$

From Eq. ([\[eq:ch4-eq-56\]](#eq:ch4-eq-56){reference-type="ref"
reference="eq:ch4-eq-56"}--[\[eq:ch4-eq-57\]](#eq:ch4-eq-57){reference-type="ref"
reference="eq:ch4-eq-57"}) it can be observed that the horizontal and
vertical displacement components, $u_{1}, u_{2}$, have a $90\,^{\circ}$
phase difference; when $u_1$ is maximum then $u_2$ is zero, and vice
versa. Due to $90\,^{\circ}$ phase difference in $u_{1}, u_{2}$ the
trajectories of the particles are ellipses, and particles at the free
surfaces moves counter-clockwise when wave travels in positive $x_1$
direction. At depth $x_{2}\approx 0.2 \Lambda$ the direction of rotation
changes (Here $\Lambda$ denotes the wavelength of *Rayleigh wave*). It
can be shown that at the free surface the normal displacement is about
1.5 times the tangential displacement.

It should be noted that in a homogeneous elastic half-space only
*Rayleigh wave* and body waves can exist. However, *Love waves* can
arise if soil layering is present. *Love waves* typically develop in
shallow surface soil layers overlying layers of stiffer materials
properties. They basically consists of *SH-waves* that are trapped by
multiple reflection within the surface layer. Exactly like *SH-waves*
*Love waves* propagates in the out-of-plane direction and they have no
vertical components of particle motion.


