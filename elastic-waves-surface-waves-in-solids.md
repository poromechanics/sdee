# Surface waves in elastic solids {#sec-surface-waves-in-solid}

So far we have discussed the homogeneous plane harmonic waves where the amplitude of wave motion remains constant in the plane of constant phase. There exists another type of wave motion for which the amplitude changes in the plane of constant phase. Consequently, one can define the plane of constant amplitude. It turns out that the two planes are normal to each other. Moreover, the plane of constant phase moves in the direction of wave propagation, therefore, the amplitude remains constant in the wave propagation direction.

Surface waves are inhomogeneous plane waves for which amplitude of disturbance exponentially decays as one moves away from the surface. However, the amplitude of motion remains constant in the wave direction of wave propagation. From an earthquake engineering viewpoint two type of surface waves are of primary importance; *Rayleigh wave* and *Love wave*.

Consider the in-plane motion of plane waves traveling in $x_1$-direction in a homogeneous elastic half-space ($x_{2} \le 0$) with free surface at $x_2 = 0$. The Motion of particles as the wave passes by can be described as,

$$
{u_1} = \left[ {{A_1}\exp \left( {{b_1}{x_2}} \right) + {A_2}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]
$$ {#eq-ch4-50}

$$
{u_2} = \left[ { - {A_1}\frac{{{b_1}}}{{i{k_R}}}\exp \left( {{b_1}{x_2}} \right) + {A_2}\frac{{i{k_R}}}{{{b_2}}}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]
$$ {#eq-ch4-51}

where $A_{1},A_{2}$ are constants to be determined, $k_{R}$ and $c_R$ are the wave number and phase-velocity of the *Rayleigh wave*, respectively. $b_{1}$ and $b_{2}$ are given by

$${b_1} = {k_R}{\left( {1 - \frac{{c_R^2}}{{c_L^2}}} \right)^{1/2}}$$

$${b_2} = {k_R}{\left( {1 - \frac{{c_R^2}}{{c_T^2}}} \right)^{1/2}}$$

Note that for $b_{1}, b_{2}$ are real valued if $c_{R}<c_{T}<c_{L}$. The mathematical expressions for $A_{1}, A_{2}$ and $c_{R}$ are obtained by satisfying the stress free conditions (i.e. $\sigma_{22}=\sigma_{12}=0$) at $x_2=0$. Subsequently, one can obtain the following relationship between the constants $A_{1}$ and $A_{2}$

$$
{A_2} =  - \frac{{2{b_1}{b_2}}}{{\left( {k_R^2 + b_2^2} \right)}}{A_1}
$$ {#eq-ch4-54}

The approximate value for $c_R$ can be given by [@Achenbach1973a]

$$
{c_R} = \frac{{0.862 + 1.14\nu }}{{1 + \nu }}{c_T}
$$ {#eq-ch4-55}

where $\nu$ denotes the Poisson's ratio of the elastic half-space.

Using @eq-ch4-54 -- @eq-ch4-55 in @eq-ch4-50 and @eq-ch4-51 following expression for displacement components can be obtained.

$$
{u_1} = {A_1}\left[ {\exp \left( {{b_1}{x_2}} \right) - \frac{{2{b_1}{b_2}}}{{k_R^2 + b_2^2}}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]
$$ {#eq-ch4-56}

$$
{u_2} = i{A_1}\left[ {\frac{{{b_1}}}{{{k_R}}}\exp \left( {{b_1}{x_2}} \right) - \frac{{2{b_1}{k_R}}}{{k_R^2 + b_2^2}}\exp \left( {{b_2}{x_2}} \right)} \right]\exp \left[ {i{k_R}\left( {{x_1} - {c_R}t} \right)} \right]
$$ {#eq-ch4-57}

From @eq-ch4-56 -- @eq-ch4-57, it can be observed that the horizontal and vertical displacement components, $u_{1}, u_{2}$, have a $90\,^{\circ}$ phase difference; when $u_1$ is maximum then $u_2$ is zero, and vice versa. Due to $90\,^{\circ}$ phase difference in $u_{1}, u_{2}$ the trajectories of the particles are ellipses, and particles at the free surfaces moves counter-clockwise when wave travels in positive $x_1$ direction. At depth $x_{2}\approx 0.2 \Lambda$ the direction of rotation changes (Here $\Lambda$ denotes the wavelength of *Rayleigh wave*). It can be shown that at the free surface the normal displacement is about 1.5 times the tangential displacement.

It should be noted that in a homogeneous elastic half-space only *Rayleigh wave* and body waves can exist. However, *Love waves* can arise if soil layering is present. *Love waves* typically develop in shallow surface soil layers overlying layers of stiffer materials properties. They basically consists of *SH-waves* that are trapped by multiple reflection within the surface layer. Exactly like *SH-waves* *Love waves* propagates in the out-of-plane direction and they have no vertical components of particle motion.
