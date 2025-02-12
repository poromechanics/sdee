# V-STFEM for Elastodynamics {#sec-elastodynamics-vstfem}

## Continuum theory of elastodynamics

Let $\Omega\subset\mathbb{R}^{n_{sd}}$ be an open and bounded region occupied by an elastic body at time $t$, where $n_{sd}$ is the number of spatial dimensions. The boundary of $\Omega$ is denoted by $\Gamma$. Let $\bar{\Omega}=\Omega\cup\Gamma$ denote the closure of $\Omega$. Further, let indices $i,j,k$ and $l$ take values from $1,\cdots,n_{sd}$, and the Einstein summation convention applies to the repeated indices $i,j,k,$ and $l$ only. Furthermore, consider the nonoverlapping partitions, $\Gamma_{i}^{g}$ and $\Gamma_{i}^{h}$, of the boundary $\Gamma$ such that

$$\begin{aligned}
\Gamma & =\Gamma_{i}^{g}\cup\Gamma_{i}^{h}, & \Gamma_{i}^{g}\cap\Gamma_{i}^{h} & =\phi, & i & =1,\cdots,n_{sd}
\end{aligned}
$$

The displacement, velocity, and stress field are denoted by $\mathbf{u}$, $mathbf{v}$, and $\sigma$, respectively. The infinitesimal strain tensor, $\varepsilon$, and stretching tensor, $\mathbf{d}$, are given by

$$
\begin{aligned}
\varepsilon:={\mathbf{e}}\left({\mathbf{u}}\right) & =\frac{1}{2}\left({{\mathbf{u}}\otimes{\nabla_{x}}}\right)+\frac{1}{2}\left({{\nabla_{x}}\otimes{\mathbf{u}}}\right), & {\varepsilon_{ij}} & ={e_{ij}}\left({\mathbf{u}}\right)=\frac{1}{2}\left({\frac{{\partial{u_{i}}}}{{\partial{x_{j}}}}+\frac{{\partial{u_{j}}}}{{\partial{x_{i}}}}}\right)
\end{aligned}
$${#eq-ch3-119}

$$
\begin{aligned}
{\mathbf{d}}:={\mathbf{e}}\left({\mathbf{v}}\right) & =\frac{1}{2}\left({{\mathbf{v}}\otimes{\nabla_{x}}}\right)+\frac{1}{2}\left({{\nabla_{x}}\otimes{\mathbf{v}}}\right), & {d_{ij}} & ={e_{ij}}\left({\mathbf{v}}\right)=\frac{1}{2}\left({\frac{{\partial{v_{i}}}}{{\partial{x_{j}}}}+\frac{{\partial{v_{j}}}}{{\partial{x_{i}}}}}\right)
\end{aligned}
$${#eq-ch3-120}

in which $\varepsilon_{ij}$, $d_{ij}$, $u_{i}$ and $v_{i}$ denote the Cartesian components of $\varepsilon$, $\mathbf{d}$, $\mathbf{u}$, and $\mathbf{v}$, respectively.

In the small strain framework, the linear elastic constitutive relationship is described as,

$$
\begin{aligned} \dot{\sigma}_{ij} & =C_{ijkl}d_{kl}, & \sigma_{ij} & =C_{ijkl}\varepsilon_{kl}
\end{aligned}
$${#eq-ch3-121}

where $\dot{\sigma}_{ij}$ denotes the first-order time derivative of the stress field, and $C_{ijkl}$ is the fourth-order elasticity tensor. In the case of an isotropic material, the elasticity tensor is expressed using the Lame parameters, $\lambda,\mu$:

$$
C_{ijkl}:=\lambda\delta_{ij}\delta_{kl}+2\mu\left(\frac{\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}}{2}\right)
$${#eq-ch3-122}

in which $\delta_{ij}$ represents the Kronecker-delta function. If $i=j$ then $\delta_{ij}=1$, otherwise $\delta_{ij}=0$.

The strong form of the initial-boundary value problem of elastodynamics
can be stated as: given the functions

$$b_{i}:\Omega\times\left[0,T\right]\rightarrow\mathbb{R},$$

$$g_{i}:\Gamma_{i}^{g}\times\left[0,T\right]\rightarrow\mathbb{R},$$

$$f_{i}^{s}:\Gamma_{i}^{h}\times\left[0,T\right]\rightarrow\mathbb{R},$$

$${u_{0}}_{i}:\Omega\rightarrow\mathbb{R},$$

$${v_{0}}_{i}:\Omega\rightarrow\mathbb{R},$$

$$\rho:\Omega\rightarrow\mathbb{R},$$

find $u_{i}:\bar{\Omega}\times\left[0,T\right]\rightarrow\mathbb{R}$, such that

$$
\begin{aligned}
\rho\frac{{{\partial^{2}}{u_{i}}}}{{\partial{t^{2}}}}-\frac{{\partial{\sigma_{ij}}}}{{\partial{x_{j}}}}-\rho{b_{i}} & =0, & \forall(\mathbf{x},t)\in\Omega\times(0,T)
\end{aligned}
$${#eq-ch3-123}

$$
\begin{aligned}
{u_{i}} & ={g_{i}}, & \forall(\mathbf{x},t)\in\Gamma_{i}^{g}\times(0,T)
\end{aligned}
$${#eq-ch3-124}

$$
\begin{aligned}
{\sigma_{ij}}{n_{j}} & =f_{i}^{s}, & \forall(\mathbf{x},t)\in\Gamma_{i}^{h}\times(0,T)
\end{aligned}
$${#eq-ch3-125}

$$
\begin{aligned}
{u_{i}}\left({{\mathbf{x}},0}\right) & ={u_{0i}}, & \forall\mathbf{x}\in\Omega
\end{aligned}
$${#eq-ch3-126}

$$\begin{aligned}
\frac{{\partial{u_{i}}\left({{\mathbf{x}},0}\right)}}{{\partial t}} & ={v_{0i}}, & \forall\mathbf{x}\in\Omega
\end{aligned}
$${#eq-ch3-127}

where $\rho$ is the mass density of the elastic body, $b_{i}$ is the body force density, $g_{i}$ is the prescribed displacement on the Dirichlet-boundary $\Gamma_{i}^{g}$, $f_{i}^{s}$ is the prescribed traction on the Neumann-boundary $\Gamma_{i}^{h}$, $u_{0i}$ is the initial value of the displacement field, and $v_{0i}$ is the initial value of the velocity field.

## Space-time notations

Let $\Omega_{h}$, the set of finite spatial elements $\Omega_{e},e=1,\cdots,n_{el}$, be the discretization of spatial domain $\Omega$, where $n_{el}$ is the total number of spatial elements in $\Omega_{h}$. Furthermore, consider a non-uniform subdivision for the time domain $\left[0,T\right]$,

$$
0=t_{0}<t_{1}<\cdots<t_{N}=T
$$

with

$$
I_{n}=(t_{n},t_{n}+1),\quad\Delta t_{n}=t_{n+1}-t_{n},\quad\Delta t=\mathop{\max}\limits_{0\leqslant n\leqslant N-1}\Delta{t_{n}}.
$$

The space-time slab $Q_{n}$ and the space-time finite element $Q_{n,e}$ are given by following expressions,

$$
\begin{aligned} Q_{n} & :=\Omega_{h}\times I_{n}, & Q_{n,e} & :=\Omega_{e},\quad e=1,\cdots,n_{el}
\end{aligned}
$$

Accordingly,

$$
{\mathbb{Q}_{h}}:=\bigcup\limits_{n=0}^{N-1}{{Q_{n}}}
$$

denotes the discretization of the entire space-time domain.

In TDG/ST/FEM, the unknown fields remains continuous in the spatial domain $\Omega_{h}$, and discontinuity in time occurs at times that belong to the finite set,

$$
{D_{t}}:=\left\{ {{t_{0}},{t_{1}},\ldots,{t_{N}}}\right\}.
$$

Therefore, at the discrete times $t\in D_{t}$ the solution have two values, and the jump discontinuity in time for some unknown scalar field $q(\mathbf{x},t)$ at time $t_{n}\in D_{t}$ is given by @eq-ch3-128

$$
{\left[\kern-0.15em \left[{q\left({\mathbf{x}}\right)}\right]\kern-0.15em \right]_{n}}=q_{n}^{+}\left({\mathbf{x}}\right)-q_{n}^{-}\left({\mathbf{x}}\right)
$${#eq-ch3-128}

where

$$
\begin{aligned}
q_{n}^{+}\left({\mathbf{x}}\right) & =\mathop{\lim}\limits_{\varepsilon\to0}q\left({{\mathbf{x}},{t_{n}}+\varepsilon}\right), & q_{n}^{-}\left({\mathbf{x}}\right) & =\mathop{\lim}\limits_{\varepsilon\to0}q\left({{\mathbf{x}},{t_{n}}-\varepsilon}\right)
\end{aligned}
$${#eq-ch3-129}

denote the right and left limits of the unknown field $q(\mathbf{x},t)$ at time $t=t_{n}$, respectively.

Let us now consider ${\wp_{l}}\left({Q_{n,e}^ {}}\right)$, the collection of all polynomials defined on $Q_{n,e}$ with a total degree of no more than $l$, and $C^{0}(\star)$, the space of piecewise continuous functions defined on domain $(\star)$. Consider also the following collection of functions:

$$
\Im_{l}^{h}:=\left\{ {\left.{{{\mathbf{u}}^{h}}}\right|{{\mathbf{u}}^{h}}\in{C^{0}}{{\left({\bigcup\limits_{n=0}^{N-1}{{Q_{n}}}}\right)}^{{n_{sd}}}},\left.{{{\mathbf{u}}^{h}}}\right|{Q_{n,e}}\in{{\left({{\wp_{l}}\left({{Q_{n,e}}}\right)}\right)}^{{n_{sd}}}}}\right\}
$${#eq-ch3-130}

where ${\left.{{{\mathbf{u}}^{h}}}\right|{Q_{n,e}}}$ is the restriction of $\mathbf{u}^{h}$ to $Q_{n,e}$. Lastly, the space of the test functions for the TDG/ST/FEM is given as

$$
{V^{h}}:=\left\{ {\left.{{{\mathbf{u}}^{h}}}\right|{{\mathbf{u}}^{h}}\in\Im_{l}^{h},u_{i}^{h}=0,\forall\left({{\mathbf{x}},t}\right)\in\Gamma_{i}^{g}\times{I_{n}},{\text{for }}i=1,\cdots,{n_{sd}}}\right\}
$${#eq-ch3-131}

In what follows, a general introduction to the two-field TDG/ST/FEM is provided, then the weak-form for the v-ST/FEM is derived by using the two-field formulation.

## Two-field TDG/ST/FEM

In the two-field formulation, both the displacement field and velocity field are taken as primary unknowns. Accordingly, the weak form should satisfy the following conditions in the weak sense:

1.  Balance of the linear momentum, @eq-ch3-123
2.  Essential and natural boundary conditions, @eq-ch3-123 and @eq-ch3-125
3.  Traction continuity in space;
4.  Continuity of the velocity field in time;

$$
{\left[\kern-0.15em \left[{{{\mathbf{v}}^{h}}\left({\mathbf{x}}\right)}\right]\kern-0.15em \right]_{n}}=0
$$

5.  Continuity of the displacement field in time;

$$
{\left[\kern-0.15em \left[{{{\mathbf{u}}^{h}}\left({\mathbf{x}}\right)}\right]\kern-0.15em \right]_{n}}=0
$$

6.  Displacement-velocity compatibility condition;

$$
\frac{{\partial{{\mathbf{u}}^{h}}}}{{\partial t}}-{{\mathbf{v}}^{h}}=0
$$

Here, Conditions (i)--(iii) are always satisfied as the weak form is derived by employing the Galerkin method [@Hughes2012]. Conditions (iv) and (v) necessary due to the time-discontinuous approximation of the displacement and velocity field, respectively. Condition (vi) is due to the independent approximation of the displacement field and velocity field.

The displacement-velocity two-field weak-form for TDG/ST/FEM in its original form as presented by @Hughes1988 is described as follows.

:::{.callout-note title="Weak form"}

Find $\mathbf{u}\in S_{u}^{h}$ and $\mathbf{v}\in S_{v}^{h}$ such that for all $\delta\mathbf{u}\in V^{h}$ and $\delta\mathbf{v}\in V^{h}$, and for all $n=1,\cdots,N-1$ @eq-ch3-132 holds.

$$
\begin{split} & \int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\delta{v_{i}}\rho\frac{{\partial{v_{i}}}}{{\partial t}}d\Omega}dt}+\int_{{\Omega_{n}}}^ {}{\delta{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)\rho{{\left[\kern-0.15em \left[{{v_{i}}\left({\mathbf{x}}\right)}\right]\kern-0.15em \right]}_{n}}d\Omega}\\
 & +\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\frac{{\partial\delta v_{i}^ {}}}{{\partial{x_{j}}}}{\sigma_{ij}}\left({{\mathbf{x}},t}\right)d\Omega}dt}-\int_{{I_{n}}}^ {}{\int_{\Gamma_{i}^{h}}^ {}{\delta{v_{i}}f_{i}^{s}ds}dt}\\
 & -\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\delta{v_{i}}\rho{b_{i}}d\Omega}dt}+\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\frac{{\partial\delta{u_{i}}}}{{\partial{x_{j}}}}{C_{ijkl}}\frac{\partial}{{\partial{x_{l}}}}\left({\frac{{\partial{u_{k}}}}{{\partial t}}-{v_{k}}}\right)d\Omega}dt}\\
 & +\int_{{\Omega_{h}}}^ {}{\frac{{\partial\delta{u_{i}}\left({{\mathbf{x}},{t_{n}}}\right)}}{{\partial{x_{j}}}}{C_{ijkl}}\frac{{\partial{{\left[\kern-0.15em \left[{{u_{k}}\left({\mathbf{x}}\right)}\right]\kern-0.15em \right]}_{n}}}}{{\partial{x_{l}}}}d\Omega}=0
\end{split}
$${#eq-ch3-132}

where $S_{u}^{h}$ and $S_{v}^{h}$ denote the collections of trial functions for the displacement field and velocity field, respectively, and given by

$$
S_{u}^{h}:=\left\{ {\left.{\mathbf{u}}\right|{\mathbf{u}}\in\Im_{l}^{h},{u_{i}}={g_{i}},\forall\left({{\mathbf{x}},t}\right)\in\Gamma_{i}^{g}\times{I_{n}},i=1,\cdots,{n_{sd}}}\right\}
$${#eq-ch3-133}

$$
S_{v}^{h}:=\left\{ {\left.{\mathbf{v}}\right|{\mathbf{v}}\in\Im_{l}^{h},{v_{i}}=\frac{{\partial{g_{i}}}}{{\partial t}},\forall\left({{\mathbf{x}},t}\right)\in\Gamma_{i}^{g}\times{I_{n}},i=1,\cdots,{n_{sd}}}\right\}
$${#eq-ch3-134}
:::

It is noteworthy that in the above-mentioned weak-form, Conditions (v) and (vi) are satisfied by using a strain-energy norm and Condition (iv) is satisfied by using a kinetic energy norm.

## Velocity based TDG/ST/FEM: v-ST/FEM

In the velocity based time-discontinuous space-time finite element method (v-ST/FEM) displacement-velocity compatibility condition is strongly enforced by computing the displacement field from the consistent time-integration of the velocity field. Accordingly, displacement field remains continuous in time while velocity field still remains discontinuous at the discrete times.

The weak-form for the v-ST/FEM can be stated as:

:::{.callout-note title="Weak form"}
Find $\mathbf{v}\in S_{v}^{h}$ such that for all $\delta\mathbf{v}\in V^{h}$, and for all $n=1,\cdots,N-1$, @eq-ch3-135 holds.

$$
\begin{split} & \int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\delta{v_{i}}\rho\frac{{\partial{v_{i}}}}{{\partial t}}d\Omega}dt}+\int_{{\Omega_{n}}}^ {}{\delta{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)\rho{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)d\Omega}\\
 & +\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\frac{{\partial\delta v_{i}^ {}}}{{\partial{x_{j}}}}{\sigma_{ij}}\left({{\mathbf{x}},t}\right)d\Omega}dt}-\int_{{I_{n}}}^ {}{\int_{\Gamma_{i}^{h}}^ {}{\delta{v_{i}}f_{i}^{s}ds}dt}\\
 & -\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\delta{v_{i}}\rho{b_{i}}d\Omega}dt}-\int_{{\Omega_{n}}}^ {}{\delta{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)\rho{v_{i}}\left({{\mathbf{x}},t_{n}^{-}}\right)d\Omega}=0
\end{split}
$${#eq-ch3-135}

:::

The displacement field is computed by

$$
{\mathbf{u}}\left({{\mathbf{x}},t}\right)={\mathbf{u}}\left({{\mathbf{x}},{t_{n}}}\right)+\int_{{t_{n}}}^{t}{{\mathbf{v}\left({{\mathbf{x}},{\tau}}\right)}d\tau},\qquad\forall\left({{\mathbf{x}},t}\right)\in{Q_{n}},\quad\forall{\mathbf{v}}\in S_{v}^{h}
$$ {#eq-ch3-136}

In the case of hyperelastic material law, the stress is computed by first computing the displacement field (see  @eq-ch3-119 and @eq-ch3-121). Alternatively, if the constitutive relationship is given in the rate form (e.g., hypoelastic form) then stress field can be obtained by time integration of the rate form. In latter case, computation of displacement field can be avoided. In the realm of small strain theory, however, these two techniques of stress recovery are equivalent, and the stress field in a space-time slab takes the form

$$
{\sigma_{ij}}\left({{\mathbf{x}},t}\right)={\sigma_{ij}}\left({{\mathbf{x}},{t_{n}}}\right)+{C_{ijkl}}{\psi_{kl}}\left({{\mathbf{v}},t}\right),\qquad\forall\left({{\mathbf{x}},t}\right)\in{Q_{n}},\forall v\in S_{v}^{h}
$$ {#eq-ch3-137}

where

$$
{\psi_{ij}}\left({{\mathbf{v}},t}\right)=\int_{{t_{n}}}^{t}{\frac{{\partial{v_{i}}}}{{\partial{x_{j}}}}d\tau,}\qquad\forall t\in{I_{n}},\forall{\mathbf{v}}\in S_{v}^{h}
$$ {#eq-ch3-138}

Subsequently, by substituting the expression for the stress given in @eq-ch3-137 into the weak-form, Eq. @eq-ch3-135, the following weak-form for v-ST/FEM is obtained;

:::{.callout-note title="Weak form"}
Find $\mathbf{v}\in S_{v}^{h}$ such that for all
$\delta\mathbf{v}\in V^{h}$, and for all $n=1,\cdots,N-1$, Eq.
[\[eq-ch3-139\]](#eq-ch3-139){reference-type="eqref"
reference="eq-ch3-139"} holds.
$$\begin{split} & \int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\delta{v_{i}}\rho\frac{{\partial{v_{i}}}}{{\partial t}}d\Omega}dt}+\int_{{\Omega_{n}}}^ {}{\delta{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)\rho{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)d\Omega}\\
 & +\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\frac{{\partial\delta v_{i}^ {}}}{{\partial{x_{j}}}}{\sigma_{ij}}\left({{\mathbf{x}},{t_{n}}}\right)d\Omega}dt}+\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\frac{{\partial\delta v_{i}^ {}}}{{\partial{x_{j}}}}{C_{ijkl}}{\psi_{kl}}d\Omega}dt}\\
 & -\int_{{I_{n}}}^ {}{\int_{\Gamma_{i}^{h}}^ {}{\delta{v_{i}}f_{i}^{s}ds}dt}-\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{\delta{v_{i}}\rho{b_{i}}d\Omega}dt}\\
 & -\int_{{\Omega_{n}}}^ {}{\delta{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)\rho{v_{i}}\left({{\mathbf{x}},t_{n}^{-}}\right)d\Omega}=0
\end{split}
$$ {#eq-ch3-139}
:::

::: {#rmk-12}
In the two-field TDG/ST/FEM, both displacement and velocity are the primary unknown. This approach usually yields a large system of linear equations which in turn increases the computation cost. In this case, special procedures are required to reduce the size of the problem. Such procedures are mainly based on the predictor-multicorrector schemes [@Bonelli2002; @Bonelli2003; @Mancuso2003; @Kunthong2005].
:::

::: {#rmk-13}
In v-ST/FEM, velocity field is the primary unknown, therefore, the resultant system of linear equations will be significantly smaller than the one in case of two-field TDG/ST/FEM. Accordingly, the computation cost of v-ST/FEM is lower than the two-field TDG/ST/FEM, however, the cost is still higher than the classical semi-discrete schemes.
:::

## Implementation of v-ST/FEM

Consider $Q_{n,e}=\Omega_{e}\times I_{n}$ denoting the space-time finite element. Let $n_{e}$ be the total number of spatial nodes in the spatial finite element $\Omega_{e}$. Let $v_{i}(\mathbf{x},t_{n}^{+})$ and $v_{i}(\mathbf{x},t_{n+1}^{-})$ be the spatial nodal velocities on the bottom and top faces of space-time slab $Q_{n}$, respectively (see @fig-ch3-17). Furthermore, time $t\in I_{n}$ is given by

$$
t=T_{1}(\theta)t_{n}+T_{2}(\theta)t_{n+1},\quad\forall\theta\in\left[-1,1\right]
$$ {#eq-ch3-140}

where

$$
\begin{aligned}
T_{1}(\theta) & =\frac{1-\theta}{2} & T_{2}(\theta) & =\frac{1+\theta}{2}
\end{aligned}
$$ {#eq-ch3-141}

![Schematic diagram of the space-time slabs $Q_{n}$ and $Q_{n+1}$](./figures/ch3-fig-17.svg){#fig-ch3-17 width="50%" .lightbox}

![Schematic diagram of the space-time finite element $Q_{n,e}$](./figures/ch3-fig-18.svg){#fig-ch3-18}

The linear trial functions for the velocity defined on $Q_{n,e}$ are give by

$$
{v_{i}}\left({{\mathbf{x}},t}\right)={}_{}^{a}{v_{iI}}{T_{a}}\left(\theta\right){N^{I}}\left({\xi,\eta}\right),\quad\,a=1,2;\quad\,I=1,\cdots,{n_{e}}
$$ {#eq-ch3-142}

where $_{}^{a}{v_{iI}}$ denotes the space-time nodal values of velocity; $a=1$ and $a=2$ correspond to the bottom and top face (i.e., temporal nodes) of the space-time element, and $I=1,\cdots,n_{e}$ denotes the spatial node of the space-time element (see @fig-ch3-18). $N^{I}(\xi,\eta)$ are the spatial shape functions defined on the local domain. Accordingly, $v_{i}(\mathbf{x},t_{n}^{+})$ and $v_{i}(\mathbf{x},t_{n+1}^{-})$ are given by @eq-ch3-143 and @eq-ch3-144, respectively.

$$
{v_{i}}\left({{\mathbf{x}},t_{n}^{+}}\right)={}_{}^{1}{v_{iI}}{N^{I}}\left({\xi,\eta}\right)
$$ {#eq-ch3-143}

$$
{v_{i}}\left({{\mathbf{x}},t_{n+1}^{-}}\right)={}_{}^{2}{v_{iI}}{N^{I}}\left({\xi,\eta}\right)
$$ {#eq-ch3-144}

The displacement field corresponding to the above-mentioned local interpolation for the velocity field (cf. @eq-ch3-144) is obtained by using @eq-ch3-136.

$$
{u_{i}}\left({{\mathbf{x}},t}\right)={u_{i}}\left({{\mathbf{x}},{t_{n}}}\right)+{\tilde{T}_{1}}\left(\theta\right){v_{i}}\left({{\mathbf{x}},{t_{n}}}\right)+{\tilde{T}_{2}}\left(\theta\right){v_{i}}\left({{\mathbf{x}},{t_{n+1}}}\right)
$$ {#eq-ch3-145}

where

$$
\begin{aligned}
{{\tilde{T}}_{1}}\left(\theta\right) & =\frac{{\Delta{t_{n}}}}{2}\left[{1-T_{1}^{2}\left(\theta\right)}\right], & {{\tilde{T}}_{2}}\left(\theta\right) & =\frac{{\Delta{t_{n}}}}{2}T_{2}^{2}\left(\theta\right)\label{eq-ch3-146}
\end{aligned}
$$

are the quadratic shape function defined on $\left[-1,1\right]$. Subsequently, $\psi_{ij}$ in @eq-ch3-138 becomes,

$$
{\psi_{ij}}\left({{\mathbf{x}},t}\right)={}^{a}{v_{iI}}{{\tilde{T}}_{a}}\frac{{\partial{N^{I}}}}{{\partial{x_{j}}}}\label{eq-ch3-147}
$$

By using the expression for $psi_{ij}$ in @eq-ch3-137 the stress $\sigma_{ij}$ inside the space-time element is described as,

$$
\begin{aligned}
{\sigma_{ij}}\left({{\mathbf{x}},t}\right) & ={\sigma_{ij}}\left({{\mathbf{x}},{t_{n}}}\right)+{{\tilde{T}}_{a}}{C_{ijkl}}\frac{{\partial{N^{I}}}}{{\partial{x_{l}}}}{}^{a}{v_{kI}}, & \forall(\mathbf{x},t)\in Q_{n,e}
\end{aligned}
$$ {#eq-ch3-148}

Lastly, the test functions for the velocity field are desribed as

$$
{\delta v_{i}}\left({{\mathbf{x}},t}\right)={}_{}^{a}{\delta v_{iI}}{T_{a}}\left(\theta\right){N^{I}}\left({\xi,\eta}\right),\quad\,a=1,2;\quad\,I=1,\cdots,{n_{e}}
$$ {#eq-ch3-149}

Using @eq-ch3-140 -- @eq-ch3-149 for the space-time finite element discretization of the weak-form which is given in Eq. @eq-ch3-139.

$$
\begin{split} & ^{a}\delta{v_{iI}}\left[{\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{{N^{I}}{T_{a}}\rho\frac{{\partial{N^{J}}{T_{b}}}}{{\partial t}}d\Omega}dt}}\right]{}^{b}{v_{iJ}}{+^{a}}\delta{v_{iI}}\left[{{\delta_{1a}}{\delta_{1b}}\int_{{\Omega_{n}}}^ {}{{N^{I}}\rho{N^{J}}d\Omega}}\right]{}^{b}{v_{iJ}}\\
 & {+^{a}}\delta v_{iI}^ {}\left\{ {\int_{{I_{n}}}^ {}{{T_{a}}\int_{{\Omega_{h}}}^ {}{\frac{{\partial{N^{I}}}}{{\partial{x_{j}}}}\sigma_{ij}^{n}d\Omega}dt}}\right\} {+^{a}}\delta v_{iI}^ {}\left[{\int_{{I_{n}}}^ {}{{T_{a}}{{\tilde{T}}_{b}}\int_{{\Omega_{h}}}^ {}{\frac{{\partial{N^{I}}}}{{\partial{x_{j}}}}{C_{ijkl}}\frac{{\partial{N^{J}}}}{{\partial{x_{l}}}}d\Omega}dt}}\right]{}^{b}{v_{kJ}}\\
 & -^{a}\delta{v_{iI}}\left\{ {\int_{{I_{n}}}^ {}{\int_{\Gamma_{i}^{h}}^ {}{{T_{a}}{N^{I}}f_{i}^{s}ds}dt}}\right\} -{}^{a}\delta{v_{iI}}\left\{ {\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{{T_{a}}{N^{I}}\rho{b_{i}}d\Omega}dt}}\right\} \\
 & -^{a}\delta{v_{iI}}\left[{{\delta_{1a}}\int_{{\Omega_{n}}}^ {}{{N^{I}}\rho{N^{J}}d\Omega}}\right]v_{iJ}^{-}=0
\end{split}
$$ {#eq-ch3-150}

Let us now use the following notation for representing the global space-time nodal vector,

$$
\left\{ {\mathbf{J}}\right\} :=\left\{ J\right\} _{i}^{a}\left(I\right)
$$

here $a=1,2$ corresponds to the temporal nodes, $I=1,\cdots$ corresponds to the $I^{th}$ spatial node, and $i=1,2$ corresponds to spatial components. Accordingly, $a=1$ corresponds to the spatial nodal values of the vector defined at the bottom space-time slab (i.e., at time $t_{n}^{+}$). Similarly, $a=2$ corresponds to the spatial nodal values of the vector at the top space-time slab (i.e., at time $t_{n+1}^{-}$).

A typical space-time finite element matrix will be denoted as:

$$
\left[{\mathbf{K}}\right]:=\left[K\right]_{ij}^{ab}\left({I,J}\right)
$$

where $a,b=1,2$ corresponds to the temporal nodes, $I=1,\cdots,n_{e}$ corresponds to the spatial node, and $i,j=1,2$ corresponds to the spatial components.

By using these notations @eq-ch3-150 can be expressed as:

$$
\left\{ {\delta v}\right\} _{i}^{a}\left(I\right)\cdot\left(\begin{gathered}\left[M\right]_{ij}^{ab}\left({I,J}\right)\cdot\left\{ v\right\} _{j}^{b}\left(J\right)+\left[K\right]_{ij}^{ab}\left({I,J}\right)\cdot\left\{ v\right\} _{j}^{b}\left(J\right)\\
+\left\{ {{J_{{\sigma^{n}}}}}\right\} _{i}^{a}\left(I\right)-\left\{ {{J_{ext}}}\right\} _{i}^{a}\left(I\right)-\left\{ {{J_{0}}}\right\} _{i}^{a}\left(I\right)
\end{gathered}
\right)=0
$$ {#eq-ch3-151}

Since above equation is true for all values of $\left\{ \delta v\right\} _{i}^{a}$ one can obtain the following the system of linear equations:

$$
\left[M\right]_{ij}^{ab}\left({I,J}\right)\cdot\left\{ v\right\} _{j}^{b}\left(J\right)+\left[K\right]_{ij}^{ab}\left({I,J}\right)\cdot\left\{ v\right\} _{j}^{b}\left(J\right)=\left\{ {{J_{ext}}}\right\} _{i}^{a}\left(I\right)+\left\{ {{J_{0}}}\right\} _{i}^{a}\left(I\right)-\left\{ {{J_{{\sigma^{n}}}}}\right\}_{i}^{a}\left(I\right)
$$ {#eq-ch3-152}

Matrix-vector form of above equation is given as,

$$
\left[{\mathbf{M}}\right]\cdot\left\{ {{\mathbf{\tilde{v}}}}\right\} +\left[{\mathbf{K}}\right]\cdot\left\{ {{\mathbf{\tilde{v}}}}\right\} =\left\{ {{{\mathbf{J}}_{ext}}}\right\} +\left\{ {{{\mathbf{J}}_{0}}}\right\} -\left\{ {{{\mathbf{J}}_{{\sigma^{n}}}}}\right\}
$$ {#eq-ch3-153}

If Rayleigh damping is used to model the material damping then above
equation becomes:

$$
\left[{\mathbf{M}}\right]\cdot\left\{ {{\mathbf{\tilde{v}}}}\right\} +\left[{\mathbf{K}}\right]\cdot\left\{ {{\mathbf{\tilde{v}}}}\right\} +\alpha\left[{{{\mathbf{M}}_{R}}}\right]\cdot\left\{ {{\mathbf{\tilde{v}}}}\right\} +\beta\left[{{{\mathbf{K}}_{R}}}\right]\cdot\left\{ {{\mathbf{\tilde{v}}}}\right\} =\left\{ {{{\mathbf{J}}_{ext}}}\right\} +\left\{ {{{\mathbf{J}}_{0}}}\right\} -\left\{ {{{\mathbf{J}}_{{\sigma^{n}}}}}\right\}
$$ {#eq-ch3-154}

where $\alpha$ and $\beta$ are the Rayleigh damping coefficients.

In @eq-ch3-153 -- @eq-ch3-154, $\left\{ {{\mathbf{\tilde{v}}}}\right\}$ denotes the space-time nodal vector of velocity, $\left[{\mathbf{M}}\right]$ denotes the space-time mass matrix, $\left[{\mathbf{K}}\right]$ denotes the space-time tangent stiffness matrix, $\left[{{{\mathbf{M}}_{R}}}\right]$ is the mass proportional Rayleigh damping matrix, and $\left[{{{\mathbf{K}}_{R}}}\right]$ is the stiffness proportional Rayleigh damping matrix. Furthermore, $\left\{ {{{\mathbf{J}}_{ext}}}\right\}$ denotes the space-time nodal vectors which contains the contribution of external body force and external boundary traction, $\left\{ {{{\mathbf{J}}_{0}}}\right\}$ contains the contribution of initial velocity, and $\left\{ {{{\mathbf{J}}_{\sigma}^{n}}}\right\}$ contains the contribution of initial stress $\sigma^{n}$. In what follows the finite element expressions for these matrices and vectors are described.

$$
\left[{\mathbf{M}}\right]:=\left[M\right]_{ij}^{ab}\left({I,J}\right)={\delta_{ij}}\int_{{I_{n}}}^ {}{{T_{a}}\frac{{\partial{T_{b}}}}{{\partial t}}\int_{{\Omega_{h}}}^ {}{{N^{I}}\rho{N^{J}}d\Omega}dt}+{\delta_{ij}}{\delta_{1a}}{\delta_{1b}}\int_{{\Omega_{n}}}^ {}{{N^{I}}\rho{N^{J}}d\Omega}
$$ {#eq-ch3-155}

$$
\left[{\mathbf{K}}\right]:=\left[K\right]_{ij}^{ab}\left({I,J}\right)=\int_{{I_{n}}}^ {}{{T_{a}}{{\tilde{T}}_{b}}\int_{{\Omega_{h}}}^ {}{\frac{{\partial{N^{I}}}}{{\partial{x_{p}}}}{C_{pijq}}\frac{{\partial{N^{J}}}}{{\partial{x_{q}}}}d\Omega}dt}
$$ {#eq-ch3-156}

$$
\left[{{{\mathbf{M}}_{R}}}\right]:=\left[{{M_{R}}}\right]_{ij}^{ab}\left({I,J}\right)={\delta_{ij}}\int_{{I_{n}}}^ {}{{T_{a}}{T_{b}}\int_{{\Omega_{h}}}^ {}{{N^{I}}{N^{J}}d\Omega}dt}
$${#eq-ch3-157}

$$
\left[{{{\mathbf{K}}_{R}}}\right]:=\left[{{K_{R}}}\right]_{ij}^{ab}\left({I,J}\right)=\int_{{I_{n}}}^ {}{{T_{a}}{T_{b}}\int_{{\Omega_{h}}}^ {}{\frac{{\partial{N^{I}}}}{{\partial{x_{p}}}}{C_{pijq}}\frac{{\partial{N^{J}}}}{{\partial{x_{q}}}}d\Omega}dt}
$$ {#eq-ch3-158}

$$
\left\{ {{{\mathbf{J}}_{ext}}}\right\} :=\left\{ {{J_{ext}}}\right\} _{i}^{a}\left(I\right)=\int_{{I_{n}}}^ {}{\int_{{\Omega_{h}}}^ {}{{N^{I}}{T_{a}}\rho{b_{i}}}d\Omega dt}+\int_{{I_{n}}}^ {}{\int_{\Gamma_{i}^{h}}^ {}{{N^{I}}{T_{a}}\rho f_{i}^{s}}dsdt}
$$ {#eq-ch3-159}

$$
\left\{ {{{\mathbf{J}}_{0}}}\right\} :=\left\{ {{J_{0}}}\right\} _{i}^{a}\left(I\right)={\delta_{1a}}{\delta_{ij}}\left[{\int_{{\Omega_{h}}}^ {}{{N^{I}}\rho{N^{J}}d\Omega}}\right]\left\{ {^{0}{v_{jJ}}}\right\}
$${#eq-ch3-160}

where, $\left\{ {^{0}{v_{jJ}}}\right\}$ denotes the initial velocity for the space-time slab $Q_{n}$, i.e., spatial nodal velocity at time $t_{n}^{-}$. Note that this velocity vector is known from the computation in the previous space-time slab $Q_{n-1}$.

$$
\left\{ {{{\mathbf{J}}_{{\sigma^{n}}}}}\right\} :=\left\{ {{J_{{\sigma^{n}}}}}\right\} _{i}^{a}\left(I\right)=\int_{{I_{n}}}^ {}{{T_{a}}\int_{{\Omega_{h}}}^ {}{\frac{{\partial{N^{I}}}}{{\partial{x_{j}}}}\sigma_{ij}^{n}d\Omega}dt}
$${#eq-ch3-161}

where $\sigma^{n}:=\sigma(\mathbf{x},t_{n})$ is the stress at time $t_{n}$ which is usually computed from the displacement field at time $t_{n}$.

This section briefly discusses the space-time matrices and vectors, however, a detailed description about the derivation of the space-time matrices and space-time vectors can be found in Appendix.

## Primary wave propagation in homogeneous linear elastic medium

A theoretical model of homogenous isotropic linear elastic media occupying a square domain, $\Omega=\left[0,L\right] \times \left[0,L\right]$, is considered in this section for a numerical analysis of v-ST/FEM. The governing equations of the problem are given by @eq-ch3-123 -- @eq-ch3-127.

Further, it is assumed that the elastic media is subjected to the periodic boundary conditions with no external body force (i.e., $b_i=0$ in @eq-ch3-123. The initial conditions for displacement and velocity field corresponding to @eq-ch3-126 and @eq-ch3-127 are given by:

$$
\begin{align}
    {{\mathbf{u}}_0}\left( {\mathbf{x}} \right) &= {{\mathbf{d}}_0}\cos \left( {k{\mathbf{x}} \cdot {\mathbf{\hat r}}} \right),&
    {{\mathbf{v}}_0}\left( {\mathbf{x}} \right) &= ck{{\mathbf{d}}_0}\sin \left( {k{\mathbf{x}} \cdot {\mathbf{\hat r}}} \right)
\end{align}
$${#eq-ch3-162}

where $k$ is the wave-number, $\hat{\mathbf{r}}$ is the direction vector of the wave propagation, $\mathbf{d}_0 \in \mathbb{R}^2$ denotes the direction of the motion of particles of the medium, and $c$ is related to the speed of the wave in the medium.

The above-mentioned initial conditions create a plane P-wave if $\hat{\mathbf{r}}$ vector is parallel to $\mathbf{d}_0$. The analytical solutions for displacement and velocity can be given by [@Achenbach1973a]:

$$
\begin{align}
    {\mathbf{u}}\left( {{\mathbf{x}},t} \right) &= {{\mathbf{d}}_0}\cos \left[ {k\left( {{\mathbf{x}} \cdot {\mathbf{\hat r}} - {c_p}t} \right)} \right]&
    {\mathbf{v}}\left( {{\mathbf{x}},t} \right) = {c_p}k{{\mathbf{d}}_0}\sin \left[ {k\left( {{\mathbf{x}} \cdot {\mathbf{\hat r}} - {c_p}t} \right)} \right]
\end{align}
$${#eq-ch3-163}

where

$${c_p} = \sqrt {\frac{{\lambda  + 2\mu }}{\rho }}$$

is the speed of the P-wave.

The values of constants and the parameters used for solving the problem are given below. @fig-ch3-19 depicts the physical dimensions of the problem. All variables have been dimensionless. Letting $\mathbf{u}^h(\mathbf{x},t)$ and $\mathbf{v}^h(\mathbf{x},t)$ denote the numerical solutions computed by employing the v-ST/FEM, the error in displacement field ${E_u}\left( {{t_n}} \right)$ and the error in velocity field ${E_v}\left( {{t_n}} \right)$ are then defined as:

$$
\begin{align}
    {E_u}\left( {{t_n}} \right) &: = {\left\| {{{\mathbf{u}}^h}\left( {{\mathbf{x}},{t_n}} \right) - {\mathbf{u}}\left( {{\mathbf{x}},{t_n}} \right)} \right\|_2},&
    {E_v}\left( {{t_n}} \right)&: = {\left\| {{{\mathbf{v}}^h}\left( {{\mathbf{x}},t_n^ - } \right) - {\mathbf{v}}\left( {{\mathbf{x}},{t_n}} \right)} \right\|_2}
\end{align}
$${#eq-ch3-164}

where $\left\| {\, \cdot \,} \right\|_2$ denotes the $L_2$ norm given by

$$
\begin{equation}
    {\left\| {\mathbf{u}} \right\|_2} = {\left[ {\int_\Omega ^{} {{\mathbf{u}}\left( {{\mathbf{x}},t} \right) \cdot {\mathbf{u}}\left( {{\mathbf{x}},t} \right)d\Omega } } \right]^{1/2}}
\end{equation}
$${#eq-ch3-165}

![(a) Schematic diagram of the spatial domain for P-wave propagation problem, (b) bilinear quadrilateral element (Quad4), (c) linear triangular element (Tria3)](./figures/ch3-fig-19.svg){#fig-ch3-19}

![Rate of convergence of the solutions in space computed at time $t=1.0$ sec: (a) displacement, and (b) velocity](./figures/ch3-fig-20.svg){#fig-ch3-20}

![Rate of convergence of the solutions in time domain at time $t=35.0$ sec: (a) displacement, and (b) velocity.](./figures/ch3-fig-21.svg){#fig-ch3-21}

The numerical experiments for determining the convergence rate of the solution in the space domain are performed on two sequences of regular linear triangular (Tria3) and bilinear quadrilateral (Quad4) meshes (see Fig. @fig-ch3-19), while keeping the time-step fixed $\Delta t=0.1$ sec.

Each sequence consists of four meshes with a decreasing mesh size. In @fig-ch3-20, the L2 norm of the errors in the displacement and velocity fields at time $t=1.0$ sec are given in relation to mesh spacing parameter $h$. Based on the convergence results, it can be stated that the v-ST/FEM formulation is nearly second-order accurate in the space for both Quad4 and Tria3 elements. In @fig-ch3-20, it is observed that the error for the triangular spatial mesh (Tria3) is less than that of the quadrilateral spatial element (Quad4) for the same mesh spacing. This can be attributed to the perfect alignment of the diagonal of the triangular elements with the characteristic lines of the wave propagation.

![(a) Exact displacement (x-component) waveforms, (b) displacement (x-component) waveforms obtained by using the v-ST/FEM, (c) exact velocity (x-component) waveforms, (d) velocity (x-component) waveforms obtained by using the v-ST/FEM.](./figures/ch3-fig-22.png){#fig-ch3-22}

Furthermore, the convergence of the solutions (displacement and velocity) in the time domain is illustrated in @fig-ch3-21. It can be seen that both displacement and velocity fields computed by the present method are third order accurate in time. Note that the results presented in @fig-ch3-21 are consistent with @eq-ch3-114.

@fig-ch3-22 illustrates the spatial variation of the displacement and velocity fields obtained by v-ST/FEM. The results are obtained at time $t=35$ seconds with linear triangular mesh of size $h=0.25$ and uniform time-step of size $\Delta t = 0.1$ sec. The results advocate the ability of v-ST/FEM to maintain the high-order accuracy of the long-term solutions, especially the displacement and velocity waveforms. This can be attributed to the very low numerical dissipation and dispersion characteristics of the present method (see also @fig-ch3-15 and @fig-ch3-16). To further emphasize the accuracy of the long-term solutions, the time histories of the computed displacement and velocity at the midpoint P1 (see @fig-ch3-19) and the corresponding relative errors are presented in @fig-ch3-23.

![Temporal variation of the (a) displacement, (b) relative error in displacement, (c) velocity, and (d) relative error in velocity, at the
midpoint P1 computed by using v-ST/FEM.](./figures/ch3-fig-23.png){#fig-ch3-23}

## Impulsive response of a fixed-free pile

In this section, we consider a pile of length $L=50$ m with a unit cross-section area having one fixed end and one end loaded by an axial impulsive force given by a step function, as shown in Fig. [6](#fig-ch3-24){reference-type="ref" reference="fig-ch3-24"}. The mass density $\rho$ is $2500.0$ $kg/m^3$ and the Young's modulus $E$ is $1.0\times 10^{10}$ $N/m^2$. The magnitude of the implusive force is $1.0 \times 10^6$ $N$. Under these circumstances, the analytical solutions for the stress and velocity fields are discontinuous in the spatial domain and given by the step functions. Furthermore, the displacements are given by the piecewise continuous linear functions [@Cormeau1991; @Li1998; @Verruijt2009].

![Geometry, boundary conditions, and impulse loading for fixed-free pile problem.](./figures/ch3-fig-24.svg){#fig-ch3-24}

To solve the problem by using v-ST/FEM, uniform linear spatial elements of size $h=0.1m$ and a uniform time-step size $\Delta t = 10^{-4}$ sec have been adopted for discretizing the space and time domain, respectively. To assess the performance of the present method, the problem is also solved with semi-discretized FEM techniques using the same mesh parameters. In the latter case, the trapezoidal rule and the HHT-$\alpha$ method (with $\alpha=-1/3$) have been used as the time-stepping algorithms. The results of the stress field, velocity field, and displacement field obtained by different schemes are compared in @fig-ch3-25, @fig-ch3-26, and @fig-ch3-27, respectively.

It is observed that the solutions for the stress and velocity fields obtained by the semi-discrete algorithms contain severe oscillations. Even worse, no improvement whatsoever is obtained when refining the mesh spacing or the time-step size. Further, in the case of the Newmark-beta method, these oscillations are present in the whole spatial domain; and thus, the accuracy of the stress and velocity fields deteriorate over time (see @fig-ch3-25 and @fig-ch3-26). The poor performance with the Newmark-beta method is due to the absence of algorithmic damping, as discussed in previous sections. The HHT-$\alpha$ method improves the solutions by attenuating higher frequencies; however, the results are not satisfactory for the selected time-step size as the oscillations are still present. It is remarkable that v-ST/FEM completely localizes the oscillations in the stress and velocity fields near the point of discontinuity and yields very accurate solutions. The localization phenomenon can be attributed to the presence of jump discontinuity in the velocity field. The jump in the velocity field adds artificial viscous damping to the system, subsequently, increasing the accuracy and stabilizing the solutions [@Hulbert1990; @Johnson1993]. Furthermore, the presence of overshooting and undershooting around the point of discontinuity is related to the famous Gibbs phenomenon in the Fourier analysis [@Olver2016].

![Stress field computed by employing the Newmark-beta method (First column), HHT-$\alpha$ (Second column), and v-ST/FEM (Third column) at various timesteps with linear spatial elements. Direction of wave propagation is denoted by the arrow, and dotted lines represent the analytical solutions.](./figures/ch3-fig-25.png){#fig-ch3-25}

![Velocity field computed by employing the Newmark-beta method (First column), HHT-$\alpha$ (Second column), and v-ST/FEM (Third column) at various time-steps with linear spatial elements. Direction of wave propagation is denoted by the arrow, and dotted lines represent the analytical solutions.](./figures/ch3-fig-26.png){#fig-ch3-26}

![Displacement field computed by employing the Newmark-beta method (First column), HHT-$\alpha$ (Second column), and v-ST/FEM (Third column) at various timesteps with linear spatial elements. Direction of wave propagation is denoted by the arrow.](./figures/ch3-fig-27.png){#fig-ch3-27}

## Dynamic plate load test (DPLT)

In this section, an attempt is made to validate v-ST/FEM by simulating the dynamic plate loading test (DPLT) using a light falling weight deflectometer (LFWD). DPLT using LFWD is a non-destructive technique for a quick assessment of the field compaction quality. In DPLT, a rigid circular loading plate of radius $0.15$ m is placed on the soil surface and subjected to an impulse load generated by the falling weight from a specified height onto the plate. Subsequently, the induced soil movements are recorded and the dynamic resilience modulus of the tested material is computed by some empirical relations. A complete description of the test can be found elsewhere [@Adam2009; @Tawfik2017].

![Geometry, boundary conditions and spatial mesh adopted for simulating DPLT using v-ST/FEM.](./figures/ch3-fig-28.png){#fig-ch3-28}

In the present study, DPLT is simulated by a two-dimensional axisymmetric v-ST/FEM model, where the center of the loading plate is positioned along the axis of symmetry. The geometry, boundary conditions, and spatial mesh of the finite element model are depicted in @fig-ch3-28. There are $5201$ linear triangular elements (Tria3) and $2688$ nodes present in the spatial mesh. The impulse load due to the falling weight is modeled by using an equivalent uniform vertical stress pulse of amplitude $100$ kPa and a time duration of $20$ ms [@Adam2009]. The stress pulse $h(t)$ acting on the loading plate is defined by a half sine wave as:

$$
h(t)=-10^{5} \sin(50 \pi t) \, \text{N/m}{}^2
$$

The total simulation time $T$ is set to $30$ ms, and the linear time elements of size $\Delta t = 1$ ms have been adopted for discretizing the time domain. Subsequently, the results computed by the proposed method are compared with the two DPLT studies available in the literature. The first study is denoted as S1 where the in-situ LFWD test was conducted by [@Tawfik2017], and the second study is denoted by S2 in which the numerical investigation was conducted by [@Adam2009] using the Boundary Element Method (BEM). In these studies (i.e., S1 and S2), the soils have different values for Young's modulus $E$ and common values for Poisson's ratio $\nu$ and mass density $\rho$, as shown in @tbl-ch3-3.

::: {#tbl-ch3-3}

|  Elastic parameters   |     Method         |  $d_{max}$  | $p_{max}$    |
| --------------------- | ------------------ | ----------- | ------------ |
| $E=65$ Mpa            |      v-ST/FEM      |  $0.33$ mm  |   $90.0$ kPa |
| $\nu = 0.212$         |    In-situ study (S1) |    $0.35$ mm   |  $100$ kPa |
| $\rho = 2000$ $kg/m^3$|                    |             |              |
| $E=32$ Mpa            |         v-ST/FEM   |   $0.67$ mm |   $90.0$ kPa |
| $\nu = 0.212$         |  BEM study (S2)    |  $0.65$ mm  |  $93.1$ kPa  |

:  List of elastic material parameters, and the results of maximum plate deflection and maximum soil-plate normal contact stress obtained by different schemes.

:::

The results of the time histories of the vertical displacement and the velocity of the loading plate, and the soil-plate contact stress are plotted in @fig-ch3-29. The maximum deflection of the plate $d_{max}$ obtained by v-ST/FEM is about $0.33$ mm and $0.67$ mm corresponding to the material parameters used in the studies S1 and S2, respectively. From @tbl-ch3-3, it is evident that these computed values for $d_{max}$ are in good agreement with those reported in studies S1 and S2. In the present study, the maximum soil-plate normal contact stress $p_{max}$ is about $90$ kPa; which lower than the values reported in studies S1 and S2 (see @tbl-ch3-3). This may be due to the use of linear triangular elements and the absence of the soil-plate interface elements in the present formulation. Furthermore, the maximum deflection of the plate occurs at time $t_d = 11$ ms, and the maximum contact stress $p_{max}$ occurs at time $t_p =10$ ms. The displacement contours at these time steps are depicted in @fig-ch3-30, where displacements are normalized with respect to the maximum deflection of the plate.

![Temporal variation of the vertical displacement of plate (top-left), the vertical velocity of plate (top-right), and soil-plate contact stress (bottom-left), and variation of normal contact stress with the plate displacement (bottom-right) obtained by v-ST/FEM.](./figures/ch3-fig-29.png){#fig-ch3-29}

![Normalized contours of the vertical displacements computed by v-ST/FEM at various time-steps.](./figures/ch3-fig-30.png){#fig-ch3-30}
