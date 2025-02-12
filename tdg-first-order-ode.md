# TDG/FEM for the first order ODE

To understand the working principle of the time-discontinuous Galerkin
finite element method (TDG/FEM) consider the following initial value
problem [^1] described by a first order ordinary differential equation
equation (ODE).

$$
\begin{split}
    \frac{{du}}{{dt}} + \lambda u - f\left( t \right)&=0 \quad t \in \left[ {0,T} \right],\\
    u(0)&=u_0
\end{split}
$${#eq-ch3-1}

where $T$ is the total time, $u:=u(t)$ is the
primary unknown, $f(t)$ represents the excitation function, and $u_0$ is
the initial condition. The analytical solution to this problems takes
the form of

$$
u\left( t \right) = {u_0}{e^{-\lambda t}} + \int_0^t {{e^{-\lambda \left( {t - \tau } \right)}}f\left( \tau  \right)d\tau }
$$

Let us now consider a non-uniform subdivision for the time domain
$\left[ {0,T} \right]$,

$$0 = {t_0} < {t_1} <  \cdots  < {t_N} = T,$$

with

$$
\begin{aligned}
  I_n &= (t_n,t_n+1), & \Delta t_n &= t_{n+1}-t_{n}, &
  \Delta t = \mathop {\max }\limits_{\left( {0 \leqslant n \leqslant N - 1} \right)} \left( {\Delta {t_n}} \right).
\end{aligned}
$$

![Discretization of time domain $\left[ 0, T\right]$ by using the
time-finite elements.](./figures/ch3-fig-1.svg){#fig-ch3-1}

In the finite element method (FEM), the time interval, $I_n$, denotes
the interior of the one-dimensional time finite element, and
$\partial {I_n} = \left\{ {{t_n},{t_{n + 1}}} \right\}$ denotes the
boundary of $I_n$ [^2]. Accordingly, the finite element discretization
of the continuous time domain can be written as (see also @fig-ch3-1),

$$
{I_h}: = \bigcup\limits_{n = 0}^{N - 1} {\left( {{I_n} \cup \partial {I_n}} \right)}
$$

In TDG/FEM, the solutions to the problem, $u^h(t)$, are considered to be
discontinuous at the boundary of the time-domain. The solutions,
however, remain continuous inside the time-finite element $I_n$, and
approximated by piecewise polynomials (e.g., Lagrange polynomials).

Therefore, the discontinuity in $u^h(t)$ occurs at times that belong to
the finite set $\left\{ {{t_0},{t_1}, \cdots ,{t_N}} \right\}$. The jump
discontinuity in time for $u^h$ is given by

$$
{\left[\kern-0.15em\left[ {{u^h}}
\right]\kern-0.15em\right]_n} = u_n^ +  - u_{n}^{-}
$${#eq-ch3-3}

where

$$
\begin{aligned}
        u_n^ +  &= \mathop {\lim }\limits_{\varepsilon  \to 0} {u^h}\left( {t + \varepsilon } \right), &
        u_n^ -  &= \mathop {\lim }\limits_{\varepsilon  \to 0} {u^h}\left( {t - \varepsilon } \right)
\end{aligned}
$${#eq-ch3-4}

are the discontinuous values of $u^h$ at time $t=t_n$. @fig-ch3-2 illustrates the concept of time discontinuity.

![Schematic diagram of time discontinuous approximation: (a) piecewise
linear interpolation, and (b) piecewise quadratic
interpolation.](./figures/ch3-fig-2.svg){#fig-ch3-2}

Let $C^0(\star)$ denote the space of piecewise continuous functions
defined on domain $(\star)$. In addition, consider
${\wp_{l}}\left( {{I_n}} \right)$, the collection of all polynomials
defined on $I_n$ with a total degree of no more than $l$. Accordingly,
the functional space required for the TDG/FEM is given by,

$$
\Im_{l}^{h}: = \left\{ {\left. {{u^h}} \right|{u^h} \in {C^0}\left( {\bigcup\limits_{n = 0}^{N - 1} {{I_n}} } \right),{{\left. {{u^h}} \right|}_{{I_n}}} \in {\wp_l}\left( {{I_n}} \right)} \right\}
$${#eq-ch3-5}

where, ${{u^h}}$, denotes the global solutions which can be
discontinuous at time steps given by the end points of the interval
$I_n$, and ${{{\left. {{u^h}} \right|}_{{I_n}}}}$ denotes the
restriction of $u_h$ to $I_n$.

## Weak forms for the TDG/FEM

For any $\delta u^h, u^h \in \Im _l^h$ consider the following
variational form corresponding to @eq-ch3-1,

$$
\delta \Pi_n := \int_{{I_n}}^{} {\delta u^h\left( {\frac{{du^h}}{{dt}} + \lambda u^h - f\left( t \right)} \right)dt}
$${#eq-ch3-6}

or

$$
\delta \Pi_n = \int_{{I_n}}^{} {\delta u^h\frac{{du^h}}{{dt}}dt}  + \int_{{I_n}}^{} {\delta u^h\left[ {\lambda u^h - f\left( t \right)} \right]dt}
$${#eq-ch3-7}

Here, once again note that the domain of integration is $I_n$, and $u^h$
is discontinuous at times $t_n$ and $t_{n+1}$.

Use of the integration by parts for the first term in @eq-ch3-7 will
will result in the following variation form for the TDG/FEM,

$$
\delta \Pi_n =- \int_{{I_n}}^{} {u^h\frac{{d\delta u^h}}{{dt}}dt}  + \left[ {\delta u^h \cdot {u^*}} \right]_{{t_n}}^{{t_{n + 1}}} + \int_{{I_n}}^{} {\delta u^h\left[ {\lambda u^h - f\left( t \right)} \right]dt}
$${#eq-ch3-8}

where, $u^*$ denotes the unique value of $u^h$ to be used at the
end-points of the interval $I_n$ and obtained by computing the
information from the neighboring time intervals. A detailed definition
concerning its concrete value will be soon introduced in the following text.

Subsequent use of the integration by parts for the first term in @eq-ch3-8 will yield another variation form for the TDG/FEM. This variation form is described below.

$$
\delta \Pi_n = \int_{{I_n}}^{} {\delta u^h\frac{{du^h}}{{dt}}dt}  + \left[ {\delta u^h\left( {{u^*} - u^h} \right)} \right]_{{t_n}}^{{t_{n + 1}}} + \int_{{I_n}}^{} {\delta u^h\left[ {\lambda u^h - f\left( t \right)} \right]dt}
$${#eq-ch3-9}

Mathematically, the formulations given in @eq-ch3-8 and @eq-ch3-9 are the TDG/FEM schemes for Eq. @eq-ch3-1 in weak and strong form, respectively. In @eq-ch3-8 due to the presence of first order time derivative of the test function $\delta u^h$, the test function should be continuous in time. In @eq-ch3-9, however, the continuity of the test functions is not required since this form does not include the time derivative of $\delta u^h$.

In @eq-ch3-8, recall that the solution $(u^h)$ is not well
defined at the $\partial I_n =\left\{ {{t_n},{t_{n + 1}}} \right \}$ as
it takes two values at each end-points of $I_n$. Hence, a unique
representative value $u^*$ has been adopted in the TDG formulation which
should be specified only at time $\partial I_n$. The value of $u^*$ at
time $t_n$ can be computed by using $u_n^+$ and $u_n^-$ (see the
footnote). [^3] The most widely used definition for $u^*$ is based on
the upwind flux treatment in time [@Eriksson1985; @Chen2006; @Cockburn2003; @Hesthaven2007] and given by,

$$
u^{*}:=\left\{ \begin{array}{cc}
u_{0} & \text{if}\quad n=0\\
u_{n}^{-} & otherwise
\end{array}\right.
$${#eq-ch3-10}

The choice of such a definition is inspired by the fact that in a
transient problem the solutions at time $t_n$ should be equal to the
value of its immediate past $t_n-\varepsilon$. Therefore, it is natural
to start the numerical procedure at $t=t_0$ with the prescribed initial
condition $u^*=u_0$, and transport this idea to the subsequent
time-slabs.

Incorporating the definition of $u^*$ in @eq-ch3-9, the TDG/FEM weak,

$$
\delta {\Pi _n} = \int_{{I_n}}^{} {\delta {u^h}\frac{{d{u^h}}}{{dt}}dt}  + \delta {u^h}\left( {{t_n}} \right){\left[\kern-0.15em\left[ {{u^h}}
\right]\kern-0.15em\right]_n} + \int_{{I_n}}^{} {\delta {u^h}\left[ {\lambda {u^h} - f\left( t \right)} \right]dt}
$${#eq-ch3-11}

Let us now focus on the TDG/FEM weak form for the initial value problem
(see Eq. @eq-ch3-1). By multiplying the residual of Eq. @eq-ch3-1 with the test function $\delta u^h$ and then integrating in the time domain $\left[ {0,T} \right]$ following weak form can be obtained.

$$
\delta \Pi : = \int_I^{} {\delta u\left[ {\frac{{du}}{{dt}} + \lambda u - f\left( t \right)} \right]dt}  = \sum\limits_{n = 0}^{N - 1} {\delta {\Pi _n}}  = 0
$${#eq-ch3-12}

Due to the weakly prescribed initial condition, the weak form given
above transforms into a local weak form defined on each time slab $I_n$.

$$
\delta {\Pi _n} = \int_{{I_n}}^{} {\delta {u^h}\frac{{d{u^h}}}{{dt}}dt}  + \delta  {u^h}\left( {{t_n}} \right){\left[\kern-0.15em\left[ {{u^h}}
\right]\kern-0.15em\right]_n} + \int_{{I_n}}^{} {\delta {u^h}\lambda {u^h}dt}  - \int_{{I_n}}^{} {\delta {u^h}f\left( t \right)dt}  = 0
$${#eq-ch3-13}

Rearranging the terms in @eq-ch3-13 and using the expression for the jump discontinuity in time (see Eq. @eq-ch3-3),

$$
\int_{{I_n}}^{} {\delta {u^h}\frac{{d{u^h}}}{{dt}}dt}  + \delta {u_n}u_n^ +  + \int_{{I_n}}^{} {\delta {u^h}\lambda {u^h}dt}  = \delta {u_n}u_n^ -  + \int_{{I_n}}^{} {\delta {u^h}f\left( t \right)dt}
$${#eq-ch3-14}

in which the first term on right hand side depicts the initial condition
for each time-slabs $I_n$. Here note that the information about $u_n^-$
is already obtained from the computation in the previous time slab
$I_{n-1}$. Further, in the beginning of the process (i.e., $n=0$) the
initial condition is incorporated in the computation according to the
definition in @eq-ch3-10. In this way, initial boundary condition becomes the part of the solution and need not to be specified in explicit terms.

## Implementation of the TDG/FEM

One dimensional time-finite elements can be employed to discretize the
TDG/FEM weak form presented in Eq.
[\[eq:ch3-eq-11\]](#eq:ch3-eq-11){reference-type="eqref"
reference="eq:ch3-eq-11"}. In this process, the solutions are locally
approximated by using the shape functions. These shape functions are
generally given by the Lagrange polynomials, and determine the degree of
accuracy. To construct a $p$-order Lagrange polynomial in time-domain,
$p+1$ nodes should be selected inside the time slab $I_n$; two nodes are
always located at the end-points $t_n$ and $t_{n+1}$, and remaining
$p-1$ nodes are distributed inside $I_n$ (see Fig.
[3](#fig:ch3-fig-3){reference-type="ref" reference="fig:ch3-fig-3"}).
The internal nodes are generally located at equidistant points in $I_n$.
[^4]

Consider the parent time-domain $I_\theta:=\left[-1,1\right]$ on which
the Lagrange polynomials $T^{(p)}_i(\theta)$ are defined. Here, $p$
denotes the degree of polynomial, and $i=1,\cdots,p+1$ corresponds to
the $p+1$ locations

$$
\left\{ {{\theta _1},{\theta_2}, \cdots ,{\theta _{p + 1}}} \right\}
$$

in $I_\theta$ with $\theta_1=-1$ and $\theta_2=+1$. Let

$$\left\{ {t_1^{(n)},t_2^{(n)}, \cdots ,t_{p}^{(n)},t_{p + 1}^{(n)}} \right\}$$

be the set of $p+1$ time-nodes in $I_n$ with $t_1^{(n)}=t_n$,
$t_2^{(n)}=t_n+1$. The internal nodes are represented by $t_i^{(n)}=t_n$
for $i=3,\cdots,p+1$ (see @fig-ch3-3).

![Conceptual diagram of (a) $p$-order time-finite element with $p+1$
local time nodes, (b) two node linear time-finite
element.](./figures/ch3-fig-3.svg){#fig-ch3-3}

Further, the relationship between $t \in I_n$ and $\theta \in I_\theta$
is given by a linear mapping,

$$t(\theta) = \frac{{\left( {1 - \theta } \right)}}{2}{t_n} + \frac{{\left( {1 + \theta } \right)}}{2}{t_{n + 1}}$${#eq-ch3-15}

The expression for the Lagrange polynomial of degree $p$ is described as
follows:

$$
T_i^{(p)} = \prod\limits_{\begin{subarray}{l}
                                j = 1 \\
                                j \ne i
                                \end{subarray}}^{p + 1}
{\frac{{\theta  - {\theta _j}}}{{{\theta_i} - {\theta _j}}}}
$${#eq-ch3-16}

Here, note that

$$
T_{i}^{\left(p\right)}\left(\theta_{j}\right)=\left\{ \begin{array}{cc}
        1 & \text{if}\quad i=j\\
        0 & \text{if}\quad i\ne j
        \end{array}\right.
$$

Accordingly, the $p$-order local approximation for the trial function
$u^h$ reads,

$$
{u^h} = T_1^{(p)}u_n^ +  + T_2^{(p)}u_{n + 1}^ -  + \sum\limits_{a = 3}^{p + 1} {T_a^{(p)}u_a^{(n)}}
$${#eq-ch3-17}

In @eq-ch3-17, $u_n^+$ and $u_{n+1}^-$ represent the discontinuous
nodal values at times $t=(t_n+\varepsilon) \in I_n$ and
$t=(t_{n+1}-\varepsilon) \in I_n$ (see also #eq-ch3-4). Further, the values of $u^h$ at the internal nodes of $I_n$ are given by $u_a^{(n)}$, for $a=3,\cdots,p+1$, which are well defined in the time-slab $I_n$. For the sake of clarity, let us denote,
$$
u^{(n)}_1=u_n^+,
$$

$$
u^{(n)}_2=u_{n+1}^-.
$$

With such convention the Eq. @eq-ch3-17 becomes

$${u^h} = \sum\limits_{a = 1}^{p + 1} {T_a^{(p)}u_a^{(n)}}$${#eq-ch3-18}

Similarly, the $p$-order local approximation for the test function
$\delta u^h$ can be given by,

$$
\delta {u^h} = \sum\limits_{a = 1}^{p + 1} {T_a^{(p)}\delta u_a^{(n)}}
$${#eq-ch3-19}

Subsequently, using the above interpolations for $u^h$ and $\delta u^h$
in the weak-form given by @eq-ch3-14, one can get the following discretized form. Henceforth, the use of summation symbol is avoided and Einstein summation convention is implied.

$$
\begin{split}
    &
    \delta u_a^{(n)}\left\{ {\left[ {\int_{{I_n}}^{} {T_a^{(p)}\frac{{dT_b^{(p)}}}{{dt}}dt} } \right]u_b^{(n)} + \left[ {\int_{{I_n}}^{} {T_a^{(p)}\lambda T_b^{(p)}dt} } \right]u_b^{(n)} - \int_{{I_n}}^{} {T_a^{(p)}f\left( t \right)dt} } \right\}
    \\
    &
    +
    \delta u_1^{(n)}\left( {u_1^{(n)} - u_n^ - } \right) = 0
\end{split}
$${#eq-ch3-20}

Further, using the Kronecker delta function
$\delta_{ab}$ in the last term of the above equation,

$$\delta u_1^{(n)}\left( {u_1^{(n)} - u_n^ - } \right) = \delta u_a^{(n)}{\delta _{a1}}{\delta _{b1}}\left( {u_b^{(n)} - u_n^ - } \right),$$

Accordingly, @eq-ch3-20 becomes,

$$
\delta u_a^{(n)}\left\{ {\left[ {\int_{{I_n}}^{} {T_a^{(p)}\frac{{dT_b^{(p)}}}{{dt}}dt}  + {\delta_{a1}}{\delta _{b1}} + \int_{{I_n}}^{} {T_a^{(p)}\lambda T_b^{(p)}dt} } \right]u_b^{(n)} - {\delta _{a1}}u_n^ -  - \int_{{I_n}}^{} {T_a^{(p)}f\left( t \right)dt} } \right\} = 0
$${#eq-ch3-21}

Since @eq-ch3-21 is true for all $\delta_a^{(n)}$, one can get the following system of $p+1$ algebraic equations.

$$
\left[ m \right]_{}^{ab}u_b^{(n)} + \lambda {\left[ c \right]^{ab}}u_b^{(n)} = {\delta _{a1}}u_n^ -  + {\left\{ J_{ext} \right\}^a}
$${#eq-ch3-22}

where $\left[ m \right]_{}^{ab}$ and $\left[ c \right]_{}^{ab}$ are the
$(p+1)\times (p+1)$ finite element matrices, and
${\left\{ f \right\}^a}$ is a vector of length $p+1$. The matrix-vector
form of above equation is depicted by

$$
            \left[ {\mathbf{m}} \right]\left\{ {{\mathbf{\tilde u}}} \right\}
            + \lambda \left[ {\mathbf{c}} \right]\left\{ {{\mathbf{\tilde u}}} \right\} = \left\{ {{\mathbf{e}}_1^p} \right\}u_n^ -  + \left\{ {\mathbf{J_{ext}}} \right\}
$${#eq-ch3-23}

The details about the terms present in @eq-ch3-21 and @eq-ch3-22 are given below.

$$
\left[ {\mathbf{m}} \right]: = \left[ m \right]_{}^{ab} = \int_{{I_n}}^{} {T_a^{(p)}\frac{{dT_b^{(p)}}}{{dt}}dt}  + {\delta_{a1}}{\delta _{b1}}
$${#eq-ch3-24}

$$
\left[ {\mathbf{c}} \right]: = {\left[ c \right]^{ab}} = \int_{{I_n}}^{} {T_a^{(p)}T_b^{(p)}dt}
$${#eq-ch3-25}

$$\left\{ {{\mathbf{e}}_1^p} \right\} := {\delta _{a1}} = \left\{ {1,0, \cdots ,0} \right\}$${#eq-ch3-26}

$$\left\{ {\mathbf{J_{ext}}} \right\}: = {\left\{ J_{ext} \right\}^a} = \int_{{I_n}}^{} {T_a^{(p)}f\left( t \right)dt}$${#eq-ch3-27}

In this way, at the beginning of the computation (i.e., n=0) the initial
condition $u_0$ is used to compute the right hand side of @eq-ch3-23, and the system of linear equations is solved to obtain $u_1^-$. This information is then used for constructing the right hand side for the next time slab.

In the subsequent sections, the numerical analysis of the TDG/FEM for
the first order ODE is discussed. Henceforth, for the sake of
simplicity, only the case of linear interpolation in time (i.e., $p=1$)
is discussed.

In case of the linear interpolation in time the shape functions are
given by

$$
\begin{aligned}
T_1(\theta):=T_1^{(1)}&=\frac{1-\theta}{2} & T_2(\theta):=T_2^{(1)}&=\frac{1+\theta}{2}.
\end{aligned}
$${#eq-ch3-28}

In this case, two linear equations in two unknowns
$u_n^+$ and $u_{n+1}^-$ can be obtained. These equations are described
below in the matrix-vector form.

$$
  \frac{1}{2}\left[\begin{array}{cc}
  1 & 1\\
  -1 & 1
  \end{array}\right]\left\{ \begin{array}{c}
  u_{n}^{+}\\
  u_{n+1}^{-}
  \end{array}\right\} +\lambda\frac{\Delta t_{n}}{6}\left[\begin{array}{cc}
  2 & 1\\
  1 & 2
  \end{array}\right]\left\{ \begin{array}{c}
  u_{n}^{+}\\
  u_{n+1}^{-}
  \end{array}\right\} =\left\{ \begin{array}{c}
  u_{n}^{-}\\
  0
  \end{array}\right\} +\left\{ \begin{array}{c}
  J_{ext}^{1}\\
  J_{ext}^{2}
  \end{array}\right\}
$${#eq-ch3-29}

where $J_{ext}^1$ and $J_{ext}^2$ are
given by

$$
\begin{aligned}
{J_{ext}^1} &= \int_{{I_n}}^{} {{T_1}f\left( t \right)dt}
\\
{J_{ext}^2} &= \int_{{I_n}}^{} {{T_2}f\left( t \right)dt}
\end{aligned}
$$

## Stability analysis of the TDG/FEM

Let us denote the exact value of the solution at any instant $t_n$ by
$u(t_n)$ and the corresponding numerical value obtained by TDG/FEM by
$u_n^-$. Let the error in $u_n^-$ be denoted by:

$$
e(t_n):=u_n^- - u(t_n)
$${#eq-ch3-30}

To motivate the appropriate notion of stability for the ease under
consideration, let us investigate the behavior of the homogeneous form
of @eq-ch3-1. The exact solutions of the homogeneous ODE is given by
[@Hughes2012 Chapter 8]

$$
u\left( {{t_n}} \right) = {u_0}{\operatorname{e} ^{ - \lambda {t_n}}}
$${#eq-ch3-31}

At time $t=t_{n+1}$ the exact solutions can be described as,

$$
            u\left( {{t_{n + 1}}} \right) = u\left( {{t_n}} \right){\operatorname{e} ^{ - \lambda \Delta {t_n}}}
$${#eq-ch3-32}

In @eq-ch3-32 one should note that the exact solutions decay in time for $\lambda>0$, and mathematically the solutions can be characterized by:

$$
\left.\begin{array}{cc}
\left|{u\left({{t_{n+1}}}\right)}\right|<\left|{u\left({{t_{n}}}\right)}\right|, & {\lambda>0}\\
{u\left({{t_{n+1}}}\right)=u\left({{t_{n}}}\right),} & {\lambda=0}
\end{array}\right\}
$${#eq-ch3-33}

In homogeneous case @eq-ch3-29 becomes

becomes

$$
\frac{1}{2}\left[\begin{array}{cc}
1 & 1\\
{-1} & 1
\end{array}\right]\left\{ \begin{array}{c}
{u_{n}^{+}}\\
{u_{n+1}^{-}}
\end{array}\right\} +\lambda\frac{{\Delta{t_{n}}}}{6}\left[\begin{array}{c}
{u_{n}^{+}}\\
{u_{n+1}^{-}}
\end{array}\right]\left\{ \begin{array}{c}
{u_{n}^{+}}\\
{u_{n+1}^{-}}
\end{array}\right\} =\left\{ \begin{array}{c}
{u_{n}^{-}}\\
0
\end{array}\right\}
$${#eq-ch3-34}

By eliminating $u_{n}^{+}$ from above
equations we can get

$$
u_{n + 1}^ -  = Au_n^{-}
$$ {#eq-ch3-35}

where,

$$
A = \frac{{6- 2\Omega }}{{{\Omega ^2} + 4\Omega  + 6}}
$${#eq-ch3-36}

is called the \emph{amplification factor} and $\Omega = \lambda \Delta t_n$.

For the stability of the TDG/FEM it is necessary that

$$
\left. {\begin{array}{cc}
{\left| {u_{n + 1}^ - } \right| < \left| {u_n^ - } \right|,}&{\lambda  > 0} \\
{u_{n + 1}^ -  = u_n^ - ,}&{\lambda  = 0}
\end{array}} \right\}
$${#eq-ch3-37}

From the definition of $A$ (see @eq-ch3-36), second condition of @eq-ch3-37 is satisfied.
The first condition is always satisfied if amplification factor satisfies the following,

$$
\left| A \right| < 1
$${#eq-ch3-38}

@fig-ch3-4 plots the variation of amplification factor with the dimensionless time $\Omega=\lambda \Delta t_n$, where it is evident that the TDG/FEM algorithm always satisfies the above-mentioned condition (cf. @eq-ch3-38). Thus, it is proved that the TDG Thus, it is proved that the TDG/FEM is an
*unconditionally stable* algorithm.

![Amplification factor for the time-discontinuous Galerkin method for
the first order ODE with linear interpolation in
time.](./figures/ch3-fig-4.svg){#fig-ch3-4}

## Convergence analysis of the TDG/FEM

TDG/FEM algorithm will be called convergent if for $t_n$ fixed,
$u_n^- \rightarrow u(t_n)$ as $\Delta t \rightarrow 0$. To establish the
convergence of an algorithm, two additional notion must be considered:
*stability* and *consistency*. In the previous section it has been shown
that the TDG/FEM algorithm is unconditionally stable. Therefore, in this
section, the consistency of the algorithm will be proven.

Rewriting the temporally discrete model problem described by @eq-ch3-35 in the form,

$$
u_{n + 1}^{-}  - Au_{n}^{-}  = 0
$${#eq-ch3-39}

Replacing $u_n^-$ and $u_{n+1}^-$ with the exact values $u(t_n)$ and $u(t_{n+1})$,
respectively, following expression is obtained

$$
u(t_{n+1}) - Au(t_n)=\Delta t_n  \cdot \tau(t_n)
$${#eq-ch3-40}

where, $\tau(t_n)$ is called the *local truncation error*.

:::{.callout-note appearance="simple"}
**Definition 1**. The algorithm defined by @eq-ch3-39 is called *consistent* if $\left| {\tau (t)} \right| \leqslant c\Delta {t^k}$, for all $t \in \left[ 0, T\right]$, where $c$ is a constant independent of $\Delta t$, and $k>0$. Moreover, $k$ is called the *order of accuracy* or *rate of convergence*.
:::

To show that TDG/FEM algorithm is consistent use Taylor expansion for
$u(t_{n+1})$ about $t_n$,

$$
u\left( {{t_{n + 1}}} \right) = u\left( {{t_n}} \right) + \Delta {t_n}\frac{{d{u_n}}}{{dt}} + \frac{{\Delta t_n^2}}{2}\frac{{{d^2}{u_n}}}{{d{t^2}}} + \frac{{\Delta t_n^3}}{6}\frac{{{d^3}{u_n}}}{{d{t^3}}} + O\left( {\Delta {t^2}} \right)
$${eq-ch3-41}

using the homogeneous form of model equation to eliminate the time derivatives in above equation one can get,

$$
u\left( {{t_{n + 1}}} \right) = \left( {1 - \Omega  + \frac{{{\Omega ^2}}}{2} - \frac{{{\Omega ^3}}}{6}} \right)u\left( {{t_n}} \right) + O\left( {\Delta {t^2}} \right)
$${#eq-ch3-42}

From @eq-ch3-36, @eq-ch3-40 and @eq-ch3-42, it follows that,

$$
\left| {\tau (t)} \right| \leqslant \frac{1}{{36}}\Delta t_n^3\qquad \forall t \in \left[ {0,T} \right]
$${#eq-ch3-43}

This completes the proof that TDG/FEM algorithm is consistent.

:::{.callout-note title="Remark"}
*Remark 1*. Therefore, the TDG/FEM algorithm is consistent and
unconditionally stable which also proves the convergence of the
algorithm according to the Lax equivalence theorem. [^6]
:::

:::{.callout-note title="Remark"}
*Remark 2*. The TDG/FEM is a third order accurate algorithm for the
linear interpolation in time. Moreover, the scheme can be classified as
a single-step algorithm.
:::
