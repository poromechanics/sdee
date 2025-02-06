# TDG/FEM for the second order ODE {#sec-ch3-sec3}

Consider a mass-spring-dashpot system as depicted in @fig-ch3-fig-5 . The governing equation of motion is described by the following second order initial value problem in time.

$$
\begin{split}\frac{{{d^{2}}u}}{{d{t^{2}}}}+2\zeta{\omega_{n}}\frac{{du}}{{dt}}+\omega_{n}^{2}u & =f\left(t\right)\quad\forall t\in\left[0,T\right]\\
u(0) & =u_{0}\\
\frac{du(0)}{dt} & =v_{0}
\end{split}
$${#eq-ch3-eq-44}

where $u:=u(t)$ is the unknown displacement, $f(t)$ is the external force acting on the system. Further, $u_{0}$ and $v_{0}$ are the prescribed initial values of the displacement and velocity, respectively. Damping ratio $\zeta$ and the natural frequency of vibration $\omega_{n}$ of the system are related to the mass $m$, stiffness of the spring $k$, and damping coefficient $c$ by:

$$
\begin{aligned}
\omega_{n} & =\sqrt{k/m}, & \zeta & =\frac{c}{2m\omega_{n}}=\frac{c}{2\sqrt{mk}}
\end{aligned}
$${#eq-ch3-eq-45}

![Schematic diagram of the mass-spring-dashpot system](figures/ch3-fig-5.pdf){#fig-ch3-fig-5}

Further, the TDG/FEM for solving the second order ODE can be arranged into two categories; the displacement-velocity based two-field TDG/FEM, and the single-field TDG/FEM. These strategies are discussed in the following sections.

## Two-field TDG/FEM

In the two-field TDG/FEM, both the displacement ($u$) and the velocity ($v$) are independently interpolated using the piecewise polynomials. [^7] The interpolation is performed such that within a time-slab $I_{n}$ the solutions remain continuous, and at the end-points (i.e., $t_{n}$ and $t_{n+1}$) the solutions are discontinuous. Accordingly, the set $\left\{ t_{1},t_{2},\cdots,t_{N}\right\}$ denotes the locations in time where discontinuity in solutions occur. Thus, the displacement and velocity are the primary unknowns in the uv-TDG/FEM.

Therefore, the critical step in solving the second order initial value problem using the uv-TDG/FEM involves recasting of @eq-ch3-eq-44 into a system of two first-order ODEs. The new system is then described by

$$
\frac{{dv}}{{dt}}+2\zeta{\omega_{n}}v+\omega_{n}^{2}u=f\left(t\right)\quad\forall t\in[0,T]
$${#eq-ch3-eq-46}

$$
\frac{{du}}{{dt}}-v=0\quad\forall t\in[0,T]
$${#eq-ch3-eq-47}

$$
u(0)={u_{0}},\quad v(0)={v_{0}}
$${#eq-ch3-eq-48}

@eq-ch3-eq-46 and @eq-ch3-eq-47 represent the first order ODE, and the first order ODE, and the TDG/FEM described in the previous section can be employed directly. However, note that these two equations cannot be solved independently due to the coupling between displacement and velocity.

Following the same procedure as described in the previous section, the weak-form of the uv-TDG/FEM can be stated as:

:::{.callout-note title="Weak form"}
*Weak-form 1*. Find $u^{h}\in\Im_{l}^{h}$ and $v^{h}\in\Im_{l}^{h}$,
such that for all $\delta u^{h}\in\Im_{l}^{h}$ and
$\delta v^{h}\in\Im_{l}^{h}$, and for all $n=0,\cdots,N-1$ @eq-ch3-eq-49 holds.

$$
\begin{split} & \int_{{I_{n}}}^ {}{\delta{v^{h}}}\left({\frac{{d{v^{h}}}}{{dt}}+2\zeta{\omega_{n}}{v^{h}}+\omega_{n}^{2}{u^{h}}-f\left(t\right)}\right)dt+\delta{v^{h}}\left({{t_{n}}}\right){\left[\kern-0.15em \left[{{v^{h}}}\right]\kern-0.15em \right]_{n}}\\
 & +\int_{{I_{n}}}^ {}{\delta{u^{h}}\left({\frac{{d{u^{h}}}}{{dt}}-{v^{h}}}\right)dt+\delta{u^{h}}\left({{t_{n}}}\right){{\left[\kern-0.15em \left[{{u^{h}}}\right]\kern-0.15em \right]}_{n}}=0}
\end{split}
$${#eq-ch3-eq-49}
:::

A careful examination of uv-TDG/FEM weak-form leads to the following remarks.

:::{#rem-3}
In the above weak-form, the presence of jump discontinuity in time for the displacement, ${{{\left[\kern-0.15em \left[{{u^{h}}}\right]\kern-0.15em \right]}_{n}}}$, and for the velocity, ${{{\left[\kern-0.15em \left[{{v^{h}}}\right]\kern-0.15em \right]}_{n}}}$, correspond to the weakly enforced initial condition for the displacement and velocity, respectively.
:::

:::{#rem-4}
Since the selection of the test functions, $\delta u^{h}$ and $\delta v^{h}$, are independent from each other @eq-ch3-eq-49 can be depicted by the combination of following two variational forms.

$$
\int_{{I_{n}}}^ {}{\delta{v^{h}}}\left({\frac{{d{v^{h}}}}{{dt}}+2\zeta{\omega_{n}}{v^{h}}+\omega_{n}^{2}{u^{h}}-f\left(t\right)}\right)dt+\delta{v^{h}}\left({{t_{n}}}\right){\left[\kern-0.15em \left[{{v^{h}}}\right]\kern-0.15em \right]_{n}}=0
$${#eq-ch3-eq-50}

$$
\int_{{I_{n}}}^ {}{\delta{u^{h}}\left({\frac{{d{u^{h}}}}{{dt}}-{v^{h}}}\right)dt+\delta{u^{h}}\left({{t_{n}}}\right){{\left[\kern-0.15em \left[{{u^{h}}}\right]\kern-0.15em \right]}_{n}}=0}
$${#eq-ch3-eq-51}

From @eq-ch3-eq-51 it follows that in two-field TDG/FEM, the
displacement-velocity compatibility relationship is satisfied in weak
form.
:::

:::{#rem-5}
It is of course possible to use the different order interpolation for the displacement and velocity in the above weak-form. Only equal order interpolations, however, yield useful and efficient algorithms [@Hulbert1992].
:::

Let us now focus on the discretization of the two-field TDG weak-form. The discretization will be performed by using the locally defined $p$-order test and trial functions of the form,

$$
\begin{aligned}
{u^{h}} & =\sum\limits_{a=1}^{p+1}{T_{a}^{\left(p\right)}u_{a}^{(n)}}, & \delta{u^{h}} & =\sum\limits_{a=1}^{p+1}{T_{a}^{\left(p\right)}\delta u_{a}^{(n)}}, & \forall t\in I_{n}
\end{aligned}
$${#eq-ch3-eq-52}

$$
\begin{aligned}
{v^{h}} & =\sum\limits_{a=1}^{p+1}{T_{a}^{\left(p\right)}v_{a}^{(n)}}, & \delta{v^{h}} & =\sum\limits_{a=1}^{p+1}{T_{a}^{\left(p\right)}\delta v_{a}^{(n)}}, & \forall t\in I_{n}
\end{aligned}
$${#eq-ch3-eq-53}

where ${T_{a}^{\left(p\right)}}$ are the $p$-order Lagrange polynomials, and given by @eq-ch3-eq-16. Besides, in above equation, following conventions have been used.

$$
\begin{aligned}
u_{1}^{(n)} & =u_{n}^{+}, & u_{2}^{(n)} & =u_{n+1}^{-}, & v_{1}^{(n)} & =v_{n}^{+}, & v_{2}^{(n)} & =v_{n+1}^{-}
\end{aligned}
$$

where, $(u_{n}^{+},u_{n+1}^{-})$ and $(v_{n}^{+},v_{n+1}^{-})$ denote the discontinuous values of the displacement and velocity, respectively.

Subsequently, using the test functions and trial function in the weak-form (@eq-ch3-eq-49) to obtain the following discretized form.

$$
\begin{split} & \delta v_{a}^{(n)}\left[{\int_{{I_{n}}}^ {}{T_{a}^{(p)}}\left({\frac{{dT_{b}^{(p)}}}{{dt}}}\right)dt+{\delta_{1a}}{\delta_{1b}}}\right]v_{b}^{(n)}+\delta v_{a}^{(n)}\left[{2\zeta{\omega_{n}}\int_{{I_{n}}}^ {}{T_{a}^{(p)}}T_{b}^{(p)}dt}\right]v_{b}^{(n)}\\
 & +\delta v_{a}^{(n)}\left[{\omega_{n}^{2}\int_{{I_{n}}}^ {}{T_{a}^{(p)}}T_{b}^{(p)}dt}\right]u_{b}^{(n)}-\delta v_{a}^{(n)}\left\{ {\int_{{I_{n}}}^ {}{T_{a}^{(p)}}f(t)dt}\right\} -\delta v_{a}^{(n)}\left\{ {{\delta_{1a}}v_{n}^{-}}\right\} \\
 & \delta u_{a}^{(n)}\left[{\int_{{I_{n}}}^ {}{T_{a}^{(p)}\frac{{dT_{b}^{(p)}}}{{dt}}dt}+{\delta_{1a}}{\delta_{1b}}}\right]u_{b}^{(n)}\\
 & -\delta u_{a}^{(n)}\left[{\int_{{I_{n}}}^ {}{T_{a}^{(p)}T_{b}^{(p)}dt}}\right]v_{b}^{(n)}-\delta u_{a}^{(n)}\left\{ {{\delta_{1a}}u_{n}^{-}}\right\} =0
\end{split}
$${#eq-ch3-eq-54}

Since @eq-ch3-eq-54 is true for all $\delta u_{a}^{(n)}$ and $\delta v_{a}^{(n)}$, one can get the following system of $2p+2$ number of algebraic equations.

$$
{\left[m\right]^{ab}}v_{b}^{(n)}+2\zeta{\omega_{n}}{\left[c\right]^{ab}}v_{b}^{(n)}+\omega_{n}^{2}{\left[c\right]^{ab}}u_{b}^{(n)}={\left\{ {{J_{ext}}}\right\} ^{a}}+{\left\{ {J_{0}^{v}}\right\} ^{a}}
$${#eq-ch3-eq-55}

$$
{\left[m\right]^{ab}}u_{b}^{(n)}-{\left[c\right]^{ab}}v_{b}^{(n)}={\left\{ {J_{0}^{u}}\right\} ^{a}}
$${#eq-ch3-eq-56}

The matrix-vector form of @eq-ch3-eq-55 and @eq-ch3-eq-56 is given by @eq-ch3-eq-57 and @eq-ch3-eq-58, respectively.

$$
\left[{\mathbf{m}}\right]\left\{ {{\mathbf{\tilde{v}}}}\right\} +2\zeta{\omega_{n}}\left[{\mathbf{c}}\right]\left\{ {{\mathbf{\tilde{v}}}}\right\} +\omega_{n}^{2}\left[{\mathbf{c}}\right]\left\{ {{\mathbf{\tilde{u}}}}\right\} =\left\{ {{{\mathbf{J}}_{ext}}}\right\} +\left\{ {{\mathbf{J}}_{0}^{v}}\right\}
$${#eq-ch3-eq-57}

$$
\left[{\mathbf{m}}\right]\left\{ {{\mathbf{\tilde{u}}}}\right\} -\left[{\mathbf{c}}\right]\left\{ {{\mathbf{\tilde{v}}}}\right\} =\left\{ {{\mathbf{J}}_{0}^{u}}\right\}
$${#eq-ch3-eq-58}

In @eq-ch3-eq-55 -- @eq-ch3-eq-58, the matrices $\left[{\mathbf{m}}\right]$,
$\left[{\mathbf{c}}\right]$, and the vector
$\left\{ {{{\mathbf{J}}_{ext}}}\right\}$ are given by @eq-ch3-eq-24, @eq-ch3-eq-25, and @eq-ch3-eq-27, respectively. The vectors, $\left\{ {{\mathbf{J}}_{0}^{v}}\right\}$ (see Eq. @eq-ch3-eq-59) and $\left\{ {{\mathbf{J}}_{0}^{u}}\right\}$ (see Eq. @eq-ch3-eq-60) correspond to the initial value of the velocity ($v_{n}^{-}$) and displacement ($u_{n}^{-}$), respectively.

$$
\left\{ {{\mathbf{J}}_{0}^{v}}\right\} :={\left\{ {J_{0}^{v}}\right\} ^{a}}={\delta_{1a}}v_{n}^{-}
$${#eq-ch3-eq-59}

$$
\left\{ {{\mathbf{J}}_{0}^{u}}\right\} :={\left\{ {J_{0}^{u}}\right\} ^{a}}={\delta_{1a}}u_{n}^{-}
$${#eq-ch3-eq-60}

:::{#rem-6}
In any time-slab $I_{n}$, there are $p+1$ unknowns for the
velocity and $p+1$ unknowns for the displacement. Consequently, there
are total $2p+2$ unknowns to be determined in each time-slab. These
unknowns are computed by solving the system of $2p+2$ equations formed
by @eq-ch3-eq-57 and @eq-ch3-eq-58. Besides, in Eq. @eq-ch3-eq-57 -- @eq-ch3-eq-58, the shape of all matrices and all vectors are $(p+1)\times(p+1)$ and $(p+1)\times(1)$.
:::

If the displacement and the velocity are linearly interpolated in time
(i.e., $p=1$) using the shape function described by @eq-ch3-eq-28, then @eq-ch3-eq-57 and @eq-ch3-eq-58 can be written as follows.

$$
\begin{split}\frac{1}{2}\left[\begin{array}{cc}
1 & 1\\
{-1} & 1
\end{array}\right]\left\{ \begin{array}{c}
{v_{n}^{+}}\\
{v_{n+1}^{-}}
\end{array}\right\} +\zeta{\omega_{n}}\frac{{\Delta t}}{3}\left[\begin{array}{cc}
2 & 1\\
1 & 2
\end{array}\right]\left\{ \begin{array}{c}
{v_{n}^{+}}\\
{v_{n+1}^{-}}
\end{array}\right\} \\
+\omega_{n}^{2}\frac{{\Delta t}}{6}\left[\begin{array}{cc}
2 & 1\\
1 & 2
\end{array}\right]\left\{ \begin{array}{c}
{u_{n}^{+}}\\
{u_{n+1}^{-}}
\end{array}\right\} =\left\{ \begin{array}{c}
{v_{n}^{-}}\\
0
\end{array}\right\} +\left\{ \begin{array}{c}
{J_{ext}^{1}}\\
{J_{ext}^{2}}
\end{array}\right\}
\end{split}
$${#eq-ch3-eq-61}

$$
\frac{1}{2}\left[\begin{array}{cc}
1 & 1\\
{-1} & 1
\end{array}\right]\left\{ \begin{array}{c}
{u_{n}^{+}}\\
{u_{n+1}^{-}}
\end{array}\right\} -\frac{{\Delta t_{n}}}{6}\left[\begin{array}{cc}
2 & 1\\
1 & 2
\end{array}\right]\left\{ \begin{array}{c}
{v_{n}^{+}}\\
{v_{n+1}^{-}}
\end{array}\right\} =\left\{ \begin{array}{c}
{u_{n}^{-}}\\
0
\end{array}\right\}
$${#eq-ch3-eq-62}

In @eq-ch3-eq-61, the expressions for $J_{ext}^{1}$ and $J_{ext}^{2}$ are identical to those given in @eq-ch3-eq-29.

## Displacement based single-field TDG/FEM

In the displacement based single-field TDG/FEM (u-TDG/FEM) only displacement is interpolated using the piecewise polynomials. The displacement remains continuous within a time-slab $I_{n}$. However, at the end-points (i.e., $t_{n}$ and $t_{n+1}$) displacement takes two different values, for example, $u_{n}^{+}$ and $u_{n}^{-}$ at $t_{n}$. In addition, the velocity is obtained by taking the time derivative of the displacement. Thus, the velocity-displacement compatibility relationship is naturally satisfied, and @eq-ch3-eq-47 is no longer required to be solved. However, note that both displacement and velocity still remain discontinuous in time.

The weak form of the displacement based single-field TDG/FEM, which is described below (see @eq-ch3-eq-63), is obtained by considering the second order ODE (@eq-ch3-eq-44).

:::{.callout-note title="Weakform"}
*Weak-form 2*. u-TDG/FEM: Find $u^{h}\in\Im_{l}^{h}$ such that for all
$\delta u^{h}\in\Im_{l}^{h}$, and for all $n=0,\cdots,N-1$ @eq-ch3-eq-63 holds.

$$
\begin{split}\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}\left({\frac{{{d^{2}}{u^{h}}}}{{d{t^{2}}}}+2\zeta{\omega_{n}}\frac{{d{u^{h}}}}{{dt}}+\omega_{n}^{2}{u^{h}}-f(t)}\right)dt}\\
+\frac{{d\delta{u^{h}}({t_{n}})}}{{dt}}{\left[\kern-0.15em \left[{\frac{{d{u^{h}}}}{{dt}}}\right]\kern-0.15em \right]_{n}}+\delta{u^{h}}({t_{n}})\omega_{n}^{2}{\left[\kern-0.15em \left[{{u^{h}}}\right]\kern-0.15em \right]_{n}}=0
\end{split}
$${#eq-ch3-eq-63}
:::

The above weak form is obtained by using the following intermediate results

$$
\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}\frac{{{d^{2}}{u^{h}}}}{{d{t^{2}}}}dt}=\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}\frac{{{d^{2}}{u^{h}}}}{{d{t^{2}}}}dt}+\left[{\frac{{d\delta{u^{h}}}}{{dt}}\left({{{\left({\frac{{du}}{{dt}}}\right)}^{*}}-\frac{{d{u^{h}}}}{{dt}}}\right)}\right]_{{t_{n}}}^{{t_{n+1}}}
$${#eq-ch3-eq-64}

$$
\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}}\omega_{n}^{2}{u^{h}}dt=\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}}\omega_{n}^{2}{u^{h}}dt+\left[{\delta{u^{h}}\omega_{n}^{2}\left({{u^{*}}-{u^{h}}}\right)}\right]_{{t_{n}}}^{{t_{n+1}}}
$${#eq-ch3-eq-65}

These intermediate results are obtained by following the procedure described in  the previous section. Further, in above equations, ${{u^{*}}}$ and ${{{\left({\frac{{du}}{{dt}}}\right)}^{*}}}$ denote the unique representative value of displacement and its first time derivative at the end points of $I_{n}$ (here, recall that at $t_{n}$ and $t_{n+1}$ both displacement and its first time derivative are discontinuous). Furthermore, by adopting the definition of representative values, which is provided in @eq-ch3-eq-10, @eq-ch3-eq-64 and @eq-ch3-eq-65 transform into

$$
\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}\frac{{{d^{2}}{u^{h}}}}{{d{t^{2}}}}dt}=\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}\frac{{{d^{2}}{u^{h}}}}{{d{t^{2}}}}dt}+\frac{{d\delta{u^{h}}\left({{t_{n}}}\right)}}{{dt}}{\left[\kern-0.15em \left[{\frac{{d{u^{h}}}}{{dt}}}\right]\kern-0.15em \right]_{n}}
$${#eq-ch3-eq-66}

$$
\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}}\omega_{n}^{2}{u^{h}}dt=\int_{{I_{n}}}^ {}{\frac{{d\delta{u^{h}}}}{{dt}}}\omega_{n}^{2}{u^{h}}dt+\delta{u^{h}}({t_{n}})\omega_{n}^{2}{\left[\kern-0.15em \left[{{u^{h}}}\right]\kern-0.15em \right]_{n}}
$${#eq-ch3-eq-67}

Let us now focus on the discretization of the u-TDG/FEM weak-form given in @eq-ch3-eq-63. Here, for the sake of clarity, the discretization will be performed by using the locally defined quadratic test and trial functions of the form 

$$
\begin{aligned}
{u^{h}} & =u_{1}^{(n)}{T_{1}}+u_{2}^{(n)}{T_{2}}+u_{3}^{(n)}{T_{3}}, & \delta{u^{h}} & =\delta u_{1}^{(n)}{T_{1}}+\delta u_{2}^{(n)}{T_{2}}+\delta u_{3}^{(n)}{T_{3}}
\end{aligned}
$${#eq-ch3-eq-68}

where, $T_{1}$, $T_{2}$, and $T_{3}$ are the quadratic shape functions described as follows. 

$$
\begin{aligned}
{T_{1}} & =\frac{1}{2}\left({{\theta^{2}}-\theta}\right), & {T_{2}} & =\frac{1}{2}\left({{\theta^{2}}+\theta}\right), & {T_{3}} & =1-{\theta^{2}}
\end{aligned}
$${#eq-ch3-eq-69}

In this case, the time derivatives of the test function and trial functions are given by following expressions.

$$
\frac{{d{u^{h}}}}{{dt}}=\frac{2}{{\Delta{t_{n}}}}\left\{ {\frac{1}{2}\left({2\theta-1}\right)u_{1}^{(n)}+\frac{1}{2}\left({2\theta+1}\right)u_{2}^{(n)}-2\theta u_{3}^{(n)}}\right\} 
$${#eq-ch3-eq-70}

$$
\frac{{{d^{2}}{u^{h}}}}{{d{t^{2}}}}=\frac{4}{{\Delta t_{n}^{2}}}\left\{ {u_{1}^{(n)}+u_{2}^{(n)}-2u_{3}^{(n)}}\right\} 
$${#eq-ch3-eq-71}

$$
\frac{{d\delta{u^{h}}}}{{dt}}=\frac{2}{{\Delta{t_{n}}}}\left\{ {\frac{1}{2}\left({2\theta-1}\right)\delta u_{1}^{(n)}+\frac{1}{2}\left({2\theta+1}\right)\delta u_{2}^{(n)}-2\theta\delta u_{3}^{(n)}}\right\}
$${#eq-ch3-eq-72}

$$
\frac{{d{u^{h}}\left({{t_{n}}}\right)}}{{dt}}=\frac{2}{{\Delta{t_{n}}}}\left\{ {-\frac{3}{2}u_{1}^{(n)}-\frac{1}{2}u_{2}^{(n)}+2u_{3}^{(n)}}\right\}
$${#eq-ch3-eq-73}

$$
\frac{{d\delta{u^{h}}\left({{t_{n}}}\right)}}{{dt}}=\frac{2}{{\Delta{t_{n}}}}\left\{ {-\frac{3}{2}\delta u_{1}^{(n)}-\frac{1}{2}\delta u_{2}^{(n)}+2\delta u_{3}^{(n)}}\right\} 
$${#eq-ch3-eq-74}

Subsequently, using the above expressions in the u-TDG/FEM weak-form following system of linear equation can be obtained.

$$
\begin{split}\frac{1}{{\Delta t_{n}^{2}}}\left[\begin{array}{ccc}
5 & {-1} & {-4}\\
7 & 5 & {-12}\\
{-12} & {-4} & {16}
\end{array}\right]\left\{ \begin{array}{c}
{u_{1}^{(n)}}\\
{u_{2}^{(n)}}\\
{u_{3}^{(n)}}
\end{array}\right\} +\frac{{2\zeta{\omega_{n}}}}{{3\Delta{t_{n}}}}\left[\begin{array}{ccc}
7 & 1 & {-8}\\
1 & 7 & {-8}\\
{-8} & {-8} & {16}
\end{array}\right]\left\{ \begin{array}{c}
{u_{1}^{(n)}}\\
{u_{2}^{(n)}}\\
{u_{3}^{(n)}}
\end{array}\right\} \\
+\frac{{\omega_{n}^{2}}}{6}\left[\begin{array}{ccc}
3 & 1 & {-4}\\
{-1} & 3 & 4\\
4 & {-4} & 0
\end{array}\right]\left\{ \begin{array}{c}
{u_{1}^{(n)}}\\
{u_{2}^{(n)}}\\
{u_{3}^{(n)}}
\end{array}\right\} =\frac{1}{{\Delta{t_{n}}}}\left\{ \begin{array}{c}
{-3v_{n}^{-}}\\
{-v_{n}^{-}}\\
{4v_{n}^{-}}
\end{array}\right\} +\left\{ \begin{array}{c}
{\omega_{n}^{2}u_{n}^{-}}\\
0\\
0
\end{array}\right\} +\left\{ \begin{array}{c}
{J_{ext}^{1}}\\
{J_{ext}^{2}}\\
{J_{ext}^{3}}
\end{array}\right\}
\end{split}
$${#eq-ch3-eq-75}

where $u_{1}^{(n)}:=u_{n}^{+}$, $u_{2}^{(n)}:=u_{n+1}^{-}$, and 

$$
\begin{split}J_{ext}^{1}=\int_{-1}^{+1}{\left({\frac{{2\theta-1}}{2}}\right)f(t)d\theta}\\
J_{ext}^{2}=\int_{-1}^{+1}{\left({\frac{{2\theta+1}}{2}}\right)f(t)d\theta}\\
J_{ext}^{3}=\int_{-1}^{+1}{\left({-2\theta}\right)f(t)d\theta}
\end{split}
$${#eq-ch3-eq-76}

:::{#rem-7}
In the case of u-TDG/FEM, the order of interpolation for the displacement should be at-least two (i.e., $p=2$). This requirement is due to the presence of second order time derivative in the weak-form (see @eq-ch3-eq-63).
:::

:::{#rem-8}
In case of quadratic interpolation, there are total three unknowns, $u_{a}^{(n)}$ for $a=1,2,3$, to be determined in each time-slab $I_{n}$. Besides, the number of unknowns in case of the u-TDG/FEM is less than that of uv-TDG/FEM.
:::

## Velocity based single-field TDG/FEM

From the displacement based single-field TDG/FEM it follows that the number of unknowns may be decreased by explicitly satisfying the displacement-velocity compatibility condition (cf. @eq-ch3-eq-47). In the u-TDG/FEM, however, the requirement of at-least quadratic time interpolation of the displacement implies that minimum three unknowns should be determined for each time-slab $I_{n}$. One of the objectives of this thesis is to further decrease the number of unknowns in a given time-slab. To achieve this goal a velocity based single-field TDG/FEM (henceforth, v-TDG/FEM) is developed.

The key idea behind the v-TDG/FEM is to treat the velocity as the only
primary unknown. In $I_{n}$, the velocity is interpolated using the
Lagrange polynomials of degree $p$; the velocity remains continuous in
$I_{n}$, but discontinuity occurs at the end-points $t_{n},t_{n+1}$.
Further, a consistent time-integration of the velocity is performed as
post-processing step to compute the displacement. Thus, the
displacement-velocity compatibility relationship is naturally satisfied,
and Eq. [\[eq:ch3-eq-47\]](#eq:ch3-eq-47){reference-type="eqref"
reference="eq:ch3-eq-47"} is no longer required to be solved. Regarding
the v-TDG/FEM, it is worth mentioning that the displacement remains
continuous throughout the time domain $\left[0,T\right]$, whereas in the
uv-TDG/FEM and u-TDG/FEM displacement is discontinuous at the discrete
times $\left\{ t_{0},t_{1},\cdots,t_{N}\right\}$.

The weak form of the v-TDG/FEM, which is described below (see Eq.
[\[eq:ch3-eq-77\]](#eq:ch3-eq-77){reference-type="ref"
reference="eq:ch3-eq-77"}), is obtained by considering the first order
ODE given by Eq.
[\[eq:ch3-eq-46\]](#eq:ch3-eq-46){reference-type="eqref"
reference="eq:ch3-eq-46"}.

::: weakform
*Weak-form 3*. Find $v^{h}\in\Im_{l}^{h}$ such that for all
$\delta v^{h}\in\Im_{l}^{h}$, and for all $n=0,\cdots,N-1$
 Eq. [\[eq:ch3-eq-77\]](#eq:ch3-eq-77){reference-type="eqref"
reference="eq:ch3-eq-77"}  holds.
$$\int_{{I_{n}}}^ {}{\delta{v^{h}}\left({\frac{{d{v^{h}}}}{{dt}}+2\zeta{\omega_{n}}{v^{h}}+\omega_{n}^{2}{u^{h}}-f(t)}\right)dt}+\delta{v^{h}}({t_{n}}){\left[\kern-0.15em \left[{{v^{h}}}\right]\kern-0.15em \right]_{n}}=0\label{eq:ch3-eq-77}$$
:::

In the above weak-form displacement is computed by using the following
relationship.
$${u^{h}}(t)=u({t_{n}})+\int_{{t_{n}}}^{t}{{v^{h}}(\tau)d\tau}\label{eq:ch3-eq-78}$$

Similar to the previous sections, discretization of the v-TDG/FEM
weak-form can be performed by employing the locally defined $p$-order
test and trial functions for the velocity.

$$\begin{aligned}
{v^{h}} & =\sum\limits_{a=1}^{p+1}{T_{a}^{(p)}v_{a}^{(n)}} & \delta{v^{h}} & =\sum\limits_{a=1}^{p+1}{T_{a}^{(p)}\delta v_{a}^{(n)}}\label{eq:ch3-eq-79}
\end{aligned}$$ where $T_{a}^{(p)}$ are the $p$-order Lagrange
polynomials (see Eq.
[\[eq:ch3-eq-16\]](#eq:ch3-eq-16){reference-type="ref"
reference="eq:ch3-eq-16"}), and $v_{1}^{(n)}=v_{n}^{+}$,
$v_{2}^{(n)}=v_{n+1}^{-}$.

Consequently, the discrete form of the displacement-velocity
compatibility relationship, which is described below, can be obtained by
using the above-mentioned trial functions for the velocity in the Eq.
[\[eq:ch3-eq-78\]](#eq:ch3-eq-78){reference-type="eqref"
reference="eq:ch3-eq-78"}.
$$u^{h}(t)=u^{h}({t_{n}})+\sum\limits_{a=1}^{p+1}{\tilde{T}_{a}^{(p)}v_{a}^{(n)}}\label{eq:ch3-eq-80}$$
where
$$\tilde{T}_{a}^{(p)}=\frac{{\Delta{t_{n}}}}{2}\int_{-1}^{+1}{T_{a}^{(p)}d\theta}\label{eq:ch3-eq-81}$$
are $p+1$ order locally defined polynomials.

::: remark
*Remark 9*. In Eq.
[\[eq:ch3-eq-80\]](#eq:ch3-eq-80){reference-type="eqref"
reference="eq:ch3-eq-80"}, by the virtue of time integration of the
$p$-order Lagrange polynomials, the displacement are described by the
$p+1$-order local piecewise polynomials. In addition,
 Eq. [\[eq:ch3-eq-80\]](#eq:ch3-eq-80){reference-type="eqref"
reference="eq:ch3-eq-80"}  is equivalent to the following form.
$${u^{h}}(t)=\sum\limits_{a=1}^{p+2}{T_{a}^{(p+1)}u_{a}^{(n)}}\label{eq:ch3-eq-82}$$
where ${T_{a}^{(p+1)}}$ are the $p+1$ order Lagrange polynomials given
by Eq. [\[eq:ch3-eq-16\]](#eq:ch3-eq-16){reference-type="eqref"
reference="eq:ch3-eq-16"}. Besides, $u_{1}^{(n)}=u_{n}$ and
$u_{2}^{(n)}=u_{n+1}$ are the continuous value of the displacement at
time $t_{n}$ and $t_{n+1}$, respectively.
:::

After using Eq. [\[eq:ch3-eq-79\]](#eq:ch3-eq-79){reference-type="eqref"
reference="eq:ch3-eq-79"} and Eq.
[\[eq:ch3-eq-80\]](#eq:ch3-eq-80){reference-type="eqref"
reference="eq:ch3-eq-80"} in the v-TDG/FEM weak-form described by Eq.
[\[eq:ch3-eq-77\]](#eq:ch3-eq-77){reference-type="eqref"
reference="eq:ch3-eq-77"} one can obtain the following discrete form.
$$\begin{split}\delta v_{a}^{(n)}\left[{\int_{{I_{n}}}^ {}{T_{a}^{(p)}\frac{{dT_{b}^{(p)}}}{{dt}}dt}+{\delta_{1a}}{\delta_{1b}}}\right]v_{b}^{(n)}+\delta v_{a}^{(n)}\left[{2\zeta{\omega_{n}}\int_{{I_{n}}}^ {}{T_{a}^{(p)}T_{b}^{(p)}dt}}\right]v_{b}^{(n)}\\
+\delta v_{a}^{(n)}\left[{\omega_{n}^{2}\int_{{I_{n}}}^ {}{T_{a}^{(p)}\tilde{T}_{b}^{(p)}dt}}\right]v_{b}^{(n)}+\delta v_{a}^{(n)}\left\{ {\int_{{I_{n}}}^ {}{T_{a}^{(p)}dt}}\right\} \omega_{n}^{2}{u_{n}}\\
-\delta v_{a}^{(n)}\left\{ {\int_{{I_{n}}}^ {}{T_{a}^{(p)}f(t)dt}}\right\} -\delta v_{a}^{(n)}\left\{ {{\delta_{1a}}v_{n}^{-}}\right\} =0
\end{split}
\label{eq:ch3-eq-83}$$

Subsequently, using the fact that Eq.
[\[eq:ch3-eq-83\]](#eq:ch3-eq-83){reference-type="eqref"
reference="eq:ch3-eq-83"} is true for all $\delta v_{a}^{(n)}$ following
system $p+1$ algebraic equations in $p+1$ unknowns can be obtained.
$${\left[m\right]^{ab}}v_{b}^{(n)}+2\zeta{\omega_{n}}{\left[c\right]^{ab}}v_{b}^{(n)}+\omega_{n}^{2}{\left[k\right]^{ab}}v_{b}^{(n)}={\left\{ {{J_{ext}}}\right\} ^{a}}-{\left\{ {J_{0}^{u}}\right\} ^{a}}+{\left\{ {J_{0}^{v}}\right\} ^{a}}\label{eq:ch3-eq-84}$$
The matrix-vector form of above equation is described by
$$\left[{\mathbf{m}}\right]\left\{ {{\mathbf{\tilde{v}}}}\right\} +2\zeta{\omega_{n}}\left[{\mathbf{c}}\right]\left\{ {{\mathbf{\tilde{v}}}}\right\} +\omega_{n}^{2}\left[{\mathbf{k}}\right]\left\{ {{\mathbf{\tilde{v}}}}\right\} =\left\{ {{{\mathbf{J}}_{ext}}}\right\} -\left\{ {{\mathbf{J}}_{0}^{u}}\right\} +\left\{ {{\mathbf{J}}_{0}^{v}}\right\} \label{eq:ch3-eq-85}$$
where the matrices $\left[{\mathbf{m}}\right]$,
$\left[{\mathbf{c}}\right]$, and the vector
$\left\{ {{{\mathbf{J}}_{ext}}}\right\}$ are given by Eq.
[\[eq:ch3-eq-24\]](#eq:ch3-eq-24){reference-type="eqref"
reference="eq:ch3-eq-24"}, Eq.
[\[eq:ch3-eq-25\]](#eq:ch3-eq-25){reference-type="eqref"
reference="eq:ch3-eq-25"}, and Eq.
[\[eq:ch3-eq-27\]](#eq:ch3-eq-27){reference-type="eqref"
reference="eq:ch3-eq-27"}, respectively. The expression for the
$\left[{\mathbf{k}}\right]$ matrix is given below.
$$\left[{\mathbf{k}}\right]:={\left[k\right]^{ab}}=\int_{{I_{n}}}^ {}{T_{a}^{(p)}\tilde{T}_{b}^{(p)}dt}\label{eq:ch3-eq-86}$$

Furthermore, in Eq.
[\[eq:ch3-eq-85\]](#eq:ch3-eq-85){reference-type="eqref"
reference="eq:ch3-eq-85"}, the element vectors,
$\left\{ {{\mathbf{J}}_{0}^{v}}\right\}$ and
$\left\{ {{\mathbf{J}}_{0}^{u}}\right\}$ correspond to the initial value
of velocity $(v_{n}^{-})$ and displacement $(u_{n})$, respectively. The
expressions for these vectors are presented as follows.
$$\left\{ {{\mathbf{J}}_{0}^{v}}\right\} :={\left\{ {J_{0}^{v}}\right\} ^{a}}={\delta_{1a}}v_{n}^{-}\label{eq:ch3-eq-87}$$
$$\left\{ {{\mathbf{J}}_{0}^{u}}\right\} :={\left\{ {J_{0}^{u}}\right\} ^{a}}=\omega_{n}^{2}{u_{n}}\int_{{I_{n}}}^ {}{T_{a}^{(p)}dt}\label{eq:ch3-eq-88}$$

Let us now consider the special case when velocity is linearly
interpolated in time by using the trial functions of the form
$${v^{h}}(t)={T_{1}}v_{n}^{+}+{T_{2}}v_{n+1}^{-}\label{eq:ch3-eq-89}$$
where $T_{1}$ and $T_{2}$ are linear shape functions which are given in
Eq. [\[eq:ch3-eq-28\]](#eq:ch3-eq-28){reference-type="eqref"
reference="eq:ch3-eq-28"}. Accordingly, Eq.
[\[eq:ch3-eq-80\]](#eq:ch3-eq-80){reference-type="eqref"
reference="eq:ch3-eq-80"} transforms into
$${u^{h}}(t)={u_{n}}+v_{n}^{+}{{\tilde{T}}_{1}}+v_{n+1}^{-}{{\tilde{T}}_{2}}\label{eq:ch3-eq-90}$$
in which, $$\begin{aligned}
{{\tilde{T}}_{1}} & =\frac{\Delta t_{n}}{2}(1-T_{1}^{2}), & {{\tilde{T}}_{2}} & =\frac{\Delta t_{n}}{2}T_{2}^{2}\label{eq:ch3-eq-91}
\end{aligned}$$ are the quadratic polynomials. Accordingly, Eq.
[\[eq:ch3-eq-85\]](#eq:ch3-eq-85){reference-type="eqref"
reference="eq:ch3-eq-85"} now reads,
$$\begin{split}\frac{1}{2}\left[\begin{array}{cc}
1 & 1\\
{-1} & 1
\end{array}\right]\left\{ \begin{array}{c}
{v_{n}^{-}}\\
{v_{n+1}^{+}}
\end{array}\right\} +\frac{2\zeta{\omega_{n}}\Delta{t_{n}}}{6}\left[\begin{array}{cc}
2 & 1\\
1 & 2
\end{array}\right]\left\{ \begin{array}{c}
{v_{n}^{-}}\\
{v_{n+1}^{+}}
\end{array}\right\} \\
+\frac{{\omega_{n}^{2}\Delta t_{n}^{2}}}{{24}}\left[\begin{array}{cc}
3 & 1\\
5 & 3
\end{array}\right]\left\{ \begin{array}{c}
{v_{n}^{-}}\\
{v_{n+1}^{+}}
\end{array}\right\} =\left\{ \begin{array}{c}
{J_{ext}^{1}}\\
{J_{ext}^{2}}
\end{array}\right\} -\frac{\omega_{n}^{2}\Delta{t_{n}}{u_{n}}}{2}\left\{ \begin{array}{c}
1\\
1
\end{array}\right\} +\left\{ \begin{array}{c}
{v_{n}^{-}}\\
0
\end{array}\right\}
\end{split}
\label{eq:ch3-eq-92}$$
