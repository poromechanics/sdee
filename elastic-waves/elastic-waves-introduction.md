# Introduction {#sec-elastic-waves-introduction}

In many fields of both engineering and physics the problems ranging from
simulation of earthquake ground motion [@Bao1998; @Bielak2003] and
dynamic soil-structure interaction [@Wolf1985; @Wolf1988], to
electromagnetic waves [@Chew1995], and quantum mechanics [@Alonso2003]
may be best represented or modeled by considering the linear wave
propagation on infinite or semi-infinite unbounded domain. Such modeling
is especially of interest in the design of earth-structures such as dam,
tunnels, embankments, hospital and residential buildings etc. against
the transient loading and vibrations caused by the high-speed trains,
road-traffic, underground explosions, and more importantly earthquake
motion. In these problems the finite dimensional structure dynamically
interacts with the adjacent unbounded soil domain, and therefore the two
domains mutually influence the dynamic responses of each other
[@Wolf1996; @Burman2012].

For the computation of the dynamic SSI problems, a surface called
*interaction horizon*  [@Wolf1996] that encloses the structure and forms
the boundary of the computational domain has to be selected. The
objective of these domain reduction techniques is twofold; decreasing
the computation burden by reducing the size of the problem and enforcing
the so called radiation condition thus prohibiting any spurious
reflections at the artificial truncated boundaries (ATB). In this way,
ATB simulate the effects of the far field on the dynamic response of
both the structure and the near field.

The structure of the chapter is as follows. In section 2 of this chapter
basic theory of wave propagation in elastic solids is given. Section 3
of this chapter reviews some of the most popular boundary conditions for
solving wave propagation problems in unbounded domains. Section 4 of
this chapter provides the derivation of viscous boundary condition first
proposed by [@Lysmer1969] and discusses its characteristics.
Subsequently, in section 5 of this chapter, viscous boundary conditions
are modified to allow seismic excitation to enter the computational
domain.

