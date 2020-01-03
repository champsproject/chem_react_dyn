# $n$-Degree-of-Freedom System-Bath Model

## Introduction and Development of the Problem
(REMOVE LD DESCRIPTION)

Isomerization is one of the most prevalent reaction when studying
atmospheric, medical, and industrial chemistry [^1], and thus have
garnered interest from both theoretical and applied
chemistry [@wieder_dissociation_1962; @mciver_structure_1972; @dugave_cistrans_2003; @donohoe_ruthenium-catalyzed_2009].
From an applied perspective, the influence of various solvents on the
rate constant of conformational isomerization has been pursued for
specific
molecules [@price_solvent_1961; @halicioglu_solvent_1969; @flom_dynamic_1986; @eberhardt_solvent_1992; @duffy_solvent_1992].
From a theoretical perspective, as is adopted in this chapter, it is
instructive to formulate a systematic approach for identifying reactive
trajectories in the spirit of transition path
sampling [@bolhuis_transition_2002] and reactive
islands [@davis_bottlenecks_1985; @davis_phase_1987; @marston_reactive_1989; @de_leon_order_1989].

**Modeling reaction in solvents**

Following the derivation of the Langevin
equations [@kramers_brownian_1940] and its generalized system dynamics
in a heat bath [@zwanzig_nonlinear_1973], the Kramers turnover of
reaction rates (escape rates from a potential well) with increasing
viscosity of the solvent (intensity of noise) has been used to obtain
reaction rate expressions in gas and condensed phase reactions by Grote,
Hynes, Pollak, Grabert,
Hänggi [@pollak_theory_1989; @pollak_transition-state_1995]. The
generalized Langevin equation of motion for a particle trapped in a
one-dimensional well with a barrier height $\mathcal{V}^{\ddagger}$ and
coupled to a medium such as gas or liquid is modeled by a bath with a
viscosity parameter. This framework has received much attention in the
literature where the system dynamics can be obtained explicitly. The
Langevin dynamics represents the motion of the particle in a bath which
is modeled using a large number of harmonic oscillators. In this set-up,
the bath coordinates are coupled with the system coordinates, and the
Hamiltonian is simply a sum of system and bath components. However, the
system dynamics of a two or more degrees of freedom coupled with bath
modes has not received a global analysis from a dynamical systems
perspective of reactions. In this direction, the first step would be to
consider a two degrees of freedom reaction with well-understood chemical
observables such as *reaction coordinate*, and where the harmonic bath
modes can be coupled to represent the reaction in a condensed phase. In
this work, we will adopt the two degrees of freedom isomerization model
of De Leon and Berne who have studied the chemical reaction dynamics
extensively using a dynamical systems
perspective [@de_leon_intramolecular_1981; @de_leon_order_1989; @de_leon_1991; @de_leon_1991-1; @de_leon_1992].
This model of isomerization describes the conformational change by the
motion of an internal angle where the isomers are represented by the
well in a double potential well separated by a barrier. The two degrees
of freedom in the model correspond to the bond that undergoes structural
change and to the bond that breaks above dissociation energy. Typically,
this isomerization barrier height is lower than the dissociation energy
of the molecule, and the activation energy for the isomerization is
imparted by molecular collisions or photoexcitation.

**Effect of solvent on reaction rates and role of phase space structures**

Traditionally, the construction of a dividing surface (DS) was focused
on critical points of the potential energy surface (PES), that is, in
the configuration space describing the molecular
system [@Komatsuzaki97; @Komatsuzaki00]. Critical points on the PES do
have significance in phase space; they are the equilibrium points for
zero momentum. But they continue to have influence for nonzero momentum
for a range of energies above the energy of the equilibrium point. The
construction of a DS separating the phase space into two parts,
reactants and products, has been a focus from the dynamical systems
point of view in recent
years [@wiggins_impenetrable_2001; @uzer_geometry_2002; @waalkens2004direct].
In phase space, that is for nonzero momentum, the role of the *saddle
point* is played by an *invariant manifold* of saddle stability type,
the normally hyperbolic invariant manifold
(NHIM) [@wiggins_geometry_1990; @wiggins2013normally]. In order to fully
appreciate the NHIM and its role in reaction rate theory, it is useful
to begin with a precursor concept the *periodic orbit dividing surface*
or PODS. For systems with two DoF described by a natural Hamiltonian,
kinetic plus potential energy, the problem of constructing the DS in
phase space was solved during the 1970s by McLafferty, Pechukas and
Pollak [@Pechukas73; @Pechukas77; @Pollak78; @Pechukas79]. They
demonstrated that the DS at a specific energy is related to an invariant
phase space structure, an unstable periodic orbit (UPO) which defines
(it is the boundary of) the bottleneck in phase space through which the
reaction occurs. The DS which intersects trajectories evolving from
reactants to products can then be shown to have the geometry of a
hemisphere in phase space whose boundary is the UPO. The same
construction can be carried out for a DS intersecting trajectories
crossing from products to reactants and these two hemispheres form a
sphere for which the UPO is the equator. Generalisation of this
construction of DS to high dimensional systems has been a central
question in reaction dynamics and has only received a satisfactory
answer in recent years
[@wiggins_impenetrable_2001; @uzer_geometry_2002; @waalkens2004direct].
The key difficulty being the high dimensional analogue of the UPO used
in the two DoF system for the construction of the DS and which is
resolved by considering the NHIM, which has the appropriate
dimensionality for anchoring the dividing surface in high dimensional
phase space [@wiggins_geometry_1990]. Normal hyperbolicity of these
invariant manifolds means that their stability, in a precise sense, is
of saddle type in the transverse direction, which implies that they
possess stable and unstable invariant manifolds that are impenetrable
barriers and mediate reactive trajectories in phase space. These
invariant manifolds of the NHIM are structurally stable, that is, stable
under perturbation [@wiggins2013normally]. For two DoF systems, the NHIM
is an unstable PO, and for an $N > 2$ DoF system at a fixed energy, the
NHIM has the topology of a $(2N-3)$-dimensional sphere and is the
equator of a $(2N-2)$-dimensional sphere which is the DS. This DS can
then be used to divide the $(2N-1)$-dimensional energy surface into two
parts, reactants and
products [@Gillilan91; @Komatsuzaki96; @Komatsuzaki97; @Komatsuzaki00; @Komatsuzaki02a].
An elementary description of the role of the NHIM in reaction dynamics
is given in [@wiggins_role_2016] along with description of their
geometry using quadratic normal form Hamiltonians. Fundamental theorems
assure the existence of the phase space structures  NHIM and its
invariant manifolds for a range of energies above that of the saddle
[@wiggins2013normally]. However, the precise extent of this range, as
well as the nature and consequences of any bifurcations of the phase
space structures that might occur as energy is increased, is not known
and is a topic of continuing
research[@Li_bifurcation_2009; @Inarrea_bifurcations_2011; @Allahem_chaotic_2012; @mauguiere_bifurcations_2013; @mackay_bifurcations_2014; @MacKay_morse_2015].

Thus, calculation of reaction rate (or flux) based on the geometry of
phase space structures requires identifying trajectories that start in
the *reactant* well, cross the dividing surface constructed from the
NHIM, and reach the *product* well. This dividing surface has been shown
to be the appropriate (locally no-recrossing) surface that reactive
trajectories must cross since the calculated reaction rates do not need
correction due to recrossings [@waalkens2004direct]. This construction
is in contrast to the "standard" transition-state-theory (TST) for
constructing the dividing surface which is only exact in gas phase
unimolecular reactions and when ergodicity of trajectories in the phase
space holds [@pollak_transition-state_1995]. As is now established, the
no recrossing (locally) property of a dividing surface is a contribution
of the phase space perspective of chemical
reactions [@uzer_geometry_2002]. While the standard TST relies on
recrossing free surface for calculating reaction flux, a dividing
surface constructed in the configuration space violates this condition
in the case of solvent, and the TST based reaction rate is not exact.
This violation of the recrossing property when DS is constructed in the
configuration space of a reaction in a high viscosity solvent also
follows from the Kramers' diffusion model, Langevin equation, of
chemical reactions [@pollak_transition-state_1995]. Thus, finding the
reactive trajectories, and the changes in the DS and NHIM due to a
solvent presents a worthwhile step towards development of the role of
phase space structures in reaction dynamics.

The geometry of unimolecular reactions dynamics has been developed using
a $N$ degrees of freedom Hamiltonian where the coordinates represent
intermolecular bonds. As a natural step in studying unimolecular
reaction dynamics in solvents, we adopt a model where the reaction
coordinates (modeled as a system Hamiltonian) are coupled with a set of
harmonic bath modes (modeled as a bath
Hamiltonian) [@berezhkovskii_activated_1992; @hershkovitz_multidimensional_1997; @reese_curvilinear_1998; @craig_chemical_2005].
This is with the intention of parametrizing the effects of a bath (or
solvent) on the system dynamics. This formulation of coupling harmonic
bath modes with system dynamics also serves as a preliminary step in
assessing the capabilities of a trajectory diagnostic called Lagrangian
descriptors (LDs) [@mancho_2013] in *realistic* (high dimensional)
chemical systems. This system-bath model will also serve as a test bed
for illustrating the use of LDs in directly computing the chemical
observables. But before doing that, we would like to identify the
reactive trajectories and develop a systematic approach for transition
path sampling (in the sense of rare event sampling) for reactions in
solvent. Thus, can we identify the qualitative changes in the
isomerization of a molecule in the presence of a solvent? In particular,
we would like to understand the influence of the coupling strength, the
solvent's viscosity, and the number of bath modes on the reactive
trajectories. We will answer this using LDs which have been shown
recently to be of use in detecting phase space structures that mediate
reactive trajectories in dissipative, time-dependent models of chemical
reactions [@junginger_uncovering_2016; @junginger_transition_2016], and
for transition path sampling in two degrees of freedom models of
chemical reactions [@patra_detecting_2018]. In this article, we will
present an approach for sampling reactive trajectories and identifying
reactive islands in high dimensional phase space of an isomerization
reaction (system) in a solvent (bath coupled to each of the system
configuration coordinates) which is following our work on quadratic
normal form Hamiltonian systems [@naik2019afinding].


### Two DoF de Leon-Berne system coupled with bath modes

In this section we describe the system-bath model for a system with two
degrees-of-freedom (DoF). While system-bath models have received a great
deal of attention in the chemistry and physics community, models with
systems having more than 1 DoF have received much less attention.

We will consider the 2 degree-of-freedom
Hamiltonian [\[eqn:ham\_db\]](#eqn:ham_db){reference-type="eqref"
reference="eqn:ham_db"} coupled with harmonic oscillators for a
system-bath model of the
form [@berezhkovskii_activated_1992; @hershkovitz_multidimensional_1997; @reese_curvilinear_1998]:

\begin{align}
\mathcal{H}(x,y,x_j,y_j,p_x,p_y,p_{x_j},p_{y_j})  =  \underbrace{\frac{p_x^2}{2 m_s} +  \frac{p_y^2}{2 m_s} + V_{DB}(x,y)}_{\text{System Hamiltonian}} + & \underbrace{\sum_{j=1}^{N_B} \frac{1}{2} \left[ \frac{p_{x_j}^2}{m_j} + \left(\omega_j x_j -\frac{c_{x,j} x}{\omega_j} \right)^2 \right]}_{\text{Coupling of the bath to $x$}} \\
& + \underbrace{\sum_{j=1}^{N_B} \frac{1}{2} \left[ \frac{p_{y_j}^2}{m_j} + \left(\omega_j y_j - \frac{c_{y,j} y}{\omega_j} \right)^2 \right]}_{\text{Coupling of the bath to $y$}} \\
%& = \mathcal{H}_{\text{system}} + \mathcal{H}_{\text{bath coupled to $x$}} + \mathcal{H}_{\text{bath coupled to $y$}}
\label{eqn:sb_ham}
\end{align}

where $x$ and $y$ denote the configuration space coordinates of the
system, $p_x$ and $p_y$ are the associated conjugate momenta, $p_{x_j}$,
$x_j$ denote the $j^{\rm th}$ bath phase space coordinates associated
with the system configuration space variable $x$, and $p_{y_j}$, $y_j$
denote the $j^{\rm th}$ bath phase space coordinates associated with the
system configuration space variable $y$. We assume that the frequencies,
$\omega_j$, are the same for each bath, and the coupling constants for
each configuration space variable to the bath are given by $c_{x,j}$ and
$c_{y,j}$. We explicitly carry out the
discretization that gives us the coupling constants $c_{x,j}$ and
$c_{y,j}$ and the frequencies $\omega_j$. 


#### Discretization of the Spectral Density: Derivation of the Parameters of the Bath {#app:discrete}

The coupling of the bath of harmonic oscillators to the configuration
space coordinates is described by a spectral density:

$$\begin{aligned}
J_x (\omega) = \frac{\pi}{2} \sum_{i=1}^{N_B} \frac{c_{x,i}^2}{\omega_i} \delta (\omega - \omega_i), \label{x-disc_spec} \\
J_y (\omega) = \frac{\pi}{2} \sum_{i=1}^{N_B}
\frac{c_{y,i}^2}{\omega_i} \delta (\omega - \omega_i),
\label{y-disc_spect}\end{aligned}$$

and these result from the discretization of a continuous Ohmic (linear)
form with an exponential cutoff:

$$\begin{aligned}
\bar{J}_x (\omega) = \eta_x \omega e^{-\frac{\omega}{\omega_{c}}}  , \label{z_1-cont_spec} \\
\bar{J}_y (\omega) = \eta_y \omega e^{-\frac{\omega}{\omega_{c}}} , \label{z_2-cont_spect}\end{aligned}$$

the discretization that gives us the coupling constants $c_{x,j}$,
$c_{y,j}$, and the frequencies $\omega_j$, and are given by
$$\omega_j = -\omega_c \log \left(\frac{j-\frac{1}{2}}{N_B} \right),
\quad j = 1, \ldots, N_B. \label{disc_freq_1}$$ and $$\begin{aligned}
c_{x,j} = \sqrt{\frac{2 \eta_{x} \omega_{c}}{\pi N_B}} \, \omega_j, \label{coupling_x} \\
c_{y,j} = \sqrt{\frac{2 \eta_{y} \omega_{c}}{\pi N_B}} \, \omega_j, \quad j = 1, \ldots, N_B. \label{coupling_y}\end{aligned}$$

In this appendix we describe a scheme for discretizing the continuous
spectral density which was given in [@craig_chemical_2005]. In the
following we will drop the subscripts $z_1$ and $z_2$ on the various
quantities for the sake of a simpler notation since we will follow the
same discretization procedure for each spectral density. The subscripts
can then be added back afterwards.

Re-establishing the notation, the discrete spectral density function is
given by:

$$J(\omega) = \frac{\pi}{2} \sum_{i=1}^{n_b} \frac{c_i^2}{\omega_i}
\delta (\omega - \omega_i), \label{sf_disc}$$

and the continuous spectral density is given by:

$$\bar{J} (\omega) = \eta \omega e^{-\frac{\omega}{\omega_c}}.
\label{sf_cont}$$

Discretization is obtained by carrying out the following steps.

1.  Require

    $$\int_{0}^{\infty} J(\omega) F(\omega) d \omega \approx
        \int_{0}^{\infty} \bar{J} (\omega) F(\omega) d \omega,
        \label{quad_1}$$

    for any integrable function $F(\omega)$

2.  Substitute ([\[sf\_disc\]](#sf_disc){reference-type="ref"
    reference="sf_disc"}) into the left-hand side of
    ([\[quad\_1\]](#quad_1){reference-type="ref" reference="quad_1"}) to
    obtain:

    $$\frac{\pi}{2} \sum_{i=1}^{n_b} \frac{c_i^2}{\omega_i} F(\omega_i).
        \label{quad_2}$$

3.  Approximate the right-hand side of
    ([\[quad\_1\]](#quad_1){reference-type="ref" reference="quad_1"}) by
    an appropriate quadrature. This is the step that we will now carry
    out in detail.

The quadrature recommended in [@craig_chemical_2005] is the midpoint
rule, after changing variables in the integral to $x =
e^{-\frac{\omega}{\omega_c}}$. The reason that they give for choosing
this quadrature is that it gives a uniform distribution of grid points
in the unit interval $0 < x_i < 1$ and therefore a logarithmic
distribution of bath frequencies $\omega_i = -\omega_c
\log x_i$. They claim that such a distribution of bath frequencies is
appropriate on physical grounds for an exponentially decaying density of
bath states.

The change of variables gives:

$$\begin{aligned}
x = e^{-\frac{\omega}{\omega_c}}  \Rightarrow  \log x =
-\frac{\omega}{\omega_c}  \Rightarrow \omega =-\omega_c \log x,
\label{cv_1} \\
dx = -\frac{1}{\omega_c} e^{-\frac{\omega}{\omega_c}} d \omega
=  -\frac{1}{\omega_c} x d \omega.\label{cv_2}\end{aligned}$$

Then ([\[sf\_cont\]](#sf_cont){reference-type="ref"
reference="sf_cont"}) becomes

$$\bar{J} (\omega) = -\eta x \omega_c \log x, \label{sf_cont_cv}$$

and using this expression, and the change of variables given in
([\[cv\_1\]](#cv_1){reference-type="ref" reference="cv_1"}) and
([\[cv\_2\]](#cv_2){reference-type="ref" reference="cv_2"}), the
right-hand side of ([\[quad\_1\]](#quad_1){reference-type="ref"
reference="quad_1"}) becomes:

$$-\int_{0}^{1} \eta \omega_c^2 ( \log x)  F(\omega (x)) dx.
\label{quad_3}$$

We discretize this integral using the midpoint rule. We partition the
unit interval into $n_b$ intervals of length $\frac{1}{n_b}$ and
evaluate the integrand at the midpoint, $x_i$, of each sub-interval:

$$x_i = \frac{i-\frac{1}{2}}{n_b}, \, i=1, \ldots , n_b,
\label{midpoint}$$

and obtain:

$$-\int_{0}^{1} \eta \omega_c^2 ( \log x)  F(\omega (x)) dx \approx
-\sum_{i=1}^{n_b} \eta \omega_c^2 \log x_i F(\omega(x_i))
\frac{i}{n_b}. \label{disc_quad_1}$$

Equating each term in the sum
([\[quad\_2\]](#quad_2){reference-type="ref" reference="quad_2"}) to
each term in the sum
([\[disc\_quad\_1\]](#disc_quad_1){reference-type="ref"
reference="disc_quad_1"}) gives:

$$\frac{\pi}{2} \frac{c_i^2}{\omega_i} = -\frac{1}{n_b} \eta
\omega_c^2 \log x_i = \frac{1}{n_b} \eta \omega_c \omega_i,$$

which gives:

$$c_i =\sqrt{\frac{2 \eta \omega_c}{\pi n_b}} \, \omega_i, \quad i=1,
\ldots, n_b. \label{coupling}$$

and from ([\[cv\_1\]](#cv_1){reference-type="ref" reference="cv_1"}) and
([\[midpoint\]](#midpoint){reference-type="ref" reference="midpoint"})
we see that:

$$\omega_i = -\omega_c \log \left(\frac{i-\frac{1}{2}}{n_b} \right),
\quad i=1, \ldots, n_b. \label{disc_freq_2}$$

[@craig_chemical_2005] choose $\omega_c = 500 {\rm cm}^{-1}$.

Thus, these quantities are given by

\begin{equation}
\omega_j = -\omega_c \log \left(\frac{j-\frac{1}{2}}{N_B} \right),
\quad j=1, \ldots, N_B \label{eqn:disc_freq_1}
\end{equation}

and

\begin{align}
c_{x,j} =  \sqrt{\frac{2 \eta_{x} \omega_{c}}{\pi N_B}} \, \omega_j, \;
c_{y,j}  = \sqrt{\frac{2 \eta_{y} \omega_{c}}{\pi N_B}} \,\omega_j, 
\quad j=1, \ldots, N_B. 
\label{eqn:coupling_coeff}
\end{align}

 We note here that the baths are coupled to
each other through the system dynamics. Hamilton's equations for the
system-bath dynamics, using
Eqn. [\[eqn:sys\_pot\]](#eqn:sys_pot){reference-type="eqref"
reference="eqn:sys_pot"} and
Eqn. [\[eqn:sb\_ham\]](#eqn:sb_ham){reference-type="eqref"
reference="eqn:sb_ham"}, are given by

\begin{align}
\dot{x} = & \frac{\partial H}{\partial p_x} =  \frac{p_x}{m_s}, \\
\dot{y} = & \frac{\partial H}{\partial p_y} = \frac{p_y}{m_s}, \\
\dot{x}_j = & \frac{\partial H}{\partial p_{x_j}} =  \frac{p_{x_j}}{m_j},  \qquad j=1, \ldots, N_B \\
\dot{y}_j = & \frac{\partial H}{\partial p_{y_j}} =  \frac{p_{y_j}}{m_j}, \qquad j=1, \ldots, N_B, \\
\dot{p}_x = &-\frac{\partial H}{\partial x} =  2 D_x \lambda \exp(-\lambda x) 
(\exp(-\lambda x) - 1) + \\ & \qquad \qquad \dfrac{\mathcal{V}^{\ddagger}}{y_w^4}\zeta 
\lambda y^2(y^2 - 2y_w^2)\exp(-\zeta \lambda x) + \sum_{j=1}^{N_B} \frac{c_{x,j}}{\omega_j} \left(\omega_j x_j -  \frac{c_{x,j} x}{\omega_j} \right), \\
\dot{p}_y = &- \frac{\partial H}{\partial y} = -4 \dfrac{\mathcal{V}^{\ddagger}}{y_w^4}y(y^2 - y_w^2)\exp(- \zeta \lambda x) + \sum_{j=1}^{N_B} \frac{c_{y,j}}{\omega_j} \left(\omega_j y_j -  \frac{c_{y,j} y}{\omega_j} \right), \\
\dot{p}_{x_j} = &- \frac{\partial H}{\partial {x_j}} = - \omega_j \left(\omega_j x_j - \frac{c_{x, j} x}{\omega_j} \right) ,  \qquad j=1, \ldots, N_B \\
\dot{p}_{y_j} = & -\frac{\partial H}{\partial {y_j}} = - \omega_j \left(\omega_j y_j - \frac{c_{y, j} y}{\omega_j} \right) , \qquad j=1, \ldots, N_B, \\
\label{eqn:Hameq_sys_bath}
\end{align}

where the bath frequencies $\omega_j$ and
coupling coefficients $c_{x,j},c_{y,j}$ are given by
Eqn. [\[eqn:disc\_freq\_1\]](#eqn:disc_freq_1){reference-type="eqref"
reference="eqn:disc_freq_1"} and
Eqn. [\[eqn:coupling\_coeff\]](#eqn:coupling_coeff){reference-type="eqref"
reference="eqn:coupling_coeff"}, respectively. The total energy will be
denoted by
$\mathcal{H}(x,y,x_j,y_j,p_x,p_y,p_{x_j},p_{y_j}) = E = E_{\rm saddle} + \Delta E$
where $j=1, 2, \ldots, N_B$ and $\Delta E$ is the excess energy with
respect to the isomerization barrier energy. In this article, we are
adopting the contraction $x_j$ to denote $x_1, x_2, \ldots, x_{N_B}$;
$p_{x_j}$ to denote $p_{x_1}, p_{x_2}, \ldots, p_{x_{N_B}}$, and so
forth.

In this form of coupling each bath mode is coupled with all the degrees
of freedom of the system, the total number of degrees of freedom are
$N_S N_B + N_S$ and as a result the dimension of phase space is
$2N_S(N_B + 1)$. Thus, it gives rise to a high dimensional Hamiltonian
model that allows us to incorporate a high number of degrees of freedom
by increasing the number of bath modes. In this study, we will first
focus on varying the number of bath modes, in particular when
$N_B = 1, 2, 4, 8, 16$ which corresponds to $8, 12, 20, 36, 68$
dimensional phase space. Even though, the number of bath modes are not
large enough to capture the effect of a solvent dynamics in the Langevin
sense, but nonetheless this formulation can be extended as long as
enough computational time is reserved to obtain trajectories and the
associated diagnostic method. However, the high dimensionality of the
system-bath model will be a setting to assess how the method of
Lagrangian descriptors performs in more realistic chemical systems. This
also gives a natural way to model the unimolecular conformational
isomerization in a solvent and allows us to vary the friction of the
environment from low-density gases to high-density
liquids [@tucker_reaction_2000].

### Model verification

\centering
\subfigure[]{\includegraphics[width=0.24\textwidth]{./figures/fig3B2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E1-500_x-px_tau100_NB1.png}}
\subfigure[]{\includegraphics[width=0.24\textwidth]{./figures/fig3B2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB1.png}}
\centering
\subfigure[]{\includegraphics[width=0.24\textwidth]{./figures/fig3C2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E1-500_x-px_tau100_NB1.png}}
\subfigure[]{\includegraphics[width=0.24\textwidth]{./figures/fig3C2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB1.png}}

## Revealing Phase Space Structures

The equilibrium points of the system-bath
model [\[eqn:Hameq\_sys\_bath\]](#eqn:Hameq_sys_bath){reference-type="eqref"
reference="eqn:Hameq_sys_bath"} are located at
$\overline{q}_s = (0,0,0_j,0_j,0,0,0_j,0_j)$ and at

$$\overline{q}_c = (x_{eq},\pm y_w,(c_{x,j}/\omega_j^2)x_{eq},(c_{y,j}/\omega_j^2)y_w,0,0,0_j,0_j)$$

where $x_{eq}$ is the $x$-coordinate of the equilibrium point in the
system model and needs to be obtained using numerical nonlinear root
solver. This shows the system coordinates $(x,y)$ of the equilibrium
points do not depend on the bath parameters and thus the isomer wells
are still at the same location. This lets us use the same surface of
section to compute Lagrangian descriptor even when increasing the bath
modes.


As noted in the article, the location of the equilibrium points in the
system-bath model are dependent on the frequencies and coupling strength
of the bath modes. Hence, we would like to analyze the stability of
these equilibrium points. The Jacobian of the system-bath Hamiltonian
vector
field [\[eqn:Hameq\_sys\_bath\]](#eqn:Hameq_sys_bath){reference-type="eqref"
reference="eqn:Hameq_sys_bath"} is

$$\begin{aligned}
\mathbb{J} =  \begin{pmatrix}
\text{\normalfont\Large\bfseries 0}& \hspace*{-\arraycolsep}\vline\hspace*{-\arraycolsep}& \mathbf{M} \\
\hline
\frac{\partial^2 V_{\rm SB}}{\partial  x_i x_j} & \hspace*{-\arraycolsep}\vline\hspace*{-\arraycolsep}& \text{\normalfont\Large\bfseries 0}
\end{pmatrix},
\quad \text{where} \; i, j = 1, 2, \ldots , N_S(N_B + 1)
\label{eqn:jacobian_sbham_vecfield}\end{aligned}$$

where each block matrix is of size $N_S(N_B + 1) \times N_S(N_B + 1)$.
The matrix $\mathbf{M}$ is given by:

$$\begin{aligned}
\mathbf{M} = \begin{pmatrix}
\frac{1}{m_s} & 0 & 0 & 0 & 0 & \cdots & \cdots & \cdots & \cdots & 0 \\
0 & \frac{1}{m_s} & 0 & 0 & 0 & \cdots & \cdots & \cdots & \cdots & 0 \\
0 & 0 & \frac{1}{m_1} & 0 & 0 & \cdots & \cdots & \cdots & \cdots & 0 \\
0 & 0 & 0 & \frac{1}{m_2} & 0 & \cdots & \cdots & \cdots & \cdots & 0 \\
0 & 0 & 0 & 0 & \ddots & \cdots & \cdots & \cdots & \cdots & 0 \\
0 & 0 & 0 & 0 & 0 & \cdots & \frac{1}{m_1} & 0 & \cdots & 0 \\
0 & 0 & 0 & 0 & 0 & \cdots & 0 & \frac{1}{m_2} & \cdots & 0 \\
0 & 0 & 0 & 0 & 0 & \cdots & \cdots & \cdots & \ddots & 0
\end{pmatrix}\end{aligned}$$

which is a diagonal matrix with mass of system and bath degrees of
freedom along the main diagonal.

The Hessian of the potential energy function for the system-bath model
is given by:

$$\begin{aligned}
\frac{\partial^2 V_{\rm SB}}{\partial  x_i x_j} = \begin{pmatrix}
- \frac{\partial^2 V_{\rm DB}}{\partial^2  x^2} - \sum\limits_{j} (\frac{c_{x,j}}{\omega_j})^2 & -\frac{\partial^2 V_{\rm DB}}{\partial y \partial x} & c_{x,1} & c_{x,2} & \cdots & c_{x,N_B} & 0 & 0 & \cdots & 0 \\
-\frac{\partial^2 V_{\rm DB}}{\partial y \partial x} & - \frac{\partial^2 V_{\rm DB}}{\partial^2  y^2} - \sum\limits_{j} (\frac{c_{y,j}}{\omega_j})^2 & 0 & 0 & \cdots & 0 & c_{y,1} & c_{y,2} & \cdots & c_{y,N_B} \\
c_{x,1} & 0 & -\omega_1^2 & 0 & \cdots & 0 & 0 & 0 & \cdots & 0 \\
c_{x,2} & 0 & 0 & -\omega_2^2 & \cdots & 0 & 0 & 0 & \cdots & 0 \\
\vdots & 0 & 0 & 0 & \ddots & 0 & 0 & 0 & \cdots & 0 \\
c_{x,N_B} & 0 & 0 & 0 & \cdots & -\omega_{N_B}^2 & 0 & 0 & \cdots & 0 \\
0 & c_{y,1} & 0 & 0 & \cdots & 0 & -\omega_1^2 & 0 & \cdots & 0 \\
0 & c_{y,2} & 0 & 0 & \cdots & 0 & 0 & -\omega_2^2 & \cdots & 0 \\
0 & \vdots & 0 & 0 & \cdots & 0 & 0 & 0 & \ddots & 0 \\
0 & c_{y,N_B} & 0 & 0 & \cdots & 0 & 0 & 0 & \cdots & -\omega_{N_B}^2 \\
\end{pmatrix}\end{aligned}$$

The Jacobian evaluated at the equilibrium points $\overline{x}$ is given
by

$$\begin{aligned}
\mathbb{J}(\overline{q}) =  \begin{pmatrix}
\text{\normalfont\Large\bfseries 0}& \hspace*{-\arraycolsep}\vline\hspace*{-\arraycolsep}& \mathbf{M} \\
\hline
\frac{\partial^2 V_{\rm SB}}{\partial  x_i x_j} (\overline{q}) & \hspace*{-\arraycolsep}\vline\hspace*{-\arraycolsep}& \text{\normalfont\Large\bfseries 0}
\end{pmatrix},
\quad \text{where} \; i, j = 1, 2, \ldots , N_S(N_B + 1)
\label{eqn:jacobian_sbham_eqpt}\end{aligned}$$

The eigenvalues are given by the solutions of the characteristic
equation:

$$\rm{det} \left( \mathbb{J}(\overline{q}) - \lambda \mathbb{I} \right) = 0$$

where $\mathbb{J}(\overline{q})$ and $\mathbb{I}$ are
$2N_S(N_B + 1) \times 2N_S(N_B + 1)$ matrices.

\centering
![Eigenvalues of the equilibrium point $\overline{q}_s$ of the
system-bath
model [\[eqn:Hameq\_sys\_bath\]](#eqn:Hameq_sys_bath){reference-type="eqref"
reference="eqn:Hameq_sys_bath"} where the real part is shown as red dots
and imaginary part is shown as blue cross. The pair of purely positive
and negative eigenvalues verify the equilibrium point is
index-1.](eigenvalues_systembath_NB64.pdf){width="0.95\linewidth"}

The linear stability analysis of these
equilibrium points gives their stability does not change from the system
model by adding the bath modes in the form of
Eqn. [\[eqn:sb\_ham\]](#eqn:sb_ham){reference-type="eqref"
reference="eqn:sb_ham"}. Next, the total energy of the index-1
equilibrium point at $\overline{q}_s$ is
$\mathcal{H}(0,0,0_j,0_j,0,0,0_j,0_j) = \epsilon_s$ which is the same as
in the system
Hamiltonian [\[eqn:ham\_db\]](#eqn:ham_db){reference-type="eqref"
reference="eqn:ham_db"}.

## Implications for Reaction Dynamics

### Detecting changes in reactive islands using Lagrangian descriptors

It has been demonstrated that the reactive islands and their
hierarchical structure can be used to search for new reactive
trajectories using the so-called shooting
method [@patra_detecting_2018]. In high dimensional molecular phase
space, this can be a non-trivial task due to the available phase space
volume and so has been phrased as "harvesting rare trajectories" and
computed using transition path sampling methods [@patra_classical_2015].
In this section, we would like to investigate if the LD contour maps can
provide a systematic approach to finding reactive trajectories and hence
building the high dimensional analogs of reactive islands in a
system-bath model.

**$N_S(N_B + 1)$-DoF: $N_B$ bath modes coupled to each of the $N_S$ system DoF**

To sample the reactive trajectories, we introduce the following
isoenergetic two-dimensional slice in the $2N_S(N_B + 1)$-dimensional
phase space of
Eqn. [\[eqn:Hameq\_sys\_bath\]](#eqn:Hameq_sys_bath){reference-type="eqref"
reference="eqn:Hameq_sys_bath"}. $$\begin{aligned}
U_{xp_x}^- &= \left\{(x,y,x_i,y_i,p_x,p_y,p_{x_i},p_{y_i}) \; | \; y = y_w, x_i = 0, y_i = 0, p_{x_i} = 0, p_{y_i} = 0, p_y(\cdot;E) < 0 \right\} 
\label{eqn:sos_Uxpx_systembath} \end{aligned}$$

### Influence of the number of bath modes

Let us consider the bath parameters $\eta_{x} = \eta_{y} = 0.01$ and
$\omega_{c} = 1.0$ and study the effect of increasing the number of bath
modes on the reactive islands of the 2 DoF isomerization.

Changes in reactive islands with increasing number of bath modes for
$\tau = 100$ which we considered to be enough integration time for the
isomerization given by the
system [\[eqn:vec\_field\_db\]](#eqn:vec_field_db){reference-type="eqref"
reference="eqn:vec_field_db"}.

\centering
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3B2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E1-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3B2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3B2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E7-500_x-px_tau100_NB64.png}}
\centering
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3C2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E1-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3C2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3C2_eta1e-1_omega2236e-3_tau100/slice1d_backward_lag_desc_systembath_E7-500_x-px_tau100_NB64.png}}

### Influence of friction in the bath modes

\centering
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3B2_eta1e-2_omega1e0/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3B2_eta1e-1_omega1e0/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3B2_eta1e0_omega1e0/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}
\centering
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3C2_eta1e-2_omega1e0_tau100/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3C2_eta1e-1_omega1e0_tau100/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}
\subfigure[]{\includegraphics[width=0.32\textwidth]{./figures/fig3C2_eta1e0_omega1e0_tau100/slice1d_backward_lag_desc_systembath_E4-500_x-px_tau100_NB64.png}}









# References
