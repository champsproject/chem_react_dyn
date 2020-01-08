---
author: " "
csl: ./elsevier-without-titles.csl
geometry: margin = 1in
output:
  html_document:
    fig_caption: yes
    fig_height: 1
  pdf_document:
    fig_caption: yes
    fig_height: 1
bibliography: ./finding_nhim_2-3dof_nonlinear.bib
---

# Introduction

It is well-known now that the paradigm of escape from a potential well
and the topology of phase space structures that mediate such escape are
used in a broad array of problems such as isomerization of molecular
clusters {% cite Komatsuzaki2001 --file finding_nhim_2and3dof_nonlinear %}, reaction rates in chemical
physics {% cite Komatsuzaki1999 WiWiJaUz2001 --file finding_nhim_2and3dof_nonlinear %}, ionization of a hydrogen atom
under electromagnetic field in atomic physics {% cite JaFaUz2000 --file finding_nhim_2and3dof_nonlinear %}, transport
of defects in solid state and semiconductor physics {% cite Eckhardt1995 --file finding_nhim_2and3dof_nonlinear %},
buckling modes in structural mechanics {% cite Collins2012 ZhViRo2018 --file finding_nhim_2and3dof_nonlinear %}, ship
motion and capsize {% cite Virgin1989 ThDe1996 NaRo2017 --file finding_nhim_2and3dof_nonlinear %}, escape and
recapture of comets and asteroids in celestial
mechanics {% cite JaRoLoMaFaUz2002 DeJuLoMaPaPrRoTh2005 Ross2003 --file finding_nhim_2and3dof_nonlinear %}, and
escape into inflation or re-collapse to singularity in
cosmology {% cite DeOliveira2002 --file finding_nhim_2and3dof_nonlinear %}. As such a method that can identify the high
dimensional phase space structures using low dimensional surface as
probes can aid in quantifying the escape rates. These low dimensional
surfaces has been shown to be of as *reactive islands* in chemical
physics and lead to insights into sampling rare transition
events {% cite patra_classical-quantum_2015 patra_detecting_2018 --file finding_nhim_2and3dof_nonlinear %}. However,
to benchmark the methodology, we first applied it to linear systems
where the closed-form analytical expression of the phase space
structures is known {% cite naik2019finding --file finding_nhim_2and3dof_nonlinear %}. As the next step, in this
chapter, we will focus on nonlinear Hamiltonian systems which have been
extensively studied as "built by hand" models of galactic dynamics and
for demonstrating quantum dynamical
tunneling {% cite barbanis_isolating_1966 brumer_variational_1976 davis_semiclassical_1979 heller_molecular_1980 waite_mode_1981 kosloff_dynamical_1981 contopoulos_simple_1985 founargiotakis_periodic_1989 barbanis_escape_1990 babyuk_hydrodynamic_2003 --file finding_nhim_2and3dof_nonlinear %}.
The nonlinear Hamiltonian systems considered here have an underlying
Hénon-Heiles type potential with the simplest form of nonlinearity, and
show regular, quasi-periodic, and chaotic trajectories along with
bifurcations of periodic orbits. A Hénon-Heiles type potential has a
well with bottlenecks connecting the region of bounded motion (trapped
region) to unbounded motion (escape off to infinity), and have
rotational symmetry. In addition, these Hénon-Heiles type potentials are
studied as first benchmark nonlinear systems in applying new phase space
transport methods to astrophysical and molecular motion. In this
chapter, we will present verification of a method that uses trajectory
diagnostic on a low dimensional surface for revealing the phase space
structures in 4 or more dimensions.

Conservative dynamics on an open potential well has received
considerable attention because the phase space structures, normally
hyperbolic invariant manifolds (NHIM) and its invariant manifolds,
explain the intricate fractal structure of ionization
rates {% cite mitchell_geometry_2003_I mitchell_geometry_2003_II mitchell_chaos-induced_2004 --file finding_nhim_2and3dof_nonlinear %}.
Furthermore, the discrepancies in observed and predicted ionization
rates in atomic systems has also been explained by accounting for the
topology of the phase space structures. These have been connected with
the breakdown of ergodic assumption that is the basis for using
ionization and dissociation rate
formulae {% cite de_leon_intramolecular_1981 --file finding_nhim_2and3dof_nonlinear %}. This rich literature on chaotic
escape of electrons from atoms sets a precedent for applying new methods
for finding NHIM and its invariant manifolds in Hamiltonian with open
potential wells
 {% cite mitchell_analysis_2004 mitchell_chaos-induced_2004 mitchell_nonlinear_2009 mitchell_structure_2007 wang_photoionization_2010 --file finding_nhim_2and3dof_nonlinear %}.

As we noted earlier, trajectory diagnostic methods which can probe phase
space to detect the high dimensional invariant manifolds have potential
to be of use in many degrees-of-freedom models. One such method is the
Lagrangian descriptors (LDs) that can reveal phase space structures by
encoding geometric property of trajectories (such as, phase space arc
length, configuration space distance or displacement, cumulative action
or kinetic energy) initialised on a two dimensional
surface {% cite madrid2detect009 mendoza2010 mancho2013 lopesino2017 --file finding_nhim_2and3dof_nonlinear %}. The
method was originally developed in the context of Lagrangian transport
in time-dependent two dimensional fluid mechanics. However, it has also
been successful in locating transition state trajectories in chemical
reactions {% cite balibrea2016lagrangian craven2017lagrangian junginger2016lagrangian --file finding_nhim_2and3dof_nonlinear %}.
Besides, also being applicable to both Hamiltonian and non-Hamiltonian
systems, as well as to systems with arbitrary time-dependence such as
stochastic and dissipative forces, and geophysical data from satellite
and numerical
simulations {% cite amism11 mendoza2014 ggmwm15 lopesino2017 ramos2018 --file finding_nhim_2and3dof_nonlinear %}.

The method of Lagrangian descriptors (LDs) is straightforward to implement
computationally and it provides a "high resolution" method for exploring
the influence of high dimensional phase space structure on trajectory
behaviour. The method of LDs takes an *opposite* approach to that of
classical Lyapunov exponent type calculations by emphasizing the initial
conditions of trajectories, rather than their advected locations that is
involved in calculating normalized rate of divergence. This is achieved
by considering a two dimensional section of the full phase space and
discretizing with a dense grid of initial conditions. Even though the
trajectories wander off in the phase space, as the initial conditions
evolve in time, there is no loss in resolution of the two dimensional
section. In contrast to inferring the phase space structures from
Poincaré sections, LD plots do not suffer from loss of resolution since
the affects of the structure are encoded in the initial conditions and
there is no need for the trajectory to return to the section. Our
objective is to clarify the use of Lagrangian descriptors as a
diagnostic on two dimensional sections of high dimensional phase space
structures. This diagnostic is also meant to be used as the preliminary
step in computing the NHIM, their stable and unstable manifolds using
other computational
means {% cite junginger2016transition bardakcioglu2018 ezra_2018 --file finding_nhim_2and3dof_nonlinear %}. In this
chapter, we will present the method's capability to detect the high
dimensional phase space structures such as the NHIM, their stable, and
unstable manifolds in the 2 DoF Barbanis system.

<!-- #region -->
# Barbanis 2 DoF Model
{#sec:model_prob_2dof}

## Model system: coupled harmonic 2 DoF Hamiltonian


As pointed out in the Introduction, our focus is to adopt a
well-understood model system which is a 2 degrees-of-freedom coupled
harmonic oscillator with the Hamiltonian

$$\begin{aligned}
\mathcal{H}(x,y,p_x,p_y) =& T(p_x, p_y) + V_{\rm B}(x,y) \\ 
=& \frac{1}{2}p_x^2 + \frac{1}{2}p_y^2 + \frac{1}{2}\omega_x^2 x^2 + \frac{1}{2}\omega_y^2 y^2 +
\delta x y^2  
\end{aligned}
\label{eqn:Hamiltonian_Barbanis}$$ 


where $\omega_x, \omega_y, \delta$ are the harmonic oscillator frequencies of the $x$ and $y$
degree-of-freedom, and the coupling strength, respectively. We will fix
the parameters as $\omega_x 
= 1.0, \omega_y = 1.1, \delta = -0.11$ in this study. The two
degrees-of-freedom potential is also referred to as *Barbanis*
potential, and has been investigated as a model of galactic
motion ({% cite contopoulos1970 barbanis_isolating_1966 --file finding_nhim_2and3dof_nonlinear %}), dynamical
tunneling and molecular spectra in physical
chemistry ({% cite heller1980 davis1981 martens1987 --file finding_nhim_2and3dof_nonlinear %}), structural
mechanics and ship capsize ({% cite ThDe1996 NaRo2017 --file finding_nhim_2and3dof_nonlinear %}).

The equilibria of the Hamiltonian vector field are located at

$$\left(-\frac{\omega_y^2}{2\delta}, \pm 
\frac{1}{\sqrt{2}}\frac{\omega_x \omega_y}{\delta}, 0, 0 \right) \qquad \text{and}  \qquad \left(0, 
0, 0, 0 \right)$$

and are at total energy $E_c = \frac{\omega_x^2 \omega_y^4}{8 \delta^2}$ and $0$ respectively. The energy of the two index-1 saddles (as defined and shown in
App. [5.2.1](#ssect:linear)) located at positive and negative
y-coordinates and positive x-coordinate for $\delta < 0$ will be
referred to as *critical energy*, $E_c$. In our discussion, we will
refer to the total energy of a trajectory or initial condition in terms
of the excess energy, $\Delta E = E_c - e$, which can be negative or
positive to denote energy below or above the critical energy. For the
parameters used in this study, the index-1 saddle equilibrium points are
located at $\left( 5.5, \pm 7.071, 0, 
0 \right)$ and have energy, $E_c = 15.125$.

The contours of the coupled harmonic 2 DoF potential energy function in \eqref{eqn:Hamiltonian_Barbanis} is shown in
