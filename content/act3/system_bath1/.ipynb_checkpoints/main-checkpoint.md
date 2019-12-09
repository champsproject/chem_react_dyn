---
author:
- 'Rafael Garcia-Meseguer'
bibliography: myBib.bib
csl: elsevier-without-titles.csl
pagetitle: Chapter
---

# Lagrangian Descriptors for a System-Bath Model

<!-- #region -->
## Introduction

Locating the Dividing Surface (DS) in chemical reactions has been the focus of many rate constant calculations both for reactions in solution and in vacuum thanks to the Transition State Theory; a classical theory developed by Wigner, Eyring, Gwynne Evans and Polanyi;[@Wigner1932;@Eyring1935;@Polanyi1935;@Eyring1941] that calculates the rate of the reaction as the equilibrium flux of reactive trajectories through that DS.

In this chapter, we take a different point of view on the two supposed assumptions of [TST](../../prologe/TST.md). By treating the locally non-recrossing criterion for the dividing surface as an approximation, one is implying that there exists some other criterion of higher priority for locating the TS. However, in our view the location of the TS is defined solely by the satisfaction of the locally non-recrossing criterion. Hence, if one finds a putative TS whose associated dividing surface suffers recrossing (leading to a reduction in the TST transmission coefficient), then the TS is, strictly speaking, in the wrong place.
It has been argued that there can be circumstances in which no dividing surface exists that would satisfy the locally non-recrossing criterion. However, this argument is based on reactive and recrossing trajectories passing through a common point. That is possible only in configuration space, because every point in phase space is unique. The original definition of the dividing surface used by Wigner was in phase space, and that is the definition we adopt here.[@MullenPeters2014] That said, there may be no problem with using such an approximate TS, as done by Grote Hynes Theory, for the computation of the overall rate constant of a single-step reaction.[@GroteHynes1980;@Hynes2002;@Hynes1985;@Hynes1982] However, as discussed in [@GarciaMeseguerTST2019], not every effect can be explained with configuration-space models. One of these effects, which we have called the inertial barrier, arises when changes in shape of a reacting solute are resisted by the solvent because of a timescale mismatch between the solute dynamics (typically ~100 fs for transit from a PES saddle point to the next local minimum, in vacuo) and solvent dynamics (typically 1 - 10 ps, and sometimes much longer, for relocation of solvent molecules in the first shell around the reacting solute).[@Bunker1972;@Maroncelli1995;@Peters1995;@WuTenMing2000;@Larregaray2005;@Kubarych2010;@Lanzani2012;@Kowalewski2014;@Carpenter2015] Another potential problem arises, for example, when one is concerned with branching to several products after a rate-determining TS. The supposed equilibrium assumption of TST is not really an assumption either, because it can be shown to be a consequence of the locally non-recrossing requirement.[@SteinfeldHase1989] Of course, finding the exact DS is not an easy task for reactions of polyatomic molecules and could even be technically impossible for reactions in solutions.

However, the construction of the DS for 2 DoF systems using the Lyapunov family of unstable POs was presented in a series of papers by Pollak, Pechukas, and Child in the late 1970's to early 1980's.[@PechukasPollak1977;@PollakPechukas1978;@PechukasPollak1979;@PollakChild1980;@Pechukas1982] The resulting Periodic Orbit Dividing Surface (PODS) is a hypersurface in phase space arising from the unstable PO. It has been shown to have the required no-recrossing properties described above. Of particular interest was the recognition that as the total energy of trajectories increased, the location of the PODS would change, and that, in general, its projection onto the PES did not need to pass through the index one saddle point that is considered to be the location of the Transition Structure in many chemical models. In the present chapter, we emphasize that it is not only energy which can cause the PODS to move, but at least in some models, its location may also be mass dependent. This feature is of particular relevance to the solvent-derived inertial barrier mentioned above.

The objective of this chapter<font color='red'>...</font>

## Development of the Problem

<font color='red'>(Introducing the model)</font>

Our model potential (see @fig:ModelFig) consists of a one-dimensional double well oscillator that represents the reactive system, coupled to a one-dimensional harmonic oscillator representing the bath. The only coupling between the two oscillators is a Lennard-Jones-like repulsion potential term. Consequently, this model is appropriate only for nonpolar systems; our concern here is the consequence of non-bonded interactions between solvent and solute, not the more commonly studied polar interactions.

![Schematic representation (left) and definitions (right) of the model system used in our study. ](figures/SB_model.png){#fig:ModelFig width=100%}

The Hamiltonian that describes the system is as follows:

\begin{equation}
H(\mathbf{x})=H(\mathbf{r},\mathbf{p})=\frac{p_1^2}{2\mu_1}+\frac{p_2^2}{2\mu_2}+\sum_{j=1}^{5} c_j r_1^{j-1}+c_6 (c_7-r_2 )^2+\frac{c_8}{(r_2-r_1 )^{12}}
\label{modelEq}
\end{equation}

where the subscripts 1 and 2 refer to the reactive system and bath oscillators respectively; $r=(r_1,r_2 )$ is the position of the two oscillators and $p=(p_1,p_2 )$ represents the conjugate momenta. The reduced mass of each oscillator is represented by $\mu$, and c are coefficients whose values are listed in @fig:ModelFig. The potential energy can be divided between $V_1 = \sum_{j=1}^{5} c_j r_1^{j-1}$ as the potential of the reactive system, $V_2=c_6 (c_7-r_2 )^2$ as the potential of the bath and $V_{int}=c_8/(r_2-r_1 )^{-12}$ as the interaction between the two; hence $V=V_1+V_2+V_{int}$. The potential of the reactant, shown in @fig:PESFig, is chosen to have a minimum at $r_1=1.0$ and a second one at $r_1=2.0$ , with respective potential energies $V_1=0.0$ and $V_1=-10$. The maximum energy is at $r_1=1.33867$ and $V_1=2.0$. The full potential, shown in @fig:PESFig, has a saddle point at $r_1=1.36561$ and $r_2=2.161769$ at $V=3.47291$. The "reactant" minimum occurs at $r_1 = 0.98779$, $r_2 = 1.80661$, $V = 0.77040$. The "product" minimum occurs at $r_1 = 1.98517$, $r_2 = 2.75642$, $V = -6.66284$.

![(Left) Reactive system's potential energy profile. (Right) Contours of the full potential energy surface. The contours are depicted in the $-7 \leq V \leq 6$ interval. ](figures/PES.png){#fig:PESFig width=100%}

## Revealing Phase Space Structures


### Method 1 - Lagrangian Descriptors based on Action Integrals

Lagrangian descriptors are calculated on a chosen phase space grid of initial conditions $\mathbf{x_0}$ at time $t = t_0$, evolving the trajectories for a fixed forward and backward integration time $\tau$. The general expression of LDs is:

\begin{equation}
M\left(\mathbf{x_0},t_0,\tau\right) = \int_{t-\tau}^{t+\tau} \mathcal{F}\left(\mathbf{x}(t);\mathbf{x}_0\right) dt
\label{LDEq}
\end{equation}

where $\mathcal{F}\left(\mathbf{x}(t);\mathbf{x}_0\right)$ is a positive and bounded scalar representing a geometrical or physical property of a trajectory with initial conditions $\mathbf{x_0}$ and initial time $t_0$; that is integrated over the time interval $[t_0-\tau, t_0+\tau]$. Because we are interested in the DS, and we know that the action is a minimum on its vicinity we will be using an action-like value of the form:

\begin{equation}
\mathcal{F}\left(\mathbf{x}(t);\mathbf{x}_0\right) = \sum_{i=1}^{N}\left(p_i \frac{dx_i}{dt}\right)^{1/2}
\label{actionLikeEQ}
\end{equation}

where N is the number of DoF, $p_i$ and $x_i$ are respectively the momenta and position of the DoF $i$. It is interesting to note that, by finding the modulus of each term separately, we could examine their effect on the LD independently, although that issue is not explored in the present paper. In this case the LD will be:

\begin{equation}
M\left(\mathbf{x_0},t_0,\tau\right) = \int_{t-\tau}^{t+\tau} \sum_{i=1}^{N}\left(p_i \frac{dx_i}{dt}\right)^{1/2} = \int_{t-\tau}^{t+\tau} \sum_{i=1}^{N}\left(v_i^2 m_i\right)^{1/2} dt
\label{LDEq2}
\end{equation}

The application of LD to higher DoF is already being studied in depth.[@NaikSWigginsPysRevE2019; @naik2019finding] However, these studies revolve around systems with a known Hamiltonian, which is not the case for many chemical processes. Also, the control of every DoF to create the required initial conditions becomes unbearable as soon as you convert the system into a full atomistic model. Therefore, if we want to calculate the LD for complex chemical systems we need a different approach that is described in the following section.
As in our previous study the mass of the reactive system ($\mu_1$) is set to a value of 1 and the reduced mass of the bath ($\mu_2$) will be given values of 0.1, 1 or 10. But in this case we will be using the LD methodology to locate the DS and the $\tau$ value will be 1 for $\mu_2=0.1, 1$ and 2 for $\mu_2=10$.

### Method 2

The main idea of this method comes from a solution that has been applied in almost every studied chemical reaction. This is, the reduction of the system to a 1 or 2 DoF problem using a Reaction Coordinate (RC). However, as in those chemical studies, the equations of motion are not reduced to a single coordinate but instead the RC is used as a measure of the location of the system at each timestep.

The method was developed with the following assumptions about the problem to be studied:

- The system's Hamiltonian is unknown.
- The system is represented as an atomistic model with Cartesian coordinates.
- We can combine the Cartesian coordinates into one or two collective variables that accurately represent the process we want to model.

Also, the following assumptions are not required but will make the computations much easier:

- The system will have a saddle point from which we will have a rough estimate of where is it.
- The DS is relatively close to that saddle point.
- It is a closed system.

As mentioned previously the definition of the initial conditions of the system can be a major problem when dealing with atoms that have three Cartesian coordinates and velocities each. However, as we want to focus our study in the reactivity of the system, the trajectories in which we are interested have something in common. They all cross a surface orthogonal to the IRC[ref] at the saddle point in their journey. So starting at that point is much easier to assign different initial velocities as a way of exploring part of the phase space.

This part of phase space, is what we will call the reactive phase space which, and by exploring it we can discriminate those phase space structures that have no direct influence in the reactive process.

The general description of the methodology is the following:

1. We define 1 or 2 RC that will measure the evolution of the system.
2. We establish a set of initial conditions for our RC in the configurational space as our initial point.
3. From this point, we integrate several trajectories with different momenta of a time length longer than $\tau$, both forward and backward in time.
4. For each trajectory of length $T$ we can calculate the LD with a value of $\tau$ for at least a number $(T-2\tau)/dt$ of *initial conditions* that have enough trajectory length before and after them.

The amount of trajectories and their length will strongly depend on the problem and the time it takes to the trajectory to explore the reactive phase space. But, in principle will follow the rule of the longer the better.

In this method, the evolution of the trajectory defines the initial conditions and not the user. Thus, avoiding the problem of defining many initial conditions. Also, although here we talk about having a single point as initial condition there is no reason for not including a second or several points in configurational space as initial conditions. This will be most useful when the DS is far from our initial condition as we will see in the System bath model.

## Implications for Reaction Dynamics

### Method 1

The images shown in @fig:LDFirst where obtained by the above mentioned grid methodology. In those figures, for each point in the plot, the momenta of the bath DoF was equal to 0. The location of the saddle point (blue line) and the PO associated to it (red line) were included in the plots to better understanding of the information obtained. Of significant relevance are the plots obtained for the position and momenta of the reactive DoF (a plots in @fig:LDFirst). In those plots we can see how the invariant manifolds converge at the PO that encloses the DS forming a crossing point that indicates its location.

![LD plots for the different reduced masses of the bath for (a) the position and momenta of the reactive DoF and (b) the position of both DoF, each with a close up of on the saddle point area. The blue straigth line in the first line of plots indicates the position of the saddle point and the red line in the close ups represents the PO enclosing the DS. ](figures/LDFirst.png){#fig:LDFirst width="96%"}

It is important to note that this plots represent a section of the full phase space in which all the kinetic energy is in the reactive DoF. Thus, all we can expect to see from the PODS is the intersection with this section. Of course this is not a big problem, as it is very easy to change the value of the momenta of the bath at which the LD is calculated. <font color='red'>(...)</font>

### Method 2

We could suppose that the 2DoF from the described system are instead 2 collective coordinates obtained from a system of multiples DoF. We then proceed to calculate the LD as described in our second method, first to check if we could replicate the results from the previous method and then to analyze the extra information we are generating with it.

In @fig:LDSecondSec we can see how the plots are very similar to those obtained in @fig:LDFirst. The main difference appears as we increase the mass of the bath where this effect reduces the ergodicity of the system. Even in those cases, there are similarities between plots in the places visited by the trajectories.

![Equivalent plots from @fig:LDFirst obtained with the second methodology ](figures/LDSecondSection.png){#fig:LDSecondSec width="100%"}

But this is just a small amount of information that we have obtained from the calculations. The full set of values that this method produces can be projected into a single plot. Obviously these projections have so many dynamical structures depicted in them that it is very difficult to see anything. But with a bit of *cleaning* we can obtain a nice picture of the dynamical structure of the system. The plots shown in @fig:LDProject where obtained by projecting in configurational space all the values obtained from the calculations. However, for each value proyected in the same point only the minimum value was represented. This is because the value of the LD is expected to be a minimum in the vicinity of the PODS. These plots show a clear picture of the behaviour of the trajectories after crossing the DS but also they have a clear definition of the projection of the PODS in configurational space.

![(a) LD plots for the different reduced masses of the bath in configurational space and (b) their respective close up. The value of LD was obtained by selecting the minimum value between all the values that share the same point in configurational space. The red line in the close ups represents the PO enclosing the DS.  ](figures/LDSecondProjection.png){#fig:LDProject width="100%"}

Moreover, if take a slice of the plots in @fig:LDProject we can visualize the effect of the DS in the value of the LD (see @fig:LDSlice). This is quite obvious in @fig:LDProject but it will be very useful in multiples degrees of freedom where there will be much more noise.

![Slice at bath position of 2.2 of the LD values along the reactive coordinate. ](figures/LDSecondSlice.png){#fig:LDSlice width="100%"}

We have shown that the DS can be easily identified by the calculation of LD, not only with the grid methodology, but also with our version of the method. It is expected that for multidimensional problems the results will not be as clear and sharp as with this simple model. After all, we are including a significant amount of DoF that will introduce a lot of noise in the calculations. However, if we can identify the DS signal in the LD calculation, then we can begin to understand the dynamical properties that make the system behave as it does, like in the former case, the coupling between the bath and the reactive system.

## References
<!-- #endregion -->

```python

```
