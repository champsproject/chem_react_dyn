---
author:
- 'Rafael Garcia-Meseguer'
bibliography: myBib.bib
csl: elsevier-without-titles.csl
pagetitle: Rafa's Chapter
---
# Chapter

## Introduction

Finding the location of Dividing Surfaces (DSs) in chemical reactions has been the focus of many rate constant calculations both for reactions in solution and in vacuum thanks to the Transition State Theory; a classical theory developed by Wigner, Eyring, Gwynne Evans and Polanyi;[@Wigner1932;@Eyring1935;@Polanyi1935;@Eyring1941] that calculates the rate of the reaction as the equilibrium flux of reactive trajectories through that DS.

Although there has been much discussion of quantum mechanical analogs of TST[@Pechukas1982] we here treat only the original classical version, because the dividing surface that will form our principal focus is incompatible with quantum mechanics. It is commonly claimed that conventional TST makes two main assumptions.[@Mahan1974] The first, called the equilibrium assumption, requires that the reactant state and TS be in thermal equilibrium. The maintenance of energetic equilibrium means that the thermalization maintaining this equilibrium is (at least) as fast as the rate at which these states are depopulated.[@TruhlarVTST2017;@GarretTruhlarGTST1979] The equilibrium condition is usually satisfied for most gas-phase bimolecular reactions and for reactions in the liquid phase, because energy exchange between solutes and solvent is usually rapid enough to maintain the equilibrium.[@JBAnderson1973;@JBAnderson1995] However, there are cases where equilibrium is not maintained, even in solution.[@EssafiHarvey2018] In addition, for unimolecular reactions of intermediates with low barriers to product formation, it is commonly the case that most trajectories coming from the reactant state will have enough energy in a product-forming reaction coordinate to cross the second barrier as soon as they reach it.[@Carpenter1985;EzraWiggins2014] The second claimed assumption, specifies that any trajectory crossing the TS dividing surface from the reactant state is on a path towards the product state and will reach it without recrossing the dividing surface prior to the product being reached.[@TruhlarVTST1980] Because the rate is calculated as the flux through the TS, any non-reactive trajectory that crosses the TS dividing surface, or reactive trajectory that crosses it more than once will increase the flux through the dividing surface, thus leading to an overestimate of the rate constant. This means that TST gives us an upper limit on the true rate constant, and that if we found a dividing surface without any recrossing then TST would give the exact value of the rate constant (subject to certain caveats[@JCP2016Maug]). In conventional TST the transition state dividing surface is located at the saddle point, which is the maximum energy point in the minimum energy path from reactants to products. However, TST is most powerful in the form of Variational Transition State Theory (VTST),[@TruhlarVTST2017;@GarretTruhlarGTST1979;@TruhlarVTST1980;@TruhlarGarretVTST1984;@GarretTruhlarJCP1979;@Keck1960;@Horiuti1938;@Wigner1937] which is a generalization of TST that removes the restriction on the dividing surface to cross the saddle point. In VTST the dividing surface is variationally optimized to minimize the rate constant, usually by finding the maximum free energy along the reaction path. Although this surface is properly located in phase space, most of the VTST calculations assume that the TS can be found in configuration space.[@Carpenter1985] TST and its modern development VTST have been extensively reviewed,[@Truhlar2008;@TruhlarVTST2017;@GarretTruhlarGTST1979;@TruhlarVTST1980;@TruhlarGarretVTST1984;@GarretTruhlarJCP1979;@Keck1960;@Horiuti1938;@Wigner1937;@TruhlarGao2002;@AllisonTruhlar1998] but here we will focus on in its application to condensed-phase reactions, its applicability and limitations when working in configuration space, and the potential of including phase space into the methodology.

In this chapter, we take a different point of view on the two supposed assumptions of TST. By treating the non-recrossing criterion for the dividing surface as an assumption, one is implying that there exists some other criterion of higher priority for locating the TS. However, in our view the location of the TS is defined solely by the satisfaction of the non-recrossing criterion. Hence, if one finds a putative TS whose associated dividing surface suffers recrossing (leading to a reduction in the TST transmission coefficient), then the TS is, strictly speaking, in the wrong place.
It has been argued that there can be circumstances in which no dividing surface exists that would satisfy the non-recrossing criterion. However, this argument is based on reactive and recrossing trajectories passing through a common point. That is possible only in configuration space, because every point in phase space is unique. The original definition of the dividing surface used by Wigner was in phase space, and that is the definition we adopt here.[MullenPeters2014] That said, there may be no problem with using such an approximate TS, as done by Grote Hynes Theory, for the computation of the overall rate constant of a single-step reaction.[@GroteHynes1980;@Hynes2002;@Hynes1985;@Hynes1982] However, as discussed in [@GarciaMeseguerTST2019], not every effect can be explained with configuration-space models. One of these effects, which we have called the inertial barrier, arises when changes in shape of a reacting solute are resisted by the solvent because of a timescale mismatch between the solute dynamics (typically ~100 fs for transit from a PES saddle point to the next local minimum, in vacuo) and solvent dynamics (typically 1 - 10 ps, and sometimes much longer, for relocation of solvent molecules in the first shell around the reacting solute).[@Bunker1972;@Maroncelli1995;@Peters1995;@WuTenMing2000;@Larregaray2005;@Kubarych2010;@Lanzani2012;@Kowalewski2014;@Carpenter2015] Another potential problem arises, for example, when one is concerned with branching to several products after a rate-determining TS. The supposed equilibrium assumption of TST is not really an assumption either, because it can be shown to be a consequence of the non-recrossing requirement.[SteinfeldHase1989] Of course, finding the exact DS is not an easy task for reactions of polyatomic molecules and could even be technically impossible for reactions in solutions.

However, the construction of the DS for 2 DoF systems using the Lyapunov family of unstable POs was presented in a series of papers by Pollak, Pechukas, and Child in the late 1970's to early 1980's.[@PechukasPollak1977;@PollakPechukas1978;@PechukasPollak1979;@PollakChild1980;@Pechukas1982] The resulting Periodic Orbit Dividing Surface (PODS) is a hypersurface in phase space arising from the unstable PO. It has been shown to have the required no-recrossing properties described above. Of particular interest was the recognition that as the total energy of trajectories increased, the location of the PODS would change, and that, in general, its projection onto the PES did not need to pass through the index one saddle point that is considered to be the location of the Transition Structure in many chemical models. In the present chapter, we emphasize that it is not only energy which can cause the PODS to move, but at least in some models, its location may also be mass dependent. This feature is of particular relevance to the solvent-derived inertial barrier mentioned above.

The objective of this chapter<font color='red'>...</font>

## Equations for the system

<font color='red'>(Introducing the model)</font>

Our model potential (see @fig:ModelFig) consists of a one-dimensional double well oscillator that represents the reactive system, coupled to a one-dimensional harmonic oscillator representing the bath. The only coupling between the two oscillators is a Lennard-Jones-like repulsion potential term. Consequently, this model is appropriate only for nonpolar systems; our concern here is the consequence of non-bonded interactions between solvent and solute, not the more commonly studied polar interactions.

![Schematic representation (left) and definitions (right) of the model system used in our study.](Figures\SB_model.png){#fig:ModelFig width=100%}

The Hamiltonian that describes the system is as follows:

\begin{equation}
H(\mathbf{x})=H(\mathbf{r},\mathbf{p})=\frac{p_1^2}{2\mu_1}+\frac{p_2^2}{2\mu_2}+\sum_{j=1}^{5} c_j r_1^{j-1}+c_6 (c_7-r_2 )^2+\frac{c_8}{(r_2-r_1 )^{12}}
\label{modelEq}
\end{equation}

where the subscripts 1 and 2 refer to the reactive system and bath oscillators respectively; $r=(r_1,r_2 )$ is the position of the two oscillators and $p=(p_1,p_2 )$ represents the conjugate momenta. The reduced mass of each oscillator is represented by $\mu$, and c are coefficients whose values are listed in @fig:ModelFig. The potential energy can be divided between $V_1 = \sum_{j=1}^{5} c_j r_1^{j-1}$ as the potential of the reactive system, $V_2=c_6 (c_7-r_2 )^2$ as the potential of the bath and $V_{int}=c_8/(r_2-r_1 )^{-12}$ as the interaction between the two; hence $V=V_1+V_2+V_{int}$. The potential of the reactant, shown in @fig:PESFig, is chosen to have a minimum at $r_1=1.0$ and a second one at $r_1=2.0$ , with respective potential energies $V_1=0.0$ and $V_1=-10$. The maximum energy is at $r_1=1.33867$ and $V_1=2.0$. The full potential, shown in @fig:PESFig, has a saddle point at $r_1=1.36561$ and $r_2=2.161769$ at $V=3.47291$. The "reactant" minimum occurs at $r_1 = 0.98779$, $r_2 = 1.80661$, $V = 0.77040$. The "product" minimum occurs at $r_1 = 1.98517$, $r_2 = 2.75642$, $V = -6.66284$.

![(Left)Reactive system's potential energy profile. (Right) Contours of the full potential energy surface. The contours are depicted in the $-7 \leq V \leq 6$ interval.](Figures\PES.png){#fig:PESFig width=100%}

## Phase Space Structures

One of the reasons for choosing this system is that periodic orbits that define the PODSs can be computed easily. All the calculations were carried out at an energy of $3.691966889$ which is slightly above the energy of the saddle point.

In order to understand the properties of the trajectories that depart from the DS we need to sample its points in phase space. The procedure, applicable to a 2 DoF Hamiltonian system, selects points on a 2D surface with fixed total energy (E), where the periodic orbit forms the one dimensional boundary of the DS. The algorithm is as described in [@JCP2016Maug; @ezra2018sampling]:

1. Locate an unstable PO.
2. Project the unstable PO into configuration space, which gives a curve in configuration space.
3. Choose points on the curve $(x_i,y_i)$ for $i=1,...,N$, where N is the desired number of points. The points are spaced uniformly according to distance along the PO.
4. For each point $(x_i,y_i )$ determine $p_{x_{max},i}$ by solving for $p_x$.

    \begin{equation}
    H(x_i,y_i,p_x,0)=\frac{p_x^2}{2\mu_x}+V(x_i,y_i)=E
    \label{PODSsampEq}
    \end{equation}

5. Note that the solution of this equation requires $E- V(x_i,y_i ) \ge 0$, and there will be two solutions, $\pm p_{x_{max},i}$.
6. For each point $(x_i,y_i)$ choose points $p_{x_j}$ for $j=1,...,K$, with $p_{x_1}=-p_{x_{max},i}$ and $p_{x_K}=p_{x_{max},i}$ and solve the equation $H(x_i,y_i,p_x,p_y)=E$ to obtain $p_y$.

The geometrical structure of the DS sampled in this manner is a one parameter family of circles. The parameter defining the family is given by the distance along the projection of the PO onto the configuration space from Steps 1-3 in the algorithm above, and the momentum-space circles are given by the following equation obtained from the Hamiltonian:

\begin{equation}
\frac{p_x^2}{2\mu_x}+\frac{p_y^2}{2\mu_y}=E-V(x_i,y_i)
\label{DSGeomEq}
\end{equation}

### Effect of the Solvent Mass

In @fig:DSCloseFig (top) we can see the projection of the calculated PODSs in configuration space. @fig:DSCloseFig also includes two approximations to the DS explained in the caption and shows how the three of them respond as $\mu_2$ changes. It can be seen that, for low reduced masses, the approximate DSs are close to the PODS. That is because the bath can rapidly adapt to the position of the reactive system. However, as $\mu_2$ increases the PODS starts to curve and to displace from the approximate DSs, moving closer to the product well.

The geometry of this one parameter family of circles depends on the nature of the projection of the PO into configuration space. In this particular case the PO projections are arcs where a configuration space point on the projection of the PO moves back and forth along the arc. This means that the endpoints of the arc are turning points with $p_x=p_y=0$, where the circles defined by eq. \ref{DSGeomEq} shrink to points. This implies that the geometry of the one parameter family of circles defines a 2D sphere (see @fig:DSCloseFig bottom).

![(Top)Close-up of the PES of the full system, near the saddle point region at different reduced masses of the bath. Each of the axis scales were weighted by the square root of its coordinate mass. The dashed red line is the intrinsic reaction coordinate (IRC). The blue line is DS if one assumes that the reaction coordinate is $r_1$. The red line is the DS projection at the saddle point, which is locally orthogonal to the IRC (It does not look orthogonal because of the choice of axis scales). The green line is the projection of the PO that defines the Dividing surface. (Bottom) Schematic representation of the DS's geometrical structure for the different reduced masses. The yellow structure represents the possible momenta depending of the location in the DS.](Figures\PODS_DSshape.png){#fig:DSCloseFig width=100%}

What the model reveals in @fig:DSCloseFig is that when $\mu_2$ is small with respect to $\mu_1$, the projections of the three dividing surfaces are quite close to each other, but as the relative reduced mass of the solvent model gets larger, the true dividing surface moves away from the ones that are rooted at the PES saddle point. The important point here is that, even if one could take the solvent into account properly in defining an IRC for a solution-phase reaction (which generally one cannot), a DS orthogonal to that IRC but still centered on the PES saddle point (the red line) would be incorrect.

But how can we be sure that these PODS are in fact the no-recrossing DS? From the sampled trajectories we can measure the time taken to reach a determined region (transit time), in this case the PES minimum identified as the product well. Then we can perform the same calculation but with trajectories starting on the DS defined only with $r_1$ (the blue line in @fig:DSCloseFig). The blue line corresponds to the common choice for solution phase reactions of assigning the transition state location to the PES saddle point, and assuming that the reaction coordinate is entirely determined by the solute. @fig:TransitFig is a representation in phase space of the transit times of trajectories that start on the true and approximate dividing surfaces with different initial $p_\perp$, the momentum normal to the dividing surface. The transit times (calculated as the time for $r_1$ to reach a value greater than that for product minimum) show brighter colors in Figure @fig:TransitFig as the transit time increases. The expected results for a DS is that trajectories starting with negative momenta normal to the dividing surface ($p_\perp<0$), i.e. directed to the reactant well, would take longer to reach the product well than those that start with positive momenta. This is clearly the case for the PODS as can be seen in Figure 6, where $p_\perp=0$ (which corresponds to the PO) provides an exact line of demarcation in the transit times. By contrast, the approximate DS shows long and short transit times on both sides of $p_1=0$. It is interesting to note that those areas where the transit times are long for $p_1>0$ or short for $p_1<0$ correspond to recrossing of trajectories, and that the amount of recrossing gets larger as $\mu_2$ increases. Thus, the approximate DS becomes a poorer and poorer choice for the transition state as the mass of the bath oscillator increases.

![A comparison of trajectory transit times from the DS to the product. (Top) Being the PODS (green DS in @fig:DSCloseFig) and (Bottom) the DS conventional definition of the TS (blue DS in @fig:DSCloseFig). The color scale goes from dark colors for short times to brighter colors for long times. The quantity $p_\perp$ is the momentum perpendicular to the dividing surface, with a positive sign being in the direction of the product.](Figures\TransitT.png){#fig:TransitFig width=100%}

The brighter colored bands visible on the reactant sides ($p_\perp <0$) in Figure @fig:TransitFig are associated with the many periodic orbits located in the reactant well. Trajectories that approach these POs can spend a long time before finally crossing over to the product well.

### Lagrangian Descriptors

#### Method 1

Lagrangian descriptors (LDs) is a scalar field diagnostic that can be used to explore the phase space structures at a given initial time. The general expression for LDs is:

\begin{equation}
M\left(\mathbf{x}_0,t_0\right)_\tau = \int_{t-\tau}^{t+\tau} \mathcal{F}\left(\mathbf{x}(t);\mathbf{x}_0\right) dt
\label{LDEq}
\end{equation} 

where initial conditions $\mathbf{x}_0 = \mathbf{x}(t_0)$ are integrated forward and backward for a fixed integration time $\tau$, and $\mathcal{F}\left(\cdot\right)$ is a positively bounded scalar representing a geometrical or physical property of a trajectory with initial conditions $\mathbf{x}_0$ and initial time $t_0$; that is integrated over the time interval $[t_0-\tau, t_0+\tau]$. Because we are interested in the DS, and we know that the action is a minimum on its vicinity we will be using an action-like value of the form:

\begin{equation}
\mathcal{F}\left(\mathbf{x}(t);\mathbf{x}_0\right) = \sum_{i=1}^{N}\left(p_i \frac{dx_i}{dt}\right)^{1/2}
\label{actionLikeEQ}
\end{equation}

where N is the number of DoF, $p_i$ and $x_i$ are respectively the momenta and position of the DoF $i$. It is interesting to note that, by finding the modulus of each term separately, we could examine their effect on the LD independently, although that issue is not explored in this chapter. In this case the LD will be:

\begin{equation}
M\left(\mathbf{x_0},t_0,\tau\right) = \int_{t_0-\tau}^{t_0+\tau} \sum_{i=1}^{N}\left(p_i \frac{dx_i}{dt}\right)^{1/2} = \int_{t_0-\tau}^{t_0+\tau} \sum_{i=1}^{N}\left(v_i^2 m_i\right)^{1/2} dt
\label{LDEq2}
\end{equation}

The application of LD to higher DoF is already being studied in depth.[@NaikSWigginsPysRevE2019; @naik2019finding] However, these studies revolve around systems with a known Hamiltonian, which is not the case for many chemical processes. Also, the control of every DoF to create the required initial conditions becomes unbearable as soon as you convert the system into a full atomistic model. Therefore, if we want to calculate the LD for complex chemical systems we need a different approach that is described in the following section.
As in our previous study the mass of the reactive system ($\mu_1$) is set to a value of 1 and the reduced mass of the bath ($\mu_2$) will be given values of 0.1, 1 or 10. But in this case we will be using the LD methodology to locate the DS and the $\tau$ value will be 1 for $\mu_2=0.1, 1$ and 2 for $\mu_2=10$.

The images shown in @fig:LDFirst where obtained by the above mentioned grid methodology. In those figures, for each point in the plot, the momenta of the bath DoF was equal to 0. The location of the saddle point (blue line) and the PO associated to it (red line) were included in the plots to better understanding of the information obtained. Of significant relevance are the plots obtained for the position and momenta of the reactive DoF (a plots in @fig:LDFirst). In those plots we can see how the invariant manifolds converge at the PO that encloses the DS forming a crossing point that indicates its location.

![LD plots for the different reduced masses of the bath for (a) the position and momenta of the reactive DoF and (b) the position of both DoF, each with a close up of on the saddle point area. The blue straigth line in the first line of plots indicates the position of the saddle point and the red line in the close ups represents the PO enclosing the DS.](Figures\LDFirst.png){#fig:LDFirst width="96%"}

It is important to note that this plots represent a section of the full phase space in which all the kinetic energy is in the reactive DoF. Thus, all we can expect to see from the PODS is the intersection with this section. Of course this is not a big problem, as it is very easy to change the value of the momenta of the bath at which the LD is calculated. <font color='red'>(...)</font>

#### Method 2

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

We could suppose that the 2DoF from the described system are instead 2 collective coordinates obtained from a system of multiples DoF. We then proceed to calculate the LD as described in our second method, first to check if we could replicate the results from the previous method and then to analyze the extra information we are generating with it.

In @fig:LDSecondSec we can see how the plots are very similar to those obtained in @fig:LDFirst. The main difference appears as we increase the mass of the bath where this effect reduces the ergodicity of the system. Even in those cases, there are similarities between plots in the places visited by the trajectories.

![Equivalent plots from @fig:LDFirst obtained with the second methodology](Figures\LDSecondSection.png){#fig:LDSecondSec width="100%"}

But this is just a small amount of information that we have obtained from the calculations. The full set of values that this method produces can be projected into a single plot. Obviously these projections have so many dynamical structures depicted in them that it is very difficult to see anything. But with a bit of *cleaning* we can obtain a nice picture of the dynamical structure of the system. The plots shown in @fig:LDProject where obtained by projecting in configurational space all the values obtained from the calculations. However, for each value proyected in the same point only the minimum value was represented. This is because the value of the LD is expected to be a minimum in the vicinity of the PODS. These plots show a clear picture of the behaviour of the trajectories after crossing the DS but also they have a clear definition of the projection of the PODS in configurational space.

![(a) LD plots for the different reduced masses of the bath in configurational space and (b) their respective close up. The value of LD was obtained by selecting the minimum value between all the values that share the same point in configurational space. The red line in the close ups represents the PO enclosing the DS. ](Figures\LDSecondProjection.png){#fig:LDProject width="100%"}

Moreover, if take a slice of the plots in @fig:LDProject we can visualize the effect of the DS in the value of the LD (see @fig:LDSlice). This is quite obvious in @fig:LDProject but it will be very useful in multiples degrees of freedom where there will be much more noise.

![Slice at bath position of 2.2 of the LD values along the reactive coordinate.](Figures\LDSecondSlice.png){#fig:LDSlice width="100%"}

We have shown that the DS can be easily identified by the calculation of LD, not only with the grid methodology, but also with our version of the method. It is expected that for multidimensional problems the results will not be as clear and sharp as with this simple model. After all, we are including a significant amount of DoF that will introduce a lot of noise in the calculations. However, if we can identify the DS signal in the LD calculation, then we can begin to understand the dynamical properties that make the system behave as it does, like in the former case, the coupling between the bath and the reactive system.

## References

```python

```
