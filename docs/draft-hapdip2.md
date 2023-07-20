---
fontsize: 11pt
geometry: margin=1.2in
title: "
  Polygenic barriers to gene flow: the role of dominance, haploid selection and
  heterogeneous genetic architectures
  "
author: 
    - name: Arthur Zwaenepoel$^{1,\ast}$, Himani Sachdeva$^2$, Christelle Fraïsse$^1$
institutes:
    - "1. University of Lille, CNRS, UMR 8198 -- Evo-Eco-Paleo, F-59000 Lille, France"
    - "2. Department of Mathematics, University of Vienna, Vienna, Austria"
header-includes: 
  - \DeclareMathSymbol{\shortminus}{\mathbin}{AMSa}{"39}
  - \newcommand{\Ex}{\mathbb{E}}
  - \newcommand{\Var}{\mathrm{Var}}
  - \newcommand{\HW}{\mathrm{HW}}
  - \newcommand{\Bin}{\mathrm{Bin}}
  - \newcommand{\Beta}{\mathrm{Beta}}
  - \newcommand{\Exp}{\mathrm{Exponential}}
  - \newcommand{\Gam}{\mathrm{Gamma}}
  - \newcommand{\Bfun}{\mathrm{B}}
  - \newcommand{\eps}{\epsilon}
  - \newcommand{\all}{A}
  - \newcommand{\erf}{\mathrm{erf}}
  - \newcommand{\erfi}{\mathrm{erfi}}
  - \newcommand{\fix}{\mathrm{fix}}
  - \newcommand{\logit}{\mathop{\mathrm{logit}}}
  - \usepackage{caption}
  - \captionsetup[figure]{labelfont=bf,font=small}
  - \usepackage{float,soul}
  - \usepackage[normalem]{ulem} 
  - \makeatletter
  - \def\fps@figure{tb} 
  - \makeatother
  - \usepackage{algorithm}
  - \usepackage{algpseudocode}
  - \usepackage{lipsum}
abstract: \lipsum[1]
---

# Introduction {#sec:intro}

When a population is subdivided across multiple habitats with different 
environmental conditions, the extent to which distinct subpopulations can
maintain locally beneficial genetic variation depends on the rate of migration
between them.
Migration between populations that maintain divergently selected alleles can
lead to maladaptive gene flow, yielding a migration load (a reduction in mean
fitness due to the influx of locally maladaptive genes) or may lead to loss
of local adaptation altogether (so-called *swamping* by gene flow) [e.g.
@lenormand2002].
While local adaptation may be driven by a few conspicuous loci (e.g. adaptive
melanism in peppermoth [@hof2016] or pocket mice [@nachman2003]), it is believed
to typically be polygenic, with alleles of different effect at many loci across
the genome responding to selection [@pritchard2010; @lecorre2012; @barghi2020;
@bomblies2022; @stankowski2022; @westram2018].
When local adaptation is polygenic, migration from a population adapted to
different environmental conditions will produce linkage disequilibria (LD)
among selected loci, and the rate at which each individual invading locally
deleterious alleles is eliminated will be affected by such associations
[@barton1983; @feder2012; @yeaman2015; @sachdeva2022].
This in turn will affect the equilibrium migration load and swamping thresholds
for the loci under selection.
Neutral variation may also come to be associated with locally selected alleles,
so that the latter constitute a 'barrier' to neutral gene flow, increasing
neutral genetic differentiation (as quantified by $F_{ST}$ for instance) beyond the
single locus neutral expectation [@bengtsson1985].

Barrier effects due to divergent local adaptation at many loci may play an
important role in the evolution of reproductive isolation, and hence speciation
[@barton2020; @nosil2012].
The colonization of a new habitat by some species will often involve selection
on polygenic traits and give rise to a subpopulation that exhibits some
divergence from its ancestors [@barton2018].
Conditional on the initial succesful establishment of such a divergent
subpopulation through polygenic adaptation, whether or not speciation ensues
depends on the migration rate and the extent to which local adaptation can be
maintained in the face of maladaptive gene flow, and to what extent the partial
reproductive isolation deriving from local adaptation may promote further
divergence and strengthen reproductive isolation, either through reinforcement,
coupling with intrinsic incompatibilities, or the establishment of additional
locally beneficial mutations [@barton2009; @bierne2011; @butlin2018;
@kulmuni2020].
In this paper, we focus on the conditions under which polygenic local
adaptation can be maintained in the face of maladaptive gene flow.

Despite mounting evidence that local adaptation is indeed often polygenic in
nature [@bomblies2022], little is known about the underlying genetic details:
How many loci are involved? What are the typical effect sizes? Are locally
beneficial alleles typically closely linked or are they spread all over the
genome? How non-additive is local adaptation? *etc.* [e.g. @yeaman2011b;
@yeaman2015; @bomblies2022].
Moreover, even if such details were known, it would remain unclear to what
extent the genetic architecture of local adaptation affects the ability of a
population to maintain reproductive isolation in the face of gene flow. 
For instance, it is not directly clear whether a heterogeneous architecture
consisting of many loci of small effect and a couple large-effect loci would
admit more adaptive differentiation in the face of gene flow than a homogeneous
architecture with the same total selective effect.
The closely related, and heavily debated, question of how much scope there is
to infer the detailed genetic architecture underlying local adaptation from
observed patterns of genomic differentiation, as for instance obtained through
so-called 'genome scans', also remains largely unanswered.
So far, most theoretical developments have assumed rather simple genetic
architectures, dealing with biallelic loci of equal additive effect (ignoring
dominance and epistasis) that are either unlinked or uniformly spread along a
block of genome [@barton1983; @fraisse2021b; @sachdeva2022]; and statistical
approaches for the inference of gene flow across the genome either make
similarly crude assumptions [@aeschbacher2017], or ignore the genetic details
of local adaptation altogether [@fraisse2021; @laetsch2022].

In a recent paper, @sachdeva2022 showed that, when the loci under selection are
unlinked, the effects of LD on equilibrium differentiation at any individual
locus in a multilocus barrier can be well described by classical (deterministic
or stochastic) single locus population genetic theory, provided that the
migration rate $m$ is substituted by an *effective* migration rate $m_e$ which
depends, essentially, only on the total trait divergence [@petry1983;
@bengtsson1985; @barton1986; @kobayashi2008], and captures the effects of
selection against the associated genetic background on gene flow at a focal
locus.
The effective migration rate for a neutral locus can furthermore serve as a
quantitative measure of reproductive isolation [@westram2022].
In her paper, @sachdeva2022 conducted a detailed study of the effects of both
drift and LD on swamping thresholds and neutral differentiation in the
mainland-island and infinite-island models of population subdivision, assuming
a haploid sexual life cycle and $L$ divergently selected loci of equal effect.
The general theoretical framework outlined in @sachdeva2022 is however readily
extended to deal with more complicated genetic architectures.
In this paper, we derive an expression for the effective migration rate for a
polygenic genetic architecture with arbitrary fitness effects and dominance
across loci (referred to as a *heterogeneous* architecture or barrier) in a
population with a haplodiplontic life cycle (which includes haplontic and diplontic
life cycles as special cases).
We use this $m_e$ to build an approximation for the marginal allele frequency
distributions at migration-selection balance in a mainland-island model, and
make use of these tools to investigate the effects of dominance and variation
among fitness effects across loci on the maintenance of polygenic local
adaptation and the genetic architecture of reproductive isolation.


# Model and Methods

## Haplodiplontic mainland-island model {#sec:model}

Here we outline a mainland-island model for a sexual population which may be
subject to selection in both the haploid and diploid phases.
We think of this model as a caricature of a bryophyte, pteridophyte or fungal
population, but as we shall see below, the model encompasses both diplontic and
haplontic life cycles as well.
Throughout, we shall assume that sexes need not be distinguished, and that
selfing is possible.
We assume a regular and synchronous alternation of generations, where an island
population of $N$ haploids (gametophytes) produces an effectively infinite pool
of gametes which unite randomly to form $N k$ diploid individuals
(sporophytes).
The $Nk$ diploids produce in turn an effectively infinite pool of haploid
spores through meiosis, of which $N$ are drawn to form the next haploid
generation.
In each generation, we assume $M$ haploid individuals on the island are
replaced by haploid individuals from a mainland population, where $M$ is
Poisson distributed with mean $Nm$.
Fitness on the island is determined by $L$ unlinked biallelic loci which are
under divergent selection relative to the mainland.
The mainland population is assumed to have a constant, but arbitrary, genetic
composition.
Unless stated otherwise, we shall assume the mainland to be fixed for the
locally deleterious allele on the island.
Fitness effects are allowed to vary arbitrarily across loci.
Denoting the alleles at locus $i$ by $A_{i,0}$ and $A_{i,1}$, we designate by
$w_{i,j}$ the relative fitness of the haploid genotype $A_{i,j}$ and $w_{i,jk}$ the
relative fitness of diploid genotype $A_{i,j}A_{i,k}$.
We suppress the index $i$ when considering a generic locus.
We assume throughout that $w_0 = 1$ and $w_1 = e^{s_1}$ for the haploid phase,
and $w_{00} = 1, w_{01} = w_{10} = e^{s_{01}}$, and $w_{11} = e^{s_{11}}$ for
the diploid phase.
Throughout, we denote the frequency of the allele with relative fitness $1$ at
locus $i$ by $p_i$, and the frequency of the alternative allele by $q_i =
1-p_i$.
Fitness is determined multiplicatively across loci, so that, for instance, the
log relative fitness of a haploid individual fixed for all the '1' alleles
(genotype $A_{1,1},A_{2,1}, \dots, A_{L,1}$) is given by $\log w = \sum_{i=1}^L s_{i,1}$.
We assume that each haploid (diploid) individual contributes gametes (spores)
to the gamete (spore) pool proportional to its fitness.
We assume symmetric mutation at a small constant rate $u$ per locus, occurring
at meiosis.

Individual-based simulations of this model are implemented in Julia [@julia]
and the code is available at
[`https://github.com/arzwa/MultilocusIsland`](https://github.com/arzwa/MultilocusIsland).
In the following sections, we build up a theoretical approximation to this
fairly general multilocus model, roughly as in @sachdeva2022, and validate the
approximations by comparing numerical results against individual-based
simulations.
We first derive the dynamics at a single locus, considering both deterministic
and stochastic models.
Next, we derive an approximation to the effective migration rate for the
multilocus model using a rather general argument based on the reproductive
value of migrant individuals.
Lastly, we approximate the allele frequency dynamics of the multilocus model by
plugging in the effective migration rate, which captures the effect of LD among
selected alleles on the dynamics at a neutral locus, in the single locus
theory.


## Single locus mainland-island model

### Deterministic dynamics {#sec:sldet}

We first consider a deterministic model for the allele frequency dynamics at a
single focal locus, ignoring the influence of the other loci as well as genetic
drift. 
As shown in detail in @sec:app1, for weak selection and migration, the
dynamics of $p$ can be described in continuous time by the nonlinear ordinary
differential equation (ODE)
\begin{equation}
   \dot{p} := \frac{dp}{dt} = -m(p-p_M) -q(s_ap + s_bpq)\ ,
   \label{eq:ode}
\end{equation}
where $p_M$ is the frequency on the mainland of the locally beneficial allele
on the island, $s_a = s_1 + s_{01}$ and $s_b = s_{11} - 2s_{01}$, the latter
being a measure of dominance (i.e. the deviation from multiplicative fitnesses,
sometimes called $\iota$ [@otto2003; @manna2011]).
Usually, $s_1, s_{01}$ and $s_{11}$ will be assumed to be negative, and $p_M$
will be assumed to be small, so that selection increases $p$, whereas migration
decreases $p$. 
When $s_1 = 0$, we obtain the standard diploid mainland-island model, which is
commonly parameterized in terms of a dominance coefficient $h$ and selection
coefficient $s$, so that $s_{01} = sh$ and $s_{11} = s$.
When $s_1 \ne 0$ (i.e. there is selection in the haploid phase), this allows us
to define an effective selection coefficient $s_e = s_1 + 2s_{11}$ and
dominance coefficient $h_e = \frac{s_1 + s_{01}}{2s_1 + s_11}$, so that, at
least when selection is sufficiently weak, the single locus dynamics for an
arbitrary haplodiplontic life cycle can be described by the standard diploid
model with these effective parameters.
The equilibria of @eq:ode are analyzed in detail in @sec:mieq.


### Diffusion approximation to the stochastic dynamics

Still considering a single focal locus, we now account for the effects of
sampling drift.
Denoting by $X_n$ the number of $A_1$ copies in the $n$th haploid generation,
and by $Y_n^{(ij)}$ the number of $A_iA_j$ genotypes in the $n$th diploid 
generation, the life cycle as outlined in @sec:model entails the following
Markov chain model:
  \begin{align}
  Y_n|X_n &\sim \HW\left(Nk, \frac{w_{h,1}}{\bar{w}_h}\frac{X_n}{N}\right)
    \label{eq:mc1} \\
  X_{n+1}|Y_n &\sim \Bin\left(N, \frac{\bar{w}_{d,1}}{\bar{w}_d}
    \left(\frac{Y^{(11)}_n + Y^{(01)}_n/2}{Nk}\right)\right)\ , \label{eq:mc2}
  \end{align}
where $w_{h,1}\ (\bar{w}_h)$ and $\bar{w}_{d,1}\ (\bar{w}_d)$ are the marginal
fitnesses of the $A_1$ allele (mean fitnesses) in the haploid and
diploid generation respectively, and where $\HW(N,p)$ refers to the multinomial
distribution with Hardy-Weinberg proportions at allele frequency $p$ for
population size $N$.
Note that one unit of time corresponds to a single *alternation* of
generations, involving two sampling stages: first we sample diploid genotypes
from the random union of fitness-weighted haploid genotypes (@eq:mc1), next we
sample haploid spores by randomly drawing gene copies from the diploid
genotypes that survive diploid selection (@eq:mc2).
This model is akin to the standard Wright-Fisher (WF) model with variable
population size, regularly alternating between $N$ and $2Nk$ gene copies.
The corresponding effective population size is hence $N_e = (N^{-1} +
(2Nk)^{-1})^{-1}$, twice the harmonic mean of the phase-specific number of gene
copies [@hein2004] (twice because our unit of time is an alternation of
generations, not a single generation). 

When evolutionary forces are sufficiently weak, diffusion theory can be applied
to approximate the equilibrium allele frequency distribution implied by the
above Markov chain by a continuous density $\phi(p)$ on the unit interval.
Specifically, we have [e.g. @felsenstein2005; @ewens2004]:
\begin{equation*}
  \phi(p) \propto V(p)^{-1} \exp\left[2\int_0^p \frac{M(x)}{V(x)}dx\right],
\end{equation*}
where for the haplodiplontic mainland-island model, the infinitesimal mean and
variance will be, respectively,
  \begin{align*}
  M(p) &= -q(s_ap + s_bpq) + u(q - p) - m(p - p_M) \\
  V(p) &= N_e^{-1} pq\ ,
  \end{align*}
where we assume mutations occur with rate $u$ for both alleles. This yields a
probability density function for the equilibrium allele frequency distribution
  \begin{equation}
  \phi(p; N_e, u, m, s) \propto p^{2N_e(u+mp_M)-1}q^{2N_e(u+mq_M)-1}e^{N_e(2s_aq + s_bq^2)},
  \label{eq:phi}
  \end{equation}
where no closed form expression is known for the normalizing constant.  This is
essentially Wright's [@wright1937] distribution, generalized to a
haplodiplontic life cycle.


## Multilocus model {#sec:ml}

### Effective migration rate

To begin constructing a useful approximation to the allele frequency dynamics
in the multilocus system, we derive an expression for the effective migration
rate $m_e$, which captures the reduction in gene flow at a focal locus due to
selection against migrant genotypes.
As shown formally in @kobayashi2008, for weak migration, the reduction in gene
flow relative to the 'raw' migration rate $m$, termed the *gene flow factor*
(gff), depends on the expected reproductive value (RV) of migrants in the
resident background (i.e. the expected long-term contribution of a migrant
individual to the neutral gene pool on the island, relative to individuals in
the resident island population).
At any time, the proportion of individuals with recent migrant ancestry on the
island is $O(m)$, so that the probability of individuals with migrant
backgrounds mating with each other to produce, for instance, F2 crosses of the
migrant and resident genotypes, is $O(m^2)$, and hence negligible for
sufficiently weak migration.
The descendants of a migrant individual will therefore most likely correspond
to F1 and subsequent backcross generations, so that to a good approximation,
the RV of a migrant individual depends on the relative fitnesses of F1, BC1,
BC2, *etc.* individuals.

Let $W_{h,k}$ and $W_{d,k}$ denote the relative fitness of an individual
derived from a $k$th generation haploid, respectively diploid, backcross of an
initial migrant individual with the resident population (i.e. $W_{d,1}$ is the
relative fitness of an F1 diploid, $W_{d,2}$ of an offspring from a F1 $\times$
resident cross (BC1 generation), *etc.*).
Assuming migration occurs in the haploid phase before selection, the gff can be
expressed as
\begin{equation}
g = \frac{m_e}{m} = \Ex[W_{h,0}\prod^\infty_{k=1} W_{d,k}W_{h,k}],
\end{equation}
where $W_{h,0}$ is the relative fitness of the haploid migrant in the resident
population [@westram2018; @barton2018; @sachdeva2022].
Note that this involves an expectation over all possible lines of descent of an
initial migrant spore.
In practice, $g$ is determined only by the first 10 backcross generations or
so, as subsequent backcross generations are essentially indistinguishable from
residents.
In order to derive a useful approximate expression for $g$, we shall make two
further important assumptions:
(1) both the resident and migrant gene pool, as well as each backcross
generation is in Hardy-Weinberg and linkage equilibrium (HWLE);
(2) the expected allele frequency at any locus in any backcross generation is
midway between that of the parents (e.g. midway the mainland and island
allele frequencies for the F1 generation), with negligible deviations from the
midpoint due to selection within F1s, BC1s, *etc.*
This last assumption becomes more plausible as local adaptation is due to
more and more loci of smaller effect.

Under these assumptions, each of the $W_{\cdot,k}$ is determined solely
by the frequencies of the selected alleles in the mainland and the island
populations at the assumed equilibrium.
This allows us to determine $\Ex[q_{i,k}]$, the expected frequency of the
locally deleterious allele at locus $i$ among $k$th generation descendants from
a migrant, in terms of the allele frequencies in the mainland and island
population.
Indeed, assumption (2) implies the recursive relation $q_{i,k} =
\frac{1}{2}(q_{i,k-1} + q_i)$, where $q_i$ is the allele frequency in the
resident population -- i.e. the average number of selected alleles carried
by a $k$th generation backcross is the mean of the number of such alleles
carried by a $k-1$th generation backcross and a resident individual.
Hence, we have $\Ex[q_{i,k}] = \frac{1}{2^k}(q_{M,i} + (2^k - 1)\Ex[q_i])$,
where $q_{M,i} = q_{0,i}$ is the mainland frequency.
Denoting the selection coefficient at locus $i$ for the haploid phase by
$s_{i1}$, we can use this to derive the expected relative fitness of a $k$th
generation haploid descendant:
\begin{align*}
\Ex[W_{h,k}] 
  &= \Ex\left[\frac{\prod_{i=1}^L(p_{k,i} + q_{k,i}e^{s_{i1}})}
                   {\prod_{i=1}^L (p_i + q_i e^{s_{i1}})}\right] \\
  &= \frac{\prod_{i=1}^L\Ex\left[1 - q_{k,i}(1 - e^{s_{i1}})\right]}
          {\prod_{i=1}^L\Ex\left[1 - q_i (1 - e^{s_{i1}})\right]} \\
  &\approx \exp\left[\sum_{i=1}^L s_{i1}\Ex[q_{k,i} - q_i]\right] \\
  &= \exp\left[2^{-k}\sum_{i=1}^L s_{i1}(q_{M,i} - \Ex[q_i])\right]
\end{align*}
where we have assumed that per-locus selection is sufficiently weak so that
terms which are $O(s^2)$ can be ignored.
For the diploid phase, a similar argument shows that for the $(k+1)$th
generation,
\begin{align*}
\Ex[W_{d,k+1}] 
  &= \exp\left[2^{-k}\sum_{i=1}^L s_{i01}(q_{M,i} - \Ex[q_i]) 
    - s_{i,b}(p_{M,i}\Ex[q_i] - \Ex[p_iq_i])\right],
\end{align*}
where $s_{i01}$ and $s_{i11}$ are the selection coefficients against
heterozygotes and homozygotes at locus $i$ respectively, and where, analogous
to the single locus model, $s_{i,b} = s_{i11} - 2s_{i01}$.
Putting everything together, the approximate gff becomes
\begin{align}
g 
 &\approx \Ex[W_{h,0}]\prod_{k=1}^\infty \left(\Ex[W_{d,k}]\Ex[W_{h,k}]\right)
 \nonumber \\
 &=\prod_{k=0}^\infty \exp\left[ 
    2^{-k}\sum_{i=1}^L s_{i,a}(q_{M,i} - \Ex[q_i]) - 
    s_{i,b}(p_{M,i}\Ex[q_i] - \Ex[p_iq_i])\right] \nonumber \\
&=\exp\left[2\sum_{i=1}^L 
    s_{i,a}(q_{M,i} - \Ex[q_i]) - s_{i,b}(p_{M,i}\Ex[q_i] - \Ex[p_iq_i])\right] 
    \label{eq:gff}
\end{align}
where, similarly, $s_{i,a} = s_{i1} + s_{i01}$.
It is worth stressing that the gff is a function of the differentiation between
the mainland and island population as well as the heterozygosity $\Ex[pq]$ on
the island, and that, although we assume migration is sufficiently rare, we do
*not* assume that alleles introduced by migrants are rare.
We shall often highlight the dependence of the gff on the allele frequencies
and heterozygosities by writing $g[\Ex[p], \Ex[pq]]$, or $g[p]$ when allele
frequencies are known deterministically. 
If we assume all loci to have the same selection coefficient (a *homogeneous
barrier*), and that the mainland is fixed for the locally deleterious allele on
the island, the gff simplifies to
\begin{equation}
  g = e^{2L(s_a \Ex[p] + s_b \Ex[pq])}
  \label{eq:eqeff}
\end{equation}
where $\Ex[p]$ and $\Ex[pq]$ are the expected beneficial allele frequency and
expected heterozygosity at any selected locus on the island.
Note that in our applications $s_a \le 0$, so that when the locally beneficial
allele is common on the island ($p \approx 1$) and $Ls_a$ is appreciable, gene
flow will indeed be reduced ($g < 1$) due to selection against migrants. 

To gain some more intuition, we can express @eq:eqeff in terms of the effective
selection coefficient $s_e$ against the invading allele, and the effective
dominance coefficient $h_e$ of the invading allele over the locally beneficial
one, as
\begin{equation}
  g = e^{-2Ls_eh_e\Ex[p]}e^{-2Ls_e(1-2h_e)\Ex[pq]}
  \label{eq:eqeff2}
\end{equation}
(see @sec:sldet). Here, the first factor is just the gff associated with a
haploid $L$-locus system with selection coefficients $s_eh_e$ [@sachdeva2022].
The second factor captures the effects of dominance and depends on the
heterozygosity $\Ex[pq]$.
Clearly, $h_e$ has opposing effects on both factors.
The immediate effect of dominance is therefore that the gff is decreased
(barrier strength increased) relative to the additive case whenever invading
alleles exhibit a dominant deleterious effect on the island ($h_e > 1/2$).
Only when heterozygosity ($\Ex[pq]$) becomes appreciable does the second factor
contribute to the increase (when $h_e > 1/2$) or decrease (when $h_e < 1/2$) of
the gff.
The implications of these observations for the maintenance of adaptive
differentiation will be explored in detail in the results section.

Two remarks are due. Firstly, the gff as derived above yields the effective
migration rate at an unlinked neutral locus.
If we wish to calculate the effective migration rate for a selected locus in
the barrier, say locus $j$, the relevant gff is
\begin{equation}
    g_j = \exp\left[2 \sum_{i\ne j}^L 
    s_{i,a}(q_{M,i} - \Ex[q_i]) - s_{i,b}(p_{M,i}\Ex[q_i] - \Ex[p_iq_i])\right] 
    \label{eq:gff2}
\end{equation}
where we have assumed that the gff at a selected locus is the same as that of a
neutral locus at the same location -- an assumption which is only expected to
work well if selection at the focal locus is sufficiently weak.
Secondly, as in the model outlined in @sec:model, we have assumed that
migration occurs at the start of the haploid phase, reflecting a process such
as spore dispersal in bryophytes, pteridophytes or Fungi.
However, it should be emphasized that, while the details of when migration
occurs in the life cycle do not matter for the single locus model as long as
selection and migration are sufficiently weak (so that the continuous-time
limit is appropriate), these details *do* matter for the effective migration
rate.
This is because, although selection *per locus* is weak ($s$ being small),
selection against migrant genotypes can be strong ($Ls$ being appreciable).
Two other cases should hence be considered. Firstly, when migration is due to
dispersal of gametes (e.g. pollen dispersal), the first generation experiencing
selection on the island will be the diploid F1 generation, so that the
appropriate gff under the same approximation is $g/\Ex[W_{h,0}]$.
Secondly, when migration occurs at the beginning of the diploid phase (e.g.
seed dispersal), the first generation experiencing selection on the island will
consist of diploid migrant individuals, so that $g\Ex[W_{d,0}]$ is the
appropriate gff, where
\begin{equation*}
  \Ex[W_{d,0}] \approx
    \frac{e^{\sum_i^L 2p_{M,i}q_{M,i}s_{i,01} + q_{M,i}^2
        s_{i,11}}}{e^{\sum_i^L2\Ex[p_iq_i]s_{i,01} + \Ex[q_i^2]s_{i,11}}}
    = \exp\left[\sum_{i=1}^Ls_{i11}(q_{M,i} - \Ex[q_i]) - 
    s_{b,i}(p_{M,i}q_{M,i} - \Ex[p_iq_i])\right]
\end{equation*}
If the haploid, diploid and gametic migration rates are $m_1, m_2$ and $m_3$
respectively, the effective migration rate will be $(m_1 + \Ex[W_{d,0}] m_2 +
\Ex[W_{h,0}]^{-1}m_3)g$.
Unless stated otherwise, in the present work, we shall always assume life
cycles in which migration is due to dispersal of haploid spores, so that
@eq:gff gives the relevant gff.


### Dynamics and equilibria for the multilocus model {#sec:dynamics}

The gff captures the effect of LD among selected loci on the rate of gene flow
from the mainland into the island at any individual locus.
The key observation is that a certain separation of time scales applies:
although selection against migrant *genotypes* can be very strong in the
polygenic case (of magnitude $Ls$, roughly), selection at any individual locus
is still assumed to be weak, so that, when linkage is weak or absent, after an
evolutionarily short period in which entire sets of alleles are efficiently
removed together, LD among selected loci becomes negligible and the standard
single locus theory should be applicable.
Hence, on the longer time scales at which migration-selection balance is
attained, the allele frequency dynamics at any individual locus should
essentially follow the single locus dynamics, but where migrant alleles are
introduced at a rate reduced by a factor equal to the gff [@sachdeva2022].

As a consequence, in the deterministic case, we expect that the effects of LD
should be well captured by substituting the effective migration rate $m_e = mg$
for $m$ in @eq:ode.
Specifically, we get a system of $L$ coupled differential equations, where for
$1 \le j \le L$,
\begin{equation}
   \dot{p_j} = -m g_j[p_{-j}]p_j - q_j(s_{a,j}p_j + s_{b,j}p_jq_j)\ ,
   \label{eq:odeme}
\end{equation}
and we have assumed the mainland to be fixed for the deleterious allele
on the island at all loci.
Here we write $g_j[p_{-j}]$ for the gff as in @eq:gff2, to highlight the
dependence of the gff at locus $j$ on the allele frequencies at the $L-1$ other
loci.
Note that in the deterministic model, the expected values in @eq:gff2
disappear, i.e. at any time $\Ex[q_j] = q_j$ and $\Ex[p_jq_j] = p_jq_j$. 
We study the equilibria of this model by numerically solving for $p$ at
stationarity ($\dot{p}_j = 0$, for $1 \le j \le L$).

As in @sachdeva2022, we can also plug in $m_e$ in the single locus diffusion
approximation to determine the equilibrium allele frequency distribution for
each locus on the island.
Specifically, we postulate that the joint distribution of allele frequencies in
the barrier factorizes as
\begin{align}
  \phi(p) = Z^{-1} \prod_{j=1}^L \phi_j(p_j|p_{-j}) = 
    Z^{-1}\prod_{j=1}^L \phi(p_j; N_e, u, mg_j[p_{-j}], s_j)
  \label{eq:mrf}
\end{align}
where $Z$ is a normalizing constant and $\phi$ was defined in @eq:phi.
@Eq:mrf can be thought of as the distribution associated with a Markov random
field over the complete graph with $L$ vertices.
The marginal allele frequency distribution at any particular locus depends on
the allele frequencies at the other $L-1$ loci.
We can compute moments of the allele frequency distribution at each locus by
solving the whole system self-consistently, that is, by assuming 
\begin{align*}
  \Ex[p_j] &= Z_j^{-1}\int p_j \phi(p_j, N_e, u, mg_j[\Ex[p_{-j}],
      \Ex[pq_{-j}]],s_j) dp_j \\
  \Ex[p_jq_j] &= Z_{j'}^{-1}\int p_j q_j \phi(p_j, N_e, u, mg_j\left[\Ex[p_{-j}],
      \Ex[pq_{-j}]\right], s_j) dp_j,
\end{align*}
where the $Z$'s are again normalizing constants, and $\Ex[pq_{-j}]$ is the
vector of expected heterozygosities at all loci excluding locus $j$ (i.e.
$(\Ex[p_1q_1], \dots \Ex[p_{j-1}q_{j-1}], \Ex[p_{j+1}q_{j+1}], \dots,
\Ex[p_Lq_L])$.
To solve this system of $2L$ nonlinear equations, we use the fixed point
iteration outlined in @sec:fp.
The numerical methods used in this paper are also implemented in the Julia
package available at 
[`https://github.com/arzwa/MultilocusIsland`](https://github.com/arzwa/MultilocusIsland).


# Results

The results section is organized as follows: we start with an analysis of the
deterministic multilocus system, contrasting the predicted barrier effect with
the single locus model and examining the effects of dominance, neglecting drift
and heterogeneity of selective effects across the barrier.
We next consider the effects of genetic drift and show how we can accurately
predict the equilibrium allele frequency distributions at individual barrier
loci with arbitrary dominance and haploid selection.
Lastly, we study the effects of heterogeneous barrier architectures on the
maintenance of local adaptation and observable patterns of equilibrium
differentiation.

## Multilocus barriers with dominance in the deterministic model

We first analyze the deterministic multilocus model for a homogeneous barrier.
We shall assume the mainland to be fixed for the locally deleterious allele (in
the island habitat) at all loci.
We only consider diploid selection here, with $s_{01} = -sh$ and $s_{11}
= -s$, where $s$ is the selection coefficient against the locally deleterious
allele, and $h=s_{01}/s_{11}$ is the dominance coefficient.
We emphasize that $h$ measures dominance *of the mainland (invading) allele
over the island (locally beneficial) allele*, so that $h=1$ corresponds to a
situation where the invading allele is fully dominant, or, equivalently, where
the allele that confers local adaptation on the island is recessive.
For the case where migration is at the haploid stage, restricting the analysis
to diploid selection incurs no loss of generality, as haploid selection then
simply amounts to a rescaling of the dominance and selection coefficients (see
methods).
The effective (diploid) dominance coefficient when there is haploid selection
with intensity $s_1$ will be $h_e = \frac{s_1+s_{01}}{2s_1+s_{11}}$, so that
the effect of haploid selection is to pull $h_e$ towards the additive case
($h=1/2$) where selection acts on each gene copy independently, as it does in a
strictly haploid model.

### The impact of dominance on barrier effects

![
Recessive local adaptation leads to stronger multilocus barriers to gene flow.
Equilibrium frequencies ($\tilde{p}$) of the locally beneficial alleles are shown for
increasing $Ls$ for \uline{(A)} the case of dominant local adaptation (recessive
migrant alleles), \uline{(B)} additive fitness effects and \uline{(C)} recessive local
adaptation.
Note that the mainland is fixed for the alternative allele, so that $\tilde{p}$
corresponds to the differentiation at equilibrium.
The thick lines show the stable equilibria for increasing $m/s$, whereas the
dotted lines show unstable equilibria.
The black dots mark the critical point beyond which swamping is obtained for
any initial condition (in (C) the approximate expression from @sec:supdet is
used).
The results for $Ls=0.01$ (gray lines) correspond to the single locus
predictions.
\uline{(D)} Swamping thresholds for different degrees of dominance (colors, see
(E) and (F)) for increasing total barrier strength $Ls$. 
\uline{(E, F)} Barrier strength $b[q] = g[q]^{-1}$ as a function of the
deleterious allele frequency $q$ on the island for different degrees of
dominance, for $Ls=0.75$ and $Ls=1.5$ respectively. The vertical dotted lines
mark the level of differentiation beyond which the barrier strength starts to
decrease when $h < 1/3$.
](/home/arthur_z/vimwiki/build/img/2023-07-17/detdom.svg){#fig:detdom}

In the homogeneous deterministic model, all loci have the same dynamics if the
initial allele frequencies are equal.
From @eq:odeme, we find that in that case, the frequency $p$ of the locally
beneficial allele at any selected locus in an $L+1$ locus system must satisfy
at equilibrium
\begin{align}
   0 &= shpq + s(1-2h)pq^2 -mg[p]p 
   \label{eq:odeq} \\
   &\text{where }\quad g[p] = e^{-2Ls(hp + (1-2h)pq)}. \nonumber
\end{align}
We can solve this numerically for the equilibrium allele frequency.
Note that if we set $g[p] = 1$, we recover the classical diploid single locus
model, for which the equilibrium behavior is well understood [@haldane1930VI;
@nagylaki1975].
We briefly recapitulate the main results for the single locus model (see also
@sec:mieq, and the gray lines in @fig:detdom).
For the case $h=0.5$ (no dominance, also referred to as codominance, or
additivity), the equilibrium frequency $\tilde{p}$ of the locally beneficial
allele decreases linearly from $1$ to $0$ as the rate of migration
approaches the strength of selection on a single allele $s/2$.
When local adaptation is due to a dominant allele (so that the invading allele
acts recessively to reduce fitness on the island, i.e. $h=0$), the migration
rate beyond which no polymorphism can be maintained is increased to $s$, while
$\tilde{p}$ is decreased relative to the additive case as long as $m <
s/4$.
When local adaptation is due to a recessive allele ($h=1$), the model has two
equilibria for the beneficial allele frequency, one stable equilibrium
$\tilde{p}_+ > 1/2$ and one unstable equilibrium $\tilde{p}_- < 1/2$ as long as
the migration rate does not exceed $s/4$. When this critical threshold is
passed, swamping occurs for any initial frequency.
Hence, for the recessive case, whether or not a polymorphism is attained
depends not only on the migration rate (which should be at most $s/4$), but
also on the history of the population: the island population cannot fix a new
recessive beneficial variant, but an established recessive variant can be
maintained upon secondary contact. Similar bistable behavior occurs for partial
recessivity as long as $h > 2/3$.

We now consider the equilibrium behavior in the multilocus case, where LD will
cause deviations from these single locus predictions.
@Fig:detdom shows the equilibrium behavior for a number of example parameter
sets.
As expected, stronger net selection against maladapted (migrant) genotypes
(larger $Ls$) increases the equilibrium frequency of the locally beneficial
allele relative to the single locus prediction, but the magnitude of this
effect depends quite strongly on dominance.
When the invading alleles are recessive ($h=0$; @fig:detdom A), gene flow is
not at all impeded when deleterious alleles are rare on the island (the gff
being near one; @fig:detdom E, F).
This is essentially because, irrespective of how many deleterious alleles a
migrant carries, if the deleterious alleles are rare on the island they will
not be found in homozygotes as long as migration is sufficiently weak, and
hence will not be 'seen' by selection.
Only once deleterious alleles are segregating at appreciable frequencies on the
island, will F1, BC1, *etc.* individuals be likely to be homozygous at several
loci, thus exposing (partially) recessive invading alleles to selection and
reducing the RV of migrants.
As a result, when invading alleles at the selected loci act recessively, a
strong genome-wide barrier effect emerges only once differentiation falls below
a critical threshold. 
The situation is clearly different when migrant alleles are dominant, as the
invading alleles will immediately express their full load in the resident
population, irrespective of the other allele at the locus, yielding efficient
selection against migrant alleles (the gff being at its minimum when migrant
alleles are rare, @fig:detdom E, F).
Any increase in the frequency of the deleterious allele on the island will
merely increase the expected relative fitness of migrants in the resident
background, and hence reduce the efficiency of selection against migrant
genotypes.
We observe a transition between these two qualitatively different types of
behaviour at intermediate values of $h$: when $h<1/3$, the barrier strength, as
measured by $g^{-1}$ [@barton1986], *increases* as the locally deleterious
allele increases in frequency on the island (and hence as differentiation
between mainland and island *decreases*), decreasing the rate of gene flow,
until a value of $q=(3h-1)/(4h-2)$ is reached (@fig:detdom, E, F).


### Effect of dominance on swamping thresholds

In the single locus model, arbitrarily small frequencies of the locally
beneficial allele can be maintained at migration-selection balance when $h <
2/3$, whereas in the case of $h > 2/3$, a sharp swamping threshold is observed
as the equilibrium frequency of the locally beneficial allele reaches some
critical frequency $p_c \le 1/2$ (@sec:mieq, gray lines in @fig:detdom).
@sachdeva2022 showed that such sharp thresholds for swamping also appear in
the absence of dominance due to LD.
LD both increases the critical migration rate ($m_c$) at which swamping occurs
and the minimum level of differentiation that can be maintained before the
swamping point is reached ($p_c$).
Our results indicate that dominance has a considerable influence on how LD
sharpens and displaces swamping thresholds (@fig:detdom D, @sec:supdet).
As $Ls$ increases, the behavior for the case where $h < 2/3$ (i.e. when local
adaptation is not strongly recessive) is roughly similar to the additive case,
with critical behavior emerging as $Ls$ surpasses some critical value
(@sec:supdet, @fig:pplot).
However, in this case, the critical migration rate at which swamping occurs is
only marginally affected by LD for moderate levels of divergence ($Ls < 2$,
say).
This is in sharp contrast with the case where local adaptation is due to
strongly recessive alleles ($h > 2/3$), where the critical migration increases
rapidly with $Ls$ (@fig:detdom D).
Importantly, the critical differentiation ($p_c$) level below which local
adaptation collapses is very different for different degrees of dominance.
In the additive case, one can show that critical behavior emerges as soon as
$Ls > 1$ (@sec:supdet), in which case $p_c = 1-1/Ls$, and hence arbitrary
differentiation can be maintained near the critical point depending on $Ls$.
For completely dominant local adaptation ($h=0$), however, $p_c$ increases from
$0$ to a maximum of $1/2$ as $Ls \rightarrow \infty$, whereas for recessive
local adaptation ($h=1$), $p_c$ increases from $1/2$ to $1$ as $Ls$ grows.
This means, in particular, that for moderate levels of divergence, $Ls > 0.75$
say, and large population sizes, one would not expect to see locally beneficial
recessives at frequencies much below $0.8$, compared to $0.5$ for the single
locus model (@fig:detdom C).
@Fig:detdom further highlights the nontrivial feedbacks between the observed
differentiation and dominance: whereas a consideration of the gff in a regime
where migrant alleles are rare would suggest that swamping thresholds and
equilibria depend roughly on $Lsh$, and not on $Ls$ and $h$ separately, this
intuition really only works well for very small rates of migration (@fig:lsh).


## Accounting for drift 

While the deterministic analysis points towards important effects of dominance
on equilibrium differentiation and thresholds for swamping when local
adaptation is polygenic, it is important to assess to what extent these carry
over to finite populations.
Indeed, @sachdeva2022 showed that, for small effective population sizes, the
sharp thresholds for swamping predicted by the deterministic multilocus theory
(when $Ls$ is appreciable) need not apply, and that the critical migration rate
may be significantly reduced.
Furthermore, for all but the largest populations, the observed frequency of a
locally adaptive allele (and hence differentiation between the mainland and
island) at migration-selection-drift balance in the model (as outlined in
@sec:model) may deviate substantially from the deterministic approximation, so
that it becomes important to understand the *distribution* of allele
frequencies on the island.

![
Genetic drift reduces the strength of a homogeneous multilocus barrier to gene
flow.
\uline{(A)}
Predicted allele frequency distributions for a single locus in the diploid
multilocus model with homogeneous selective effects for different values of
$N_e$ and $h$ (dominance coefficient of the invading alleles).  Lines show the
numerical approximations based on the diffusion theory, dots show results from
individual-based simulations, based on taking a sample every 10 generations for
50000 generations after an initial 10000 generations to reach equilibrium.
In these simulations, $s = 0.02$, $Ls = 1$ and $u/s = 0.005$.
\uline{(B)}
The expected frequency of the locally beneficial allele ($\Ex[p]$)
is shown as a function of $N_es$ and $m/s$ for different values of $Ls$ (from
\uline{left} to \uline{right}, $Ls = 0.5, 1, 1.5, 2$) and $h$ (from \uline{top}
to \uline{bottom}, $h=0, 0.5,
1$), computed using the numerical approximation.
All results assume $s=0.02$ and $u/s=0.005$. 
](/home/arthur_z/vimwiki/build/img/2023-06-17/drift-homo-ep-af.svg){#fig:drifthm}

Similar to @sachdeva2022, we find that substituting $m_e$ for $m$ in the
single-locus diffusion theory (@sec:dynamics) yields a remarkably accurate
approximation to individual-based simulations of the model.
Indeed, even in parameter regimes where the approximation is expected to break
down ($Ls$ appreciable with $L$ small and $s$ large, small population size) we
obtain good predictions (@fig:Lsdom).
Not only can we reliably obtain the expected frequency of alleles on the
island, we also obtain very good predictions for the entire allele frequency
distribution (@fig:drifthm A).
The sharp swamping thresholds observed in the deterministic model, in
particular with recessive local adaptation ($h=1$), correspond to strongly
bimodal allele frequency distributions in the stochastic model.
As described in more detail in @sec:init, this may render our numerical
approaches sensitive to the assumed initial state of the island population.
This sensitivity is itself biologically relevant, corresponding to assumptions
on the (recent) history of the populations considered at equilibrium.
Throughout, we shall assume a scenario of secondary contact, so that the island
starts as fixed for the locally beneficial allele at all loci. 

For a homogeneous genetic architecture, the swamping threshold at any
individual locus depends on the total barrier strength $Ls$, the dominance
coefficient $h$, and the strength of selection per locus relative to drift
$N_es$.
Concomitantly, for any given migration rate $m$, these three key parameters
will determine jointly whether significant adaptive differentiation is to be
expected at a particular locus, as we show in @fig:drifthm (B).
Unsurprisingly, the general consequence of genetic drift is to reduce the
barrier effect, hence decreasing the expected differentiation at equilibrium
(@fig:drifthm B; @fig:drift).
Swamping thresholds are both decreased and made less sharp by drift, and the
detailed behavior depends on the dominance coefficient.
Whereas for dominant local adaptation the effect of drift on the critical $m/s$
value decreases markedly for $N_es > 20$ (@fig:drifthm B, top row), for the
additive case this happens, roughly, for $N_es > 10$ (@fig:drifthm B, middle
row).
For substantial $N_es$ (e.g. $N_es=16$ in @fig:drift), we now see clearly that
for dominant local adaptation (invading alleles are (partially) recessive,
$h<0.5$), two phases can be distinguished as we consider increasing migration
rates.
As noted above, the barrier strength initially increases with migration as
deleterious alleles increase in frequency on the island, exposing more invading
recessives to selection.
However, above a certain allele frequency (determined by $h$), the barrier
strength starts to decline again with increasing $m$ (@fig:drift, $h=0$ and
$h=0.25$ panels).
As in the deterministic analysis above, this behavior is in sharp contrast with
the case of recessive local adaptation (@fig:drifthm B, bottom row), where the
positive feedback between the gff and the frequency of the locally deleterious
allele ($q$) gives rise to sharp thresholds for swamping even when drift is
quite strong (e.g. $N_es = 4$).


## Selection in both the haploid and diploid phase


We now consider in more detail the accuracy of our approximations when
selection acts both in the haploid and diploid phase.
We parameterize the general model in such a way that we can investigate, for a
given total barrier strength 
(corresponding to the relative fitness of a migrant individual in an otherwise
perfectly adapted island population),
the effects of the relative strength of selection in the diploid and the
haploid phase and the degree of dominance ($h$) in the diploid phase.
To this end, we assume
\begin{align}
   s_1 = -(1-\tau)s & & s_{01} = -h\tau s & & s_{11} = -\tau s,
   \label{eq:hapdipmodel}
\end{align}
where $0 \le \tau \le 1$ measures the relative strength of selection in the
diploid phase (if one assumes selection to operate continuously with constant
intensity throughout the life cycle, this can be interpreted as the relative
length of the diploid phase).
Similar models have appeared in the study of life cycle modifiers [e.g.
@otto1994; @scott2017].
Recall furthermore that we assume a regular alternation of $N$ haploid and $Nk$
diploid individuals.
We find that we can indeed accurately predict equilibrium allele frequencies
for haplodiplontic life cycles with selection in both phases using the
diffusion theory (@fig:domtau).
We emphasize that the diffusion approximation depends only on $N, k, s, h$ and
$\tau$ through $N_es_e = N_e(2-\tau)s$ and $h_e = s_e^{-1}(1 - \tau(1-h))$,
where $N_e = (N^{-1} + (2Nk)^{-1})^{-1}$, showing that, at least for weak
selection, life cycle details can be accounted for by means of a set of
suitable effective parameters (@fig:hapdipeff).

![
Multilocus migration-selection balance and swamping with dominance and
haploid selection.
For a given total barrier strength ($Ls$, see main text), we vary the relative
strength of selection in the diploid and haploid phase ($\tau$, colors) and the
degree of dominance (columns) in the diploid phase.
Specifically, we assume $s_1 = -(1-\tau)s, s_{01} = -h\tau s$ and $s_{11} =
-\tau s$. Hence, $\tau = 1$ corresponds to a diplontic life cycle (or at least,
absence of selection in the haploid phase), whereas $\tau=0$ corresponds to a
haplontic life cycle.
Other parameters were as follows: $Ls =0.8, L=40, N_es=8, k=5, u=s/100$.
The dots show results from individual-based simulations, whereas the lines were
computed using the numerical approximation based on diffusion theory. 
](/home/arthur_z/vimwiki/build/img/2023-04-19/domtau.svg){#fig:domtau}

The results displayed in @fig:domtau suggest furthermore that, for a given
total strength of selection $Ls$, predominantly haploid populations should be able
to maintain more adaptive variation, and exhibit stronger reproductive
isolation, irrespective of the degree of dominance in the diploid phase.
Although we have shown that recessive local adaptation in diploids leads to
stronger barriers, increasing the relative strength of haploid selection, and
thereby pulling the effective dominance coefficient closer to $0.5$, does not
lead to weaker barriers when local adaptation acts recessively in the diploid
phase.
This is because the effective selection coefficient in this model is
$(2-\tau)s$, so that the strength of selection per gene in haploids is twice
that in diploids in the absence of dominance.
The relevance of these observations for the evolution and maintenance of
haplodiplontic life cycles is however not very clear, as a life cycle modifier
need not keep the overall strength of selection constant [@scott2017].
Lastly, we note that the phase of the life cycle in which migration occurs can
have an impact on the ability to maintain adaptive differentiation.
When there is selection in both phases, and migration occurs in the diploid
phase, maladapted diploid genotypes are directly exposed to selection on the
island, whereas this is not the case when migration occurs in the haploid
phase.
In the latter case, the first diploid generation exposed to selection is an F1
cross between the mainland and island, so that there is no selection on the
homozygous diploid effect.
We would therefore expect that migration in the diploid phase generally leads
to stronger barriers to gene flow when selection acts in both phases, and this
is indeed observed (@fig:dipmig).


## Heterogeneous genetic architectures

We now depart from the unrealistic assumption of equal effects and dominance
coefficients across the polygenic barrier.
In @sec:ml, we developed the multilocus theory for potentially heterogenous
(unlinked) genetic architectures, where the selection coefficients $s_1, s_{01}$
and $s_{11}$ can vary arbitrarily across loci, and we verify that we do indeed
obtain accurate predictions also in this setting (@fig:het).
This allows us to address in more detail a number of questions pertaining to
the genetic architecture of local adaptation at migration-selection balance,
accounting for both LD and genetic drift.
The scenario we envisage is the following: at some point in time, the island is
colonized by individuals from the mainland and rapid adaptation to the local
conditions from the standing genetic variation ensues (by driving $L$ locally
beneficial alleles to high frequencies on the island).
After this initial idealized phase of directional selection on standing
variation, we assume a situation of secondary contact, where, on average, $Nm$
haploid individuals on the island are replaced by mainland individuals in each
generation, and the island evolves to an equilibrium state.
We then ask what sort of loci can contribute to local adaptation and how the
resulting partial reproductive isolation is manifested in observable allelic
differentiation.
We also consider to what extent selective interference among loci results in
departures from single locus predictions. 
We emphasize that we do not explicitly consider the buildup of divergence
between the mainland and island population, nor some plausible model for what
sort of standing variation is the source of the initial polygenic response, but
merely ask how the genetic architecture of local adaptation affects observable
patterns of adaptive differentiation.
All results discussed in the remaining sections are obtained using the
numerical approximations based on diffusion theory.

### Effect of an additive polygenic barrier on a focal selected locus

![
Polygenic barriers can protect recessive alleles from swamping.
\uline{(A)} Expected equilibrium frequency of a recessive locally beneficial
allele ($h=1$) with selective effect $s_f$ in the presence of an additive
polygenic background of strength $Ls$.
The heavy lines show the expected allele frequency at the focal selected locus
for increasing strength of migration and different strengths of selection
against the additive background ($Ls$, colors), whereas the transparent lines
show the same for a locus in the additive background.
The black line shows the single locus prediction (i.e. in the absence of the
additive polygenic barrier).
The vertical dotted lines mark the $m/s_f$ values shown in (B) and (C).
\uline{(B, C)} The associated allele frequency distributions at equilibrium for
the focal locus when selection (at the focal locus, $s_f$) is six,
respectively four, times as strong as migration ($m$). 
\uline{(D, E, F)} As in (A, B, C) but for the case of a dominant locally
adaptive variant ($h=0$).
We assume $N_e = 1000, L=50$ and $u/s=0.005$ throughout.
](/home/arthur_z/vimwiki/build/img/2023-07-17/focal.svg){#fig:focal}

Firstly, we ask to what extent the equilibrium differentiation at a single
focal locus is affected by LD in the presence of a polygenic barrier of some
strength.
In particular, single locus theory indicates a strong effect of dominance on
the ability to maintain local adaptation for a given migration rate.
For instance, as already remarked by @haldane1930VI, a single recessive allele
conferring local adaptation is rather sensitive to swamping unless its
selective effect is relatively large (at least $4m$ in an infinite population).
However, in the more plausible case where local adaptation has at least a
partly polygenic basis, it becomes important to understand to what extent LD
with other barrier loci assists in maintaining adaptive differentiation and
increasing swamping thresholds at such loci.
Indeed, recent findings such as those of @stankowski2022 emphasize that readily
discovered large-effect loci may often be associated with a polygenic
background that is less clearly associated with divergent selection.
In @fig:focal, we consider the cases of a recessive ($h=1$) and dominant
($h=0$) locally adaptive variant of relatively large homozygous effect ($s_f =
0.05$) in the presence of an additive multilocus barrier with total selective
effect $Ls$.
As $Ls$ becomes appreciable, equilibrium frequencies at the focal locus start
to deviate substantially from the single locus prediction whenever the
migration rate is below the swamping threshold of the loci composing the
additive polygenic 'background' barrier (@fig:focal A, D).
In the case with recessive local adaptation, LD also strongly affects the
critical migration rate at the focal locus (@fig:focal A), and once the
swamping threshold for the additive background loci exceeds that of the focal
locus (in this example at $Ls \ge 1$), the latter is essentially increased to
that of the additive background loci.
Clearly, Haldane's (1930) remark that, for a recessive, "if $s$ is not much
greater than $4m$, a chance fluctuation may push the population past the point
of unstable equilibrium, and the recessives be finally eliminated" does not
straightforwardly apply in the polygenic setting, where LD can have a
considerable protective effect and render recessives more robust to
swamping, even when selected loci are unlinked (@fig:focal B & C).


### Effect of variation in fitness effects across a polygenic barrier

![
Variation among selective effects weakens polygenic barriers.
\uline{(A)} The boxplots
show the mean per-locus differentiation ($\bar{\Delta} = \sum_{i}^L
\Ex[p_i]/L$) across the $L$-locus barrier ($L=100$), for 50 replicate
simulations of an additive polygenic barrier where selection coefficients are
distributed according to a $\Gam(\kappa, \kappa/\bar{s})$ distribution, with
$\Ex[s] = \bar{s} = 0.01$ and six different values of $\kappa$
(note that $\mathrm{Var}[s] = \bar{s}/\kappa^2$).
The solid horizontal line shows the predicted equilibrium differentiation per
locus for a homogeneous barrier of strength $L\bar{s}$, whereas the dashed line
shows the single locus prediction for a diploid locus with $s = \bar{s}$.
The colored lines show the average differentiation across the barrier predicted
using single locus theory $\sum_i^L{\Ex[p_i|s_i]}/L$, averaged over replicate
simulations.
\uline{(A, inset)} Density functions for the six different Gamma distributions
used in (A).
\uline{(B)} Expected beneficial allele frequencies across the barrier in a single
simulation replicate for each of the six assumed distributions, sorted by
allele frequency, assuming $m/\bar{s} = 0.1$ (horizontal lines as in (A);
colors as in (A) and (B)). Other parameters are $N_e\bar{s} = 10$,
$u/\bar{s} = 0.005$.
\uline{(C)} As in (A) but for $m/\bar{s} = 0.4$. 
](/home/arthur_z/vimwiki/build/img/2023-07-14/gammas.svg){#fig:gammas}

We now turn to fully heterogeneous barriers, where all $L$ loci have
different fitness effects.
We first consider the case with variable selection coefficients, assuming
no dominance.
In @fig:gammas (A), we show the average expected differentiation per selected
locus ($\bar{\Delta}$) in polygenic barriers with selection coefficients
sampled from a Gamma distribution ($L=100, L\bar{s}=1$).
When migration is weak relative to selection (roughly $m/\bar{s} < 1/4$),
increasing the variance in fitness effects, while keeping $\Ex[s] = \bar{s}$
constant, yields on average weaker barriers, with lower equilibrium
differentiation than a homogeneous barrier of strength $L\bar{s}$.
Often the average per-locus differentiation is even lower than that of a single
locus selected against with intensity $\bar{s}$ (@fig:gammas A, $m/\bar{s} =
0.05, 0.1, 0.2$).
One should be careful, however, in the interpretation of $\bar{\Delta}$.
As shown in @fig:gammas (B, C), differentiation across loci in the barrier
often shows a strongly sigmoidal pattern, especially when $\Var[s]$ is large,
where most loci are either strongly differentiated or not at all, and with
rather few loci having $\Ex[p]$ near $\bar{\Delta}$.
This entails that empirically, instead of detecting $L$ selected loci with an
average differentiation of $\bar{\Delta}$, we are more likely to observe about
$L\bar{\Delta}$ strongly differentiated loci.
This highlights the difference between the potential and realized genetic
architecture of local adaptation (see also below): of the $L$ loci under
divergent selection, only about $L\bar{\Delta}$ contribute to local adaptation
and reproductive isolation at migration-selection balance.

As the migration rate becomes larger, a shift occurs, and heterogeneous
barriers tend to yield higher equilibrium differentiation than a homogeneous
one with the same average effect.
This is largely driven by a subset of loci with relatively large $s$ that
resist swamping (@fig:gammas C).
The effect of increasing the variance appears less substantial in this regime.
Given that the distribution of selection coefficients is generally believed to
be at least somewhat leptokurtic, these results suggest that heterogeneity in
selection coefficients has important consequences for observable
differentiation at migration selection-balance, which would not adequately be
captured by substituting an average selection coefficient in either single
locus or multilocus theory.
Keeping selection coefficients fixed while sampling dominance coefficients from
a Beta distribution with mean $1/2$ and increasing variance has less dramatic
effects, although we do see systematic increases and decreases in equilibrium
differentiation depending on whether the migration rate exceeds the swamping
threshold for recessives (which are associated with higher equilibrium
frequencies) or not (@fig:betas).

![
\uline{(A)} Deviation of predicted allele frequencies for loci in heterogeneous
polygenic barriers when accounting for LD (multilocus) from predictions based
on single locus diffusion theory. 
The \uline{rows} show results for different total strengths of selection (different
number of loci $L\bar{s}$ with $\bar{s} = 0.01$), whereas the \uline{columns} show
results for increasing rates of migration relative to selection ($m/\bar{s}$).
We assume the $s_i$ to be exponentially distributed with mean $\bar{s}$ and
dominance coefficients are sampled uniformly from the $[0,1]$ interval.
Each dot is associated with a single locus in an $L$-locus barrier, and is
colored according to its dominance coefficient (yellow for locally beneficial
recessives ($h=1$), purple for dominants ($h=0$)).
Each plot shows results for 1000 such loci, subssampled from a total of
$150000/L$ simulations of $L$-locus barriers.
\uline{(B)} Monte Carlo approximation to the marginal distribution of the selection and
dominance coefficient conditional on observing a divergent allele on the island
(i.e. $f(s_i|X_i=1)$ and $f(h_i|X_i=1)$, see @eq:msbdfe). The distribution graphed
in gray shows $f_\text{DFE}$, i.e. the marginal distribution of the selection
and dominance coefficient for a random locus in the $L$-locus barrier (not
considering migration).
We assumed $N_e\bar{s} = 20$ and $u/\bar{s} = 0.005$ for all results.
\label{fig:diffdetail}
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe1.svg)

If we compare the expected differentiation at any particular locus in a
heterogeneous polygenic barrier to the single locus prediction for that locus,
we find that depending on $L\bar{s}, N_e\bar{s}$ and $m$, the single locus
prediction may be off by quite a lot.
Assuming selection and dominance coefficients to be independently distributed
according to an Exponential and Uniform distribution respectively (see
@sec:dfe); we find that, when $L\bar{s}$ is not small ($>0.5$, say), for large
parts of parameter space, single locus predictions underestimate the extent of
differentiation by more than 10% on average across loci, and depending on the
extent of drift and the total strength of selection, predictions can be off by
40% (@fig:randdiff).
Considering in more detail how different loci across the barrier behave at
migration-selection balance in the polygenic regime, we find that, as expected,
the error of the single locus prediction is largest for locally beneficial
recessives, where it can be as large as 90% (@fig:diffdetail A).
The error of the single locus model for dominant variants is considerably less.
As expected, the excess differentiation relative to the single locus
predictions increases markedly with $L\bar{s}$.
@Fig:diffdetail further illustrates clearly how many, predominantly partially
recessive, alleles are protected from swamping in a polygenic setting roughly
when $L\bar{s} > 1$, although this will, however, also depend on the variance
of the $s_i$ (cfr. @fig:gammas).

### The realized genetic architecture of local adaptation at migration-selection balance

This does not, however, imply that recessives necessarily contribute more to
local adaptation at migration-selection balance than dominant alleles do.
Although strongly selected recessives will be associated with strong
differentiation, weakly selected recessive alleles will be much more prone
to swamping than partially dominant ones.
One way to quantify how these two phenomena interact to yield the *realized*
genetic architecture of local adaptation (related to the concept of *adaptive
architecture*, as defined in @barghi2020)
is by considering the conditional probability density for the selection and
dominance coefficient at a locus, given that a divergent allele is observed on
the island for that locus, i.e.
\begin{align}
  f(s_i,h_i|X_i=1) 
  &= \frac{\Pr\{X_i=1|s_i,h_i\}f_{\mathrm{DFE}}(s_i,h_i)}{\Pr\{X_i=1\}} 
  \propto \int_\mathcal{B} \Ex[p_i|s_i,h_i,B]f_{\mathrm{DFE}}(s_i,h_i,B)dB
  \label{eq:msbdfe}
\end{align}
where $f_\text{DFE}$ denotes the joint density of the selection and dominance
coefficient in the $L$-locus barrier, $X_i$ is an indicator random variable (equalling
1 when a locally beneficial allele is observed at locus $i$ and zero
otherwise), $B$ is a shorthand for the selection and dominance coefficients at
the $L-1$ other loci ('$B$' for background), and we integrate over the set of
all possible such backgrounds $\mathcal{B}$.
Note that $f_\text{DFE}$ is equivalent to $f(s_i,h_i|X_i=1)$ in the absence of
migration.
We can characterize this conditional probability density using a Monte
Carlo approach by sampling random $L$-locus genetic architectures from a DFE
model and calculating for each $(s_i,h_i)$ pair in the barrier the expected
beneficial allele frequency $\Ex[p_i|s_i,h_i,B]$ as a weight.
The weighted sample will be distributed according to $f$.
@Fig:diffdetail (B) shows approximations to the marginal distributions
$f(s_i|X_i=1)$ and $f(h_i|X_i=1)$ obtained in this way for the heterogeneous
barrier model assumed in the preceding section.

As expected, we find that as migration rates go up (colors in @fig:diffdetail
B), the distribution of selection coefficients in the barrier at
migration-selection balance shifts towards higher values of $s$, and that this
effect becomes weaker with increasing $L\bar{s}$, which increases the extent by
which small-effect alleles are protected from swamping by LD.
Notably, we observe that recessives contribute *less* to differentiation than
dominants do when migration is sufficiently strong, despite the fact that,
conditional on no swamping, equilibrium frequencies of recessives are most
affected by LD (@fig:diffdetail A). 
This is most strongly observed when $L\bar{s}$ is not large (top row in
@fig:diffdetail).
When $L\bar{s} = 1.5$ for instance, the depression in the conditional density
at $h=1$ becomes very slight even for relatively large migration rates (bottom
row in @fig:diffdetail).
We observe a similar shift in the distribution of dominance coefficients when
$h$ is Beta distributed with mean $2/3$ instead of uniformly on the unit
interval (@fig:dfe1).
It is noteworthy that, despite $s$ and $h$ being independent at each locus in
the barrier, migration-selection balance induces a correlation between $s$
and $h$ in the distribution conditional on observed divergence, with variants
of relatively large effect observed at equilibrium being more likely to act
recessively than variants of small effect (@fig:dfe1j, see also @fig:dfecomp,
top row).
The correlation is negligible for small migration rates, but as the strength of
migration increases so that swamping effects become relevant, the correlation
coefficient can become as large as $0.25$, depending on $L\bar{s}$.

![
The realized genetic architecture of local adaptation for different DFE models
for increasing strength of migration.
Contour plots for the joint density of $h$ and $s$ conditional on observing a
divergent allele on the island (see @eq:msbdfe) are shown for the three DFE
models (\uline{rows}) for increasing rates of migration (\uline{columns}).
Values of $\bar{\Delta}$ in the lower right corner denote the mean expected
differentiation per locus.
We assume $L\bar{s}=0.8, \bar{s}=0.01, N_e\bar{s}=20, u/\bar{s}=0.005$ and
exponentially distributed selection coefficients, and parameterize the DFE
models so that $\Ex[h] = 2/3$, assuming $\alpha=2, \beta=1$ for the independent
model, $a=7.2, b=1.2, \sigma=1$ for the logistic model and $K=50$ for the CK94
model (see @sec:dfe for details on the different DFE models considered here).
The densities are approximated using a Monte Carlo approach, simulating 500
replicate $L$ locus genetic architectures from the assumed DFE model, fitting a
kernel density estimate to the sample so obtained.
](/home/arthur_z/vimwiki/build/img/2023-07-17/dfecomp.svg){#fig:dfecomp}


The simple DFE model assumed above where selection and dominance coefficients
are independent is almost certainly inadequate.
Theoretical and empirical work has indicated that, on the one hand, a
correlation between selective effect and degree of dominance can be expected in
the standing genetic variation that forms the basis for a polygenic selection
response, with large-effect alleles more likely to act recessively
[@caballero1994; @zhang2004; @agrawal2011].
On the other hand, it is well appreciated that during the process of
adaptation, different loci enjoy different probabilities of rising to high
frequencies, with dominant beneficial alleles having higher establishment
probabilities than (partially) recessive ones with the same homozygous effect
(Haldane's sieve; @haldane1927, @turner1981). 
These two aspects interact when adaptation is from standing variation, as
the on average higher initial frequency of partially recessive alleles
increases the fixation probability, whereas its recessivity decreases it
[@orr2001].
To examine how the realized genetic architecture of local adaptation depends on
such assumptions, we consider two alternative, admittedly *ad hoc*, DFE models,
outlined in @sec:dfe.
Both models assume Gamma distributed selection coefficients and incorporate a
positive correlation between $s$ and $h$, so that alleles of large effect tend
to be more recessive (recall once more that $h$ in our case is the dominance
coefficient of the invading allele, so $h=1$ corresponds to recessive local
adaptation).
We keep the average dominance coefficient fixed to $2/3$ for each model.
In contrast with the independent model, we find that for the models that
incorporate such a correlation between $s$ and $h$, recessives are typically
more likely to contribute to the realized differentiation at equilibrium
(@fig:dfecomp, @fig:dfe1, @fig:dfe2, @fig:dfe3).
When we make the opposite assumption that locally beneficial alleles tend to be
dominant (corresponding, for instance, to the case with a strong Haldane's sieve
effect during adaptation), we find a somewhat less dramatic shift in the joint
density as migration rates go up (@fig:dfe4). The distribution also shifts
towards higher selection coefficients, but somewhat less so than in the model
with the opposite correlation.
Swamping of partially recessive alleles of small effect further shifts the
distribution towards smaller values of $h$. 
Clearly, these examples show how correlations between $s$ and $h$ among the
loci under divergent selection can have a rather important influence on the
realized genetic architecture at migration-selection balance (i.e. on which
loci actually contribute to adaptive differentiation), driving up the
relative contribution of recessives in one case but not in the other.


# Discussion

Speciation, in essence, amounts to the buildup and maintenance of linkage
disequilibria [@felsenstein1981], with different sets of alleles maintained in
different subpopulations.
Local adaptation due to heterogeneous selection (as well as drift and mutation)
across different subpopulations contributes to the buildup of LD, whereas
migration and recombination (including hybridization) break up patterns of LD.
Populations are reproductively isolated to the extent that selection
counteracts this breaking up of associations through recombination.
In this paper, we have focused on reproductive isolation through the
maintenance of local adaptation in the face of gene flow, assuming a scenario
of secondary contact.
Of key importance is that, in the presence of selection, the rate at which
patterns of LD are broken up due to migration and recombination is not
independent of the magnitude of LD.
Selection against migrant genotypes eliminates sets of alleles jointly, leading
to stronger reproductive isolation due to divergent selection when local
adaptation is polygenic.
The effects of polygenic selection against migrant genotypes can be quantified
by the gene flow factor $g$, or the associated effective migration rate $m_e =
mg$, which accounts for the (potentially quite strong) selection against
migrant genotypes and the rapid break up of LD among surviving invading alleles
[@sachdeva2022].

We derived an expression for the effective migration rate in a diploid or
haplodiplontic mainland-island model with heterogeneous selective effects
across loci and showed how it can be used together with classical single locus
population genetic theory to yield accurate predictions of equilibrium allele
frequencies on the island.
Our results show how the maintenance of adaptive differentiation in the face of
gene flow depends jointly on the extent of LD, drift, dominance and variation
in selective effects across the set of selected loci, and that details about
the sexual life cycle can readily be incorporated through a set of effective
parameters.
The general succes of the approach indicates two important features of
polygenic migration-selection balance.
Firstly, it suggests that the 'separation of time scales' argument that is at
the root of the approach indeed works, and does so beyond the haploid case with
homogeneous selective effects [@sachdeva2022].
At least in the unlinked (and probably also weakly linked) case, strong
selection against multilocus genotypes occurs only in the first couple of
generations after a migrant arrives, and the long term fate of a migrant allele
is unaffected by LD conditional on having survived these initial generations.
As a consequence, the effects of LD are well described by the usual single
locus dynamics, but with a reduced migration rate.
Secondly, it indicates that our rather crude approximation to the expected
reproductive value of a migrant individual on the island (which assumes HWLE
within the island population, that migrants only cross with residents, and that
in each such cross the proportion of migrant alleles is exactly halved) is an
adequate estimator of the gene flow factor.

The approach enables us to study the model at population genetic equilibrium
using efficient numerical methods, and in particular to examine the
relationship between the genetic architecture of local adaptation and the
ability to maintain adaptive differentiation in the face of gene flow.
Our analyses for homogeneous genetic architectures indicate that when there is
selection the diploid phase, dominance can have a considerable impact on both
the extent of adaptive differentiation and the ability to maintain it.
Single locus theory can be somewhat misleading in this regard, in that it tends
to highlight the relative precariousness of local adaptation due to recessive
alleles, whereas in the multilocus setting (partially) recessive variants may
produce a much stronger barrier to gene flow than dominant variants with the
same homozygous effect. The reason for this is that for recessive local
adaptation, the dominant invading alleles are immediately exposed to selection,
whereas for dominant local adaptation, the recessive invading alleles can
introgress easily as long as differentiation is low. In the extreme case of
maximal differentiation and completely dominant local adaptation, gene flow
will be unimpeded whatever the extent of LD when migration is at the haploid
stage. 
Our results show that the feedbacks between the level of differentiation and
the strength of selection againts migrants (as measured by the gff, which
depends on the extent of differentiation) are strongly affected by the degree
of dominance.
It should be emphasized, however, that all our results assume a mainland-island
model of migration. The effects of dominance may turn out more subtle in more
symmetric models of population subdivision [e.g. @burger2013].

The assumption of a homogeneous architecture with some shared dominance
coefficient across the entire barrier does not, however, appear very realistic.
We hence directed our attention to the more biologically relevant situation
of a heterogeneous genetic architecture, where both dominance and selection
coefficients vary across loci.
The extent of adaptive differentiation that can be maintained in the polygenic
regime was shown to be rather strongly affected by the variance of selective
effects in the barrier, but much less so by variation in dominance
coefficients when the latter are independent of the former.
When migration is not too strong, heterogeneity in selection coefficients is
likely to weaken the barrier to gene flow generated by local adaptation.
relative to a homogeneous barrier with the same average selection intensity
per locus.
We have further shown that correlations between selection and dominance
coefficients can have a rather strong influence on the realized genetic
architecture of local adaptation at migration-selection balance.

Throughout, we have ignored physical linkage of the loci under selection.
Accounting for (tight) linkage using an approach like ours, based on plugging
in a suitable effective migration rate in single locus theory, is likely not
straightforward.
Associations between tightly linked loci will be broken up by recombination at
rates comparable to or slower than their elimination by selection, so that it
renders the separation of time scales argument on which the approach is based
inappropriate.
Having a theory that can take into account tight linkage is however crucial to
understand the transition from so-called divergence hitchhiking, where some
physical linkage is required for divergent selection to generate an appreciable
barrier effect, to genome hitchhiking, where genome-wide LD leads to a
relatively strong barrier everywhere across the genome [@feder2012].

We have assumed that fitness is determined multiplicatively across loci and
that the population is subject to directional selection, considering a history
where a rapid polygenic selection response has driven allele frequencies at $L$
loci up near fixation.
Thus, we effectively assume that the locally adapted population is not yet
close to a fitness optimum, so that we can ignore stabilizing selection.
When there is abundant standing variation, a polygenic selection response may
however only involve subtle changes in allele frequencies [@hayward2022], and
there may be considerable genetic redundancy [@yeaman2015], leading to a
scenario that is quite different from the one assumed in this paper.
The extent of reproductive isolation that can be maintained when these aspects
of polygenic adaptation become important remains unclear and likely requires
different approaches.
More generally, our focus on the maintenance of polygenic local adaptation and
the reproductive isolation it causes provides only half of the picture, as we
have both ignored the initial polygenic response, and the further building up
of divergence in the face of gene flow.
We have considered how a *given* genetic architecture underlying local adaptation
results in observable patterns of adaptive differentiation at equilibrium, but
remain ignorant about just what sort of genetic variation is likely the source
of local adaptation, neither have we considered how the resulting barrier to
gene flow promotes further divergence.
All these are important topics deserving further study if we are to understand
how populations can remain locally adapted when subjected to maladaptive gene
flow and what processes drive speciation.

# Acknowledgements

This work was supported by the European Research Council (ERC) (\hl{some
number/name}, to CF).

# References
<div id="refs"></div>
\clearpage
\setcounter{page}{1}


\clearpage
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thealgorithm}{S\arabic{algorithm}}
\setcounter{figure}{0}
\setcounter{section}{0}
\setcounter{equation}{0}
\setcounter{algorithm}{0}

# Supplementary figures

\clearpage

![Critical swamping thresholds for the multilocus model. 
Equilibria of the multilocus system correspond to the zeros of $f(p) = hq +
(1-2h)q^2 - m_e/s$. 
Examples for $f(p)$ in the case with dominant local adaptation (top row),
additive local adaptation (middle row) and recessive local adaptation (bottom
row) near the critical point.
The stable equilibrium is indicated by a filled dot, the unstable by an
unfilled dot.
When there is bistability, i.e. both a stable and unstable equilibrium, the
critical migration rate at which the two equilibria collide and cease to exist
corresponds to the value of $m$ for which $f(p)$ reaches its maximum in the
critical point, so that both $f(p) = 0$ and $f'(p)=0$ are satisfied.
](/home/arthur_z/vimwiki/build/img/2023-05-02/stab.svg){#fig:mlstab}

![
Critical equilibrium differentiation ($p_c$, the frequency of the locally
beneficial allele on the island just before swamping) and critical migration rate
($m_c$) for intermediate dominance ($0\le h \le 1$) and low to appreciable
divergence ($Ls \le 1.5$). The solid white line marks the region of parameter
space where the system exhibits bistability (i.e. a sharp swamping threshold at
a critical differentiation level $p_c > 0$). The dashed lines mark $h=1/3$ and
$h=2/3$. For $h>2/3$, bistability occurs for all $m$. For $0 < h < 1/3$, the
minimum $Ls$ for which bistable behavior is observed increases, with increasing
$h$, after which it quickly falls.
](/home/arthur_z/vimwiki/build/img/2023-05-02/phase2.svg){#fig:pplot}

![
Equilibrium differentiation and swamping thresholds for the deterministic
multilocus model, comparing different degrees of dominance on the basis of
$Lsh$. The dashed line shows results for $h=1$ (recessive local adaptation),
whereas the solid line shows results for $h=0.5$.
](/home/arthur_z/vimwiki/build/img/2023-06-28/Lsh.svg){#fig:lsh}


![Comparison of the multilocus diffusion approximation (gray line) against
individual-based simulations (black dots). $Ls =1$ and $N_es = 5$ for all
plots, while $L$ is varied across columns ($L\in [5,10,25,50]$) and $h$ varies
over rows ($h\in[0,0.5,1]$).  We assumed $k=5$ diploids per haploid individual
and set $N = N_e/2k + N_e$ so that the desired $N_e = N_es/(Ls/L)$ is obtained.
Allele frequencies for the individual-based simulations are obtained by
simulating for 110000 generations, sampling every 5 generations after
discarding the first 60000, and averaging across loci. For each $L$ we simulate
$n$ replicates so that $nL = 50$. The mutation rate was set to $u=0.005s$. 
](/home/arthur_z/vimwiki/build/img/2023-05-03/Ls-dominance.svg){#fig:Lsdom}


![
Predicted allele frequency distributions for a single locus in the diploid
multilocus model with homogeneous selective effects for different values of
$N_e$ and $h$ (dominance coefficient of the invading alleles).
Lines show the numerical approximations based on the diffusion theory, dots
show results from individual-based simulations, based on taking a sample every
10 generations for 50000 generations after an initial 10000 generations to
reach equilibrium.
Other parameter settings are $Ls = 0.8, L=40, m/s=0.2, u=0.005s$.
](/home/arthur_z/vimwiki/build/img/2023-05-03/domdriftdist.svg){#fig:driftdist}

![
Effect of drift on equilibrium differentiation and swamping thresholds for
a range of dominance values $h$, ranging from overdominant local adaptation
$h=-1$ (hybrids have an advantage), to underdominant local adaptation $h=2$
(hybrids perform worse than mainland individuals on the island). All results
use $L=40, Ls=0.8, k=5$.
](/home/arthur_z/vimwiki/build/img/2023-05-03/domdrift.svg){#fig:drift}

![Effect of genetic drift, dominance and total barrier strength on equilibrium
adaptive differentiation for a homogeneous polygenic barrier in a diploid
population.
As in @fig:drifthm, but highlighting the $Ls$ by $N_es$ interaction for several
values of $m/s$.
](/home/arthur_z/vimwiki/build/img/2023-06-17/drift-homo-ep.svg){#fig:drifthm1}

![Effect of genetic drift, dominance and total barrier strength on equilibrium
adaptive differentiation for a homogeneous polygenic barrier in a diploid
population.
As in @fig:drifthm, but highlighting the $Ls$ by $m/s$ interaction for several
values of $N_es$.
](/home/arthur_z/vimwiki/build/img/2023-06-17/drift-homo-ep2.svg){#fig:drifthm3}

![
Effective parameters accurately describe equilibrium dynamics for
haplodiplontic populations when selection is sufficiently weak.
The line shows the numerical prediction of the locally beneficial allele
frequency on the island for increasing strength of migration relative to
*effective* selection. The dots show results from individual based simulations
with different degrees of haploid vs. diploid selection and different relative
sizes of the haploid and diploid population, keeping $N_e, s_e$ and $h_e$
however constant.
Simulation results are based on 110000 generations, where we sampled every 10th
generation after discarding the first 10000 generations.
](/home/arthur_z/vimwiki/build/img/2023-06-17/hapdipeff.svg){#fig:hapdipeff}

![
Migration in the diploid phase of a haplodiplontic life cycle with selection in
both phases leads to stronger barriers to gene flow.
The lines show predictions from the multilocus diffusion theory, whereas the
dots show results from individual-based simulations (taking a sample every
fifth generation during 20000 generations after discarding the first 5000
generations). The migration rates in the haploid an diploid stage are $m_1$ and
$m_2$ respectively.
](/home/arthur_z/vimwiki/build/img/2023-06-30/dipmigk2.svg){#fig:dipmig}

![
Predicted equilibrium allele frequencies for increasing migration rates (left)
and frequency distributions (right) for six loci in a $L$-locus multilocus
barrier in a diploid system, where $s \sim \Exp(\bar{s}=0.02)$ and $h \sim
\Beta(1,1)$. Lines show predictions from the multilocus diffusion
approximation, whereas dots show results from individual-based simulations
(simulating for 200000 generations after an initial 10000, sampling every 10th
generation). The frequency distributions are shown for $m/\bar{s} = 0.2$.
\label{fig:het}
](/home/arthur_z/vimwiki/build/img/2023-05-04/exhet2.svg){width=80%}

![
As in @fig:gammas, but now keeping the selection coefficient fixed at $\bar{s}
= 0.01$ and using randomly sampled dominance coefficients, from a symmetric
Beta distribution with parameter $\alpha$. Again, $L=100, N_es = 10,
u/s=0.005$. In (C) we assumed $m/s = 0.3$.
](/home/arthur_z/vimwiki/build/img/2023-06-28/betas.svg){#fig:betas}

![Average difference in predicted allele frequency for the single locus vs.
multilocus model. We assume $L=100$, $s_i \sim \Exp(\bar{s})$ and $h_i \sim
\Beta(1,1)$ for $i=1,\dots,L$. We show $\frac{1}{L}\sum_{i}^L |\Ex[p_{i,L}] -
\Ex[p_{i,\text{single}}]|$ where $p_{i,L}$ and $p_{i,\text{single}}$ are the
equilibrium frequency of the locally beneficial allele at locus $i$ in the
multilocus model and single locus model respectively.
The results are averaged across 10 random $L$-locus barriers. We show results
for different strengths of genetic drift ($N_e\bar{s}$).
Note that values of $L\bar{s}$ range from 0.5 to 2 ($y$-axis). 
](/home/arthur_z/vimwiki/build/img/2023-06-21/driftdiff.svg){#fig:randdiff}

![
As in @fig:diffdetail, but with $h \sim \Beta(2,1)$ (so that $\Ex[h] = 2/3$).
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe1b.svg){#fig:dfe1}

![
As in @fig:diffdetail, but for the logistic regression model (with $\bar{s} = 1,
\kappa=1, a=7.2, b=1.2, \Ex[h] \approx 2/3, \sigma=1$; see @sec:dfe).
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe2b.svg){#fig:dfe2}

![
As in @fig:diffdetail, but for the CK94 model (with $\bar{s}=0.01, \kappa=1$
and $\Ex[h] = 2/3$, yielding $K=50$; see @sec:dfe).
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe3b.svg){#fig:dfe3}

![
Monte Carlo approximation to the joint probability distribution of $s$ and $h$
conditional on observing a divergent allele on the island (@eq:msbdfe) for the
DFE model with independent selection and dominance coefficients (see
@fig:diffdetail and @sec:dfe). Deep blue designates regions of low probability
density, bright yellow regions of high probability density. The estimated
correlation $\rho$ between $s$ and $h$ is shown in the upper right corner.
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe1joint.svg){#fig:dfe1j}

![
As in @fig:dfe1j, but with $h \sim \Beta(2,1)$ (see @fig:dfe1 and @sec:dfe).
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe1bjoint.svg){#fig:dfe1bj}

![
As in @fig:dfe1j, but for the logistic model (see @fig:dfe2 and @sec:dfe).
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe2joint.svg){#fig:dfe2j}

![
As in @fig:dfe1j, but for the CK94 model (see @fig:dfe3 and @sec:dfe).
](/home/arthur_z/vimwiki/build/img/2023-06-25/dfe3joint.svg){#fig:dfe3j}

![
As in @fig:dfecomp, but for the CK94$^\ast$ model with a negative correlation
between $s$ and $h$, see @sec:dfe.
](/home/arthur_z/vimwiki/build/img/2023-06-30/dfe4joint.svg){#fig:dfe4}


\clearpage

# Appendix

## Single locus allele frequency dynamics for weak selection \label{sec:app1}

Consider a single locus in a population of organisms with a haplodiplontic life
cycle. We assume a finite number $n$ of alleles exist at the locus. Let $p_i$,
$i \in [0..n-1]$, denote the frequency of the $i$th allele ($\all_i$).
Ignoring mutation and migration for now, the change in allele frequency of an
allele $A_i$ througout the life cycle is assumed to take the form:
\begin{align}
  \underbrace{
  p_i \underset{\text{haploid selection}}{\longrightarrow}
  p_i^\ast \underset{\text{gametogenesis}}{\longrightarrow}
  p_i^\ast }_{\text{haploid (gametophytic) phase}} 
  \underset{\text{syngamy}}{\longrightarrow}
  \underbrace{p_i^\ast \underset{\text{diploid selection}}{\longrightarrow}
  p_i' \underset{\text{meiosis}}
  {\longrightarrow}}_{\text{diploid (sporophytic) phase}} 
  p_i',
\end{align}
where we have assumed that gametogenesis, syngamy and spore formation do not
affect the allele frequencies.
Generally, migration could take place at any stage in the life cycle, for
instance right after meiosis (e.g. dispersal of meiospores in bryophytes and
Fungi), at the end of the haploid phase (e.g. gamete dispersal in algae), early
in the diploid phase (e.g. seed dispersal in spermatophytes) or at the level of
adult diploids (e.g. migration in animals).

Let the relative haploid fitness of a haploid individual carrying allele $i$ be
$1 + \eps s_i$, defined as the relative contribution to the diploid
(sporophytic) generation within the haploid (gametophytic) generation.
Similarly, we let $1+\eps s_{ij}$ denote the relative fitness of a diploid
individual with genotype $\all_i \all_j$.
The allele frequency change over a single generation is determined by the
following dynamical system:
\begin{align}
  p_i^\ast &= \frac{w_{g,i}}{\bar{w}_{g}(p)} p_i \nonumber \\
  p_i' &= \frac{w_{s,i}(p^\ast)}{\bar{w}_s(p^\ast)}p_i^\ast \qquad 1 \le i \le n,
  \label{eq:ds}
\end{align}
where $p = (p_1, p_2, \dots, p_n)$ (and similarly for $p^\ast$). The marginal
fitnesses $w_{g,i}$ and $w_{s,i}$ associated with allele $\all_i$ in the
gametophytic and sporophytic phase respectively are
\begin{align*}
  w_{g,i} &= 1 + \eps s_i \\
  w_{s,i}(p) &= \sum_j (1+\eps s_{ij})p_j = 1 + \eps \sum_j s_{ij}p_j := 1 +
  \eps \bar{s}_{s,i}.
\end{align*}
The mean fitnesses in the gametophytic and sporophytic phases are
\begin{align*}
  \bar{w}_g(p) &= \sum_i (1+\eps s_i)p_i := 1 + \eps \bar{s}_g \\
  \bar{w}_s(p) &= \sum_i \sum_j (1+\eps s_{ij})p_ip_j := 1 + \eps \bar{s_s}
\end{align*}
The allele frequency change over a single alternation of generations for the
dynamical system defined in @eq:ds has the usual form
  \begin{equation}
  \Delta p_i = \frac{\bar{w}_i - \bar{w}}{\bar{w}} p_i
  \label{eq:change}
  \end{equation}
Where, from @eq:ds, we have 
\begin{equation}
  \bar{w} = \Big(1 + \epsilon \sum_j \sum_k
  \frac{(1+\epsilon s_j)p_j (1+\epsilon s_k)p_k}{(1+\eps \bar{s}_g)^2} s_{jk}\Big) 
  (1+\eps \bar{s}_g) = 
  1 + \eps \bar{s}_g + \frac{\eps \bar{s}_s}{(1+\epsilon \bar{s}_g)^2} +
    O(\epsilon^2)
  \label{eq:wm}
\end{equation}
and
  $$
  \bar{w_i} = \Big(1 + \epsilon
    \frac{\sum_j s_{ij}(1+\eps s_j)p_j}{1+\eps \sum_k s_k p_k}\Big)
    (1+\eps s_i)
  = 1 + \eps s_i
  + \frac{\eps \bar{s}_{s,i}}{(1 + \eps \bar{s}_g)}+
    O(\epsilon^2),
  $$
so that
  \begin{equation}
  \bar{w_i} - \bar{w}
  = \epsilon (s_i - \bar{s}_g)
  + \frac{\epsilon}{(1 + \epsilon \bar{s}_g)^2} (\bar{s}_{s,i} - \bar{s}_s)
  + O(\epsilon^2).
  \label{eq:diff}
  \end{equation}
Assuming the intensity of selection per generation is weak (all $s$ are small)
and that $\eps$ measures the generation time, we obtain a continuous-time model
of allele frequency change by considering the per-generation change in allele
frequency and taking the limit as $\eps$ goes to zero, specifically, plugging
@eq:diff and @eq:wm in @eq:change, we get
\begin{align}
  \dot{p_i} = \lim_{\eps \rightarrow 0}
     \frac{\Delta p_i}{\epsilon} &= 
  \big[(s_i - \bar{s}_g) + (\bar{s}_{s,i} - \bar{s}_s)\big]p_i \nonumber \\
  &= (\bar{s}_i - \bar{s})p_i,
\end{align}
where
\begin{align}
    \bar{s}_i &= s_i + \bar{s}_{s,i} = s_i + \sum_{j}s_{ij}p_j \nonumber \\
    \bar{s} &= \bar{s}_g + \bar{s}_{s} = \sum_i s_i p_i + \sum_i \sum_j s_{ij}
    p_i p_j.
    \label{eq:malth}
\end{align}
This has the same form as the classical diploid or haploid continuous-time
model of allele frequency change[^fn] but with marginal and mean Malthusian
fitnesses given by $\bar{s}_i$ and $\bar{s}$ respectively.

As an example, consider the biallelic case with alleles $\all_0$ and $\all_1$,
so that $s_0 = s_{00} = 0$, $p_0 = p$ and $p_1=q$. Note that we shall always
assume $s_{ij} = s_{ji}$.  We have from @eq:malth $\bar{s}_1 = s_1 + s_{01}p +
s_{11}q$ and $\bar{s} = s_1q + 2s_{01}pq + s_{11}q^2$. Some algebra shows that
we can write the ODE for the frequency of the selected allele ($\all_1$) as
  \begin{equation}
  \dot{q} = (\bar{s}_1 - \bar{s})q = pq(r_1 + r_2q)
  \label{eq:biall}
  \end{equation}
where $r_1 = s_1 + s_{01}$ and $r_2 = s_{11} - 2s_{01}$.
As expected, this is the same dynamical law as for the strictly diploid model,
in which case the dynamics of the selected allele are given by @eq:biall but
with $r_1 = s_{01}$.
This enables us to identify a pair of 'effective' selection coefficients,
  \begin{align}
  s_{01}^\ast = s_1 + s_{01} \nonumber \\
  s_{11}^\ast = 2s_1 + s_{11},
  \end{align}
so that, for weak selection, a diploid biallelic model with parameters
$s_{01}^\ast$ and $s_{11}^\ast$ yields the same allele frequency dynamics[^fit]
as a haplodiplontic model with parameters $s_1, s_{01}$ and $s_{11}$.

[^fn]: The same result can be obtained in a slightly less cumbersome manner by
first noting that, under the assumption of weak selection, allele frequency
changes within a single alternation of generations are negligible, so that
$w_{s,i}(p^\ast) = w_{s,i}(p)$ and $\bar{w}_s(p^\ast) = \bar{w}_s(p)$ in
@eq:ds.

[^fit]: Note that a diploid model with these effective parameters does *not*
yield the same mean fitness (and hence genetic load) as a haplodiplontic model
in the original parameterization if we define mean fitness in the continuous
time model as $\sum_i e^{s_i} p_i \big(\sum_j e^{s_{ij}} p_j\big)$.


## Equilibrium structure of the mainland-island model \label{sec:mieq}

We describe the equilibrium structure of the haplodiplontic single-locus
deterministic mainland-island model for the biallelic case. The dynamics are
given by the ODE
  \begin{align}
    \dot{q} = - \dot{p} 
    &= m\Delta q + pq(s_a + s_b q) \label{eq:ode1} \\
    &= m\Delta q + pq(s_1 + s_{01} + (s_{11} - 2s_{01}) q),
  \end{align}
where $q$ is the frequency of the locally selected allele $\all_1$, and $p=1-q$
is the frequency of the allele with relative fitness of 1 on the island when
homozygous.

When $m=0$ (no migration), there will be an admissible fixed point when either
of the following conditions holds
\begin{align}
  s_{01} &> -s_1 \text{ and } s_{01} > s_1 + s_{11} \\
  s_{01} &< -s_1 \text{ and } s_{01} < s_1 + s_{11} ,
\end{align}
i.e. when there is *ploidally antagonistic selection*, diploid over- or
underdominance, or both.
The fixed point is obtained at
  \begin{equation}
  \tilde{p} = \frac{s_a + s_b}{s_a} = \frac{s_1 + s_{11} - s_{01}}{s_{11} - 2s_{01}}
  \end{equation}
This will correspond to a stable polymorphism whenever $s_{01} > 0$. This case
was first analyzed in a discrete-time model by @scudo1967.

Now consider $m>0$. We shall assume that $\Delta q = q_{\text{mainland}} - q =
1-q = p$, i.e. the mainland is fixed for the locally selected allele.  To
describe the equilibrium behavior, it is helpful to factor the dynamical law as
  \begin{equation}
  \dot{q} = mp\left(1+ \frac{s_a}{m}q + \frac{s_b}{m}q^2\right)
  \label{eq:quad}
  \end{equation}
Linear stability at a fixed point $\tilde{q}$ is determined by
  \begin{align}
  \frac{d\dot{q}}{dq}\big|_{\tilde{q}} 
    &= (s_a - m) + 2(s_b-s_a)\tilde{q} -3s_b \tilde{q}^2
  \end{align}
If $s_b = 0$, we have an effectively haploid model (i.e. *genic selection*),
and will have a stable polymorphic equilibrium at $\tilde{q} = -m/s_a$ whenever
$m < -s_a$, and a stable boundary equilibrium at $\tilde{q} = 1$ when $m >
-s_a$.
When $s_b \ne 0$, polymorphic equilibria, when they exist, will correspond to
the roots of the quadratic expression in parentheses in @eq:quad. These are
  $$q_-, q_+ = \frac{-s_a/m \pm \sqrt{(s_a/m)^2 - 4s_b/m}}{2s_b/m}$$
We have the following biologically relevant equilibria:
       
![Equilibrium and stability behavior of the single-locus biallelic
haplodiplontic mainland-island model. The dark gray zone indicates the
parameter region where there is a single protected stable polymorphic
equilibrium. The light gray zone shows the parameter region where there is both
a stable (but unprotected) and an unstable polymorphic equilibrium.
$s_1$, $s_{01}$ and $s_{11}$ are the haploid and diploid selection coefficients
for the invading allele ($\all_1$) on the island.
\label{fig:stab}](/home/arthur_z/vimwiki/build/notes/img/quadstab.pdf){width=40%}

i. $\tilde{q}=1$ (*swamping*) is always a stable equilibrium when $m > -(s_a + s_b)$.
ii. When $0 < m < -(s_a + s_b)$ there is always a single stable polymorphic
   equilibrium at $q_-$ (dark gray zone in @fig:stab) and $q_+$ will not lie in
   $[0,1]$.
iii. When $-(s_a + s_b) < m < s_b$ and $4s_b/m < (s_a/m)^2 \iff 4m < s_a^2/s_b$, there is,
   besides the stable boundary equilibrium at $\tilde{q}=1$, an unstable
   (repelling) equilibrium at $q_+$, and a stable polymorphic equilibrium at
   $q_-$ (light gray zone in @fig:stab).

The relation between the key parameters $s_a/m$ and $s_b/m$ and the equilibrium
behavior of the system when $m>0$ is illustrated in @fig:stab.

When condition (iii) holds, sharp thresholds for swamping are observed, in
which case there is a certain critical allele frequency $p_c$ below which no
local adaptation cannot be maintained whatever the migration rate.
We can ask for which degree of dominance such sharp thresholds for swamping can
possibly be observed. 
From @eq:quad we see that at an equilibrium which does not correspond to $p=0$,
the condition $f(q) = m + s_a q + s_bq^2 = 0$ holds.
A sufficient condition for observing a sharp threshold is that $f$ obtains a
maximum for some $q < 1$, hence that $f'(1) > 0$ where $f'(q) = s_a + 2s_bq$.
This will be the case whenever $s_a + 2s_b > 0$.
In the diploid case with the usual parameterization where $s_a = -sh$ and $s_b
= -s(1-2h)$, this shows that critical behavior is expected as soon as $h >
2/3$.


## Fixed point iteration algorithm {#sec:fp}

Our approximations yield a system of equations for the expected allele
frequencies and heterozygosities on the island which are coupled through the
gff.
To calculate expected allele frequencies and allele frequency distributions at
equilibrium, we solve the system self-consistently by performing a fixed point
iteration.
In words: for a given initial set of allele frequencies and heterozygosities,
we calculate the gff at each locus using @eq:gff2; using these gff values, we
next calculate expected allele frequencies and heterozygosities at each locus
using numerical quadrature.
This process is repeated until convergence.
The algorithm is more formally outlined in algorithm \autoref{alg:fp}.

\begin{algorithm}[t]                                                                                      
\caption{Fixed point iteration for calculating the expected allele frequency
and expected heterozygosity on the island.}\label{alg:fp}
\begin{algorithmic}[1]
\Require Initialization $p^{(0)} = (p_1^{(0)}, \dots, p_L^{(0)})$, tolerance $\epsilon$
\State $(pq)^{(0)} \leftarrow (p_1^{(0)}q_1^{(0)}, \dots, p_L^{(0)}q_L^{(0)})$
\State $n \leftarrow 1, \Delta \leftarrow \infty$
\While {$\Delta > \epsilon$}
\For{$j=1,\dots,L$}
\State $m_{e,j}^{(n)} \leftarrow \exp\left[ \sum_{i\ne j} s_{a,i} q_i^{(n-1)} + s_{b,i} (p_iq_i)^{(n-1)} \right]$
\State $p_j^{(n)} \leftarrow \int_0^1 p \phi(p; N_e, u, m_{e,j}^{(n)},s_j) dp$
\State $(p_jq_j)^{(n)} \leftarrow \int_0^1 p(1-p) \phi(p; N_e, u,
m_{e,j}^{(n)},s_j) dp$
\EndFor
\State $\Delta \leftarrow \sum_j (p_j^{(n)} - p_j^{(n-1)})^2$
\State $n \leftarrow n+1$
\EndWhile
\State \Return $p^{(n)}, (pq)^{(n)}$
\end{algorithmic}
\end{algorithm}


## Swamping thresholds for the deterministic multilocus model {#sec:supdet}

We now take a closer look at the equilibria of @eq:odeq and their critical
behavior.
Clearly, $p=0$ is always a solution of @eq:odeq, and it will correspond to a locally stable
equilibrium (i.e. swamping) whenever $m/s > 1-h$.
Other equilibria, when they exist, are given by the zeros of the function
\begin{equation}
  f(p) = hq + (1-2h)q^2 - \frac{m}{s}g[p] \label{eq:eq1}
\end{equation}
Note that $g[p] > 0$ and we assume $s > 0$, so that for any fixed $h$, as $m$
increases, there will indeed be a critical migration rate beyond which $f(p) <
0$, from which point onwards the only stable equilibrium will be $p=0$.
At the critical point, the equilibrium allele frequency will satisfy the
additional constraint $f'(p) = 0$ (see @fig:mlstab)[^crit], 
  \begin{align}
  f'(p) &= h + 2(1-2h)q - 2Lm(1 - 3h - 2(1-2h)q)g[p] = 0 
  \label{eq:eq2}
  \end{align}
We can solve @eq:eq1 for $g[p]$, and then plug in $g[p]$ in @eq:eq2.
This yields a cubic polynomial in $p$ which can be solved for the allele
frequency $p_c$ at the critical point:
\begin{equation}
    0=Ls ((1-h)^{2}- 2p^{3} (1-2h)^2 + p^{2} (14 h^{2} - 17 h + 5) - p (7 h^{2} -
    11 h + 4)) -1 + \frac{3}{2} h + (1 - 2 h)p 
\end{equation}
We can then plug $p_c$ into @eq:eq1 and solve for $m_c/s$.
While the general expressions yield not much insight, we can focus on a number
of special cases.

[^crit]: This appears to hold for arbitrary $h$, at least I've convinced myself
of this using some graphs. Haven't proved that it does though.

Firstly, in the additive case ($h=0.5$), a critical point different from $1-h$
appears when $Ls > 1$, in which case the equilibrium frequency at the critical
point will be $p_c = 1-1/Ls$. The corresponding critical migration rate is
  $$\frac{m_c}{s} = \frac{e^{Ls-1}}{2Ls}.$$
In the case where local adaptation is due to dominant alleles ($h=0$), we have
again critical behavior as soon as $Ls>1$, with the swamping threshold
occurring at $m/s=1$ otherwise.
In this case, we find
\begin{align*}
    p_c &= \frac{3}{4} - \frac{\sqrt{L s(Ls + 8)}}{4 L s} < \frac{1}{2}, 
  & \frac{m_c}{s} = 
    \left(\frac{1}{4} + \frac{\sqrt{L s \left(L s + 8\right)}}{4Ls}\right)^{2} e^{\frac{L s}{4}
    + \frac{\sqrt{L s \left(L s + 8\right)}}{4} - 1}.
\end{align*}
In contrast with the additive case (where as $Ls$ increases, arbitrary
equilibrium differentiation can be maintained near the critical point),
equilibrium differentiation will be below $0.5$ near $m_c$ when $h=0$.
Lastly, for recessive local adaptation ($h=1$), we have bistable critical
behavior for all $Ls > 0$. The equilibrium frequency at the critical point is
always larger than $1/2$ and is given by the zeros of the cubic polynomial
  $$4 Ls p^{3} - 4 Ls p^{2} + 2 p - 1 = 0$$
for which we have no simple expressions. A fair approximation for $Ls < 1.5$ is
given by
\begin{align*}
  p_c &\approx \frac{1}{2} + \frac{Ls}{4}, 
  & \frac{m_c}{s} \approx \left(\frac{1}{4} - \frac{(L s)^2}{16}  \right) 
  e^{\left(\frac{Ls}{2}\right)^{3} +  \left(\frac{L s}{\sqrt 2}\right)^2 +
    \frac{Ls}{2}}
\end{align*}
The swamping threshold is seen to increase strongly with increasing $Ls$.

## Sensitivity to initial conditions {#sec:init}

![
Different apparent equilibria depending on initial conditions.
In the left plot, the yellow line ($p_0 = 1$) indicates the expected allele
frequencies as determined using the fixed point iteration of algorithm
\autoref{alg:fp}, starting with $p^{(0)} = (1,1,\dots,1)$ (i.e. secondary
contact, maximal initial differentiation), whereas the green line assumes
$p^{(0)} = (0,0,\dots,0)$ (no initial differentiation). The dots show
results from individual-based simulations with the same initial conditions
(50000 generation, keeping the last 25000 and subsampling every 5 generations).
The three plots on the right show the evolution of the fixed point iteration
for different initial initial conditions $p^{(0)} = (p_0, p_0, \dots, p_0)$ for
three values of $m$.  Note the bifurcation of the dynamical system defined by
the algorithm: for $m/s = 0.2$ and $m/s=0.4$ there is a single globally stable
fixed point, whereas for $m/s = 0.3$, there are two locally stable fixed
points.
\label{fig:bifurcation}
](/home/arthur_z/vimwiki/build/img/2023-04-19/example-bif-fp.svg){width=70%}

Although the equilibrium allele frequency distribution should be independent of
the initial condition (the individual-based model can be thought of as an
ergodic Markov chain on the space of $N$ $L$-locus genotypes), for appreciable
$Ls$, the observed allele frequency distribution in any finite-time simulation
can depend strongly on the initial conditions (that is, as $Ls$ increases,
stochastic jumps between the different modes of the $L$-dimensional joint
allele frequency distribution become increasingly less likely, and occur on
time scales that are neither biologically relevant nor computationally
feasible).
This is especially true in the strongly recessive case ($h > 2/3$) and when LD
is substantial.
This is similar to the behavior in the deterministic model: when $Ls$ or $h$ is
sufficiently large, and sharp swamping thresholds appear, the system is
bistable, and the polymorphic equilibrium cannot be reached when the initial
condition corresponds to a state of no or little differentiation.

The fixed point iteration will in that case converge to an expectation computed
near one of the modes of the allele frequency distribution (see
@fig:bifurcation for an illustration).
In other words, considering the fixed point iteration outlined in algorithm
\autoref{alg:fp} as a discrete dynamical system, and treating $m/s$ as a
bifurcation parameter, two bifurcation points occur succesively, as shown in
@fig:bifurcation. For small $m/s$, a single globally stable polymorphic
equilibrium is obtained. After the first bifurcation point, this equilibrium
ceases to be globally stable, and a second locally stable equilibrium
corresponding to almost no differentiation appears. After the second
bifurcation point the lower equilibrium becomes globally stable.
The region of parameter space where the two stable equilibria coexist
corresponds to the situation where the assumption of population genetic
equilibrium becomes questionable, where the state of the population after a
large but finite time of evolution depends strongly on the detailed history of
the population. 


## Distribution of fitness effects (DFE) models \label{sec:dfe}

![
Joint distributions for an example of each of the three DFE models outlined in
@sec:dfe. The joint probability density is shown on a logarithmic scale, with
yellow marking high density and (deep) blue low density.
The marginal density for the selection coefficient is a Gamma distribution with
$\kappa = 2/3, 1, 3/2$ in the top, middle and bottom row respectively (see
@sec:dfe for the relevant definitions).
(A) Independent selection and dominance coefficients. (B) The logistic model
with $a=9.2$ and $b=2$. (C) The CK94 model with $\bar{h} = 1/3$. \label{fig:dfejoint}
](/home/arthur_z/vimwiki/build/img/2023-06-22/dfejoint.svg){width=70%}

![
Marginal distributions of $s$ and $h$ for the three example DFE model shown in
@fig:dfejoint. \label{fig:dfemarg}
](/home/arthur_z/vimwiki/build/img/2023-06-22/dfemarg.svg){width=70%}

### Independent selection and dominance coefficients

For the independent model, we assume, for $i=1,\dots,L$,
\begin{align*}
    s_i &\sim \Gam(\kappa, \lambda) \\
    h_i &\sim \Beta(\alpha, \beta),
\end{align*}
where $\kappa$ is the shape parameter of the Gamma distribution, and $\lambda$
the rate parameter (i.e. the Gamma distribution with density $f(s) =
\Gamma(\kappa)^{-1}\lambda^{\kappa} s^{\kappa-1}e^{-\lambda s}$).  The mean is
$\kappa/\lambda$, and smaller values of $\kappa$ yield a more leptokurtic
distribution. For $\kappa = 1$, this reduces to the Exponential distribution
with rate $\lambda$. See @fig:dfejoint (A) and @fig:dfemarg for examples.

### Logistic model

In the logistic model, we assume, for $i=1,\dots,L$,
\begin{align*}
    s_i &\sim \Gam(\kappa, \lambda) \\
    \logit h_i |s_i &\sim \text{Normal}(a + b\log s_i,
    \sigma^2),
\end{align*}
where $\logit h = \log\frac{h}{1-h}$ is the logit transform.
In other words, we assume $h$ to be distributed according to a linear
regression on $\log s$ with slope $a$ and intercept $b$, on a logit scale.
The marginal density of $h$ is then
  $$f(h) = \frac{1}{h(1-h)}\int_0^\infty \mathrm{N}[\logit h ; a +
     b\log(s),\sigma] \mathrm{G}(s; \kappa, \lambda)ds\ ,$$
where $\mathrm{N}(\cdot; \mu, \sigma)$ and $\mathrm{G}(\cdot;\kappa,\lambda)$
denote the density functions for the Normal distribution (with mean $\mu$ and
standard deviation $\sigma$) and Gamma distribution respectively.
Instead of setting the $a$ and $b$ parameters directly, we parameterize the
regression by choosing two reference pointst, $s^{(1)}$ and $s^{(2)}$, together
with their respective expected dominance coefficients $h^{(1)} =
\Ex[h|s^{(1)}]$ and $h^{(2)}= \Ex[h|s^{(2)}]$ using
\begin{align*}
    b &= \frac{\logit h^{(2)} - \logit h^{(1)}}{\log s^{(2)} - \log s^{(1)}} \\
    a &= \logit h^{(1)} - b\log s^{(1)}
\end{align*}
This model is also illustrated in @fig:dfejoint (B) and @fig:dfemarg.

### Model after @caballero1994 

In the model of @caballero1994 (see also @zhang2004 and discussion in
@agrawal2011), referred to as CK94, we assume
\begin{align*}
    s_i &\sim \Gam(\kappa, \lambda) \\
    h_i^\ast | s_i &\sim \text{Uniform}(0,e^{-Ks}) \\
    h_i &= 1 - h_i^\ast
\end{align*}
It should be noted that this distribution is supposed to be a reasonable model
for dominance coefficients of deleterious mutations at mutation-stabilizing
selection equilibrium, where mutations of large effect segregating at
appreciable frequencies tend to be recessive.
In our case, we regard this as the distribution of dominance coefficients of
*mutant* alleles on the *mainland* that constitutes the standing variation
which is the source of locally adaptive alleles during the initial polygenic
response (which we do not explicitly model).
These are hence the dominance coefficients of the locally *beneficial* alleles
on the *island* (assuming dominance coefficients to be constant across
environments and genetic backgrounds).
The $h_i$ as we defined them are however the dominance coefficients of the
invading wild-type alleles from the mainland over the locally beneficial ones,
so that we use $h_i = 1-h_i^\ast$ where the $h_i^\ast$ are distributed
according to the CK94 model.
We set the $K$ parameter so that $\Ex[h^\ast] = \bar{h}$ for some $\bar{h}$, i.e.
  $$K =  \lambda \left((2\bar{h})^{-\frac{1}{\kappa}} - 1 \right).$$
The marginal density for $h^\ast$ is
\begin{align*}
  f(h) &= 
  \int_0^{-\frac{\log h}{K}} \frac{\lambda^\kappa}{\Gamma(\kappa)}
    \lambda^\kappa s^{\kappa-1} e^{-(\lambda - K)s} ds \\
    &= \left(\frac{\lambda}{\lambda - K}\right)^\kappa
  \int_0^{-\frac{\log h}{K}} \frac{(\lambda - K)^\kappa}{\Gamma(\kappa)}
    \lambda^\kappa s^{\kappa-1} e^{-(\lambda - K)s} ds 
    = \left(\frac{\lambda}{\lambda - K}\right)^\kappa
  \int_0^{-\frac{\log h}{K}} \mathrm{G}(s;\kappa,\lambda - K)ds \\
    &= \left(\frac{\lambda}{\lambda - K}\right)^\kappa
       \frac{\gamma\left(\kappa, -(\lambda - K) \frac{\log h}{K}\right)}{\Gamma(\kappa)}
\end{align*}
where $\gamma$ is the lower incomplete gamma function.
This model is also illustrated in @fig:dfejoint (C) and @fig:dfemarg.
To study the effect of a negative correlation between $s$ and $h$ (i.e. where
strongly selected locally beneficial alleles tend to be dominant), we use the
same model, but with $h_i = h_i^\ast$. We refer to this model as CK94$^\ast$
(see @fig:dfe4).


## Two-locus haplodiplontic mainland-island model \label{sec:twolocus}

We here consider the deterministic dynamics of a two-locus haplodiplontic
model with mainland-island migration in continuous-time.
We assume a locus $A$ with alleles $A_0$ and $A_1$ and a linked locus $B$ with
alleles $B_0$ and $B_1$, with recombination between the two loci occurring at
rate $r$.
We assume arbitrary dominance and no epistasis, with relative Malthusian
fitnesses of all possible two-locus genotypes in the two phases given by the
following tables
\begin{equation}
\begin{matrix}
\text{Haploid} & & \\
& B_0 & B_1 \\
A_0 & 0     & \beta_1  \\
A_1 & \alpha_1 & \alpha_1 + \beta_1 
\end{matrix} \qquad \qquad 
\begin{matrix}
\text{Diploid} & & & \\
& B_0B_0 & B_0B_1 & B_1B_1 \\
A_0A_0 & 0     & \beta_{01} & \beta_{11} \\
A_0A_1 & \alpha_{01} & \alpha_{01} + \beta_{01} & \alpha_{01} + \beta_{11} \\
A_0A_1 & \alpha_{11} & \alpha_{11} + \beta_{01} & \alpha_{11} + \beta_{11}
\end{matrix}
\end{equation}
Similar to our notation for the single-locus case, we let $\alpha_a = \alpha_1
+ \alpha_{01}$ and $\alpha_b = \alpha_{11} - 2\alpha_{01}$ (and similarly for
$\beta_a$ and $\beta_b$).
Let $x_{ij}$ and $y_{ij}$ denote the frequency of the $A_iB_j$ haplotype on the
island and mainland respectively. The two-locus dynamics in continuous-time are
given by
\begin{equation}
\dot{x}_{ij} = m(y_{ij} - x_{ij}) + (\omega_{ij} - \bar{\omega})x_{ij} - \eta_{ij} rD 
\label{eq:ode-haplotype}
\end{equation}
where $D = x_{00}x_{11} - x_{01}x_{10}$ is the usual measure of linkage
disequilibrium, $\omega_{ij}$ the marginal fitness of the $A_iB_j$ haplotype,
$\bar{\omega}$ the mean Malthusian fitness and $\eta_{ij} = 1$ when $i=j$ and
$-1$ otherwise.
Defining
\begin{align*}
  Q_A &= \alpha_a + \alpha_bq_A &
  Q_B &= \beta_a + \beta_bq_B \\
  P_A &= \alpha_a + \alpha_bp_A & 
  P_B &= \beta_a + \beta_bp_B,
\end{align*}
one can find that $\bar{\omega} = q_AQ_A + q_BQ_B$ and write the marginal
fitnesses as
\begin{align}
  \omega_{00} - \bar{\omega} &= -q_AQ_A - q_BQ_B \nonumber \\
  \omega_{01} - \bar{\omega} &= \beta_a -q_AQ_A - q_BP_B \nonumber \\
  \omega_{10} - \bar{\omega} &= \alpha_a -q_AP_A - q_BQ_B \nonumber \\
  \omega_{11} - \bar{\omega} &= \alpha_a + \beta_a -q_AP_A - q_BP_B.
  \label{eq:twolocus-marg}
\end{align}
Let $p_A = x_{00} + x_{01}$ be the allele frequency of $A_0$ and 
$p_B = x_{00} + x_{10}$ that of $B_0$, and let $q_A = 1-p_A$ and $q_B=1-p_B$.
We shall assume that the mainland is fixed for haplotype $A_1B_1$ (i.e.
$y_{11}=1$).
Using @eq:twolocus-marg with @eq:ode-haplotype, one can derive the dynamics for
$p_A, p_B$ and $D$:
\begin{align}
  \dot{p}_A &= -mp_A - p_Aq_AQ_A - Q_BD \nonumber \\
  \dot{p}_B &= -mp_B - p_Bq_BQ_B - Q_AD \nonumber \\
  \dot{D} &= m(p_Ap_B - D) + (Q_A(p_A-q_A) + Q_B(p_B - q_B))D - rD.
  \label{eq:twolocus-ode}
\end{align}

We can use @eq:twolocus-ode to derive the effective migration rate at a neutral
locus linked to a barrier locus maintained at migration-selection balance.
Let $A$ be the selected locus, and $B$ the linked neutral locus.
The dynamics of the system are given by @eq:twolocus-ode, where $Q_B = 0$.
Assuming $r$ is sufficiently large (linkage is sufficiently weak) so that $D$
equilibrates much faster than the allele frequencies, we can solve the system
for $D$ at equilibrium to find
  \begin{equation}
  \tilde{D} = \frac{mp_Ap_B}{m+r-Q_A(p_A - q_A)}
  \end{equation}
We can plug this into the ODE for the neutral locus to find
  \begin{equation}
  \dot{p}_B = -m\left(1 + \frac{p_AQ_A}{m+r-Q_A(p_A-q_A)}\right)p_B
  \label{eq:qle}
  \end{equation}
Suggesting that the effective migration rate under the stated assumptions
should be
  \begin{equation}
  m_e = m\left(1 + \frac{p_AQ_A}{m+r-Q_A(p_A-q_A)}\right) 
    \approx m\left(1 + \frac{\alpha_a p_A + \alpha_b p_A q_A}{r}\right),
  \end{equation}
where the approximation holds well for weak linkage, which we assumed when we
derived @eq:qle.




