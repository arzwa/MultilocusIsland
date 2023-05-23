---
fontsize: 11pt
geometry: margin=1.2in
title: Maintenance of polygenic local adaptation and reproductive isolation in diploid and haplodiplontic populations
author: Arthur Zwaenepoel, Himani Sachdeva, Christelle Fraïsse
header-includes: 
  - \DeclareMathSymbol{\shortminus}{\mathbin}{AMSa}{"39}
  - \newcommand{\Ex}{\mathbb{E}}
  - \newcommand{\HW}{\mathrm{HW}}
  - \newcommand{\Bin}{\mathrm{Bin}}
  - \newcommand{\Beta}{\mathrm{Beta}}
  - \newcommand{\Exp}{\mathrm{Exponential}}
  - \newcommand{\Bfun}{\mathrm{B}}
  - \newcommand{\eps}{\epsilon}
  - \newcommand{\all}{A}
  - \newcommand{\erf}{\mathrm{erf}}
  - \newcommand{\erfi}{\mathrm{erfi}}
  - \newcommand{\fix}{\mathrm{fix}}
  - \usepackage{caption}
  - \captionsetup[figure]{labelfont=bf,font=small}
  - \usepackage{float,soul}
  - \makeatletter
  - \def\fps@figure{tb} 
  - \makeatother
  - \usepackage{algorithm}
  - \usepackage{algpseudocode}
  - \usepackage{lipsum}
abstract: \lipsum[1]
---

# Introduction

When a population is subdivided across multiple habitats with different 
environmental conditions, the extent to which distinct subpopulations can
maintain locally beneficial genetic variants depends on the rate of migration
between these subpopulations.
Migration between populations that maintain divergently selected alleles can
lead to maladaptive gene flow, yielding a migration load (a reduction in mean
fitness due to the influx of locally maladaptive genes) or may lead to loss
of local adaptation altogether (so-called *swamping* by gene flow).
While local adaptation may be due to a few conspicuous loci (e.g. industrial
melanism in *Biston*; @hof2016), it is believed to typically be polygenic, with
alleles of different effects responding to selection in the local environment
at many loci spread across the genome [@bomblies2022; @stankowski2022;
@westram2018].
When local adaptation is polygenic, migration from a population adapted to
different environmental conditions will lead to linkage disequilibria (LD)
among selected alleles (i.e. migrant alleles will tend to be found together in
the genome), and the rate at which the invading locally deleterious alleles are
eliminated by selection will be affected by such associations [@feder2012;
@yeaman2015; @sachdeva2022].
This in turn will affect the equilibrium migration load and swamping thresholds
for the loci under selection.
Neutral variation may also come to be associated with locally selected alleles,
so that the latter constitute a 'barrier' to neutral gene flow, increasing
neutral genetic differentiation (as quantified by $F_{ST}$ for instance) beyond the
single locus neutral expectation.
Such barrier effects due to divergent local adaptation may play an important
role in the evolution of reproductive isolation, and hence speciation
[@nosil2012]. 

Despite mounting evidence that local adaptation is often polygenic in nature,
little is known about the underlying genetic details: How many loci are
involved? What are the typical effect sizes? Are locally beneficial alleles
typically closely linked (forming so-called 'genomic islands'), or are they
spread all over the genome? How non-additive is local adaptation? *etc.*
[e.g. @yeaman2015; @bomblies2022].
Moreover, even assuming such details to be known, it is unclear to what extent
the genetic architecture of local adaptation affects the ability of a
population to maintain reproductive isolation in the face of gene flow.
Conversely, it remains unclear just to what extent one could hope to infer the
detailed genetic architecture underlying local adaptation from observed
patterns of genomic differentiation, as for instance obtained through so-called
'genome scans'.
So far, most theoretical developments have assumed rather simple genetic
architectures, dealing with biallelic loci of equal additive effect (ignoring
dominance and epistasis) that are either unlinked or uniformly spread along a
block of genome; and statistical approaches for the inference of gene flow
across the genome either make similarly crude assumptions [@aeschbacher2017],
or ignore the genetic details of local adaptation altogether [@fraisse2021;
@laetsch2022].


In a recent paper, @sachdeva2022 showed that, when the loci under selection are
unlinked, the effects of LD on migration-selection balance at any individual
locus in a multilocus barrier can be well described by classical (deterministic
or stochastic) single locus population genetic theory, provided that the
migration rate $m$ is substituted by a suitably defined (selective) *effective*
migration rate $m_e$ [@petry1983; @bengtsson1985; @kobayashi2008], which
captures the effects of selection against the associated genetic background.
The extent to which *neutral* gene flow is reduced by LD with alleles under
selection can be similarly quantified by the (neutral) effective migration rate
[@petry1983; @barton1986; @kobayashi2008], which can serve as a quantitative
measure of reproductive isolation [@westram2022].
In her paper, @sachdeva2022 conducted a detailed study of the effects of both
drift and LD on swamping thresholds and neutral differentiation in the
mailand-island and infinite-island models of population subdivision, assuming
a haploid sexual life cycle and $L$ divergently selected loci of equal effect.
Here we extend some of this work to investigate the effects of dominance and
variation among selective effects on the maintenance of polygenic local
adaptation and reproductive isolation in populations with a diploid or
haplodiplontic life cycle, focusing on the mainland-island model.

All sexual organisms have a life cycle with a distinct haploid and diploid
stage, with remarkable diversity in the relative duration and development of
the two phases.
The existence of a diploid phase entails that interactions between different
alleles at a locus (genetic dominance) will generally affect phenotypic traits
and hence fitness.
The effects of dominance on equilibrium frequencies and swamping thresholds in
the single locus mainland-island model are well known [@haldane1930VI;
@nagylaki1975].
Consider a diploid model (i.e. a model with no selection in the haploid phase),
where the relative fitnesses of genotypes $A_0A_0:A_0A_1:A_1A_1$ are given by
$1:1-hs:1-s$ and where we assume a proportion $m$ of the population is replaced
by migrants of the $A_1A_1$ genotype each generation.
For the case $h=0.5$ (no dominance, also referred to as codominance, or
additivity), the equilibrium frequency $\tilde{p}$ of the locally beneficial
$A_0$ allele decreases linearly from $1$ to $0$ as the rate of migration
approaches the strength of selection on a single allele $\frac{s}{2}$.
When local adaptation is due to a dominant allele (so that the invading allele
acts recessively to reduce fitness on the island, i.e. $h=0$), the migration
rate beyond which no polymorphism can be maintained is increased to $s$, while
$\tilde{p}$ is decreased relative to the additive case as long as $m <
s/4$.
When local adaptation is due to a recessive allele ($h=1$), the model has two
equilibria for the beneficial allele frequency, one stable equilibrium
$\tilde{p}_+ > 1/2$ and one unstable equilibrium $\tilde{p}_- < 1/2$.
The two equilibria collide at the critical migration threshold of $m_c = s/4$,
beyond which swamping occurs for any initial frequency.
Hence, for the recessive case, whether or not a polymorphism is attained
depends not only on the migration rate (which should be at most $s/4$), but
also on the history of the population: the island population cannot fix a new
recessive beneficial variant, but an established recessive variant can be
maintained upon secondary contact.
All these phenomena were first noted by @haldane1930VI.
A main goal of the present paper is to understand how this behavior is affected
when we consider the polygenic case.

A second question we aim to address is to what extent LD among alleles with
different fitness and dominance effects protects weakly selected alleles in the
barrier from swamping and determines variation in the levels of differentiation
maintained across the barrier at equilibrium.
Despite the many recent large-scale genomic studies, the genetic architecture
of local adaptation has remained rather elusive, and it appears relevant to
ask, for a given genetic architecture of local adaptation, what sort of signals
one is likely to observe in empirical data, and to what extent selective
interference between different loci blurs the signatures of individual loci.
[\hl{linkage to a barrier locus? connect to} @yeaman2015?]

We start by outlining a single locus mainland-island model for a haplodiplontic
life cycle which in the weak selection continuous-time limit encompasses both
the haplontic, diplontic and haplodiplontic cases.
We then extend this to the multilocus setting by deriving an approximation to
the effective migration rate based on the expected reproductive value of a
migrant individual in the island population.
We study the thresholds for swamping in the deterministic multilocus model with
$L$ divergently selected loci harboring alleles of equal effect, and find that
the effect of increasing the total strength of selection against migrants
($Ls$) on swamping thresholds and equilibrium differentiation depends rather
strongly on the dominance coefficient.
We then study the effects of drift on the level of differentiation maintained
at equilibrium and consider heterogeneous barrier architectures, where both
dominance and the intensity of selection varies across the loci under
selection.
Throughout, we highlight the remarkable accuracy of the heuristic approach
first outlined in @sachdeva2022, and illustrate how it extends in a rather
straightforward way to multilocus models with dominance, haploid selection and
heterogeneous selective effects.


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
The mainland population is assumed to have a constant, but arbitrary, genetic
composition.
Unless stated otherwise, we shall assume the mainland to be fixed for the
locally deleterious allele on the island.
Fitness on the island is determined by $L$ unlinked biallelic loci which are
under divergent selection relative to the mainland.
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
and the code is available at \hl{GitHub}.
In the following sections, we build up a theoretical approximation to this
fairly general multilocus model, roughly as in @sachdeva2022.
We first derive the dynamics at a single locus, considering both deterministic
and stochastic models.
Next, we derive an approximation to the effective migration rate in the
multilocus system using a rather general argument based on the reproductive
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
where $s_a = s_1 + s_{01}$ and $s_b = s_{11} - 2s_{01}$, the latter being a
measure of dominance (i.e. the deviation from multiplicative fitnesses,
sometimes called $\iota$ [@otto2003; @manna2011]).
Usually, $s_1, s_{01}$ and $s_{11}$ will be assumed to be negative, and $p_M$
will be assumed to be small, so that selection increases $p$, whereas migration
decreases $p$.
As expected, this is the same dynamical law as for a strictly diploid model, in
which case $s_a = s_{01}$.
This enables us to identify a pair of selection coefficients $s_{01}^\ast =
s_1+s_{01}$ and $s_{11}^\ast = 2s_1 + s_{11}$.
When $s_{11}^\ast \ne 0$, we can hence describe the haplodiplontic system as
the familiar diploid model with some effective degree of dominance $h_e =
s_{01}^\ast/s_{11}^\ast = \frac{s_1 + s_{01}}{2s_1 + s_{11}}$ and an effective
selection coefficient $s_e = s_{11}^\ast = 2s_1 + s_{11}$, at least when
selection is sufficiently weak so that allele frequencies do not change
appreciably within any one alternation of generations.
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
This Markov chain model is readily extended to include mutation and migration.
This model is akin to the standard Wright-Fisher (WF) model for a diploid
or haploid population.
Indeed, for the neutral case, it is easily seen that the model corresponds to a
haploid WF model with variable population size, regularly alternating between
$N$ and $2Nk$ gene copies.
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
  V(p) &= N_e^{-1} pq\ .
  \end{align*}
This yields a probability density function for the equilibrium allele frequency
distribution
  \begin{equation}
  \phi(p; N_e, u, m, s) \propto p^{2N_e(u+mp_M)-1}q^{2N_e(u+mq_M)-1}e^{N_e(2s_aq + s_bq^2)},
  \end{equation}
where no closed form expression is known for the normalizing constant.
This is essentially Wright's [@wright1937] distribution for a general
haplodiplontic life cycle.

## Multilocus model {#sec:ml}

### Effective migration rate

To begin constructing a useful approximation to the allele frequency dynamics
in the multilocus system, we derive an appropriate expression for the effective
migration rate $m_e$, which captures the reduction in gene flow at a focal
locus due to selection against migrant genotypes.
As shown formally in @kobayashi2008, for weak migration, the reduction in gene
flow relative to the 'raw' migration rate $m$, termed the *gene flow factor*
(gff), depends on the expected reproductive value of migrants in the resident
background.
At any time, the proportion of individuals with recent migrant ancestry on the
island is $O(m)$, so that the probability of individuals with migrant
backgrounds mating with each other to produce, for instance, F2 crosses of the
migrant and resident genotypes, is $O(m^2)$, and hence negligable for
sufficiently weak migration.
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
In order to derive a useful approximate expression for $g$, we shall make two
further important assumptions:
(1) both the resident and migrant gene pools are are in Hardy-Weinberg and
linkage equilibrium;
(2) variation and selection among offspring of a $k$th generation migrant
descendant $\times$ resident cross (i.e. segregation variance within F1s, BC1s,
*etc.* and selection thereupon) is negligable.
This last assumption becomes more plausible as local adaptation is due to
more and more loci of smaller effect.

Under these assumptions, each of the $W_{\cdot,k}$ is determined solely
by the frequencies of the selected alleles in the mainland and the island
populations at the assumed equilibrium.
This allows us to determine $\Ex[q_{i,k}]$, the expected frequency of the
locally deleterious allele at locus $i$ among $k$th generation descendants from
a migrant, in terms of the allele frequencies in the mainland and island
population.
Indeed, assumption (3) implies the recursive relation $q_{i,k} =
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
  &= \exp\left[2^{-k}\sum_{i=1}^L s_{i01}(q_{M,i} - \Ex[q_i]) - s_{i,b}(p_{M,i}\Ex[q_i] -
  \Ex[p_iq_i])\right],
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
It is worth stressing that the gff is a function of the actual differentiation
between the mainland and island population, and that although we assume
migration is sufficiently rare, we do *not* assume that alleles of the type
introduced by migrants are rare.
We shall often highlight the dependence of the gff on the allele frequencies by
writing $g[p]$. 
If we assume all loci to have the same selection coefficients, and that the
mainland is fixed for the locally deleterious allele on the island, the gff
simplifies to
\begin{equation}
  g = e^{2L(s_a \Ex[p] + s_b \Ex[pq])}
  \label{eq:eqeff}
\end{equation}
where $\Ex[p]$ and $\Ex[pq]$ are the expected beneficial allele frequency and
expected heterozygosity at any selected locus on the island.
Note that in our applications $s_a < 0$, so that when the locally beneficial
allele is common ($p \approx 1$) and $Ls_a$ is appreciable, gene flow will
indeed be reduced due to selection against migrants. 

To gain some intuition for this result, we can express @eq:eqeff in terms of
the effective selection coefficient $s_e$ against the invading allele, and the
effective dominance coefficient $h_e$ of the invading allele over the locally
beneficial one, as
\begin{equation}
  g = e^{-2Ls_eh_e\Ex[p]}e^{-2Ls_e(1-2h_e)\Ex[pq]}
  \label{eq:eqeff}
\end{equation}
(see @sec:sldet). Here, the first factor is just the gff associated with a
haploid $L$-locus system with selection coefficients $s_eh_e$ [@sachdeva2022].
The second factor captures the effects of dominance when heterozygosity is
appreciable.
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
dispersal of gametes, the first generation experiencing selection on the island
will be the diploid F1 generation, so that the appropriate gff under the same
approximation is $g/\Ex[W_{h,0}]$. Secondly, when migration occurs at the
beginning of the diploid phase (e.g. seed dispersal), the first generation
experiencing selection on the island will consist of diploid migrant
individuals, so that $g\Ex[W_{d,0}]$ is the appropriate gff, where
\begin{equation*}
\Ex[W_{d,0}] = \exp\left[\sum_{i=1}^Ls_{i11}(q_M - \Ex[q]) + s_{b,i}(\Ex[pq] - p_Mq_M)\right]
\end{equation*}
In the present work, we shall always assume migration is due to dispersal of
haploid spores, so that @eq:gff gives the relevant gff.

### Dynamics and equilibria for the multilocus model {#sec:dynamics}

The gff captures the effect of LD among selected loci on the rate of gene flow
from the mainland into the island at any individual locus.
The key observation is that a certain separation of time scales applies:
although selection against migrant *genotypes* can be very strong in the
polygenic case (of magnitude $Ls$, roughly), selection at any individual locus
is still assumed to be weak, so that after an evolutionarily short period in
which entire sets of alleles are efficiently removed together, LD among
selected loci quickly becomes negligable and the standard single locus theory
should be appliccable.
Hence, on the longer time scales at which migration-selection balance is
attained, the allele frequency dynamics at any individual locus should
essentially follow the single locus dynamics, but where migrant alleles are
introduced at a reduced rate [@sachdeva2022].

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
dependence of the gff at locus $j$ on the $L-1$ other loci.
We study the equilibria of this model by numerically solving for $p$ at
stationarity ($\dot{p}_j = 0$, for $1 \le j \le L$).

As in @sachdeva2022, we can also plug in $m_e$ in the single locus diffusion
approximation to determine the equilibrium allele frequency distribution for
each locus on the island.
Specifically, we postulate that the joint distribution of allele frequencies in
the barrier factorizes as
\begin{align*}
\phi(p) = Z^{-1} \prod_{j=1}^L \phi_j(p_j|p_{-j}) = 
    Z^{-1}\prod_{j=1}^L \phi(p_j; N_e, u, mg_j[p_{-j}], s_j)
\end{align*}
where $Z$ is a normalizing constant
(this can be thought of as a Markov random field over the complete graph with
$L$ vertices).
We can compute the expected allele frequencies by solving the system
self-consistently, assuming $\Ex[p_j] = Z_j^{-1}\int p_j
\phi_j(p_j|\Ex[p_{-j}]) dp_j$.
Specifically, assuming the mainland to be fixed for the deleterious allele on
the island at all loci, we solve the nonlinear system of $2L$ equations
\begin{align*}
\Ex[p_j] &= Z^{-1}\int p^{2N_eu} q^{2N_e(u + mg_j[\Ex[p_{-j}]]) - 1} \psi_j(p) dp \\
\Ex[p_jq_j] &= Z^{-1}\int p^{2N_eu} q^{2N_e(u + mg_j[\Ex[p_{-j}]])} \psi_j(p) dp,
\end{align*}
where
\begin{align*}
g_j[\Ex[p_{-j}]] &= 
    e^{\sum_{i\ne j}^L s_{a,i} \Ex[q_i] + s_{b,i} \Ex[p_iq_i]} \\
\psi_j(p_j) &= e^{N_e(2s_{j,a}q_j + s_{b,j}q_j^2)},
\end{align*}
for $\Ex[p_j]$ and $\Ex[p_jq_j]$.
To do so, we use the fixed point iteration outlined in @sec:fp.


# Results

## Barrier effect and swamping thresholds in the deterministic model

![Equilibrium differentiation and swamping thresholds for the deterministic
multilocus model. Equilibrium frequencies for increasing $Ls$ are shown for (A)
the case of dominant local adaptation (recessive migrant alleles), (B) additive
fitness effects and (C) recessive local adaptation. The thick lines show the
stable equilibria for increasing $m/s$, whereas the dotted lines show unstable
equilibria. The black dots mark the critical point beyond which swamping is
obtained for any initial condition (in (C) the approximate expression discussed
in the main text is used). The results for $Ls=0.01$ (gray lines) correspond to
the single locus predictions. (D) Swamping thresholds for different degrees of
dominance for increasing total barrier strength $Ls$. 
](/home/arthur_z/vimwiki/build/img/2023-05-02/detdom.svg){#fig:detdom}

We first analyze the deterministic multilocus model for a homogeneous barrier,
where the fitness and dominance effects of all loci are the same.
We shall assume the mainland to be fixed for the locally deleterious allele at
all loci.
We only consider diploid selection here, with $s_{01} = -sh$ and $s_{11}
= -s$, where $s$ is the selection coefficient against the locally deleterious
allele, and $h=s_{01}/s_{11}$ is the dominance coefficient.
Note that $h$ measures dominance of the mainland allele over the island allele,
so that $h=1$ corresponds to a situation where the invading allele is fully
dominant, or, equivalently, where the allele that confers local adaptation on
the island is recessive.
For the case where migration is at the haploid stage, restricting the analysis
to diploid selection incurs no loss of generality, as haploid selection then
simply amounts to a rescaling of the dominance and selection coefficients (see
methods).
The effective (diploid) dominance coefficient when there is haploid selection
with intensity $s_1$ will be $h_e = \frac{s_1+s_{01}}{2s_1+s_{11}}$, so that
the effect of haploid selection is to pull $h_e$ towards the additive case
($h=1/2$).

In the homogeneous deterministic model, all loci have the same dynamics if the
initial allele frequencies are equal.
From @eq:odeme, we find that in that case, the frequency $p$ of the locally
beneficial allele at any individual selected locus in an $L+1$ locus system
must satisfy at equilibrium
\begin{align}
   0 &= shpq + s(1-2h)pq^2 -mg[p]p 
   \label{eq:odeq} \\
   &\text{where }\quad g[p] = e^{-2Ls(hp + (1-2h)pq)}. \nonumber
\end{align}
We can solve this numerically for the equilibrium allele frequency.
@Fig:detdom shows the equilibrium behavior for a number of example parameter
sets.
As expected, increasing the total strength of selection ($Ls$) increases the
equilibrium frequency of the locally beneficial allele with respect to the
single locus prediction, but the magnitude of this effect depends quite
strongly on dominance.
In particular the effect is much more pronounced when local adaptation is due
to recessive variants ($h=1$).
In this case, the gff is at its minimum, and hence the impact of LD on
migration-selection balance is strongest, when the deleterious allele is
rare (@fig:gff).
On the other hand, when the invading alleles are recessive ($h=0$), gene flow
is not at all impeded when the deleterious allele is rare (the gff being near
one).
When $h<1/3$, the barrier strength, as measured by $g^{-1}$ [@barton1986],
*increases* as the deleterious allele increases in frequency on the island (and
hence as differentiation between mainland and island *decreases*),
decreasing the rate of gene flow, until a value of $q=(3h-1)/(4h-2)$ is reached
(@fig:detdom, @fig:gff).
This is essentially because in the latter case, irrespective of how many
deleterious alleles a migrant carries, if the deleterious alleles are rare on
the island they will not be found in homozygotes as long as migration is
sufficiently weak, and hence will not be 'seen' by selection.
On the other hand, when the deleterious alleles are segregating at appreciable
frequencies on the island, F1, BC1, *etc.* individuals will be more likely to 
be homozygous at several loci, thus exposing invading alleles to selection and
reducing the RV of migrants.
As a result, when invading alleles at the selected loci act recessively, a
strong genome-wide barrier effect emerges only once differentiation falls below
a critical threshold. 
The situation is clearly different when migrant alleles are dominant, as the
invading alleles will immediately express their full load in the resident
population, irrespective of the other allele at the locus, yielding efficient
selection against migrant alleles.
When the frequency of the deleterious allele increases on the island (and
differentiation decreases), this will merely increase the expected relative
fitness of migrants in the resident background, and hence reduce the efficiency
of selection against migrant genotypes.

In the single locus model, arbitrary small frequencies of the locally
beneficial allele can be maintained at migration-selection balance when $h <
2/3$, whereas in the case of $h > 2/3$, a sharp swamping threshold is observed
as the equilibrium frequency of the locally beneficial allele reaches some
critical frequency $p_c \le 1/2$.
@sachdeva2022 showed that such sharp thresholds for swamping also appear in
the absence of dominance due to LD.
LD both increases the critical migration rate ($m_c$) at which swamping occurs
and the minimum level of differentiation that can be maintained before the
swamping point is reached ($p_c$).
Our results indicate that dominance has a considerable influence on how LD
sharpens and displaces swamping thresholds (@fig:detdom, @sec:supdet).
As $Ls$ increases, the behavior for the case where $h < 2/3$ (i.e. when local
adaptation is not strongly recessive) is roughly similar to the additive case,
with critical behavior emerging as $Ls$ surpasses some critical value
(@sec:supdet, @fig:pplot).
However, the critical migration rate at which swamping occurs is only
marginally affected by LD in the latter case for moderate levels of divergence
($Ls < 2$, say).
This is in sharp contrast with the case where local adaptation is due to
strongly recessive alleles, where the critical migration increases rapidly with
$Ls$.
Importantly, the critical differentiation ($p_c$) level below which local
adaptation collapses is very different for different degrees of dominance.
In the additive case, one can show that $p_c = 1-1/Ls$ (@sec:supdet), so that
arbitrary differentiation can be obtained at the critical point depending on
$Ls$.
For completely dominant local adaptation, however, $p_c$ increases from $0$ to
$1/2$ as $Ls \rightarrow \infty$, whereas for recessive local adaptation, $p_c$
increases from $1/2$ to $1$ as $Ls$ grows.
This means, in particular, that for moderate levels of divergence, $Ls > 0.75$
say, and large population sizes, one would not expect to see locally beneficial
recessives at frequencies much below $0.8$, compared to $0.5$ for the single
locus model.
\hl{As Himani noted before, it may be worthwhile to highlight that although a
consideration of the gff in the regime where deleterious alleles are rare may
suggest a strong barrier, it may be that swamping thresholds are hardly
affected because of what happens to the gff as differentiation decreases...}

## Accounting for drift and comparison to individual-based simulations

While the deterministic analysis points towards important effects of dominance
on equilibrium differentiation and swamping thresholds, it is important to
assess whether we expect such effects to be as pronounced in finite
populations.
Indeed, @sachdeva2022 showed that, for small populations, the sharp thresholds
for swamping predicted by the deterministic multilocus theory when $Ls$ is
appreciable need not apply, and that the critical migration rate may be
significantly reduced.
Furthermore, for all but the largest populations, the actually observed
frequency of a locally adaptive allele (and hence differentiation between the
mainland and island) at migration-selection-drift balance in the model as
outlined in @sec:model may deviate substantially from the deterministic
approximation, so that it becomes important to understand the *distribution* of
allele frequencies.
More pragmatically, it is hard to evaluate the accuracy of our approximations
against simulations, as any individual-based simulation will necessarily
exhibit the effects of genetic drift.

![
Predicted allele frequency distributions for the diploid multilocus model with
homogeneous selective effects for different values of $N_e$ and $h$ (dominance
coefficient of the invading alleles). Lines show the numerical approximations
based on the diffusion theory, dots show results from individual-based
simulations, based on taking a sample every 10 generations for 50000
generations after an initial 10000 generations to reach equilibrium.
Other parameter settings are $Ls = 0.8, L=40, m/s=0.2, u=0.005s$.
](/home/arthur_z/vimwiki/build/img/2023-05-03/domdriftdist.svg){#fig:driftdist}

Similar to @sachdeva2022, we find that the heuristic multilocus diffusion
approximation obtained by substituting $m_e$ for $m$ in the single-locus
diffusion theory yields remarkably accurate predictions to the multilocus model
outlined in @sec:model.
Indeed, even in parameter regimes where the approximation is expected to break
down ($Ls$ appreciable with $L$ small and $s$ large, small population size) we
obtain good predictions (@fig:Lsdom).
Not only can we reliably obtain the expected frequency of alleles on the
island, we also obtain very good predictions for the entire allele frequency
distribution as observed in individual-based simulations (@fig:driftdist).
The sharp swamping thresholds observed in the deterministic model, in
particular with recessive local adaptation ($h=1$), correspond to strongly
bimodal allele frequency distributions in the model with drift.
As described in more detail in @sec:init, this may render our numerical
approaches sensitive to the initilization of the island population.
Throughout, we shall assume a scenario of secondary contact, so that the island
starts as fixed for the locally beneficial allele at all loci.

The general effects of genetic drift are as expected, in that the barrier
effect is considerably reduced as $N_e$ decreases (@fig:drift).
The characteristic effects of dominance hence only become apparent for
appreciable $N_e$.
For large population sizes, we now see clearly that for dominant local
adaptation (invading alleles are (partially) recessive, $h<0.5$), two phases
can be distinguished as the migration rate increases.
As noted above, the barrier strength increases initially as migration
increases, so that as long as the deleterious allele is below some frequency
determined by $h$, the equilibrium differentiation declines less fast with
increasing $m$. Once a threshold frequency is reached, the barrier strength
declines again as $m$ increases (@fig:drift $h=0$ and $h=0.25$ panels).
\hl{perhaps not worth to spell out.}

![Multilocus migration-selection balance and swamping with dominance and
haploid selection.
For a given total barrier strength (see main text), we vary the relative
strength of selection in the diploid and haploid phase ($\tau$) and the degree
of dominance in the diploid phase.
Specifically, we assume $s_1 = -(1-\tau)s, s_{01} = -h\tau s$ and $s_{11} =
-\tau s$.
Other parameters were as follows: $Ls =0.8, L=40, N_es=8, k=5, u=s/100$.
](/home/arthur_z/vimwiki/build/img/2023-04-19/domtau.svg){#fig:domtau}

We now consider in more detail to what extent our approximations apply to the
case with both selection in the haploid and diploid phase.
In @fig:domtau, we parameterize the model in such a way that we can
investigate, for a given total barrier strength (corresponding to the fitness
difference between an island population fixed for the locally beneficial and an
island fixed for the locally deleterious allele) the effects of the relative
strength of selection in the diploid and the haploid phase and the degree of
dominance in the diploid phase.
To this end, we assume
\begin{align*}
   s_1 = -(1-\tau)s & & s_{01} = -h\tau s & & s_{11} = -\tau s,
\end{align*}
where $0 \le \tau \le 1$ measures the relative strength of selection in the
diploid phase (if one assumes selection to operate continuously with constant
intensity throughout the life cycle, this can be interpreted as the relative
length of the diploid phase).
Similar models have appeared in the study of life cycle modifiers, see for
instance, @otto1994 and @scott2017.
Using this parameterization, we find that when selection is not too strong, we
can indeed accurately predict equilibrium frequencies for haplodiplontic life
cycles where selection occurs in both the haploid and diploid phase.
Furthermore, these results suggest that, for a given total strength of
selection, predominantly haploid populations should be able to maintain more
adaptive variation, and exhibit stronger reproductive isolation, irrespective
of the degree of dominance in the diploid phase.
Although haploid selection pulls the effective dominance coefficient towards
$0.5$ (so that one might expect, based on our results for diploids above, that
a diploid phase with recessive local adaptation would enable more adaptive
differentiation), the effective selection coefficient in this model is
$(2-\tau)s$, so that the strength of selection per gene in haploids is twice
that in diploids in the absence of dominance. 
The relevance of these observations for the evolution and maintenance of
haplodiplontic life cycles is however not very clear, as a life cycle modifier
need not keep the overall strength of selection constant [@scott2017].

## Heterogeneous genetic architectures

![
Predicted equilibrium allele frequencies for increasing migration rates (left)
and frequency distributions (right) for six loci in a $L$-locus multilocus
barrier in a diploid system, where $s \sim \Exp(\bar{s}=0.02)$ and $h \sim
\Beta(1,1)$. Lines show predictions from the multilocus diffusion
approximation, whereas dots show results from individual-based simulations
(simulating for 200000 generations after an initial 10000, sampling every 10th
generation).
\label{fig:het}
](/home/arthur_z/vimwiki/build/img/2023-05-04/exhet2.svg){width=80%}

In @sec:ml, we developed the multilocus theory for potentially heterogenous
unlinked architectures, where the $s_1, s_{01}$ and $s_{11}$ can vary
arbitrarily across loci.
We verify that we obtain accurate predictions also in this setting, using
simulations with randomly sampled selection and dominance coefficients
(@fig:het).


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
The barrier strength $b=g^{-1}$ (top row), and its derivative with respect to
the deleterious allele frequency $q$ (bottom row), for $Ls=0.5, 1, 1.5$
(columns) and different degrees of (effective) dominance ($h$, colors).
](/home/arthur_z/vimwiki/build/img/2023-05-03/barrier.svg){#fig:gff}

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
Effect of drift on equilibrium differentiation and swamping thresholds for
a range of dominance values $h$, ranging from overdominant local adaptation
$h=-1$ (hybrids have an advantage), to underdominant local adaptation $h=2$
(hybrids perform worse than mainland individuals on the island). All results
use $L=40, Ls=0.8, k=5$.
](/home/arthur_z/vimwiki/build/img/2023-05-03/domdrift.svg){#fig:drift}

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
    \approx m\left(1 + \frac{p_AQ_A}{r}\right),
  \end{equation}
where the approximation holds well for weak linkage, which we assumed when we
derived @eq:qle.

## Fixed point iteration algorithm  {#sec:fp}

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





