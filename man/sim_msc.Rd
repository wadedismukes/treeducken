% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/deprecated_sim_functions.R
\name{sim_msc}
\alias{sim_msc}
\alias{sim_multispecies_coal}
\title{Simulate multispecies coalescent on a species tree}
\usage{
sim_msc(
  species_tree,
  ne,
  num_sampled_individuals,
  num_genes,
  rescale = TRUE,
  mutation_rate = 1L,
  generation_time = 1L
)

sim_multispecies_coal(
  species_tree,
  ne,
  num_sampled_individuals,
  num_genes,
  rescale = TRUE,
  mutation_rate = 1
)
}
\arguments{
\item{species_tree}{input species tree of class "phylo"}

\item{ne}{Effective population size}

\item{num_sampled_individuals}{number of individuals sampled within each lineage}

\item{num_genes}{number of genes to simulate within each locus}

\item{rescale}{Rescale the tree into coalescent units (otherwise assumes it is in those units)}

\item{mutation_rate}{The rate of mutation per generation}

\item{generation_time}{The number of time units per generation}
}
\value{
A list of coalescent trees
}
\description{
Simulates the multispecies coalescent on a species tree.
}
\details{
This a multispecies coalescent simulator with two usage options.
The function can rescale the given tree into coalescent units given the `mutation_rate`, `ne`, and the `generation_time`.
These result in a tree with coalescent times in units of expected number of mutations per site.
The generation_time parameter default is 1 time unit per generation if the units of the tree are in millions of years
The mutation_rate parameter is by default set to 1 mutations per site per generation (which is nonsensical).
Rescale is set to true by default.

If rescale is set to false the tree is assumed to be in coalescent units and `ne` is used as the population
genetic parameter theta.
}
\examples{
# first simulate a species tree
mu <- 0.5
lambda <- 1.0
nt <- 6
tr <- sim_stBD(sbr = lambda, sdr = mu, numbsim = 1, n_tips = nt)
# for a locus tree with 100 genes sampled per locus tree
gentrees <- sim_msc(tr[[1]],
                    ne = 10000,
                    mutation_rate = 1e-9,
                    generation_time = 1e-6,
                    num_sampled_individuals = 1,
                    num_genes = 100)

}
\references{
Bruce Rannala and Ziheng Yang (2003) Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences From Multiple Loci Genetics August 1, 2003 vol. 164 no. 4 1645-1656
Mallo D, de Oliveira Martins L, Posada D (2015) SimPhy: Phylogenomic Simulation of Gene, Locus and Species Trees. Syst. Biol. doi: http://dx.doi.org/10.1093/sysbio/syv082
}
\seealso{
sim_ltBD, sim_stBD, sim_stBD_t
}
