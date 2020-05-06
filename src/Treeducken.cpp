#include <iostream>
#include <string>
#include "SpeciesTree.h"
#include "Simulator.h"
#include <sstream>
#include <RcppArmadillo.h>

using namespace Rcpp;


Rcpp::List bdsim_species_tree(double sbr,
                        double sdr,
                        int numbsim,
                        int n_tips){

    Simulator *phySimulator;
    List multiphy(numbsim);
    for(int i = 0; i < numbsim; i++){
        phySimulator = new Simulator(n_tips,
                                     sbr,
                                     sdr,
                                     1);
        phySimulator->simSpeciesTree();

        List phy = List::create(Named("edge") = phySimulator->getSpeciesEdges(),
                                Named("edge.length") = phySimulator->getSpeciesEdgeLengths(),
                                Named("Nnode") = phySimulator->getSpeciesNnodes(),
                                Named("tip.label") = phySimulator->getSpeciesTipNames(),
                                Named("root.edge") = phySimulator->getSpeciesTreeRootEdge());;
        phy.attr("class") = "phylo";
        multiphy[i] = phy;
    }




    delete phySimulator;
    return multiphy;
}

Rcpp::List sim_bdsimple_species_tree(double sbr,
                                     double sdr,
                                     int numbsim,
                                     double timeToSimTo){
    Simulator *phySimulator;
    List multiphy(numbsim);
    for(int i = 0; i < numbsim; i++){
        phySimulator = new Simulator(1,
                                     sbr,
                                     sdr,
                                     1);
        phySimulator->simSpeciesTreeTime();

        List phy = List::create(Named("edge") = phySimulator->getSpeciesEdges(),
                                Named("edge.length") = phySimulator->getSpeciesEdgeLengths(),
                                Named("Nnode") = phySimulator->getSpeciesNnodes(),
                                Named("tip.label") = phySimulator->getSpeciesTipNames(),
                                Named("root.edge") = phySimulator->getSpeciesTreeRootEdge());;
        phy.attr("class") = "phylo";
        multiphy[i] = phy;
    }

    delete phySimulator;
    return multiphy;
}

Rcpp::List sim_locus_tree(SpeciesTree* species_tree,
                          double gbr,
                          double gdr,
                          double lgtr,
                          int numbsim){
    Rcpp::List multiphy;
    int ntax = species_tree->getNumExtant();
    double lambda = 0.0;
    double mu = 0.0;
    double rho = 0.0;
    unsigned numLociToSim = numbsim;
    for(int i = 0; i < numbsim; i++){
        Simulator *phySimulator = new Simulator( ntax,
                                     lambda,
                                     mu,
                                     rho,
                                     numLociToSim,
                                     gbr,
                                     gdr,
                                     lgtr);
        phySimulator->setSpeciesTree(species_tree);
        phySimulator->simLocusTree();
        List phy = List::create(Named("edge") = phySimulator->getLocusEdges(),
                                Named("edge.length") = phySimulator->getLocusEdgeLengths(),
                                Named("Nnode") = phySimulator->getLocusNnodes(),
                                Named("tip.label") = phySimulator->getLocusTipNames(),
                                Named("root.edge") = phySimulator->getLocusTreeRootEdge());

        phy.attr("class") = "phylo";

        multiphy.push_back(std::move(phy));
        delete phySimulator;
    }



    return multiphy;
}

Rcpp::List sim_locus_tree_gene_tree(SpeciesTree* species_tree,
                                    double gbr,
                                    double gdr,
                                    double lgtr,
                                    int numLoci,
                                    double theta,
                                    int samples_per_lineage,
                                    int numGenesPerLocus){
    Rcpp::List multiphy;
    int ntax = species_tree->getNumExtant();
    double lambda = 0.0;
    double mu = 0.0;
    double rho = 1.0;
    unsigned numLociToSim = numLoci;
    double genTime = 1.0;
    double ts = 1.0;
    bool sout = false;
    double og = 0.0;
    for(int i = 0; i < numLoci; i++){
        Simulator *phySimulator = new Simulator(ntax,
                                                lambda,
                                                mu,
                                                rho,
                                                numLociToSim,
                                                gbr,
                                                gdr,
                                                lgtr,
                                                samples_per_lineage,
                                                theta,
                                                genTime,
                                                numGenesPerLocus,
                                                og,
                                                ts,
                                                sout);
        phySimulator->setSpeciesTree(species_tree);
        phySimulator->simLocusTree();

        List phyGenesPerLoc(numGenesPerLocus);
        for(int j=0; j<numGenesPerLocus; j++){
            phySimulator->simGeneTree(j);

            List phyGene = List::create(Named("edge") = phySimulator->getGeneEdges(j),
                         _("edge.length") = phySimulator->getGeneEdgeLengths(j),
                         _("Nnode") = phySimulator->getGeneNnodes(j),
                         _("tip.label") = phySimulator->getGeneTipNames(j),
                         _("root.edge") = phySimulator->getGeneTreeRootEdge(j));
            phyGene.attr("class") = "phylo";
            phyGenesPerLoc[j] = phyGene;
        }
        List phyLoc = List::create(Named("edge") = phySimulator->getLocusEdges(),
                                   Named("edge.length") = phySimulator->getLocusEdgeLengths(),
                                   Named("Nnode") = phySimulator->getLocusNnodes(),
                                   Named("tip.label") = phySimulator->getLocusTipNames(),
                                   Named("root.edge") = phySimulator->getLocusTreeRootEdge());
        phyLoc.attr("class") = "phylo";
        List locusGeneSet = List::create(Named("locus.tree") = phyLoc,
                                         Named("gene.trees") = phyGenesPerLoc);

        multiphy.push_back(std::move(locusGeneSet));

        delete phySimulator;
    }



    return multiphy;
}


Rcpp::List sim_genetree_msc(SpeciesTree* species_tree,
                            double popsize,
                            int samples_per_lineage,
                            int numbsim){
    return sim_locus_tree_gene_tree(species_tree,
                             0.0,
                             0.0,
                             0.0,
                             1,
                             samples_per_lineage,
                             popsize,
                             numbsim);
    // this one is a wrapper for above function with locus tree parameters set to 0
}