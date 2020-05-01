#include <iostream>
#include <string>
#include "SpeciesTree.h"
#include "Simulator.h"
#include "Engine.h"
#include <sstream>
#include <RcppArmadillo.h>

using namespace Rcpp;

std::vector<std::string> split(std::string strToSplit, char delimeter)
{
    std::stringstream ss(strToSplit);
    std::string item;
    std::vector<std::string> splittedStrings;
    while (std::getline(ss, item, delimeter))
    {
        splittedStrings.push_back(item);
    }
    return splittedStrings;
}


int run_treeducken(std::string params_file) {
    std::string setFileName;
    int mt;
    Engine *phyEngine;
    if(params_file.empty()){
        return 0;
    }
    else{
        std::vector<std::string> params_vector;

        params_vector = split(params_file, ' ');

        std::string outName;
        std::string stn;
        int nt = 100, r = 100, nloc = 0, ipp = 0, ne = 0, ngen = 0;
        double sbr = 0.5, sdr = 0.2, gbr = 0.0, gdr = 0.0, lgtr = 0.0;
        double ts = -1, og = 0.0;
        bool sout = 1;
        bool mst = 0;
        for (int i = 0; i < params_vector.size(); i++){
            if(params_vector[i] == "-sbr"){
                sbr = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-sdr"){
                sdr = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-gbr"){
                gbr = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-gdr"){
                gdr = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-lgtr"){
                lgtr = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-ne"){
                ne = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-ipp"){
                ipp = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-nt"){
                nt = std::atoi(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-sc"){
                ts = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-nl"){
                nloc = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-ng"){
                ngen = std::atoi(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-r"){
                r = std::atoi(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-o"){
                outName = params_vector[i+1];
            }
            else if(params_vector[i] == "-istnw"){
                stn = params_vector[i+1];
            }
            else if(params_vector[i] == "-og"){
                og = std::atof(params_vector[i+1].c_str());
            }
            else if(params_vector[i] == "-sout"){
                sout = std::atoi(params_vector[i+1].c_str());
            }
            else {
                //TODO: need error handling
                // Rcout << "\n############################ !!! ###########################\n";
                // Rcout << "\n\n\tPerhaps you mis-typed something, here are the \n\tavailable options:\n";
                // printHelp();
                // Rcout << "\n############################ !!! ###########################\n";
            }
        }
        if(stn != ""){
            mt = 4;
            // Rcout << "Species tree is set. Simulating only locus and gene trees...\n";
            if(nloc > 0){
                if (gbr <= 0.0){
                    if(gbr < 0.0){
                        // std::cerr << "Gene birth rate is a negative number, no loci or gene trees will be simulated.\n";
                        // std::cerr << "You have input a tree with no other parameters set to simulate within.\n";
                    }
                    else{
                        if(ne >  0 && ipp > 0 && ipp <= ne){
                            //   Rcout << "Gene birth rate is 0.0, locus trees will match species trees.\n";
                        }
                        else{
                            //    std::cerr << "Gene tree parameters are incorrectly specified. Only simulating species and locus trees\n";
                            //    std::cerr << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                        }
                    }
                    // printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);

                }
                else if (ne <= 0 || ipp <= 0 || ipp > ne){
                    //    std::cerr << "Gene tree parameters are incorrectly specified. Only simulating species and locus trees\n";
                    //    std::cerr << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                }
                else{
                    //    Rcout << "Simulating locus and gene trees on input species tree.\n";
                    // printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);
                }
            }
        }
        else{
            if(sbr <= 0.0){
                //    std::cerr << "Species birth rate of 0.0 is invalid, exiting...\n";
                return 0;
            }

            if(nloc > 0){
                if (gbr <= 0.0){
                    if(gbr < 0.0){
                        mt = 1;
                        //         Rcout << "Gene birth rate is a negative number, no loci or gene trees will be simulated.\n";
                    }
                    else{
                        if(ne >  0 && ipp > 0 && ipp <= ne){
                            mt = 3;
                            //            Rcout << "Gene birth rate is 0.0, locus trees will match species trees.\n";
                        }
                        else{
                            //            Rcout << "Gene tree parameters are incorrectly specified. Only simulating species and locus trees\n";
                            //            Rcout << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                            mt = 2;
                        }
                    }
                    //printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);

                }
                else if (ne <= 0 || ipp <= 0 || ipp > ne){
                    mt = 2;
                    //   Rcout << "Gene tree parameters are incorrectly specified.\n";
                    //    Rcout << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                    //printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);
                }
                else{
                    mt = 3;
                    //    Rcout << "Simulating sets of three trees.\n";
                    //printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);
                }
            }
            else{
                //    Rcout << "Number of loci to simulate is set to 0." << std::endl;
                mt = 1;
                //printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);
                if(mst)
                    mt = 5;
            }
        }

        phyEngine = new Engine(outName,
                               mt,
                               sbr,
                               sdr,
                               gbr,
                               gdr,
                               lgtr,
                               ipp,
                               ne,
                               1.0,
                               ts,
                               r,
                               nt,
                               nloc,
                               ngen,
                               og,
                               sout,
                               mst);
        if(stn != ""){
            phyEngine->setInputSpeciesTree(stn);
            phyEngine->doRunSpTreeSet();
        }
        else
            phyEngine->doRunRun();


    }
    delete phyEngine;

    return 0;
}


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
                                    double popsize,
                                    int samples_per_lineage,
                                    int numGenesPerLocus){
    Rcpp::List multiphy;
    int ntax = species_tree->getNumExtant();
    Rcout << ntax << std::endl;
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
                                                popsize,
                                                genTime,
                                                numGenesPerLocus,
                                                og,
                                                ts,
                                                sout);
        phySimulator->setSpeciesTree(species_tree);
        bool good = phySimulator->simLocusTree();
        Rcout << good << std::endl;

        List phyGenesPerLoc(numGenesPerLocus);
        for(int j=0; j<numGenesPerLocus; j++){
            phySimulator->simGeneTree();

            List phyGene = List::create(Named("edge") = phySimulator->getGeneEdges(),
                         _("edge.length") = phySimulator->getGeneEdgeLengths(),
                         _("Nnode") = phySimulator->getGeneNnodes(),
                         _("tip.label") = phySimulator->getGeneTipNames(),
                         _("root.edge") = phySimulator->getGeneTreeRootEdge());
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