#include <iostream>
#include <string>
#include "SpeciesTree.h"
#include "Simulator.h"
#include "Engine.h"
#include <sstream>
#include <RcppArmadillo.h>

using namespace Rcpp;

void printHelp();
void printSettings();
void printVersion();
//
// void printHelp(){
//     Rcout << "\tHere are the available options that you can change (default values are in []):\n";
//     Rcout << "\t\t-i    : input settings file \n";
//     Rcout << "\t\t-o    : output file name prefix [= ""]\n";
//     Rcout << "\t\t-nt   : number of extant taxa [= 100]\n";
//     Rcout << "\t\t-r    : number of replicates [= 10]\n";
//     Rcout << "\t\t-nl : number of loci to simulate [= 0]\n";
//     Rcout << "\t\t-sc   : tree scale [= 1.0]\n";
//     Rcout << "\t\t-sd1  : seed 1 (use this if you only pass in one seed) \n";
//     Rcout << "\t\t-sd2  : seed 2 \n";
//     Rcout << "\t\t-sbr  : species birth rate [= 0.5]\n";
//     Rcout << "\t\t-sdr  : species death rate [= 0.2]\n";
//     Rcout << "\t\t-gbr  : gene birth rate [= 0.0]\n";
//     Rcout << "\t\t-gdr  : gene death rate [= 0.0]\n";
//     Rcout << "\t\t-lgtr : gene transfer rate [= 0.0]\n";
//     Rcout << "\t\t-ipp  : individuals to sample per locus [= 0]\n";
//     Rcout << "\t\t-ne   : effective population size per locus [= 0] \n";
//     Rcout << "\t\t-ng   : number of genes to simulate per locus [= 0] \n";
//     Rcout << "\t\t-og   : fraction of tree to use as length of branch between outgroup [=0.0] \n" ;
//     Rcout << "\t\t-istnw  : input species tree (newick format) [=""] \n";
//     Rcout << "\t\t-sc     : tree scale [=1.0] \n";
//     Rcout << "\t\t-sout   : turn off standard output (improves runtime) \n";
//     Rcout << "\t\t-mst    : Moran species tree ";
// }
//
void printSettings(std::string of, int nt, int r, int nloc, int ts, double sbr, double sdr,
                   double gbr, double gdr, double lgtr, int ipp, int ne, int ngen, double og,
                   std::string stn, bool mst){
    Rcout << "\t\toutput file name prefix         = " << of << "\n";
    Rcout << "\t\tNumber of extant taxa           = " << nt << "\n";
    Rcout << "\t\tNumber of replicates            = " << r << "\n";
    Rcout << "\t\tNumber of loci to simulate      = " << nloc << "\n";
    Rcout << "\t\tNumber of gene trees to simulate     = " << ngen << "\n";
    Rcout << "\t\tSpecies birth rate              = " << sbr << "\n";
    Rcout << "\t\tSpecies death rate              = " << sdr << "\n";
    Rcout << "\t\tGene birth rate                 = " << gbr << "\n";
    Rcout << "\t\tGene death rate                 = " << gdr << "\n";
    Rcout << "\t\tGene transfer rate              = " << lgtr << "\n";
    Rcout << "\t\tIndividuals to sample per locus = " << ipp << "\n";
    Rcout << "\t\tEffective pop size per locus    = " << ne << "\n";
    Rcout << "\t\tTree fraction to set outgroup   = " << og << "\n";
    Rcout << "\t\tSpecies tree input as newick    = " << stn << "\n";
    Rcout << "\t\tTree scale                      = " << ts  << "\n";
    Rcout << "\t\tMoran process species tree      = " << mst << "\n";
}
//
// void printVersion(){
//     Rcout << "############################################################\n";
//     Rcout << "####\ttreeducken, version 0.1 \t\t\t####\n";
//     Rcout << "####\t" << GIT_HASH << "\t####\n";
//     Rcout << "############################################################\n";
// }
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
                        exit(1);
                    }
                    else{
                        if(ne >  0 && ipp > 0 && ipp <= ne){
                            //   Rcout << "Gene birth rate is 0.0, locus trees will match species trees.\n";
                        }
                        else{
                            //    std::cerr << "Gene tree parameters are incorrectly specified. Only simulating species and locus trees\n";
                            //    std::cerr << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                            exit(1);
                        }
                    }
                    // printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);

                }
                else if (ne <= 0 || ipp <= 0 || ipp > ne){
                    //    std::cerr << "Gene tree parameters are incorrectly specified. Only simulating species and locus trees\n";
                    //    std::cerr << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                    exit(1);
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
        printSettings(outName,
                      nt,
                      r,
                      nloc,
                      ts,
                      sbr,
                      sdr,
                      gbr,
                      gdr,
                      lgtr,
                      ipp,
                      ne,
                      ngen,
                      og,
                      stn,
                      mst);

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

    std::vector<SpeciesTree*> speciesTrees;
    speciesTrees.resize(numbsim);
    Simulator *phySimulator;


    List multiphy(numbsim);
    for(int i = 0; i < numbsim; i++){
        phySimulator = new Simulator(n_tips,
                                     sbr,
                                     sdr,
                                     1);
        phySimulator->simSpeciesTree();
        speciesTrees[i] = phySimulator->getSpeciesTree();

        NumericMatrix edges_rmat = speciesTrees[i]->getEdges();


        List phy = List::create(Named("edge") = edges_rmat,
                                Named("edge.length") = speciesTrees[i]->getEdgeLengths(),
                                Named("Nnode") = speciesTrees[i]->getNnodes(),
                                Named("tip.label") = speciesTrees[i]->getTipNames(),
                                Named("root.edge") = speciesTrees[i]->getRoot()->getDeathTime() -
                                    speciesTrees[i]->getRoot()->getBirthTime());
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
    Simulator *phySimulator = nullptr;
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
        Rcout << "********" << std::endl;

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