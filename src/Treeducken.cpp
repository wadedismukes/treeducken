#include <iostream>
#include <string>
#include "SpeciesTree.h"
#include "Simulator.h"
#include "Engine.h"

void printHelp();
void printSettings();
void printVersion();
//
// void printHelp(){
//     std::cout << "\tHere are the available options that you can change (default values are in []):\n";
//     std::cout << "\t\t-i    : input settings file \n";
//     std::cout << "\t\t-o    : output file name prefix [= ""]\n";
//     std::cout << "\t\t-nt   : number of extant taxa [= 100]\n";
//     std::cout << "\t\t-r    : number of replicates [= 10]\n";
//     std::cout << "\t\t-nl : number of loci to simulate [= 0]\n";
//     std::cout << "\t\t-sc   : tree scale [= 1.0]\n";
//     std::cout << "\t\t-sd1  : seed 1 (use this if you only pass in one seed) \n";
//     std::cout << "\t\t-sd2  : seed 2 \n";
//     std::cout << "\t\t-sbr  : species birth rate [= 0.5]\n";
//     std::cout << "\t\t-sdr  : species death rate [= 0.2]\n";
//     std::cout << "\t\t-gbr  : gene birth rate [= 0.0]\n";
//     std::cout << "\t\t-gdr  : gene death rate [= 0.0]\n";
//     std::cout << "\t\t-lgtr : gene transfer rate [= 0.0]\n";
//     std::cout << "\t\t-ipp  : individuals to sample per locus [= 0]\n";
//     std::cout << "\t\t-ne   : effective population size per locus [= 0] \n";
//     std::cout << "\t\t-ng   : number of genes to simulate per locus [= 0] \n";
//     std::cout << "\t\t-og   : fraction of tree to use as length of branch between outgroup [=0.0] \n" ;
//     std::cout << "\t\t-istnw  : input species tree (newick format) [=""] \n";
//     std::cout << "\t\t-sc     : tree scale [=1.0] \n";
//     std::cout << "\t\t-sout   : turn off standard output (improves runtime) \n";
//     std::cout << "\t\t-mst    : Moran species tree ";
// }
//
// void printSettings(std::string of, int nt, int r, int nloc, int ts, double sbr, double sdr,
//                    double gbr, double gdr, double lgtr, int ipp, int ne, int ngen, double og,
//                    std::string stn, bool mst){
//     std::cout << "\t\toutput file name prefix         = " << of << "\n";
//     std::cout << "\t\tNumber of extant taxa           = " << nt << "\n";
//     std::cout << "\t\tNumber of replicates            = " << r << "\n";
//     std::cout << "\t\tNumber of loci to simulate      = " << nloc << "\n";
//     std::cout << "\t\tNumber of gene trees to simulate     = " << ngen << "\n";
//     std::cout << "\t\tSpecies birth rate              = " << sbr << "\n";
//     std::cout << "\t\tSpecies death rate              = " << sdr << "\n";
//     std::cout << "\t\tGene birth rate                 = " << gbr << "\n";
//     std::cout << "\t\tGene death rate                 = " << gdr << "\n";
//     std::cout << "\t\tGene transfer rate              = " << lgtr << "\n";
//     std::cout << "\t\tIndividuals to sample per locus = " << ipp << "\n";
//     std::cout << "\t\tEffective pop size per locus    = " << ne << "\n";
//     std::cout << "\t\tTree fraction to set outgroup   = " << og << "\n";
//     std::cout << "\t\tSpecies tree input as newick    = " << stn << "\n";
//     std::cout << "\t\tTree scale                      = " << ts  << "\n";
//     std::cout << "\t\tMoran process species tree      = " << mst << "\n";
// }
//
// void printVersion(){
//     std::cout << "############################################################\n";
//     std::cout << "####\ttreeducken, version 0.1 \t\t\t####\n";
//     std::cout << "####\t" << GIT_HASH << "\t####\n";
//     std::cout << "############################################################\n";
// }

int run_treeducken(std::string params_file) {
    std::string setFileName;
    int mt;
    Engine *phyEngine;
    if(params_file.empty()){
        return 0;
    }
    else{
        std::vector<std::string> params_vector;
        int startIndex = 0;
        int  endIndex = 0;
        while( (endIndex = params_file.find(' ', startIndex)) < params_file.size() )
        {

            std::string val = params_file.substr(startIndex, endIndex - startIndex);
            params_vector.push_back(val);
            startIndex = endIndex + 1;

        }
        if(startIndex < params_file.size())
        {
            std::string val = params_file.substr(startIndex);
            params_vector.push_back(val);
        }

        std::string outName;
        std::string stn;
        int nt = 100, r = 10, nloc = 10, ipp = 0, ne = 0, sd1 = 0, sd2 = 0, ngen = 0;
        double sbr = 0.5, sdr = 0.2, gbr = 0.0, gdr = 0.0, lgtr = 0.0, ts = -1, og = 0.0;
        bool sout = 1;
        bool mst = 0;
        for (int i = 0; i < params_vector.size(); i++){
            if(params_vector[i] == "-sbr"){
                sbr = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-sdr"){
                sdr = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-gbr"){
                gbr = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-gdr"){
                gdr = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-lgtr"){
                lgtr = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-ne"){
                ne = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-ipp"){
                ipp = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-nt"){
                nt = std::atoi(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-sc"){
                ts = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-nl"){
                nloc = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-ng"){
                ngen = std::atoi(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-r"){
                r = std::atoi(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-o"){
                outName = params_vector[i+1];
                i++;
            }
            else if(params_vector[i] == "-istnw"){
                stn = params_vector[i+1];
                i++;
            }
            else if(params_vector[i] == "-og"){
                og = std::atof(params_vector[i+1].c_str());
                i++;
            }
            else if(params_vector[i] == "-sout"){
                sout = std::atoi(params_vector[i+1].c_str());
                i++;
            }
            else {
                //TODO: need error handling
                // std::cout << "\n############################ !!! ###########################\n";
                // std::cout << "\n\n\tPerhaps you mis-typed something, here are the \n\tavailable options:\n";
                // printHelp();
                // std::cout << "\n############################ !!! ###########################\n";
            }
        }
        if(stn != ""){
            mt = 4;
            // std::cout << "Species tree is set. Simulating only locus and gene trees...\n";
            if(nloc > 0){
                if (gbr <= 0.0){
                    if(gbr < 0.0){
                        // std::cerr << "Gene birth rate is a negative number, no loci or gene trees will be simulated.\n";
                        // std::cerr << "You have input a tree with no other parameters set to simulate within.\n";
                        exit(1);
                    }
                    else{
                        if(ne >  0 && ipp > 0 && ipp <= ne){
                            //   std::cout << "Gene birth rate is 0.0, locus trees will match species trees.\n";
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
                    //    std::cout << "Simulating locus and gene trees on input species tree.\n";
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
                        //         std::cout << "Gene birth rate is a negative number, no loci or gene trees will be simulated.\n";
                    }
                    else{
                        if(ne >  0 && ipp > 0 && ipp <= ne){
                            mt = 3;
                            //            std::cout << "Gene birth rate is 0.0, locus trees will match species trees.\n";
                        }
                        else{
                            //            std::cout << "Gene tree parameters are incorrectly specified. Only simulating species and locus trees\n";
                            //            std::cout << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                            mt = 2;
                        }
                    }
                    //printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);

                }
                else if (ne <= 0 || ipp <= 0 || ipp > ne){
                    mt = 2;
                    //   std::cout << "Gene tree parameters are incorrectly specified.\n";
                    //    std::cout << "Population size and individuals per population must both be positive integers and individuals per population must be less than or equal to the population size.\n";
                    //printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);
                }
                else{
                    mt = 3;
                    //    std::cout << "Simulating sets of three trees.\n";
                    //printSettings(outName, nt, r, nloc, ts, sbr, sdr, gbr, gdr, lgtr, ipp, ne, ngen, og, stn, mst);
                }
            }
            else{
                //    std::cout << "Number of loci to simulate is set to 0." << std::endl;
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