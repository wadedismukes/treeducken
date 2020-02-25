//
//  Simulator.hpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/9/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef Simulator_h
#define Simulator_h
#include "GeneTree.h"
#include "SymbiontTree.h"
#include <set>
#include <map>
#include <RcppArmadillo.h>

class Simulator
{
    protected:
        double      currentSimTime;
        unsigned    simType;
        unsigned    numTaxaToSim, gsaStop;
        unsigned    numLoci;
        unsigned    numGenes;
        double      speciationRate, extinctionRate;
        double      samplingRate;
        double      treeScale;
        double      geneBirthRate, geneDeathRate, transferRate;
        double      propTransfer, propDuplicate;
        unsigned    indPerPop;
        unsigned    popSize;
        double      generationTime;
        double      outgroupFrac;
        bool        printSOUT;
        std::vector<SpeciesTree*>   gsaTrees;
        SpeciesTree*    spTree;
        LocusTree*      lociTree;
        std::vector<LocusTree*> locusTrees;
        GeneTree*       geneTree;
        std::vector<std::vector<GeneTree*> > geneTrees;
        // symbiont tree stuff
        SymbiontTree*   symbiontTree;
        double      cospeciationRate;
        double      timeToSim;
        int         hostLimit;
        arma::umat   assocMat;

    public:
        // Simulating species tree only
        Simulator(unsigned numTaxaToSim,
                  double speciationRate,
                  double extinctionRate,
                  double rho);
        // Simulating species and locus tree
        Simulator(unsigned numTaxaToSim,
                  double speciationRate,
                  double extinctionRate,
                  double rho,
                  unsigned numLociToSim,
                  double geneBirthRate,
                  double geneDeathRate,
                  double transferRate);
        // Simulating species and locus tree with proportion of transfer (e.g. hybridization, linkage)
        // // CAN THIS BE DELETED??
        //.
        Simulator(unsigned numTaxaToSim,
                  double speciationRate,
                double extinctionRate,
                double rho,
                unsigned numLociToSim,
                double geneBirthRate,
                double geneDeathRate,
                double transferRate,
                double propTransfer);
        // Simulating species and locus trees with one gene tree per locus tree
        Simulator(unsigned numTaxaToSim,
                double speciationRate,
                double extinctionRate,
                double rho,
                unsigned numLociToSim,
                double geneBirthRate,
                double geneDeathRate,
                double transferRate,
                unsigned indPerPop,
                unsigned popSize,
                double genTime,
                int ng,
                double og,
                double ts,
                bool sout);
        //
        Simulator(double timeToSimTo,
                  double hostSpeciationRate,
                  double hostExtinctionRate,
                  double symbSpeciationRate,
                  double symbExtinctionRate,
                  double switchingRate,
                  double cospeciationRate,
                  double rho,
                  int hostLimit);
        ~Simulator();

        void    setSpeciesTree(SpeciesTree *st) { spTree = st; }
        bool    gsaBDSim();
        bool    bdsaBDSim();
        bool    pairedBDPSim();
        bool    coalescentSim();
        bool    simSpeciesTree();
        bool    simLocusTree();
        bool    simSpeciesLociTrees();
        bool    simThreeTree();
        bool    simLocusGeneTrees();
        bool    simHostSymbSpeciesTreePair();
        bool    gsaCheckStop();
        void    initializeSim();
        void    processGSASim();
        void    prepGSATreeForReconstruction();
        void    processSpTreeSim();
        void    graftOutgroup(Tree *tr, double trDepth);
        double  calcSpeciesTreeDepth();
        double  calcExtantSpeciesTreeDepth();
        double  calcLocusTreeDepth(int i);
        int     findNumberTransfers();
        double  findTMRCAGeneTree(int i, int j);
        std::string    printSpeciesTreeNewick();
        std::string    printExtSpeciesTreeNewick();
        std::string    printLocusTreeNewick(int i);
        std::string    printGeneTreeNewick(int i, int j);
        std::string    printExtantGeneTreeNewick(int i, int j);
        std::set<double, std::greater<double> > getEpochs();
        //SpeciesTree*    getSpeciesTree() {SpeciesTree* spec_tree = new SpeciesTree(*spTree); return spec_tree;}
        SpeciesTree*    getSpeciesTree() { return spTree; }
        LocusTree*      getLocusTree() {return lociTree;}
        SymbiontTree*   getSymbiontTree() {return symbiontTree;}
        double          getTimeToSim() {return timeToSim; }

        NumericMatrix   getSymbiontEdges() { return symbiontTree->getEdges(); }
        NumericMatrix   getSpeciesEdges() { return spTree->getEdges(); }
        NumericMatrix   getLocusEdges() { return lociTree->getEdges(); }

        std::vector<double>    getSymbiontEdgeLengths() { return symbiontTree->getEdgeLengths(); }
        std::vector<double>    getSpeciesEdgeLengths() { return spTree->getEdgeLengths(); }
        std::vector<double>    getLocusEdgeLengths() { return lociTree->getEdgeLengths(); }

        int    getSymbiontNnodes() { return symbiontTree->getNnodes(); }
        int    getSpeciesNnodes() { return spTree->getNnodes(); }
        int    getLocusNnodes() { return lociTree->getNnodes(); }

        std::vector<std::string> getSpeciesTipNames() { return spTree->getTipNames(); }
        std::vector<std::string> getSymbiontTipNames() { return symbiontTree->getTipNames(); }
        std::vector<std::string> getLocusTipNames() { return lociTree->getTipNames(); }

        double    getSpeciesTreeRootEdge();
        double    getLocusTreeRootEdge();
        double    getSymbiontTreeRootEdge();
        arma::umat    getAssociationMatrix() { return assocMat; }
        arma::umat    cophyloEvent(double eventTime, arma::umat assocMat);
        arma::umat    cophyloERMEvent(double eventTime, arma::umat assocMat);
        arma::umat    cospeciationEvent(double eventTime, arma::umat assocMat);
        //List  getEventHistory() { return updateIndices(events); }
       //List  updateIndices(List evs);
};


#endif /* Simulator_h */
