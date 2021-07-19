#ifndef Simulator_h
#define Simulator_h
#include "GeneTree.h"
#include "SymbiontTree.h"
#include <set>
#include <map>
#include <RcppArmadillo.h>

class Simulator
{
    private:
        double      currentSimTime;
        unsigned    numTaxaToSim, gsaStop;
        unsigned    numLoci;
        unsigned    numGenes;
        double      speciationRate, extinctionRate;
        double      samplingRate;
        double      geneBirthRate, geneDeathRate, transferRate;
        double      propTransfer, propDuplicate;
        double      dispersalRate, extirpationRate;
        unsigned    indPerPop;
        double      popSize;
        double      generationTime;
        bool        host_switch_mode;
        std::vector<std::shared_ptr<SpeciesTree>>   gsaTrees;
        std::shared_ptr<SpeciesTree>    spTree;
        std::shared_ptr<LocusTree>      lociTree;
        std::vector<std::shared_ptr<LocusTree>> locusTrees;
        std::shared_ptr<GeneTree>       geneTree;
        std::vector<std::shared_ptr<GeneTree>> geneTrees;
        // symbiont tree varibles
        std::shared_ptr<SymbiontTree>   symbiontTree;
        double      cospeciationRate;
        double      timeToSim;
        int         hostLimit;
        arma::umat   assocMat;
        std::string  transferType;

        Rcpp::IntegerVector inOrderVecOfHostIndx;
        Rcpp::IntegerVector inOrderVecOfSymbIndx;
        Rcpp::CharacterVector inOrderVecOfEvent;
        Rcpp::NumericVector inOrderVecOfEventTimes;

    public:
        struct SpeciesSimTime
        {   
            inline SpeciesSimTime(double speciesBirthRate,
                    double speciesDeathRate,
                    double time):
                        sbr(speciesBirthRate),
                        sdr(speciesDeathRate),
                        t(time) {}
            double sbr;
            double sdr;
            double t;
            double currentSimTime = 0;

        };
        struct SpeciesSimTips
        {
            inline SpeciesSimTips(double speciesBirthRate,
            double speciesDeathRate,
            unsigned gsaStop,
            unsigned targetTips):
                        sbr(speciesBirthRate),
                        sdr(speciesDeathRate),
                        gsaStop(gsaStop),
                        numTaxaToSim(targetTips) {}
            double sbr;
            double sdr;
            unsigned gsaStop;
            unsigned numTaxaToSim;
     
        };

   
        struct LocusSim {
                inline LocusSim(std::shared_ptr<SpeciesTree> species_tree,
                        double gbr,
                        double gdr,
                        double lgtr,
                        std::string trans_type):
                            spTree(species_tree),
                            geneBirthRate(gbr),
                            geneDeathRate(gdr),
                            transferRate(lgtr),
                            transferType(trans_type) {}
                std::shared_ptr<SpeciesTree> spTree;
                double geneBirthRate;
                double geneDeathRate;
                double transferRate;
                std::string transferType;

        };  
        struct LocusAndGeneMSCSim {
                inline LocusAndGeneMSCSim(std::shared_ptr<SpeciesTree> species_tree,
                                        double gbr,
                                        double gdr,
                                        double lgtr,
                                        double popsize,
                                        int numLoci,
                                        int samples_per_lineage,
                                        int numbsim):
                    spTree(species_tree),
                    geneBirthRate(gbr),
                    geneDeathRate(gdr),
                    transferRate(lgtr),
                    popsize(popsize),
                    numLoci(numLoci),
                    samples_per_lineage(samples_per_lineage),
                    numbsim(numbsim) {}
                std::shared_ptr<SpeciesTree> spTree;
                double geneBirthRate;
                double geneDeathRate;
                double transferRate;
                double popsize;
                unsigned numLoci;
                int samples_per_lineage;
                int numbsim;

        };
        struct CophySim {
            inline CophySim(double hostbr,
                                double hostdr,
                                double symbbr,
                                double symbdr,
                                double switchrate,
                                double cosprate,
                                double timeToSimTo,
                                int host_limit,
                                bool hsMode): 
                                    hostBirthRate(hostbr),
                                    hostDeathRate(hostdr),
                                    symbDeathRate(symbdr),
                                    symbBirthRate(symbbr),
                                    switchingRate(switchrate),
                                    cospeciationRate(cosprate),
                                    timeToSimTo(timeToSimTo),
                                    hostLimit(host_limit),
                                    hsMode(hsMode) {}
            double hostBirthRate;
            double hostDeathRate;
            double symbDeathRate;
            double symbBirthRate;
            double switchingRate;
            double cospeciationRate;
            double timeToSimTo;
            unsigned hostLimit;
            bool hsMode;
        };
        struct CophySimAna {

            inline CophySimAna(double hostbr,
                                double hostdr,
                                double symbbr,
                                double symbdr,
                                double switchrate,
                                double cosprate,
                                double timeToSimTo,
                                int host_limit,
                                bool hsMode,
                                double symbdispersal,
                                double symbextirpation): 
                                    hostBirthRate(hostbr),
                                    hostDeathRate(hostdr),
                                    symbDeathRate(symbdr),
                                    symbBirthRate(symbbr),
                                    switchingRate(switchrate),
                                    cospeciationRate(cosprate),
                                    timeToSimTo(timeToSimTo),
                                    hostLimit(host_limit),
                                    hsMode(hsMode),
                                    symbDispRate(symbdispersal),
                                    symbExtRate(symbextirpation) {}
            double hostBirthRate;
            double hostDeathRate;
            double symbDeathRate;
            double symbBirthRate;
            double switchingRate;
            double cospeciationRate;
            double timeToSimTo;
            unsigned hostLimit;
            bool hsMode;
            double symbDispRate;
            double symbExtRate;
                        
        };

        Simulator(const SpeciesSimTips &spSimTips);
        Simulator(const SpeciesSimTime &spSimTime);
        Simulator(const LocusSim &locSim);
        Simulator(const LocusAndGeneMSCSim &lgMSC);
        Simulator(const CophySim &cophy);
        Simulator(const CophySimAna &cophyAna);


        ~Simulator();
        void    setGSAStop(int g) { gsaStop = g; }
        void    setSpeciesTree(std::shared_ptr<SpeciesTree> st) { spTree = st; }
        void    setLocusTree(std::shared_ptr<LocusTree> lt) { lociTree = lt; }

        bool    gsaBDSim();
        bool    bdsaBDSim();
        bool    bdSimpleSim();
        bool    pairedBDPSim();
        bool    pairedBDPSimAna();
        bool    coalescentSim();
        bool    simSpeciesTree();
        bool    simSpeciesTreeTime();
        bool    simLocusTree();
        bool    simGeneTree(int j);
        bool    simHostSymbSpeciesTreePair();
        bool    simHostSymbSpeciesTreePairWithAnagenesis();
        void    initializeSim();
        void    processGSASim();
        void    prepGSATreeForReconstruction();
        void    processSpTreeSim();
        double  calcSpeciesTreeDepth();
        double  calcExtantSpeciesTreeDepth();
        double  calcLocusTreeDepth(int i);
        int     findNumberTransfers();
        std::set<double, std::greater<double> > getEpochs();
        //SpeciesTree*    getSpeciesTree() {SpeciesTree* spec_tree = new SpeciesTree(*spTree); return spec_tree;}
        std::shared_ptr<SpeciesTree>    getSpeciesTree() { return spTree; }
        std::shared_ptr<LocusTree>      getLocusTree() {return lociTree;}
        std::shared_ptr<SymbiontTree>   getSymbiontTree() {return symbiontTree;}
        std::shared_ptr<GeneTree>       getGeneTree() {return geneTree; }
        double          getTimeToSim() {return timeToSim; }
        void            setTimeToSim(double tts) {timeToSim = tts; }
        NumericMatrix   getSymbiontEdges() { return symbiontTree->getEdges(); }
        NumericMatrix   getSpeciesEdges() { return spTree->getEdges(); }
        NumericMatrix   getLocusEdges() { return lociTree->getEdges(); }
        NumericMatrix   getGeneEdges(int j) { return geneTrees[j]->getGeneEdges(); }

        std::vector<double>    getSymbiontEdgeLengths() { return symbiontTree->getEdgeLengths(); }
        std::vector<double>    getSpeciesEdgeLengths() { return spTree->getEdgeLengths(); }
        std::vector<double>    getLocusEdgeLengths() { return lociTree->getEdgeLengths(); }
        std::vector<double>    getGeneEdgeLengths(int j) { return geneTrees[j]->getEdgeLengths(); }

        int    getSymbiontNnodes() { return symbiontTree->getNnodes(); }
        int    getSpeciesNnodes() { return spTree->getNnodes(); }
        int    getLocusNnodes() { return lociTree->getNnodes(); }
        int    getGeneNnodes(int j) { return geneTrees[j]->getNnodes(); }

        std::vector<std::string> getSpeciesTipNames() { return spTree->getTipNames(); }
        std::vector<std::string> getSymbiontTipNames() { return symbiontTree->getTipNames(); }
        std::vector<std::string> getLocusTipNames() { return lociTree->getTipNames(); }
        std::vector<std::string> getGeneTipNames(int j) { return geneTrees[j]->getTipNames(); }

        std::vector<std::string> getLocusTreeNodeLabels() { return lociTree->getNodeLabels(); }

        double    getSpeciesTreeRootEdge();
        double    getLocusTreeRootEdge();
        double    getSymbiontTreeRootEdge();
        double    getGeneTreeRootEdge(int j);
        arma::umat    hostLimitCheck(arma::umat assocMat, int hostLimit);
        arma::umat    getAssociationMatrix() { return assocMat; }
        arma::umat    cophyloEvent(double eventTime, arma::umat assocMat);
        arma::umat    cophyloERMEvent(double eventTime, arma::umat assocMat);
        arma::umat    cospeciationEvent(double eventTime, arma::umat assocMat);
        arma::umat    symbiontTreeEvent(double eventTime, arma::umat assocMat);
        Rcpp::DataFrame createEventDF();
        void      updateEventIndices();
        void      updateEventVector(int h, int s, int e, double time);
        void    clearEventDFVecs();
        void    initializeEventVector();
        Rcpp::CharacterVector  getExtantHostNames(std::vector<std::string> hostNames);
        Rcpp::CharacterVector  getExtantSymbNames(std::vector<std::string> symbNames);
        // anagenetic functions
        double    getTimeToAnaEvent(double dispRate, double extRate, arma::umat assocMat);
        arma::umat symbiontDispersalEvent(int symbInd, arma::umat assocMat);
        arma::umat symbiontExtirpationEvent(int symbInd, arma::umat assocMat);
        arma::umat anageneticEvent(double dispersalRate, double extirpationRate, double currTime, arma::umat assocMat);

};

extern Rcpp::List bdsim_species_tree(double sbr,
                                     double sdr,
                                     int numbsim,
                                     int n_tips,
                                     int gsa_stop);

extern Rcpp::List sim_bdsimple_species_tree(double sbr,
                                            double sdr,
                                            int numbsim,
                                            double timeToSimTo);

extern Rcpp::List sim_locus_tree(std::shared_ptr<SpeciesTree> species_tree,
                                 double gbr,
                                 double gdr,
                                 double lgtr,
                                 int numLoci,
                                 std::string trans_type);

extern Rcpp::List sim_host_symb_treepair(double hostbr,
                                         double hostdr,
                                         double symbbr,
                                         double symbdr,
                                         double switchrate,
                                         double cosprate,
                                         double timeToSimTo,
                                         int host_limit,
                                         int numbsim,
                                         bool hsMode);

extern Rcpp::List sim_host_symb_treepair_ana(double hostbr,
                                            double hostdr,
                                            double symbbr,
                                            double symbdr,
                                            double symbdispersal,
                                            double symbextirpation,
                                            double switchrate,
                                            double cosprate,
                                            double timeToSimTo,
                                            int host_limit,
                                            int numbsim,
                                            bool hsMode);

extern Rcpp::List sim_locus_tree_gene_tree(std::shared_ptr<SpeciesTree> species_tree,
                                           double gbr,
                                           double gdr,
                                           double lgtr,
                                           int numLoci,
                                           double popsize,
                                           int samples_per_lineage,
                                           int numGenesPerLocus);

extern Rcpp::List sim_genetree_msc(std::shared_ptr<SpeciesTree> species_tree,
                                   double popsize,
                                   int samples_per_lineage,
                                   int numbsim);

#endif /* Simulator_h */
