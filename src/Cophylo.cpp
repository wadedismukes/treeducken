#include <Rcpp.h>
#include "Engine.h"
#include "SymbiontTree.h"
using namespace Rcpp;


Rcpp::List sim_host_symb_treepair(double hostbr,
                                  double hostdr,
                                  double symbbr,
                                  double symbdr,
                                  double switchRate,
                                  double cospeciationRate,
                                  double timeToSimTo,
                                  int numbsim){
    std::vector<SpeciesTree*> hostTrees;
    hostTrees.resize(numbsim);

    std::vector<SymbiontTree*> symbTrees;
    symbTrees.resize(numbsim);
    double rho = 1.0;
    Simulator *phySimulator = nullptr;
    Rcpp::List multiphy;
    Rcpp::List hostSymbPair;
    unsigned numLociToSim = numbsim;
    for(int i = 0; i < numbsim; i++){
        // TODO: write new constructor for simulator
        Simulator *phySimulator = new Simulator( timeToSimTo,
                                                 hostbr,
                                                 hostdr,
                                                 symbbr,
                                                 symbdr,
                                                 switchRate,
                                                 cospeciationRate,
                                                 rho);
        phySimulator->simHostSymbSpeciesTreePair();
        symbTrees[i] = phySimulator->getSymbiontTree();
        hostTrees[i] = phySimulator->getSpeciesTree();





        NumericMatrix edges_rmat = symbTrees[i]->getEdges();

        List phySymb = List::create(Named("edge") = edges_rmat,
                                Named("edge.length") = symbTrees[i]->getEdgeLengths(),
                                Named("Nnode") = symbTrees[i]->getNnodes(),
                                Named("tip.label") = symbTrees[i]->getTipNames(),
                                Named("root.edge") = symbTrees[i]->getRoot()->getDeathTime() -
                                    symbTrees[i]->getRoot()->getBirthTime());
        phySymb.attr("class") = "phylo";

        edges_rmat = hostTrees[i]->getEdges();

        List phyHost = List::create(Named("edge") = edges_rmat,
                                    Named("edge.length") = hostTrees[i]->getEdgeLengths(),
                                    Named("Nnode") = hostTrees[i]->getNnodes(),
                                    Named("tip.label") = hostTrees[i]->getTipNames(),
                                    Named("root.edge") = hostTrees[i]->getRoot()->getDeathTime() -
                                        hostTrees[i]->getRoot()->getBirthTime());
        phyHost.attr("class") = "phylo";

        hostSymbPair = List::create(Named("host_tree") = phyHost,
                                    Named("symb_tree") = phySymb);

        multiphy[i] = hostSymbPair;
    }



    delete phySimulator;
    return multiphy;
}