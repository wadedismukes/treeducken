#include <RcppArmadillo.h>
#include "SymbiontTree.h"
#include "Simulator.h"
using namespace Rcpp;

Rcpp::List sim_host_symb_treepair_ana(double hostbr,
                                      double hostdr,
                                      double symbbr,
                                      double symbdr,
                                      double symb_dispersal,
                                      double symb_extirpation,
                                      double switchRate,
                                      double cospeciationRate,
                                      double timeToSimTo,
                                      int host_limit,
                                      int numbsim,
                                      std::string hsMode) {
    double rho = 1.0;
    Rcpp::List multiphy;
    Rcpp::List hostSymbPair;
    for(int i = 0; i < numbsim; i++){
        Simulator phySimulator(Simulator::CophySimAna(hostbr, 
                                                hostdr, 
                                                symbbr, 
                                                symbdr, 
                                                switchRate, 
                                                cospeciationRate, 
                                                timeToSimTo, 
                                                host_limit, 
                                                hsMode, 
                                                symb_dispersal, 
                                                symb_extirpation));
        phySimulator.simHostSymbSpeciesTreePairWithAnagenesis();





        List phyHost = List::create(Named("edge") = phySimulator.getSpeciesEdges(),
                                    Named("edge.length") = phySimulator.getSpeciesEdgeLengths(),
                                    Named("Nnode") = phySimulator.getSpeciesNnodes(),
                                    Named("tip.label") = phySimulator.getSpeciesTipNames(),
                                    Named("root.edge") = phySimulator.getSpeciesTreeRootEdge());
        phyHost.attr("class") = "phylo";


        List phySymb = List::create(Named("edge") = phySimulator.getSymbiontEdges(),
                                    Named("edge.length") = phySimulator.getSymbiontEdgeLengths(),
                                    Named("Nnode") = phySimulator.getSymbiontNnodes(),
                                    Named("tip.label") = phySimulator.getSymbiontTipNames(),
                                    Named("root.edge") = phySimulator.getSymbiontTreeRootEdge());
        phySymb.attr("class") = "phylo";
        Rcpp::NumericMatrix assocMat = Rcpp::wrap(phySimulator.getAssociationMatrix());
        assocMat = Rcpp::transpose(assocMat);
        Rcpp::CharacterVector hostNames = phySimulator.getExtantHostNames(phySimulator.getSpeciesTipNames());
        Rcpp::CharacterVector symbNames = phySimulator.getExtantSymbNames(phySimulator.getSymbiontTipNames());
        Rcpp::rownames(assocMat) = hostNames;
        Rcpp::colnames(assocMat) = symbNames;
        hostSymbPair = List::create(Named("host_tree") = phyHost,
                                    Named("symb_tree") = phySymb,
                                    Named("association_mat") = assocMat,
                                    Named("event_history") = phySimulator.createEventDF());
        hostSymbPair.attr("class") = "cophy";
        multiphy.push_back(hostSymbPair);
    }


    multiphy.attr("class") = "multiCophy";
    return multiphy;
}

Rcpp::List sim_host_symb_treepair(double hostbr,
                                  double hostdr,
                                  double symbbr,
                                  double symbdr,
                                  double switchRate,
                                  double cospeciationRate,
                                  double timeToSimTo,
                                  int host_limit,
                                  int numbsim,
                                  std::string hsMode,
                                  bool mutualism){

    double rho = 1.0;
    Rcpp::List multiphy;
    Rcpp::List hostSymbPair;
    for(int i = 0; i < numbsim; i++){
        Simulator phySimulator(Simulator::CophySim(hostbr, 
                                                   hostdr, 
                                                   symbbr, 
                                                   symbdr, 
                                                   switchRate, cospeciationRate,
                                                   timeToSimTo, 
                                                   host_limit, 
                                                   hsMode, 
                                                   mutualism));
        phySimulator.simHostSymbSpeciesTreePair();






        List phyHost = List::create(Named("edge") = phySimulator.getSpeciesEdges(),
                                    Named("edge.length") = phySimulator.getSpeciesEdgeLengths(),
                                    Named("Nnode") = phySimulator.getSpeciesNnodes(),
                                    Named("tip.label") = phySimulator.getSpeciesTipNames(),
                                    Named("root.edge") = phySimulator.getSpeciesTreeRootEdge());
        phyHost.attr("class") = "phylo";


        List phySymb = List::create(Named("edge") = phySimulator.getSymbiontEdges(),
                                Named("edge.length") = phySimulator.getSymbiontEdgeLengths(),
                                Named("Nnode") = phySimulator.getSymbiontNnodes(),
                                Named("tip.label") = phySimulator.getSymbiontTipNames(),
                                Named("root.edge") = phySimulator.getSymbiontTreeRootEdge());
        phySymb.attr("class") = "phylo";
        Rcpp::NumericMatrix assocMat = Rcpp::wrap(phySimulator.getAssociationMatrix());
        assocMat = Rcpp::transpose(assocMat);
        Rcpp::CharacterVector hostNames = phySimulator.getExtantHostNames(phySimulator.getSpeciesTipNames());
        Rcpp::CharacterVector symbNames = phySimulator.getExtantSymbNames(phySimulator.getSymbiontTipNames());
        Rcpp::rownames(assocMat) = hostNames;
        Rcpp::colnames(assocMat) = symbNames;
        hostSymbPair = List::create(Named("host_tree") = phyHost,
                                    Named("symb_tree") = phySymb,
                                    Named("association_mat") = assocMat,
                                    Named("event_history") = phySimulator.createEventDF());
        hostSymbPair.attr("class") = "cophy";
        multiphy.push_back(hostSymbPair);
    }


    multiphy.attr("class") = "multiCophy";

    return multiphy;
}