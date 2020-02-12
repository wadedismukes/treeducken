//
// Created by Dismukes, Wade T [EEOBS] on 10/31/19.
//

#include "SymbiontTree.h"
SymbiontTree::SymbiontTree(int numTaxa,
                           double ct,
                           double br,
                           double dr,
                           double hostExpanRate,
                           int K) : Tree(numTaxa, 0.0){
    currentTime = 0.0;
    numTaxa = numTaxa;
    symbSpecRate = br;
    symbExtRate = dr;
    hostExpanRate = hostExpanRate;
    numExpansions = 0;
    hostLimit = K;

}


SymbiontTree::~SymbiontTree(){

}

double SymbiontTree::getTimeToNextEvent(double hostSpecRate,
                                        double hostExtRate,
                                        double cospeciaRate,
                                        int numHosts){
    double sumrt_host =  (hostSpecRate + hostExtRate) * numHosts;
    double sumrt_symb = (symbSpecRate + symbExtRate + hostExpanRate) * numExtant;
    double sumrt_both = cospeciaRate * this->getNumHostSymbPairs();
    double returnTime = -log(unif_rand()) / (sumrt_host + sumrt_symb + sumrt_both);
    return returnTime;
}

void SymbiontTree::setSymbTreeInfoSpeciation(int indxToFind, int indxToReplace){
    for(std::vector<Node*>::iterator s = extantNodes.begin(); s != extantNodes.end(); ++s){
        std::vector<int> hostsOfS = (*s)->getHosts();
        for(std::vector<int>::iterator it=hostsOfS.begin(); it != hostsOfS.end(); ++it){
            if((*it) == indxToFind)
                (*it) = indxToReplace;
        }
    }
}

void SymbiontTree::setSymbTreeInfoExtinction(int deadIndx){
    std::vector<int> toBeExtincted;
    for(int i = 0; i < extantNodes.size(); i++){
        std::vector<int> hostsOf = extantNodes[i]->getHosts();
        for(int j = 0; j < hostsOf.size(); j++){
            if(hostsOf[j] == deadIndx){
                std::swap(hostsOf.back(),j);
                hostsOf.pop_back();
            }
        }
        if(hostsOf.empty())
            toBeExtincted.push_back(std::move(i));
    }
    if(!(toBeExtincted.empty())){
        for(int i = 0; i < toBeExtincted.size(); i++){
            this->lineageDeathEvent(toBeExtincted[i]);
        }
    }
}