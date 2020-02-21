//
// Created by Dismukes, Wade T [EEOBS] on 10/31/19.
//

#include "SymbiontTree.h"
SymbiontTree::SymbiontTree(int nt,
                           double ct,
                           double br,
                           double dr,
                           double her,
                           int K) : Tree(nt, 0.0){
    currentTime = 0.0;
    numTaxa = nt;
    symbSpecRate = br;
    symbExtRate = dr;
    hostExpanRate = her;
    numExpansions = 0;
    hostLimit = K;
    root->addHost(0);
    std::vector<int> initialHosts;
    initialHosts.resize(1);
    initialHosts[0] = 0;
    symbHostMap.insert(std::pair<int,std::vector<int>> (0,initialHosts));
}

SymbiontTree::SymbiontTree(const SymbiontTree& symbionttree, unsigned numTaxa) : Tree(numTaxa) {
    nodes = symbionttree.nodes;
    extantNodes = symbionttree.extantNodes;
    root = symbionttree.root;
    symbSpecRate = symbionttree.symbSpecRate;
    symbExtRate = symbionttree.symbExtRate;
    extantRoot = symbionttree.extantRoot;
    currentTime = symbionttree.currentTime;
    numNodes = symbionttree.numNodes;
    numTotalTips = symbionttree.numTotalTips;
    numExtant = symbionttree.numExtant;
    numExtinct = symbionttree.numExtinct;
}

SymbiontTree::~SymbiontTree(){

}


double SymbiontTree::getTimeToNextJointEvent(double hostSpecRate,
                                        double hostExtRate,
                                        double cospeciaRate,
                                        int numHosts){
    double sumrt_host =  (hostSpecRate + hostExtRate) * numHosts;
    double sumrt_symb = (symbSpecRate + symbExtRate + hostExpanRate) * numExtant;
    double sumrt_both = cospeciaRate * this->getNumHostSymbPairs();
    // TODO: make map for host->symbionts
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

std::vector<int> SymbiontTree::getSymbsOnHost(int hostIndx){
    std::vector<int> symbs = symbHostMap[hostIndx];

    return symbs;
}

void SymbiontTree::lineageBirthEvent(unsigned indx){
    Node *sis, *right;
    right = new Node();
    sis = new Node();
    setNewLineageInfo(indx, right, sis);
}

void SymbiontTree::lineageDeathEvent(unsigned indx){
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(true);
    extantNodes[indx]->setIsExtinct(true);
    extantNodes.erase(extantNodes.begin() + indx);
    numExtinct += 1;
    numExtant = (int) extantNodes.size();
}

void SymbiontTree::setNewLineageInfo(int indx, Node*r, Node*l){
    extantNodes[indx]->setLdes(l);
    extantNodes[indx]->setRdes(r);
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsTip(false);
    extantNodes[indx]->setIsExtant(false);

    r->setLdes(NULL);
    r->setRdes(NULL);
    r->setSib(l);
    r->setAnc(extantNodes[indx]);
    r->setBirthTime(currentTime);
    r->setIsTip(true);
    r->setIsExtant(true);
    r->setIsExtinct(false);
    //r->setHosts(extantNodes[indx]->getHosts());
    //std::vector<int> hostCheckVec = r->getHosts();

    l->setLdes(NULL);
    l->setRdes(NULL);
    l->setSib(r);
    l->setAnc(extantNodes[indx]);
    l->setBirthTime(currentTime);
    l->setIsTip(true);
    l->setIsExtinct(false);
    l->setIsExtant(true);
    //l->setHosts(extantNodes[indx]->getHosts());

    extantNodes.erase(extantNodes.begin() + indx);
    extantNodes.push_back(std::move(r));
    extantNodes.push_back(std::move(l));
    nodes.push_back(std::move(r));
    nodes.push_back(std::move(l));
    numNodes = (int) nodes.size();
    numExtant = (int) extantNodes.size();
    r->setIndx(numNodes - 2);
    l->setIndx(numNodes - 1);
}

arma::mat SymbiontTree::ermJointEvent(double ct, arma::mat assocMat){
    currentTime = ct;
    this->setCurrentTime(ct);

    // pick a row at random
    int nodeInd = unif_rand()*(numExtant - 1);

    arma::rowvec rvec = assocMat.row(nodeInd);
    assocMat.shed_row(nodeInd);

    // which event
    double relBr = symbSpecRate / (symbExtRate + symbSpecRate + hostExpanRate);
    double relDr = relBr + (symbExtRate / (symbExtRate + symbSpecRate + hostExpanRate));
    double dec = unif_rand();
    if(dec < relBr){
        // its a birth
        this->lineageBirthEvent(nodeInd);
        assocMat.resize(numExtant, assocMat.n_cols);
        assocMat(numExtant - 2, arma::span::all) = rvec;
        assocMat(numExtant - 1, arma::span::all) = rvec;
    }
    else if(dec < relDr)
        this->lineageDeathEvent(nodeInd);

    else{
        int hostInd = unif_rand() * assocMat.n_cols;
        this->hostExpansionEvent(nodeInd, hostInd);
        assocMat.resize(numExtant, assocMat.n_cols);
        assocMat(numExtant - 2, arma::span::all) = rvec;
        rvec(hostInd) = 1;
        assocMat(numExtant - 1, arma::span::all) = rvec;
    }
    return assocMat;
}

void SymbiontTree::hostExpansionEvent(int indx, int hostIndx){
    Node *sis, *right;
    right = new Node();
    sis = new Node();
    this->setNewLineageInfoExpan(indx, right, sis, hostIndx);
}

void SymbiontTree::setNewLineageInfoExpan(int indx, Node* r, Node* l, int hostIndx){
    extantNodes[indx]->setLdes(l);
    extantNodes[indx]->setRdes(r);
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsTip(false);
    extantNodes[indx]->setIsExtant(false);

    r->setLdes(NULL);
    r->setRdes(NULL);
    r->setSib(l);
    r->setAnc(extantNodes[indx]);
    r->setBirthTime(currentTime);
    r->setIsTip(true);
    r->setIsExtant(true);
    r->setIsExtinct(false);
    //r->setHosts(extantNodes[indx]->getHosts());

    l->setLdes(NULL);
    l->setRdes(NULL);
    l->setSib(r);
    l->setAnc(extantNodes[indx]);
    l->setBirthTime(currentTime);
    l->setIsTip(true);
    l->setIsExtinct(false);
    l->setIsExtant(true);
    //l->setHosts(extantNodes[indx]->getHosts());
    // if(unif_rand() < 0.5){
    //     l->addHost(hostIndx);
    // }
    // else{
    //     r->addHost(hostIndx);
    // }
    extantNodes.erase(extantNodes.begin() + indx);
    extantNodes.push_back(std::move(r));
    extantNodes.push_back(std::move(l));
    nodes.push_back(std::move(r));
    nodes.push_back(std::move(l));
    numExtant = (int) extantNodes.size();
    r->setIndx(numExtant - 2);
    l->setIndx(numExtant - 1);
}

void SymbiontTree::setTreeTipNames(){
    unsigned nodeIndx = numExtant + numExtinct;
    unsigned tipIt = 0;
    recTipNamer(this->getRoot(), nodeIndx, tipIt);
}

void SymbiontTree::recTipNamer(Node *p, unsigned &nodeIndx, unsigned &tipIndx){
    if(p != NULL){
        std::stringstream tn;
        if(p->getIsTip()){
            tipIndx++;
            p->setIndx(tipIndx);
            if(p->getIsExtinct()){
                tn << p->getIndex();
                std::string name = "X" + tn.str();
                p->setName(name);

            }
            else{
                tn << p->getIndex();
                std::string name = "S" + tn.str();
                p->setName(name);
            }
        }
        else{
            nodeIndx++;
            p->setIndx(nodeIndx);
            recTipNamer(p->getLdes(), nodeIndx, tipIndx);
            recTipNamer(p->getRdes(), nodeIndx, tipIndx);

        }
    }
}

void SymbiontTree::setBranchLengths(){
    double bl;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        bl = (*it)->getDeathTime() - (*it)->getBirthTime();
        branchLengths.push_back(std::move(bl));
        (*it)->setBranchLength(bl);
    }
}

void SymbiontTree::setPresentTime(double currentT){
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        (*it)->setDeathTime(currentT);
        (*it)->setIsExtant(true);
    }
    this->setBranchLengths();
    this->setTreeTipNames();
}


void SymbiontTree::updateCurrentMap(int oldHostIndx, int newHostIndx){
    std::map<int, std::vector<int>>::iterator found;
    found = symbHostMap.find(oldHostIndx);
    if (found != symbHostMap.end()) {
        // Swap value from oldKey to newKey, note that a default constructed value
        // is created by operator[] if 'm' does not contain newKey.
        std::swap(symbHostMap[newHostIndx], found->second);
        // Erase old key-value from map
        symbHostMap.erase(found);
    }

}



void SymbiontTree::cospeciationMapUpdate(int oldHostIndx,
                                         int numNodesHost,
                                         int oldSymbIndx){

    std::vector<int> symbsOnHost = symbHostMap[oldHostIndx];
    std::vector<int> leftHostSymbiontsValues;
    std::vector<int> rightHostSymbiontsValues;
    std::vector<Node*> nodesForUpdating = this->getNodes();
    for(int i = 0; i < symbsOnHost.size(); i++){
        std::vector<int> hostsInSymb = nodesForUpdating[symbsOnHost[i]]->getHosts();
        if(oldSymbIndx == symbsOnHost[i]){
            leftHostSymbiontsValues.push_back(this->getNodesSize() - 1);
            rightHostSymbiontsValues.push_back(this->getNodesSize() - 2);
        }
        else{
            double which = unif_rand();
            if(which < 0.5){
                leftHostSymbiontsValues.push_back(symbsOnHost[i]);
                for(int i=0; i < hostsInSymb.size(); i++){
                    if(hostsInSymb[i] == oldHostIndx)
                        hostsInSymb[i] = numNodesHost - 1;
                }
            }
            else{
                rightHostSymbiontsValues.push_back(symbsOnHost[i]);
                for(int i=0; i < hostsInSymb.size(); i++){
                    if(hostsInSymb[i] == oldHostIndx)
                        hostsInSymb[i] = numNodesHost - 2;
                }
            }

        }
        nodesForUpdating[symbsOnHost[i]]->setHosts(hostsInSymb);
    }
    symbHostMap[numNodesHost - 2] = rightHostSymbiontsValues;
    symbHostMap[numNodesHost - 1] = leftHostSymbiontsValues;
    // need to update nodes[i].hosts
    std::map<int, std::vector<int>>::iterator it = symbHostMap.find(oldHostIndx);
    symbHostMap.erase(it);


}

void SymbiontTree::updateHostsInNodes(){
    for(int i=0; i < nodes.size(); i++){

    }
}

int SymbiontTree::getExtantIndxFromNodes(int nodesIndx){
    int count = 0;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it ){
        if((*it)->getIndex() == nodesIndx)
            break;
        count++;
    }
    return count;
}