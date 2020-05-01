#include "SpeciesTree.h"
#include <iostream>

using namespace Rcpp;


SpeciesTree::SpeciesTree(unsigned numTaxa, double ct, double br, double dr) : Tree(numTaxa, 0.0){
    currentTime = 0.0;
    extantStop = numTaxa;
    speciationRate = br;
    extinctionRate = dr;

}

SpeciesTree::SpeciesTree(unsigned numTaxa) : Tree(numTaxa){
    extantStop = numTaxa;
}

SpeciesTree::SpeciesTree(SEXP rtree) : Tree(rtree){
  speciationRate = 0.0;
  extinctionRate = 0.0;
}

SpeciesTree::SpeciesTree(const SpeciesTree& speciestree, unsigned numTaxa) : Tree(numTaxa) {
  extantStop = numTaxa;
  nodes = speciestree.nodes;
  extantNodes = speciestree.extantNodes;
  root = speciestree.root;
  speciationRate = speciestree.speciationRate;
  extinctionRate = speciestree.extinctionRate;
  extantStop = speciestree.extantStop;
  extantRoot = speciestree.extantRoot;
  currentTime = speciestree.currentTime;
  numNodes = speciestree.numNodes;
  numTotalTips = speciestree.numTotalTips;
  numExtant = speciestree.numExtant;
  numExtinct = speciestree.numExtinct;
}


SpeciesTree::~SpeciesTree(){

}

double SpeciesTree::getTimeToNextEvent(){
    double sumrt = speciationRate + extinctionRate;
    double returnTime = 0.0;
    returnTime = -log(unif_rand()) / (double(numExtant) * sumrt);
    return returnTime;
}

void SpeciesTree::lineageBirthEvent(unsigned indx){
    Node *sis, *right;
    right = new Node();
    sis = new Node();
    setNewLineageInfo(indx, right, sis);
}

void SpeciesTree::lineageDeathEvent(unsigned int indx){
    extantNodes[indx]->setDeathTime(currentTime);
    extantNodes[indx]->setIsExtant(false);
    extantNodes[indx]->setIsTip(true);
    extantNodes[indx]->setIsExtinct(true);
    extantNodes.erase(extantNodes.begin() + indx);
    numExtinct += 1;
    numExtant = (int) extantNodes.size();
}

void SpeciesTree::ermEvent(double cTime){
    currentTime = cTime;
    int nodeInd = unif_rand()*(numExtant - 1);
    double relBr = speciationRate / (speciationRate + extinctionRate);
    bool isBirth = (unif_rand() < relBr ? true : false);
    if(isBirth)
        lineageBirthEvent(nodeInd);
    else
        lineageDeathEvent(nodeInd);
}

void SpeciesTree::setNewLineageInfo(unsigned int indx, Node *r, Node *l){
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

    l->setLdes(NULL);
    l->setRdes(NULL);
    l->setSib(r);
    l->setAnc(extantNodes[indx]);
    l->setBirthTime(currentTime);
    l->setIsTip(true);
    l->setIsExtinct(false);
    l->setIsExtant(true);

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

void SpeciesTree::setBranchLengths(){
    double bl;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
      bl = (*it)->getDeathTime() - (*it)->getBirthTime();
      branchLengths.push_back(std::move(bl));
      (*it)->setBranchLength(bl);
    }
}

void SpeciesTree::setPresentTime(double currentT){
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        (*it)->setDeathTime(currentT);
        (*it)->setIsExtant(true);
    }
    this->setBranchLengths();
    this->setTreeTipNames();
}

void SpeciesTree::setTreeInfo(){
  //  double trDepth = this->getTreeDepth();
    std::set<double> deathTimes;
    std::vector<Node*>::iterator it = nodes.begin();
    (*it)->setBirthTime(0.0);
    (*it)->setDeathTime((*it)->getBranchLength() + (*it)->getBirthTime());
    (*it)->setIndx(0);
    ++it;
    for(; it != nodes.end(); ++it){
        (*it)->setBirthTime((*it)->getAnc()->getDeathTime());
        (*it)->setDeathTime((*it)->getBranchLength() + (*it)->getBirthTime());
        deathTimes.insert(deathTimes.begin(),(*it)->getBranchLength() + (*it)->getBirthTime());
        (*it)->setIndx((int)std::distance(nodes.begin(), it));
    }
    it = nodes.begin();
    std::set<double>::iterator set_iter = deathTimes.end();
    --set_iter;
    double currentTime = *(set_iter);
    for(; it != nodes.end(); ++it){
        if((*it)->getIsTip()){
            if(std::abs((*it)->getDeathTime() - currentTime) < 0.1){
                (*it)->setIsExtant(true);
                (*it)->setIsExtinct(false);
                (*it)->setDeathTime(currentTime);
                numTaxa++;
                extantNodes.push_back(std::move(*it));
            }
            else{
                (*it)->setIsExtant(false);
                (*it)->setIsExtinct(true);
            }
        }
    }
    return;
}

void SpeciesTree::setTreeTipNames(){
  unsigned nodeIndx = numExtant + numExtinct;
  unsigned tipIt = 0;
  std::stringstream tn;

  for(int i=0; i < nodes.size(); i++){
    if(nodes[i]->getIsTip()){
      tipIt++;
      nodes[i]->setIndx(tipIt);
      if(nodes[i]->getIsExtant()){
        tn << nodes[i]->getIndex();
        std::string name = "H" + tn.str();
        nodes[i]->setName(name);

      }
      else{
        tn << nodes[i]->getIndex();
        std::string name = "X" + tn.str();
        nodes[i]->setName(name);
      }
    }
    else{
      nodeIndx++;
      nodes[i]->setIndx(nodeIndx);
    }
    tn.clear();
    tn.str(std::string());
  }
}


void SpeciesTree::recTipNamer(Node *p, unsigned &nodeIndx, unsigned &tipIndx){
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
        std::string name = "H" + tn.str();
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

void SpeciesTree::recGetNewickTree(Node *p, std::stringstream &ss){
    if(p != NULL){
        if( p->getRdes() == NULL)
            ss <<  p->getName();
        else{
            ss << "(";
            recGetNewickTree(p->getLdes(), ss);
            ss << "[&index=" << p->getLdes()->getIndex() << "]" << ":" << p->getLdes()->getBranchLength();
            ss << ",";
            recGetNewickTree(p->getRdes(), ss);
            ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
            ss << ")";        }
    }
}



std::string SpeciesTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string spTree = ss.str();
    return spTree;
}

std::string SpeciesTree::printExtNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string spTree = ss.str();
    return spTree;
}


void SpeciesTree::setGSATipTreeFlags(){
    zeroAllFlags();
    numTotalTips = 0;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        if((*it)->getIsTip()){
            numTotalTips++;
            (*it)->setFlag(1);

        }
        else{
            (*it)->setFlag(2);
        }
    }
    setSampleFromFlags();
}


//void SpeciesTree::setSampleFromFlags(){
//    int flag;
//    Node *q = NULL;
//    for(std::vector<Node*>::iterator p=nodes.begin(); p!=nodes.end(); p++){
//        if((*p)->getIsTip()){
//            flag = (*p)->getFlag();
//            q = (*p);
//            if(flag == 1){
//                while(q->getIsRoot() == false && flag < 2){
//                    q = q->getAnc();
//                    flag = q->getFlag();
//                    flag++;
//                    q->setFlag(flag);
//
//                }
//            }
//        }
//    }
//}
void SpeciesTree::popNodes(){
    nodes.clear();
    extantNodes.clear();

    recPopNodes(this->getRoot());

   // int indx;
    // for(std::vector<Node*>::iterator p=nodes.begin(); p!=nodes.end(); p++){
    //     indx = (int) (p - nodes.begin());
    //     (*p)->setIndx(indx);
    // }
}

void SpeciesTree::recPopNodes(Node *p){
    if(p != nullptr){
        if(p->getIsTip()){
            if(p->getIsExtant()){
                extantNodes.push_back(std::move(p));
                nodes.push_back(std::move(p));
            }
            else{
                nodes.push_back(std::move(p));
            }
        }
        else{
            nodes.push_back(std::move(p));
            recPopNodes(p->getLdes());
            recPopNodes(p->getRdes());
        }
    }
}

void SpeciesTree::reconstructTreeFromGSASim(Node *oRoot){
    Node *n = new Node();
    unsigned tipCounter = 0;
    unsigned intNodeCounter = extantStop;
    reconstructLineageFromGSASim(n, oRoot, tipCounter, intNodeCounter);
    delete n;
}

void SpeciesTree::reconstructLineageFromGSASim(Node *currN, Node *prevN, unsigned &tipCounter, unsigned &intNodeCounter){
    Node *p;
    bool rootN = prevN->getIsRoot();
    double brlen = prevN->getBranchLength();
    int oFlag = prevN->getFlag();
    if(prevN->getIsTip() && oFlag == 1){
        // need to recalculate branchlength
        Node *prevAnc = prevN->getAnc();
        int ancFlag = prevAnc->getFlag();
        if(ancFlag == 1){
            brlen += prevAnc->getBranchLength();
            while(!prevAnc->getIsRoot() && ancFlag < 2){
                prevAnc = prevAnc->getAnc();
                ancFlag = prevAnc->getFlag();
                if(ancFlag == 1)
                    brlen += prevAnc->getBranchLength();
            }
        }

        p = new Node();
        tipCounter++;
        p->setIndx(tipCounter);
        p->setBranchLength(brlen);
        p->setIsTip(true);
        p->setBirthTime(prevN->getBirthTime());
        p->setDeathTime(prevN->getDeathTime());
        p->setIsExtant(prevN->getIsExtant());
        p->setIsExtinct(prevN->getIsExtinct());
        p->setAnc(currN);
        if(currN->getLdes() == NULL)
            currN->setLdes(p);
        else if(currN->getRdes() == NULL)
            currN->setRdes(p);
        else{
            stop("ERROR: Problem adding a tip to the tree!");
        }

    }
    else{
        if(oFlag > 1){
            Node *s1 = new Node();
            intNodeCounter++;
            s1->setIndx(intNodeCounter);
            if(prevN->getLdes()->getFlag() > 0)
                reconstructLineageFromGSASim(s1, prevN->getLdes(), tipCounter, intNodeCounter);
            if(prevN->getRdes()->getFlag() > 0)
                reconstructLineageFromGSASim(s1, prevN->getRdes(), tipCounter, intNodeCounter);


            if(rootN == false){
                Node *prevAnc = prevN->getAnc();
                int ancFlag = prevAnc->getFlag();
                if(ancFlag == 1){
                    brlen += prevAnc->getBranchLength();
                    while(!(prevAnc)->getIsRoot() && ancFlag < 2){
                        prevAnc = prevAnc->getAnc();
                        ancFlag = prevAnc->getFlag();
                        if(ancFlag == 1)
                            brlen += prevAnc->getBranchLength();
                    }
                }

                if(currN != NULL){
                    s1->setBranchLength(brlen);
                    s1->setBirthTime(prevN->getBirthTime());
                    s1->setDeathTime(prevN->getDeathTime());
                    s1->setAnc(currN);
                    if(currN->getLdes() == NULL)
                        currN->setLdes(s1);
                    else if(currN->getRdes() == NULL)
                        currN->setRdes(s1);
                    else{
                        stop("ERROR: Probem adding an internal node to the tree");
                    }
                }
                else{
                    s1->setAsRoot(true);
                    setRoot(s1);
                    s1->setBranchLength(brlen);
                    s1->setBirthTime(prevN->getBirthTime());
                    s1->setDeathTime(prevN->getDeathTime());
                }

            }
            else{
                s1->setAsRoot(true);
                setRoot(s1);
                s1->setBranchLength(0.0);
                s1->setBirthTime(prevN->getBirthTime());
                s1->setDeathTime(prevN->getDeathTime());
            }

        }
        else if(oFlag == 1){
            if(prevN->getRdes()->getFlag() == 0 && prevN->getLdes()->getFlag() > 0)
                reconstructLineageFromGSASim(currN, prevN->getLdes(), tipCounter, intNodeCounter);
            else
                reconstructLineageFromGSASim(currN, prevN->getRdes(), tipCounter, intNodeCounter);
        }
    }
}


std::map<int,int> SpeciesTree::makeIndxMap(){
  std::map<int,int> indxMap;
  for(int i=0; i < nodes.size(); i++){
    int rIndx = nodes[i]->getIndex();
    int tdckenIndx = i;
    indxMap.insert(std::pair<int,int>(rIndx, tdckenIndx));
  }
  return indxMap;
}


std::map<int,double> SpeciesTree::getBirthTimesFromNodes(){
    int indx;
    double birthTime;
    std::map<int,double> birthTimeMap;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        indx = (*it)->getIndex();
        birthTime = (*it)->getBirthTime();
        birthTimeMap.insert(std::pair<int,double>(indx, birthTime));
    }
    return birthTimeMap;
}

std::map<int,double> SpeciesTree::getDeathTimesFromNodes(){
    int indx;
    double deathTime;
    std::map<int,double> deathTimeMap;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if(!((*it)->getIsExtant())){
            indx = (*it)->getIndex();
            deathTime = (*it)->getDeathTime();
            Rcout << "* " << deathTime << std::endl;

            deathTimeMap.insert(std::pair<int,double>(indx, deathTime));
        }
    }
    return deathTimeMap;
}

std::pair<int,int> SpeciesTree::preorderTraversalStep(int indx){
    std::pair<int,int> sibs;
    sibs.first = nodes[indx]->getLdes()->getIndex();
    sibs.second = nodes[indx]->getRdes()->getIndex();
    return sibs;
}

int SpeciesTree::postOrderTraversalStep(int index){
    int d;
    d = nodes[index]->getAnc()->getIndex();
    return d;
}

bool SpeciesTree::macroEvent(int indx){
    bool isSpec;
    Node* n = nodes[indx];

    Rcout << n->getIsTip() << ":SP:" << n->getIndex() << std::endl;
    if(n->getIsTip())
        isSpec = false;
    else
        isSpec = true;
    return isSpec;
}


int SpeciesTree::findLastToGoExtinct(double EventTime){
  int indxExtinct = -1;
  double epsi = std::numeric_limits<double>::epsilon();
  bool is_near = false;
  for(int i=0; i < nodes.size(); i++){
    if(nodes[i]->getIsTip() && nodes[i]->getIsExtinct()){
      double scale = std::max(abs(EventTime), abs(nodes[i]->getDeathTime()));
      is_near = abs(nodes[i]->getDeathTime() - EventTime) <= scale *(2*epsi);
      if(is_near){
        indxExtinct = i;
        break;
      }
    }
  }

  return indxExtinct;
}