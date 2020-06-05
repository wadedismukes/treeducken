//
//  GeneTree.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 12/20/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "GeneTree.h"
#include <iostream>
#include <cmath>
#include <Rcpp.h>

GeneTree::GeneTree(unsigned nt, unsigned ipp, double ne, double genTime) : Tree(nt){
    numTaxa = nt;
    individualsPerPop = ipp;
    popSize = ne;
    generationTime = genTime;
    delete root;
}

GeneTree::~GeneTree(){
    // for(std::vector<Node*>::iterator p=nodes.begin(); p != nodes.end(); ++p){
    //     if((*p) != nullptr){
    //         delete (*p);
    //         (*p) = nullptr;
    //     }
    // }
    //wclearNodes(extantRoot);
}


//TODO:  go back and speed this up by removing push_back calls
void GeneTree::initializeTree(std::vector< std::vector<int> > extantLociInd, double presentTime){
    for(std::vector<Node*>::iterator p=nodes.begin(); p != nodes.end(); ++p){
        delete (*p);
        (*p) = nullptr;
    }
    nodes.clear();
    for(std::vector<Node*>::iterator p=extantNodes.begin(); p != extantNodes.end(); ++p){
        delete (*p);
        (*p) = nullptr;
    }
    extantNodes.clear();
    Node *p;
    int numLociInPresnt = extantLociInd[0].size();
    for(int i = 0; i < numLociInPresnt; i++){
        for(int j = 0; j < individualsPerPop; j++){
            p = new Node();
            p->setDeathTime(presentTime);
            p->setLindx(extantLociInd[0][i]);
            p->setLdes(NULL);
            p->setRdes(NULL);
            p->setAnc(NULL);
            p->setIsExtant(true);
            p->setIsTip(true);
            p->setIsExtinct(false);
            extantNodes.push_back(p);
            nodes.push_back(p);
            p->setIndx((int) nodes.size());
        }
    }
}

double GeneTree::getCoalTime(int n){
    double ct;
    double lambda = (double)(n * (n - 1)) / (popSize) ;
    ct = -log(unif_rand()) / (lambda);
    return ct;
}

bool GeneTree::censorCoalescentProcess(double startTime, double stopTime, int contempSpeciesIndx, int ancSpIndx, bool chck){
    int leftInd, rightInd;
    int leftIndExtN, rightIndExtN;
    int extIndx;
    Node *l, *r;
    Node *n;
    double t = startTime;
    bool allCoalesced = false;
    // search extantNodes for members with Lindx = contempSpecisIndx
    std::vector<int> indInExtNodes;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        if((*it)->getLindx() == contempSpeciesIndx){
            extIndx = std::distance(extantNodes.begin(), it);
            indInExtNodes.push_back(extIndx);
        }
    }
  //  std::cout << contempSpeciesIndx << "   ($)$)%    " << indInExtNodes.size() << std::endl;
    if(indInExtNodes.size() > 1){
        while(t > stopTime){
            t -= getCoalTime(indInExtNodes.size());
            // std::cout << t << std::endl;
            if(t < stopTime){
                if(chck){
                    t = stopTime;
                    allCoalesced = true;
                    break;
                }
                else{
                    t = stopTime;
                    allCoalesced = false;
                    break;
                }
            }
            rightInd = unif_rand() * (indInExtNodes.size() - 1);
            iter_swap(indInExtNodes.begin() + rightInd, indInExtNodes.begin());
            rightIndExtN = indInExtNodes[0];
            r = extantNodes[rightIndExtN];

            std::reverse(indInExtNodes.begin(), indInExtNodes.end());
            leftInd = unif_rand() * (indInExtNodes.size() - 2);
            iter_swap(indInExtNodes.begin() + leftInd, indInExtNodes.begin());
            leftIndExtN = indInExtNodes[0];
            l = extantNodes[leftIndExtN];
            std::reverse(indInExtNodes.begin(), indInExtNodes.end());

            n = coalescentEvent(t, l, r);
            //iter_swap(extantNodes.begin() + rightIndExtN, extantNodes.end() - 1);
            //iter_swap(extantNodes.begin() + leftIndExtN, extantNodes.end() - 2);
            if(leftIndExtN > rightIndExtN){
                extantNodes.erase(extantNodes.begin() + leftIndExtN);
                extantNodes.erase(extantNodes.begin() + rightIndExtN);
            }
            else{
                extantNodes.erase(extantNodes.begin() + rightIndExtN);
                extantNodes.erase(extantNodes.begin() + leftIndExtN);
            }
            indInExtNodes.clear();
            //indInExtNodes.erase(indInExtNodes.begin(), indInExtNodes.begin() + 2);
            extantNodes.insert(extantNodes.begin(), n);

            // if(!(indInExtNodes.empty())){
            //     for(std::vector<int>::iterator it = indInExtNodes.begin(); it != indInExtNodes.end(); ++it){
            //         (*it) += 1;
            //     }
            // }
            for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
                if((*it)->getLindx() == contempSpeciesIndx){
                    extIndx = std::distance(extantNodes.begin(), it);
                    indInExtNodes.push_back(extIndx);
                }
            }
//            extantNodes.erase(extantNodes.end() - 2, extantNodes.end());
         //  std::cout << "size of extantNodes " << extantNodes.size() << std::endl;
            if(indInExtNodes.size() == 1){
                allCoalesced = true;
                break;
            }
        }
    }
    else if (indInExtNodes.size() == 1){
        t = stopTime;
        allCoalesced = true;
        extantNodes[indInExtNodes[0]]->setLindx(ancSpIndx);
    }
    else{
        allCoalesced = true;
    }

    if(allCoalesced == true){
        for(int i = 0; i < indInExtNodes.size(); ++i){
            extantNodes[indInExtNodes[i]]->setLindx(ancSpIndx);
        }
    }

    indInExtNodes.clear();
    return allCoalesced;
}

Node* GeneTree::coalescentEvent(double t, Node *p, Node *q){
    Node *n = new Node();
    n->setDeathTime(t);
    n->setLdes(p);
    n->setRdes(q);
    n->setIsExtant(false);
    n->setIsTip(false);
    n->setIsExtinct(false);
    n->setLindx(p->getLindx());
    nodes.push_back(n);
    n->setIndx((int) nodes.size());
    p->setBirthTime(t);
    p->setAnc(n);
   // p->setSib(q);

    q->setBirthTime(t);
    q->setAnc(n);
    // q->setSib(p);


    return n;
}

std::multimap<int, double> GeneTree::rescaleTimes(std::multimap<int, double> timeMap){
    std::multimap<int, double> rescaledTimeMap;
    std::pair<int, double> p;
    for(std::multimap<int, double>::iterator it = timeMap.begin(); it != timeMap.end(); ++it){
        p.first = (*it).first;
        p.second = ((*it).second);
        rescaledTimeMap.insert(p);
    }

    return rescaledTimeMap;

}



void GeneTree::rootCoalescentProcess(double startTime){
    int leftInd, rightInd;
    Node *l, *r;
    Node *n;
    double t = startTime;
    // search extantNodes for members with Lindx = contempSpecisIndx
    std::vector<int> indInExtNodes;
    for(std::vector<Node*>::iterator it = extantNodes.begin(); it != extantNodes.end(); ++it){
        (*it)->setLindx(0);
    }
    while(extantNodes.size() > 1){
        t -= getCoalTime(extantNodes.size());

        rightInd = unif_rand() * (extantNodes.size() - 1);
        r = extantNodes[rightInd];
        extantNodes.erase(extantNodes.begin() + rightInd);

        leftInd = unif_rand() * (extantNodes.size() - 1);
        l = extantNodes[leftInd];
        extantNodes.erase(extantNodes.begin() + leftInd);

        n = coalescentEvent(t, l, r);
        extantNodes.push_back(n);
    }
    extantNodes[0]->setAsRoot(true);
    extantNodes[0]->setBirthTime(n->getDeathTime());
    this->setRoot(extantNodes[0]);
    this->setBranchLengths();
}

void GeneTree::recursiveRescaleTimes(Node* r, double add){
    if(r != NULL){
        if( r->getRdes() == NULL){
            r->setBirthTime(r->getBirthTime() + add);
            r->setDeathTime(r->getDeathTime() + add);
        }
        else{

            r->getLdes()->setBirthTime(r->getLdes()->getBirthTime() + add);
            r->getLdes()->setDeathTime(r->getLdes()->getDeathTime() + add);
            recursiveRescaleTimes(r->getLdes(), add);

            r->getRdes()->setBirthTime(r->getRdes()->getBirthTime() + add);
            r->getRdes()->setDeathTime(r->getRdes()->getDeathTime() + add);
            recursiveRescaleTimes(r->getRdes(), add);

        }
    }
}

void GeneTree::setBranchLengths(){
    double brlen;
    numExtant = 0;
    numExtinct = 0;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        brlen = (*it)->getDeathTime() - (*it)->getBirthTime();
        (*it)->setBranchLength(brlen);

        if((*it)->getIsTip()){
          if((*it)->getIsExtant())
            numExtant++;
          else
            numExtinct++;
          branchLengths.push_back(std::move(brlen));
        }
        else if((*it)->getIsRoot()){
          branchLengths.emplace(branchLengths.begin(), brlen);
        }
        else{
          branchLengths.push_back(std::move(brlen));
        }
    }
    this->setTreeTipNames();
}

void GeneTree::addExtinctSpecies(double bt, int indx){
    Node *p;
    for(int i = 0; i < individualsPerPop; i++){
        p = new Node();
        p->setDeathTime(bt);
        p->setLindx(indx);
        p->setLdes(NULL);
        p->setRdes(NULL);
        p->setAnc(NULL);
        p->setIsExtant(false);
        p->setIsTip(true);
        p->setIsExtinct(true);
        extantNodes.push_back(p);
        nodes.push_back(p);
        p->setIndx((int) nodes.size() + 1);

    }
    //delete p;
}


void GeneTree::setIndicesBySpecies(std::map<int, int> spToLocusMap){
    numExtant = 0;
    numExtinct = 0;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if((*it)->getIsTip()){
            // indx = (*it)->getLindx();
            // spIndx = spToLocusMap.find(indx)->second;
            //(*it)->setIndx((*it)->getLindx());

            if((*it)->getIsExtant())
              numExtant++;
            else
              numExtinct++;
        }
    }
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
      int notTipCount = numExtant + numExtinct + 1;
      if(!((*it)->getIsTip())){
        (*it)->setIndx(notTipCount);
        notTipCount++;
      }
    }
}



void GeneTree::setTreeTipNames(){
    int indNumber = 0;
    std::stringstream tn;
    std::string name;
    int locusIndxCounter = 0;

    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        if((*it)->getIsTip()){
            tn << locusIndxCounter + 1;
            name = tn.str();
            tn.clear();
            tn.str(std::string());
            indNumber++;
            tn << indNumber;
            name += "_" + tn.str();
            tn.clear();
            tn.str(std::string());
            (*it)->setName(name);
            if(indNumber == individualsPerPop){
              indNumber = 0;
              locusIndxCounter++;
            }
        }

    }

}

void GeneTree::reindexForR(){
  int intNodeCount = numExtant + numExtinct + 1;
  int tipCount = 1;
  for(int i = nodes.size() - 1; i > -1; i--){
    if(nodes[i]->getIsTip()){
      nodes[i]->setIndx(tipCount);
      tipCount++;
    }
    else{
      nodes[i]->setIndx(intNodeCount);
      intNodeCount++;
    }
  }
}


NumericMatrix GeneTree::getGeneEdges(){
  this->GeneTree::reindexForR();
  int numRows = (int) nodes.size() - 1;
  NumericMatrix edgeMat(numRows, 2);
  for(int i=0; i < nodes.size()-1; i++){
    if(!(nodes[i]->getIsRoot())){

      NumericMatrix::Row row = edgeMat(i, _);

      row[0] = nodes[i]->getAnc()->getIndex();
      row[1] = nodes[i]->getIndex();
    }
  }
  return edgeMat;
}

