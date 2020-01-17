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

GeneTree::GeneTree(unsigned nt, unsigned ipp, unsigned ne, double genTime) : Tree(nt){
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
            p->setIndx(extantLociInd[0][i]);
            p->setLdes(NULL);
            p->setRdes(NULL);
            p->setAnc(NULL);
            p->setIsExtant(true);
            p->setIsTip(true);
            p->setIsExtinct(false);
            if(extantLociInd[0][i] == -1){
                this->setOutgroup(p);
                p->setName("OUT");
            }
            else{
                extantNodes.push_back(p);
            }
            nodes.push_back(p);
        }
    }
}

double GeneTree::getCoalTime(int n){
    double ct;
    double lambda = (double)(n * (n - 1)) / (2 * popSize) ;
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

            leftInd = (0.1 + unif_rand()) *(indInExtNodes.size() - 1);
            iter_swap(indInExtNodes.begin() + leftInd, indInExtNodes.begin() + 1);
            leftIndExtN = indInExtNodes[1];
            l = extantNodes[leftIndExtN];

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
        // std::cout << "this shouldn't even be happnening" << std::endl;
        allCoalesced = true;
    }
    
    if(allCoalesced == true){
        for(int i = 0; i < indInExtNodes.size(); ++i){
          //  std::cout << "**************" << std::endl;
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
    n->setIndx(p->getIndex());
    nodes.push_back(n);
    
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



void GeneTree::rootCoalescentProcess(double startTime, double ogf){
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
    if(ogf == 0.0){
        extantNodes[0]->setAsRoot(true);
        this->setRoot(extantNodes[0]);
    }
    else{
        Node *nRoot = new Node(); 
        t -= getCoalTime(2);
        extantNodes[0]->setBirthTime(t);
        nRoot->setBirthTime(t);
        nRoot->setLdes(extantNodes[0]);
        nRoot->setRdes(this->getOutgroup());
        nRoot->setDeathTime(extantNodes[0]->getBirthTime());
        nRoot->setAsRoot(true);
        
        this->setRoot(nRoot);
        this->getOutgroup()->setBirthTime(extantNodes[0]->getBirthTime());
        this->getOutgroup()->setAnc(nRoot);
        extantNodes[0]->setAnc(nRoot);
        nodes.push_back(nRoot);

    }


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
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        brlen = (*it)->getDeathTime() - (*it)->getBirthTime();
        (*it)->setBranchLength(std::abs(brlen));
    }
}

void GeneTree::addExtinctSpecies(double bt, int indx){
    Node *p;
    for(int i = 0; i < individualsPerPop; i++){
        p = new Node();
        p->setDeathTime(bt);
        p->setIndx(indx);
        p->setLindx(indx);
        p->setLdes(NULL);
        p->setRdes(NULL);
        p->setAnc(NULL);
        p->setIsExtant(false);
        p->setIsTip(true);
        p->setIsExtinct(true);
        extantNodes.push_back(p);
        nodes.push_back(p);

    }
    //delete p;
}


void GeneTree::setIndicesBySpecies(std::map<int, int> spToLocusMap){
    int indx;
    int spIndx;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        if((*it)->getIsTip()){
            indx = (*it)->getIndex();
            spIndx = spToLocusMap.find(indx)->second;
            (*it)->setIndx(spIndx);
        }
        else{
            indx = (*it)->getLindx();
            spIndx = spToLocusMap.find(indx)->second;
            (*it)->setIndx(spIndx);
        }
    }
    this->setTreeTipNames();
}

std::string GeneTree::printNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string geneTreeString = ss.str();
    return geneTreeString;
}

std::string GeneTree::printExtantNewickTree(){
    std::stringstream ss;
    recGetNewickTree(this->getRoot(), ss);
    ss << ";";
    std::string geneTreeString = ss.str();
    return geneTreeString;
}


void GeneTree::recGetExtNewickTree(Node *p, std::stringstream &ss, double brlen){ 
    if(p->getRdes() == NULL){
        if(p->getIsExtant())
            ss << p->getName();
    }
    else{
    // if(p != NULL){
        int flag = p->getFlag();
        if(flag == 2){
            ss << "(";
            recGetExtNewickTree(p->getRdes(), ss, brlen);
            ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
            ss << ",";
            recGetExtNewickTree(p->getLdes(), ss, brlen);
            ss << "[&index=" << p->getLdes()->getIndex() << "]" << ":" << p->getLdes()->getBranchLength();
            ss << ")";
        }
        else{
            // if(p->getRdes() == NULL){
            //     ss << p->getName();
            // }
            // else{
            if(p->getLdes()->getIsExtinct()){
                recGetExtNewickTree(p->getRdes(), ss, brlen);
            }
            if(p->getRdes()->getIsExtinct()){
                recGetExtNewickTree(p->getLdes(), ss, brlen);
            }
            //}
        }
   }
}

void GeneTree::recGetNewickTree(Node *p, std::stringstream &ss){
    if(p != NULL){
        if( p->getRdes() == NULL)
            ss << p->getName();
        else{
            ss << "(";
            recGetNewickTree(p->getRdes(), ss);
            ss << "[&index=" << p->getRdes()->getIndex() << "]" << ":" << p->getRdes()->getBranchLength();
            ss << ",";
            recGetNewickTree(p->getLdes(), ss);
            ss << "[&index=" << p->getLdes()->getIndex() << "]" << ":" << p->getLdes()->getBranchLength();
            ss << ")";
        }
    }
}

void GeneTree::setTreeTipNames(){
    int indNumber = 0;
    std::stringstream tn;
    std::string name;
    int indx;
    for(std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
        if((*it)->getIsTip()){
            indx  = (*it)->getIndex();
            tn << indx;
            name = tn.str();
            tn.clear();
            tn.str(std::string());
            indNumber++;
            tn << indNumber;
            name += "_" + tn.str();
            tn.clear();
            tn.str(std::string());
            if((*it) == this->getOutgroup())
                (*it)->setName("OUT");
            else
                (*it)->setName(name);
            if(indNumber == individualsPerPop)
                indNumber = 0;
        }

    }

}
