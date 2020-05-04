//
//  GeneTree.hpp
//  treeducken
//
//  Created by Dismukes, Wade T [EEOBS] on 12/20/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef GeneTree_h
#define GeneTree_h

#include "LocusTree.h"
#include <algorithm>

class GeneTree : public Tree {
    private:
        unsigned individualsPerPop;
        unsigned popSize;
        double   generationTime; // specified in generations per year

    public:
                    GeneTree(unsigned nt, unsigned ipp, unsigned ne, double genTime);
        virtual     ~GeneTree();
        double      getCoalTime(int n); // what do you need to determine this?
        Node*       coalescentEvent(double t, Node *p, Node *q);
        bool        censorCoalescentProcess(double startTime, double stopTime, int contempSpIndx, int newSpIndx, bool chck);
        void        initializeTree(std::vector< std::vector<int> > extantLociIndx, double presentTime);
        std::multimap<int,double> rescaleTimes(std::multimap<int, double> timeMap);
        void        rootCoalescentProcess(double startTime, double ogf);
        void        recursiveRescaleTimes(Node *r, double add);
        void        setBranchLengths();
        void        setIndicesBySpecies(std::map<int,int> spToLocusMap);
        std::string printNewickTree();
        std::string printExtantNewickTree();
        void        recGetNewickTree(Node *r, std::stringstream &ss);
        void        recGetExtNewickTree(Node *r, std::stringstream &ss, double brlen);
        void        setTreeTipNames();
        void        addExtinctSpecies(double bt, int indx);
        NumericMatrix        getGeneEdges();
        void        reindexForR();
};

#endif /* GeneTree_hpp */
