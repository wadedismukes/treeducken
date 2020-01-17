//
//  LocusTree.h
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/13/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#ifndef LocusTree_h
#define LocusTree_h

#include "SpeciesTree.h"
#include <algorithm>
#include <set>

class LocusTree : public Tree 
{
    private:
        double geneBirthRate, geneDeathRate, transferRate;
        double currentTime;
        double stopTime;
        unsigned numTaxa;
        unsigned numTransfers;

    public:
        LocusTree(unsigned nt, double stop, double gbr, double gdr, double lgtr);
        virtual         ~LocusTree();
        virtual double  getTimeToNextEvent();
        virtual void    lineageBirthEvent(unsigned indx);
        virtual void    lineageDeathEvent(unsigned indx);
        virtual void    setNewLineageInfo(int indx, Node *r, Node *s);
        void    lineageTransferEvent(int indx, bool randTrans);
        void    ermEvent(double ct);

        int     speciationEvent(int indx, double time, std::pair<int,int> sibs);
        void    extinctionEvent(int indx, double time);
        void    setNewIndices(int indx, std::pair<int,int> sibs, int count);
        std::string   printNewickTree();
        void    setTreeTipNames();
        void    recTipNamer(Node *p, unsigned &copyNumber);
        void    recGetNewickTree(Node *r, std::stringstream &ss);
        void    setBranchLengths();
        void    setPresentTime(double currentT);
        void    setStopTime(double st) {stopTime = st; currentTime = 0;}
        double  getCurrentTime() { return currentTime; }
        void    setCurrentTime(double ct) {currentTime = ct; }
        int     getNumberTransfers();
        int     chooseRecipientSpeciesID(Node *d);
        std::map<int,double>     getBirthTimesFromNodes();
        std::set<int>            getExtLociIndx();
        std::set<int>            getCoalBounds();
        std::multimap<int,double>     getDeathTimesFromNodes();
        std::multimap<int,double>     getDeathTimesFromExtinctNodes();
        std::map<int,int>             getLocusToSpeciesMap();
        std::vector< std::vector<int> >     getExtantLoci(std::set<double, std::greater<double> > epochSet);
        std::vector< std::string >    printSubTrees();
        int     postOrderTraversalStep(int indx);
    
    
    

};
#endif /* LocusTree_h*/
