//
// Created by Dismukes, Wade T [EEOBS] on 10/31/19.
//

#ifndef SRC_SYMBIONTTREE_H
#define SRC_SYMBIONTTREE_H

#include "Tree.h"
#include <sstream>
#include <map>
#include <set>


class SymbiontTree : Tree {

private:
    double symbSpecRate, symbExtRate, hostExpanRate;
    double currentTime;
    double stopTime;
    unsigned numTaxa;
    unsigned numExpansions;
  //  SpeciesTree* hostTree; needed to add the header for Species Tree 

public:
    SymbiontTree(unsigned nt, double stop, double ssr, double ser, double her);
    virtual         ~SymbiontTree();
    virtual double  getTimeToNextEvent();
    virtual void    lineageBirthEvent(unsigned indx);
    virtual void    lineageDeathEvent(unsigned indx);
    virtual void    setNewLineageInfo(int indx, Node *r, Node *s);
    void            hostExpansionEvent(int indx);
    void            ermEvent(double ct);

    int             hostSpeciationEvent(int indx, double time, std::pair<int,int> sibs);
    void            hostExtinctionEvent(int indx, double time);
    void            setNewIndices(int indx, std::pair<int,int> sibs, int count);
    std::string     printNewickTree();
    void            setTreeTipNames();
    void            recTipNamer(Node *p, unsigned &copyNumber);
    void            recGetNewickTree(Node *r, std::stringstream &ss);
    void            setBranchLengths();
    void            setPresentTime(double currentT);
    void            setStopTime(double st) {stopTime = st; currentTime = 0;}
    double          getCurrentTime() { return currentTime; }
    void            setCurrentTime(double ct) {currentTime = ct; }
    int             getNumberExpansoins();
    std::map<int,double>     getBirthTimesFromNodes();
    std::set<int>            getExtSymbIndx();
    std::set<int>            getCoalBounds();
    std::multimap<int,double>     getDeathTimesFromNodes();
    std::multimap<int,double>     getDeathTimesFromExtinctNodes();
    std::map<int,int>             getSymbToHostsMap();
    std::vector< std::vector<int> >     getExtantLoci(std::set<double, std::greater<double> > epochSet);
    std::vector< std::string >    printSubTrees();
    int     postOrderTraversalStep(int indx);




};


#endif //SRC_SYMBIONTTREE_H
