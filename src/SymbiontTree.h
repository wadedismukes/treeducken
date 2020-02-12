//
// Created by Dismukes, Wade T [EEOBS] on 10/31/19.
//

#ifndef SRC_SYMBIONTTREE_H
#define SRC_SYMBIONTTREE_H

#include "Tree.h"
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include "SpeciesTree.h"

class SymbiontTree : public Tree {

    private:
      double symbSpecRate, symbExtRate, hostExpanRate;
      double currentTime;
      double stopTime;
      unsigned numTaxa;
      unsigned numExpansions;
      unsigned numPairs;
      unsigned hostLimit;
      std::multimap<int,int> symbHostMap; // keys are symb indices

    public:
      SymbiontTree(int nt,
                   double currSimTime,
                   double symbsr,
                   double symber,
                   double hostExpanRate,
                   int hostLimit);
      virtual         ~SymbiontTree();
      virtual double  getTimeToNextEvent(double hostSpecRate,
                                         double hostExtRate,
                                         double cospeciaRate,
                                         int numHosts);
      virtual void    lineageBirthEvent(unsigned indx);
      virtual void    lineageDeathEvent(unsigned indx);
      virtual void    setNewLineageInfo(int indx, Node *r, Node *s);
      void            hostExpansionEvent(int indx);
      void            ermEvent(double ct);

      void            setSymbTreeInfoSpeciation(int ancIndx, int desIndx);
      void            setSymbTreeInfoExtinction(int deadIndx);

      std::string     printNewickTree();
      void            setTreeTipNames();
      void            recTipNamer(Node *p, unsigned &copyNumber);
      void            recGetNewickTree(Node *r, std::stringstream &ss);
      void            setBranchLengths();
      void            setPresentTime(double currentT);
      void            setStopTime(double st) {stopTime = st; currentTime = 0;}
      double          getCurrentTime() { return currentTime; }
      void            setCurrentTime(double ct) {currentTime = ct; }
      int             getNumberExpansions();
      int             getNumHostSymbPairs() { return symbHostMap.size(); }
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
