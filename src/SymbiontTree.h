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
      unsigned hostLimit;
      std::map<int,std::vector<int>> symbHostMap; // keys are symb indices

    public:
      SymbiontTree(int nt,
                   double currSimTime,
                   double symbsr,
                   double symber,
                   double hostExpanRate,
                   int hostLimit);
      SymbiontTree(const SymbiontTree& symbionttree,
                   unsigned numTaxa);

      virtual         ~SymbiontTree();
      double  getTimeToNextJointEvent(double hostSpecRate,
                                         double hostExtRate,
                                         double cospeciaRate,
                                         arma::umat assocMat);
      virtual void    lineageBirthEvent(unsigned indx);
      virtual void    lineageDeathEvent(unsigned indx);
      virtual void    setNewLineageInfo(int indx, Node *r, Node *s);
      void            setNewLineageInfoExpan(int indx,
                                             Node* r,
                                             Node* s,
                                             int hostIndx);
      void            hostExpansionEvent(int indx, int hostIndx);
      arma::umat       ermJointEvent(double ct, arma::umat assocMat);

      void            setSymbTreeInfoSpeciation(int ancIndx, int desIndx);
      void            setSymbTreeInfoExtinction(int deadIndx);

      //std::string     printNewickTree();
      void            setTreeTipNames();
      void            recTipNamer(Node *p, unsigned &extinctCount, unsigned &tipCount);

//      void            recGetNewickTree(Node *r, std::stringstream &ss);
      void            setBranchLengths();
      void            setPresentTime(double currentT);
      void            setStopTime(double st) { stopTime = st; currentTime = 0;}
      double          getCurrentTime() { return currentTime; }
      void            setCurrentTime(double ct) { currentTime = ct; }
      int             getNumberExpansions() {return numExpansions; }
      int             getNumHostSymbPairs() { return symbHostMap.size(); }

      std::vector<int>  getSymbsOnHost(int hostIndx);
      void            updateCurrentMap(int oldHostIndx, int newHostIndx);
      int             getExtantIndxFromNodes(int extantNodesIndx);
      void            cospeciationMapUpdate(int oldHostIndx,
                                            int numHosts,
                                            int symbIndx);
      void            updateHostsInNodes();
      int             getNodesIndxFromExtantIndx(int i) { return extantNodes[i]->getIndex(); }

};


#endif //SRC_SYMBIONTTREE_H
