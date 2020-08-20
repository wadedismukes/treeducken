#ifndef SpeciesTree_h
#define SpeciesTree_h

#include "Tree.h"
#include <sstream>
#include <map>
#include <set>
#include <RcppArmadillo.h>

using namespace Rcpp;

class SpeciesTree : public Tree
{
    private:

        double        speciationRate, extinctionRate;
        unsigned      extantStop;

    public:
                      SpeciesTree(unsigned numTaxa, double curTime, double specRate, double extRate);
                      SpeciesTree(unsigned numTaxa);
                      SpeciesTree(SEXP rtree);
                      SpeciesTree(const SpeciesTree& speciestree, unsigned numTaxa);
        virtual       ~SpeciesTree();

        SpeciesTree*  clone() const { return new SpeciesTree(*this); }
        void          setSpeciationRate(double sr) {speciationRate = sr; }
        void          setExtinctionRate(double er) {extinctionRate = er; }
        void          setCurrentTime(double et) { currentTime = et; }
        // tree-building functions
        virtual double        getTimeToNextEvent();
        double                getTimeToNextEventMoran();
        virtual void          lineageBirthEvent(unsigned indx);
        virtual void          lineageDeathEvent(unsigned indx);
        void          ermEvent(double curTime);
        void          setNewLineageInfo(unsigned indx, Node *r, Node *l);

        // set node parameters across tree
        void          setBranchLengths();
        void          setPresentTime(double currentT);
        void          setTreeTipNames();
        void          recTipNamer(Node *p, unsigned &extinctCount, unsigned &tipCount);

        // simulation functions
        void          setGSATipTreeFlags();
        void          reconstructTreeFromGSASim(Node *oRoot);
        void          setTreeInfo();
        void          popNodes();
        void          recPopNodes(Node *p);
        void          reconstructLineageFromGSASim(Node *currN, Node *prevN, unsigned &tipCounter, unsigned &intNodeCounter);
      //  void          setSampleFromFlags();
        std::map<int,int>           makeIndxMap();
        std::map<int, std::string>  makeTipMap();
        std::map<int,double>        getBirthTimesFromNodes();
        std::map<int,double>        getDeathTimesFromNodes();
        double                      getCurrentTimeFromExtant() {return extantNodes[0]->getDeathTime();}
        double                      getCurrentTime();
        bool                        getIsExtantFromIndx(int indx) { return nodes[indx]->getIsExtant(); }
        bool                        macroEvent(int indx);

        std::pair<int, int>         preorderTraversalStep(int index);
        int                         postOrderTraversalStep(int index);
        int           findLastToGoExtinct(double eventTime);
        int           getNodesIndxFromExtantIndx(int extanIndx) {return extantNodes[extanIndx]->getIndex(); }


};


#endif /* SpeciesTree_h */
