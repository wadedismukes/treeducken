//
//  Simulator.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 11/9/17.
//  Copyright © 2017 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Simulator.h"
#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace arma;

Simulator::Simulator(unsigned nt, double lambda, double mu, double rho)
{
    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    simType = 1;
    currentSimTime = 0.0;
    numTaxaToSim = nt;
    gsaStop = 100*nt;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;

    treeScale = -1;
    numLoci = 0;
    numGenes = 0;
    geneBirthRate = 0.0;
    geneDeathRate = 0.0;
    transferRate = 0.0;
    propTransfer = 0.0;
    indPerPop = 0;
    popSize = 0;

}


Simulator::Simulator(unsigned ntax,
                     double lambda,
                     double mu,
                     double rho,
                     unsigned numLociToSim,
                     double gbr,
                     double gdr,
                     double lgtr)
{
    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    simType = 2;
    currentSimTime = 0.0;
    numTaxaToSim = ntax;
    gsaStop = 100*ntax;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;

    numLoci = numLociToSim;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    propTransfer = 0.0;
    indPerPop = 0;
    popSize = 0;


}

Simulator::Simulator(unsigned ntax,
                     double lambda,
                     double mu,
                     double rho,
                     unsigned numLociToSim,
                     double gbr,
                     double gdr,
                     double lgtr,
                     unsigned ipp,
                     unsigned Ne,
                     double genTime,
                     int ng,
                     double og,
                     double ts,
                     bool sout)
{
    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    simType = 3;
    currentSimTime = 0.0;
    numTaxaToSim = ntax;
    gsaStop = 100*ntax;
    speciationRate = lambda;
    extinctionRate = mu;
    samplingRate = rho;
    numLoci = numLociToSim;
    numGenes = ng;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    propTransfer = 0.0;
    indPerPop = ipp;
    popSize = Ne;
    printSOUT = sout;
    generationTime = genTime;
    outgroupFrac = og;
    geneTrees.resize(numLoci);
    treeScale = ts;
}

Simulator::Simulator(double stopTime,
          double hostSpeciationRate,
          double hostExtinctionRate,
          double symbSpeciationRate,
          double symbExtinctionRate,
          double switchingRate,
          double csr,
          double rho,
          int hl){

    speciationRate = hostSpeciationRate;
    extinctionRate = hostExtinctionRate;
    samplingRate = rho;
    cospeciationRate = csr;
    geneBirthRate = symbSpeciationRate;
    geneDeathRate = symbExtinctionRate;
    transferRate = switchingRate;
    timeToSim = stopTime;

    hostLimit = hl;

    spTree = nullptr;
    geneTree = nullptr;
    lociTree = nullptr;
    symbiontTree = nullptr;

}



Simulator::~Simulator(){
    for(std::vector<SpeciesTree*>::iterator p=gsaTrees.begin(); p != gsaTrees.end(); ++p){
        delete (*p);
    }
    gsaTrees.clear();
    int i = 0;
    for(std::vector<LocusTree*>::iterator p=locusTrees.begin(); p != locusTrees.end(); ++p){
        delete (*p);
        for(std::vector<GeneTree*>::iterator q=geneTrees[i].begin(); q != geneTrees[i].end(); ++q){
            delete (*q);
        }
        geneTrees[i].clear();
        ++i;
    }
    locusTrees.clear();


}

void Simulator::initializeSim(){
    spTree = new SpeciesTree(numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
}


/*
 Below is the machinery to use GSA sampling (Hartmann 2010) to simulate a species tree.
 Much of this code is modified from FossilGen (written by Tracy Heath)
 */
bool Simulator::gsaBDSim(){
    double timeIntv, sampTime;
    bool treeComplete = false;
    SpeciesTree st =  SpeciesTree(numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
    spTree = &st;
    double eventTime;

    while(gsaCheckStop()){
        eventTime = spTree->getTimeToNextEvent();
        currentSimTime += eventTime;
        spTree->ermEvent(currentSimTime);
        if(spTree->getNumExtant() < 1){
            treeComplete = false;
            return treeComplete;
        }
        else if(spTree->getNumExtant() == numTaxaToSim){
            timeIntv = spTree->getTimeToNextEvent();
            sampTime = (unif_rand() * timeIntv) + currentSimTime;
            spTree->setPresentTime(sampTime);
            processGSASim();
        }

    }
    unsigned gsaRandomTreeID = unif_rand() * (gsaTrees.size() - 1);
    // delete spTree;
    spTree = gsaTrees[gsaRandomTreeID];
    processSpTreeSim();
    spTree->setBranchLengths();
    spTree->setTreeTipNames();
    currentSimTime = spTree->getCurrentTimeFromExtant();
    if(treeScale > 0.0){
      spTree->scaleTree(treeScale, currentSimTime);
      currentSimTime = treeScale;
    }

    treeComplete = true;

    return treeComplete;
}


bool Simulator::gsaCheckStop(){

  bool keepSimulating = true;

  if(spTree->getNumExtant() >= gsaStop){
      keepSimulating = false;
  }

  return keepSimulating;
}

void Simulator::processGSASim(){
    SpeciesTree *tt = new SpeciesTree(numTaxaToSim + spTree->getNumExtinct());
    this->prepGSATreeForReconstruction();
    Node *simRoot = spTree->getRoot();
    tt->setRoot(simRoot);
    tt->reconstructTreeFromGSASim(simRoot);
    gsaTrees.push_back(std::move(tt));
}

void Simulator::processSpTreeSim(){
    spTree->setSpeciationRate(speciationRate);
    spTree->setExtinctionRate(extinctionRate);
    spTree->popNodes();
    spTree->setNumExtant();
    spTree->setNumExtinct();

}

void Simulator::prepGSATreeForReconstruction(){
    spTree->setGSATipTreeFlags();
}


bool Simulator::simSpeciesTree(){
    bool good = false;
    while(!good){
        good = gsaBDSim();
    }
    if(outgroupFrac > 0.0)
        this->graftOutgroup(spTree, spTree->getTreeDepth());
    return good;
}

bool Simulator::simHostSymbSpeciesTreePair(){
  bool good = false;
  while(!good){
    good = pairedBDPSim();
  }
  return good;
}


bool Simulator::pairedBDPSim(){
  bool treePairGood = false;
  currentSimTime = 0.0;
  double stopTime = this->getTimeToSim();

  spTree = new SpeciesTree(1, currentSimTime, speciationRate, extinctionRate);
  symbiontTree = new SymbiontTree(1,
                                  currentSimTime,
                                  geneBirthRate,
                                  geneDeathRate,
                                  transferRate,
                                  hostLimit);

  double eventTime;
  // TODO: the tree is not initialzied with links between hosts and parasites
  //symbiontTree->setSymbTreeInfo(spTree);
  assocMat = ones<umat>(1,1);
  while(currentSimTime < stopTime){

    eventTime = symbiontTree->getTimeToNextJointEvent(speciationRate,
                                                 extinctionRate,
                                                 cospeciationRate,
                                                 assocMat);
    currentSimTime += eventTime;
    assocMat = this->cophyloEvent(currentSimTime, assocMat);

    if(spTree->getNumExtant() < 1){
      treePairGood = false;
      return treePairGood;
    }
    else if(symbiontTree->getNumExtant() < 1){
      treePairGood = false;
      return treePairGood;
    }
  }
  treePairGood = true;
 // symbiontTree->reindexForR();
 //spTree->reindexForR();
  symbiontTree->setPresentTime(currentSimTime);
  spTree->setPresentTime(currentSimTime);

  return treePairGood;
}


arma::umat Simulator::cophyloEvent(double eventTime, arma::umat assocMat){
  double hostEvent = speciationRate + extinctionRate;
  double symbEvent = geneBirthRate + geneDeathRate + transferRate;
  double cospecEvent = cospeciationRate;
  // which tree then ermEvent whichever or cospeciation event
  double hostEventProb = hostEvent / (hostEvent + symbEvent + cospecEvent);
  double symbEventProb = symbEvent / (hostEvent + symbEvent + cospecEvent);
  symbEventProb += hostEventProb;
  double whichEvent = unif_rand();
  if(whichEvent < hostEventProb){
    assocMat = this->cophyloERMEvent(eventTime, assocMat);
  }
  else if(whichEvent < symbEventProb){
    assocMat = symbiontTree->ermJointEvent(eventTime, assocMat);
  }
  else{

    assocMat = this->cospeciationEvent(eventTime, assocMat);
  }
  return assocMat;
}


arma::umat Simulator::cophyloERMEvent(double eventTime, arma::umat assocMat){
  int numExtantHosts = spTree->getNumExtant();
  int nodeInd = unif_rand()*(numExtantHosts);
  double relBr = speciationRate / (speciationRate + extinctionRate);
  bool isBirth = (unif_rand() < relBr ? true : false);
  spTree->setCurrentTime(eventTime);
  symbiontTree->setCurrentTime(eventTime);
  int numExtantSymbs = symbiontTree->getNumExtant();

  arma::ucolvec cvec = assocMat.col(nodeInd);

  // List rowOfEvents;

  assocMat.shed_col(nodeInd);

  if(isBirth){
    spTree->lineageBirthEvent(nodeInd);
    numExtantHosts = spTree->getNumExtant();

    assocMat.resize(numExtantSymbs, numExtantHosts);

    assocMat(span::all, numExtantHosts-2) = cvec;
    assocMat(span::all, numExtantHosts-1) = cvec;

   // sort symbs on new hosts
   for(int i = 0; i < numExtantSymbs; i++){
      if(cvec[i] == 1 && i != nodeInd){
      // shuffle might be better here?
        arma::umat rr = zeros<umat>(1,2);
        rr(0,0) = 1;
        arma::umat shuffledRandomRow = shuffle(rr);
        assocMat(i, span(numExtantHosts-2,numExtantHosts-1)) = shuffledRandomRow;
      }
    }
  }
  else{

    spTree->lineageDeathEvent(nodeInd);

    numExtantHosts = spTree->getNumExtant();
    // check rows for 0's
    //
    arma::uvec hostless = zeros<uvec>(numExtantSymbs);
    for(int i = assocMat.n_rows - 1; i > -1; i--){

      if(sum(assocMat.row(i)) < 1){
        symbiontTree->lineageDeathEvent(i);
        hostless(i) = 1;
      }
    }

    uvec toBeDeleted = find(hostless);
    assocMat.shed_rows(toBeDeleted);
  }

  return assocMat;
}


arma::umat Simulator::cospeciationEvent(double eventTime, arma::umat assocMat){
  // draw index of host

  spTree->setCurrentTime(eventTime);
  symbiontTree->setCurrentTime(eventTime);
  // issue is that here you can choose hosts without symbionts!
  int numExtantHosts = spTree->getNumExtant();
  std::vector<int> hostIndices;
  hostIndices.resize(numExtantHosts);

  for(int i = 0; i < numExtantHosts; i++){
    if(sum(assocMat.col(i)) >= 0.5)
      hostIndices[i] = i;
  }

  int indxOfHost = unif_rand() * hostIndices.size(); //col of assocMat
 // Rcout << "hostIndx = " << indxOfHost << std::endl;

  arma::ucolvec cvec = assocMat.col(hostIndices[indxOfHost]);
  arma::uvec symbIndices = find(cvec);
  //assocMat.print();
  int indxOfSymb = unif_rand() * symbIndices.size();

  arma::urowvec rvec = assocMat.row(symbIndices[indxOfSymb]);

  spTree->lineageBirthEvent(hostIndices[indxOfHost]);
  symbiontTree->lineageBirthEvent(symbIndices[indxOfSymb]);

  numExtantHosts = spTree->getNumExtant();
  int numExtantSymbs = symbiontTree->getNumExtant();

  assocMat.shed_col(hostIndices[indxOfHost]);
  assocMat.shed_row(symbIndices[indxOfSymb]);

  assocMat.resize(numExtantSymbs, numExtantHosts);

  assocMat.submat(numExtantSymbs-2,
                  numExtantHosts-2,
                  numExtantSymbs-1,
                  numExtantHosts-1) = eye<umat>(2,2);

  //List rowOfEvents;
  // rowOfEvents["Host Index"] = spTree->getNodesIndxFromExtantIndx(numExtantHosts - 2);
  // rowOfEvents["Symb Index"] = symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 2);
  // rowOfEvents["Event Type"] = 1;
  // rowOfEvents["Time"] = eventTime;
  // events.push_back(rowOfEvents);
  //
  // rowOfEvents.erase(rowOfEvents.begin(), rowOfEvents.end());
  // rowOfEvents["Host Index"]  = spTree->getNodesIndxFromExtantIndx(numExtantHosts - 1);
  // rowOfEvents["Symb Index"] = symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1);
  // rowOfEvents["Event Type"]  = 1;
  // rowOfEvents["Time"] = eventTime;
  // events.push_back(rowOfEvents);
  // rowOfEvents.erase(rowOfEvents.begin(), rowOfEvents.end());

 //  assocMat.print();

  // if only cospeciation is on don't do this
  // loop through cvec and make little vecs
  if(speciationRate > 0.0 || transferRate > 0.0){
    for(int i = 0; i < cvec.size() - 1; i++){
      if(cvec[i] == 1){
        arma::umat rr = ones<umat>(1,2);
        int randOne = unif_rand() * 2;
        if(randOne == 0){
          rr(0, 0) = 0;
         // rowOfEvents["Host Index"] = spTree->getNodesIndxFromExtantIndx(numExtantHosts - 1);

        }
        else{
          rr(0, 1) = 0;
         // rowOfEvents["Host Index"] = spTree->getNodesIndxFromExtantIndx(numExtantHosts - 2);
        }

//
//         rowOfEvents["Symb Index"] = symbiontTree->getNodesIndxFromExtantIndx(i);
//         rowOfEvents["Event Type"] = 3;
//         rowOfEvents["Time"] = eventTime;

        assocMat(i, span(numExtantHosts-2,numExtantHosts-1)) = rr;

      //  events.push_back(rowOfEvents);
    //    rowOfEvents.erase(rowOfEvents.begin(), rowOfEvents.end());

      }
    }
  }

  if(geneBirthRate > 0.0 || transferRate > 0.0){
    // loop through rvec
    for(int i = 0; i < rvec.size() - 1; i++){
      if(rvec[i] == 1){
        arma::umat rr = ones<umat>(1,2);
    //    rowOfEvents["Host Index"] = spTree->getNodesIndxFromExtantIndx(i);

        int randOne = unif_rand() * 2;
        if(randOne == 0){
          rr(0,0) = 0;
    //      rowOfEvents["Symb Index"] = symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1);

        }
        else{
          rr(0, 1) = 0;
    //      rowOfEvents["Symb Index"] = symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 2);
        }
 //       rowOfEvents["Event Type"] = 3;
        // rowOfEvents["Time"] = eventTime;
        // events.push_back(rowOfEvents);
        // rowOfEvents.erase(rowOfEvents.begin(), rowOfEvents.end());

        assocMat(span(numExtantSymbs-2,numExtantSymbs-1),i) = rr.t();
      }
    }
  }
  return assocMat;
}

// List Simulator::updateIndices(List evnts){
//
//   // List eventRow;
//   // for(int i = 0; i < evnts.size(); i++){
//   //   eventRow = evnts[i];
//   //   int oldHostInd = eventRow["Host Index"];
//   //   int oldSymbInd = eventRow["Symb Index"];
//   //   int newHostInd = spTree->getIndexFromNodes(oldHostInd);
//   //   int newSymbInd = symbiontTree->getIndexFromNodes(oldSymbInd);
//   //   eventRow["Host Index"] = newHostInd;
//   //   eventRow["Symb Index"] = newSymbInd;
//   //   evnts[i] = eventRow;
//   // }
//
//   return evnts;
// }

std::string Simulator::printExtSpeciesTreeNewick(){
    SpeciesTree *tt = new SpeciesTree(numTaxaToSim);
    spTree->getRootFromFlags(false);
    if(outgroupFrac > 0.0){
        tt->setOutgroup(spTree->getOutgroup());
        tt->setRoot(spTree->getOutgroup()->getAnc());
    }
    else{
        tt->setRoot(spTree->getExtantRoot());
    }
    tt->setExtantRoot(tt->getRoot());
    tt->reconstructTreeFromSim(spTree->getRoot());
    std::string newickTree = tt->printExtNewickTree();
    delete tt;
    tt = nullptr;
    return newickTree;
}

std::string Simulator::printSpeciesTreeNewick(){
    return spTree->printNewickTree();
}

bool Simulator::bdsaBDSim(){

    bool treesComplete = false;
    double stopTime = spTree->getCurrentTimeFromExtant();

    double eventTime;
    bool isSpeciation;
    lociTree = new LocusTree(numTaxaToSim, currentSimTime, geneBirthRate, geneDeathRate, transferRate);
    Rcout << "end poooo" << std::endl;

    std::map<int,double> speciesBirthTimes = spTree->getBirthTimesFromNodes();
    std::map<int,double> speciesDeathTimes = spTree->getDeathTimesFromNodes();
    std::set<int> contempSpecies;
    std::pair<int, int> sibs;

    Node* spRoot = spTree->getRoot();
    lociTree->setStopTime(stopTime);
    currentSimTime = 0;
    if(!(contempSpecies.empty()))
        contempSpecies.clear();
    contempSpecies.insert(spRoot->getIndex());

    while(currentSimTime < stopTime){
        eventTime = lociTree->getTimeToNextEvent();
        currentSimTime += eventTime;
        for(std::set<int>::iterator it = contempSpecies.begin(); it != contempSpecies.end();){
            if(currentSimTime > speciesDeathTimes[(*it)]){
                isSpeciation = spTree->macroEvent((*it));
                if(isSpeciation){
                    sibs = spTree->preorderTraversalStep(*it);
                    lociTree->speciationEvent((*it), speciesDeathTimes[(*it)], sibs);
                    it = contempSpecies.erase(it);
                    it = contempSpecies.insert( it, sibs.second);
                    ++it;
                    it = contempSpecies.insert( it, sibs.first);
                }
                else{
                    if(!(spTree->getIsExtantFromIndx(*it))){
                        lociTree->extinctionEvent(*it, speciesDeathTimes[(*it)]);
                        it = contempSpecies.erase(it);
                    }
                    else{
                        ++it;
                    }
                }
            }
            else{
                ++it;
            }

            if(lociTree->getNumExtant() < 1){
                treesComplete = false;
                return treesComplete;
            }

        }

        if(currentSimTime >= stopTime){
            currentSimTime = stopTime;
            lociTree->setCurrentTime(stopTime);
        }
        else{
            lociTree->ermEvent(currentSimTime);
        }



    }
    lociTree->setPresentTime(currentSimTime);
    treesComplete = true;

    return treesComplete;
}

bool Simulator::simSpeciesLociTrees(){
    bool good = false;
    bool spGood = false;
    for(int i = 0; i < numLoci; i++){
        while(!good){
            while(!spGood){
                spGood = gsaBDSim();
            }
            if(outgroupFrac > 0.0)
                this->graftOutgroup(spTree, spTree->getTreeDepth());
            if(printSOUT)
                std::cout << "Simulating loci #" <<  i + 1 << std::endl;
            good = bdsaBDSim();
        }
        if(outgroupFrac > 0.0)
            this->graftOutgroup(lociTree, lociTree->getTreeDepth());

        locusTrees.push_back(lociTree);

        good = false;
    }
    return good;
}

bool Simulator::simLocusTree(){
  bool good = false;
  Rcout << "start poooo" << std::endl;

  while(!good){
    good = bdsaBDSim();

  }
  return good;
}

std::string Simulator::printLocusTreeNewick(int i){
    std::string newickTree;
    std::vector<LocusTree*>::iterator it = locusTrees.begin();
    std::advance(it, i);
    newickTree = (*it)->printNewickTree();
    return newickTree;
}


std::set<double, std::greater<double> > Simulator::getEpochs(){
    std::set<double, std::greater<double> > epochs;
    std::vector<Node*> lociTreeNodes = lociTree->getNodes();
    for(std::vector<Node*>::iterator it = lociTreeNodes.begin(); it != lociTreeNodes.end(); ++it){
        if(!((*it)->getIsExtinct())){
            if((*it)->getIsTip())
                epochs.insert((*it)->getDeathTime());
            epochs.insert((*it)->getBirthTime());
        }
        else
            epochs.insert((*it)->getDeathTime());
    }
    return epochs;
}


bool Simulator::coalescentSim(){
    bool treeGood = false;
    geneTree = new GeneTree(numTaxaToSim, indPerPop, popSize, generationTime);

    std::map<int,int> spToLo;

    int ancIndx;
    int epochCount = 0;

    double stopTime, stopTimeEpoch, stopTimeLoci;
    bool allCoalesced = false, deathCheck = false;
    bool is_ext;

    std::set<double, std::greater<double> > epochs = getEpochs();
    int numEpochs = (int) epochs.size();
    std::set<int> extinctFolks = lociTree->getExtLociIndx();
    std::set<int> coalescentBounds = lociTree->getCoalBounds();
    std::vector< std::vector<int> > contempLoci = lociTree->getExtantLoci(epochs);
    std::map<int, double> stopTimes = lociTree->getBirthTimesFromNodes();
    geneTree->initializeTree(contempLoci, *(epochs.begin()));
    if(outgroupFrac != 0.0)
        contempLoci[0].pop_back();
    std::set<int>::iterator extFolksIt;

    for(std::set<double, std::greater<double> >::iterator epIter = epochs.begin(); epIter != epochs.end(); ++epIter){
        currentSimTime = *epIter;
        if(epochCount != numEpochs - 1){
            epIter = std::next(epIter, 1);
            stopTimeEpoch = *epIter;
            for(int j = 0; j < contempLoci[epochCount].size(); ++j){
                extFolksIt = extinctFolks.find(contempLoci[epochCount][j]);
                is_ext = (extFolksIt != extinctFolks.end());
                if(is_ext){
                    geneTree->addExtinctSpecies(currentSimTime, contempLoci[epochCount][j]);
                    extinctFolks.erase(extFolksIt);
                }
                stopTimeLoci = stopTimes[contempLoci[epochCount][j]];

                if(stopTimeLoci > stopTimeEpoch){
                    stopTime = stopTimeLoci;
                    deathCheck = true;
                }
                else{
                    stopTime = stopTimeEpoch;
                    deathCheck = false;
                }

                ancIndx = lociTree->postOrderTraversalStep(contempLoci[epochCount][j]);
                allCoalesced = geneTree->censorCoalescentProcess(currentSimTime, stopTime, contempLoci[epochCount][j], ancIndx, deathCheck);


                // if all coalesced remove that loci from the matrix of loci
                if(allCoalesced){
                    int check = contempLoci[epochCount][j];
                    for(int k = epochCount + 1; k < numEpochs; k++){
                        for(int m = 0; m < contempLoci[k].size(); ++m){
                            if(contempLoci[k][m] == check){
                                contempLoci[k].erase(contempLoci[k].begin() + m);
                                break;
                            }
                        }
                    }
                }
                allCoalesced = false;
                is_ext = false;
            }
            epIter = std::prev(epIter, 1);
        }
        else{
            // finish coalescing
            geneTree->rootCoalescentProcess(currentSimTime, outgroupFrac);
            treeGood = true;
        }
        epochCount++;
    }

    spToLo = lociTree->getLocusToSpeciesMap();
    geneTree->setIndicesBySpecies(spToLo);
    return treeGood;
}

bool Simulator::simThreeTree(){
    bool gGood = false;
    bool spGood = false;
    bool loGood = false;
    while(!spGood){
        spGood = gsaBDSim();

    }
    for(int i = 0; i < numLoci; i++){
        while(!loGood){
            if(printSOUT)
                Rcout << "Simulating loci # " <<  i + 1 << std::endl;
            loGood = bdsaBDSim();
        }
        if(outgroupFrac > 0.0){
            this->graftOutgroup(lociTree, lociTree->getTreeDepth());
        }
        for(int j = 0; j < numGenes; j++){
            while(!gGood){
                if(printSOUT)
                    Rcout << "Simulating gene # " <<  j + 1 << " of loci # " << i + 1 << std::endl;
                gGood = coalescentSim();
            }
            geneTrees[i].push_back(geneTree);

            gGood = false;
        }
        locusTrees.push_back(lociTree);
        loGood = false;
    }
    if(outgroupFrac > 0.0)
        this->graftOutgroup(spTree, spTree->getTreeDepth());

    return gGood;
}


std::string Simulator::printGeneTreeNewick(int i, int j){
    std::string newickTree;
    newickTree = geneTrees[i][j]->printNewickTree();
    return newickTree;
}

std::string Simulator::printExtantGeneTreeNewick(int i, int j){
    std::string newickTree;
    GeneTree *tt = new GeneTree(numTaxaToSim, indPerPop, popSize, generationTime);

    geneTrees[i][j]->getRootFromFlags(true);


    if(outgroupFrac > 0.0){
        tt->setOutgroup(geneTrees[i][j]->getOutgroup());
        tt->setRoot(geneTrees[i][j]->getOutgroup()->getAnc());
    }
    else
        tt->setRoot(geneTrees[i][j]->getExtantRoot());
    tt->setExtantRoot(geneTrees[i][j]->getExtantRoot());
    tt->reconstructTreeFromSim(geneTrees[i][j]->getRoot());
    newickTree = tt->printNewickTree();
    delete tt;
    tt = nullptr;
    return newickTree;
}


void Simulator::graftOutgroup(Tree *tr, double trDepth){
    Node *rootN = new Node();
    Node *currRoot = tr->getRoot();
    rootN->setBirthTime(currRoot->getBirthTime());
    Node *outgroupN = new Node();
    tr->rescaleTreeByOutgroupFrac(outgroupFrac, trDepth);
    double tipTime = tr->getEndTime();
    tr->setNewRootInfo(rootN, outgroupN, currRoot, tipTime);
}

bool Simulator::simLocusGeneTrees(){
    bool loGood = false;
    bool gGood = false;
    for(int i = 0; i < numLoci; i++){
        while(!loGood){
            if(printSOUT)
                std::cout << "Simulating loci # " <<  i + 1 << std::endl;
            loGood = bdsaBDSim();
        }
        for(int j = 0; j < numGenes; j++){
            while(!gGood){
                if(printSOUT)
                    std::cout << "Simulating gene # " <<  j + 1 << " of loci # " << i + 1 << std::endl;
                gGood = coalescentSim();
            }
            geneTrees[i].push_back(geneTree);

            gGood = false;
        }
        locusTrees.push_back(lociTree);
        loGood = false;
    }
    return gGood;
}

double Simulator::calcSpeciesTreeDepth(){
    return spTree->getTreeDepth();
}

double Simulator::calcExtantSpeciesTreeDepth(){
    SpeciesTree *tt = new SpeciesTree(numTaxaToSim);
    spTree->getRootFromFlags(false);
    tt->setRoot(spTree->getExtantRoot());
    tt->setExtantRoot(tt->getRoot());
    tt->reconstructTreeFromSim(tt->getExtantRoot());
    double extTreeDepth = tt->getTreeDepth();
    delete tt;
    tt = nullptr;
    return extTreeDepth;
}

double Simulator::calcLocusTreeDepth(int i){
    return locusTrees[i]->getTreeDepth();
}

int Simulator::findNumberTransfers(){
    int numTrans = 0;
    for(int i = 0; i < locusTrees.size(); i++){
        numTrans += locusTrees[i]->getNumberTransfers();
    }
    return numTrans;
}

double Simulator::findTMRCAGeneTree(int i, int j){
    return geneTrees[i][j]->getTreeDepth();
}

double Simulator::getSpeciesTreeRootEdge(){
  return spTree->getRoot()->getDeathTime() - spTree->getRoot()->getBirthTime();
}

double Simulator::getLocusTreeRootEdge(){
  return lociTree->getRoot()->getDeathTime() - lociTree->getRoot()->getBirthTime();
}

double Simulator::getSymbiontTreeRootEdge(){
  return symbiontTree->getRoot()->getDeathTime() - symbiontTree->getRoot()->getBirthTime();
}

