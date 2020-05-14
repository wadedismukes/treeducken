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
    popSize = 0.0;

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
                     double Ne,
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
    geneTrees.resize(ng);
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
        for(std::vector<GeneTree*>::iterator q=geneTrees.begin(); q != geneTrees.end(); ++q){
            delete (*q);
        }
        geneTrees.clear();
        ++i;
    }
    locusTrees.clear();


}

void Simulator::initializeSim(){
    spTree = new SpeciesTree(numTaxaToSim, currentSimTime, speciationRate, extinctionRate);
}


void Simulator::initializeEventVector(){
  inOrderVecOfHostIndx.push_back(0);
  inOrderVecOfSymbIndx.push_back(0);
  inOrderVecOfEvent.push_back("I");
  inOrderVecOfEventTimes.push_back(0.0);
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
    return good;
}


bool Simulator::simSpeciesTreeTime(){
  bool good = false;
  while(!good){
    good = bdSimpleSim();
  }
  return good;
}


bool Simulator::bdSimpleSim(){
  bool treeComplete = false;
  currentSimTime = 0.0;
  double stopTime = this->getTimeToSim();
  double eventTime;

  spTree = new SpeciesTree(1, currentSimTime, speciationRate, extinctionRate);
  while(currentSimTime < stopTime){

    eventTime = spTree->getTimeToNextEvent();
    currentSimTime += eventTime;

    if(currentSimTime >= stopTime){
      currentSimTime = stopTime;
      spTree->setPresentTime(currentSimTime);
    }
    else{
      spTree->ermEvent(currentSimTime);
    }

    if(spTree->getNumExtant() < 1){
      delete spTree;
      treeComplete = false;
      return treeComplete;
    }
  }

  if(spTree->getNumExtant() <= 1){
    delete spTree;
    treeComplete = false;
    return treeComplete;
  }

  treeComplete = true;
  currentSimTime = stopTime;


  spTree->setPresentTime(currentSimTime);
  return treeComplete;
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
  this->initializeEventVector();
  assocMat = ones<umat>(1,1);
  while(currentSimTime < stopTime){
    eventTime = symbiontTree->getTimeToNextJointEvent(speciationRate,
                                                 extinctionRate,
                                                 cospeciationRate,
                                                 assocMat);
    currentSimTime += eventTime;
    if(currentSimTime >= stopTime){
      currentSimTime = stopTime;
      spTree->setPresentTime(stopTime);
      symbiontTree->setPresentTime(stopTime);
    }
    else{
      assocMat = this->cophyloEvent(currentSimTime, assocMat);
    }

    if(spTree->getNumExtant() < 1 || symbiontTree->getNumExtant() < 1 ||
       assocMat.n_rows < 1 || assocMat.n_cols < 1){
      treePairGood = false;
      this->clearEventDFVecs();
      delete spTree;
      delete symbiontTree;
      return treePairGood;
    }
  }
  if(spTree->getNumExtant() <= 1 || symbiontTree->getNumExtant() <= 1){
    treePairGood = false;
    this->clearEventDFVecs();
    delete spTree;
    delete symbiontTree;
    return treePairGood;
  }
  treePairGood = true;
  currentSimTime = stopTime;

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
    // assocMat = symbiontTree->ermJointEvent(eventTime, assocMat);
    assocMat = this->symbiontTreeEvent(eventTime, assocMat);
  }
  else{

    assocMat = this->cospeciationEvent(eventTime, assocMat);
  }
  return assocMat;
}


Rcpp::DataFrame Simulator::createEventDF(){
  this->updateEventIndices();
  DataFrame df = DataFrame::create(Named("Host Index") = inOrderVecOfHostIndx,
                                   Named("Symbiont Index") = inOrderVecOfSymbIndx,
                                   Named("Event Type") = inOrderVecOfEvent,
                                   Named("Event Time") = inOrderVecOfEventTimes);
  return df;
}


void Simulator::clearEventDFVecs(){
  inOrderVecOfHostIndx.erase(inOrderVecOfHostIndx.begin(),
                             inOrderVecOfHostIndx.end());
  inOrderVecOfSymbIndx.erase(inOrderVecOfSymbIndx.begin(),
                             inOrderVecOfSymbIndx.end());
  inOrderVecOfEvent.erase(inOrderVecOfEvent.begin(),
                          inOrderVecOfEvent.end());
  inOrderVecOfEventTimes.erase(inOrderVecOfEventTimes.begin(),
                               inOrderVecOfEventTimes.end());
}
void Simulator::updateEventIndices(){

  for(int i = 0; i < inOrderVecOfHostIndx.size(); i++){
    int oldHostIndx = inOrderVecOfHostIndx(i);
    int oldSymbIndx = inOrderVecOfSymbIndx(i);
    int newHostIndx = spTree->getIndexFromNodes(oldHostIndx);

    int newSymbIndx = symbiontTree->getIndexFromNodes(oldSymbIndx);
    // inOrderVecOfHostIndx.erase(i);
    // inOrderVecOfHostIndx.insert(i,newHostIndx);
    inOrderVecOfHostIndx(i) = newHostIndx;
    inOrderVecOfSymbIndx(i) = newSymbIndx;
    // inOrderVecOfSymbIndx.erase(i);
    // inOrderVecOfSymbIndx.insert(i,newSymbIndx);
  }
}


void Simulator::updateEventVector(int h, int s, int e, double time){
  inOrderVecOfHostIndx.push_back(std::move(h));
  inOrderVecOfSymbIndx.push_back(std::move(s));
  switch(e) {
    case 0:
      inOrderVecOfEvent.push_back(std::move("SL"));
      break;
    case 1:
      inOrderVecOfEvent.push_back(std::move("HL"));
      break;
    case 2:
      inOrderVecOfEvent.push_back(std::move("SG"));
      break;
    case 3:
      inOrderVecOfEvent.push_back("HG");
      break;
    case 4:
      inOrderVecOfEvent.push_back("AG");
      break;
    case 5:
      inOrderVecOfEvent.push_back("AL");
      break;
    case 6:
      inOrderVecOfEvent.push_back("C");
      break;
    default:
      Rcout << "not sure what happened there folks." << std::endl;
  }
  inOrderVecOfEventTimes.push_back(std::move(time));
}




arma::umat Simulator::symbiontTreeEvent(double eventTime, arma::umat assocMat){
  int numExtantSymbs = symbiontTree->getNumTips();
  int nodeInd = unif_rand()*(numExtantSymbs - 1);
  double relBr = geneBirthRate / (geneBirthRate
                                    + geneDeathRate
                                    + transferRate);
  double relDr = relBr + (geneDeathRate / (geneBirthRate
                                            + geneDeathRate
                                            + transferRate));
  double decid = unif_rand();
  spTree->setCurrentTime(eventTime);
  symbiontTree->setCurrentTime(eventTime);
  int numExtantHosts = spTree->getNumExtant();

  arma::urowvec rvec = assocMat.row(nodeInd);

  // List rowOfEvents;
  assocMat.shed_row(nodeInd);

  if(decid < relBr){
    updateEventVector(spTree->getNodesIndxFromExtantIndx(assocMat.n_cols - 1),
                      symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                      2,
                      eventTime);
    symbiontTree->lineageBirthEvent(nodeInd);
    numExtantSymbs = symbiontTree->getNumExtant();


    assocMat.resize(numExtantSymbs, numExtantHosts);

    assocMat(numExtantSymbs-2, span::all) = rvec;
    assocMat(numExtantSymbs-1, span::all) = rvec;

    // sort symbs on new hosts
    for(int i = 0; i < rvec.n_cols; i++){
      if(rvec(i) == 1){
        // shuffle might be better here?
        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                          4,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                          4,
                          eventTime);
      }
    }
  }
  else if(decid < relDr){
    //assocMat.print();
    updateEventVector(spTree->getNodesIndxFromExtantIndx(assocMat.n_cols - 1),
                      symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                      0,
                      eventTime);
    numExtantSymbs = symbiontTree->getNumExtant();
    // check rows for 0's
    // arma::uvec symbless = zeros<uvec>(numExtantHosts);
    // for(int i = assocMat.n_cols - 1; i != -1; i--){
    //   if(!(any(assocMat.col(i)))){
    //     updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
    //                       symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
    //                       0,
    //                       eventTime);
    //     spTree->lineageDeathEvent(i);
    //
    //     symbless(i) = 1;
    //   }
    // }
    for(int i = 0; i< rvec.n_cols; i++){
      if(rvec[i] == 1){
        // this may be a waste of code
            updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                              symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                              4,
                              eventTime);
      }
    }
    symbiontTree->lineageDeathEvent(nodeInd);
//
//     uvec toBeDeleted = find(symbless);
//     assocMat.shed_rows(toBeDeleted);
  }
  else{
    updateEventVector(spTree->getNodesIndxFromExtantIndx(assocMat.n_cols - 1),
                      symbiontTree->getNodesIndxFromExtantIndx(nodeInd),
                      2,
                      eventTime);
    if(sum(rvec) != assocMat.n_cols){
      std::vector<int> hostIndices;
      for(int i = 0; i < rvec.n_cols; i++){
        if(rvec(i) < 1)
          hostIndices.push_back(std::move(i));
      }
      //uvec hostIndices = find(rvec < 1);
      int hostInd = unif_rand() * (hostIndices.size());
      symbiontTree->lineageBirthEvent(nodeInd);
      numExtantSymbs = symbiontTree->getNumTips();

      assocMat.resize(numExtantSymbs, numExtantHosts);

      assocMat(numExtantSymbs-2, span::all) = rvec;
      rvec(hostIndices[hostInd]) = 1;
      assocMat(numExtantSymbs-1, span::all) = rvec;

      // sort symbs on new hosts
      for(int i = 0; i < rvec.n_cols; i++){
        if(rvec(i) == 1){
          updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                            symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                            4,
                            eventTime);
            updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                              symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                              4,
                              eventTime);
        }
      }
    }
    else{
      symbiontTree->lineageBirthEvent(nodeInd);
      numExtantSymbs = symbiontTree->getNumExtant();


      assocMat.resize(numExtantSymbs, numExtantHosts);

      assocMat(numExtantSymbs-2, span::all) = rvec;
      assocMat(numExtantSymbs-1, span::all) = rvec;

      // sort symbs on new hosts
      for(int i = 0; i < rvec.n_cols; i++){
        if(rvec(i) == 1){
          // shuffle might be better here?
          updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                            symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                            4,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                            symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                            4,
                            eventTime);
        }
      }

    }
  }
  return assocMat;
}




arma::umat Simulator::cophyloERMEvent(double eventTime, arma::umat assocMat){
  int numExtantHosts = spTree->getNumExtant();
  int nodeInd = unif_rand()*(numExtantHosts - 1);
  double relBr = speciationRate / (speciationRate + extinctionRate);
  bool isBirth = (unif_rand() < relBr ? true : false);
  spTree->setCurrentTime(eventTime);
  symbiontTree->setCurrentTime(eventTime);
  int numExtantSymbs = symbiontTree->getNumExtant();

  arma::ucolvec cvec = assocMat.col(nodeInd);
  // List rowOfEvents;

  assocMat.shed_col(nodeInd);

  if(isBirth){
    updateEventVector(spTree->getNodesIndxFromExtantIndx(nodeInd),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      3,
                      eventTime);
    spTree->lineageBirthEvent(nodeInd);
    numExtantHosts = spTree->getNumExtant();


    assocMat.resize(numExtantSymbs, numExtantHosts);
    updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      4,
                      eventTime);
    updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      4,
                      eventTime);
    //  assocMat(span::all, numExtantHosts-2) = cvec;
  //  assocMat(span::all, numExtantHosts-1) = cvec;

   // sort symbs on new hosts
   for(int i = 0; i < cvec.n_elem; i++){
      if(cvec[i] == 1){
      // shuffle might be better here?
        arma::umat rr = randi<umat>(1,2, distr_param(0,1));
        if(rr(0,0) == 0 && rr(0,1) == 1){
          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            5,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

        }
        else if(rr(0,0) == 1 && rr(0,1) == 0){

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

            updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            5,
                            eventTime);

        }
        else if(rr(0,0) == 1 && rr(0,1) == 1){
          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);
        }
        else{
          rr.replace(0,1);
          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

          updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                            symbiontTree->getNodesIndxFromExtantIndx(i),
                            4,
                            eventTime);

        }
        assocMat(i, span(numExtantHosts-2,numExtantHosts-1)) = rr;
      }
      else{
        arma::umat rr = zeros<umat>(1,2);
        assocMat(i, span(numExtantHosts-2,numExtantHosts-1)) = rr;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          5,
                          eventTime);
      }
    }
  }
  else{
    updateEventVector(spTree->getNodesIndxFromExtantIndx(nodeInd),
                      symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs - 1),
                      1,
                      eventTime);
    //assocMat.print();
    numExtantSymbs = symbiontTree->getNumExtant();
    // check rows for 0's
    arma::uvec hostless = zeros<uvec>(numExtantSymbs);
    for(int i = assocMat.n_rows - 1; i != -1; i--){
      if(!(any(assocMat.row(i)))){
        updateEventVector(spTree->getNodesIndxFromExtantIndx(nodeInd),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          0,
                          eventTime);
        symbiontTree->lineageDeathEvent(i);

        hostless(i) = 1;
      }
    }
    spTree->lineageDeathEvent(nodeInd);

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

  for(int i = 0; i < numExtantHosts; i++){
    if(sum(assocMat.col(i)) >= 0.5)
      hostIndices.push_back(std::move(i));
  }

  int indxOfHost = unif_rand() * (hostIndices.size() - 1); //col of assocMat

  arma::ucolvec cvec = assocMat.col(hostIndices[indxOfHost]);

  arma::uvec symbIndices = find(cvec);
 // assocMat.print();
  int indxOfSymb = unif_rand() * (symbIndices.size() - 1);

  arma::urowvec rvec = assocMat.row(symbIndices[indxOfSymb]);
  updateEventVector(spTree->getNodesIndxFromExtantIndx(hostIndices[indxOfHost]),
                    symbiontTree->getNodesIndxFromExtantIndx(symbIndices[indxOfSymb]),
                    6,
                    eventTime);
  spTree->lineageBirthEvent(hostIndices[indxOfHost]);
  symbiontTree->lineageBirthEvent(symbIndices[indxOfSymb]);


  numExtantHosts = spTree->getNumExtant();
  int numExtantSymbs = symbiontTree->getNumExtant();

  assocMat.shed_col(hostIndices[indxOfHost]);
  assocMat.shed_row(symbIndices[indxOfSymb]);

  cvec.shed_row(symbIndices[indxOfSymb]);
  rvec.shed_col(hostIndices[indxOfHost]);

  assocMat.resize(numExtantSymbs, numExtantHosts);

  assocMat.submat(numExtantSymbs-2,
                  numExtantHosts-2,
                  numExtantSymbs-1,
                  numExtantHosts-1) = eye<umat>(2,2);
  updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                    symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                    4,
                    eventTime);

  updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                    symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                    4,
                    eventTime);
  // loop through cvec and make little vecs

  for(int i = 0; i < cvec.n_rows; i++){
    if(cvec[i] == 1){

      arma::umat rr = ones<umat>(1,2);
      int randOne = unif_rand() * 2;
      if(randOne == 0){
        rr(0, 0) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          5,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          4,
                          eventTime);
      }
      else{
        rr(0, 1) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          4,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                          symbiontTree->getNodesIndxFromExtantIndx(i),
                          5,
                          eventTime);
      }
      assocMat(i, span(numExtantHosts-2,numExtantHosts-1)) = rr;
    }
    else{
      arma::umat rr = zeros<umat>(1,2);
      assocMat(i, span(numExtantHosts-2,numExtantHosts-1)) = rr;
      updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-2),
                        symbiontTree->getNodesIndxFromExtantIndx(i),
                        5,
                        eventTime);
      updateEventVector(spTree->getNodesIndxFromExtantIndx(numExtantHosts-1),
                        symbiontTree->getNodesIndxFromExtantIndx(i),
                        5,
                        eventTime);
    }
  }



  // loop through rvec
  for(int i = 0; i < rvec.n_elem; i++){
    if(rvec[i] == 1){
      arma::umat rr = ones<umat>(1,2);

      int randOne = unif_rand() * 2;
      if(randOne == 0){
        rr(0,0) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                          5,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                          4,
                          eventTime);
      }
      else{
        rr(0, 1) = 0;
        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-2),
                          4,
                          eventTime);

        updateEventVector(spTree->getNodesIndxFromExtantIndx(i),
                          symbiontTree->getNodesIndxFromExtantIndx(numExtantSymbs-1),
                          5,
                          eventTime);
      }


      assocMat(span(numExtantSymbs-2,numExtantSymbs-1),i) = rr.t();
    }
    else{
      arma::umat rr = zeros<umat>(1,2);
      assocMat(span(numExtantSymbs-2,numExtantSymbs-1), i) = rr.t();
    }
  }
  return assocMat;
}




bool Simulator::bdsaBDSim(){
    bool treesComplete = false;
    double stopTime = spTree->getCurrentTimeFromExtant();
    double eventTime;
    bool isSpeciation;
    lociTree = new LocusTree(numTaxaToSim,
                             currentSimTime,
                             geneBirthRate,
                             geneDeathRate,
                             transferRate);
    std::map<int, int> rToTdckenIndxMap = spTree->makeIndxMap();
    spTree->switchIndicesFirstToSecond(rToTdckenIndxMap);


    std::map<int,double> speciesBirthTimes = spTree->getBirthTimesFromNodes();
    std::map<int,double> speciesDeathTimes = spTree->getDeathTimesFromNodes();
    std::set<int> contempSpecies;
    std::pair<int, int> sibs;

    Node* spRoot = spTree->getRoot();
    lociTree->setStopTime(stopTime);
    currentSimTime -= lociTree->getTimeToNextEvent();
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
        if(geneBirthRate > 0.0 || geneDeathRate > 0.0 || transferRate > 0.0){
            lociTree->ermEvent(currentSimTime);
        }




    }
    lociTree->setPresentTime(currentSimTime);
    treesComplete = true;

    return treesComplete;
}

bool Simulator::simLocusTree(){
  bool good = false;

  while(!good){
    good = bdsaBDSim();
  }
  return good;
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
            geneTree->rootCoalescentProcess(currentSimTime);
            treeGood = true;
        }
        epochCount++;
    }

    //spToLo = lociTree->getLocusToSpeciesMap();
    //geneTree->setIndicesBySpecies(spToLo);
    //geneTree->setBranchLengths();
    return treeGood;
}


bool Simulator::simGeneTree(int j){
  bool gGood = false;
  RNGScope scope;

  while(!gGood){
    gGood = coalescentSim();
  }
  geneTrees[j] = geneTree;
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


double Simulator::getSpeciesTreeRootEdge(){
  return spTree->getRoot()->getDeathTime() - spTree->getRoot()->getBirthTime();
}

double Simulator::getLocusTreeRootEdge(){
  return lociTree->getRoot()->getDeathTime() - lociTree->getRoot()->getBirthTime();
}

double Simulator::getSymbiontTreeRootEdge(){
  return symbiontTree->getRoot()->getDeathTime() - symbiontTree->getRoot()->getBirthTime();
}


double Simulator::getGeneTreeRootEdge(int j){
  return geneTrees[j]->getRoot()->getBranchLength();
}