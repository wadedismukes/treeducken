//
//  Engine.cpp
//  multiTree
//
//  Created by Dismukes, Wade T [EEOBS] on 1/28/18.
//  Copyright © 2018 Dismukes, Wade T [EEOBS]. All rights reserved.
//

#include "Engine.h"
#include "Rcpp.h"

using namespace Rcpp;


Engine::Engine(std::string of,
                 int mt,
                 double sbr,
                 double sdr,
                 double gbr,
                 double gdr,
                 double lgtr,
                 int ipp,
                 int popsize,
                 double genTime,
                 double ts,
                 int reps,
                 int ntax,
                 int nloci,
                 int ngen,
                 double og,
                 bool sout,
                 bool mst){
    outfilename = of;
    inputSpTree = "";
    simType = mt;
    spBirthRate = sbr;
    spDeathRate = sdr;
    geneBirthRate = gbr;
    geneDeathRate = gdr;
    transferRate = lgtr;
    doScaleTree = false;
    treescale = ts;
    seedset = 0;
    printOutputToScreen = sout;
    individidualsPerPop = ipp;
    populationSize = popsize;
    generationTime = genTime;
    numSpeciesTrees = reps;
    proportionToSample = 1.0;
    numTaxa = ntax;
    numLoci = nloci;
    numGenes = ngen;
    outgroupFrac = og;
    // if(sd1 > 0 && sd2 > 0)
    //     rando.setSeed(sd1, sd2);
    // else
    //     rando.setSeed();
    // seedType gs1, gs2;
    // rando.getSeed(gs1, gs2);
    //Rcout << "\nSeeds = {" << gs1 << ", " << gs2 << "}" << std::endl;
}


Engine::~Engine(){
    for(std::vector<TreeInfo*>::iterator p=simSpeciesTrees.begin(); p != simSpeciesTrees.end(); ++p){
        delete (*p);
    }
    simSpeciesTrees.clear();
}


void Engine::doRunRun(){
    // double TS = 0.0;
    TreeInfo *ti = nullptr;
    for(int k = 0; k < numSpeciesTrees; k++){

        Simulator *treesim = new Simulator(numTaxa,
                                           spBirthRate,
                                           spDeathRate,
                                           1.0,
                                           numLoci,
                                           geneBirthRate,
                                           geneDeathRate,
                                           transferRate,
                                           individidualsPerPop,
                                           populationSize,
                                           generationTime,
                                           numGenes,
                                           outgroupFrac,
                                           treescale,
                                           printOutputToScreen);
        if(printOutputToScreen)
            Rcout << "Simulating species tree replicate # " << k + 1 << std::endl;

        switch(simType){
            case 1:
                treesim->simSpeciesTree();
                break;
            case 2:
                treesim->simSpeciesLociTrees();
                break;
            case 3:
                treesim->simThreeTree();
                break;
            case 4:
                treesim->simLocusGeneTrees();
            default:
                treesim->simSpeciesTree();
                break;
        }
        // TODO: find seg fault between here and Rcout statement
        ti =  new TreeInfo(k, numLoci);
        ti->setWholeTreeStringInfo(treesim->printSpeciesTreeNewick());
        ti->setExtTreeStringInfo(treesim->printExtSpeciesTreeNewick());
        ti->setSpeciesTreeDepth(treesim->calcSpeciesTreeDepth());
        ti->setExtSpeciesTreeDepth(treesim->calcExtantSpeciesTreeDepth());
        ti->setNumberTransfers(treesim->findNumberTransfers());
        for(int i = 0; i < numLoci; i++){
            ti->setLocusTreeByIndx(i, treesim->printLocusTreeNewick(i));

            if(simType == 3){
                for(int j = 0; j < numGenes; j++){
                    // Rcout << "gene index " << j << "?????????" << std::endl;
                    // Rcout << treesim->printGeneTreeNewick(i, j) << std::endl;
                    // Rcout << treesim->printExtantGeneTreeNewick(i, j) << std::endl;

                    ti->setGeneTreeByIndx(i, j, treesim->printGeneTreeNewick(i, j));
                    ti->setExtantGeneTreeByIndx(i, j, treesim->printExtantGeneTreeNewick(i, j));
                }
            }
        }
        simSpeciesTrees.push_back(ti);
        delete treesim;
        treesim = nullptr;
    }
    this->writeTreeFiles();
    return;
}


void Engine::writeTreeFiles(){

    for(std::vector<TreeInfo *>::iterator p = simSpeciesTrees.begin(); p != simSpeciesTrees.end(); p++){
        int d = (int) std::distance(simSpeciesTrees.begin(), p);
        (*p)->writeTreeStatsFile(d, outfilename);
        (*p)->writeWholeTreeFileInfo(d, outfilename);
        (*p)->writeExtantTreeFileInfo(d, outfilename);
        for(int i = 0; i < numLoci; i++){
            (*p)->writeLocusTreeFileInfoByIndx(d, i, outfilename);
            if(simType == 3)
                // for(int j = 0; j < numGenes; j++){
                //     (*p)->writeGeneTreeFileInfoByIndx(d, i, j, outfilename);
                // }
                (*p)->writeExtGeneTreeFileInfo(d, i, numGenes, outfilename);
        }
    }
}


TreeInfo* Engine::findTreeByIndx(int i){
    TreeInfo *tf = 0;
    int count = 0;
    for(std::vector<TreeInfo*>::iterator it = simSpeciesTrees.begin(); it != simSpeciesTrees.end(); ++it){
        if(count == i){
            tf = (*it);
            break;
        }
        else
            count++;
    }
    return tf;
}


void Engine::calcAverageRootAgeSpeciesTrees(){
    std::ofstream out;
    out.open("Average_root_depths_spTree.out");
    double sumRH = 0.0;
    for(std::vector<TreeInfo*>::iterator p = simSpeciesTrees.begin(); p != simSpeciesTrees.end(); ++p){
        sumRH += (*p)->getSpeciesTreeDepth();
        out << (*p)->getSpeciesTreeDepth() << "\n";
    }
    out.close();

    Rcout << "\n #### Average Root Age of Species Trees is " << (double) sumRH / numSpeciesTrees << " #####" << std::endl;


}


// these newick reading functions are heavily modified by those used by Paul Lewis
// https://phylogeny.uconn.edu/phylogenetic-software-development-tutorial/build-a-tree-from-a-newick-description/#
unsigned int Engine::countNewickLeaves(const std::string spTreeStr){
    std::regex taxonregexpr ("(\\w*\\.?[\\w|\\s|\\.]?\\w*)\\:");
    std::sregex_iterator it1(spTreeStr.begin(), spTreeStr.end(), taxonregexpr);
    std::sregex_iterator it2;
    return (unsigned) std::distance(it1, it2);
}

std::string Engine::stripCommentsFromNewickTree(std::string spTreeStr){
    std::string commentlessNewick;
    std::regex commentregexpr ("\\[.*?\\]");
    commentlessNewick = std::regex_replace(spTreeStr,commentregexpr,std::string(""));
    return commentlessNewick;
}

std::string Engine::formatTipNamesFromNewickTree(std::string spTreeStr){
    std::string formattedNewick;
    std::regex taxonregexpr ("(\\w*\\.?[\\w|\\s|\\.]?\\w*)\\:");

    formattedNewick = std::regex_replace(spTreeStr,taxonregexpr, "'$1':");
    return formattedNewick;
}

SpeciesTree* Engine::buildTreeFromNewick(const std::string spTreeStr){
    SpeciesTree* spTree = nullptr;
    Node* currNode = nullptr;
    Node* prevNode = nullptr;
    std::string commentlessSpTreeStr = spTreeStr;
    commentlessSpTreeStr = stripCommentsFromNewickTree(commentlessSpTreeStr);
    numTaxa = countNewickLeaves(commentlessSpTreeStr);
    commentlessSpTreeStr = formatTipNamesFromNewickTree(commentlessSpTreeStr);
    spTree = new SpeciesTree(numTaxa);
    currNode = new Node();
    spTree->setRoot(currNode);
    currNode->setAsRoot(true);

    enum {
        Prev_Tok_LParen		= 0x01,	// previous token was a left parenthesis ('(')
        Prev_Tok_RParen		= 0x02,	// previous token was a right parenthesis (')')
        Prev_Tok_Colon		= 0x04,	// previous token was a colon (':')
        Prev_Tok_Comma		= 0x08,	// previous token was a comma (',')
        Prev_Tok_Name		= 0x10,	// previous token was a node name (e.g. '2', 'P._articulata')
        Prev_Tok_EdgeLen	= 0x20	// previous token was an edge length (e.g. '0.1', '1.7e-3')
    };
    unsigned previous = Prev_Tok_LParen;

    // Some useful flag combinations
    unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
    unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    unsigned Comma_Valid  = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    unsigned Colon_Valid  = (Prev_Tok_RParen | Prev_Tok_Name);
    unsigned Name_Valid   = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);
    std::string::const_iterator newickStart = commentlessSpTreeStr.begin();
    std::string::const_iterator it = newickStart;

    for(; it != commentlessSpTreeStr.end(); ++it){
        char ch = (*it);
        if(iswspace(ch))
            continue;
        switch(ch){
            case ';': {
                Rcout << "Species tree read in successfully.\n";
                break;
            }
            case '(': {
                if(!(previous & LParen_Valid)){
                    Rcerr << "Your newick tree is not formatted properly. Exiting...\n";
                    Rcerr << "A left parenthetical is not in the right place maybe...\n";
                    exit(1);
                }
                prevNode = currNode;
                currNode =  new Node();
                currNode->setAnc(prevNode);
                prevNode->setLdes(currNode);
                previous = Prev_Tok_LParen;
                break;
            }
            case ':': {
                if(!(previous & Colon_Valid)){
                    Rcerr << "Your newick tree is not formatted properly. Exiting...\n";
                    Rcerr << "A colon appears to be in the wrong place...\n";
                    exit(1);
                }
                previous = Prev_Tok_Colon;
                break;
            }
            case ')': {
                if(!(previous & RParen_Valid)){
                    Rcerr << "Your newick tree is not formatted properly. Exiting...\n";
                    Rcerr << "A right parenthetical is not in the right place maybe...\n";
                    exit(1);
                }
                //prevNode = currNode;
                prevNode = nullptr;
                currNode = currNode->getAnc();
                previous = Prev_Tok_RParen;
                break;
            }
            case ',': {
                if(!(previous & Comma_Valid)){
                    Rcerr << "Your newick tree is not formatted properly. Exiting...\n";
                    Rcerr << "A comma is not in the right place maybe...\n";
                    exit(1);
                }
                prevNode = currNode;
                currNode = new Node();
                prevNode->setSib(currNode);
                currNode->setSib(prevNode);
                currNode->setAnc(prevNode->getAnc());
                prevNode->getAnc()->setRdes(currNode);

                previous = Prev_Tok_Comma;
                break;
            }
            case '\'': {
                std::string tipname = "";
                for (++it; it != commentlessSpTreeStr.end(); ++it){
                    ch = *it;
                    if (ch == '\''){
                        break;
                    }
                    else if (iswspace(ch))
                        tipname += ' ';
                    else
                        tipname += ch;
                }
                if(previous == Prev_Tok_RParen){
                    //prevNode->setIsTip(true);
                //    prevNode->setName(tipname);
                }
                else if(previous == Prev_Tok_Comma){
                    currNode->setIsTip(true);
                    currNode->setName(tipname);
                }
                else if(previous == Prev_Tok_LParen){
                    currNode->setName(tipname);
                    currNode->setIsTip(true);
                }

                previous = Prev_Tok_Name;
                break;
            }
            default: {
                // branch length
                if(previous == Prev_Tok_Colon){
                    std::string::const_iterator jit = it;
                    for (; it != commentlessSpTreeStr.end(); ++it){
                        ch = *it;
                        if (ch == ',' || ch == ')' || iswspace(ch)){
                            --it;
                            break;
                        }
                        bool valid = (ch =='e' || ch == 'E' || ch =='.' || ch == '-' || ch == '+' || isdigit(ch));
                        if (!valid){
                            Rcerr << "Invalid branch length character in tree description\n";
                            exit(1);
                        }
                        std::string edge_length_str = std::string(jit,it+1);
                        currNode->setBranchLength(atof(edge_length_str.c_str()));
                        if (currNode->getBranchLength() < 1.e-10)
                            currNode->setBranchLength(1.e-10);
                        previous = Prev_Tok_EdgeLen;
                    }
                }
                else{
                    // Get the node name
                    std::string tipname = "";
                    for (; it != commentlessSpTreeStr.end(); ++it){
                        ch = *it;
                        if (ch == '('){
                            Rcerr << "Didn't expect a left parenthesis here! Check your newick tree...\n";
                            Rcerr << "Exiting... :(\n";
                            exit(1);
                        }
                        if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')'){
                            --it;
                            break;
                        }
                        tipname += ch;
                    }

                    // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                    if (!(previous & Name_Valid)){
                       stop("Unexpected placement of name of tip. Exiting...\n");
                    }
                    currNode->setIsTip(true);
                    currNode->setName(tipname);
                    previous = Prev_Tok_Name;
                }
                if(it == commentlessSpTreeStr.end())
                    break;
            }

        }

    }
    spTree->popNodes();
    spTree->setTreeInfo();
    return spTree;
}

void Engine::doRunSpTreeSet(){

    Rcout << "Setting species tree to this newick tree: " << inputSpTree << std::endl;

    TreeInfo *ti = nullptr;
    Simulator *treesim = new Simulator(numTaxa,
                                        spBirthRate,
                                        spDeathRate,
                                        1.0,
                                        numLoci,
                                        geneBirthRate,
                                        geneDeathRate,
                                        transferRate,
                                        individidualsPerPop,
                                        populationSize,
                                        generationTime,
                                        numGenes,
                                        outgroupFrac,
                                        treescale,
                                        printOutputToScreen);


    treesim->setSpeciesTree(this->buildTreeFromNewick(inputSpTree));
    treesim->simLocusGeneTrees();


    ti =  new TreeInfo(0, numLoci);
    ti->setWholeTreeStringInfo(treesim->printSpeciesTreeNewick());
    ti->setSpeciesTreeDepth(treesim->calcSpeciesTreeDepth());
    ti->setExtSpeciesTreeDepth(treesim->calcExtantSpeciesTreeDepth());
    ti->setNumberTransfers(treesim->findNumberTransfers());
    for(int i = 0; i < numLoci; i++){
        ti->setLocusTreeByIndx(i, treesim->printLocusTreeNewick(i));
        if(simType == 3){
            for(int j = 0; j < numGenes; j++){
                ti->setGeneTreeByIndx(i, j, treesim->printGeneTreeNewick(i, j));
                ti->setExtantGeneTreeByIndx(i, j, treesim->printExtantGeneTreeNewick(i, j));

            }
        }
    }
    simSpeciesTrees.push_back(ti);
    delete treesim;
    treesim = nullptr;

    this->writeTreeFiles();

}


/*
    TreeInfo functions to write tree information to file in various file formats
                                                                */
TreeInfo::TreeInfo(int idx, int nl){
    locusTrees.resize(nl);
    geneTrees.resize(nl);
    extGeneTrees.resize(nl);
    spTreeLength = 0.0;
    extSpTreeLength = 0.0;
    spTreeDepth = 0.0;
    extSpTreeDepth = 0.0;
    spTreeNess = 0.0;
    spAveTipLen = 0.0;
    loTreeLength = 0.0;
    loTreeDepth = 0.0;
    loTreeNess = 0.0;
    loAveTipLen = 0.0;
    aveTMRCAGeneTree = 0.0;
    numTransfers = 0;
}

TreeInfo::~TreeInfo(){
    gsaTrees.clear();
    locusTrees.clear();
    geneTrees.clear();
    extGeneTrees.clear();
    speciesTree.clear();
}

void TreeInfo::writeTreeStatsFile(int spIndx, std::string ofp){
    std::string path = "";
    std::string fn = ofp;
    std::stringstream tn;
    tn << spIndx;
    fn += "_" + tn.str() + ".sp.tre.stats.txt";
    path += fn;
    std::ofstream out(path);
    out << "Species Tree Statistics\n";

    tn.clear();
    tn.str(std::string());
    tn << getSpeciesTreeDepth();
    out << "Tree depth\t" << tn.str() << std::endl;

    tn.clear();
    tn.str(std::string());

    tn << getExtSpeciesTreeDepth();
    out << "Extant Tree depth\t" << tn.str() << std::endl;

    tn.clear();
    tn.str(std::string());
    // tn << getSpeciesAveTipLen();
    // out << "Average branch length: " << tn.str() << std::endl;
    tn << getNumberTransfers();
    out << "Transfers\t" << tn.str() << std::endl;

    tn.clear();
    tn.str(std::string());
    //TODO: add the statistics for locus and gene trees
}

void TreeInfo::writeWholeTreeFileInfo(int spIndx, std::string ofp){
    std::string path = "";

    std::string fn = ofp;
    std::stringstream tn;

    tn << spIndx;
    //path += tn.str();
    // path +=  "/";


    fn += "_" + tn.str() + ".sp.full.tre";
    path += fn;
    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree wholeT_" << spIndx << " = ";
    out << getWholeSpeciesTree() << "\n";
    out << "end;";

}

void TreeInfo::writeExtantTreeFileInfo(int spIndx, std::string ofp){
    std::string path = "";

    std::string fn = ofp;
    std::stringstream tn;

    tn << spIndx;
    //path += tn.str();
    // path +=  "/";


    fn += "_" + tn.str() + ".sp.tre";
    path += fn;
    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree extT_" << spIndx << " = ";
    out << getExtantSpeciesTree() << "\n";
    out << "end;";
}

void TreeInfo::writeLocusTreeFileInfoByIndx(int spIndx, int indx, std::string ofp){
    std::string path = "";

    std::string fn = ofp;
    std::stringstream tn;

    tn << spIndx;
    // path += tn.str();
    // path += "/";


    fn += "_" + tn.str();

    tn.clear();
    tn.str(std::string());

    tn << indx;

    fn += "_" + tn.str() + ".loc.tre";
    path += fn;

    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree locT_" << indx << " = ";
    out << getLocusTreeByIndx(indx) << "\n";
    out << "end;";
}

void TreeInfo::writeGeneTreeFileInfoByIndx(int spIndx, int Lindx, int indx, std::string ofp){
    std::string path = "";

    std::string fn = ofp;
    std::stringstream tn;

    tn << spIndx;
    //path += tn.str();
   // path += "/";


    fn += "_" + tn.str();

    tn.clear();
    tn.str(std::string());

    tn << Lindx;

    fn += "_" + tn.str();
    tn.clear();
    tn.str(std::string());

    tn << indx;
    fn += "_" + tn.str() + ".gen.tre";
    path += fn;

    std::ofstream out(path);
    out << "#NEXUS\nbegin trees;\n    tree geneT_" << indx << " = ";
    out << getGeneTreeByIndx(Lindx, indx) << "\n";
    out << "tree extGeneT_" << indx << " = ";
    out << getExtGeneTreeByIndx(Lindx, indx) << "\n";
    out << "end;";

}

void TreeInfo::writeExtGeneTreeFileInfo(int spIndx, int Lindx, int numGenes, std::string ofp){
    std::string path = "";

    std::string fn = ofp;
    std::stringstream tn;

    tn << spIndx;
    //path += tn.str();
    // path += "/";
    fn += "genetrees_";

    fn += tn.str();

    tn.clear();
    tn.str(std::string());

    tn << Lindx;

    fn += "_" + tn.str() + ".tre";;
    tn.clear();
    tn.str(std::string());

    // tn << i;
    // fn += "_" + tn.str() + ".gen.tre";
    path += fn;

    std::ofstream out(path);
    // out << "#NEXUS\nbegin trees;\n";
    // for(int i = 0; i < numGenes; i++){
    //         out << "tree geneT_" << indx << " = ";
    //         out << getGeneTreeByIndx(Lindx, indx) << "\n";
    // }
    // out << "end;";


    out << "#NEXUS\nbegin trees;\n";
    for(int i = 0; i < numGenes; i++){
        out << "tree extGeneT_" << i << " = ";
        out << getExtGeneTreeByIndx(Lindx, i) << "\n";
    }
    out << "end;";

}

void TreeInfo::writeGeneTreeFileInfo(int spIndx, int Lindx, int numGenes, std::string ofp){
    std::string path = "";

    std::string fn = ofp;
    std::stringstream tn;

    tn << spIndx;
    //path += tn.str();
    // path += "/";
    fn += "genetrees_";

    fn += tn.str();

    tn.clear();
    tn.str(std::string());

    tn << Lindx;

    fn += "_" + tn.str() + ".tre";;
    tn.clear();
    tn.str(std::string());

    // tn << i;
    // fn += "_" + tn.str() + ".gen.tre";
    path += fn;

    std::ofstream out(path);
    // out << "#NEXUS\nbegin trees;\n";
    // for(int i = 0; i < numGenes; i++){
    //         out << "tree geneT_" << indx << " = ";
    //         out << getGeneTreeByIndx(Lindx, indx) << "\n";
    // }
    // out << "end;";


    out << "#NEXUS\nbegin trees;\n";
    for(int i = 0; i < numGenes; i++){
        out << "tree geneT_" << i << " = ";
        out << getExtGeneTreeByIndx(Lindx, i) << "\n";
    }
    out << "end;";

}
