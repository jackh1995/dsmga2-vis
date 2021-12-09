/*
 * dsmga2.h
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */

#ifndef _DSMGA2_H_
#define _DSMGA2_H_

/* ----------------------------- SECTION IMPORT ----------------------------- */
#include "chromosome.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"
#include "statistics.h"
#include "trimatrix.h"
#include <list>
#include <unordered_map>
#include <set>
/* -------------------------------------------------------------------------- */

class DSMGA2 {
public:
    DSMGA2(int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff);
    ~DSMGA2();
    void oneRun();
    int doIt();
    // DSM
    void selection();
    void tournamentSelection();
    void buildGraph();
    // Mixing
    void mixing();
    void findMask(Chromosome &ch, vector<int> &mask, int startNode);
    void findMaskWithBound(Chromosome &ch, vector<int> &mask, int startNode, int bound);
    void buildGraph_sizecheck();
    void restrictedMixing(Chromosome &);
    bool restrictedMixing(Chromosome &ch, vector<int> &mask);
    void backMixing(Chromosome &source, vector<int> &mask, Chromosome &des);
    void backMixingE(Chromosome &source, vector<int> &mask, Chromosome &des);
    // Termination
    bool shouldTerminate();
    bool foundOptima();
    int getGeneration() const { return generation; }
    bool isInP(const Chromosome &) const;
    // Traversal
    void genOrderN();
    void genOrderELL();
    void showStatistics();
    bool isSteadyState();
    static unordered_map<unsigned long, unsigned> ch_count;
    static set<Chromosome> population_set;
    const Chromosome* const getPopulation(void) const;
    int const getBestIndex(void) const;

private:
    // Population
    int ell;
    int nCurrent; 
    int selectionPressure;
    // Termination
    int maxGen, maxFe, repeat;
    // Statistics
    int generation;
    int bestIndex;
    Statistics stFitness;
    // Population hashing (ch ->  fitness)
    unordered_map<unsigned long, double> pHash;
    // Heap memory allocation
    vector<int> *masks;
    int *selectionIndex;
    int *orderN; 
    int *orderELL;
    Chromosome *population;
    FastCounting *fastCounting;
    // DSM
    TriMatrix<double> graph;
    TriMatrix<double> graph_size;
    // DSM
    double computeMI(double, double, double, double) const;
    void buildFastCounting();
    int countXOR(int, int) const;
    int countOne(int) const;
    // Mixing
    bool EQ;
    size_t findSize(Chromosome &, vector<int> &) const;
    size_t findSize(Chromosome &, vector<int> &, Chromosome &) const;
    // Print
    void countPop(void);
    void printPop(void);
};

#endif