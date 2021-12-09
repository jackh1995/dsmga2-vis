/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

/* ----------------------------- SECTION IMPORT ----------------------------- */
#include "dsmga2.h"
#include "chromosome.h"
#include "fastcounting.h"
#include "statistics.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <vector>
#include <unordered_set>
using namespace std;
/* -------------------------------------------------------------------------- */

/* ----------------------------- SECTION DSMGA2 ----------------------------- */
DSMGA2::DSMGA2(int _ell, int _nInitial, int _maxGan, int _maxFe, int fffff) :
    // Stack memory allocation
    ell(_ell),
    nCurrent((_nInitial / 2) * 2),
    selectionPressure(2),
    maxGen(_maxGan),
    maxFe(_maxFe),
    repeat((ell > 50) ? ell / 50 : 1),
    generation(0),
    bestIndex(-1),
    // Heap memory allocation
    masks(new vector<int>[ell]),
    selectionIndex(new int[nCurrent]),
    orderN(new int[nCurrent]),
    orderELL(new int[ell]),
    population(new Chromosome[nCurrent]),
    fastCounting(new FastCounting[ell]) {
    // Global initialization
    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0,
    Chromosome::lsnfe = 0,
    Chromosome::rmnfe = 0,
    Chromosome::bmnfe = 0,
    Chromosome::bmeq = 0,
    Chromosome::bmgt = 0,
    Chromosome::bmfail = 0,
    Chromosome::initnfe = 0,
    Chromosome::hit = false,
    // DSM initialization
    graph.init(ell);
    graph_size.init(ell);
    // FastCounting initialization
    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);
    // Population hashing (ch ->  fitness)
    pHash.clear();
    for (int i = 0; i < nCurrent; ++i) {
        population[i].initR(ell);
        double f = population[i].getFitness(Chromosome::initnfe);
        pHash[population[i].getKey()] = f;
    }
    // Greedy hill climbing
    if (GHC) {
        for (int i = 0; i < nCurrent; i++)
            population[i].GHC();
    }
}

DSMGA2::~DSMGA2() {
    delete[] masks;
    delete[] orderN;
    delete[] orderELL;
    delete[] selectionIndex;
    delete[] population;
    delete[] fastCounting;
}

int DSMGA2::doIt () {
    generation = 0;
    while (!shouldTerminate()) {
        if (verbose == POPULATION || verbose == RM) {
            countPop();
            printPop();
        }
        oneRun();
    }
    return generation;
}

void DSMGA2::oneRun() {
    if (CACHE) {Chromosome::cache.clear();}
    mixing();
    double max = -INF;
    stFitness.reset();
    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record(fitness);
    }

    if (verbose != NO) {showStatistics();}
    ++generation;
}

bool DSMGA2::shouldTerminate() {
    bool termination = false;
    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }
    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }
    if (population[0].getMaxFitness() - EPSILON <= stFitness.getMax())
        termination = true;
    if (stFitness.getMax() - EPSILON <= stFitness.getMean())
        termination = true;
    return termination;
}

bool DSMGA2::foundOptima() {
    return (population[0].getMaxFitness() - EPSILON <= stFitness.getMax());
}

void DSMGA2::showStatistics() {
    printf("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n",
           generation, stFitness.getMax(), stFitness.getMean(),
           stFitness.getMin());
    fflush(NULL);
}

void DSMGA2::buildFastCounting() {
    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++) {
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }
        }
    } else {
        for (int i = 0; i < nCurrent; i++) {
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
        }
    }
}

int DSMGA2::countOne(int x) const {
    int n = 0;
    for (int i = 0; i < fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;
        val = fastCounting[x].gene[i];
        n += myBD.countOne(val);
    }
    return n;
}

int DSMGA2::countXOR(int x, int y) const {
    int n = 0;
    for (int i = 0; i < fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;
        val = fastCounting[x].gene[i];
        val ^= fastCounting[y].gene[i];
        n += myBD.countOne(val);
    }
    return n;
}

void DSMGA2::findMaskWithBound(Chromosome &ch, vector<int> &result, int startNode, int bound) {
    result.clear();
    DLLA rest(ell);
    for (int i = 0; i < ell; i++) {
        if (orderELL[i] == startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }
    double *connection = new double[ell];
    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
        pair<double, double> p = graph_size(startNode, *iter);
        int i = ch.getVal(startNode);
        int j = ch.getVal(*iter);
        // p00 or p11
        if (i == j) 
            connection[*iter] = p.first;
        // p01 or p10
        else 
            connection[*iter] = p.second;
    }
    bound--;
    while (!rest.isEmpty() && bound > 0) {
        bound--;
        double max = -INF;
        int index = -1;
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }
        rest.erase(index);
        result.push_back(index);
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
            pair<double, double> p = graph_size(index, *iter);
            int i = ch.getVal(index);
            int j = ch.getVal(*iter);
            // p00 or p11
            if (i == j) 
                connection[*iter] += p.first;
            // p01 or p10
            else 
                connection[*iter] += p.second;
        }
    }
    delete[] connection;
}

void DSMGA2::findMask(Chromosome &ch, vector<int> &result, int startNode) {
    result.clear();
    DLLA rest(ell);
    genOrderELL();
    for (int i = 0; i < ell; i++) {
        if (orderELL[i] == startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }
    double *connection = new double[ell];
    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
        pair<double, double> p = graph(startNode, *iter);
        int i = ch.getVal(startNode);
        int j = ch.getVal(*iter);
        // p00 or p11
        if (i == j) 
            connection[*iter] = p.first;
        // p01 or p10
        else 
            connection[*iter] = p.second;
    }
    while (!rest.isEmpty()) {
        double max = -INF;
        int index = -1;
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }
        rest.erase(index);
        result.push_back(index);
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
            pair<double, double> p = graph(index, *iter);
            int i = ch.getVal(index);
            int j = ch.getVal(*iter);
            // p00 or p11
            if (i == j) 
                connection[*iter] += p.first;
            // p01 or p10
            else 
                connection[*iter] += p.second;
        }
    }
    delete[] connection;
}

size_t DSMGA2::findSize(Chromosome &ch, vector<int> &mask) const {
    DLLA candidate(nCurrent);
    for (int i = 0; i < nCurrent; ++i)
        candidate.insert(i);
    size_t size = 0;
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        int allele = ch.getVal(*it);
        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);
            if (candidate.isEmpty())
                break;
        }
        if (candidate.isEmpty())
            break;
        ++size;
    }
    return size;
}

void DSMGA2::restrictedMixing(Chromosome &ch) {
    
    int startNode = myRand.uniformInt(0, ell - 1);
    // size_2e
    vector<int> mask;
    findMask(ch, mask, startNode);
    size_t size = findSize(ch, mask);
    // size_1e
    vector<int> mask_size;
    findMaskWithBound(ch, mask_size, startNode, size);
    size_t size_original = findSize(ch, mask_size);
    // Trim
    if (size > size_original) {size = size_original;}
    if (size > (size_t)ell / 2) { size = ell / 2;}
    while (mask.size() > size) {mask.pop_back();}
    // Restricted mixing
    bool taken = restrictedMixing(ch, mask);
    // Back mixing
    EQ = true;
    if (taken) {
        // FIXME
        Chromosome::bmeq=0;
        Chromosome::bmgt=0;
        Chromosome::bmfail=0;
        genOrderN();
        for (int i = 0; i < nCurrent; ++i) {
            if (EQ) {backMixingE(ch, mask, population[orderN[i]]);}
            else {backMixing(ch, mask, population[orderN[i]]);}
        }
        
        if (verbose == RM) {
            cout << ' ' << GREEN << Chromosome::bmgt << RESET << '/' \
            << YELLOW << Chromosome::bmeq << RESET << '/'  << \
            RED << Chromosome::bmfail << RESET << endl;
        }
        
    }
}

void DSMGA2::backMixing(Chromosome &source, vector<int> &mask, Chromosome &des) {
    Chromosome trial(ell);
    trial = des;
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));
    if (trial.getFitness(Chromosome::bmnfe) > des.getFitness(Chromosome::bmnfe)) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        ++Chromosome::bmgt;
        return;
    }
    ++Chromosome::bmfail;
}

void DSMGA2::backMixingE(Chromosome &source, vector<int> &mask, Chromosome &des) {
    Chromosome trial(ell);
    trial = des;
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));
    // Strong bm
    if (trial.getFitness(Chromosome::bmnfe) > des.getFitness(Chromosome::bmnfe)) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        EQ = false;
        des = trial;
        ++Chromosome::bmgt;
        return;
    }
    // Weak bm
    if (trial.getFitness(Chromosome::bmnfe) >= des.getFitness(Chromosome::bmnfe) - EPSILON) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        ++Chromosome::bmeq;
        return;
    }
    ++Chromosome::bmfail;
}

bool DSMGA2::restrictedMixing(Chromosome &ch, vector<int> &mask) {
    bool taken(false);
    size_t lastUB(0);
    Chromosome trial(ell);
    trial = ch;
    
    for (size_t ub = 1; ub <= mask.size(); ++ub) {
        size_t size = 1;
        Chromosome trial(ell);
        trial = ch;
	    vector<int> takenMask;
        for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
			takenMask.push_back(*it);
            trial.flip(*it);
            ++size;
            if (size > ub) break;
        }
        if (isInP(trial)) break;
        if (trial.getFitness(Chromosome::rmnfe) >= ch.getFitness(Chromosome::rmnfe) - EPSILON) {
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();
            taken = true;
            ch = trial;
        }
        if (taken) {
            lastUB = ub;
            break;
        }
    }
    
    if (verbose == RM)
        if (mask.size()) {
            unordered_set<int> success;
            unordered_set<int> fail;
            for (size_t idx=0; idx!=mask.size(); ++idx) {
                if (idx < lastUB)
                    success.insert(mask[idx]);
                else
                    fail.insert(mask[idx]);
            }
            for (size_t idx=0; idx!=size_t(ell); ++idx) {
                if (success.count(idx))
                    cout << GREEN << idx%10 << RESET;
                else if (fail.count(idx))
                    cout << RED << idx%10 << RESET;
                else
                    cout << idx%10;
            } 
            cout << ' ' << lastUB << '/' << mask.size();
            
            if (!taken)
                cout << endl;
        }
    
    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }
    return taken;
}

void DSMGA2::mixing() {
    
    if (SELECTION) selection();
    buildFastCounting();
    buildGraph();
    buildGraph_sizecheck();

    //for (int i=0; i<ell; ++i)
    //    findClique(i, masks[i]); // replaced by findBoundedMask in restrictedMixing
    for (int k = 0; k < repeat; ++k) {
        genOrderN();
        for (int i = 0; i < nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
            if (Chromosome::hit)
                break;
        }
        if (Chromosome::hit)
            break;
    }
}

inline bool DSMGA2::isInP(const Chromosome &ch) const {
    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent - 1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell - 1);
}

void DSMGA2::buildGraph() {
    int *one = new int[ell];
    for (int i = 0; i < ell; ++i) {
        one[i] = countOne(i);
    }
    for (int i = 0; i < ell; ++i) {
        for (int j = i + 1; j < ell; ++j) {
            int n00, n01, n10, n11;
            int nX = countXOR(i, j);
            n11 = (one[i] + one[j] - nX) / 2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;
            double p00 = (double)n00 / (double)nCurrent;
            double p01 = (double)n01 / (double)nCurrent;
            double p10 = (double)n10 / (double)nCurrent;
            double p11 = (double)n11 / (double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;
            double linkage = computeMI(p00, p01, p10, p11);
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00 * log(p00 / p_0 / p0_);
            if (p11 > EPSILON)
                linkage00 += p11 * log(p11 / p_1 / p1_);
            if (p01 > EPSILON)
                linkage01 += p01 * log(p01 / p0_ / p_1);
            if (p10 > EPSILON)
                linkage01 += p10 * log(p10 / p1_ / p_0);
            if (Chromosome::nfe < 0) {
                pair<double, double> p(linkage, linkage);
                graph.write(i, j, p);
            } else {
                pair<double, double> p(linkage00, linkage01);
                graph.write(i, j, p);
            }
        }
    }
    delete[] one;
}

void DSMGA2::buildGraph_sizecheck() {
    int *one = new int[ell];
    for (int i = 0; i < ell; ++i) {
        one[i] = countOne(i);
    }
    for (int i = 0; i < ell; ++i) {
        for (int j = i + 1; j < ell; ++j) {
            int n00, n01, n10, n11;
            int nX = countXOR(i, j);
            n11 = (one[i] + one[j] - nX) / 2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;
            double p00 = (double)n00 / (double)nCurrent;
            double p01 = (double)n01 / (double)nCurrent;
            double p10 = (double)n10 / (double)nCurrent;
            double p11 = (double)n11 / (double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;
            double linkage = computeMI(p00, p01, p10, p11);
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00 * log(p00 / p_0 / p0_);
            if (p11 > EPSILON)
                linkage00 += p11 * log(p11 / p_1 / p1_);
            if (p01 > EPSILON)
                linkage01 += p01 * log(p01 / p0_ / p_1);
            if (p10 > EPSILON)
                linkage01 += p10 * log(p10 / p1_ / p_0);
            pair<double, double> p(linkage, linkage);
            graph_size.write(i, j, p);
        }
    }

    delete[] one;
}

double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {
    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1 - p0;
    double q1 = 1 - q0;
    double join = 0.0;
    if (a00 > EPSILON)
        join += a00 * log(a00);
    if (a01 > EPSILON)
        join += a01 * log(a01);
    if (a10 > EPSILON)
        join += a10 * log(a10);
    if (a11 > EPSILON)
        join += a11 * log(a11);
    double p = 0.0;
    if (p0 > EPSILON)
        p += p0 * log(p0);
    if (p1 > EPSILON)
        p += p1 * log(p1);
    double q = 0.0;
    if (q0 > EPSILON)
        q += q0 * log(q0);
    if (q1 > EPSILON)
        q += q1 * log(q1);
    return -p - q + join;
}

void DSMGA2::selection() {
    tournamentSelection();
}

void DSMGA2::tournamentSelection() {
    int i, j;
    int randArray[selectionPressure * nCurrent];
    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray(randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);
    for (i = 0; i < nCurrent; i++) {
        int winner = 0;
        double winnerFitness = -INF;
        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();
            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }
        }
        selectionIndex[i] = winner;
    }
}

void DSMGA2::countPop(void) {
    ch_count.clear();
    for (int idx=0; idx!=nCurrent; ++idx) {
        // has seen
        if (ch_count.count(population[idx].getKey())) {
            ++ch_count[population[idx].getKey()];
        }
        else {
            ch_count[population[idx].getKey()] = 1;
        }
    } cout << endl;
}

void DSMGA2::printPop(void) {
    population_set.clear();
    for (int idx=0; idx!=nCurrent; ++idx) {
        population_set.insert(population[idx]);
    }
    set<Chromosome>::iterator it;
    for (it = population_set.begin(); it != population_set.end(); ++it) {
        cout << *it << ' ' << it->getCount() << endl;
    }
}
/* -------------------------------------------------------------------------- */

const Chromosome* const DSMGA2::getPopulation(void) const {
    return population;
}

int const DSMGA2::getBestIndex(void) const {
    return bestIndex;
}