/***************************************************************************
 *   Copyright (C) 2015 by TEIL                                            *
 ***************************************************************************/

#include <cstdio>
#include <cstring>
#include "spin.h"
#include "chromosome.h"
#include "nk-wa.h"
#include "sat.h"
#include <iostream>
#include <bitset>
#include "dsmga2.h"
using namespace std;

#define TRAP_K 5

Chromosome::Chromosome () {
    length = 0;
    lengthLong = 0;
    gene = NULL;
    evaluated = false;
}

Chromosome::Chromosome (int n_length) {
    gene = NULL;
    init (n_length);
}

// FIXMe
Chromosome::Chromosome(const Chromosome &c) {
    length = c.length;
    evaluated = c.evaluated;
    fitness = c.fitness;
    lengthLong = c.lengthLong;
    key = c.key;
    gene = new unsigned long[c.lengthLong];
    memcpy(gene, c.gene, sizeof(long) * lengthLong);
}

Chromosome::~Chromosome () {
    if (gene != NULL) delete []gene;
}

void Chromosome::init (int _length) {
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];
    gene[lengthLong-1] = 0;

    evaluated = false;
}

void Chromosome::init0 (int _length) {
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];

    for (int i=0; i<lengthLong; ++i)
        gene[i] = 0;

    key = 0;
    evaluated = false;
}

void Chromosome::initR (int _length) {
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];
    gene[lengthLong-1] = 0;

    key = 0;
    for (int i=0; i<length; ++i) {

        int val = myRand.flip();
        setValF(i, val);
        if (val == 1)
            key ^= zKey[i];
    }

    evaluated = false;
}

double Chromosome::getFitness (int& counter) {
    if (evaluated)
        return fitness;
    else {
        fitness = evaluate(counter);
        if (!hit && fitness > getMaxFitness()) {
            hit = true;
        }
        return fitness;
    }
}

double Chromosome::getFitness () const {
    if (evaluated)
        return fitness;
    else {
        throw runtime_error("Error: Chromosome not evaluated");
    }
}

int Chromosome::getLength () const { return length; }
int Chromosome::getLength () { return length; }
int Chromosome::getLengthLong () const { return lengthLong; }
int Chromosome::getLengthLong () { return lengthLong; }
unsigned long * Chromosome::getGene() { return gene; }    
unsigned long * Chromosome::getGene() const { return gene; }

bool Chromosome::isEvaluated () const {
    return evaluated;
}

bool Chromosome::hasSeen() const {
    unordered_map<unsigned long, double>::iterator it = cache.find(key);
    if (it != cache.end())
        return true;
    return false;
}

double Chromosome::evaluate (int& counter) {
    if (CACHE)
        if (hasSeen()) {
            evaluated = true;
            return cache[key];
        }
    
    if (!Chromosome::hit) {
        ++nfe;
        ++counter;
    }
        
    evaluated = true;
    double accum = 0.0;

    switch (function) {
        case ONEMAX:
            accum = oneMax();
            break;
        case MKTRAP:
            accum = mkTrap(1, 0.8);
            break;
        case CYCTRAP:
            accum = cycTrap(1, 0.8);
            break;
        case FTRAP:
            accum = fTrap();
            break;
        case SPINGLASS:
            accum = spinGlass();
            break;
        case NK:
            accum = nkFitness();
            break;
        case SAT:
            accum = satFitness();
            break;
        case MAXCUT:
            accum = maxcutFitness();
            break;
        // NOTE [Add new problems] 2021-12-06 15:28:52
        case LEADINGONES:
            accum = leadingOnes();
            break;
        default:
            accum = mkTrap(1, 0.8);
            break;
    }

    if (CACHE)
        cache[key]=accum;

    return accum;
}


double Chromosome::maxcutFitness() const {
    char *x = new char[length];

    for (int i=0; i!=length; ++i){
        x[i] = getVal(i);
    }

    double result = evaluateMAXCUT(x, &myMAXCUT);
    delete []x;
    return result;
}

double Chromosome::spinGlass () const {

    int *x = new int[length];
    double result;

    for (int i=0; i<length; i++)
        if (getVal(i) == 1)
            x[i] = 1;
        else
            x[i] = -1;

    result = evaluateSPIN(x, &mySpinGlassParams);

    delete []x;

    return result;
}

double Chromosome::nkFitness() const {
    char *x = new char[length];

    for ( int i = 0; i < length; ++i) {
        x[i] = (char) getVal(i);
    }

    double result = evaluateNKProblem(x, &nkwa);
    //double result = evaluateNKWAProblem(x, &nkwa);
    delete []x;
    return result;
}

// OneMax
double Chromosome::oneMax () const {

    double result = 0;

    for (int i = 0; i < length; ++i)
        result += getVal(i);

    return result;
}

bool Chromosome::operator== (const Chromosome& c) const {
    
    // return getKey() == c.getKey();
    
    if (length != c.length)
        return false;

    for (int i=0; i<lengthLong; i++)
        if (gene[i] != c.gene[i])
            return false;

    return true;
}

Chromosome& Chromosome::operator= (const Chromosome& c) {

    if (length != c.length) {
        length = c.length;
        init (length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;
    lengthLong = c.lengthLong;
    key = c.key;

    memcpy(gene, c.gene, sizeof(long) * lengthLong);

    return *this;
}

/* -------------------------- SECTION TEST PROBLEMS ------------------------- */
double Chromosome::trap (int unitary, double fHigh, double fLow, int trapK) const {
    if (unitary > trapK)
        return 0;

    if (unitary == trapK)
        return fHigh;
    else
        return fLow - unitary * fLow / (trapK-1);
}
double Chromosome::fTrap() const {
    double result = 0.0;
    for (int i=0; i<length/6; ++i) {
        int u=0;
        for (int j=0; j<6; ++j)
            u += getVal(i*6+j);
        if (u==0)
            result += 1.0;
        else if (u==1)
            result += 0.0;
        else if (u==2)
            result += 0.4;
        else if (u==3)
            result += 0.8;
        else if (u==4)
            result += 0.4;
        else if (u==5)
            result += 0.0;
        else // u == 6
            result += 1.0;
    }
    return result;
}
double Chromosome::cycTrap(double fHigh, double fLow) const {
    int i, j;
    int u;
    int TRAP_M = length / (TRAP_K-1);
    if (length % (TRAP_K-1) != 0)
        outputErrMsg ("TRAP_k doesn't divide length for Cyclic Setting");
    double result = 0;
    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        int idx = i * TRAP_K - i;
        for (j = 0; j < TRAP_K; j++) {
            int pos = idx + j;
            if (pos == length)
                pos = 0;
            else if (pos > length)
                outputErrMsg ("CYCLIC BUG");
            //
            u += getVal(pos);
        }
        result += trap (u, fHigh, fLow, TRAP_K);
    }
    return result;
}
double Chromosome::mkTrap (double fHigh, double fLow) const {
    int i, j;
    int u;
    int TRAP_M = length / TRAP_K;
    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");
    double result = 0;
    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        result += trap (u, fHigh, fLow, TRAP_K);
    }
    return result;
}
// NOTE [Add new problem] 2021-12-06 15:28:27
double Chromosome::leadingOnes() const {
    double numOnes(0.0);
    for (size_t idx=0; idx!=size_t(length); ++idx) {
        if (getVal(idx))
            ++numOnes;
        else
            break;
    }
    return numOnes;
}
/* -------------------------------------------------------------------------- */

double Chromosome::getMaxFitness () const {
    double maxF;
    switch (function) {
        case ONEMAX:
            maxF = length;
            break;
        case MKTRAP:
            maxF = length/TRAP_K;
            break;
        case FTRAP:
            maxF = length/6;
            break;
        case CYCTRAP:
            maxF =  length/(TRAP_K - 1);
            break;
        case SPINGLASS:
            maxF = mySpinGlassParams.opt;
            break;
        case NK:
            maxF = nkwa.maxF;
            break;
        case SAT:
            maxF = 0;
            break;
        case MAXCUT:
            maxF = myMAXCUT.opt;
            break;
        // NOTE [Add new problem] 2021-12-06 15:28:40
        case LEADINGONES:
            maxF = length;
            break;
        default:
            // Never converge
            maxF = INF;
    }
    return maxF - EPSILON;
}

// contribute to lsnfe
bool Chromosome::tryFlipping(int index) {
    // int oldNFE = nfe;
    double oldF = getFitness(Chromosome::lsnfe);
    flip(index);
    if (getFitness(Chromosome::lsnfe) - EPSILON <= oldF) {
        flip(index);
        evaluated = true;
        fitness = oldF;
        // lsnfe += nfe - oldNFE;
        // nfe = oldNFE;
        return false;
    } else {
        // lsnfe += nfe - oldNFE;
        // nfe = oldNFE;
        return true;
    }
}

bool Chromosome::GHC() {

    int* order = new int [length];
    myRand.uniformArray(order, length, 0, length-1);

    bool flag = false;
    for (int i=0; i<length; ++i) {
        if (tryFlipping(order[i])) flag = true;
    }

    delete []order;
    return flag;

}

double Chromosome::satFitness() const {
    int *x = new int[length];

    for ( int i = 0; i < length; ++i) {
        x[i] = getVal(i);
    }

    double result = evaluateSAT(x, &mySAT);
    delete []x;
    return result;
}

ostream& operator<<(ostream& os, const Chromosome& ch){
    for (int i=0; i!=ch.getLength(); ++i) {
        if (ch.getVal(i)) {
            cout << GREEN << 1 << RESET;
        } else {
            cout << RED << 0 << RESET;
        }
    } 
    printf(" %6.3f", ch.getFitness());
    return os;
}

bool Chromosome::operator>(const Chromosome& c) const{
    if (getKey() == c.getKey()) {
        return false;
    }    
    if (getFitness() > c.getFitness()) {
        return true;
    } else if (getFitness() == c.getFitness() && getCount() >= c.getCount()) {
        return true;
    } else {
        return false;
    }
}

bool Chromosome::operator<(const Chromosome& c) const{
    if (getKey() == c.getKey()) {
        return false;
    }
    return *this > c;
    if (getFitness() < c.getFitness()) {
        return true;
    } else if (getFitness() == c.getFitness() && getCount() <= c.getCount()) {
        return true;
    } else {
        return false;
    }
}

unsigned int Chromosome::getCount(void) const {
    return DSMGA2::ch_count[getKey()];
}