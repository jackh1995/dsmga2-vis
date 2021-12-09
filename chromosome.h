/***************************************************************************
 *   Copyright (C) 2011 by TEIL                                        *
 *                                                                         *
 ***************************************************************************/
#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include <unordered_map>
#include "global.h"
#include "nk-wa.h"

using namespace std;

class Chromosome {

public:

    static enum Function {
        ONEMAX      = 0,
        MKTRAP      = 1,
        FTRAP       = 2,
        CYCTRAP     = 3,
        NK          = 4,
        SPINGLASS   = 5,
        SAT         = 6,
        MAXCUT      = 7,
        // NOTE [Add new problems]
        LEADINGONES = 8
    } function;


    Chromosome ();
    Chromosome (int n_ell);
    Chromosome(const Chromosome &);

    ~Chromosome ();

    bool hasSeen() const;

    bool GHC();
    void steepestDescent();

    void init (int _ell);
    void init0 (int _ell);
    void initR (int _ell);

    bool tryFlipping (int index);

    int getVal (int index) const {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        if ( (gene[q] & (1lu << r)) == 0 )
            return 0;
        else
            return 1;
    }

    void setVal (int index, int val) {

        assert (index >= 0 && index < length);

        if (getVal(index) == val) return;

        setValF(index, val);
        key ^= zKey[index];
    }

    unsigned long getKey () const {
        return key;
    }


    void setValF (int index, int val) {

        assert (index >= 0 && index < length);
        //outputErrMsg ("Index overrange in Chromosome::operator[]");

        int q = quotientLong(index);
        int r = remainderLong(index);

        if (val == 1)
            gene[q] |= (1lu<<r);
        else
            gene[q] &= ~(1lu<<r);

        evaluated = false;
    }

    void flip (int index) {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        gene[q] ^= (1lu<<r);
        key ^= zKey[index];

        evaluated = false;
    }

    /** real evaluator */
    double evaluate (int& counter);

    bool isEvaluated () const;

    bool operator== (const Chromosome & c) const;
    Chromosome& operator= (const Chromosome & c);

    double getFitness (int& counter);
    double getFitness () const;
    int getLength ();
    int getLength () const;
    int getLengthLong ();
    int getLengthLong () const;
    unsigned long * getGene();
    unsigned long * getGene() const;
    double trap (int u, double high, double low, int trapK) const;
    double oneMax () const;
    double mkTrap (double high, double low) const;
    double cycTrap(double fHigh, double fLow) const;
    double fTrap () const;
    // NOTE [Add new problem] 2021-12-06 15:23:14
    double leadingOnes() const;
    double spinGlass () const;
    double nkFitness() const;
    double satFitness() const;
    double maxcutFitness() const;
    void setLength ();
    double getMaxFitness () const;

public:
    static int nfe;
    static int lsnfe;
    static int initnfe;
    static int rmnfe;
    static int bmnfe;
    static int bmeq;
    static int bmgt;
    static int bmfail;
    static bool hit;
    static unordered_map<unsigned long, double> cache;
    friend ostream& operator<<(ostream& os, const Chromosome& ch);
    bool operator>(const Chromosome&) const;
    bool operator<(const Chromosome&) const;
    unsigned int getCount(void) const;

private:
    unsigned long *gene;
    int length;
    int lengthLong;
    double fitness;
    bool evaluated;
    unsigned long key;
};

#endif
