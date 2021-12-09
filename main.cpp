/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/
#include <math.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#include "chromosome.h"
using namespace std;

int main (int argc, char *argv[]) {
    // Default
    Problem problem;
    if (argc != 9 && argc != 4 && argc != 5) {
        printf ("DSMGA2 [problem] [ell] [pop] [maxGen] [maxFe] [repeat] [verbose] [randseed]\n");
        printf ("Functions: \n");
        printf ("    0. ONEMAX     : Any ell\n");
        printf ("    1. MK         : Any ell\n");
        printf ("    2. FTRAP      : Any ell\n");
        printf ("    3. CYC        : Any ell\n");
        printf ("    4. NK         : 25, 50, 100, 200, 400\n");
        printf ("    5. SPIN       : 36, 100, 196, 400, 784\n");
        printf ("    6. SAT        : 20, 50, 75, 100, 200\n");
        printf ("    7. MAXCUT     : 50, 75, 100\n");
        printf ("    8. LEADINGONES: Any ell\n");
        return -1;
    }
    auto it = Problem_table.find(argv[1]);
    if (it != Problem_table.end()) {
        problem = it->second;
    } else { throw runtime_error("Error: Problem not defined"); }
    int ell = atoi (argv[2]);
    int nInitial = atoi (argv[3]);
    int maxGen = 200;
    int maxFe = -1;
    int repeat = 1;
    verbose = NO;
    if (argc == 5) {
        verbose = Verbosity(atoi(argv[4]));
    }
    int rand_seed = -1;
    // Load the problem
    switch(problem) {
        char filename[200];
        char opt_filename[200];
        FILE *fp;
        int instance_id;
        // NOTE [Add new problem] 2021-12-06 15:37:33
        case onemax: case mktrap: case cyctrap: case ftrap: case leadingones:
            break;
        case nk: 
            instance_id = 1;
            sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, 5, instance_id);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            fp = fopen(filename, "r");
            loadNKWAProblem(fp, &nkwa);
            break;
        case spin:
            instance_id = 1;
            sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, instance_id);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            loadSPIN(filename, &mySpinGlassParams); 
            break;
        case maxsat:
            instance_id = 1;
            sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf", ell, ell, instance_id);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            loadSAT(filename, &mySAT);
            break;
        case maxcut:
            instance_id = 1;
            sprintf(filename, "./maxcut/w05_%d/w05_%d.%d", ell, ell, instance_id);
            sprintf(opt_filename, "./maxcut/g_w05_%d/g_w05_%d.%d", ell, ell, instance_id);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            if (SHOW_BISECTION) printf("Loading groundtruth: %s\n", opt_filename);
            loadMAXCUT(filename, opt_filename, &myMAXCUT);
            break;
        default:
            throw runtime_error("Error: Problem not defined");
    }
    // Customized
    if (argc != 4 && argc != 5) {
        maxGen = atoi (argv[4]);
        maxFe = atoi (argv[5]);
        repeat = atoi (argv[6]);
        verbose = Verbosity(atoi (argv[7]));
        rand_seed = atoi (argv[8]);
    }
    // Set random seed
    if (rand_seed != -1)
        myRand.seed((unsigned long)rand_seed);
    // Start evolving
    Statistics stGen, stFE, stINITFE, stLSFE, stRMFE, stBMFE;
    int usedGen;
    int failNum = 0;
    for (int i = 0; i < repeat; ++i) {
        DSMGA2 ga (ell, nInitial, maxGen, maxFe, problem);
        usedGen = ga.doIt();
        if (!ga.foundOptima()) {
            failNum++;
            if (repeat > 1)
                printf ("-");
        } else {
            stGen.record (usedGen);
            stFE.record (Chromosome::nfe);
            stINITFE.record(Chromosome::initnfe);
            stLSFE.record (Chromosome::lsnfe);
            stRMFE.record (Chromosome::rmnfe);
            stBMFE.record (Chromosome::bmnfe);
            if (repeat > 1) {
                printf ("+");
            }
        }
        // Print the best seen chromosome if repeat == 1
        if (repeat == 1) {
            cout << endl;
            cout << "SUMMARY: ";
            if (ga.foundOptima())
                cout << GREEN << "SUCCESS" << RESET << endl;
            else
                cout << RED << "FAIL" << RESET << endl;
            cout << ga.getPopulation()[ga.getBestIndex()];
        }
        fflush (NULL);
    }
    
    cout << fixed << setprecision(2) << endl;
    cout << "GEN=[" << BOLDCYAN << stGen.getMean() << RESET << "] ";
    cout << "NFE=[" << BOLDCYAN << stFE.getMean() << RESET << "] ";
    cout << "INITNFE=[" << BOLDCYAN << stINITFE.getMean() << RESET << "] ";
    cout << "LSNFE=[" << BOLDCYAN << stLSFE.getMean() << RESET << "] ";
    cout << "RMNFE=[" << BOLDCYAN << stRMFE.getMean() << RESET << "] ";
    cout << "BMNFE=[" << BOLDCYAN << stBMFE.getMean() << RESET << "] ";
    cout << "SUCCESS=[" << BOLDCYAN << repeat-failNum << '/' << repeat << RESET << "]\n";
    
    if (problem == 4) freeNKWAProblem(&nkwa);
    return EXIT_SUCCESS;
}