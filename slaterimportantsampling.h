#ifndef SLATERIMPORTANTSAMPLING_H
#define SLATERIMPORTANTSAMPLING_H

#include <iostream>
#include "slaterenergy.h"
#include "wavefunction.h"
#include <string>
#include "../../home/ole/Desktop/include/armadillo"
#include <fstream>
#include <cmath>

using namespace std;
using namespace arma;
class SlaterImportantSampling
{
public:
    SlaterImportantSampling(int d, int npar, int cha, int maxvar, int numcyc, double stepl, double hhh, double hhh2, double mora, int num_at, const mat &atom_R);
    void ImpSamp(double *cumulative_e, double *cumulative_e2, double beta, double alpha, int my_rank, bool isminimizing, int totalranks);
    double NormalDist(long int *idum);
    mat SetTempR(const mat &r, int at_nr);

private:
    int number_particles, dimension, charge;
    int cycles, variate, variate2, accept, i, j, k;
    int J, I;
    double alpha, betaperm, Jastrow_old, Jastrow_new, wfnew, wfold;
    double energy, energy2, delta_e;
    //double wfnew, wfold, wfold_up, wfold_down, wfnew_up, wfnew_down;
    vec wfold_up, wfold_down, wfnew_up, wfnew_down;

    long int idum;
    string name;
    double greensfunction;
    double timestep;
    int max_variations, number_cycles;
    double step_lenght;
    double D;
    double h, h2;
    double temporary, temporary2;
    double thermalization;
    int number_atoms, atom_nr;
    mat R;
    double RadPos;
};

#endif // SLATERIMPORTANTSAMPLING_H
