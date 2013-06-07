#ifndef MINIMAZATION_H
#define MINIMAZATION_H

#include <iostream>
#include <string>
#include "slaterimportantsampling.h"
#include "../../home/ole/Desktop/include/armadillo"

using namespace std;
using namespace arma;

class minimazation
{
public:
    minimazation(double alfha, double beta, int number_variables, int numpar, int charg, double steplen, double hhhh, double hhhh2, int dimdim, int myran, int totsiz);
    vec ConjugateGratient(double inputR);
    vec returnEnergy(double beta, double alpha, int xaxa, const mat &Rmat);

private:
    int max_variations, cycles, thermalization;
    double alpha, beta;
    int dimension, number_particles, charge;
    double h, h2, step_lenght;
    int xyzdimension;
    int my_rank;
    double timestep;
    int size;
    double Rpos;
    mat R;
    mat R_temp;
};

#endif // MINIMAZATION_H
