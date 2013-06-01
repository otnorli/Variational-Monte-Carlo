#ifndef SLATERENERGY_H
#define SLATERENERGY_H

#include <iostream>
#include "wavefunction.h"
#include <string>
#include "../../home/ole/Desktop/include/armadillo"

using namespace std;
using namespace arma;

class SlaterEnergy
{
public:
    SlaterEnergy(const mat &atom_R, double Radi);
    double local_energy(const mat &r, double alpha, double wfold, int dim, int number_particles, int charge, double beta, int number_atoms);

private:
    int i,j,k;
    double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12, r_single_particle;
    double h;
    double h2;
    mat R;
};

#endif // SLATERENERGY_H
