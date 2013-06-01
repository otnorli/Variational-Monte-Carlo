#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <string>
#include <iostream>
#include <cmath>
#include "../../home/ole/Desktop/include/armadillo"
#include "../../home/ole/Desktop/cppLibrary/lib.h"

using namespace std;
using namespace arma;

class Wavefunction
{
public:
    Wavefunction(int dimens, int number_particles, int num_at, double Radi);
    double psi_orbital_slater(const mat &r, const int &particle_number, const int &Orbital_number, double alpha, int atom_nr);
    double slater_det(const mat &r, int number_particle, double alpha, double beta);
    mat slater_reduced_det(const mat &r, int K, double alpha, const int &number_particles, int atom_nr);
    double Jastrow_Part(const mat &r, int number_particles, double beta);
    double Laplacian(const mat &r, const int number_particle, const int &Orbital_number, double alpha, int atom_nr);
    rowvec3 Gradient(const mat &r, const int number_particle, const int &Orbital_number, double alpha, int atom_nr);
    double Laplace_Ratio(const mat &r, const int number_particle, double alpha, int atom_nr);
    mat Gradient_Radio(const mat &r, double alpha, int number_particles, int at_nr);
    mat Jastrow_Gradient_Ratio(const mat &r, const int &number_particles, const double &beta);
    double Jastrow_Laplace_Ratio(const mat &r, const double &beta, const int &number_particles);
    double Kinetic_Energy_Combo(const mat &r, const double &beta, const int &number_particles, const double &alpha, int atom_nr);
    double seta(int atom1, int atom2, int number_particles);

private:
    int dimension;
    mat Spin_Up_Slater_det;
    mat Spin_Down_Slater_det;
    mat Spin_Up_Slater_det_Inv;
    mat Spin_Down_Slater_det_Inv;
    mat R;
    int number_atoms;

};

#endif // WAVEFUNCTION_H
