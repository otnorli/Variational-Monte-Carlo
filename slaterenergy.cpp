#include "slaterenergy.h"

SlaterEnergy::SlaterEnergy(const mat &atom_R, double Radi)
{
    h = 0.001;
    h2 = 1000000;
    R = atom_R;
    R(1,0) = Radi;
}

double SlaterEnergy::local_energy(const mat &r, double alpha, double wfold, int dim, int number_particles, int charge, double beta, int number_atoms)
{
    Wavefunction WaveFunksjon(dim, number_particles, number_atoms, R(1,0));

    int atom_nr;
    e_kinetic = 0;

    //Kinetisk energi beregning, H*psi
    //e_kinetic = -0.5 * WaveFunksjon.Laplace_Ratio(r, number_particles, alpha);
    //number_atoms = 1;
    for (atom_nr = 0; atom_nr < number_atoms; atom_nr++)
    {
        e_kinetic -= 0.5 * WaveFunksjon.Kinetic_Energy_Combo(r, beta, number_particles, alpha, atom_nr);
    }

    //Beregner potensiell energi
    e_potential = 0;

    //Bidrag fra elektron-proton virkning
    for (int Y = 0; Y < number_atoms; Y++)
    {
        for (i=0;i<(number_atoms*number_particles);i++){
            r_single_particle = 0;
            for (j=0;j<dim;j++){
                r_single_particle += (r(i,j)-R(Y,j))*(r(i,j)-R(Y,j));
            }
            //e_potential -= charge/sqrt(r_single_particle);
            e_potential -= charge/sqrt(r_single_particle);
        }
    }

    //Bidrag fra elektron-elektron virkning
    for (i=0;i<(number_atoms*number_particles)-1;i++){
        for (j=i+1;j<(number_atoms*number_particles);j++){
            r_12 = 0;
            for (k=0;k<dim;k++){
                r_12 += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }
            e_potential += 1/sqrt(r_12);
        }
    }

    //Bidrag fra atom-atom virkning
    r_12 = 0;
    for (i=0; i<number_atoms-1; i++)
    {
        for (j=i+1; j<number_atoms; j++)
        {
            for (k=0; k<dim;k++)
            {
                r_12 += (R(i,k)-R(j,k)) * (R(i,k)-R(j,k));
            }
        }
    }

    r_12 = sqrt(r_12);
    //e_potential += charge/r_12;
    e_potential += charge*charge/r_12;

    e_local = e_potential + e_kinetic;
    return e_local;
}
