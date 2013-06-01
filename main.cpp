#include <iostream>
#include <cmath>
#include <fstream>
#include "../../home/ole/Desktop/include/armadillo"
#include "slaterimportantsampling.h"
#include "minimazation.h"

using namespace std;
using namespace arma;

int main()
{
    //E0 = 2.9, Helium
    //E0 = 14.7, Beryllium

    //Initialiserer
    double h = 0.001;
    double h2 = 1000000;
    double beta;
    int i, number_cycles, max_variations, thermalization, charge, dimension, number_particles, number_atoms;
    double step_lenght, alpha;
    int number_variableZ;
    double *cumulative_e, *cumulative_e2;
    double *Renergy;
    mat R;
    //vec variable_list(2);
    //ofstream out;
    //out.open("data.txt");

    //Spesifiserer startdata
    dimension = 3;
    number_particles = 4;
    charge = number_particles; //Hvis du ikke har ioner
    max_variations = 1;
    number_cycles  = 200000;
    thermalization =  10000;
    step_lenght = 0.05;
    //beta = 0.94714228;
    //alpha = 4.62;

    beta = 0.87731;
    alpha = 3.92879;

    number_variableZ = 3;
    number_atoms = 2;

    R = zeros(number_atoms, dimension);
    R(1,0) = 4.5625;

    //vec variable_list(number_variableZ);
    //variable_list(0) = alpha;
    //variable_list(1) = beta;
    //variable_list(2) = R(1,0);

    cumulative_e = new double[max_variations*max_variations];
    cumulative_e2 = new double[max_variations*max_variations];
    Renergy = new double[100];

    //minimerer
    minimazation Mini(alpha, beta, number_variableZ, number_particles, charge, step_lenght, h, h2, dimension, 0, 0);

    //variable_list = Mini.ConjugateGratient();
    //alpha = variable_list(0);
    //beta = variable_list(1);
    //R(1,0) = variable_list(2);

    for (int T = 1; T < 101; T++)
    {
        R(1,0) = 3.53 + (double) T / 100;

        //Importantsampling med Slater determinant som bølgefunksjon
        SlaterImportantSampling SIS(dimension, number_particles, charge, max_variations, number_cycles, step_lenght, h, h2, thermalization, number_atoms, R);

        SIS.ImpSamp(cumulative_e, cumulative_e2, beta, alpha, 0, false, 1);
        Renergy[T] = cumulative_e[0];
    }

    //Printer litt resultater
    //out << "Alfa var: " << alpha << endl;
    //out << "Beta var: " << beta << endl;
    //out << "R var: " << R(1,0) << endl;
    //out << "Energien: " << endl;
    for (i=0;i<(max_variations*max_variations);i++)
    {
    //    out << cumulative_e[i] << " " << cumulative_e2[i] << endl;
    }

    cout << endl;
    for (int T = 1; T < 101; T++)
    {
        cout << Renergy[T] << " ";
    }

    //Ferdig med programmet
    //out.close();
    cout << "Program complete." << endl;
    return 0;
}
