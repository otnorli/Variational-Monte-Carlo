#include "wavefunction.h"

//Denne koden er SPESIFISERT til (1 atom) eller (2 atomer i et molekyl)

Wavefunction::Wavefunction(int dimens, int number_particles, int num_at, double Radi)
{
    dimension = dimens;
    int number_particless = number_particles/2;
    Spin_Up_Slater_det = zeros(number_particless ,number_particless);
    Spin_Down_Slater_det = zeros(number_particless ,number_particless);
    Spin_Up_Slater_det_Inv = zeros(number_particless, number_particless);
    Spin_Down_Slater_det_Inv = zeros(number_particless ,number_particless);
    number_atoms = num_at;
    R = zeros(number_atoms, dimension);

    //Insert startposisjoner her, for enkelthetsskyld
    R(0,0) = Radi;
}

//Generell orbitalfunksjon, returnerer verdien for gitt r og gitt orbitalnummer
double Wavefunction::psi_orbital_slater(const mat &r, const int &particle_number, const int &Orbital_number, double alpha, int atom_nr)
{
    double Sum_In_Exponential=0;
    double wavefunc_value=0;
    int k;

    for (k=0; k<dimension; k++){
        Sum_In_Exponential += (r(particle_number, k) - R(atom_nr, k))*(r(particle_number, k) - R(atom_nr, k));
    }

    Sum_In_Exponential = sqrt(Sum_In_Exponential);
    Sum_In_Exponential *= alpha;

    //For første orbital er dette prøvebølgefunksjonen
    if (Orbital_number == 0)
    {
        wavefunc_value = exp(-Sum_In_Exponential);
    }

    //For andre orbital er dette prøvebølgefunksjonen
    else if (Orbital_number == 1)
    {
        Sum_In_Exponential /= 2;
        wavefunc_value = (1 - Sum_In_Exponential) * exp(-Sum_In_Exponential);

        //wavefunc_value = (Sum_In_Exponential-2)*exp(-0.5*Sum_In_Exponential);
    }

    //For tredje orbital er dette prøvebølgefunksjonen
    else if (Orbital_number == 2)
    {
        wavefunc_value = (r(particle_number, 0)-R(atom_nr,0)) * exp(-Sum_In_Exponential/2);
    }

    //For fjerde orbital er dette prøvebølgefunksjonen
    else if (Orbital_number == 3)
    {
        wavefunc_value = (r(particle_number, 1)-R(atom_nr,1)) * exp(-Sum_In_Exponential/2);
    }

    //For femte orbital er dette prøvebølgefunksjonen
    else if (Orbital_number == 4)
    {
        wavefunc_value = (r(particle_number, 2)-R(atom_nr, 2)) * exp(-Sum_In_Exponential/2);
    }

    return wavefunc_value;
}

//Beregning av redusert slater determinant, etter man deler opp etter spin
mat Wavefunction::slater_reduced_det(const mat &r, int K, double alpha, const int &number_particles, int atom_nr)
{
    //Input er 2 atomer og finner determinanten
    mat Reduced_Slater_Matrix = zeros<mat>(number_particles/2,number_particles/2);

    K = (int) K*number_particles/2;

    int half_p = (int) number_particles/2;

    for (int i=0; i<half_p; i++)
    {
        for (int j=K; j<half_p+K; j++)
        {
            Reduced_Slater_Matrix(j-K, i) = psi_orbital_slater(r, j, i, alpha, atom_nr);
        }
    }
    return Reduced_Slater_Matrix;
}

//Beregning av jastrow faktoren
double Wavefunction::Jastrow_Part(const mat &r, int number_particles, double beta)
{
    double r_12;
    double Spin_variable;
    double Jastrow_Part_Argument = 0;
    double Jastrow_P=0;
    int i,j,k;

    for (i=0;i<number_particles-1;i++)
    {
        for (j=i+1;j<number_particles;j++)
        {
            r_12 = 0;
            for (k=0;k<dimension;k++){
                r_12 += (r(i,k)-r(j,k))*(r(i,k)-r(j,k));
            }

            r_12 = sqrt(r_12); //Avstanden mellom 2 atomer

            //Variabelen avhenger av spin opp eller spin ned
            //Vi velger spin opp til å være de første halvparten av atomene, resten ned
/*
            if (((i < number_particles/2) && (j < number_particles/2)) || ((i > number_particles/2) && (j > number_particles/2)))
            {
                Spin_variable = 1.0/4;
            }

            else
            {
                Spin_variable = 1.0/2;
            }
*/
            Spin_variable = seta(i,j, number_particles);

            //Summerer opp argumentene
            Jastrow_Part_Argument += Spin_variable*r_12/(1+beta*r_12);
            //Jastrow_Part_Argument /= number_atoms; // <-- number atoms variabelen brukes her av TEKNISKE grunner
        }
    }
    Jastrow_P = exp(Jastrow_Part_Argument);
    return Jastrow_P;
}

//Beregning av dobbeltderiverte til orbitalene
double Wavefunction::Laplacian(const mat &r, const int number_particle, const int &Orbital_number, double alpha, int atom_nr)
{
    double Sum_I_Eksponenten = 0;
    double Laplacen = 0;

    for (int k = 0; k < dimension; k++)
    {
        Sum_I_Eksponenten += (r(number_particle, k) - R(atom_nr, k)) * (r(number_particle, k) - R(atom_nr,k));
    }

    Sum_I_Eksponenten = sqrt(Sum_I_Eksponenten);

    //Velger orbital

    if (Orbital_number == 0)
    {
        Laplacen = alpha * (alpha * Sum_I_Eksponenten - 2) * exp(-alpha*Sum_I_Eksponenten) / Sum_I_Eksponenten;
    }

    else if (Orbital_number == 1)
    {
        //Laplacen = alpha*(alpha*Sum_I_Eksponenten-8)*(alpha*Sum_I_Eksponenten-2)*exp(-0.5*alpha*Sum_I_Eksponenten)/(4*Sum_I_Eksponenten);
        Laplacen = -alpha*(alpha*Sum_I_Eksponenten-8)*(alpha*Sum_I_Eksponenten-2)*exp(-0.5*alpha*Sum_I_Eksponenten)/(8*Sum_I_Eksponenten);
    }

    else if (Orbital_number == 2)
    {
        Laplacen = alpha * (r(number_particle, 0)-R(atom_nr, 0)) * exp(-alpha*Sum_I_Eksponenten/2) * (alpha*Sum_I_Eksponenten - 8) / (4 * Sum_I_Eksponenten);
    }

    else if (Orbital_number == 3)
    {
        Laplacen = alpha * (r(number_particle, 1)-R(atom_nr,1)) * exp(-alpha*Sum_I_Eksponenten/2) * (alpha*Sum_I_Eksponenten - 8) / (4 * Sum_I_Eksponenten);
    }

    else if (Orbital_number == 4)
    {
        Laplacen = alpha * (r(number_particle, 2)-R(atom_nr,2)) * exp(-alpha*Sum_I_Eksponenten/2) * (alpha*Sum_I_Eksponenten - 8) / (4 * Sum_I_Eksponenten);
    }

    return Laplacen;
}

//Beregning av enkeltderiverte til orbitalene
rowvec3 Wavefunction::Gradient(const mat &r, const int number_particle, const int &Orbital_number, double alpha, int atom_nr)
{
    double Sum_I_Eksponenten = 0;
    rowvec3 Gradienten;
    Gradienten(0)=0; Gradienten(1)=0; Gradienten(2)=0;

    for (int k = 0; k < dimension; k++)
    {
        Sum_I_Eksponenten += (r(number_particle, k) - R(atom_nr, k)) * (r(number_particle, k) - R(atom_nr, k));
    }

    Sum_I_Eksponenten = sqrt(Sum_I_Eksponenten);

    if (Orbital_number == 0)
    {
        Gradienten = -(alpha * (r.row(number_particle)-R.row(atom_nr)) / Sum_I_Eksponenten) * exp(-alpha*Sum_I_Eksponenten);
    }

    if (Orbital_number == 1)
    {
        Gradienten = (alpha*(r.row(number_particle) - R.row(atom_nr))*(alpha*Sum_I_Eksponenten-4) / (4*Sum_I_Eksponenten)) * exp(-alpha*Sum_I_Eksponenten/2);
    }

    if (Orbital_number == 2)
    {
        Gradienten(0) = -alpha*(r(number_particle, 0)-R(atom_nr, 0)) * (r(number_particle ,0)-R(atom_nr, 0)) + 2*Sum_I_Eksponenten;
        Gradienten(1) = -alpha*(r(number_particle, 0)-R(atom_nr,0)) * (r(number_particle ,1)-R(atom_nr,1));
        Gradienten(2) = -alpha*(r(number_particle, 0)-R(atom_nr,0)) * (r(number_particle ,2)-R(atom_nr,2));
        Gradienten = Gradienten * exp(-alpha*Sum_I_Eksponenten / 2) / (2*Sum_I_Eksponenten);
    }

    if (Orbital_number == 3)
    {
        Gradienten(0) = -alpha*(r(number_particle, 1)-R(atom_nr,1)) * (r(number_particle ,0)-R(atom_nr, 0));
        Gradienten(1) = -alpha*(r(number_particle, 1)-R(atom_nr,1)) * (r(number_particle ,1)-R(atom_nr,1)) + 2*Sum_I_Eksponenten;
        Gradienten(2) = -alpha*(r(number_particle, 1)-R(atom_nr,1)) * (r(number_particle ,2)-R(atom_nr,2));
        Gradienten = Gradienten * exp(-alpha*Sum_I_Eksponenten / 2) / (2*Sum_I_Eksponenten);
    }

    if (Orbital_number == 4)
    {
        Gradienten(0) = -alpha*(r(number_particle, 2)-R(atom_nr,2)) * (r(number_particle ,0)-R(atom_nr, 0));
        Gradienten(1) = -alpha*(r(number_particle, 2)-R(atom_nr,2)) * (r(number_particle ,1)-R(atom_nr,1));
        Gradienten(2) = -alpha*(r(number_particle, 1)-R(atom_nr,1)) * (r(number_particle ,2)-R(atom_nr,2)) + 2*Sum_I_Eksponenten;
        Gradienten = Gradienten * exp(-alpha*Sum_I_Eksponenten / 2) / (2*Sum_I_Eksponenten);
    }

    return Gradienten;
}


double Wavefunction::Laplace_Ratio(const mat &r, const int number_particles, double alpha, int atom_nr)
{
    double Ratio = 0;

    //Slater up biten
    Spin_Up_Slater_det = slater_reduced_det(r, 0, alpha, number_particles, atom_nr);
    Spin_Up_Slater_det_Inv = inv(Spin_Up_Slater_det);

    //Slater ned biten
    Spin_Down_Slater_det = slater_reduced_det(r, 1, alpha, number_particles, atom_nr);
    Spin_Down_Slater_det_Inv = inv(Spin_Down_Slater_det);

    for (int i = 0; i < number_particles/2; i++){
        for (int j = 0; j < number_particles/2; j++){
            Ratio += Laplacian(r, i, j, alpha, atom_nr) * Spin_Up_Slater_det_Inv(j,i)
                    + Laplacian(r, (i + number_particles/2), j, alpha, atom_nr) * Spin_Down_Slater_det_Inv(j,i);
        }
    }
    return Ratio;
}


mat Wavefunction::Gradient_Radio(const mat &r, double alpha, int number_particles, int at_nr)
{
    //Slater up biten
    Spin_Up_Slater_det = slater_reduced_det(r, 0, alpha, number_particles, at_nr);
    Spin_Up_Slater_det_Inv = inv(Spin_Up_Slater_det);

    //Slater ned biten
    Spin_Down_Slater_det = slater_reduced_det(r, 1, alpha, number_particles, at_nr);
    Spin_Down_Slater_det_Inv = inv(Spin_Down_Slater_det);

    mat gradientradio = zeros(number_particles, dimension);

    for (int i = 0; i < number_particles/2; i++)
    {
        for (int j=0; j<number_particles/2; j++)
        {
            gradientradio.row(i) += Gradient(r, i, j, alpha, at_nr) * Spin_Up_Slater_det_Inv(j,i);
            gradientradio.row(i+number_particles/2) += Gradient(r, i + number_particles/2, j, alpha, at_nr) * Spin_Down_Slater_det_Inv(j,i);
        }
    }

    return gradientradio;
}


mat Wavefunction::Jastrow_Gradient_Ratio(const mat &r, const int &number_particles, const double &beta)
{
    mat Jastrow_Gradient_Radio = zeros(number_particles, dimension);
    rowvec r_12_vec(3);
    rowvec temp_r(3);
    double r_12, arg, Spin_variable;
    double bb;

    for (uint i = 0; i < number_particles; i++)
    {
        for (uint j = 0; j < i; j++)
        {

                temp_r(0)=0;
                temp_r(1) =0;
                temp_r(2)=0;

                r_12_vec = r.row(i) - r.row(j);

                r_12 = 0;
                for (uint k = 0; k < dimension; k++)
                {
                    r_12 += r_12_vec(k) * r_12_vec(k);
                }
                r_12 = sqrt(r_12);

                //Finner spin variablen

                /*
                if (((i < number_particles/2) && (j < number_particles/2)) || ((i > number_particles/2) && (j > number_particles/2)))
                {
                    Spin_variable = 1.0/4;
                }

                else
                {
                    Spin_variable = 1.0/2;
                }
*/
                Spin_variable = seta(i,j,number_particles);

                bb = 1+beta*r_12;
                //bb = number_atoms * bb; //number_atoms brukes av TEKNISKE grunner
                arg = Spin_variable / (bb*bb);
                Jastrow_Gradient_Radio.row(i) += arg * r_12_vec / r_12;
        }

        for (uint j = i+1; j < number_particles; j++)
        {
            r_12_vec = r.row(i) - r.row(j);

            r_12 = 0;
            for (uint k = 0; k < dimension; k++)
            {
                r_12 += r_12_vec(k) * r_12_vec(k);
            }
            r_12 = sqrt(r_12);

            //Finner spin variablen

            /*
            if (((i < number_particles/2) && (j < number_particles/2)) || ((i > number_particles/2) && (j > number_particles/2)))
            {
                Spin_variable = 1.0/4;
            }
            else
            {
                Spin_variable = 1.0/2;
            }
            */

            Spin_variable = seta(i,j,number_particles);

                bb = 1+beta*r_12;
                //bb = number_atoms * bb; //number_atoms brukes av TEKNISKE grunner
                arg = Spin_variable / (bb*bb);
                Jastrow_Gradient_Radio.row(i) += r_12_vec*arg/r_12;
            }
    }
    return Jastrow_Gradient_Radio;
}


double Wavefunction::Jastrow_Laplace_Ratio(const mat &r, const double &beta, const int &number_particles)
{
    mat Jastrow_Gradient_Radio = zeros(number_particles, dimension);
    Jastrow_Gradient_Radio = Jastrow_Gradient_Ratio(r, number_particles, beta);
    double Laplace_Radio=0, Grad_Radio;
    rowvec r_12_vec(3);
    double r_12, Spin_variable;
    double bb;

    //Tar vekselvirkningene for alle untatt i med seg selv
    for (uint i = 0; i < number_particles; i++)
    {
        for (uint j = 0; j < i; j++)
        {
            r_12 = 0;
            r_12_vec = r.row(i) - r.row(j);
            for (uint k=0; k < dimension; k++)
            {
                r_12 += r_12_vec(k) * r_12_vec(k);
            }
            r_12 = sqrt(r_12);

            /*
            if (((i < number_particles/2) && (j < number_particles/2)) || ((i > number_particles/2) && (j > number_particles/2)))
            {
                Spin_variable = 1.0/4;
            }

            else
            {
                Spin_variable = 1.0/2;
            }
            */

            Spin_variable = seta(i,j,number_particles);

            bb = (1+beta*r_12);
            //bb *= number_atoms; //number_atoms brukes av TEKNISKE grunner

            Laplace_Radio += 2*Spin_variable/(r_12*bb*bb);
            Laplace_Radio -= 2*beta*Spin_variable/(bb*bb*bb);
        }
        for (uint j = i+1; j < number_particles; j++)
        {
            r_12 = 0;
            r_12_vec = r.row(i) - r.row(j);
            for (uint k=0; k < dimension; k++)
            {
                r_12 += r_12_vec(k) * r_12_vec(k);
            }
            r_12 = sqrt(r_12);

            /*
            if (((i < number_particles/2) && (j < number_particles/2)) || ((i > number_particles/2) && (j > number_particles/2)))
            {
                Spin_variable = 1.0/4;
            }

            else
            {
                Spin_variable = 1.0/2;
            }
            */

            Spin_variable = seta(i,j,number_particles);

            bb = (1+beta*r_12);
            //bb *= number_atoms; //number_atoms brukes av TEKNISKE grunner

            Laplace_Radio += 2*Spin_variable/(r_12*bb*bb);
            Laplace_Radio -= 2*beta*Spin_variable/(bb*bb*bb);
        }

        //Tar med faktoren fra gradienten
        Grad_Radio = 0;

        for (uint k = 0; k < dimension; k++)
        {
            Grad_Radio += Jastrow_Gradient_Radio(i, k) * Jastrow_Gradient_Radio(i, k);
        }

        Laplace_Radio += Grad_Radio;
    }

    return Laplace_Radio;
}


double Wavefunction::Kinetic_Energy_Combo(const mat &r, const double &beta, const int &number_particles, const double &alpha, int atom_nr)
{
    double EK = 0;
    mat Slat = zeros<mat>(number_particles, dimension);
    mat Jast = zeros<mat>(number_particles*number_atoms, dimension);

    mat temp_r = zeros<mat>(number_particles, dimension);

    for (int F = 0; F < number_particles; F++)
    {
        temp_r.row(F) = r.row(F+atom_nr*number_particles);
    }

    Slat = Gradient_Radio(temp_r, alpha, number_particles, atom_nr);
    Jast = Jastrow_Gradient_Ratio(r, (number_atoms*number_particles), beta);

    //Tar her halve laplacen
    for (uint i = 0; i < number_particles; i++)
    {
        EK += Slat(i,0)*Jast(i+atom_nr*number_particles,0)
                + Slat(i,1)*Jast(i+atom_nr*number_particles,1)
                + Slat(i,2)*Jast(i+atom_nr*number_particles,2);
    }
    EK *= 2;

    EK += Laplace_Ratio(temp_r, number_particles, alpha, atom_nr);

    //Legger til laplacen bare en gang
    if (atom_nr == 0)
    {
        EK += Jastrow_Laplace_Ratio(r, beta, number_atoms*number_particles);
    }

    return EK;
}

double Wavefunction::seta(int atom1, int atom2, int number_particles)
{
    double Spin_variable;

    if (atom1 >= number_particles/number_atoms)
    {
        atom1 = (int) atom1 - number_particles/number_atoms;
    }

    if (atom2 >= number_particles/number_atoms)
    {
        atom2 = (int) atom2 - number_particles/number_atoms;
    }

    if (((atom1 < number_particles/(2*number_atoms)) && (atom2 < number_particles/(2*number_atoms))) || ((atom1 >= number_particles/(2*number_atoms)) && (atom2 >= number_particles/(2*number_atoms))))
    {
        Spin_variable = 1.0/4;
    }

    else
    {
        Spin_variable = 1.0/2;
    }

    return Spin_variable;
}

