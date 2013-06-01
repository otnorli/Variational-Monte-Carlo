#include "slaterimportantsampling.h"

SlaterImportantSampling::SlaterImportantSampling(int d, int npar, int cha, int maxvar, int numcyc, double stepl, double hhh, double hhh2, double mora, int num_at, const mat &atom_R)
{
    idum = -1;
    D = 0.5;
    dimension = d;
    number_particles = npar;
    charge = cha;
    max_variations = maxvar;
    number_cycles = numcyc;
    step_lenght = stepl;
    h = hhh;
    h2 = hhh2;
    thermalization = mora;
    timestep = step_lenght;
    number_atoms = num_at;
    R = atom_R;
    RadPos = R(1,0);
}

void SlaterImportantSampling::ImpSamp(double *cumulative_e, double *cumulative_e2, double beta, double alpha, int my_rank, bool isminimizing, int totalranks)
{
    ofstream out;
    out.open("data.txt");
    //ostringstream a;
    //a << "data" << my_rank << ".txt";
    //string s = a.str();
    //out.open(s.c_str());
    Wavefunction WaveFunc(dimension, number_particles, number_atoms, RadPos);


    ofstream oga;
    ostringstream b;
    b << "posisjon.txt";
    string ss = b.str();
    oga.open(ss.c_str());

    //============================================0
    //Midlertidig
    number_particles = 4;
    number_atoms = 2;
    //============================================0

    int N = number_particles*number_atoms;
    SlaterEnergy Energi(R, RadPos);

    //mat r_old = zeros<mat>(number_particles, dimension);
    //mat r_new = zeros<mat>(number_particles, dimension);
    //mat qforce_old = zeros<mat>(number_particles, dimension);
    //mat qforce_new = zeros<mat>(number_particles, dimension);

    mat r_old = zeros<mat>(N, dimension);
    mat r_new = zeros<mat>(N, dimension);
    mat qforce_old = zeros<mat>(N, dimension);
    mat qforce_new = zeros<mat>(N, dimension);
    mat r_temp = zeros<mat>(number_particles, dimension);
    mat q_temp = zeros<mat>(number_particles, dimension);

    wfold_up = vec(2);
    wfold_down = vec(2);
    wfnew_up = vec(2);
    wfnew_down = vec(2);

    int par_num_check;

    //Varierer alfa her
    for (variate2 = 0; variate2 < max_variations; variate2++){
        //Varierer beta her
    for (variate = 0; variate < max_variations; variate++){

        //Initialiserer verdier
        energy = 0;
        energy2 = 0;
        accept = 0;
        delta_e = 0;

        //Initial trial position, skal variere avhengig av hvilken prossessor vi er på

        //Kjører noen runder i random funksjonen
        for (i=0; i<my_rank;i++)
        {
            r_old(0,0) = NormalDist(&idum);
        }

        //Initialiserer posisjoner
        for (i=0; i<N;i++)
        {
            atom_nr = (int) i/number_particles;
            for (j=0; j<dimension; j++)
            {
                //cout << i << " " << j << endl;
                r_old(i,j) = R(atom_nr,j) + NormalDist(&idum) * sqrt(timestep);
            }
        }


        //Finner bølgefunksjonsverdien i starten
        wfold = 1;
        for (atom_nr = 0; atom_nr < number_atoms; atom_nr++)
        {
            r_temp = SetTempR(r_old, atom_nr);
            wfold_up(atom_nr) = det(WaveFunc.slater_reduced_det(r_temp, 0, alpha, number_particles, atom_nr));
            wfold_down(atom_nr) = det(WaveFunc.slater_reduced_det(r_temp, 1, alpha, number_particles, atom_nr));

            wfold *= wfold_up(atom_nr) * wfold_down(atom_nr);
        }
        Jastrow_old = WaveFunc.Jastrow_Part(r_old, N, beta);
        wfold *= Jastrow_old;

        //atom_nr = 0;
        //Finner bølgefunksjonen i starten
        //wfold_up = det(WaveFunc.slater_reduced_det(r_old, 0, alpha, number_particles, atom_nr));
        //wfold_down = det(WaveFunc.slater_reduced_det(r_old, 1, alpha, number_particles, atom_nr));
        //Jastrow_old = WaveFunc.Jastrow_Part(r_old, number_particles, beta);
        //wfold = wfold_up * wfold_down * Jastrow_old;

        //Finner kvantekraften i starten
        qforce_old = 2*WaveFunc.Jastrow_Gradient_Ratio(r_old, N, beta);
        //qforce_old = WaveFunc.Jastrow_Gradient_Ratio(r_old, N, beta);
        for (atom_nr = 0; atom_nr < number_atoms; atom_nr++)
        {
            r_temp = SetTempR(r_old, atom_nr);
            q_temp = WaveFunc.Gradient_Radio(r_temp, alpha, number_particles, atom_nr);

            for (k = 0; k < number_particles; k++)
            {
                qforce_old.row(k+atom_nr*number_particles) += 2*q_temp.row(k);
                //qforce_old.row(k+atom_nr*number_particles) += q_temp.row(k);
            }
        }
        //qforce_old = (WaveFunc.Jastrow_Gradient_Ratio(r_old, number_particles, beta) + WaveFunc.Gradient_Radio(r_old, alpha, number_particles, atom_nr));

        //Analytisk derivert den over her

        for (cycles=0; cycles <= (thermalization+number_cycles); cycles++){

            //cout << cycles << endl;

            //New position
            for (i=0; i<N; i++)
            {
                //Bestem atom nr.
                atom_nr = (int) i/number_particles;
                par_num_check = i;
                while (par_num_check >= number_particles)
                {
                    par_num_check -= number_particles;
                }

                //Oppdaterer posisjonen til partikkel i
                for (j=0; j<dimension; j++)
                {
                    r_new(i,j) = r_old(i,j) + NormalDist(&idum) * sqrt(timestep) + qforce_old(i,j) * timestep * D;
                }

                //Bare en partikkel av gangen, setter de andre posisjonene til r_old
                for (k=0; k<N;k++)
                {
                    if (k != i)
                    {
                        for (j=0; j<dimension; j++)
                        {
                            r_new(k,j) = r_old(k,j);
                        }
                    }
                }



                //Oppdaterer bølgefunksjonen med nye posisjonen.
                //Avgjør om partikkelen har spin opp eller ned
                //og endrer dermed bare en av redusert slater determinantene

                //if (i < number_particles/2)
                if (par_num_check < number_particles/2)
                {
                    //wfnew_up = det(WaveFunc.slater_reduced_det(r_new, 0, alpha, number_particles, atom_nr));
                    //wfnew_down = wfold_down;

                    wfnew_up = wfold_up;
                    wfnew_down = wfold_down;

                    r_temp = SetTempR(r_new, atom_nr);
                    wfnew_up(atom_nr) = det(WaveFunc.slater_reduced_det(r_temp, 0, alpha, number_particles, atom_nr));
                }

                else
                {
                    //wfnew_down = det(WaveFunc.slater_reduced_det(r_new, 1, alpha, number_particles, atom_nr));
                    //wfnew_up = wfold_up;

                    wfnew_up = wfold_up;
                    wfnew_down = wfold_down;

                    r_temp = SetTempR(r_new, atom_nr);
                    wfnew_down(atom_nr) = det(WaveFunc.slater_reduced_det(r_temp, 1, alpha, number_particles, atom_nr));
                }

                //Jastrow_new = WaveFunc.Jastrow_Part(r_new, number_particles, beta);
                Jastrow_new = WaveFunc.Jastrow_Part(r_new, N, beta);

                wfnew = 1;
                for (k=0; k<number_atoms; k++)
                {
                    wfnew *= wfnew_up(k);
                    wfnew *= wfnew_down(k);
                }

                //wfnew = wfnew_up * wfnew_down * Jastrow_new;
                wfnew *= Jastrow_new;

                //Oppdaterer kvantekraften
                //qforce_new = (WaveFunc.Jastrow_Gradient_Ratio(r_new, number_particles, beta) + WaveFunc.Gradient_Radio(r_new, alpha, number_particles, atom_nr));
                qforce_new = 2*WaveFunc.Jastrow_Gradient_Ratio(r_new, N, beta);
                //qforce_new = WaveFunc.Jastrow_Gradient_Ratio(r_new, N, beta);
                for (k = 0; k < number_atoms; k++)
                {
                    r_temp = SetTempR(r_new, k);
                    q_temp = WaveFunc.Gradient_Radio(r_temp, alpha, number_particles, k);
                    for (int K = 0; K < number_particles; K++)
                    {
                        qforce_new.row(K+k*number_particles) += 2*q_temp.row(K);
                        //qforce_new.row(K+k*number_particles) += q_temp.row(K);
                    }
                }

                //Bruker green funksjoner
                greensfunction = 0.0;

                for (j=0; j<dimension; j++)
                {
                    greensfunction += (qforce_old(i,j) + qforce_new(i,j)) *
                            (D * timestep * 0.5 * (qforce_old(i,j) - qforce_new(i,j)) - r_new(i,j) + r_old(i,j));
                }
                greensfunction = exp(0.5*greensfunction);

                //Sjekker om vi skal godta nytt steg
                if (ran2(&idum) <= greensfunction*wfnew*wfnew/wfold/wfold)
                {
                    if (cycles > thermalization)
                    {
                        accept+=1;
                    }

                    //r_old = r_new;
                    //qforce_old = qforce_new;

                    for (j=0; j<dimension; j++){
                        r_old(i,j) = r_new(i,j);
                        qforce_old(i,j) = qforce_new(i,j);
                    }
                    wfold_up = wfnew_up;
                    wfold_down = wfnew_down;
                    Jastrow_old = Jastrow_new;
                    wfold = wfnew;
                }
            }

            if (cycles > thermalization)
            {
                //Local Energy
                delta_e = Energi.local_energy(r_old, alpha, wfold, dimension, number_particles, charge, beta, number_atoms);
                energy += delta_e;
                energy2 += delta_e*delta_e;

               // if (isminimizing == false)
                //{
                    out << delta_e << " " << delta_e*delta_e << endl;
                    oga << r_old(0,0) << " " << r_old(0,1) << " " << r_old(5,0) << " " << r_old(5,1) << endl;
                //}
            }
        }

        if (isminimizing == false)
        {
            cout << "alpha = "<< alpha << " beta = " << beta;
            cout << " step length = " << timestep;
            cout << " R = " << R(1,0);
            cout << " Accepted steps: " << (double) accept/(number_cycles*number_particles*number_atoms) << endl;
            cout << "Energy = " << energy/number_cycles << endl;
        }

        //Oppdaterer energien
        cumulative_e[variate2*max_variations + variate] = energy/number_cycles;
        cumulative_e2[variate2*max_variations + variate] = energy2/number_cycles;

    }
    }

}

double SlaterImportantSampling::NormalDist(long *idum)
{
    double u = ran0(idum) * 2 - 1;
    double v = ran0(idum) * 2 - 1;
    double r = u*u + v*v;
    if (r==0 || r > 1){
        return NormalDist(idum);}
    else{
        double c = sqrt(-2 * log(r)/r);
        return u*c;
    }
}

mat SlaterImportantSampling::SetTempR(const mat &r, int at_nr)
{
    mat temp = zeros<mat>(number_particles, dimension);
    for (int F = 0; F < number_particles; F++)
    {
        temp.row(F) = r.row(F+at_nr*number_particles);
    }
    return temp;
}

