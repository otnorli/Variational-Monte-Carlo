#include "minimazation.h"

minimazation::minimazation(double alfha, double betna, int number_variables, int numpar, int charg, double steplen, double hhhh, double hhhh2, int dimdim, int myran, int totsiz)
{
    alpha = alfha;
    dimension = number_variables; //Antall variabler vi varierer er dimensjonen her, ikke x y z som vanlig :O
    number_particles = numpar;
    charge = charg;
    step_lenght = 0.1;
    timestep = steplen;
    h = hhhh;
    h2 = hhhh2;
    max_variations = 1;
    thermalization = 10000;
    cycles = 50000;
    xyzdimension = dimdim;
    beta = betna;
    my_rank = myran;
    size = totsiz;
}

vec minimazation::ConjugateGratient()
{
    double number_atoms = 2;
    double fraction = 0.5;
    R = zeros(number_atoms, dimension);
    R_temp = zeros(number_atoms, dimension);
    int Itermax, i, j, k, m;
    //Rpos = 5.0;

    //Starter ved i=500 for runde 2, trenger flere iterasjoner
    Rpos = 5.00;

    /////////////////

    R(1,0) = Rpos;
    double tempSGA;
    double forskjell;
    double scaling;
    double tolerance = 0.0001;
    vec x_old(dimension);
    m = 10;
    Itermax = 1000;
    vec alpha_liste(Itermax);
    vec beta_liste(Itermax);

    alpha_liste(0) = alpha;
    beta_liste(0) = beta;

    mat SGA_upp = zeros<mat>(m, dimension);
    mat SGA_org = zeros<mat>(m, dimension);
    mat w_i = zeros<mat>(m, dimension);
    vec m_tilde(dimension);
    vec SGA_org_sum(dimension);
    vec SGA_upp_sum(dimension);

    vec SGA_tot_sum(dimension);

    x_old(0) = alpha;
    x_old(1) = beta;
    x_old(2) = Rpos;

    for (i=50; i<Itermax; i++)
    {
        m_tilde(0) = 0;
        m_tilde(1) = 0;
        m_tilde(2) = 0;

        SGA_org_sum(0) = 0;
        SGA_org_sum(1) = 0;
        SGA_org_sum(2) = 0;

        SGA_upp_sum(0) = 0;
        SGA_upp_sum(1) = 0;
        SGA_upp_sum(2) = 0;

        SGA_tot_sum(0) = 0;
        SGA_tot_sum(1) = 0;
        SGA_tot_sum(2) = 0;

        //Vi bruker m walkere og snitter over dette. Trenger da ikke så mange steg per walker
        for (j=0; j<m; j++)
        {
            //cout << "walker nr: " << j << endl;

            //Variansen eller energien
            //0 = energien,
            //1 = variansen
            R(1,0) = x_old(2);
            R_temp(1,0) = x_old(2) - step_lenght;

            //tempSGA = returnEnergy(x_old(1), x_old(0), j, R)(1);

            //tempSGA = returnEnergy(x_old(1), x_old(0), j)(1);

            //Energien i et punkt
            //SGA_org(j,0) = tempSGA;
            //SGA_org(j,1) = tempSGA;
            //SGA_org(j,2) = tempSGA;

            SGA_org(j,0) = returnEnergy(x_old(1), x_old(0) - step_lenght, j, R)(1);
            SGA_org(j,1) = returnEnergy(x_old(1) - step_lenght, x_old(0), j, R)(1);
            //SGA_org(j,1) = returnEnergy(x_old(1), x_old(0), j, R)(0);
            SGA_org(j,2) = returnEnergy(x_old(1), x_old(0), j, R_temp)(0);

            R_temp(1,0) = x_old(2) + step_lenght;

            //Energien i et punkt veldig nær, for deriverte
            SGA_upp(j,0) = returnEnergy(x_old(1), x_old(0) + step_lenght, j, R)(1);
            SGA_upp(j,1) = returnEnergy(x_old(1) + step_lenght, x_old(0), j, R)(1);
            SGA_upp(j,2) = returnEnergy(x_old(1), x_old(0), j, R_temp)(0);

            //Minimere variansen
            //SGA_upp(j,0) = returnEnergy(x_old(1), x_old(0)+step_lenght, j)(1);
            //SGA_upp(j,1) = returnEnergy(x_old(1)+step_lenght, x_old(0), j)(1);

            //Vektene for stokastic gradient method
            w_i(j,0) = SGA_upp(j, 0)*SGA_upp(j,0)/SGA_org(j,0)/SGA_org(j,0);
            w_i(j,1) = SGA_upp(j, 1)*SGA_upp(j,1)/SGA_org(j,1)/SGA_org(j,1);
            w_i(j,2) = SGA_upp(j, 2)*SGA_upp(j,2)/SGA_org(j,2)/SGA_org(j,2);

            //m tilde, summen av alle vektene
            m_tilde(0) += w_i(j,0);
            m_tilde(1) += w_i(j,1);
            m_tilde(2) += w_i(j,2);

            //verdien ganger vekten, summerer opp
            SGA_org_sum(0) += w_i(j,0) * SGA_org(j,0);
            SGA_org_sum(1) += w_i(j,1) * SGA_org(j,1);
            SGA_org_sum(2) += w_i(j,2) * SGA_org(j,2);

            SGA_upp_sum(0) += w_i(j,0) * SGA_upp(j,0);
            SGA_upp_sum(1) += w_i(j,1) * SGA_upp(j,1);
            SGA_upp_sum(2) += w_i(j,2) * SGA_upp(j,2);
        }

        cout << SGA_upp_sum(2) << " " << SGA_org_sum(2) << " ";

        for (k=0; k<dimension; k++)
        {
            SGA_tot_sum(k) = SGA_upp_sum(k) - SGA_org_sum(k);
            SGA_tot_sum(k) = SGA_tot_sum(k) / (m_tilde(k)*step_lenght*2);


            //scaling = 100+i;
            scaling = (double) i;///10000.0;


            forskjell = (double) SGA_tot_sum(k)/scaling;

            if (k==2)
            {
                forskjell *= 3;//5000;
            }

            if (k==0)
            {
                forskjell *= 50;
            }

            if (k==1)
            {
                forskjell *= 1000;
            }

            x_old(k) = x_old(k) - forskjell;

            if (x_old(k) < 0)
            {
                cout << "Look out!";
                x_old(k) = 0;
            }

            cout << x_old(k) << " ";
        }
        cout << endl;

        alpha_liste(i) = x_old(0);
        beta_liste(i) = x_old(1);

        if (abs(SGA_tot_sum(1)) == 0 && abs(SGA_tot_sum(0)) == 0 && (abs(SGA_tot_sum(3)) == 0))
        {
            i = Itermax-1;
        }
    }

    cout << "alpha_liste : ";
    for (i=1; i<Itermax; i++)
    {
        cout << alpha_liste(i) << " ";
    }
    cout << endl << "beta_liste :";

    for (i=1; i<Itermax; i++)
    {
        cout << beta_liste(i) << " ";
    }
    cout << endl;

    return x_old;
}


vec minimazation::returnEnergy(double beta, double alpha, int xaxa, const mat &Rmat)
{
    double *cumulative_e, *cumulative_e2;
    cumulative_e = new double[1];
    cumulative_e2 = new double[1];

    cumulative_e[0] = 0;
    cumulative_e2[0] = 0;

    SlaterImportantSampling SIS(xyzdimension, number_particles, charge, max_variations, cycles, timestep, h, h2, 2, thermalization, Rmat);
    SIS.ImpSamp(cumulative_e, cumulative_e2, beta, alpha, xaxa, true, size);

    vec Values(2);

    //The energy
    Values(0) = cumulative_e[0];

    //The variance
    Values(1) = sqrt((cumulative_e2[0] - cumulative_e[0]*cumulative_e[0])/cycles);

    //If estimating -RMS/E use this
    //Values(1) = -Values(1)/Values(0);

    return Values;
}
