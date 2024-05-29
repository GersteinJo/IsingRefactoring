#include<iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <numeric>
#include <fstream>
#include <filesystem>

class ising_grid
{
private:

	int L, N_mc, N_sample, N_therm; // grid params
	int num_T, num_h;
    float h=0., h_max = 8.;//magnetic field range
	float T = 0., T_max = 4.; //temperature range
	float* temperature; // temperature array
    float* field;
	float* s; //average spin array
    float* Gamma; // sum of spin-spin interaction
	float** grid; // grig for calculating the average

	//randomizer
	long long seed;
	std::mt19937 generator;
	std::uniform_real_distribution<> distribution;

public:
	ising_grid(int L, int N_mc, int N_sample, int N_therm, int num_T, int num_h) :
		L(L), N_mc(N_mc), N_sample(N_sample), N_therm(N_therm), num_T(num_T), num_h(num_h)
	{
		grid = new float* [L];
		for (int i = 0; i < L; i++) grid[i] = new float[L];
		for (int i = 0; i < L; i++)
		{
			for (int j = 0; j < L; j++)
			{
				grid[i][j] = 1;
			}
		}

		temperature = new float[num_T];
        for (int i = 0; i < num_T; i++) temperature[i] = T_max / num_T * (i + 1); // for T = 0 we know the polarisation

        field = new float[num_h];
        for (int i = 0; i < num_h; i++) field[i] = h_max / num_h * (i + 1); 

		s = new float[num_T];
        Gamma = new float[num_T];

		seed = std::chrono::steady_clock::now().time_since_epoch().count();
		generator = std::mt19937(seed);
		distribution = std::uniform_real_distribution<>(0, 1);
	}


    void visualize()
    {
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < L; j++)
            {
                grid[i][j] >0 ? std::cout<<"^ ": std::cout << "v ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

	void ising_simulate_temp(float h = 0)
	{
		float** grid_1 = new float* [L];
		for (int i = 0; i < L; i++) grid_1[i] = new float[L];

        for (int k = 0; k < num_T; k++)
        {
            int dim = (N_mc - N_therm) / N_sample;
            float* p_mean = new float[dim];
            float* gamma = new float[dim];

            for (int t = 0; t < dim; t++) p_mean[t] = 0.0;

            for (int a = 0; a < N_mc; a++)
            {
                int p = 0, g = 0;
                for (int i = 0; i < L; i++)
                {
                    for (int j = 0; j < L; j++)
                    {
                        int i_, i_i, j_, j_j;
                        i + 1 == L ? i_ = 0 : i_ = i + 1;
                        j + 1 == L ? j_ = 0 : j_ = j + 1;
                        i == 0 ? i_i = L - 1 : i_i = i - 1;
                        j == 0 ? j_j = L - 1 : j_j = j - 1;
                        float delta_E = 2 * grid[i][j] * (grid[i_][j] + grid[i_i][j] + grid[i][j_] + grid[i][j_j]-h);

                        if (distribution(generator) < std::fmin(std::exp(-delta_E / temperature[k]), 1))
                        {
                            grid_1[i][j] = -grid[i][j];
                        }
                        else
                            grid_1[i][j] = grid[i][j];
                    }
                }

                grid = grid_1;

                for (int i = 0; i < L; i++)
                {
                    for (int j = 0; j < L; j++)
                    {
                        p += grid[i][j];

                        for (int i_ = 0; i_ < L; i_++)
                        {
                            for (int j_ = 0; j_ < L; j_++)
                            {
                                g += grid[i][j] * grid[i_][j_];
                            }
                        }
                    }
                }



                if ((a - N_therm) % N_sample == 0 && (a - N_therm) >= 0)
                {
                    int n = (a - N_therm) / N_sample;
                    p_mean[n] = p;
                    gamma[n] = g;
                }
            }
            //visualize();

            float p_sum = 0, g_sum = 0;

            for (int t = 0; t < dim; t++)
            {
                p_sum += p_mean[t];
                g_sum += gamma[t];
            }

            s[k] = p_sum / (float)(N_mc - N_therm) * (float)N_sample;
            Gamma[k] = g_sum / (float)(N_mc - N_therm) * (float)N_sample;
            
        }
	}

    float* spin()
    {
        return s;
    }

    double* chi()
    {
        double* chi = new double[num_T];
        //here we insert eq 2.10
        for (int i = 0; i < num_T; i++)
        {
            chi[i] = 1 / (double)temperature[i] * ((double)Gamma[i] - (double)s[i] * (double)s[i]); // we devide by L^2 in the output
//            if (chi[i] == 0.) std::cout<<Gamma[i]<<" "<<s[i]<<"\n";
        }
        return chi;
    }

    float* t()
    {
        return temperature;
    }
};

int main()
{
    int num_T = 80, num_h = 80;
    int L = 20, N_mc = 1200, N_sample = 20, N_therm = 600;
    ising_grid my_grid(L, N_mc, N_sample, N_therm, num_T, num_h);
    my_grid.ising_simulate_temp();
    float* polarisation = my_grid.spin();
    double* chi = my_grid.chi();
    //for (int i = 0; i < num_T; i++) std::cout << polarisation[i] << "\n";


    std::fstream my_file_stream("chi_20x20.txt", std::ios::out);
    for (int i = 0; i < num_T; i++) my_file_stream << chi[i]/L/L << "\n";

    my_file_stream.close();

    //my_grid.ising_simulate_magn();
    //float* magnetization = my_grid.magnetization();

    //std::fstream my_file_stream("polarisation_20x20 from H temp 0.5.txt", std::ios::out);
    //for (int i = 0; i < num_h; i++) my_file_stream << magnetization[i] / L / L << "\n";
    //my_file_stream.close();

}