// Original Compile: g++  -Wall  -lm -lgsl -lgslcblas -c  "%f"
// Original Build: g++ RFIM1.cc -Wall   -lm -lgsl -lgslcblas -o RFIM1.x 
// Original run: ./RFIM1.x-L 32  -g 1.0  -T 1.0   -r 2 -s 615 -n seg -M 1000
// Compile: g++  -Wall  -lm -lgsl -lgslcblas -c  "%f"
// Build: g++ "%f" -Wall   -lm -lgsl -lgslcblas -o "%e" 
// Run: "./%e"  -k 17 -0 100  -1 100 -y 0

/* Code: 2D Random Field Ising Model
 * Input: (default values)
 *        system size "-L": L = 16 (length), W = L (width)
 *        temperature "-T" T = 2
 *        variance of hi "-g": s2hi = 1
 *        configuration output MC sweeps interval "-M": dMCs = 100.2
 *        total MC sweeps: total_MCs = 1000 
 *        repeat_num "-r": number of repetitions repeated_num = 1 
 *          (each repetition will re-run the simulation and export new outputs)
 *        NAME "-n": name = "" ; any string with length < 8 to help you name the output csv files;
 *        seed "-s": seed=getpid() (usually 1 or 2); 
 * ### theory of fractal 
 * Output:
 *        configurations: config_L***_T***_r**_NAME.csv
 *                      1st column: MC step;
 *                      Each row except 1st column: Spin configurations
 *        Magnetizations: mag_L***_T***_r**_NAME.csv
 *                      1st column: MC step;
 *                      2nd column: average magnetization (M/N);
 *        Energy: e__L***_T***_r**_NAME.csv
 *                      1st column: MC step;
 *                      2nd column: average energy (E/N);
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <term.h>
#include <ncurses.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <unistd.h>
#include <cstdio>
#include<fstream>

using namespace std;

#define L_DEF 1000
#define L_MAX 10000
#define T_DEF 2.0
#define S_DEF 1.0
#define NAME_LEN 8
#define DIM 2
#define JCON -1

gsl_rng *rnd = gsl_rng_alloc(gsl_rng_mt19937); // global random number setup

class RFIM2d {
  private:
    int L, N;     // system size N=L^2
    double s_hi;  // 
    int *spins;   // length is N; mememory is allocate later
    double *his;  // static internal random field 
    double T;     // temperature
    double E;     // energy
    int M;        // magnetization
    
  public:
    int check_E();
    int check_M();
    int get_L() { return(L);}
    int get_N() { return(N);}
    int get_M() { return(M);}
    int *get_spins() {return(spins);}
    double get_E() { return(E); }
    double get_T() { return(T); }
    int initialize();  // initialization with default values
    int *neighbors(int iv);
    int neighbor_check(int iv);
    void set_L (int Lp = L_DEF) {L = Lp; N=L*L;} // set system size
    void set_s_hi (double sp = S_DEF) {s_hi = sp;}
    void set_T (double Tp = T_DEF) {T = Tp;}
    int set_spins ();
    int set_his ();
    double calc_E();
    int calc_M();
    int print_info(int config);
    int update();
};

double RFIM2d::calc_E()
{
    double energy_s = 0.0;
    double energy_h = 0.0;
    double energy = 0.0;
    int *ng; 
    for (int i = 0; i < N; i++)
    {
        ng = neighbors(i);
        for (int j = 0; j < 2*DIM; j++)
        {
            energy_s = energy_s + JCON * spins[i] * spins[ng[j]];   
        }
        energy_h = energy_h + his[i] * spins[i];
    }
    energy = energy_s / 2.0 + energy_h;
    
    return(energy);
}

int RFIM2d::calc_M()
{
    int mag = 0; 
    for (int i = 0; i < N; i++)
    {
        mag = mag + spins[i];
    }
    return(mag);
}

int RFIM2d::check_E()
{
    double energy = calc_E();
    if((energy - E) < 1e-12)
    {
        cout << "  Energy is correct: E = calc_energy = " << energy << endl;
    }
    else
    {
        cout << "  Energy is wrong: E = " << E << endl;
        cout << "  But --------calc_E = " << energy << endl;
        E = energy;
    }
    return(0);
}

int RFIM2d::check_M()
{
    int mag = calc_M();
    if(M == mag)
    {
        cout << "  Magnetization is correct: M = calc_M = " << M << endl;
    }
    else
    {
        cout << "  Magnetization is wrong: M = " << M << endl;
        cout << "  But -------------- calc_M = " << mag << endl;
        M = mag;
    }
    return(0);
  
}

int * RFIM2d::neighbors(int iv){
    int x, ld=1;  
    int direction=-1, n=0;
    int d=0;
    static int neighs[DIM*2];
    
    for(d=0; d<DIM; d++){
      for(direction=-1; direction<2; direction += 2){
        ld=1;
        if(d>0) ld *= L;
        x = (iv/ld)%L;
        neighs[n++]=(iv+((x+L+direction)%L - x)*ld);
      }
    }
    return(neighs);
}



int RFIM2d::neighbor_check(int iv)
{
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }
    int *ng = neighbors(iv);
    cout << "   Neighbors of site " << iv << endl;
    for (int i = 0; i < 2*DIM; i++)
    {
        cout << "   " << ng[i] << endl;
    }
    return(1);
}


// initialization with default values
int RFIM2d::initialize() 
{
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }
    set_spins();
    set_his();
    
    E = calc_E();
    M = calc_M();
    return(1);
}

int RFIM2d::set_spins () 
{  
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }    
    spins = new int[N];
    for(int i=0; i<N; i++)
    {
        if(i < N/2) { spins[i] = -1; }
        else { spins[i] = 1; }
    }
    return(1);
}

int RFIM2d::set_his()
{
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }
    his = new double[N];
    for(int i=0; i<N; i++)
    {
        his[i] = gsl_ran_gaussian(rnd, s_hi);
        //his[i] = 0.0; //del
    }
    return(1);
}

int RFIM2d::print_info(int config=0)
{
    int i; 
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        System not set up yet!\n", L, N);
        return(0);
    }
    printf("*** Simulation System Information ***\n");
    printf("   L=%d; N=%d;\n   T=%.2f; hi_sigma=%.2f\n   E=%.2f; M=%d;\n"\
    , L, N, T, s_hi, E, M);
    if(config != 0)
    {
        printf("   ** Configurations **\n");
        for(i=0; i<N; i++)
        {
            printf("    site %-4d: %-4d, %-.4f\n", i, spins[i], his[i]);
        }
    }
    return(1);
}

int RFIM2d::update()
{
    int site;
    site = N*gsl_rng_uniform(rnd);
    //printf("site:%d\n", site);//del
    double dE = 0.0;
    int *ng = neighbors(site);
    for (int i=0; i < 2*DIM; i++)
    {
        dE = dE - 2.0 * JCON * spins[site] * spins[ng[i]];
    }
    dE = dE - 2.0 * his[site] * spins[site];  //??????
    if(dE < 0.0 || gsl_rng_uniform(rnd) < exp(-(dE + gsl_ran_gaussian(rnd, sqrt(T))*spins[site])/T))
    {
        spins[site] = -spins[site];
        E = E + dE;
        M = M + 2*spins[site];
    }
    return(1);
}




/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, RFIM2d *age, double *total_MCs, int *repeat_num, int *seed, char *name)
{
  	int i;
  	*seed = getpid();
  	for (i = 1; i < argc; i++)
    {   //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-')
        {
      			switch (argv[i][1])
            {
      				case 'L':       age->set_L( atoi(argv[++i]) );
        			break;
      				case 's':       *seed = atoi(argv[++i]);
        			break;
              case 'g':       age->set_s_hi( atof(argv[++i]) );
              break;
      				case 'r':       *repeat_num = atoi(argv[++i]);
        			break;
      				case 'T':       age->set_T( atof(argv[++i]) );
        			break;
      				case 'n':       strncpy(name, argv[++i], 7);
        			break;
              case 'M':       *total_MCs = atof(argv[++i]);
              break;
              case 'd':       age->set_L();
                              age->set_T();
                              age->set_s_hi();
              break;
      				default:
        			fprintf(stderr,"\nError:  Incorrect option %s\n",argv[i]);
        			fprintf(stderr,"\n\
              Available options: \n\
              -L = L (Length; N=L^2; default=16)\n\
              -g = variance for hi (default=16)\n\
              -n = name (default='')\
              -r = repeat number (default= 1)\n\
              -s = seed: random seed  (default=getpid())\n\
              -t = total MC sweeps (default=1000)\
              \n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}

int main (int argc, char *argv[]) {
  /************* System Setup ****************/
    int repeat_num = 1, iter = 0, seed, *spins; // iter=0:repeat_num-1
    RFIM2d age;  // simulation main object
    char name[NAME_LEN + 1];
    for(int i=0; i <= NAME_LEN; i++) {name[i] = '\0';}
    double total_MCs = 100.0, total_rnd_steps = 0.0, total_rnd_steps_ran=0.0, steps_vis = 0.0;
    
    Commandlineparse(argc, argv, &age, &total_MCs, &repeat_num, &seed, name);
    gsl_rng_set(rnd, (unsigned long int)seed); // set up random seed before initialize()
    total_rnd_steps = total_MCs * age.get_N();
    age.initialize();
    printf("total_MC=%.2f; repeat_num=%d; seed=%d; name=%s (len=%d)\n", total_MCs, repeat_num, seed, name, (int)strlen(name));
    // Test input:
    age.print_info();
    //age.neighbor_check(15);
    age.check_E();
    age.check_M();
    ofstream out; 
    FILE *pipe = popen("gnuplot", "w");

    fprintf(pipe, "set xrange [-0.5:31.5]\n"); //telling gnuplot the range for axes
    fprintf(pipe, "set yrange [-0.5:31.5]\n");
    fprintf(pipe, "set palette defined (-1 \"blue\", 1 \"black\") \n");
    for (iter = 0; iter < repeat_num; iter++)
    {
        /*********** System Initialization *************/
        total_rnd_steps_ran = 0.0;
        age.initialize();
        /**************** Simulation *******************/
        while (total_rnd_steps_ran < total_rnd_steps)
        {
            age.update();
            total_rnd_steps_ran += 1.0;
            steps_vis += 1.0; 
            if(steps_vis > 10)
            {
                //printf("Open coe: %d\n", out.open("spins.txt")); 
                spins = age.get_spins();
                for (int i = 0; i < age.get_N(); i++)
                {
                    out<< i % age.get_L() <<"\t"<< i/age.get_L() << "\t" << spins[i] << endl;
                }
                out.close();
                out.clear();
                fprintf(pipe, "plot 'spins.txt' using 1:2:3 with image pixels; \n");

                //cout << "plot: " << iter << endl;
                steps_vis = 0.0;
            }
            //if(total_rnd_steps_ran-1.0<1e-12) { sleep(1);}
        }
        
        //age.check_E();
        //age.check_M();

    }
    
    /**************** Output ***********************/
    getchar(); //This line keeps the gnuplot window open after the code runs through.
    return(0);
}







/* Neighbor test
    int *de;
    int ss[] = {1, 6, 9, 10};
    for(i=0; i<4; i++)
    {
        de = age.neighbors(ss[i]);
        for(j=0; j<2*DIM; j++)
        {
            printf(" %-3d:%-3d\n", ss[i], de[j]);
        }
    }
    */
