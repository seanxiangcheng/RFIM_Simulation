#include "mtrand.h" //random number generator 
#include <cstdio>
#include<iostream>
#include<fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main() {

double end=5000; //upper bound for x and y coordinate
double begin=0;//lower bound
double length_data=(end-begin);
double rando_x; //randomly generated x coordinate
double rando_y; //randomly generated y coordinate

ofstream out; 
///////////////////////Initialization of the random number generator
unsigned long init[4] = {0x123, 0x234, 0x345, 0x456}, length = 4;
MTRand drand; // double in [0, 1) generator, already init

FILE *pipe = popen("gnuplot", "w"); //this opens gnuplot

fprintf(pipe, "set xrange [0:5000]\n"); //telling gnuplot the range for axes
fprintf(pipe, "set yrange [0:5000]\n");


for (int j=0;j<5000;j++) //loops through 5000 time steps
{

out.open("gp_rand_test.txt"); 

for (int i = 0; i < 1000; ++i) //This loop generates the set of randomly placed points
{
rando_x=(drand()*length_data+begin);
rando_y=(drand()*length_data+begin);
out<<rando_x<<"\t"<<rando_y<<endl;
}

out.close();
out.clear();


fprintf(pipe, "plot 'gp_rand_test.txt' using 1:2\n"); //plotting the .txt file 



}

getchar(); //This line keeps the gnuplot window open after the code runs through.
return 0;

} 
