//
//  SBS.cpp
//
//
//  Created by Fyodor Kondrashov on 28.06.19.
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;

#include<map>
#include <string>
#include <iomanip>
using std::setw;

#include <random>
#include <algorithm>
#include <vector>
#include <chrono>
using namespace chrono;

#include <boost/random.hpp>

#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/skew_normal.hpp>

#include <cstdlib>
using std::atoi;
using std::atof;
#include <ctime>
using std::time;
#include <boost/random/variate_generator.hpp>
using boost::variate_generator;

#define GENERATIONS 30000
#define MAXPOPSIZE 5000
#define GENOMESIZE 5000
#define MU 0.0000012

//USAGE: popsize, num_muts
//RESULT: introduces ave(lambda) poisson dist mutations to pop
void mutation (int, int);

// USAGE: old population size, new population size
//RESULT: creates newpop from pop and swaps them (pop becomes newpop, newpop becomes pop)
void reproduction (int, int);

// USAGE: current population size
//RESULT: creates newpop from pop and swaps them (pop becomes newpop, newpop becomes pop)
void update_fitness (int);

//USAGE: current popsize, number of alleles per genome (numhetero). and number of homozygous sites. heterozygosity
//RESULT:calculates heterozygosity and other stats about the pop (# of alleles per genome)
// also updates fitness (fit) by removing the impact of fixed alleles (necessary to keep
// the same distribution of s in muts)
int get_stats (int, double*, double*, double*, double*);

//initialize random seed

/* random_device MAY BE CRYPTOGRAPHICALLY STRONG */
random_device rd;
/* SEED mt19937 and mt19937_64 ENGINE WITH STRONG RANDOM NUMBER */
mt19937 mt(rd());
/* CREATE RANDOM DISTRIBUTION OBJECTS FOR INTS AND DOUBLES */
uniform_real_distribution<double> real_dist(0, 1.0);

double **muts;
double *fit;

std::map<int, std::map<int, int> > pop;
std::map<int, std::map<int, int> > newpop;
std::map<int, int> haplo;


int main (int argc, char** argv)
{
    
    // These are variable definitions
    long i, gener, popsize;
    int num_muts, fixed, size;
    double *numhetero, *numhomo, *heteroz, *het_pq, rand, h_ave, h_std, mu, sigma, alpha;
    char* dt;

    if (not argv[5]){cout << "Usage: ./RNS mu sigma alpha h_average outputfilename\n Remember that sigma must be positive\n"; return 0;}
    
    mu    = atof(argv[1]); //-0.01;
    sigma = atof(argv[2]); //0.00001;
    alpha = atof(argv[3]); //-3;
    h_ave = atof(argv[4]); //0.1;
    h_std = 0.3;
    popsize = MAXPOPSIZE;

    ofstream outfile;
    outfile.open (argv[5]);

    std::cout.precision(10);
    
    //do not forget to assign memory to all the variables
    
    fixed = 0;
    
    numhetero = new double [1];
    numhomo   = new double [1];
    heteroz   = new double [1];
    het_pq    = new double [1];
    
    numhomo[0]=numhetero[0]=heteroz[0]=het_pq[0]=NULL;
    
    //initialize poisson_distribution variables - redo later if popsize changes
    
    double lambda = MU * 2 * popsize * GENOMESIZE; //do not forget about POPSIZE
    
    poisson_distribution<int> pdist(lambda);
    variate_generator< mt19937, poisson_distribution<int> > rvt(mt, pdist);
    
    // initialize the population, fitness and mutations array (muts array is for s and h
    
    // initialize fitness array
    fit = new double [MAXPOPSIZE];
    for(i = 0; i < MAXPOPSIZE; i++)
    {fit[i]=100;}
    
    // DISTRIBUTION OF S AND SH TO MUTS ARRAY
    
    auto skew_norm_dist = boost::math::skew_normal_distribution<double>(mu, sigma, alpha);
    std::normal_distribution<> norm(h_ave,h_std);

    // initialize mutation array [0] for s and [1] for sh
    muts = new double* [GENOMESIZE];
    for(i = 0; i < GENOMESIZE; i++)
    {   muts[i] = new double [2];
        
        muts[i][0]=boost::math::quantile(skew_norm_dist, real_dist(mt));
             if (muts[i][0] < -0.1){muts[i][0] = -0.1;}
        else if (muts[i][0] > 0.1){muts[i][0] = -0.1;}
        
        muts[i][1] = norm(mt);
        while ((muts[i][1] >= 1) or (muts[i][1] <= 0)){muts[i][1] = norm(mt);} //make this a truncated normal

        muts[i][1] = muts[i][1] * muts[i][0];
    }
    
    for (i = 0; i < 10; i++)
    {
        
        time_t now = time(0);
        dt = ctime(&now);
        outfile << "Run " << i << " start time : " << dt << endl;

        for (gener = 0; gener < GENERATIONS; gener++)
        {
        
            num_muts  = rvt();
            
            mutation(popsize, num_muts);
        
            reproduction(popsize, popsize);
        
            update_fitness(popsize);
        
            fixed += get_stats (popsize, numhetero, numhomo, heteroz, het_pq);
        
            // print results
            outfile << " generation\t" << gener << "\tnum_muts\t" << num_muts;
            outfile << "\tfixed\t" << fixed << "\tnumhomo\t" << numhomo[0] << "\tnumhetero\t" << numhetero[0];
            outfile <<  "\theterozygosity_2pq\t" << heteroz[0] << "\theterozygosity_1_p2_q2\t" << het_pq[0] << endl;
        
    }
        time_t now1 = time(0);
        dt = ctime(&now1);
        outfile << "Run " << i << "time at the end: " << dt << endl;
    }
    
    outfile.close();
    return 0;
    
}

int get_stats (int popsize, double* numhetero, double* numhomo, double* heteroz, double* het_pq)
{
    int i, fixed;
    double q, temp, counthomo, counthet;

    std::map<int, std::map<int, float> > homo_het;

    numhetero[0]=numhomo[0]=heteroz[0]=het_pq[0]=fixed=q=0;
    
    //count the number of sites per genome and make homo_het
    for(auto const &ent1 : pop)
    {counthet=counthomo=0;
        for (auto const &ent2 : ent1.second)
        {
                 if (ent2.second == 2) {homo_het[ent2.first][0]++;counthomo++;} //  if 11
            else if (ent2.second == 1) {homo_het[ent2.first][1]++; counthet++;} //  if 01
        }
        numhomo[0]   += counthomo/GENOMESIZE;
        numhetero[0] += counthet/GENOMESIZE;
    }
    
    // count fixed sites and get them back to 00
    for(auto const &genes : homo_het)
    {
        if (homo_het[genes.first][0] == popsize)
        { fixed++;
            for(i = 0; i < popsize; i++) {pop[i].erase(genes.first);} //return to 00
        }
    }
    
    // calculate heterozygosity
    //first as a fraction of heterozygous sites
    for(auto const &genes : homo_het)
    { if (homo_het[genes.first][1]){heteroz[0] += homo_het[genes.first][1]/popsize;}}
    
    heteroz[0] /= GENOMESIZE;

    // then as 1 minus homozygous sites
    for(auto const &ent1 : homo_het)
    {
            q = (2 * homo_het[ent1.first][0] + homo_het[ent1.first][1]) / (2 * popsize);
            het_pq[0] += 1 - pow(q,2) - pow((1-q),2);
    }
    
    het_pq[0] /= GENOMESIZE;
    
    homo_het.clear();
    
    return fixed;
}

void update_fitness (int popsize)
{

    int i;
    
    for(i = 0; i < popsize; i++)
    {fit[i]=100;}
    
    
    for(auto const &ent1 : pop)
    {
        for (auto const &ent2 : ent1.second)
        {
            if (ent2.second == 2) // if 11
            {fit[ent1.first] += fit[ent1.first]*muts[ent2.first][0];} //  s
            else // if 01
            {fit[ent1.first] += fit[ent1.first]*muts[ent2.first][1];} //  sh
        }
    }
    
}

void reproduction (int popsize, int newpopsize)
{
    bool selected;
    double maxfit, selection;
    int i, mother, father, haplotype;
    
    mother=father=haplotype=selection=0;
    
    //calculate the maximum fit value in the population
    
    maxfit = -999; //start with very low fitness
    for (i=0; i < popsize; i++)
    {if (fit[i] > maxfit){maxfit=fit[i];}}
    
    if (maxfit < 0) {maxfit=0;}
    
    //populate the new population
    
    newpop.clear();
    
    for (i=0; i < newpopsize; i++)
    {
        selected = 0;
        while (not selected) //select mother
        {   mother = real_dist(mt) * popsize;
            selection = real_dist(mt) * maxfit;
            
            if (fit[mother] >= selection){selected = 1;}
        }
        
        selected = 0;
        while (not selected) //select father
        {   father = real_dist(mt) * popsize;
            selection = real_dist(mt) * maxfit;
            
            if (fit[father] >= selection and mother != father){selected = 1;}
        }
        
        // give half of mother
        haplo = pop[mother];
        for(auto const &ent1 : haplo)
        {
                if (ent1.second == 2) // mother is 11
                    {newpop[i][ent1.first] = 1;} //baby gets a 1
                else if (ent1.second == 1) //mother is 01
                {haplotype = real_dist(mt)*2;
                if (haplotype == 1) // 0 -> 1
                    {newpop[i][ent1.first] = 1;}
                }
        }
        
        haplo.clear();
        
        //give half of father
        haplo=pop[father];
        for(auto const &ent1 : haplo)
        {
                if (ent1.second == 2) // father is 11
                {newpop[i][ent1.first]++;}//x0 -> x1
                else if (ent1.second == 1) // father is 01
                {haplotype = real_dist(mt)*2;
                    if (haplotype == 1) // 0 -> 1
                    { newpop[i][ent1.first]++;} //x0 -> x1
                }
        }
    }
    pop.clear();
    pop = newpop;
    newpop.clear();
}

void mutation (int popsize, int num_muts)
{
    int i, indiv, gene, parent;
    
    //number of mutations in a generation - from poisson dist in main
    
    //this function also automatically updates the fitness array
     
    for (i=0; i < num_muts; i++)
    {
        indiv   = real_dist(mt) * popsize;
        gene    = real_dist(mt) * GENOMESIZE;
        
             if (pop[indiv][gene] == 2) // 11 -> 01 mutation
                {pop[indiv][gene]=1;
                 fit[indiv] += fit[indiv]*(1/(1+muts[gene][0]) - 1); // r=1/(1+y) -1, revert to 00 fitness
                 fit[indiv] += fit[indiv]*muts[gene][1]; // add impact of sh
                }
            else if (pop[indiv][gene] == 1) // 01 -> 00 or 11 mutation
                {parent = real_dist(mt) * 2;
                    if (parent == 1)  // 01 -> 00
                      {pop[indiv][gene]=0;
                      fit[indiv] += fit[indiv]*(1/(1+muts[gene][1]) - 1);
                      }
                    else  // 01 -> 11
                      {pop[indiv][gene]=2;
                       fit[indiv] += fit[indiv]*(1/(1+muts[gene][1]) - 1); // r=1/(1+y) -1, revert to 00 fitness from 01 fitness
                       fit[indiv] += fit[indiv]*muts[gene][0]; // add impact of s
                      }
                }
            else // 00 -> 01 mutation
                {pop[indiv][gene]=1;
                fit[indiv] += fit[indiv]*muts[gene][1];}
        }
}

/*   cout << "pop after reproduction: " << endl;
 
 
 for(auto  &ent1 : pop)
 {
 for (auto  &ent2 : ent1.second)
 {
 cout << ent1.first << " " << ent2.first << " " << ent2.second << endl;
 }
 } */

/*    cout << "fit after update: " << endl;
 
 for (i = 0; i < popsize; i++)
 {  cout << i << " " << fit[i] << endl;}
 
 */
