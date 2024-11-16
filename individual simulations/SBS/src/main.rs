use rand::prelude::*;
use rand::Rng;
use rand_distr::{Poisson, Distribution, SkewNormal, Normal, Gamma};
use std::collections::HashMap;
use std::process;
use std::fs::File;
use std::io::prelude::*;
use std::env;

use chrono::{Local, Datelike, Timelike};

fn main() {

    // Collect the command-line arguments into a vector of strings
    let args: Vec<String> = env::args().collect();

    // Ensure the correct number of arguments is provided
    if args.len() != 4 {
        eprintln!("Usage: {} <arg1> <arg2> <arg3>", args[0]);
        std::process::exit(1);
    }

    //these variables define the gamma distribution and fraction of lethals
    let lethals:     f64 = args[1].parse().expect("Failed to parse arg1 as f64");
    let gamma_alpha: f64 = args[2].parse().expect("Failed to parse arg2 as f64");
    let gamma_beta:  f64 = args[3].parse().expect("Failed to parse arg3 as f64");

    
    // Convert the floats to strings for use in filenames
    let arg1_str = format!("{:.2}", lethals);
    let arg2_str = format!("{:.10}", gamma_alpha);
    let arg3_str = format!("{:.10}", gamma_beta);

    // Create filenames using the float values
    let filename1 = format!("SBS_fitness_before_{}_{}_{}.txt",  arg1_str, arg2_str, arg3_str);
    let filename2 = format!("SBS_frequency_before{}_{}_{}.txt", arg1_str, arg2_str, arg3_str);
    let filename3 = format!("SBS_fitness_after{}_{}_{}.txt",    arg1_str, arg2_str, arg3_str);
    let filename4 = format!("SBS_frequency_after{}_{}_{}.txt",  arg1_str, arg2_str, arg3_str);
    let filename5 = format!("SBS_fitness_after_bneck{}_{}_{}.txt",    arg1_str, arg2_str, arg3_str);
    let filename6 = format!("SBS_frequency_after_bneck{}_{}_{}.txt",  arg1_str, arg2_str, arg3_str);
    let filename7 = format!("SBS_{}_{}_{}.txt",                 arg1_str, arg2_str, arg3_str);
 

    // Create the files
    let mut file_fit_before        = File::create(&filename1).unwrap();
    let mut file_freq_before       = File::create(&filename2).unwrap();
    let mut file_fit_after         = File::create(&filename3).unwrap();
    let mut file_freq_after        = File::create(&filename4).unwrap();
    let mut file_fit_after_bneck   = File::create(&filename5).unwrap();
    let mut file_freq_after_bneck  = File::create(&filename6).unwrap();
    let mut outfile                = File::create(&filename7).unwrap();


for runs in 0..100
{

writeln!(outfile, "gen\tpopsize\tfixed\ttotal_1\thet2pq\tave_homo\thet1p2q2\tsample_total_1\tsample_het2pq\tsample_ave_homo\tsample_het1p2q2").expect("Failed to write to outfile");

let mut popsize: usize = 8000;
let genomesize:  usize = 300000;
let mut_rate           = 0.00000012;
let sample             = 7;

let expect_muts = mut_rate * genomesize as f64 * 2.0 * 10000.0; //this is the lambda parameter in Poisson

//DO NOT FORGET ABOUT POPSIZIE - AS POPSIZE CHANGES, SO MUST THE EXPECTED NUMBER OF MUTATIONS.

//this tuple carries the variables needed to define s and sh, in this order: my, sigma, alpha and beta.
// beta is used for the exponential to define sh, the rest are for the skewed normal distribution
let skew: (f64, f64, f64, f64) = (-0.0003, 0.0002, -4.0, 1000.0/10.0);

/*
lethals: 0.0: alpha=1.00000000e+00,  beta=5.92687785e-05
lethals: 0.1: alpha=1.00000000e+00 , beta=4.36764895e-05
lethals: 0.2: alpha=1.00000000e+00 , beta=2.77869816e-05
lethals: 0.3: alpha=1.00000000e+00 , beta=1.13052239e-05
*/


let mut fixed: usize   = 0;

let mut total_1   = 0.0;
let mut ave_hetero= 0.0;
let mut ave_homo= 0.0;
let mut het_1p2q2 = 0.0;

let mut sample_total_1   = 0.0;
let mut sample_ave_hetero= 0.0;
let mut sample_ave_homo= 0.0;
let mut sample_het_1p2q2 = 0.0;

// define different random distributions
let skew_normal = SkewNormal::new(skew.0, skew.1, skew.2).unwrap();
let gamma = Gamma::new(gamma_alpha, gamma_beta).unwrap();
let normal      =     Normal::new(-0.01, 0.005).unwrap();
let poisson = Poisson::new(expect_muts).unwrap();
let mut rng = rand::thread_rng();

/*define the population of popsize individuals but when the length of the genome u16 - this is 65536 sites with bool alleles (empty,0,1) */
let mut pop: Vec<HashMap<u32, u8>> = Vec::with_capacity(popsize);
pop.resize(popsize, HashMap::new());

/*for each individual in the population this array keeps track of its fitness */
let mut fitness: Vec<f64> = Vec::with_capacity(popsize);

/*for each site in the genome this array keeps track of s and sh - this is used to define fitness*/
let mut selection: Vec<(f64, f64)> = vec![(0.0, 0.0); genomesize];

selection_distribution (&mut selection, &genomesize, &skew.3, &skew_normal, &normal, &gamma, &lethals);

//print_selection(&selection);

popsize = 2000;
fitness.resize(popsize, 1.0);

let mut newpopsize = popsize;

//burnin
for i in 0..30000
{

let num_muts = poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0;

mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);

(pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
= reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);

writeln!(outfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, newpopsize, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero,sample_ave_homo, sample_het_1p2q2).expect("Failed to write to outfile");

update_fitness (& pop, &selection, &mut fitness, &popsize);
}

let formatted_runs = format!("{}", runs);
file_fit_before.write_all(b"\t").unwrap();
file_fit_before.write_all(formatted_runs.as_bytes()).unwrap();
file_fit_before.write_all(b"\tbefore\n").unwrap();

file_freq_before.write_all(formatted_runs.as_bytes()).unwrap();
file_freq_before.write_all(b"\tbefore\n").unwrap();

print_stats(&pop, &selection, &popsize, &genomesize, &sample, &mut file_fit_before, &mut file_freq_before); 


// Copy the state of the population before the complex demographic changes
let pop_saved: Vec<HashMap<u32, u8>> = pop.iter().cloned().collect();

//growth
for i in 0..1000
{

newpopsize = popsize + 6;

let num_muts = poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0;

mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);

(pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
= reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);

writeln!(outfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, newpopsize, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero,sample_ave_homo, sample_het_1p2q2).expect("Failed to write to outfile");

popsize = newpopsize;

update_fitness (&pop, &selection, &mut fitness, &newpopsize);
}

//stable
for i in 0..200
{

let num_muts = poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0;

mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);

(pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
= reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);

writeln!(outfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, newpopsize, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero,sample_ave_homo, sample_het_1p2q2).expect("Failed to write to outfile");

update_fitness (& pop, &selection, &mut fitness, &popsize);
}

//decline
for i in 0..565
{

newpopsize = popsize - 14;

let num_muts = poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0;

mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);

(pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
= reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);

writeln!(outfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, newpopsize, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero,sample_ave_homo, sample_het_1p2q2).expect("Failed to write to outfile");

popsize = newpopsize;

update_fitness (& pop, &selection, &mut fitness, &popsize);
}

let formatted_runs = format!("{}", runs);
file_fit_after.write_all(b"\t").unwrap();
file_fit_after.write_all(formatted_runs.as_bytes()).unwrap();
file_fit_after.write_all(b"\tafter\n").unwrap();

file_freq_after.write_all(b"\t").unwrap();
file_freq_after.write_all(formatted_runs.as_bytes()).unwrap();
file_freq_after.write_all(b"\tafter\n").unwrap();

print_stats(&pop, &selection, &popsize, &genomesize, &sample, &mut file_fit_after, &mut file_freq_after); 


// Restart with the population saved before the demographic increase
pop.clear();
pop.extend(pop_saved.iter().cloned());
popsize = 2000;
update_fitness (& pop, &selection, &mut fitness, &popsize);

//bottleneck
for i in 0..136
{

newpopsize = popsize - 14;

let num_muts = poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0;

mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);

(pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
= reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);

writeln!(outfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, newpopsize, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero,sample_ave_homo, sample_het_1p2q2).expect("Failed to write to outfile");

popsize = newpopsize;

update_fitness (& pop, &selection, &mut fitness, &popsize);
}

let formatted_runs = format!("{}", runs);
file_fit_after_bneck.write_all(b"\t").unwrap();
file_fit_after_bneck.write_all(formatted_runs.as_bytes()).unwrap();
file_fit_after_bneck.write_all(b"\tafter\n").unwrap();

file_freq_after_bneck.write_all(b"\t").unwrap();
file_freq_after_bneck.write_all(formatted_runs.as_bytes()).unwrap();
file_freq_after_bneck.write_all(b"\tafter\n").unwrap();

print_stats(&pop, &selection, &popsize, &genomesize, &sample, &mut file_fit_after_bneck, &mut file_freq_after_bneck); 


}
}


fn selection_distribution(selection: &mut Vec<(f64, f64)>, genomesize: &usize, beta: &f64, skew_normal: &SkewNormal<f64>, normal: &Normal<f64>, gamma: &Gamma<f64>, lethals: &f64) 
{

let non_lethals = (*genomesize as f64 * (1.0 - lethals)).round() as usize;

    for i in 0..non_lethals
    {

//        let mut s = skew_normal.sample(&mut rand::thread_rng());

        let s = -1.0 * gamma.sample(&mut rand::thread_rng());

//        while (s < -1.0 || s > 0.001)
//        {s = skew_normal.sample(&mut rand::thread_rng());}
   
 //    s *= 10.0; // must be different in RNS and SBS

//          if s >  1.0 {println!("s is outside the range: {}", s); s =  1.0;}
//     else if s < -1.0 {println!("s is outside the range: {}", s); s = -1.0;}

     selection[i].0 = s * 10.0; //selection coefficient of 11 - this is the homozygous state
//     selection[i].1 = 1.0 / (1.0 + f64::exp(-beta * selection[i].0)); //norm(mt); sigmoidal

//     selection[i].1 = selection[i].1.max(0.0).min(1.0);

//     selection[i].1 *= selection[i].0; //get h into sh form.

//     selection[i].0 = -0.000;
     selection[i].1 = selection[i].0 * 0.5;
 
    }

   for i in non_lethals..*genomesize
    {
  //   let mut s = normal.sample(&mut rand::thread_rng());

//     while (s < -1.0 || s > 0.000)
//     {s = normal.sample(&mut rand::thread_rng());}

 //    s *= 10.0; // must be different in RNS and SBS

//          if s >  1.0 {println!("s is outside the range: {}", s); s =  1.0;}
//     else if s < -1.0 {println!("s is outside the range: {}", s); s = -1.0;}

     selection[i].0 = -0.1; //selection coefficient of 11 - this is the homozygous state
     selection[i].1 = 0.0; //get h into sh form.
    // selection[i].1 = 1.0 / (1.0 + f64::exp(-beta * selection[i].0)); //norm(mt); sigmoidal

//     selection[i].1 = selection[i].1.max(0.0).min(1.0);

//     selection[i].1 *= selection[i].0; //get h into sh form.

    }
}
 
 
fn print_stats(pop: &Vec<HashMap<u32, u8>>, selection: &Vec<(f64,f64)>, popsize: &usize, genomesize: &usize, sample: &usize, file_fit: &mut File, file_freq: &mut File) 
{

    let mut q:        Vec<f64> = vec![0.0; *genomesize];
    let mut sample_q: Vec<f64> = vec![0.0; *genomesize];

    for i in 0..*popsize 
       {
        if let Some(pop_gene) =  pop.get(i) 
         {

            let mut fitness = 1.0;

            for j in 0..*genomesize
                {
                if let Some(gene) = pop_gene.get(&(j as u32))
                    {

                    if *gene == 2 
                    { fitness *= f64::exp(selection[*gene as usize].0);
                        q[j] += 2.0;

                    if i < *sample
                    {sample_q[j] +=  2.0;}

                    } // if 11 - s 
                    else if *gene == 1 
                    { fitness *= f64::exp(selection[*gene as usize].1);
                        q[j] += 1.0;

                        if i < *sample
                        {sample_q[j] +=  1.0;}

                    } // 01 - sh
                    
                    }                
                } 

                let formatted_fitness = format!("{}", fitness);
                file_fit.write_all(formatted_fitness.as_bytes()).unwrap();
                file_fit.write_all(b"\n").unwrap();
          }
    }

for i in 0..*genomesize
    {
        let formatted_freq = format!("{}", q[i]);
        let formatted_sample_freq = format!("{}", sample_q[i]);
        file_freq.write_all(formatted_freq.as_bytes()).unwrap();
        file_freq.write_all(b"\t").unwrap();
        file_freq.write_all(formatted_sample_freq.as_bytes()).unwrap();
        file_freq.write_all(b"\n").unwrap();
    }

}

fn update_fitness(pop: &Vec<HashMap<u32, u8>>, selection: &Vec<(f64, f64)>, fitness: &mut Vec<f64>, popsize: &usize)
{

 fitness.resize(*popsize, 1.0); // for (_, fit) in fitness.iter_mut().enumerate() {*fit = 1.0;} //makes all fitness elements equal to 1

    let mut maxfit = 0.0;

    for i in 0..*popsize
    {

        fitness[i] = 1.0;

        if let Some(gene_map) = pop.get(i) 
        {
            for (gene, &value) in gene_map.iter() 
            {
                if value == 2 
                { fitness[i] *= f64::exp(selection[*gene as usize].0);} // if 11 - s 
             else if value == 1 
                { fitness[i] *= f64::exp(selection[*gene as usize].1);} // 01 - sh

                if fitness[i] < 0.0 {fitness[i] = 0.0;}
            }

            if fitness[i] > maxfit {maxfit = fitness[i];}
        }
    }

    for i in 0..*popsize
    {
        fitness[i] /= maxfit;
    }
    
}

fn mutation(pop: &mut Vec<HashMap<u32, u8>>, fitness: &mut Vec<f64>, selection: &Vec<(f64, f64)>,
num_muts: &usize, popsize: &usize, genomesize: &usize, rng: &mut impl Rng)
{
       
     let mut mutation_list: Vec<usize> = Vec::new();

       //itterate across all new mutations being added to the population
       for _ in 0..*num_muts  
       {
         //getting the random individuals, sites and nucleotides we want to chage
         let mut_ind       = rng.gen_range(0..*popsize   );  // generate a random individual from the distribution
         let mut_site: u32 = rng.gen_range(0..*genomesize) as u32;  // generate a random site from the distribution

         if !mutation_list.contains(&mut_ind) {mutation_list.push(mut_ind);}
                   
         //change the population adding new mutations
         if let Some(pop_genome) = pop.get_mut(mut_ind)  /*individual has at least one mutant site*/ 
          {
            if let Some(value) = pop_genome.get(&mut_site) 
            {
                if *value == 2     {pop_genome.insert(mut_site, 1);}                   //11 -> 01
                else if *value == 1
                {
                let random = rng.gen_range(1..=2);                                      // generate a random 1 or 2
                    if random == 1 {pop_genome.insert(mut_site, 2);}                   // 01 -> 11
                    else           {pop_genome.remove(&mut_site);}                     // 01 -> 00   
                }                                                                      // 10 -> 11 or 00            
            } /*the site is not 00 - has at least one mutation */
            else                    {pop_genome.insert(mut_site, 1);}                  //00 -> 01
          }
          else                      {pop[mut_ind].insert(mut_site, 1);}                //00  -> 01
        }        


     let mut update_maxfit: bool = false;
     let mut maxfit = 1.0;

        for i in mutation_list.iter() 
        { 
            fitness[*i] = 1.0;

            if let Some(pop_genome) = pop.get(*i)  /*individual has at least one mutant site*/ 
            {
                for (gene, &value) in pop_genome.iter() 
                {
                    if value == 2 
                    { fitness[*i] *= f64::exp(selection[*gene as usize].0);} // if 11 - s 
                 else if value == 1 
                    { fitness[*i] *= f64::exp(selection[*gene as usize].1);} // 01 - sh
    
                    if fitness[*i] < 0.0 {fitness[*i] = 0.0;}
                }                              
            }

            if fitness[*i] > maxfit 
             {
               update_maxfit = true;
               maxfit = fitness[*i];
             }  
        }

        if update_maxfit
        {
            for i in 0..*popsize
            {fitness[i] /= maxfit;}   
        }
}


fn reproduction (pop: &mut Vec<HashMap<u32, u8>>, fitness: &Vec<f64>, popsize: &usize, newpopsize: &usize, genomesize: &usize, sample: &usize, rng: &mut impl Rng)
 -> (Vec<HashMap<u32, u8>>, usize, f64, f64, f64, f64, f64, f64, f64, f64)
{

    let mut newpop: Vec<HashMap<u32, u8>> = Vec::new();
    
    let mut q: Vec<f64> = vec![0.0; *genomesize];
    let mut ave_hetero: f64 = 0.0;
    let mut ave_homo:   f64 = 0.0;
    let mut total_1:    f64 = 0.0;

    let mut sample_q: Vec<f64> = vec![0.0; *genomesize];
    let mut sample_ave_hetero: f64 = 0.0;
    let mut sample_ave_homo:   f64 = 0.0;
    let mut sample_total_1:    f64 = 0.0;

    for i in 0..*newpopsize
    {
        let mut selected: bool = false;
        let mut mother = 0;
        let mut father = 0;

        let mut hetero: usize = 0;
        let mut homo:   usize = 0;

        let mut sample_hetero: usize = 0;
        let mut sample_homo:   usize = 0;

        while !selected
        {
         mother = rng.gen_range(0..*popsize);  //generate a random mother from the population
         let selection = rng.gen::<f64>();

         if fitness[mother] >= selection {selected = true;} 
        }

        selected = false;
        while !selected
        {
         father = rng.gen_range(0..*popsize);  //generate a random mother from the population
         let selection = rng.gen::<f64>();

         if (fitness[father] >= selection) && (mother != father){selected = true;} 
        }

        if let Some(mother_genome) = pop.get(mother) //pass the mutations from mother to child
        {

            let mut child_genome: HashMap<u32, u8> = HashMap::new();

            for (gene, value) in mother_genome.iter() 
            {
                let go = if *value == 1 { rng.gen::<u32>() % 2 == 0 } else { true }; //mother is 10 random, decide if you are 0 or 1

                if go // if value == 2 or you are passing 1 from when value == 1
                {
                    child_genome.insert(*gene, 1);
                    hetero += 1;
                    q[*gene as usize] += 1.0;

                    if i < *sample 
                    {   sample_hetero += 1;
                        sample_q[*gene as usize] += 1.0;
                    }
                } //add one hetero
            }
            newpop.insert(i, child_genome);
        }  

        if let Some(father_genome) = pop.get(father) // pass the mutations from father to child
        {
            if let Some(child_genome) = newpop.get_mut(i) 
            {
                for (gene, value) in father_genome.iter() 
                {
                    let go = if *value == 1 { rng.gen::<u32>() % 2 == 0 } else { true }; //father is 10 random, decide if you are 0 or 1

                    if go // if value == 2 or you are passing 1 from when value == 1
                    {
                        if let Some(_) = child_genome.get(&gene) 
                        {
                            child_genome.insert(*gene, 2);
                            hetero -= 1;
                            homo += 1;
                            q[*gene as usize] += 1.0;

                            if i < *sample 
                            {   sample_hetero -= 1;
                                sample_homo +=1;
                                sample_q[*gene as usize] += 1.0;
                            }        
                        } //the child is already 01, add homo and q, substract hetero
                        else
                        {
                            child_genome.insert(*gene, 1);
                            hetero += 1;

                            q[*gene as usize] += 1.0;
                            if i < *sample 
                            {   sample_hetero += 1;
                                sample_q[*gene as usize] += 1.0;
                            }
                        } //the child is 00, add hetero
                    }
                }
            }
            else
            {
                let mut child_genome: HashMap<u32, u8> = HashMap::new();

                for (gene, value) in father_genome.iter() 
                {
                    let go = if *value == 1 { rng.gen::<u32>() % 2 == 0 } else { true }; //father is 10 random, decide if you are 0 or 1

                    if go // if value == 2 or you are passing 1 from when value == 1
                        {
                            child_genome.insert(*gene, 1);
                            hetero += 1;
                            q[*gene as usize] += 1.0;

                            if i < *sample 
                            {   sample_hetero += 1;
                                sample_q[*gene as usize] += 1.0;
                            }
                        } //the child is 00, add hetero
                }
                newpop.insert(i, child_genome);
            }
        }  

         ave_hetero += hetero as f64;
         ave_homo   += homo as f64;

         sample_ave_hetero += sample_hetero as f64;
         sample_ave_homo   += sample_homo as f64;
        }

        let mut fixed: usize = 0;

        let mut het_1p2q2: f64 = 0.0;
        let mut sample_het_1p2q2: f64 = 0.0;

//IF FIXED for population and sample
for i in 0..*genomesize
{
    if q[i] as usize == *newpopsize * 2 // 11 was fixed, so we remove this from the population, keeping the distribution of new mutations the same
    {
        fixed += 1;
        q[i] = 0.0;
        ave_homo -= *newpopsize as f64;
        for j in 0..*newpopsize
        {newpop[j].remove(&(i as u32));}
    }
    else
    {
     q[i] /= (*newpopsize as f64) * 2.0;

     het_1p2q2 += 1.0 - q[i].powi(2) - (1.0 - q[i]).powi(2);
    }

    if sample_q[i] as usize == *sample * 2 // 11 was fixed in the sample
      {
        sample_q[i] = 0.0;
        sample_ave_homo -= *sample as f64;
      }
     else
      {
        sample_q[i] /= (*sample as f64) * 2.0;
        sample_het_1p2q2 += 1.0 - sample_q[i].powi(2) - (1.0 - sample_q[i]).powi(2);
      } 
}


het_1p2q2 /= *genomesize as f64;
sample_het_1p2q2 /= *genomesize as f64;

total_1 = (ave_hetero + ave_homo * 2.0) / (*newpopsize as f64);
sample_total_1 = (sample_ave_hetero + sample_ave_homo * 2.0) / (*sample as f64);

ave_hetero /= *newpopsize as f64 * *genomesize as f64;
ave_homo /= *newpopsize as f64 * *genomesize as f64;
sample_ave_hetero /= *sample as f64 * *genomesize as f64;
sample_ave_homo /= *sample as f64 * *genomesize as f64;

return (newpop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2);
}


/*
//after running this function maxfit == 1
fn normalize_fitness(fitness: &mut Vec<f64>, popsize: &usize) 
{
    let mut maxfit = -1000.0;  //start with very low fitness

    for i in 0..*popsize
    {
    let fit = fitness[i];
    if fit > maxfit {maxfit = fit;}
    }

    for i in 0..*popsize
    { fitness[i] /= maxfit; }
}

*/

//functions for debugging or printout

fn print_time()
{
            // Get the current local time
            let local_time = Local::now();

            // Extract components (year, month, day, hour, minute, second)
            let year = local_time.year();
            let month = local_time.month();
            let day = local_time.day();
            let hour = local_time.hour();
            let minute = local_time.minute();
            let second = local_time.second();
              
            println!("{:02}:{:02}:{:02}",hour, minute, second);
        
}

/*
fn print_fitness(fitness: &Vec<f64>)
 {
    for &num in fitness {
        println!("{}", num);
    }
}


fn print_selection(selection: &Vec<(f64,f64)>)
 {
    for &(x, y) in selection {
        println!("{}\t{}", x, y);
    }
}


fn print_pop(pop: &Vec<HashMap<u32, u8>>, popsize: &usize, genomesize: &u32) 
{
        for i in 0..*popsize {
        print!("Individual: {}\t", i);
       if let Some(pop_gene) =  pop.get(i) 
       {
        for j in 0..*genomesize
            {
            if let Some(gene) = pop_gene.get(&j) {print!("{}", gene);}
               else{print!("{}", 0);}
             } 
          }
         else  
         {
            for j in 0..*genomesize
            {
            print!("{}", 0);
            }
        }
        println!("");
        }
}
*/