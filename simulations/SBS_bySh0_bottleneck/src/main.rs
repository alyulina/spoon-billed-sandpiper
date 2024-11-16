//use rand::prelude::*;
use rand::Rng;
use rand_distr::{Poisson, Distribution, SkewNormal, Normal};
use std::collections::HashMap;
use std::process;
use std::fs::File;
use std::io::prelude::*;

use chrono::{Local, Datelike, Timelike};

fn main() {

for runs in 0..1
{

print_time();

println!("s\ttotal_1\thet2pq\thomo\tsample_total_1\tsample_het2pq\tsample_homo\ttotal_1\thet2pq\thomo\tsample_total_1\tsample_het2pq\tsample_homo");

let nastya_s = vec![
/*0.00010000000000000,
0.00008799225435691,
0.00007742636826811,
0.00006812920690580,
0.00005994842503189,
0.00005274997063703,
0.00004641588833613,
0.00004084238652675,
0.00003593813663805,
0.00003162277660168,
0.00002782559402207,
0.00002448436746822,
0.00002154434690032,
0.00001895735652406,
0.00001668100537200,
0.00001467799267622,
0.00001291549665015,
0.00001136463666386,
0.00001000000000000,*/
0.00000000000000000,
-0.00001000000000000,
-0.00001138098301710,
-0.00001295267744355,
-0.00001474142020111,
-0.00001677718529567,
-0.00001909408609248,
-0.00002173094695456,
-0.00002473195382353,
-0.00002814739464454,
-0.00003203450204251,
-0.00003645841237071,
-0.00004149325720215,
-0.00004722340555418,
-0.00005374487766218,
-0.00006116695399294,
-0.00006961400646015,
-0.00007922758252753,
-0.00009016877712317,
-0.00010262093211115,
-0.00011679270855561,
-0.00013292158325925,
-0.00015127782816797,
-0.00017216903932435,
-0.00019594529126210,
-0.00022300500321348,
-0.00025380161543011,
-0.00028885118749228,
-0.00032874104593189,
-0.00037413962607748,
-0.00042580767304122,
-0.00048461098954332,
-0.00055153494418929,
-0.00062770098331558,
-0.00071438542309320,
-0.00081304083678879,
-0.00092532039557025,
-0.00105310557073619,
-0.00119853766157626,
-0.00136405367717549,
-0.00155242717343479,
-0.00176681472961469,
-0.00201080884321080,
-0.00228849812952179,
-0.00260453583467540,
-0.00296421781018708,
-0.00337357125567263,
-0.00383945571677891,
-0.00436967803075706,
-0.00497312314582438,
-0.00565990300645773,
-0.00644152599949322,
-0.00733108980044445,
-0.00834350085156978,
-0.00949572414948803,
-0.01080706752803955,
-0.01229950520012777,
-0.01399804598013949,
-0.01593115235725613,
-0.01813121744207746,
-0.02063510777876421,
-0.02348478111861564,
-0.02672798950712899,
-0.03041907946618789,
-0.03461990268005186,
-0.03940085244553465,
-0.04484204325419174,
-0.05103465327280476,
-0.05808245221814101,
-0.06610354022862211,
-0.07523232687121760,
-0.08562178344582902,
-0.09744600632908472,
-0.11090313431155939,
-0.12621866881430852,
-0.14364925262166944,
-0.16348697045064187,
-0.18606424342159841,
-0.21175939944708896,
-0.24100301288187434,
-0.27428511966786967,
-0.31216342887834397,
-0.35527266826243264,
-0.40433522039347730,
-0.46017322765137969,
-0.52372236888247359,
-0.59604753859271320,
-0.67836069141083022,
-0.77204115084153024,
-0.87865872262302647,
-1.00000000000000000];

for step in &nastya_s 
{ 
    
let mut popsize: usize    = 8000;
let mut genomesize: usize = 300000;
let mut mut_rate = 0.00000012;
let mut newpopsize: usize = 0;
let sample = 7;

let expect_muts = mut_rate * genomesize as f64 * 2.0 * 10000.0; //this is the lambda parameter in Poisson

//DO NOT FORGET ABOUT POPSIZIE - AS POPSIZE CHANGES, SO MUST THE EXPECTED NUMBER OF MUTATIONS.

//this tuple carries the variables needed to define s and sh, in this order: my, sigma, alpha and beta.
// beta is used for the exponential to define sh, the rest are for the skewed normal distribution

//mu=-0.0000546277, sigma=0.0000100000, alpha=0.0, beta=3000
let mut skew: (f64, f64, f64, f64) = (-0.0027, 0.0010, -0.0, 3000.0/1.0);

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
//let skew_normal = SkewNormal::new(skew.0, skew.1, skew.2).unwrap();
//let normal      =     Normal::new(-0.05, 0.1).unwrap();
//let normal      =     Normal::new(-0.0008, 0.0005).unwrap();
let poisson = Poisson::new(expect_muts).unwrap();
let mut rng = rand::thread_rng();

/*define the population of popsize individuals but when the length of the genome u16 - this is 65536 sites with bool alleles (empty,0,1) */
let mut pop: Vec<HashMap<u32, u8>> = Vec::with_capacity(popsize);
pop.resize(popsize, HashMap::new());

/*for each individual in the population this array keeps track of its fitness */
let mut fitness: Vec<f64> = Vec::with_capacity(popsize);

/*for each site in the genome this array keeps track of s and sh - this is used to define fitness*/
let mut selection: Vec<(f64, f64)> = vec![(0.0, 0.0); genomesize];


//let s = -0.1 + (step as f64 * 0.001);
let s = step * 10.0;
let h = s * 0.0;

/*for each site in the genome this array keeps track of s and sh - this is used to define fitness*/
let mut selection: Vec<(f64, f64)> = vec![(s, h); genomesize];

popsize = 2000;
fitness.resize(popsize, 1.0);

newpopsize = popsize;

//burnin
for i in 0..30000
{

    let num_muts = (poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0);

    mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);
    
    (pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
    = reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);
    
    update_fitness (& pop, &selection, &mut fitness, &popsize);
}

print!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t", s, total_1, ave_hetero, ave_homo, sample_total_1, sample_ave_hetero, sample_ave_homo);

print_stats(&pop, &selection, &popsize, &genomesize, &sample); 

/*
//growth
for i in 0..1000
{
    newpopsize = popsize + 6;

    let num_muts = poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0;
    
    mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);
    
    (pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
    = reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);
       
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
    
    update_fitness (& pop, &selection, &mut fitness, &popsize);    
}
*/

//decline
for i in 0..136 //565
{

    newpopsize = popsize - 14;

    let num_muts = poisson.sample(&mut rand::thread_rng()) * popsize as f64/10000.0;
    
    mutation (&mut pop, &mut fitness, &selection, &(num_muts as usize), &popsize, &genomesize, &mut rng);
    
    (pop, fixed, total_1, ave_hetero, ave_homo, het_1p2q2, sample_total_1, sample_ave_hetero, sample_ave_homo, sample_het_1p2q2)
    = reproduction (&mut pop, &fitness, &popsize, &newpopsize, &genomesize, &sample, &mut rng);

    popsize = newpopsize;
    
    update_fitness (& pop, &selection, &mut fitness, &popsize);
}

println!("{}\t{}\t{}\t{}\t{}\t{}", total_1, ave_hetero, ave_homo, sample_total_1, sample_ave_hetero, sample_ave_homo);

print_stats(&pop, &selection, &popsize, &genomesize, &sample); 

}
}
}

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
        
            println!("Hello, world!");
        
            println!(
                "{:02}:{:02}:{:02}",
                hour, minute, second
            );
        
}

fn print_stats(pop: &Vec<HashMap<u32, u8>>, selection: &Vec<(f64,f64)>, popsize: &usize, genomesize: &usize, sample: &usize) 
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

                    if (i < *sample)
                    {sample_q[j] +=  2.0;}

                    } // if 11 - s 
                    else if *gene == 1 
                    { fitness *= f64::exp(selection[*gene as usize].1);
                        q[j] += 1.0;

                        if (i < *sample)
                        {sample_q[j] +=  1.0;}

                    } // 01 - sh
                    
                    }                
                } 
          }
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
