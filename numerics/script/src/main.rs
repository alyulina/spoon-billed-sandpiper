//! Compute fixed parameters or run a CLI to ask for novel parameters.

const N: f64 = 2e4; // 2e3; //  
const U: f64 = 1.2e-8; // 1.2e-7; // 
const MAX_ITERS: usize = 1000_0000;
const TOL: f64 = 1e-10;

fn main() {
    // Toogle to run parameters
    if false {
        // Run parameters
        use rayon::prelude::*;
        println!("Running parameters");
        // Computing
        let results = PARAMETERS.into_par_iter().map(|(p, shape, scale)| {
            let alpha = shape;
            let beta = 1. / scale;
            expectation::mixture(p, alpha, beta)
        }).collect::<Vec<_>>();
        // Reporting
        for ((p, shape, scale), result) in PARAMETERS.iter().zip(results.iter()) {
            println!("{p}, {shape}, {scale}, {result:.7}");
        }
    }
    // Toogle to precompute values
    if false {
        use rayon::prelude::*;
        println!("Precomputing values");
        let h = 0.5;
        // Computing
        let results = VALUES.into_par_iter().map(|(s, _)| expectation::conditional(s, h) ).collect::<Vec<_>>();
        // Reporting
        for ((s, _), result) in VALUES.iter().zip(results.iter()) {
            println!("({s}, {result:.8}),");
        }
    }
    // Toogle to compute plot values for LaTeX report
    if false {
        use rayon::prelude::*;
        println!("Computing plot values");
        
        let h = 0.0;
        println!("\nh = {h}\n");
        // Computing
        let results = PLOT_S.into_par_iter().map(|s| expectation::conditional(s, h) ).collect::<Vec<_>>();
        // Reporting
        for (s, result) in PLOT_S.iter().zip(results.iter()) {
            println!("{s} {result:.8}");
        }

        let h = 0.5;
        println!("\nh = {h}\n");
        // Computing
        let results = PLOT_S.into_par_iter().map(|s| expectation::conditional(s, h) ).collect::<Vec<_>>();
        // Reporting
        for (s, result) in PLOT_S.iter().zip(results.iter()) {
            println!("{s} {result:.8}");
        }

        println!("\nVarying alpha/beta\n");
        let alpha = 1.0;
        // Computing
        let results = PLOT_GAMMA.into_par_iter().map(|beta_inv| expectation::gamma(alpha, 1./beta_inv) ).collect::<Vec<_>>();
        // Reporting
        for (beta_inv, result) in PLOT_GAMMA.iter().zip(results.iter()) {
            println!("{beta_inv} {result:.8}");
        }
    }
    // Toogle to run CLI
    if true {
        // Cli
        cli::run_cli();
    }
}

mod cli {
    use clap::{Parser, Subcommand};

    /// Approximate, up to five digits, the expected polymorphisms.
    #[derive(Parser)]
    #[command(version, about, long_about = None)]
    #[command(propagate_version = true)]
    struct Cli {
        #[command(subcommand)]
        command: Commands,
    }

    #[derive(Subcommand)]
    enum Commands {
        /// The selection is a mixture of a point mass and a gamma distribution
        Mixture {
            /// Mixture probability for the negative selection
            p: f64,
            /// Shape of the gamma distribution
            shape: f64,
            /// Scale of the gamma distribution (parameter beta in Wikipedia)
            scale: f64,
        },
        /// The selection and dominance are given directly
        Precise {
            /// Selection coefficient
            s: f64,
            /// Dominand coefficient
            h: f64,
        },
    }

    pub fn run_cli() {
        let cli = Cli::parse();
        let n = crate::N;
        let u = crate::U;

        match &cli.command {
            Commands::Mixture { p, shape, scale } => {
                // Feedback
                println!("Parameters ");
                println!("N = {n}\nU = {u}\np = {p}\n𝛼 = {shape}\n𝛽 = {scale}");
                // Computing
                let result = crate::expectation::mixture(*p, *shape, *scale);
                // Reporting
                println!("Result");
                println!("{result:.10}");
            }
            Commands::Precise { s, h } => {
                // Feedback
                println!("Parameters ");
                println!("N = {n}\nU = {u}\ns = {s}\nh = {h}");
                // Computing
                let result = crate::expectation::conditional(*s, *h);
            	let upper_bound = crate::expectation::s_derivative_conditional(*s, *h); // This is an upper bound because neutrality (s = 0) should give the largest expecte polymorphisms. Therefore, we think the expected polymorphisms are increasing for negative selections.
            	if upper_bound < result {
            		println!("Careful! Unexpected result found accoring to population genetics");
            	}
                // Reporting
                println!("Result");
                println!("{:.10}", crate::expectation::conditional_old(*s, *h));
                /*
                println!("Debugging");
                println!("{:.10}", result);
                println!("Upper bound = {:.10}", upper_bound);
                println!("Upper bound = {:.10}", crate::expectation::s_derivative_conditional_old(*s, *h));
                println!("C(s, h) = {:.10}", crate::function::normalizing_constant(*s, *h));
                println!("C(s, h) = {:.10}", crate::function::normalizing_constant_old(*s, *h));
                */
            }
        }
    }
}

///
mod expectation {

    use accurate::traits::*;
    use statrs::distribution::ContinuousCDF;

    /// Expected polymorphisms of a gamma distributed selection and a fixed deleterious selection.
    pub fn mixture(p: f64, alpha: f64, beta: f64) -> f64 {
        let value_at_mass = 0.00008501514709555435; // Value of s = -1e-2, h = 0.0
        p * value_at_mass + (1. - p) * gamma(alpha, beta)
    }

    /// Expected polymorphisms when the selection `s` is such that -s follows a Gamma distribution with shape alpha and scale beta.
    pub fn gamma(alpha: f64, beta: f64) -> f64 {
        let expected_poly = crate::init_interpolation();

        assert!(
            alpha <= beta * beta * 1e-7,
            "Parameters of Gamma distribution have too much mass outside [0, 0.1]"
        );

        let gamma = statrs::distribution::Gamma::new(alpha, beta)
            .expect("Gamma distribution is not defined for these parameters");

        // Find a step that suitable for expected_poly and gamma
        let mean = alpha / beta;
        let mut s_min = mean;
        while gamma.cdf(s_min) > 1e-7 {
            s_min /= 1.1;
        }
        let mut s = s_min; // this is the start
        let delta = s_min.min(1e-8); // this is the step

        let mut total = accurate::sum::OnlineExactSum::zero();
        let mut probability = gamma.cdf(s + delta) - gamma.cdf(s);
        let mut value = expected_poly
            .clamped_sample(-s)
            .expect("This is a bug. There were no enough percomputed points. Please report");

        // while the integral has some value remaining
        while (1. - gamma.cdf(s)) * value > 1e-8 {
            total += probability * value;
            s += delta;
            probability = gamma.cdf(s + delta) - gamma.cdf(s);
            value = expected_poly
                .clamped_sample(-s)
                .expect("This is a bug. There were no enough percomputed points. Please report");
        }

        total.sum()
    }

    pub fn conditional(s: f64, h: f64) -> f64 {
        let mut total = 0.;
        let n = 100_000;
        let n_inv = 1. / n as f64;
        // [0, 1]
        //
        let grid = (0..n).map(|i| i as f64 * n_inv);
        total += grid.map(|x| {
            crate::function::conditional_integrand(s, h, x)
                }).map(|v| v * n_inv).sum::<f64>();

        total /= crate::function::normalizing_constant(s, h);

        total
    }
    
    /// Conditional expectation given both s and h.
    pub fn conditional_old(s: f64, h: f64) -> f64 {
        assert!(s <= 0.0);
        assert!(s >= - crate::N * 1e-6, "be careful!! population size and selection are not in a controlled range.");

        let eps: f64 = (1e-20 / s.abs().powi(2)).min(1e-5);
        let n: i32 = {
            let mut n = 1;
            let mut x = eps;
            while x < 0.5 {
                n += 1;
                x = (eps * x).sqrt() + x;
            }
            n
        };
        let mut total: accurate::sum::OnlineExactSum<f64> = accurate::sum::OnlineExactSum::zero();

        // Interval [0, 1/2]
        //
        // grid  = { eps, sqrt(eps * eps) + eps, ..., 0.5 }
        //
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                x
            }
        };
        let grid = (2..=n).map(|_| generator());
        //
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                x
            }
        };
        let mut values = (2..=n)
            .map(|_| generator())
            .map(|x| crate::function::conditional_integrand(s, h, x));
        //
        let mut prev = eps;
        for next in grid {
            let value = values.next().expect("this is a bug! please report");
            assert!(next > prev);
            total += (next - prev) * value; // Never fails since the iterators are of the same length
            prev = next;
        }

        // Interval [1/2, 1]
        //
        // We generate a grid from 1 to 0.5
        // grid = {1 - eps, 1 - sqrt(eps * eps) - eps, ..., 0.5}
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                1.0 - x
            }
        };
        let grid = (2..=n)
            .map(|_| generator());
        //
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                1.0 - x
            }
        };
        let mut values = (2..=n)
            .map(|_| generator())
            .map(|x| crate::function::conditional_integrand(s, h, x));
        //
        let mut prev = 1.0 - eps;
        for next in grid {
            let value = values.next().expect("this is a bug! please report");
            assert!(next < prev);
            total += (next - prev).abs() * value; // Never fails since the iterators are of the same length
            prev = next;
        }

        //
        total.sum() / crate::function::normalizing_constant(s, h)
    }

    /// Conditional expectation given both s and h.
    pub fn s_derivative_conditional(s: f64, h: f64) -> f64 {
        assert!(s <= 0.0);
        assert!(s >= - crate::N * 1e-6, "be careful!! population size and selection are not in a controlled range.");

        let mut total = 0.;
        let n = 10_000_000;
        let n_inv = 1. / n as f64;
        // [0, 1]
        //
        let grid = (1..n).map(|i| i as f64 * n_inv);
        total += grid.map(|x| {
            crate::function::s_derivative_conditional_integrand(s, h, x)
                }).map(|v| v * n_inv).sum::<f64>();

        total /= crate::function::s_derivative_normalizing_constant(s, h);

        total
    }

    /// Conditional expectation given both s and h.
    pub fn s_derivative_conditional_old(s: f64, h: f64) -> f64 {
        assert!(s <= 0.0);
        assert!(s >= - crate::N * 1e-6, "be careful!! population size and selection are not in a controlled range.");

        let eps: f64 = (1e-20 / s.abs().powi(2)).min(1e-5);
        let n: i32 = {
            let mut n = 1;
            let mut x = eps;
            while x < 0.5 {
                n += 1;
                x = (eps * x).sqrt() + x;
            }
            n
        };
        let mut total: accurate::sum::OnlineExactSum<f64> = accurate::sum::OnlineExactSum::zero();

        // Interval [0, 1/2]
        //
        // grid  = { eps, sqrt(eps * eps) + eps, ..., 0.5 }
        //
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                x
            }
        };
        let grid = (2..=n).map(|_| generator());
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                x
            }
        };
        let mut values = (2..=n)
            .map(|_| generator())
            .map(|x| crate::function::s_derivative_conditional_integrand(s, h, x));

        let mut prev = eps;
        for next in grid {
            let value = values.next().expect("this is a bug! please report");
            // dbg!(next, value);
            total += (next - prev) * value; // Never fails since the iterators are of the same length
            prev = next;
        }

        // Interval [1/2, 1]
        //
        // We generate a grid from 1 to 0.5
        // grid = {1 - eps, 1 - sqrt(eps * eps) - eps, ..., 0.5}
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                1.0 - x
            }
        };
        let grid = (2..=n)
            .map(|_| generator());
        let mut generator = {
            let mut x = eps;
            move || {
                x = ((eps * x).sqrt() + x).min(0.5);
                1.0 - x
            }
        };
        let mut values = (2..=n)
            .map(|_| generator())
            .map(|x| crate::function::s_derivative_conditional_integrand(s, h, x));

        let mut prev = 1.0 - eps;
        for next in grid {
            let value = values.next().expect("this is a bug! please report");
            // dbg!(next, value);
            total += (next - prev).abs() * value; // Never fails since the iterators are of the same length
            prev = next;
        }

        total.sum() / crate::function::s_derivative_normalizing_constant_old(s, h)
    }

}

/// Functions used for the numerical integration.
mod function {

    use crate::{MAX_ITERS, N, TOL, U};
    use gkquad::{single::Integrator, Tolerance};

    /// Conditional integrand.
    ///
    /// # Bounds
    ///
    /// Bounds in `[0, 1]`.
    /// Evaluation `1`
    pub fn conditional_integrand(s: f64, h: f64, x: f64) -> f64 {
        2. * (4. * N * U * (x * (1. - x)).ln() + 2. * N * s * x * (x + 2. * h * (1. - x))).exp()
    }

    pub fn s_derivative_conditional_integrand(s: f64, h: f64, x: f64) -> f64 {
        x * (x + 2. * h * (1. - x)) * 2. * (4. * N * U * (x * (1. - x)).ln() + 2. * N * s * x * (x + 2. * h * (1. - x))).exp()
    }

    /// Computes the normalizing constant for the distribution in `[0, 1]` propostional to
    /// `exp( 2 N s x (x + 2 h (1 - x)) )  x^{4 N U - 1}  (1 - x)^{4 N U - 1}`.
    /// In other words, it computes
    /// `int_{[0, 1]} exp( 2 N s x (x + 2 h (1 - x)) )  x^{4 N U - 1}  (1 - x)^{4 N U - 1}`.
    pub fn normalizing_constant(s: f64, h: f64) -> f64 {
        let mut total = 0.;
        // [0, 1/2]
        //
        let n = 10_000_000;
        let n_inv = 1. / n as f64;
        let grid = (1..n).map(|i| i as f64 * n_inv * 0.5);
        total += grid.map(|x| {
            normalizing_constant_integrand(s, h, x) 
                    - x.powf(4. * N * U - 1.)
                }).map(|v| v * n_inv * 0.5).sum::<f64>();

        // [1/2, 1]
        //
        let grid = (0..(n - 1)).map(|i| i as f64 * n_inv * 0.5 + 0.5);
        total += grid.map(|x| {
            normalizing_constant_integrand(s, h, x)
                    - (2. * N * s).exp() * (1. - x).powf(4. * N * U - 1.)
                }).map(|v| v * n_inv * 0.5).sum::<f64>();

        // Renormalize
        //
        total += (1. + (2. * N * s).exp()) * (0.5_f64).powf(4. * N * U) / (4. * N * U);

        // Report
        total
    }
    pub fn normalizing_constant_old(s: f64, h: f64) -> f64 {
        Integrator::new(|x| {
                normalizing_constant_integrand(s, h, x) 
                    - x.powf(4. * N * U - 1.)
            })
            .max_iters(MAX_ITERS)
            .tolerance(Tolerance::AbsAndRel(TOL, TOL))
            .run((0.0)..0.5)
            .estimate()
            .expect("This is a bug. Computing the normalizing constant failed for the desired accuracy. Please report")
            + Integrator::new(|x| {
                normalizing_constant_integrand(s, h, x)
                    - (2. * N * s).exp() * (1. - x).powf(4. * N * U - 1.)
            })
            .max_iters(MAX_ITERS)
            .tolerance(Tolerance::AbsAndRel(TOL, TOL))
            .run((0.5)..1.0)
            .estimate()
            .expect("This is a bug. Computing the normalizing constant failed for the desired accuracy. Please report")
            + (1. + (2. * N * s).exp()) * (0.5_f64).powf(4. * N * U) / (4. * N * U)
    }

    fn normalizing_constant_integrand(s: f64, h: f64, x: f64) -> f64 {
        ( ( 4. * N * U - 1. ) * (x * (1. - x)).ln() + 2.0 * N * s * x * (x + 2. * h * (1. - x))).exp()
    }

    /// Computes the derivative on `s` of the normalizing constant (divided by 2 N)
    /// for the distribution in `[0, 1]` propostional to
    /// `exp( 2 N s x (x + 2 h (1 - x)) )  x^{4 N U - 1}  (1 - x)^{4 N U - 1}`.
    pub fn s_derivative_normalizing_constant(s: f64, h: f64) -> f64 {
        let mut total = 0.;
        let n = 10_000_000;
        let n_inv = 1. / n as f64;
        // [0, 1/2]
        //
        let grid = (1..n).map(|i| i as f64 * n_inv * 0.5);
        total += grid.map(|x| {
            s_derivative_normalizing_constant_integrand(s, h, x) 
                }).map(|v| v * n_inv * 0.5).sum::<f64>();

        // [1/2, 1]
        //
        let grid = (0..(n - 1)).map(|i| i as f64 * n_inv * 0.5 + 0.5);
        total += grid.map(|x| {
            s_derivative_normalizing_constant_integrand(s, h, x)
                    - (2. * N * s).exp() * (1. - x).powf(4. * N * U - 1.)
                }).map(|v| v * n_inv * 0.5).sum::<f64>();

        // Renormalize
        //
        total += (2. * N * s).exp() * (0.5_f64).powf(4. * N * U) / (4. * N * U);

        // Report
        total
    }
    pub fn s_derivative_normalizing_constant_old(s: f64, h: f64) -> f64 {
        Integrator::new(|x| {
                s_derivative_normalizing_constant_integrand(s, h, x) 
            })
            .max_iters(MAX_ITERS)
            .tolerance(Tolerance::AbsAndRel(TOL, TOL))
            .run((0.0)..0.5)
            .estimate()
            .expect("This is a bug. Computing the normalizing constant failed for the desired accuracy. Please report")
            + Integrator::new(|x| {
                s_derivative_normalizing_constant_integrand(s, h, x)
                    - (2. * N * s).exp() * (1. - x).powf(4. * N * U - 1.)
            })
            .max_iters(MAX_ITERS)
            .tolerance(Tolerance::AbsAndRel(TOL, TOL))
            .run((0.5)..1.0)
            .estimate()
            .expect("This is a bug. Computing the normalizing constant failed for the desired accuracy. Please report")
            + (2. * N * s).exp() * (0.5_f64).powf(4. * N * U) / (4. * N * U)
    }


    fn s_derivative_normalizing_constant_integrand(s: f64, h: f64, x: f64) -> f64 {
        x * (x + 2. * h * (1. - x)) * ( ( 4. * N * U - 1. ) * (x * (1. - x)).ln() + 2.0 * N * s * x * (x + 2. * h * (1. - x))).exp()
    }

}

/// Initialize interpolation of previously computed values.
fn init_interpolation() -> splines::Spline<f64, f64> {
    splines::Spline::from_vec(Vec::from_iter(
        VALUES
            .into_iter()
            .map(|(x, f)| splines::Key::new(x, f, splines::Interpolation::Linear)),
    ))
}

/// Previously computed values `(s, E(poly | s))` where
/// `s` is the selection, and
/// `E(poly | s)` is the expected polymorphisms given s and dominance `h` is 0.5.
const VALUES: [(f64, f64); 102] = [
(-0.1, 0.00000048),
(-0.09, 0.00000053),
(-0.08, 0.00000060),
(-0.07, 0.00000069),
(-0.06, 0.00000080),
(-0.05, 0.00000096),
(-0.04, 0.00000120),
(-0.03, 0.00000160),
(-0.02, 0.00000240),
(-0.01, 0.00000480),
(-0.0099, 0.00000485),
(-0.0098, 0.00000490),
(-0.0097, 0.00000495),
(-0.0096, 0.00000500),
(-0.0095, 0.00000505),
(-0.0094, 0.00000511),
(-0.0093, 0.00000516),
(-0.0092, 0.00000522),
(-0.0091, 0.00000527),
(-0.009, 0.00000533),
(-0.008, 0.00000600),
(-0.007, 0.00000686),
(-0.006, 0.00000800),
(-0.005, 0.00000960),
(-0.004, 0.00001200),
(-0.003, 0.00001600),
(-0.003, 0.00001600),
(-0.0029, 0.00001655),
(-0.0028, 0.00001714),
(-0.0027, 0.00001778),
(-0.0026, 0.00001846),
(-0.0025, 0.00001920),
(-0.0024, 0.00002000),
(-0.0023, 0.00002087),
(-0.0022, 0.00002182),
(-0.0021, 0.00002286),
(-0.002, 0.00002400),
(-0.0019, 0.00002526),
(-0.0018, 0.00002667),
(-0.0017, 0.00002823),
(-0.0016, 0.00003000),
(-0.0015, 0.00003200),
(-0.0014, 0.00003428),
(-0.0013, 0.00003692),
(-0.0012, 0.00004000),
(-0.0011, 0.00004363),
(-0.001, 0.00004800),
(-0.0009, 0.00005333),
(-0.0008, 0.00006000),
(-0.0007, 0.00006857),
(-0.0006, 0.00007999),
(-0.0005, 0.00009599),
(-0.0004, 0.00011998),
(-0.0003, 0.00015997),
(-0.0002, 0.00023977),
(-0.0001, 0.00046241),
(-0.00009, 0.00050457),
(-0.00008, 0.00055251),
(-0.00007, 0.00060648),
(-0.00006, 0.00066615),
(-0.00005, 0.00073016),
(-0.00004, 0.00079566),
(-0.00003, 0.00085786),
(-0.00002, 0.00091025),
(-0.00001, 0.00094562),
(-0.000009, 0.00094797),
(-0.000008, 0.00095009),
(-0.000007, 0.00095196),
(-0.000006, 0.00095360),
(-0.000005, 0.00095499),
(-0.000004, 0.00095613),
(-0.000003, 0.00095702),
(-0.000002, 0.00095765),
(-0.000001, 0.00095803),
(-0.0000009, 0.00095806),
(-0.0000008, 0.00095808),
(-0.0000007, 0.00095810),
(-0.0000006, 0.00095811),
(-0.0000005, 0.00095813),
(-0.0000004, 0.00095814),
(-0.0000003, 0.00095815),
(-0.0000002, 0.00095816),
(-0.0000001, 0.00095816),
(-0.00000009, 0.00095816),
(-0.00000008, 0.00095816),
(-0.00000007, 0.00095816),
(-0.00000006, 0.00095816),
(-0.00000005, 0.00095815),
(-0.00000004, 0.00095815),
(-0.00000003, 0.00095814),
(-0.00000002, 0.00095814),
(-0.00000001, 0.00095814),
(-0.000000009, 0.00095814),
(-0.000000008, 0.00095814),
(-0.000000007, 0.00095814),
(-0.000000006, 0.00095814),
(-0.000000005, 0.00095814),
(-0.000000004, 0.00095814),
(-0.000000003, 0.00095814),
(-0.000000002, 0.00095814),
(-0.000000001, 0.00095814),
(0.0, 0.00095814),
];

/// Parameters (p, shape, scale) for mixture distributed selection
const PARAMETERS: [(f64, f64, f64); 134] = [
    (0., 300.0000000000001, 1.541467657920456e-7),
    (0., 300.0000000000001, 1.541467657920456e-7),
    (0., 290.8506518295917, 1.584893192461114e-7),
    (0., 225.56161320737237, 2.059863942616875e-7),
    (0., 182.5276007813775, 2.5118864315095823e-7),
    (0., 169.59347117570735, 2.708036442160055e-7),
    (0., 127.51258982610176, 3.652004813843986e-7),
    (0., 115.93811031632079, 3.9810717055349687e-7),
    (0., 95.87315155141827, 4.836048047249234e-7),
    (0., 72.81762821148743, 6.30957344480193e-7),
    (0., 72.08434242404265, 6.376784776269807e-7),
    (0., 54.19820188053225, 8.62049856020774e-7),
    (0., 46.27986268057486, 9.999999999999997e-7),
    (0., 40.750112830372316, 1.1416356918193267e-6),
    (0., 30.638870628004046, 1.5174645544368143e-6),
    (0., 29.267338838171042, 1.584893192461114e-6),
    (0., 23.036510285681892, 2.0426117609592174e-6),
    (0., 18.58524958360993, 2.5118864315095823e-6),
    (0., 17.320508075688775, 2.7076390188324098e-6),
    (0., 13.022805810412262, 3.6300647024173715e-6),
    (0., 11.84596516073017, 3.981071705534969e-6),
    (0., 9.791483623609768, 4.887475719264104e-6),
    (0., 7.549938308749734, 6.30957344480193e-6),
    (0., 7.3619428061166206, 6.487995060251171e-6),
    (0., 5.535238985626912, 8.80608239687674e-6),
    (0., 4.874872287742436, 1e-5),
    (0., 4.161791450287818, 1.1941732589246017e-5),
    (0., 3.151143795792919, 1.5848931924611138e-5),
    (0., 3.1291346445318977, 1.5980727030317896e-5),
    (0., 2.3527088612123075, 2.2084851528262447e-5),
    (0., 2.0796912080337586, 2.5118864315095822e-5),
    (0., 1.7689360204744258, 3.0532846214531185e-5),
    (0., 1.386785662694208, 3.9810717055349695e-5),
    (0., 1.3300135414628025, 4.2040057592132016e-5),
    (0., 1.0, 5.9268778545316244e-5),
    (0.1, 300.0000000000001, 1.2221738250849407e-7),
    (0.1, 300.0000000000001, 1.2221738250849407e-7),
    (0.1, 225.56161320737237, 1.5745614479235895e-7),
    (0.1, 224.1010946598654, 1.584893192461114e-7),
    (0.1, 169.59347117570732, 2.144465552856113e-7),
    (0.1, 143.34770274920493, 2.5118864315095823e-7),
    (0.1, 127.51258982610176, 2.8579681663925166e-7),
    (0.1, 95.87315155141827, 3.7423151438881107e-7),
    (0.1, 90.0283006949983, 3.981071705534969e-7),
    (0.1, 72.08434242404266, 5.083976144022012e-7),
    (0.1, 57.05958587182932, 6.30957344480193e-7),
    (0.1, 54.19820188053225, 6.666962666444516e-7),
    (0.1, 40.750112830372316, 8.930699466512873e-7),
    (0.1, 36.22876263209082, 1e-6),
    (0.1, 30.63887062800405, 1.1993468186386521e-6),
    (0.1, 23.036510285681896, 1.5714416165845784e-6),
    (0.1, 22.837445924637542, 1.584893192461114e-6),
    (0.1, 17.320508075688775, 2.1300562536844447e-6),
    (0.1, 14.603973433793314, 2.5118864315095823e-6),
    (0.1, 13.02280581041226, 2.8419398537103077e-6),
    (0.1, 9.791483623609768, 3.7781975987967553e-6),
    (0.1, 9.284429639760333, 3.981071705534969e-6),
    (0.1, 7.36194280611662, 5.109136821319283e-6),
    (0.1, 5.947353627019855, 6.30957344480193e-6),
    (0.1, 5.535238985626913, 6.827122654861981e-6),
    (0.1, 4.161791450287818, 9.18627914562949e-6),
    (0.1, 3.8238639591838774, 1e-5),
    (0.1, 3.129134644531898, 1.2482821607816822e-5),
    (0.1, 2.475753088277439, 1.584893192461114e-5),
    (0.1, 2.3527088612123075, 1.6813321581146697e-5),
    (0.1, 1.7689360204744258, 2.294015464650834e-5),
    (0.1, 1.6201389072592205, 2.5118864315095822e-5),
    (0.1, 1.3300135414628023, 3.168060528558289e-5),
    (0.1, 1.0756784679359217, 3.9810717055349695e-5),
    (0.1, 1.0, 4.3676489470296766e-5),
    (0.2, 256.2807710325777, 1e-7),
    (0.2, 256.2807710325777, 1e-7),
    (0.2, 225.56161320737237, 1.1364696122304818e-7),
    (0.2, 169.59347117570735, 1.5123953747595992e-7),
    (0.2, 161.96443250931463, 1.584893192461114e-7),
    (0.2, 127.51258982610179, 2.0055942186273017e-7),
    (0.2, 102.14031138870588, 2.5118864315095823e-7),
    (0.2, 95.87315155141827, 2.6824058482269253e-7),
    (0.2, 72.08434242404265, 3.5497809161196566e-7),
    (0.2, 64.36502993948046, 3.981071705534969e-7),
    (0.2, 54.19820188053225, 4.7318691333941774e-7),
    (0.2, 40.750112830372316, 6.294673535524286e-7),
    (0.2, 40.65347306971581, 6.30957344480193e-7),
    (0.2, 30.63887062800405, 8.365144655244301e-7),
    (0.2, 25.639846151723688, 1e-6),
    (0.2, 23.036510285681896, 1.1178082096049572e-6),
    (0.2, 17.320508075688775, 1.4822120906896504e-6),
    (0.2, 16.194502664483952, 1.584893192461114e-6),
    (0.2, 13.022805810412258, 1.9836209570038384e-6),
    (0.2, 10.253734393489943, 2.5118864315095823e-6),
    (0.2, 9.791483623609768, 2.6384100777932896e-6),
    (0.2, 7.3619428061166206, 3.52472795994203e-6),
    (0.2, 6.510512597897407, 3.981071705534969e-6),
    (0.2, 5.535238985626911, 4.726606709275465e-6),
    (0.2, 4.161791450287818, 6.2733360853971586e-6),
    (0.2, 4.137573074316612, 6.30957344480193e-6),
    (0.2, 3.129134644531898, 8.45758051600518e-6),
    (0.2, 2.643460952791566, 1e-5),
    (0.2, 2.3527088612123075, 1.1347301376863671e-5),
    (0.2, 1.7689360204744258, 1.5164806245905173e-5),
    (0.2, 1.6917548381801009, 1.584893192461114e-5),
    (0.2, 1.3300135414628025, 2.056177337299613e-5),
    (0.2, 1.0917428030640703, 2.511886431509582e-5),
    (0.2, 1.0, 2.7786981581076333e-5),
    (0.3, 131.73101632587696, 1e-7),
    (0.3, 131.73101632587696, 1e-7),
    (0.3, 127.51258982610177, 1.036031327225367e-7),
    (0.3, 95.87315155141827, 1.3560818222650752e-7),
    (0.3, 82.91681924936498, 1.584893192461114e-7),
    (0.3, 72.08434242404265, 1.8379256714765547e-7),
    (0.3, 54.19820188053225, 2.4049898594786135e-7),
    (0.3, 52.162340819646396, 2.5118864315095823e-7),
    (0.3, 40.750112830372316, 3.206455565491242e-7),
    (0.3, 32.987505571800696, 3.981071705534969e-7),
    (0.3, 30.638870628004042, 4.3005407274758674e-7),
    (0.3, 23.036510285681896, 5.621415666214232e-7),
    (0.3, 20.647798282350223, 6.30957344480193e-7),
    (0.3, 17.320508075688775, 7.51704200002037e-7),
    (0.3, 13.023131078147005, 1e-6),
    (0.3, 13.02280581041226, 1.0000253540123269e-6),
    (0.3, 9.791483623609768, 1.3084998663927643e-6),
    (0.3, 8.090646325469281, 1.584893192461114e-6),
    (0.3, 7.3619428061166206, 1.7381939942229566e-6),
    (0.3, 5.535238985626912, 2.27819286184552e-6),
    (0.3, 5.01313317811982, 2.5118864315095823e-6),
    (0.3, 4.161791450287818, 2.992023861936028e-6),
    (0.3, 3.129134644531898, 3.941117170750688e-6),
    (0.3, 3.0963385271032418, 3.981071705534969e-6),
    (0.3, 2.3527088612123075, 5.1142938718884704e-6),
    (0.3, 1.883291180781875, 6.30957344480193e-6),
    (0.3, 1.7689360204744258, 6.694418211172288e-6),
    (0.3, 1.3300135414628025, 8.679788297386442e-6),
    (0.3, 1.1403134896010028, 9.999999999999999e-6),
    (0.3, 1.0, 1.1305223895756521e-5),
];

const PLOT_S: [f64; 110] = [
    -0.01, -0.0099, -0.0098, -0.0097, -0.0096, -0.0095, -0.0094, -0.0093, -0.0092, -0.0091, -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.0029, -0.0028, -0.0027, -0.0026, -0.0025, -0.0024, -0.0023, -0.0022, -0.0021, -0.002, -0.0019, -0.0018, -0.0017, -0.0016, -0.0015, -0.0014, -0.0013, -0.0012, -0.0011, -0.001, -0.0009, -0.0008, -0.0007, -0.0006, -0.0005, -0.0004, -0.0003, -0.00029, -0.00028, -0.00027, -0.00026, -0.00025, -0.00024, -0.00023, -0.00022, -0.00021, -0.0002, -0.00019, -0.00018, -0.00017, -0.00016, -0.00015, -0.00014, -0.00013, -0.00012, -0.00011, -0.0001, -9e-5, -8e-5, -7e-5, -6e-5, -5e-5, -4e-5, -3e-5, -2e-5, -1e-5, -9e-6, -8e-6, -7e-6, -6e-6, -5e-6, -4e-6, -3e-6, -2e-6, -1e-6, -9e-7, -8e-7, -7e-7, -6e-7, -5e-7, -4e-7, -3e-7, -2e-7, -1e-7, -9e-8, -8e-8, -7e-8, -6e-8, -5e-8, -4e-8, -3e-8, -2e-8, -1e-8, -9e-9, -8e-9, -7e-9, -6e-9, -5e-9, -4e-9, -3e-9, -2e-9, -1e-9, 0.0
];

const PLOT_GAMMA: [f64; 8] = [
                0.000046, 
                0.000048, 
                0.00005, 
                0.000052, 
                0.000054, 
                0.000056,
                0.000058, 
                0.00006];