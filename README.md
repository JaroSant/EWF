#EWF

An efficient simulator for exact Wright-Fisher diffusion and diffusion bridge paths, accounting for a wide class of selective regimes (genic, diploid and arbitrary polynomial selection), and the presence/absence of mutation.

*Dependencies*

EWF has been tested out on Ubuntu 20.04, and requires the following:

- g++ compiler (tested on version 9.4.0)
- libconfig library (tested on version 1.7.3), available from http://hyperrealm.github.io/libconfig
- boost library (tested on version 1.78.0), available from https://boost.org
- R (tested on version 4.2.1 using RStudio version 2022.07.1+554), available from https://www.r-project.org/

*Compilation*

Please ensure that the compiler and linker flags in 'Makefile' point towards where the 'libconfig' and 'boost' libraries are on your platform! 

*Configuration files*

The underlying Wright-Fisher diffusion/diffusion bridge can be configured via the 'config.cfg' file, where the mutation and selection parameters can be suitably altered. 

The configuration setup for simulating draws from the law of a _diffusion_ are found in 'configDiffusion.cfg' which allows for the start points, start times and sample times to be modified (as well as number of samples to generate and mesh size if the truncated transition density is desired). If multiple simulation setups are desired, the corresponding setup inputs need to be entered as an array. Precise instructions on input syntax can be found in the file itself.

The configuration setup for simulating draws from the law of a _diffusion bridge_ are found in 'configBridge.cfg' which allows for the start/end points and times, sampling times, number of bridges to simulate, etc. to be modified. Please see the details within the configuration file for exact instructions with regards to input syntax. The number of simulations and mesh sizes (for the truncated transition density) can also be modified.

*Running the program*

When in the root directory run 'run.sh'. This first compiles the program by invoking the makefile, and subsequently calls the program using './main horses' where the second argument invokes the demo described below. The program can be run as a diffusion or diffusion bridge simulator by changing the program invocation to simply './main', whence the program asks whether the user desires to simulate draws from a diffusion law or from a diffusion bridge law, whether they wish to condition on non-absorption and further offers the option of computing a truncation to the transition density.

*Output files*

Output for the diffusion simulator is saved using the format 'YYYY-MM-DD-HH-mmAbsDiffusionSamplesX%T%S%.txt' for the samples generated (where 'Abs' is either 'Conditioned' (if absorption is not allowed at the boundaries) or 'Unconditioned' (if absorption at the boundaries is allowed), 'T' denotes the start time and 'S' the sampling time), and 'YYYY-MM-DD-HH-mmAbsDiffusionDensityX%T%S%.txt' for the truncated transition density.

A similar system is in place for the diffusion bridge simulator, where the output is saved as 'YYYY-MM-DD-HH-mmAbsBridgeSamplesX%Z%T1%T2%S%.txt' with 'X' denoting the start point, 'Z' the end point, 'T1' the start time, 'T2' the end time and 'S' the sampling time. A similar setup is in place for the truncated transition density.
