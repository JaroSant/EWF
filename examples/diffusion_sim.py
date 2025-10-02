import EWF_pybind as EWF
import numpy as np
import matplotlib.pyplot as plt 

# Simple python script to construct Wright--Fisher class and run diffusion simulator

# Define Wright--Fisher diffusion class parameters
changepoints = np.array([0.0])
mutation_vector = np.array([[1.0, 0.75]])
non_neutral = False
sigma = np.array([0.0])
selectionSetup = 0
dominance_parameter = np.array([0.0])
selectionPolynomialDegree = 1
selectionCoefficients = np.array([[]])

# Create and initialise WrightFisher class
WF = EWF.WrightFisher(changepoints, mutation_vector, non_neutral, sigma, selectionSetup, dominance_parameter, selectionPolynomialDegree, selectionCoefficients)

# Define simulation parameters
nSim = 10000
x = 0.5
startT = 0.0
endT = 0.5
Absorption = False
Filename_sim = "EWF_diffusion_sim.txt"
verbose = True

# Run simulator
WF.DiffusionRunner(nSim, x, startT, endT, Absorption, Filename_sim, verbose)

# Define parameters for pointwise transition density evaluation
meshSize = 100
Filename_eva = "EWF_diffusion_eva.txt"

# Run transition density evaluator
WF.DiffusionDensityCalculator(meshSize, x, startT, endT, Absorption, Filename_eva, verbose)

# Load in data, create histogram and generate plot
data = np.loadtxt(Filename_sim)
density = np.loadtxt(Filename_eva)
Filename_png = "EWF_diffusion_new.png"

fig, ax = plt.subplots(1,1)
counts = plt.hist(data, bins=100, density=True)
ax.plot(density[density[:, 1] <= max(counts[0]), 0], density[density[:, 1] <= max(counts[0]), 1], linewidth = 2)
ax.set_xlabel("Sample draws")
ax.set_ylabel("Density")
ax.set_title("Histogram of sample draws (blue), truncated density (orange)")
plt.savefig(Filename_png)
plt.close()