import EWF_pybind as EWF
import numpy as np
import matplotlib.pyplot as plt

# Simple python script to construct Wright--Fisher class and run diffusion simulator

# Define Wright--Fisher diffusion class parameters
mutation_vector = np.array([0.5, 0.5])
non_neutral = False
sigma = 0.0
selectionSetup = 0
dominance_parameter = 0.0
selectionPolynomialDegree = 1
selectionCoefficients = np.array([])

# Create and initialise WrightFisher class
WF = EWF.WrightFisher(mutation_vector, non_neutral, sigma, selectionSetup, dominance_parameter, selectionPolynomialDegree, selectionCoefficients)

# Define simulation parameters
nSim = 10000
x = 0.5
z = 0.5
startT = 0.0
endT = 0.5
sampleT = 0.25
Absorption = False
Filename_sim = "EWF_diffusion_bridge_sim.txt"

# Run simulator
WF.BridgeDiffusionRunner(nSim, x, z, startT, endT, sampleT, Absorption, Filename_sim)

# Define parameters for pointwise transition density evaluation
meshSize = 100
Filename_eva = "EWF_diffusion_bridge_eva.txt"

# Run transition density evaluator
WF.BridgeDiffusionDensityCalculator(meshSize, x, z, startT, endT, sampleT, Absorption, Filename_eva)

# Load in data, create histogram and generate plot
data = np.loadtxt(Filename_sim)
density = np.loadtxt(Filename_eva)
Filename_png = "EWF_diffusion_bridge.png"

fig, ax = plt.subplots(1,1)
counts = plt.hist(data, bins=50, density=True)
ax.plot(density[density[:, 1] <= max(counts[0]), 0], density[density[:, 1] <= max(counts[0]), 1], linewidth = 2)
ax.set_xlabel("Sample draws")
ax.set_ylabel("Density")
ax.set_title("Histogram of sample draws (blue), truncated density (orange)")
plt.savefig(Filename_png)
plt.close()