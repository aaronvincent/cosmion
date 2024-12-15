#!/usr/bin/env python
# coding: utf-8

"""
Purpose:
    This example script processes dark matter (DM) particle collision data from cosmion.
    It reads position and velocity data from multiple simulation files, calculates
    the energy transfer during each collision, and computes the luminosity profile
    as a function of radial distance within the star. The results, including energy
    transport and associated errors, are saved to an output file.

Usage:
    python script.py <folder_number> <start_file> <end_file> <output_file> <dm_mass>

Arguments:
    folder_number : str   -> Name of the 'run#/' folder containing the input files.
    start_file    : int   -> Index of the first 'positions#.dat' file to process.
    end_file      : int   -> Index of the last 'positions#.dat' file to process.
    output_file   : str   -> Name of the output file to save the processed data.
    dm_mass       : float -> Mass of the dark matter particle in GeV.

Example:
    python script.py 1 0 10 output.txt 0.1

    This will process files from 'run1/positions0.dat' to 'run1/positions10.dat' and
    save the results to 'output.txt', assuming a dark matter mass of 0.1 GeV.

Output:
    The output file contains the following columns:
    1. Radial bins (fraction of stellar radius)
    2. Change in luminosity per radial bin (dL/dr) in erg/s/cm
    3. Error in dL/dr
    4. Cumulative luminosity (L) in erg/s
    5. Error in cumulative luminosity (L)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp

# Constants
starR = 69.57e9     # Stellar radius in cm (Sun: 69.57e9 cm, dwarf: 7.2795e9 cm, Earth: 6.378e8 cm)
starM = 1.989e33    # Stellar mass in g (Sun: 1.989e33 g, dwarf: 1.989e31 g, Earth: 5.972e27 g)
mp = 1.67262e-24    # Proton mass in g
n = 1e-15 * starM / mp  # Number of dark matter (DM) particles
GeV = 1.78266e-24   # Conversion factor for GeV to g
mdm = float(sys.argv[5]) * GeV  # Dark matter particle mass in g (read from command-line argument)

# Read the file names and output file name from the command line.
fileFolder = sys.argv[1]       # Folder name for 'run#/'
fileStart = int(sys.argv[2])   # First index of 'positions#.dat' file
fileEnd = int(sys.argv[3])     # Last index of 'positions#.dat' file
fileRange = np.arange(fileStart, fileEnd + 1, 1)  # Range of file indices

# Create a list of file names based on the folder and indices.
names = []
for suffix in fileRange:
    names.append(f'run{fileFolder}/positions{suffix}.dat')

filename = sys.argv[4]         # Output file name

# Initialize simulation time and bins
t = 0                          # Total simulation time in seconds
bins = np.linspace(0, 1, 300)  # Radial bins (fractions of stellar radius)
dr = (bins[1] - bins[0]) * starR  # Radial step size in cm

# Initialize arrays to store energy and related calculations per bin
binE = np.zeros(len(bins))     # Energy per radial bin in erg
binE2 = np.zeros(len(bins))    # Energy squared per radial bin in erg^2
binEtrap = np.zeros(len(bins)) # Integral of energy per bin
err_E = np.zeros(len(bins))    # Error in energy
N = np.zeros(len(bins))        # Number of collisions per radial bin

# Loop through each file in the list of file names
for name in names:
    # Read the data from the file, using the first 11 columns
    data = np.genfromtxt(name, usecols=range(11))
    
    # Extract position and velocity data
    x1 = data[1:, 0]
    x2 = data[1:, 1]
    x3 = data[1:, 2]
    vin1 = data[1:, 3]
    vin2 = data[1:, 4]
    vin3 = data[1:, 5]
    vout1 = data[1:, 6]
    vout2 = data[1:, 7]
    vout3 = data[1:, 8]
    dt = data[1:, 9]          # Time since last collision in seconds
    out = data[1:, 10]        # Flag indicating if the particle is outside the star (1 for outside, 0 for inside)
    
    # Compute scalar radius (distance from center) in cm
    r = np.sqrt(x1**2 + x2**2 + x3**2)
    fracR = r / starR         # Radius as a fraction of stellar radius
    
    # Compute speeds before and after collision
    vin = np.sqrt(vin1**2 + vin2**2 + vin3**2)
    vout = np.sqrt(vout1**2 + vout2**2 + vout3**2)
    
    # Calculate energy transfer for each collision in erg
    E = mdm * (vout**2 - vin**2) / 2

    # Update the total simulation time
    t += sum(dt)
    
    # Distribute energy transfer into radial bins
    binplace = np.digitize(fracR, bins)
    for i in range(len(bins)):
        places = np.where(binplace == i)[0]
        for j in places:
            if out[j] == 0:   # Only consider particles inside the star
                binE[i] += E[j]
                binE2[i] += E[j]**2
                N[i] += 1
    print('done ', name)

# Calculate luminosity and energy transport
binL = binE * n / t                # Total luminosity per radial bin in erg/s
dLdr = -binL / dr                  # Change in luminosity per radial bin
L = sp.cumulative_trapezoid(binL, initial=0)  # Cumulative luminosity as a function of radius

# Calculate error on energy transport
err_dEdr = np.sqrt(binE2 - binE**2 / N)
err_dLdr = err_dEdr * n / t / dr

# Calculate errors for cumulative luminosity
for i in range(len(bins)):
    binEtrap[i] = np.trapz(binE[:i+1])
    err_E[i] = np.sqrt(sum(binE2[:i+1]) - binEtrap[i]**2 / N[i])
err_L = err_E * n / t

# Print the total number of collisions
print('N =', sum(N))

# Prepare output data and save to file
reduced = np.array([bins, dLdr, err_dLdr, L, err_L])
np.savetxt(filename, reduced)
