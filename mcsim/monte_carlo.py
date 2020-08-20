"""
<<<<<<< HEAD
Functions for running a Monte Carlo Simulation
"""

import math
import random
import os

from .energy import *
from .coordinates import *

import math
import os
import random

import matplotlib.pyplot as plt

def calculate_LJ(r_ij):
    """
    The LJ interaction energy between two particles.
    
    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced unites.
    
    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.
        
    Returns
    -------
    pairwise_energy : floa
        The pairwise Lennard Jones interaction energy in reduced units.
    
    """
    
    r6_term = math.pow(1/r_ij,6)
    r12_term = math.pow(r6_term,2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    
    return pairwise_energy



def calculate_distance(coord1,coord2,box_length=None):
    """
    Calculate the distance between two 3D coordinates.
    
    Parameters
    ----------
    coord1, coord2 : list
        The atomic coordinates [x, y, z]
    
    box_length : float, optional
        The box length. This function assumes box is a cube.
    
    Returns
    -------
    distance : float
        The distance between the two atoms.
    """
    
    distance = 0
    vector = [0,0,0]
    
    for i in range(3):
        vector[i] = coord1[i] -coord2[i]
        if box_length is None:
            pass
        else:
            if vector[i] > box_length/2:
                vector[i] -= box_length
            elif vector[i] < -box_length/2:
                vector[i] += box_length
        dim_dist = vector[i] ** 2
        distance += dim_dist
        
    distance = math.sqrt(distance)
    
    return distance



def calculate_tail_correction(cutoff, box_length, num_atoms):
    """
    Calculate the tail correction.
    
    Parameters
    ----------
    cutoff : float
        The curoff distance.
        
    box_length : float
        The length of the cell.
       
    num_atoms : int
        Number of atoms in a given system.
       
    Returns
    -------
    tail_co_LJ : float
        A float number that shows the value of tail correction energy for the given system.
    """
    
    tail_co_LJ = 0
    
    coeff = 0
    
    r3 = math.pow(1/cutoff,3)
    r9 = math.pow(r3,3)
    
    coeff = 8 * math.pi * (num_atoms ** 2)/(3 * (box_length ** 3))

    tail_co_LJ = coeff * (r9/3 - r3)
    
    return tail_co_LJ



def calculate_total_energy(coordinates, cutoff=3, box_length=None):
    """
    Calculate the total Lennard Jones energy of a system of particles.
    
    Parameters
    ----------
    coordinates : list
        Nested list containing particle coordinates.
        
    cutoff : float
        A criteria distance for intermolecular interaction truncation
    
    box_length : float, optional
        The box length. This function assumes box is a cube.
        
    Returns
    -------
    total_energy : float
        The total pairwise Lennard Jones energy of the system of particles.
    """
    
    total_energy = 0
    
    num_atoms = len(coordinates)
    
    for i in range(num_atoms):
        for j in range(i+1,num_atoms):
            
            
            dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length)
            
            if dist_ij < cutoff:
                interaction_energy = calculate_LJ(dist_ij)
                total_energy += interaction_energy
            
    return total_energy



def read_xyz(filepath):
    """
    Reads coordinates from an xyz file.
    
    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.
       
    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates
    """
    
    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())
        coordinates = f.readlines()
    
    atomic_coordinates = []
    
    for atom in coordinates:
        split_atoms = atom.split()
        
        float_coords = []
        
        # We split this way to get rid of the atom label.
        for coord in split_atoms[1:]:
            float_coords.append(float(coord))
            
        atomic_coordinates.append(float_coords)

        
    
    return atomic_coordinates, box_length

def run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement, freq=1000):
    """
    Run a Monte Carlo simulation with the specified parameters. 
    """

    # Reporting information
    steps = []
    energies = []
    all_coordinates = []

    # Calculated quantities
    beta = 1/reduced_temperature
    num_particles = len(coordinates)

    # Calculated based on simulation inputs
    total_energy = calculate_total_energy(coordinates, box_length, cutoff)
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff)


    for step in range(num_steps):
        
        # 1. Randomly pick one of num_particles particles
        random_particle = random.randrange(num_particles)
        
        # 2. Calculate the interaction energy of the selected particle with the system. Store this value.
        current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        
        # 3. Generate a random x, y, z displacement range (-max_displacement, max_displacement) - uniform distribution

    
    return atomic_coordinates, box_length



def accept_or_reject(delta_e, beta):
    """
    Accept or reject based on change in energy and temperature.
    
    Parameters
    ----------
    delta_e : float
       The change of the system's energy.
       
    beta : float
        1 over the temperature.
       
    Returns
    -------
    accept : bool
        Accept the move or not.
    
    """
    if delta_e <= 0:
        accept = True
        
    else:
        rand_num = random.random()
        p_acc = math.exp(-beta * delta_e)
        
        if rand_num < p_acc:
            accept = True
        else:
            accept = False
            
    return accept



def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of the particles with its environment (all other particles in the system)
    
    Parameters
    ----------
    coordinates : list
        The coordinates for all the particles within the system.
        
    i_particle : int
        The particle index for which to calculate the energy.
        
    box_length : float
        The length of the simulation box.
        
    cutoff : float
        The simulation cutoff. Beyond this distance, interactions are not calculated.
        
    Returns
    -------
    e_total : float
        The pairwise interaction energy with the i_th particle with all other particles in the system.
        
    """
    e_total = 0
    
    num_atoms = len(coordinates)
    
    for i in range(num_atoms):
        
        # only consider the interactions with particles that is not the i_particle
        if i != i_particle:
            
            r_ij = calculate_distance(coordinates[i_particle], coordinates[i], box_length)
            
            if r_ij < cutoff:
            
                e_pair = calculate_LJ(r_ij)
                e_total += e_pair
    
    return e_total



def initial_config(box_volume, num_particles):
    """
    Generate the initial configuration randomly by random()
    
    Parameters
    ----------
    box_volume : float
        The volume of the simulation box. 
        Here, the simulation box are regarded as a cubic box so that the box_length can be calculated directly.
        
    num_particles : int
        The number of particles within this system.
        
    Returns
    -------
    box_length : float
        The length of the cubic simulation box.
    
    coordinates : list
        The list that containing the coordinates of the atoms in the generated configuration.
    
    """
    
    # Assume that the simulation box is a cubic box and calculate the box_length
    box_length = box_volume ** (1/3)
    
    # Create a new empty list to store the generated coordinates.
    coordinates = []
    
    # Generate num_particles of coordinates
    for i in range(num_particles):
        
        x_rand = random.uniform(0, box_length)
        y_rand = random.uniform(0, box_length)
        z_rand = random.uniform(0, box_length) 
        
        coordinate = [x_rand, y_rand, z_rand]
        
        coordinates.append(coordinate)
    
    return box_length, coordinates



def run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement = 0.1, freq = 1000):
    """
    Run a Monte Carlo simulation with the specified parameters.
    """
    
    # Reporting information
    steps = []
    energies = []
    beta = 1 / reduced_temperature
    num_particles = len(coordinates)


    # Calculate based on the inputs
    try:
        total_energy = calculate_total_energy(coordinates, cutoff, box_length)
    except ZeroDivisionError:
        raise ZeroDivisionError("Infinite energy calculated - particles overlapping! Halting MC simulation.")

    total_energy += calculate_tail_correction(cutoff, box_length, num_particles)

    for step in range(num_steps):
        
        # 1. Randomly pick one particle in the num_particles particles.
        random_particle = random.randrange(0,num_particles)
        
        # 2. Calculate the interaction energy of the selected particles with the system and store this value.
        current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        
        # 3. Generate a random displacement in x, y, z directions with range (-max_displacement, max_displacement).

        x_rand = random.uniform(-max_displacement, max_displacement)
        y_rand = random.uniform(-max_displacement, max_displacement)
        z_rand = random.uniform(-max_displacement, max_displacement)
       

        # 4. Modify the coordinate of the selected particle by generated displacement.

        coordinates[random_particle][0] += x_rand
        coordinates[random_particle][1] += y_rand
        coordinates[random_particle][2] += z_rand
        
        # 5. Calculate the new interaction energy of the new particle and store this value.
        proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        
        # 6. Calculate energy change and decide if this move is accepted.

        delta_energy = proposed_energy - current_energy
        
        accept = accept_or_reject(delta_energy, beta)
        

        # 7. If accepted, keep movement. Else, revert to the old position.
        if accept == True:
            total_energy += delta_energy
        else:
            # if rejected, roll back to the origin coordinates of the selected particle.

            coordinates[random_particle][0] -= x_rand
            coordinates[random_particle][1] -= y_rand
            coordinates[random_particle][2] -= z_rand

        # 8. Print the energy and store the coordinates at certain intervals

        # 8. Print the energy and store the coordinates at certain intervals.

        if step % freq == 0:
            print(step, total_energy/num_particles)
            steps.append(step)
            energies.append(total_energy/num_particles)



    return coordinates

