import math
import random

from .coordinates import calculate_distance

def calculate_LJ(r_ij):
    """
    The LJ interaction energy between two particles.
    
    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.
    
    Parameters
    ----------
    r_ij : float
        The distance between the particles in reduced units.
        
    Returns
    -------
    pairwise_energy : float
        The pairwise Lennard Jones interaction energy in reduced units.

    """
    
    r6_term = math.pow(1/r_ij, 6)
    r12_term = math.pow(r6_term, 2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    
    return pairwise_energy

def calculate_tail_correction(num_particles, box_length, cutoff):
    """
    Calculate the long range tail correction
    """
    
    const1 = (8 * math.pi * num_particles ** 2) / (3 * box_length ** 3)
    const2 = (1/3) * (1 / cutoff)**9 - (1 / cutoff) **3
    
    return const1 * const2

def calculate_total_energy(coordinates, box_length, cutoff):
    """
    Calculate the total Lennard Jones energy of a system of particles.
    
    Parameters
    ----------
    coordinates : list
        Nested list containing particle coordinates.
    
    Returns
    -------
    total_energy : float
        The total pairwise Lennard Jones energy of the system of particles.
    """
    
    total_energy = 0
    
    num_atoms = len(coordinates)
    
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            
            #print(f'Comparings atom number {i} with atom number {j}')
            
            dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length=box_length)
            
            if dist_ij < cutoff:
                interaction_energy = calculate_LJ(dist_ij)
                total_energy += interaction_energy
            
    return total_energy

def accept_or_reject(delta_e, beta):
    """
    Accept or reject based on change in energy and temperature.
    """
    
    if delta_e <= 0:
        accept = True
    else:
        random_number = random.random()
        p_acc = math.exp(-beta*delta_e)
        
        if random_number < p_acc:
            accept = True
        else:
            accept = False
    
    return accept

def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system)
    
    
    Parameters
    ----------
    coordinates : list
        The coordinates for all the particles in the system.
        
    i_particle : int
        The particle index for which to calculate the energy.
    
    box_length : float
        The length of the simulation box.
    
    cutoff : float
        The simulation cutoff. Beyond this distance, interactions are not calculated
    
    Returns
    -------
    e_total : float
        The pairwise interaction energy of the i-th particle with all other particles in the system.
    
    """
    
    e_total = 0
    
    i_position = coordinates[i_particle]
    
    num_atoms = len(coordinates)
    
    for j_particle in range(num_atoms):
        if i_particle != j_particle:
            j_position = coordinates[j_particle]
            rij = calculate_distance(i_position, j_position, box_length)

            if rij < cutoff:
                e_pair = calculate_LJ(rij)
                e_total += e_pair
    
    return e_total