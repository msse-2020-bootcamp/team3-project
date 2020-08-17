"""
Functions related to coordinates.
"""

import math, random

def generate_random_coordinates(num_atoms, density):
    """
    Generate the initial configuration randomly by random()
    
    Parameters
    ----------

    num_atoms : int
        The number of atoms within this system.

    density : float
        The particle density of the simulation box. 
        Here, the simulation box are regarded as a cubic box so that the box_length can be calculated directly.    
    
    Returns
    -------
    box_length : float
        The length of the cubic simulation box.
    
    coordinates : list
        The list that containing the coordinates of the atoms in the generated configuration.
    
    """
    
    # Assume that the simulation box is a cubic box and calculate the box_length
    box_volume = num_atoms / density
    box_length = box_volume ** (1/3)
    
    # Create a new empty list to store the generated coordinates.
    coordinates = []
    
    # Generate num_particles of coordinates
    for i in range(num_atoms):
        
        x_rand = random.uniform(0, box_length)
        y_rand = random.uniform(0, box_length)
        z_rand = random.uniform(0, box_length) 
        
        coordinate = [x_rand, y_rand, z_rand]
        
        coordinates.append(coordinate)
    
    return box_length, coordinates
    
    

def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.
    
    Parameters
    ----------
    coord1, coord2: list
        The atomic coordinates
    
    Returns
    -------
    distance: float
        The distance between the two points.
    """
    
    distance = 0
    for i in range(3):
        dim_dist = (coord1[i] - coord2[i]) 
        if box_length:
            dim_dist = dim_dist - box_length * round(dim_dist / box_length)
        
        dim_dist = dim_dist**2
        distance += dim_dist
    
    distance = math.sqrt(distance)
    return distance