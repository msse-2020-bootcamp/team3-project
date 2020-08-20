# include <iostream>
# include <cmath>
# include <array>
# include <vector>
# include <random>
# include <chrono>
# include <fstream>
# include <utility>
# define PI 3.1415926

typedef std::array<double, 3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;

std::default_random_engine re;

double random_double(double lower_bound, double upper_bound)
{
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    return dist(re);
}

double random_integer(int lower_bound, int upper_bound)
{           
   //dist will return [a,b] but we want [a,b)
   std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
   return dist(re);
} 

double calculate_LJ(double r_ij)
{
    /*
    The LJ interaction energy between two particles.
    
    Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced unites.
    
    Parameters
    ----------
    r_ij : double
        The distance between the particles in reduced units.
        
    Returns
    -------
    pairwise_energy : double
        The pairwise Lennard Jones interaction energy in reduced units.
    */

    double r6_term, r12_term, pairwise_energy;

    r6_term = pow(1.0/r_ij, 6);
    r12_term = pow(r6_term, 2);

    pairwise_energy = 4.0 * (r12_term - r6_term);

    return pairwise_energy;
}

double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length=100)
{
    /*
    Calculate the distance between two 3D coordinates.
    
    Parameters
    ----------
    coord1, coord2 : AtomCoord
        The atomic coordinates [x, y, z]
    
    box_length : double
        The box length. This function assumes box is a cube.
    
    Returns
    -------
    distance : double
        The distance between the two atoms.
    */

    double distance = 0;
    double dim_dist;

    for(int i=0; i<3; i++)
    {
        dim_dist = coord1[i] - coord2[i];
        dim_dist = dim_dist - box_length * round(dim_dist / box_length);
        dim_dist = pow(dim_dist, 2);
        distance += dim_dist;
    }

    distance = sqrt(distance);
    return distance;
}

double calculate_tail_correction(int num_particles, double box_length, double cutoff)
{

    double const1, const2;
    const1 = (8.0 * PI * pow(num_particles, 2)) / (3.0 * pow(box_length, 3));
    const2 = (1.0/3.0) * pow(1.0/cutoff, 9) - pow(1.0/cutoff, 3);

    return (const1 * const2);
}

double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff)
{
    /*
    Calculate the total Lennard Jones energy of a system of particles.
    
    Parameters
    ----------
    coordinates : Coordinates
        Nested list containing particle coordinates.
    
    Returns
    -------
    total_energy : double
        The total pairwise Lennard Jones energy of the system of particles.
    */

    double total_energy = 0.0;
    double dist_ij, interaction_energy;
    int num_atoms = coordinates.size();

    for (int i=0; i < num_atoms; i++)
    {
        for (int j = i+1; j < num_atoms; j++)
        {
            dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length);
            
            if(dist_ij < cutoff)
            {
                interaction_energy = calculate_LJ(dist_ij);
                total_energy += interaction_energy;
            }
        }
    }
    return total_energy;
}

bool accept_or_reject(double delta_e, double beta)
{
    /*
    Accept or reject based on change in energy and temperature.
    */
    
    bool accept;
    double random_number, p_acc;

    if (delta_e <= 0)
    {
        accept = true;
    }
    else
    {
        random_number = random_double(0.0, 1.0);
        p_acc = exp(-beta * delta_e);

        if (random_number < p_acc)
        {
            accept = true;
        }
        else
        {
            accept = false;
        }
    }
    return accept;
}