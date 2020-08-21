# include <iostream>
# include <cmath>
# include <array>
# include <vector>
# include <random>// for random numbers
# include <chrono>// for generating the seed
# include <fstream>
# include <utility>
# include <stdexcept>

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

double calculate_pair_energy(Coordinates coordinates, int i_particle, double box_length, double cutoff)
{
    /*
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
        
    */

    double e_total = 0;
    
    int num_atoms = coordinates.size();
    double r_ij;
    double e_pair;
    
    for (int i = 0; i < num_atoms; i++)
        {
        // only consider the interactions with particles that is not the i_particle
        if (i != i_particle)
            {
            r_ij = calculate_distance(coordinates[i_particle], coordinates[i], box_length);

            if (r_ij < cutoff)
                {
                e_pair = calculate_LJ(r_ij);
                e_total += e_pair;
                }
            }
        }
    return e_total;
}


std::pair<Coordinates, double> read_xyz(std::string file_path)
{
    // Opens up a file stream for input
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if(!infile.is_open())
    {   
        throw std::runtime_error("File path in read_xyz does not exist!");
    }
    
    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;
    
    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;
    
    // now the number of atoms
    infile >> num_atoms;
    
    // Uncomment to help troubleshoot
    //std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;
    
    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;
    
    for(int i = 0; i < num_atoms; i++)
    {   
        AtomCoord coord;
        
        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];
        
        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}


Coordinates run_simulation(Coordinates coordinates, double box_length, double cutoff, double reduced_temperature, int num_steps, double max_displacement = 0.1, int freq = 1000)
{
    /*
    Run a Monte Carlo simulation with the specified parameters.
    */
    
    // Reporting information
    std::vector<int> steps; //= []
    std::vector<double> energies; //= []
    double beta = 1 / reduced_temperature;
    int num_particles = coordinates.size();


    // Calculate based on the inputs
    /*try {
        std::cout << convert_F_to_C(-607.8) << std::endl;
        }
    catch (std::exception & ex)
    {
        std::cout << "Program encountered an error!" << std::endl;
        std::cout << ex.what() << std::endl;
        return 3;
    }*/

    /*try
    {
        double total_energy = calculate_total_energy(coordinates, cutoff, box_length);
    }
    catch (std::exception & ex)
    {
        std::cout << "Infinite energy calculated - particles overlapping! Halting MC simulation." << std::endl;
        std::cout << ex.what() << std::endl;
        return 1;
    }*/
    //int num_particles, double box_length, double cutoff
    double total_energy;
    total_energy = calculate_total_energy(coordinates, cutoff, box_length);
    //std::cout << total_energy << std::endl;
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);
    //std::cout << total_energy << std::endl;

    for (int step = 0; step < num_steps; step++)
        {
        //1. Randomly pick one particle in the num_particles particles.
        int random_particle = random_integer(0,num_particles);
        
        //2. Calculate the interaction energy of the selected particles with the system and store this value.
        double current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);
        
        //3. Generate a random displacement in x, y, z directions with range (-max_displacement, max_displacement).

        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);
       

        //4. Modify the coordinate of the selected particle by generated displacement.

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;
        
        //5. Calculate the new interaction energy of the new particle and store this value.
        double proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);
        
        //6. Calculate energy change and decide if this move is accepted.

        double delta_energy = proposed_energy - current_energy;
        
        bool accept = accept_or_reject(delta_energy, beta);
        

        //7. If accepted, keep movement. Else, revert to the old position.
        if (accept == true)
            {
            //std::cout << total_energy << std::endl;
            total_energy += delta_energy;
            
            }
        else
            {
            // if rejected, roll back to the origin coordinates of the selected particle.

            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
            }
        
        //8. Print the energy and store the coordinates at certain intervals.

        if (step % freq == 0)
            {
                double avg_energy = total_energy/num_particles;
                std::cout << step << ":" << avg_energy << std::endl;
                steps.push_back(step);
                energies.push_back(avg_energy);
            }
        }

    return coordinates;
}

int main(void)
{
    // set initial data
    Coordinates atom_coordinates;
    double box_length;
    std::pair<Coordinates, double> read_pair;
    read_pair = read_xyz("lj_sample_config_periodic1.txt");
    atom_coordinates = read_pair.first;
    box_length = read_pair.second;
    double cutoff = 1.0;
    double reduced_temperature = 0.9;
    int num_steps = 1000;
    double max_displacement = 0.1;
    int freq = 100;

    // call the run_mc function
    Coordinates atom_coordinates_new;
    //std::cout << atom_coordinates[0][1] << std::endl;
    atom_coordinates_new = run_simulation(atom_coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement, freq);
    //double ss;
    //double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff)
    //ss = calculate_total_energy(atom_coordinates, box_length, cutoff);
    //std::cout << ss << std::endl;

    return 0;
}

