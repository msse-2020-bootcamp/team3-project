# include <iostream>
# include <cmath>
# include <array>
# include <vector>
# include <random>
# include <chrono>
# include <fstream>
# include <utility>
# include <time.h>
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
    AtomCoord vec = {0, 0, 0};
    double dim_dist;

    for(int i=0; i<3; i++)
    {
        vec[i] = coord1[i] - coord2[i];
        if (vec[i] > box_length/2)
        {
            vec[i] -= box_length;
        }
        else if (vec[i] < - box_length/2)
        {
            vec[i] += box_length;
        }
        dim_dist = pow(vec[i], 2);
        distance += dim_dist;
    }

    distance = sqrt(distance);
    return distance;
}

double calculate_tail_correction(int num_atoms, double box_length, double cutoff)
{
    /*
    Calculate the tail correction.
    
    Parameters
    ----------
    num_atoms : int
        Number of atoms in a given system.

    cutoff : float
        The cutoff distance.   

    box_length : float
        The length of the cell.
       
    Returns
    -------
    tail_co_LJ : float
        A float number that shows the value of tail correction energy for the given system.
    */
    double tail_co_LJ, coeff;
    double r3 = pow(1/cutoff, 3);
    double r9 = pow(r3, 3);
    
    coeff = 8.0 * PI * pow(num_atoms, 2)/(3.0 * pow(box_length, 3));
    tail_co_LJ = coeff * (r9 / 3.0 - r3);

    return tail_co_LJ;
}

double calculate_total_energy(Coordinates coordinates, double box_length=100, double cutoff=3)
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
    
    */

    double e_total = 0;
    AtomCoord i_position = coordinates[i_particle];
    AtomCoord j_position;
    int num_atoms = coordinates.size();
    double rij, e_pair;

    for (int j_particle=0; j_particle < num_atoms; j_particle++)
    {
        if (i_particle != j_particle)
        {
            j_position = coordinates[j_particle];
            rij = calculate_distance(i_position, j_position, box_length);

            if (rij < cutoff)
            {
                e_pair = calculate_LJ(rij);
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

std::pair< std::vector<int>, std::vector<double> > run_simulation(Coordinates coordinates, double box_length, double cutoff, double reduced_temperature, int num_steps, double max_displacement, int freq=1000)
{
    
    // Reporting information
    std::vector<int> steps;
    std::vector<double> energies;
    double beta = 1 / reduced_temperature;
    int num_particles = coordinates.size();
    int random_particle;
    double current_energy, proposed_energy, delta_energy;
    double x_rand, y_rand, z_rand;
    bool accept;

    // Calculate based on the inputs
    double total_energy = calculate_total_energy(coordinates, box_length, cutoff);
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);

    for (int step=0; step < num_steps; step++)
    {
        random_particle = random_integer(0, num_particles);
        current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);
        x_rand = random_double(-max_displacement, max_displacement);
        y_rand = random_double(-max_displacement, max_displacement);
        z_rand = random_double(-max_displacement, max_displacement);

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;

        proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        delta_energy = proposed_energy - current_energy;

        accept = accept_or_reject(delta_energy, beta);

        if (accept == true)
        {
            total_energy += delta_energy;
        }
        else
        {
            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
        }

        if (step % freq == 0)
        {
            steps.push_back(step);
            energies.push_back(total_energy / num_particles);
        }
    }
    return std::make_pair(steps, energies);

}

int main(void)
{
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());
    //time_t start, end, due;
    //start = clock();
    std::pair<Coordinates, double> xyz_info = read_xyz('..',"lj_sample_configurations/lj_sample_config_periodic1.txt");
    Coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;

    std::pair< std::vector<int>, std::vector<double> > sim_info = run_simulation(coords, box_length, 3, 0.9, 10000, 0.1);
    //start = clock();
    std::vector<int> simulation_steps = sim_info.first;
    std::vector<double> simulation_energies = sim_info.second;
    std::cout << simulation_energies[0] << std::endl;
    //due = end - start;
    //std::cout << due/1000000.0 << std::endl;
    return 0;
}