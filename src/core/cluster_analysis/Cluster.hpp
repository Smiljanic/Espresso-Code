#ifndef CLUSTER_ANALYSIS_CLUSTER_HPP
#define CLUSTER_ANALYSIS_CLUSTER_HPP

#include <vector>
#include <valarray>
#include "Vector.hpp"
#include <map>
#include <string>

#include "interaction_data.hpp"
#include "particle_data.hpp"


namespace ClusterAnalysis {

/** @brief Represents a single cluster of particles */
/** @brief Represents a single cluster of particles */
class Cluster {
  public: 
    /** @brief Ids of the particles in the cluster */
    std::vector<int> particles;
    /** @brief add a particle to the cluster */
    void add_particle(const Particle& p) {
         particles.push_back(p.p.identity);
    }
    /** @brief Ids of the particles in the spherical 
    *   shell defined by r_min and r_max */
    std::vector<int> particles_ids_in_spherical_shell(
	double r_min, double r_max);

    /** @brief Ids of all particles in the corresponding
    *   spherical shells*/
    //std::map<int,std::shared_ptr<Cluster>>  particles_ids_in_all_spherical_shells();

    
    /** @brief Calculate the center of mass of the cluster */
    Vector3d center_of_mass();
    /** @brief Longest distance between any combination of two particles */
    double longest_distance();
    /** @brief Calculate radius of gyration of the cluster */
    double radius_of_gyration();
    /** @brief Calculate the fractal dimension 
    *  N(r) via r^d, where N(r) counts the number of particles in a sphere
    *  of radius n, and d denotes the fractal dimension.
    *  The fitting is done by the Gnu Scientific Library. 
    *  @param dr:   increment for when constructing the discrete version of N(r)
    *  @param mean_sq_residual:  Mean square residual returned by the fit */
    double fractal_dimension(double dr, double& mean_sq_residual);
    /** @biref: Calculate the vector between the farthest particles in the cluster,
    *  which defines the main axis of the cluster */
    double max_radius();
};


} // namespace ClusterAnalysis



#endif

