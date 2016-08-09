/*
  Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "cluster_analysis.hpp"
#include "interaction_data.hpp"
#include "utils.hpp"
#include "config.hpp"
#include "errorhandling.hpp"

void ClusterStructure::clear() {
 clusters.clear();
 cluster_id.clear();
 cluster_identities.clear();
}

// Analyze the cluster structure of the given particles
void ClusterStructure::analyze_pair()
{
  // clear data structs
  clear();
  
  // Iterate over pairs
  for (int i=0;i<=max_seen_particle;i++) {
    if (! local_particles[i]) continue;
    for (int j=i+1;j<=max_seen_particle;j++) {
      if (! local_particles[j]) continue;
      add_pair(*local_particles[i],*local_particles[j]); // maybe no *
    }
  }
}


//MILENA-HACK: missing part
void ClusterStructure::analyze_energy()
{
  // clear data structs
  clear();
 
 printf("ANALYZE ENERGY------"); 
  // Iterate over pairs
  for (int i=0;i<=max_seen_particle;i++) {
    if (! local_particles[i]) continue;
    for (int j=i+1;j<=max_seen_particle;j++) {
      if (! local_particles[j]) continue;
      add_pair(*local_particles[i],*local_particles[j]); // maybe no *
    }
  }
}


// ---------------begin of gemetrical analysis-------------------

std::vector<double>  Cluster::center_of_mass(Particle& p) 
{
  std::vector<double> com; //initialized com
  double temp[3] = {0,0,0}; //initialized temporary array triplet storing sum of particle positions of particle
  for (auto const& it : particles) { //iterate over all particles within a cluster
    int pid = particles[it]; //ID of the indexed particle from (vector) particles
    for (int i=0; i<3; i++){ 
      temp[i] += local_particles[pid]->r.p[i]; //ad to the temporary array current particle positions
    }
  }
  for (int i=0; i<3; i++) {
    com[i] = temp[i]*(1.0/particles.size()); //divide by nuber of particles in aggregate
  }
  return com;
}


double Cluster::largest_distance(Particle& p)
{
  double ld = 0.0; //the longest distance
  double ld_vec[3] ={0,0,0}; //longest distance vector
  double position[3] = {0,0,0}; //position of current particle
//calculate com  
  double com[3]; //center of mass
  double temp[3] = {0,0,0}; // intermediate variable storing sums of positions for coordiantes x,y,z
  for (auto const& it : particles) { //iterate over list of particles within an aggregate
    int pid = particles[it]; //take ID of the indexed particle from (vector) particles
    for (int i=0; i<3; i++){ 
      temp[i] += local_particles[pid]->r.p[i]; //add current particle position to the temp
    }
  } 
  for (int i=0; i<3; i++) {
    com[i] = temp[i]*(1.0/particles.size()); //divide sum ov all positions with number of particles
  } 
//compare the distance of each particle from the c_o_m to get the longest    
  for (auto const& it2 : particles) { //again iterate over particles within an aggregate
    int pid = particles[it2]; //ID of the indexed particle from (vector) particles
    for (int i=0; i<3; i++){ 
     ld_vec[i] = com[i]-local_particles[pid]->r.p[i]; //calculate relative particle position to the com
    }
    if ((sqrlen(ld_vec))>ld) //compare that distance with the longest distance
      ld=sqrlen(ld_vec); //save bigger value as longest distance - ld
  }
  return ld;
}


double Cluster::radius_of_gyration(Particle& p)
{
  double rg2;
  int cluster_size = particles.size();
  double position[3] = {0,0,0};
  double com[3];
// routine to calculate com 
  double temp[3] = {0,0,0}; //initialized position of particle
  for (auto const& it : particles) {
    int pid = particles[it]; //ID of the indexed particle from (vector) particles
    for (int i=0; i<3; i++){ 
      temp[i] += local_particles[pid]->r.p[i];
    }
  } 
  for (int i=0; i<3; i++) {
    com[i] = temp[i]*(1.0/particles.size()); 
  } 
//compare the distance of each particle from the c_o_m to get the longest    
  double current[3];
  double current_modul;
  double current2;
  for (auto const& it3 : particles) {
    int pid = particles[it3]; //ID of the indexed particle from (vector) particles
// calculate distance between com and pid and store in variable current  
    get_mi_vector(current, com, local_particles[pid]->r.p);
// calculate square length of this distance  
    current2 += sqrlen(current)*sqrlen(current);
  }      
// divide with number of particles 
  rg2 = current2/particles.size(); 
//return square root of it
  return sqrt(rg2);
}


double Cluster::fractal_dimension(Particle& p)
{
  double df = 3.0;    //maximum df for spheres
  double ppos[3];     //current particle position
  double p_to_com[3]; //vector of the particle to the com
  double distance;    //distance of the particle from the center of the mass of the agglomerate
  int pid;            // particle ID
//  calculate com of the aggregate
  double com[3];
  double temp[3] = {0,0,0}; //initialized position of particle
  for (auto const& it : particles) {
    int pid = particles[it]; //ID of the indexed particle from (vector) particles
    for (int i=0; i<3; i++){ 
      temp[i] += local_particles[pid]->r.p[i];
    }   
  }   
  for (int i=0; i<3; i++) {
    com[i] = temp[i]*(1.0/particles.size()); 
  } 

  std::vector<double> distances; //list of distances , same size as particles
  std::vector<double> diameters; //all diameters=(radii*2) of circles around the com of aggregate
  std::vector<int> pcounts; // numbers of particles within given diameters
  int rad = 0;
  int k = 0;  
// iterate over particles within an aggregate 
  for (auto const& it : particles) {
    diameters.push_back(rad*2.0); //diameters are taken as doubled counter rad=0,1,2,3,..,particles.size()
// get particle's ID
    pid = particles[it];
    for (int i=0; i<3; i++){ 
// get particle position
      ppos[i] = local_particles[pid]->r.p[i];
// calculate particle vector positions from the COM
      p_to_com[i] = com[i]-ppos[i]; 
    }
//calculate particle distance from the COM 
    distance = sqrlen(p_to_com);
    distances.push_back(distance); //add distance from the current particle to the com in the distances vectors
    rad+=1;
  }
//now calculate pcounts for all diameters or iterate over distances or do while loop; k is initialized as 0
  for (auto const& co : diameters ) { //iterate over diameters
    double diam = diameters[co];
    int pcount=0;
    for (auto const& it : distances) { //go over distances
      double dist = distances[it]; //ID of the indexed particle from (vector) particles
      if (dist<diam)
      {
        pcount+=1; //count number of particcles within given diameter
      }     
    pcounts.push_back(pcount);
    }
  }   
	  
//calculate Df using linear regression on the logarithms of diameters [__std::vector<double> diameters__] and num of particles [__std::vector<int> pcounts__] within the diameters
  std::vector<double> log_diameters;
  std::vector<double> log_pcounts;
  for (auto const& co : diameters ) { 
    log_diameters[co]=log(diameters[co]); //save the logarithms of diameters and num of particles --> do it in a more fashionable way : maybe with map
    log_pcounts[co]=log(pcounts[co]);
  }

#ifdef GSL
  df=3.0;
  double c1, c2, c3, c4, c5, c6;
  if (diameters.size() > 1) : 
    gsl_fit_linear(x,1,y,1,c1,c2,c3,c4,c5,c6);
    gsl_fit_linear (&diameters.front(), 1, &pcounts.front(), 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);  
  return df;

#else
  runtimeErrorMsg()<< "GSL (gnu scientific library) is not found.";
#endif
}
// -----------------end of gemetrical analysis-------------------



void ClusterStructure::add_pair(Particle& p1, Particle& p2) {
// * check, if there's a neighbor
 //   * No: Then go on to the next particle
 // * Yes: Then if
 //   * One of them belongs to a cluster, give the other one the same cluster
 //     id.
 //   * None of them belongs to a cluster: Give them both a new cluster id
 //   * Both belong to different clusters: Mark the clusters as identical 
 //   * so that they can be put together later
  if (! nc) {
    runtimeErrorMsg() << "No cluster criterion defined"; 
    return;
  }
  if (nc->are_neighbors(p1,p2)) {
     if // None belongs to a cluster
     ((cluster_id.find(p1.p.identity)==cluster_id.end()) && (cluster_id.find(p2.p.identity)==cluster_id.end()))
    {  
      // Both particles belong to the same, new cluster
      int cid=get_next_free_cluster_id();

      // assign the 
      cluster_id[p1.p.identity]=cid;
      cluster_id[p2.p.identity]=cid;
    }
    else if // j belongs to a cluster but i doesn't
      ((cluster_id.find(p1.p.identity)==cluster_id.end()) 
      && 
      (cluster_id.find(p2.p.identity) !=cluster_id.end()))
    {
     // Give p1 the same cluster id as p2
     cluster_id[p1.p.identity]=find_id_for(cluster_id[p2.p.identity]);
    }
    else if // i belongs to a cluster but j doesn't
      ((cluster_id.find(p2.p.identity)==cluster_id.end()) 
      && 
      (cluster_id.find(p1.p.identity) !=cluster_id.end()))
    {
     // give p2 the cluster id from p1
     cluster_id[p2.p.identity]=find_id_for(cluster_id[p1.p.identity]);
    }
    else if // Both belong to different clusters
      (cluster_id[p1.p.identity] != cluster_id[p2.p.identity])
    {
     // Clusters of p1 and p2 are one and the same. Add an identity to the list
     // The lower number must be inserted as first value of tjhe pair
     // because the substituions later have to be done in ascending order
     if (cluster_id[p1.p.identity]<cluster_id[p2.p.identity])
     {
       cluster_identities[find_id_for(cluster_id[p2.p.identity])] =find_id_for(cluster_id[p1.p.identity]);
     }
     else
     {
       cluster_identities[find_id_for(cluster_id[p1.p.identity])] =find_id_for(cluster_id[p2.p.identity]);
     }
     
     // Connected clusters will be merged later
    }
    // The case for both particles being in the same cluster does not need to be
    // treated.
  }
}

void ClusterStructure::merge_clusters() {
  // Relabel particles according to the cluster identities map
  // Also create empty cluster objects for the final cluster id
  for (auto it : cluster_id) { 
    // particle id is in it.first and cluster id in it.second
    // We change the cluster id according to the cluster identities
    // map
    int cid=find_id_for(it.second);
    it.second =cid;
    // Empty cluster object
    if (clusters.find(cid)==clusters.end()) {
      clusters[cid]=Cluster();
    }
  }

  
  // Now fill the cluster objects with particle ids
  // Iterate over particles, fill in the cluster map 
  // to each cluster particle the corresponding cluster id 
  for (auto it : cluster_id) {
    clusters[it.second].particles.push_back(it.first);
  }
}


int ClusterStructure::find_id_for(int x)
{
 while (cluster_identities.find(x)!=cluster_identities.end())
 {
  x =cluster_identities[x];
 }
 return x;
}

int ClusterStructure::get_next_free_cluster_id(){
  //iterate over cluster_id'
  int max_seen_cluster = 1;
  for (auto it : clusters){
    int cid=it.first;
    if (max_seen_cluster <= cid ) {
      max_seen_cluster=cid+1;
    }
  }
  return max_seen_cluster;
}


ClusterStructure cluster_structure;

ClusterStructure& cluster_analysis() {
  return cluster_structure;
}

void ClusterStructure::set_criterion(NeighborCriterion* c) {
  if (nc)
  {
    delete nc;
    nc=0;
  }
  nc=c;
}

   

 

