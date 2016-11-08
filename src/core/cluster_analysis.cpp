#include "cluster_analysis.hpp"
#include "interaction_data.hpp"
#include <algorithm>
#include "config.hpp"
#include <gsl/gsl_fit.h>

void ClusterStructure::clear() {
 clusters.clear();
 cluster_id.clear();
 cluster_identities.clear();
}

inline bool ClusterStructure::part_of_cluster(const Particle& p)
{
 if (cluster_id.find(p.p.identity)==cluster_id.end()) return false;
 return true;
}

// Analyze the cluster structure of the given particles
void ClusterStructure::analyze_pair()
{
 // printf("Came to the analyze_pair part!\n");
  // clear data structs
  clear();
  
  // Iterate over pairs
  for (int i=0;i<=max_seen_particle;i++) {
    if (! local_particles[i]) continue;
    for (int j=i+1;j<=max_seen_particle;j++) {
      if (! local_particles[j]) continue;
      add_pair(*local_particles[i],*local_particles[j]); 
    }
  }
  merge_clusters();
}




void ClusterStructure::analyze_bonds() {

  //printf("Came to the analyze_bond part!\n");
clear();
for (int i=0;i<=max_seen_particle;i++) {
  if (local_particles[i]) {
    auto p=local_particles[i];
    int j=0;
    while (j<p->bl.n) {
      int bond_type=p->bl.e[j];
      int partners =bonded_ia_params[bond_type].num;
      if (partners!=1) {
        j+=1+partners;
        continue;
      }
      // We are only here if bond has one partner
      add_pair(*p,*(local_particles[p->bl.e[j+1]]));
      j+=2; // Type id + one partner
    }
  }
}
merge_clusters();
}



void ClusterStructure::add_pair(Particle& p1, Particle& p2) {
// * check, if there's a neighbor
 //   * No: Then go on to the next particle
 // * Yes: Then if
 //   * One of them belongs to a cluster, give the other one the same cluster
 //     id.
 //   * None of them belongs to a cluster: Give them both a new cluster id
 //   * Both belong to different clusters: Mark the clusters as identical 
 //   * so that they can be put together later
  //printf("Came to the add_pair part!\n");
  if (! nc) {
    runtimeErrorMsg() << "No cluster criterion defined"; 
    return;
  }
  if (nc->are_neighbors(p1,p2)) {
     
     if // None belongs to a cluster
     ((!part_of_cluster(p1)) && (!part_of_cluster(p2)))
    {  
      // Both particles belong to the same, new cluster
      int cid=get_next_free_cluster_id();

      // assign the cluster_ids 
      cluster_id[p1.p.identity]=cid;
      cluster_id[p2.p.identity]=cid;
    }
    else if // p2 belongs to a cluster but p1 doesn't
      (part_of_cluster(p2) && !part_of_cluster(p1))
    {
     // Give p1 the same cluster id as p2
     cluster_id[p1.p.identity]=find_id_for(cluster_id.at(p2.p.identity));
    }
    else if // i belongs to a cluster but j doesn't
      (part_of_cluster(p1) && !part_of_cluster(p2))
    {
     // give p2 the cluster id from p1
     cluster_id[p2.p.identity]=find_id_for(cluster_id.at(p1.p.identity));
    }
    else if // Both belong to different clusters
      (part_of_cluster(p1) && part_of_cluster(p2) &&
       cluster_id.at(p1.p.identity)!=cluster_id.at(p2.p.identity))
    {
     // Clusters of p1 and p2 are one and the same. Add an identity to the list
     // The higher number must be inserted as first value of tjhe pair
     // because the substituions later have to be done in descending order
     int cid1=find_id_for(cluster_id.at(p1.p.identity));
     int cid2=find_id_for(cluster_id.at(p2.p.identity));
     if (cid1>cid2)
     {
       cluster_identities[cid1] =cid2;
     }
     else if (cid1<cid2) 
     {
       cluster_identities[cid2]=cid1;
     }
     // else do nothing. The clusters are already noted for merging.
     // Connected clusters will be merged later
    }
    // The case for both particles being in the same cluster does not need to be
    // treated.
  }
}

void ClusterStructure::merge_clusters() {
  // Relabel particles according to the cluster identities map
  // Also create empty cluster objects for the final cluster id
  //printf("Came to MERGE_CLUSTERS !\n"); 
  std::map<int,int> to_be_changed;
  for (auto it : cluster_id) { 
    // particle id is in it.first and cluster id in it.second
    // We change the cluster id according to the cluster identities
    // map
    int cid=find_id_for(it.second);
    // We note the list of changes here, so we don't modify the map
    // while iterating
    to_be_changed[it.first]=cid;
    // Empty cluster object
    if (clusters.find(cid)==clusters.end()) {
      clusters[cid]=Cluster();
    }
  }
  
  // Now act on the changes marked in above iteration
  for (auto it : to_be_changed) {
    cluster_id[it.first]=it.second;
    //printf("MERGE_CLUSTERS %d\n", it); 
   
  }
  
  // Now fill the cluster objects with particle ids
  // Iterate over particles, fill in the cluster map 
  // to each cluster particle the corresponding cluster id 
  for (auto it : cluster_id) {
    clusters[it.second].particles.push_back(it.first);
  }

  // Sort particles ids in the clusters
  for (auto c : clusters) {
    std::sort(c.second.particles.begin(),c.second.particles.end());
  }
   
}



 // Geometry analysis

//Center of mass of an aggregate
std::vector<double>  Cluster::calculate_cluster_center_of_mass() 
{
 std::vector<double> com; //initialized com
  for (int i=0; i<3; i++) {
    com.push_back(0.0);
  }
  //printf("initialized center of mass is: %f, %f, %f\n", com[0], com[1], com[2]);

  // due to the periodic boundary conditions, positions have to be folded 
  // Instead using fold_coordinate() from grid.hpp, the position of the first
  // particle of the cluster is taken as reference, and for the other particles 
  // distance is calculated with get_mi_vector(reference, current part), added to 
  // the reference and finally divided with num of part. in cluster 

  double reference_position[3] = {0.0};
  double relative_to_reference[3] = {0.0};
  double sum_of_distances[3] = {0.0};

  // accessing first particle of an aggregate
  for (int i=0; i<3; i++)
    reference_position[i] = local_particles[particles[0]]->r.p[i];
//  printf("the reference particle is: %d at %f , %f, %f\n", local_particles[particles[0]]->p.identity, reference_position[0], reference_position[1],reference_position[2]); 
  for (int it : particles)  //iterate over all particles within a cluster
  {
//    printf("%d: [%f, %f, %f]\n", local_particles[particles[it]]->p.identity, local_particles[it]->r.p[0], local_particles[it]->r.p[1], local_particles[it]->r.p[2]);
    get_mi_vector(relative_to_reference, local_particles[it]->r.p, reference_position); //add current particle positions
//    printf("relative_to_reference is for particle %d:  %f %f %f\n", it, relative_to_reference[0],relative_to_reference[1],relative_to_reference[2]);
    for (int i=0; i<3; i++)
    {
    sum_of_distances[i] += relative_to_reference[i]+reference_position[i];
    }
//    printf("sum of distances %d:  %f %f %f\n", it, sum_of_distances[0], sum_of_distances[1], sum_of_distances[2]);

  }
  for (int i = 0; i < 3; i ++) {
    com[i] = fmod((sum_of_distances[i])*(1.0/particles.size()), box_l[i]);    //divide by number of particles in aggregate
      if (com[i] < 0.0) {
        com[i] = box_l[i]+com[i];
      if (com[i] > box_l[i]) {
        com[i] = com[i] - box_l[i];
      }
      }
   }
//  printf("**********************************************************\n");
//  printf("Cluster center of mass is: [%f,%f,%f].\n", com[0], com[1], com[2]);
//  printf("**********************************************************\n");
  return com;
}



//Longest distance
double Cluster::calculate_longest_distance()
{
  double relative_distance[3]={0.0, 0.0, 0.0};
  double partA[3], partB[3];
  
  double ld = 0.0; //the longest distance
  for (int it2a : particles) { //iterate over particles within an aggregate
//    printf ("it2 is: %d\n", it2);
//    printf ("particle it2 is: %d at the [%f %f %f]\n", local_particles[it2]->p.identity, local_particles[it2]->r.p[0], local_particles[it2]->r.p[1], local_particles[it2]->r.p[2]);
    for (int it2b : particles) { 
      for (int i=0; i!=3; ++i) {
        partA[i]=local_particles[it2a]->r.p[i];
        partB[i]=local_particles[it2b]->r.p[i];
      }
// caluclate vector  between particle A and B  
    get_mi_vector(relative_distance, partA, partB); //add current particle positions
// printf("relative Distance is: %f, %f, %f\n", relative_distance[0], relative_distance[1], relative_distance[2]);
       
    if (ld < (sqrt(sqrlen(relative_distance))))  //compare that distance with the longest distance
      ld=sqrt(sqrlen(relative_distance)); //save bigger value as longest distance - ld

//  printf("current distance and ld are: %f : %f\n", sqrt(sqrlen(relative_distance)), ld);
    } 
  }
//  printf("last ld is: %f\n", ld);
// printf("*****************************\n");
//  printf("The longest distance is: %f\n", ld);
//  printf("*****************************\n");

  return ld;
}


//Radius of gyration
double Cluster::calculate_radius_of_gyration()
{
  double distance[3]={0.0, 0.0, 0.0};
  double distance2 = 0.0;
  double rg2 = 0.0; 
  double rg =0.0;;

//  calculate com of the aggregate
  std::vector<double> com; //center of mass
// get vector center of mass
  com = calculate_cluster_center_of_mass();  
// get array center of mass
  double *comarray = &com[0];
  for (auto const& it3 : particles) {
// calculate distance between com and pid and store in variable current  
    get_mi_vector(distance, comarray, local_particles[it3]->r.p);
// calculate square length of this distance  
    distance2 += sqrlen(distance);
  }   
 
// divide with number of particles 
  rg2 = distance2/particles.size(); 
// return square root of it
  rg = sqrt(rg2);
//  printf("*****************************\n");
//  printf("The radius of gyration is: %f.\n", rg);
//  printf("*****************************\n");
  return rg;
}



// Fractal dimension
double Cluster::calculate_fractal_dimension()
{
  double df = 3.0;    //maximum df for spheres
  double relative_to_com[3]; //vector of the particle to the com
  double distance;    //distance of the particle from the center of the mass of the agglomerate
  int pid;            // particle ID
// calculate com of the aggregate
  std::vector<double> com; //center of mass
// double com; //center of mass
  com = calculate_cluster_center_of_mass();  
  double *comarray = &com[0];
  std::vector<double> distances; //list of distances 
  std::vector<double> distances_smaller_than_rad; //list of distances 
  std::vector<double> diameters; //all diameters=(radii*2) of circles around the com of aggregate
  std::vector<double> pcounts; // numbers of particles within given diameters
  
  int cluster_size = particles.size();

// calculate Df using linear regression on the logarithms of diameters [__std::vector<double> diameters__] and num of particles [__std::vector<int> pcounts__] within the diameters
  std::vector<double> log_diameters;
  std::vector<double> log_pcounts;

// calculate relative distance for each particle to the center of mass and store it into vector distances 
  for (auto const& it3 : particles) {
// calculate particle vector positions from the COM
    get_mi_vector(relative_to_com, comarray, local_particles[it3]->r.p); 
// calculate particle distance from the COM 
    distance = sqrlen(relative_to_com);
// printf("Particles distance is %f\n",distance );
    distances.push_back(sqrt(distance)); //add distance from the current particle to the com in the distances vectors
  }
  
  double rad = 0.0;
  int k = 0;
  int pcount = 0;  
// iterate over particles within an aggregate 
  while (k < particles.size()) 
  { 
    rad += 1;  //increase the radius for sigma=1
    for (int i = 0; i < distances.size(); i ++) {
     if (distances[i] < rad) {
      distances_smaller_than_rad.push_back(distances[i]);
// printf ("distances are %f\n", distances[i]);
     }
    }
    k = distances_smaller_than_rad.size(); 
    distances_smaller_than_rad.clear();
// printf("Number of particles within given rad %f:  %d\n", rad, k );
    if (k > 0) 
    {
      pcounts.push_back(k); //append number of particles wihin given diameter
      diameters.push_back(rad*2.0); //diameters are taken as doubled counter rad=0,1,2,3,..,particles.size()
        
      log_pcounts.push_back(log(k)); //append log of pcounts
      log_diameters.push_back(log(rad*2.0)); //append log of diameters
    }
}

#ifdef GSL
// usage: Function: int gsl_fit_linear (const double * x, const size_t xstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * sumsq) 
  int n=pcounts.size();
  double c0, c1, cov00, cov01, cov11, sumsq;
  if (diameters.size() > 1) 
  {
   gsl_fit_linear (&log_diameters.front(), 1, &log_pcounts.front(), 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);  
  }

#else
  runtimeErrorMsg()<< "GSL (gnu scientific library) is not found.";
#endif

 
//  printf("c0	c1	cov00	cov01	cov11	sumsq\n");
//  printf("%f	%f	%f	%f	%f	%f \n", c0, c1, cov00, cov01, cov11, sumsq);

  return (c1); 

}



int ClusterStructure::find_id_for(int x)
{
 int tmp=x;
 while (cluster_identities.find(tmp)!=cluster_identities.end())
 {
  tmp =cluster_identities[tmp];
 }
 return tmp;
}

int ClusterStructure::get_next_free_cluster_id(){
  //iterate over cluster_id's
  int max_seen_cluster = 0;
  for (auto it : cluster_id){
    int cid=it.second;
    if (max_seen_cluster < cid ) {
      max_seen_cluster=cid;
    }
  }
  return max_seen_cluster+1;
}


ClusterStructure cluster_structure;

ClusterStructure& cluster_analysis() {
  return cluster_structure;
}


//!!!!!!!!!!!!!!!!!!
// Geometry analysis performed on all aggregates
// com of all aggregates
std::vector<std::vector<double> > ClusterStructure::calculate_centers_of_masses()
{
  std::vector<std::vector<double> > coms;
  printf("COMS are:\n");
  for (auto it: clusters) {
    std::vector<double> single_com;
    single_com = it.second.calculate_cluster_center_of_mass();
    coms.push_back(single_com);
  // coms.push_back(it.second.calculate_cluster_center_of_mass());
  //  printf("[%f, %f, %f]\n", single_com[0], single_com[1], single_com[2]);
  printf("**********************************************************\n");
  printf("Cluster centers of masses is: [%f,%f,%f].\n", single_com[0], single_com[1], single_com[2]);
  printf("**********************************************************\n");
  }
  return coms;
}

 
// ld of all aggregates
std::vector<double> ClusterStructure::calculate_longest_distances() 
{
 std::vector<double> lds;

 double temp;
 for (auto  it: clusters) {
   temp = it.second.calculate_longest_distance();
   lds.push_back(temp);
 } 
return lds;
}


// rg of all aggregates
std::vector<double> ClusterStructure::calculate_radii_of_gyration() 
{
 std::vector<double> rgs;

 for (auto  it: clusters) {
   rgs.push_back(it.second.calculate_radius_of_gyration());
 } 
return rgs;
}


// df of all aggregates
std::vector<double> ClusterStructure::calculate_fractals_dimensions()
{
 std::vector<double> dfs;
 for (auto  it: clusters) {
   dfs.push_back(it.second.calculate_fractal_dimension());
 }
return dfs;
}

//!!!!!!!!!!!!!!!!!



void ClusterStructure::set_criterion(NeighborCriterion* c) {
  if (nc)
  {
    delete nc;
    nc=0;
  }
  nc=c;
}

   
