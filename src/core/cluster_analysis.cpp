#include "cluster_analysis.hpp"
#include "interaction_data.hpp"
#include <algorithm>
#include "config.hpp"

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
  calculate_cluster_center_of_mass();
  //calculate_fractal_dimension();
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
std::vector<double>  Cluster::calculate_cluster_center_of_mass() 
{
  double zeros[] = {0.0,0.0,0.0};
  std::vector<double> com; //initialized com
  for (int i=0; i<3; i++)
    com.push_back(0.0);   
   
  for (auto const& it : particles)  //iterate over all particles within a cluster
  {
    for (int i=0; i<3; i++)
    { 
      com[i] += local_particles[it]->r.p[i]; //add current particle positions
    }
  }
  for (int i=0; i<3; i++) {
    com[i] = com[i]*(1.0/particles.size()); //divide by number of particles in aggregate
  }
  return com;
}


double Cluster::calculate_fractal_dimension()
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
//usage: Function: int gsl_fit_linear (const double * x, const size_t xstride, const double * y, const size_t ystride, size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * sumsq) 
  df=3.0;
  double c0, c1, cov00, cov01, cov11, sumsq;
  if (diameters.size() > 1) : 
    gsl_fit_linear(x,1,y,1,c0, c1, cov00, cov01, cov11, sumsq);
    gsl_fit_linear (&diameters.front(), 1, &pcounts.front(), 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);  
  return df; 

#else
  runtimeErrorMsg()<< "GSL (gnu scientific library) is not found.";
#endif
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

void ClusterStructure::set_criterion(NeighborCriterion* c) {
  if (nc)
  {
    delete nc;
    nc=0;
  }
  nc=c;
}

   

 

