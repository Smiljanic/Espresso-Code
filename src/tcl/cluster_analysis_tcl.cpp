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
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "cluster_analysis.hpp"
#include "parser.hpp"
#include <sstream>
#include <cstdio>
int tclcommand_cluster_analysis(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  // If no argumens are given, print status
  if (argc==1) {
      if (cluster_analysis().get_criterion()) {
        Tcl_AppendResult(interp, cluster_analysis().get_criterion()->name().c_str(), (char*) NULL);
      } else {
        Tcl_AppendResult(interp, "off", (char*) NULL);
      }        
      return TCL_OK;
    }

  argc--; argv++;

  // Otherwise, we set parameters
  if (ARG0_IS_S("off")) {
    cluster_analysis().set_criterion(NULL);
    return TCL_OK;
  }
  if (ARG0_IS_S("distance")) {
      if (argc != 2) {
      	Tcl_AppendResult(interp, "The distnace criterion needs a distance as argument.", (char*) NULL);
      	return TCL_ERROR;
      }
      double d;
      if (!ARG_IS_D(1,d)) {
        	Tcl_AppendResult(interp, "Need a distance as 1st arg.", (char*) NULL);
        	return TCL_ERROR;
      }
      cluster_analysis().set_criterion(new DistanceCriterion(d));
      argc -= 2; argv += 2;
    }
  else if (ARG0_IS_S("bond")) {
      if (argc != 2) {
      	Tcl_AppendResult(interp, "The bond criterion needs a bond type as argument.", (char*) NULL);
      	return TCL_ERROR;
      }
      int b;
      if (!ARG_IS_I(1,b)) {
        	Tcl_AppendResult(interp, "Need a bond type as 1st arg.", (char*) NULL);
        	return TCL_ERROR;
      }
      cluster_analysis().set_criterion(new BondCriterion(b));
      argc -= 2; argv += 2;
    }
    else if (ARG0_IS_S("analyze_pair")) {
      cluster_analysis().analyze_pair();
      argc -= 1; argv += 1;
    }
    else if (ARG0_IS_S("analyze_bond")) {
      cluster_analysis().analyze_bonds();
      argc -= 1; argv += 1;
    }

    else if (ARG0_IS_S("print")) {
      std::stringstream res;
      for (auto it : cluster_analysis().clusters) {
        res << "{"<<it.first<<" {";
        Cluster cluster = it.second;
        for (int pid : cluster.particles) {
          res << pid<<" " << "[" << local_particles[pid]->r.p[0] << "," << local_particles[pid]->r.p[1] << "," << local_particles[pid]->r.p[2] << "] ";
        }
        res << "} } ";
      }
      argc -= 1; argv += 1;
    	Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }

//print particles within given clusters

    else if (ARG0_IS_S("particles-within-clusters")) {
      std::stringstream res;
      res << "#N	cluster_ID	particle_ID		pos[x]		pos[y]		pos[z]\n";
      int N = 0;
      for (auto it : cluster_analysis().clusters) {
//        res << N <<"	";
        Cluster cluster = it.second;
        N  = cluster.size_of_cluster();
        for (int pid : cluster.particles) {
         // if (!(local_particles[pid]->p.isVirtual)) {
            res << N << "	" << it.first << "	" << pid << "	" <<  local_particles[pid]->r.p[0] << "	" << local_particles[pid]->r.p[1] << "	" << local_particles[pid]->r.p[2] <<"\n";
          }
        //}
      }
      argc -= 1; argv += 1;
    	Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }



//postprocessing  properties

    else if (ARG0_IS_S("geometry-all")) {
//      std::cout << "\nClusters geometry analysis: \n";
      std::stringstream res;
      int N = 0;
      double ld, rg, df;
      res << "#N		com[x]		com[y]		com[z]		ld		rg		df\n";
      for (auto it : cluster_analysis().clusters) {
//        res << it.first << " {";
        Cluster& cluster = it.second;
        std::vector<double> com = cluster.calculate_cluster_center_of_mass();
        ld = cluster.calculate_longest_distance();
        rg = cluster.calculate_radius_of_gyration();
        df = cluster.calculate_fractal_dimension();
//      N = sizeof(cluster); //returns 24 as size of an int=8 , pointer to begin of vector, to end, and end of reserved memory
        N  = cluster.size_of_cluster();
        res << N << "		" << com[0]  << "		" << com[1] << "		" << com[2]  << "		" << ld << "		"  <<  rg  <<"		" << df <<"\n";
       }

      argc -= 1; argv += 1;
        Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }



 ///coms
     else if (ARG0_IS_S("com")) {
      std::cout << "\nCluster center of mass: \n";
      std::stringstream res;
      for (auto it : cluster_analysis().clusters) {
        res << it.first << " {";
        Cluster& cluster = it.second;
        std::vector<double> com = cluster.calculate_cluster_center_of_mass();
        for (int pid : cluster.particles) 
          res << pid<<" " << "[" << local_particles[pid]->r.p[0] << "," << local_particles[pid]->r.p[1] << "," << local_particles[pid]->r.p[2] << "] ";
      //  res << com[0] << "," << com[1] << "," << com[2] << "}\n";
      }
      argc -= 1; argv += 1;
        Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }


  

    else if (ARG0_IS_S("ld")) {
      std::cout << "\nCluster longest distance: \n";
      std::stringstream res;
      for (auto it : cluster_analysis().clusters) {
        res << it.first << " {";
        Cluster& cluster = it.second;
        for (int pid : cluster.particles) 
          res << pid<<" " << "[" << local_particles[pid]->r.p[0] << "," << local_particles[pid]->r.p[1] << "," << local_particles[pid]->r.p[2] << "] ";
        double ld = cluster.calculate_longest_distance();
        for (int pid : cluster.particles) 
          res << pid<<" " << "[" << local_particles[pid]->r.p[0] << "," << local_particles[pid]->r.p[1] << "," << local_particles[pid]->r.p[2] << "] ";
        res << ld << "}\n";
      }
      argc -= 1; argv += 1;
        Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }


    else if (ARG0_IS_S("rg")) {
      std::cout << "Cluster radius of gyration: \n";
      std::stringstream res;
      for (auto it : cluster_analysis().clusters) {
        res << it.first << " {";
        Cluster& cluster = it.second;
        double rg;
        rg = cluster.calculate_radius_of_gyration();
        res << rg << "}\n";
      }
      argc -= 1; argv += 1;
    	Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }




    else if (ARG0_IS_S("df")) {
      std::cout << "Cluster fractal dimensions: \n";
      std::stringstream res;
      for (auto it : cluster_analysis().clusters) {
        res << it.first << " {";
        Cluster& cluster = it.second;
        double df;
        df = cluster.calculate_fractal_dimension();
        res << df << "}\n";
      }
      argc -= 1; argv += 1;
    	Tcl_AppendResult(interp, res.str().c_str(), (char*) NULL);
      return TCL_OK;
    }


    else {
    	Tcl_AppendResult(interp, "Unknown argument.", (char*) NULL);
	    return TCL_ERROR;
    }
  return gather_runtime_errors(interp,TCL_OK);
}
