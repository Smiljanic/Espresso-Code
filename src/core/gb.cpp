/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file gb.cpp
 *
 *  Implementation of \ref gb.hpp
 */
#include "gb.hpp"
#include "communication.hpp"

#ifdef GAY_BERNE

int gay_berne_set_params(int part_type_a, int part_type_b,
			 double eps1, double eps2, double eps3, 
                         double sig1, double sig2, double sig3, 
                         double cut,
			 double mu, double nu)
{
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data) return ES_ERROR;

  data->GB_eps1    = eps1;
  data->GB_eps2    = eps2;
  data->GB_eps3    = eps3;
  data->GB_sig1    = sig1;
  data->GB_sig2    = sig2;
  data->GB_sig3    = sig3;
  data->GB_cut    = cut;
  data->GB_mu     = mu;
  data->GB_nu     = nu;
 
  /* Calculate dependent parameters */

  data->GB_chi1 = ((data->GB_sig1*data->GB_sig2) - 1) / ((data->GB_sig1*data->GB_sig3) + 1);
  data->GB_chi2 = (pow(data->GB_sig1,(1/data->GB_mu))-1)/(pow(data->GB_sig3,(1/data->GB_mu))+1);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif
