#ifndef SPATIAL_H
#define SPATIAL_H

#include "fourlevel.h"

std::vector<fourlevel_state> spatial_solve
(
    fourlevel_state, 
    boundary_conditions, // for z-grid
    filenames // for monitoring
);

#endif // SPATIAL_H
