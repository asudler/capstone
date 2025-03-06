#ifndef SPATIAL_H
#define SPATIAL_H

#include "/home/asudler/git/capstone/fourlevel/fourlevel.h"

std::vector<fourlevel_state> spatial_solve
(
    fourlevel_state, 
    boundary_conditions, // for z-grid
    filenames // for monitoring
);

// the next two functions need to be implemented
void spatial_print
(
    std::vector<fourlevel_state>,
    boundary_conditions,
    filenames
);
void spatial_print_lite
(
    std::vector<fourlevel_state>,
    boundary_conditions,
    filenames
);
std::vector<fourlevel_state> spatial_read
(
    std::string,
    boundary_conditions,
    filenames
);

void rotate(fourlevel_state state, filenames files);
void unrotate(fourlevel_state state, filenames files);

#endif // SPATIAL_H
