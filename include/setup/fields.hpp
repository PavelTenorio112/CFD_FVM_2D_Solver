#pragma once
#include<vector>
#include<array>
namespace f
{
    struct Fields
    {
        std::vector <std::array<double, 2> > triangles_centroids_velocities_before;
        std::vector <std::array<double, 2> > triangles_centroids_velocities_after;
        std::vector <std::array <std::array<double, 2>, 2> > triangles_centroids_gradient;
        
        std::vector <std::array<double, 2> > edges_midpoints_velocities;
        std::vector <std::array <std::array<double, 2>, 2> > edges_midpoints_gradient;

    };
}