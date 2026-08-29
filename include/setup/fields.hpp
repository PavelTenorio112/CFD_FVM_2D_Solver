#pragma once
#include<vector>
#include<array>
namespace f
{
    struct Fields
    {
        std::vector <std::array<double, 2> > domain_triangles_centroids_velocities_before;
        std::vector <std::array<double, 2> > domain_triangles_centroids_velocities_after;
        std::vector <std::array<double, 2> >  edges_midpoints_velocities;
    };
}