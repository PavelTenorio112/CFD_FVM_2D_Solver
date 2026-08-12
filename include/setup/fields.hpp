#pragma once
#include"include/geometric_preprocess/types.hpp"

namespace fds
{
    struct Fields
    {
        std::vector <double> triangles_centroids_vx_before;
        std::vector <double> triangles_centroids_vy_before;
        std::vector <double> triangles_centroids_vx_after;
        std::vector <double> triangles_centroids_vy_after;

        std::vector <double> edges_vx;
        std::vector <double> edges_vy;
    };
}