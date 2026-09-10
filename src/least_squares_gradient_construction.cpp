#include"include/least_squares_gradient_construction.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include<iostream>
namespace lsgc
{
    void least_squares_gradient_construction(t::MeshInfo& mesh_info, t::DomainTriangles& domain_triangles, f::Fields& fields)
    {
        for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
        {
            auto[vix, viy] = fields.triangles_centroids_velocities_before[i];

            double sum_dy2 = 0.0,
                   sum_dx2 = 0.0, 
                   sum_dxdy = 0.0, 
                   sum_dxdvx = 0.0, 
                   sum_dxdvy = 0.0, 
                   sum_dydvx = 0.0, 
                   sum_dydvy = 0.0;

            for(int j = 0; j < 3; ++j)
            {
                int neighbor_triangle_ID = domain_triangles.neighbor_triangles_IDs[i][j];
                auto[vkx, vky] = fields.triangles_centroids_velocities_before[neighbor_triangle_ID];
                auto[dx, dy] = domain_triangles.from_self_centroid_to_neighbor_triangles_centroids_vectors[i][j];
                double dvx = vkx - vix,
                       dvy = vky - viy;

                sum_dx2 = sum_dx2 + dx * dx;
                sum_dy2 = sum_dy2 + dy * dy;
                sum_dxdy = sum_dxdy + dx * dy;
                sum_dxdvx = sum_dxdvx + dvx * dx;
                sum_dydvx = sum_dydvx + dvx * dy;
                sum_dxdvy = sum_dxdvy + dvy * dx;
                sum_dydvy = sum_dydvy + dvy * dy;
            }
            double denominator = sum_dx2 * sum_dy2 - sum_dxdy * sum_dxdy;

            fields.triangles_centroids_gradient[i][0][0] = (sum_dy2 * sum_dxdvx - sum_dxdy * sum_dydvx) / denominator;
            fields.triangles_centroids_gradient[i][0][1] = (sum_dx2 * sum_dydvx - sum_dxdy * sum_dxdvx) / denominator;
            fields.triangles_centroids_gradient[i][1][0] = (sum_dy2 * sum_dxdvy - sum_dxdy * sum_dydvy) / denominator; 
            fields.triangles_centroids_gradient[i][1][1] = (sum_dx2 * sum_dydvy - sum_dxdy * sum_dxdvy) / denominator; 
        }
        return;
    }
}



