#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/boundary_conditions.hpp"
#include "iostream"

namespace bc
{
    void boundary_conditions(t::MeshInfo& mesh_info, t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles, f::Fields& fields, sp::SimulationParameters& simulation_parameters)
    {
        for(int i = mesh_info.domain_triangles_number; i < mesh_info.triangles_number; ++i)
        {
            int type = ghost_triangles.types[i - mesh_info.domain_triangles_number];
            int inner_triangle_ID = ghost_triangles.inner_triangles_IDs[i - mesh_info.domain_triangles_number];
            if(type == 1)
            {
                fields.domain_triangles_centroids_velocities_before[i][0] = 2 * simulation_parameters.vx_inlet - fields.domain_triangles_centroids_velocities_before[inner_triangle_ID][0];
                fields.domain_triangles_centroids_velocities_before[i][1] = 2 * simulation_parameters.vy_inlet - fields.domain_triangles_centroids_velocities_before[inner_triangle_ID][1];
            }
            else if(type == 2)
            {
                fields.domain_triangles_centroids_velocities_before[i][0] = fields.domain_triangles_centroids_velocities_before[inner_triangle_ID][0];
                fields.domain_triangles_centroids_velocities_before[i][1] = fields.domain_triangles_centroids_velocities_before[inner_triangle_ID][1];
            }
            /* No-Slip condition (for walls)
            else if(type == 3)
            {
                fields.domain_triangles_centroids_velocities_before[i][0] =  - fields.domain_triangles_centroids_velocities_before[inner_triangle_ID][0];
                fields.domain_triangles_centroids_velocities_before[i][1] =  - fields.domain_triangles_centroids_velocities_before[inner_triangle_ID][1];
            }
            */
            //free-slip condition (for walls)
            else if(type == 3)
            {
                double nijx = 0.0, nijy = 0.0;
                auto [vix, viy] = fields.domain_triangles_centroids_velocities_before[inner_triangle_ID];
                for(int j = 0; j < 3; ++j)
                {
                    int neighbor_triangle_ID = domain_triangles.neighbor_triangles_IDs[inner_triangle_ID][j];
                    if(neighbor_triangle_ID == i)
                    {
                        nijx = domain_triangles.unit_normal_vectors[inner_triangle_ID][j][0];
                        nijy = domain_triangles.unit_normal_vectors[inner_triangle_ID][j][1];
                        break;
                    }
                }
                fields.domain_triangles_centroids_velocities_before[i][0] = vix - 2 * (vix * nijx + viy * nijy) * nijx;
                fields.domain_triangles_centroids_velocities_before[i][1] = viy - 2 * (vix * nijx + viy * nijy) * nijy;
            }
            
        }
        return;
    }
}
