#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/boundary_conditions.hpp"
#include "iostream"

namespace bc
{
    void boundary_conditions(t::MeshData& mesh, fds::Fields& fields, sp::SimulationParameters& simulation_parameters)
    {
        int i, j, type, inner_triangle_ID, neighbor_triangle_id;
        double vix = 0.0, viy = 0.0, nijx = 0.0, nijy = 0.0;
        for(i = mesh.domain_triangles_number; i < mesh.triangles_number; ++i)
        {
            type = mesh.ghost_triangles[i - mesh.domain_triangles_number].type;
            inner_triangle_ID = mesh.ghost_triangles[i - mesh.domain_triangles_number].interior_triangle_id;

            if(type == 1)
            {
                fields.triangles_centroids_vx_before[i] = 2 * simulation_parameters.vx_inlet - fields.triangles_centroids_vx_before[inner_triangle_ID];
                fields.triangles_centroids_vy_before[i] = 2 * simulation_parameters.vy_inlet - fields.triangles_centroids_vy_before[inner_triangle_ID];
            }
            else if(type == 2)
            {
                fields.triangles_centroids_vx_before[i] = fields.triangles_centroids_vx_before[inner_triangle_ID];
                fields.triangles_centroids_vy_before[i] = fields.triangles_centroids_vy_before[inner_triangle_ID];
            }
            /* No-Slip condition (for walls)
            else if(type == 3)
            {
                fields.triangles_centroids_vx_before[i] = - fields.triangles_centroids_vx_before[inner_triangle_ID];
                fields.triangles_centroids_vy_before[i] = - fields.triangles_centroids_vy_before[inner_triangle_ID];
            }
            */
            else if(type == 3) // Free-Slip condition (for walls)
            {
                vix = fields.triangles_centroids_vx_before[inner_triangle_ID];
                viy = fields.triangles_centroids_vy_before[inner_triangle_ID];
                for(j = 0; j < 3; ++j)
                {
                    neighbor_triangle_id = mesh.domain_triangles[inner_triangle_ID].neighbor_triangles_ids[j];
                    if(neighbor_triangle_id == i)
                    {
                        nijx = mesh.domain_triangles[inner_triangle_ID].normal_unitary_vectors[j][0];
                        nijy = mesh.domain_triangles[inner_triangle_ID].normal_unitary_vectors[j][1];

                        break;
                    }
                }
                fields.triangles_centroids_vx_before[i] = vix - 2 * (vix * nijx + viy * nijy) * nijx;
                fields.triangles_centroids_vy_before[i] = viy - 2 * (vix * nijx + viy * nijy) * nijy;
            }
        }
        return;
    }
}
