#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/velocity_boundary_conditions.hpp"
#include "iostream"

namespace vbc
{
    void velocity_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles, f::Fields& fields, sp::SimulationParameters& simulation_parameters)
    {
        for(int i = mesh_info.domain_triangles_number; i < mesh_info.triangles_number; ++i)
        {
            int type = ghost_triangles.types[i - mesh_info.domain_triangles_number],
                inner_triangle_ID = ghost_triangles.inner_triangles_IDs[i - mesh_info.domain_triangles_number];
            double vx_inlet = simulation_parameters.vx_inlet,
                   vy_inlet = simulation_parameters.vy_inlet,
                   ghost_triangle_vx = 0.0, 
                   ghost_triangle_vy = 0.0;
            auto [inner_triangle_vx, inner_triangle_vy] = fields.triangles_centroids_velocities_before[inner_triangle_ID];
            


            if(type == 1) // Inlet
            {
                ghost_triangle_vx = 2 * vx_inlet - inner_triangle_vx;
                ghost_triangle_vy = 2 * vy_inlet - inner_triangle_vy;
            }
            else if(type == 2) // Outlet
            {
                ghost_triangle_vx = inner_triangle_vx;
                ghost_triangle_vy = inner_triangle_vy;
            }
            else if(type == 3) // Walls
            {
                // no-slip condition
                ghost_triangle_vx = - inner_triangle_vx;
                ghost_triangle_vy = - inner_triangle_vy;

                /*// free- slip condition
                double nijx = 0.0, nijy = 0.0;
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
                ghost_triangle_vx = inner_triangle_vx - 2 * (inner_triangle_vx * nijx + inner_triangle_vy * nijy) * nijx;
                ghost_triangle_vy = inner_triangle_vy - 2 * (inner_triangle_vx * nijx + inner_triangle_vy * nijy) * nijy;
                */
            }
            fields.triangles_centroids_velocities_before[i][0] = ghost_triangle_vx;
            fields.triangles_centroids_velocities_before[i][1] = ghost_triangle_vy;
        }
        return;
    }
}
