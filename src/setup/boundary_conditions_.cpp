#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/boundary_conditions.hpp"

namespace bc
{
    void boundary_conditions(t::MeshData &mesh, fds::Fields &fields, sp::SimulationParameters &simulation_parameters)
    {
        int i, type, inner_triangle_ID;;
        for(i=mesh.domain_triangles_number; i<mesh.domain_triangles_number + mesh.ghost_triangles_number; ++i)
        {
            type = mesh.ghost_triangles[i].type;
            inner_triangle_ID = mesh.ghost_triangles[i].interior_triangle_id;

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
            else if(type == 3)
            {
                fields.triangles_centroids_vx_before[i] = - fields.triangles_centroids_vx_before[inner_triangle_ID];
                fields.triangles_centroids_vy_before[i] = - fields.triangles_centroids_vy_before[inner_triangle_ID];
            }
        }
        return;
    }
}
