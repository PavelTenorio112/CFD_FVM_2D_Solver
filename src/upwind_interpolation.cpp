#include "include/upwind_interpolation.hpp"
#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/simulation_parameters_reader.hpp"

namespace ui
{
    void upwind_interpolation(t::MeshInfo &mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles, f::Fields &fields)
    {
        for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
        {
            double vix = fields.triangles_centroids_velocities_before[i][0],
                   viy = fields.triangles_centroids_velocities_before[i][1];

            for(int j = 0; j < 3; ++j)
            {
                int edge_ID = domain_triangles.edges_IDs[i][j],
                    neighbor_triangle_ID = domain_triangles.neighbor_triangles_IDs[i][j],
                    type = edges.types[edge_ID];

                auto[vkx, vky] = fields.triangles_centroids_velocities_before[neighbor_triangle_ID];

                double vx_average = (vix + vkx) / 2,
                       vy_average = (viy + vky) / 2;

                if(type == 1 or type == 2 or type == 3)
                {
                    fields.edges_midpoints_velocities[edge_ID][0] = vx_average;
                    fields.edges_midpoints_velocities[edge_ID][1] = vy_average;
                }
                else
                {
                    double nijx = domain_triangles.unit_normal_vectors[i][j][0],
                           nijy = domain_triangles.unit_normal_vectors[i][j][1],
                           average_velocity_flux = (vx_average * nijx) + (vy_average * nijy);

                    if(average_velocity_flux >= 0)
                    {
                        fields.edges_midpoints_velocities[edge_ID][0] = vix;
                        fields.edges_midpoints_velocities[edge_ID][1] = viy;
                    }
                    else
                    {
                        fields.edges_midpoints_velocities[edge_ID][0] = vkx;
                        fields.edges_midpoints_velocities[edge_ID][1] = vky;
                    }
                }

            }
        }
        return;
    }
}