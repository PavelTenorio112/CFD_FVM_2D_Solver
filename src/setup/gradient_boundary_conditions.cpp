#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/gradient_boundary_conditions.hpp"
#include <cmath>

namespace gbc
{
    void gradient_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, [[maybe_unused]] t::GhostTriangles& ghost_triangles, f::Fields& fields)
    {
        for(int i = mesh_info.domain_triangles_number; i < mesh_info.triangles_number; ++i)
        {
            int inner_triangle_ID = ghost_triangles.inner_triangles_IDs[i - mesh_info.domain_triangles_number],
                edge_ID = ghost_triangles.boundary_edges_IDs[i - mesh_info.domain_triangles_number];
            double rijx = 0.0, rijy = 0.0;

            auto[Dxvxi, Dyvxi] = fields.triangles_centroids_gradient[inner_triangle_ID][0];
            auto[Dxvyi, Dyvyi] = fields.triangles_centroids_gradient[inner_triangle_ID][1];
            auto[vix, viy] = fields.triangles_centroids_velocities_before[inner_triangle_ID];
            auto[vkx, vky] = fields.triangles_centroids_velocities_before[i];
            


            for(int j = 0; j < 3; ++j)
            {
                if(domain_triangles.edges_IDs[inner_triangle_ID][j] == edge_ID)
                {
                    rijx = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[inner_triangle_ID][j][0];
                    rijy = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[inner_triangle_ID][j][1];
                    break;
                }
            }

            double rij = std::sqrt(rijx * rijx + rijy * rijy),
                   unit_rijx = rijx / rij,
                   unit_rijy = rijy / rij,
                   
                   Dxvxk = (unit_rijx / rij) * (vkx - vix) - Dxvxi,
                   Dyvxk = (unit_rijy / rij) * (vkx - vix) - Dyvxi,
                   Dxvyk = (unit_rijx / rij) * (vky - viy) - Dxvyi,
                   Dyvyk = (unit_rijy / rij) * (vky - viy) - Dyvyi;
                   


            fields.triangles_centroids_gradient[i][0][0] = Dxvxk;
            fields.triangles_centroids_gradient[i][0][1] = Dyvxk;
            fields.triangles_centroids_gradient[i][1][0] = Dxvyk;
            fields.triangles_centroids_gradient[i][1][1] = Dyvyk;
        }
        return;
    }
}
