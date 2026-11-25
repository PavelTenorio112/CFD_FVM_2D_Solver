#include "interpolation_schemes.hpp"
#include "geometric_preprocess/types.hpp"
#include "setup/fields.hpp"
#include <cmath>
#include <algorithm>

namespace is
{
    void first_order_upwind_interpolation(t::MeshInfo &mesh_info, t::Edges& edges, [[maybe_unused]] t::DomainTriangles& domain_triangles, f::Fields &fields)
    {
        for(int j = 0; j < mesh_info.edges_number; ++j)
        {
            auto[triangle_0_ID, triangle_1_ID] = edges.triangles_IDs[j];
            auto[v0x, v0y] = fields.triangles_centroids_velocities_before[triangle_0_ID];
            auto[v1x, v1y] = fields.triangles_centroids_velocities_before[triangle_1_ID];
            double vx_average = (v0x + v1x) / 2;
            double vy_average = (v0y + v1y) / 2;
            auto[nx, ny] = edges.unit_normal_vectors[j];
            double flux = vx_average * nx + vy_average * ny;
            if(flux >= 0)
            {
                fields.edges_midpoints_velocities[j][0] = v0x;
                fields.edges_midpoints_velocities[j][1] = v0y;
            }
            else
            {
                fields.edges_midpoints_velocities[j][0] = v1x;
                fields.edges_midpoints_velocities[j][1] = v1y;
            }
        }
        return;
    }
    void weights_based_linear_interpolation_for_velocity_gradient(t::MeshInfo& mesh_info, t::Edges& edges, f::Fields &fields)
    {   
        for(int i = 0; i < mesh_info.edges_number; ++i)
        {
            auto[triangle_0_ID , triangle_1_ID] = edges.triangles_IDs[i];
            auto[geometric_weight_0, geometric_weight_1] = edges.geometric_weights[i];
            
            fields.edges_midpoints_gradient[i][0][0] = fields.triangles_centroids_gradient[triangle_0_ID][0][0] * geometric_weight_0 + fields.triangles_centroids_gradient[triangle_1_ID][0][0] * geometric_weight_1; 
            fields.edges_midpoints_gradient[i][0][1] = fields.triangles_centroids_gradient[triangle_0_ID][0][1] * geometric_weight_0 + fields.triangles_centroids_gradient[triangle_1_ID][0][1] * geometric_weight_1; 
            fields.edges_midpoints_gradient[i][1][0] = fields.triangles_centroids_gradient[triangle_0_ID][1][0] * geometric_weight_0 + fields.triangles_centroids_gradient[triangle_1_ID][1][0] * geometric_weight_1; 
            fields.edges_midpoints_gradient[i][1][1] = fields.triangles_centroids_gradient[triangle_0_ID][1][1] * geometric_weight_0 + fields.triangles_centroids_gradient[triangle_1_ID][1][1] * geometric_weight_1; 
        }
        return;
    }
 
    void orthogonality_correction_for_velocity_gradient_interpolation(t::MeshInfo& mesh_info, t::DomainTriangles& domain_triangles, f::Fields& fields)
    {
        for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
        {
            auto[vix, viy] = fields.triangles_centroids_velocities_before[i];
            for(int j = 0; j < 3; ++j)
            {
                auto[rikx, riky] = domain_triangles.from_self_centroid_to_neighbor_triangles_centroids_vectors[i][j];
                double rik = hypot(rikx, riky);
                double unit_rikx = rikx / rik;
                double unit_riky = riky / rik;
                int neighbor_triangle_ID = domain_triangles.neighbor_triangles_IDs[i][j];
                int edge_ID = domain_triangles.edges_IDs[i][j];
                auto[vkx, vky] = fields.triangles_centroids_velocities_before[neighbor_triangle_ID];
                auto[Dxvjx, Dyvjx] = fields.edges_midpoints_gradient[edge_ID][0];
                auto[Dxvjy, Dyvjy] = fields.edges_midpoints_gradient[edge_ID][1];
                fields.edges_midpoints_gradient[edge_ID][0][0] = ((vkx - vix) / rik) * unit_rikx + Dxvjx - (Dxvjx * unit_rikx + Dyvjx * unit_riky) * unit_rikx;
                fields.edges_midpoints_gradient[edge_ID][0][1] = ((vkx - vix) / rik) * unit_riky + Dyvjx - (Dxvjx * unit_rikx + Dyvjx * unit_riky) * unit_riky;
                fields.edges_midpoints_gradient[edge_ID][1][0] = ((vky - viy) / rik) * unit_rikx + Dxvjy - (Dxvjy * unit_rikx + Dyvjy * unit_riky) * unit_rikx;
                fields.edges_midpoints_gradient[edge_ID][1][1] = ((vky - viy) / rik) * unit_riky + Dyvjy - (Dxvjy * unit_rikx + Dyvjy * unit_riky) * unit_riky;
            }
        }
        return;
    }
}
