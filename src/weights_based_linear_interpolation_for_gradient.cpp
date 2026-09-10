#include"include/weights_based_linear_interpolation_for_gradient.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include <iostream>
namespace wblifg
{
    void weights_based_linear_interpolation_for_gradients(t::MeshInfo& mesh_info, t::Edges& edges, f::Fields &fields)
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
}