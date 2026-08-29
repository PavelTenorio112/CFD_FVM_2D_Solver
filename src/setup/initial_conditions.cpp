#include"include/setup/initial_conditions.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"

namespace ic
{
    void initial_conditions(t::MeshInfo &mesh_info, f::Fields &fields)
    {
        fields.domain_triangles_centroids_velocities_before.resize(mesh_info.triangles_number);
        fields.domain_triangles_centroids_velocities_after.resize(mesh_info.triangles_number);
        
        fields.edges_midpoints_velocities.resize(mesh_info.edges_number);

        //Quiescent initial condition
        for(int i = 0; i < mesh_info.domain_triangles_number; ++i) 
        {
            fields.domain_triangles_centroids_velocities_before[i][0] = 0.0;
            fields.domain_triangles_centroids_velocities_before[i][1] = 0.0;
        }
        return;
    }
}

