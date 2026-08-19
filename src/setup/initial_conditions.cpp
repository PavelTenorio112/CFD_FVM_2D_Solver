#include"include/setup/initial_conditions.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"

namespace ic
{
    void initial_conditions(t::MeshData &mesh, fds::Fields &fields)
    {
        int i;
        fields.triangles_centroids_vx_before.reserve(mesh.triangles_number);
        fields.triangles_centroids_vx_after.reserve(mesh.triangles_number);
        fields.triangles_centroids_vy_before.reserve(mesh.triangles_number);
        fields.triangles_centroids_vy_after.reserve(mesh.triangles_number);

        fields.edges_vx.reserve(mesh.edges_number);
        fields.edges_vy.reserve(mesh.edges_number);

        for(i=0; i<mesh.domain_triangles_number; ++i) 
        {
            fields.triangles_centroids_vx_before[i] = 0;
            fields.triangles_centroids_vy_before[i] = 0;
        }
        return;
    }
}

