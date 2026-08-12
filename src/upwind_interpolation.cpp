#include "include/upwind_interpolation.hpp"
#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"

namespace ui
{
    void upwind_interpolation(t::MeshData &mesh, fds::Fields &fields)
    {
        int i, j, triangle_0_ID, triangle_1_ID, edge_ID;
        double self_flux, neighbor_flux, average, vxi, vyi, vxk, vyk, nxij, nyij;
        for(i=0; i<mesh.domain_triangles_number; ++i)
        {
            vxi = fields.triangles_centroids_vx_before[i];
            vyi = fields.triangles_centroids_vy_before[i];

            for(j=0; j<3; ++j)
            {
                edge_ID = mesh.domain_triangles[i].edges_ids[j];

                vxk = fields.triangles_centroids_vx_before[mesh.domain_triangles[i].neighbor_triangles_ids[j]];
                vyk = fields.triangles_centroids_vy_before[mesh.domain_triangles[i].neighbor_triangles_ids[j]];

                nxij = mesh.domain_triangles[i].normal_unitary_vectors[j][0];
                nyij = mesh.domain_triangles[i].normal_unitary_vectors[j][1];

                self_flux = vxi * nxij + vyi * nyij;
                neighbor_flux = vxk * nxij + vyk * nyij;
                average = (self_flux + neighbor_flux) / 2;
                if(average > 0)
                {
                    fields.edges_vx[edge_ID] = vxi;
                    fields.edges_vy[edge_ID] = vyi;
                }
                else
                {
                    fields.edges_vx[edge_ID] = vxk;
                    fields.edges_vy[edge_ID] = vyk;
                }
            }
        }
        return;
    }
}