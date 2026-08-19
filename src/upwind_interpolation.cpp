#include "include/upwind_interpolation.hpp"
#include "include/geometric_preprocess/types.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/simulation_parameters_reader.hpp"

namespace ui
{
    void upwind_interpolation(const t::MeshData &mesh, fds::Fields &fields, const sp::SimulationParameters& simulation_parameters)
    {
        int i, j, edge_id, type, neighbor_triangle_id, ghost_triangle_id;
        double average_velocity_flux, vxi, vyi, vxk, vyk, nxij, nyij, vx_average, vy_average;
        for(i=0; i<mesh.domain_triangles_number; ++i)
        {
            vxi = fields.triangles_centroids_vx_before[i];
            vyi = fields.triangles_centroids_vy_before[i];

            for(j=0; j<3; ++j)
            {
                edge_id = mesh.domain_triangles[i].edges_ids[j];
                type = mesh.edges[edge_id].type;
                neighbor_triangle_id = mesh.domain_triangles[i].neighbor_triangles_ids[j];
                if(type == 1)
                {
                    fields.edges_vx[edge_id] = simulation_parameters.vx_inlet;
                    fields.edges_vy[edge_id] = simulation_parameters.vy_inlet;
                }
                else if(type == 2)
                {
                    fields.edges_vx[edge_id] = - fields.triangles_centroids_vx_before[i];
                    fields.edges_vy[edge_id] = - fields.triangles_centroids_vy_before[i];
                }
                else if(type == 3)
                {
                    ghost_triangle_id = mesh.domain_triangles[i].neighbor_triangles_ids[j];
                    
                    fields.edges_vx[edge_id] = (fields.triangles_centroids_vx_before[i] + fields.triangles_centroids_vx_before[ghost_triangle_id]) / 2;
                    fields.edges_vy[edge_id] = (fields.triangles_centroids_vy_before[i] + fields.triangles_centroids_vy_before[ghost_triangle_id]) / 2;
                }
                else
                {
                    nxij = mesh.domain_triangles[i].normal_unitary_vectors[j][0];
                    nyij = mesh.domain_triangles[i].normal_unitary_vectors[j][1];
                
                    vxk = fields.triangles_centroids_vx_before[neighbor_triangle_id];
                    vyk = fields.triangles_centroids_vy_before[neighbor_triangle_id];

                    vx_average = (vxi + vxk) / 2;
                    vy_average = (vyi + vyk) / 2;

                    average_velocity_flux = (vx_average * nxij) + (vy_average * nyij);

                    if(average_velocity_flux >= 0)
                    {
                        fields.edges_vx[edge_id] = vxi;
                        fields.edges_vy[edge_id] = vyi;
                    }
                    else
                    {
                        fields.edges_vx[edge_id] = vxk;
                        fields.edges_vy[edge_id] = vyk;
                    }
                }
            }
        }
        return;
    }
}