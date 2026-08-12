#pragma once
#include"inviscid_burgers.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/boundary_conditions.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
#include"include/simulation_parameters.hpp"
#include"include/upwind_interpolation.hpp"

namespace ib
{
    void inviscid_burguers(t::MeshData &mesh, fds::Fields &fields, sp::SimulationParameters &simulation_parameters)
    {
        int n, i, j, N = simulation_parameters.T / simulation_parameters.dt;
        double sumx = 0, sumy = 0, vxijn, vyijn, nxij, nyij, area;
        for(n=0; n<N; ++n)
        {
            bc::boundary_conditions(mesh, fields, simulation_parameters);
            ui::upwind_interpolation(mesh, fields);
            for(i=0; i<mesh.domain_triangles_number; ++i)
            {
                sumx = 0;
                sumy = 0;
                for(j=0; j<3; ++j)
                {
                    vxijn = fields.edges_vx[mesh.domain_triangles[i].edges_ids[j]];
                    vxijn = fields.edges_vy[mesh.domain_triangles[i].edges_ids[j]];
                    nxij = mesh.domain_triangles[i].normal_unitary_vectors[j][0];
                    nyij = mesh.domain_triangles[i].normal_unitary_vectors[j][1];
                    
                    sumx = sumx - (vxijn * nxij + vyijn * nyij) * vxijn;
                    sumy = sumy - (vxijn * nxij + vyijn * nyij) * vyijn;
                }
                fields.triangles_centroids_vx_after[i] = fields.triangles_centroids_vx_before[i] + (mesh.domain_triangles[i].area / simulation_parameters.dt) * sumx;
                fields.triangles_centroids_vy_after[i] = fields.triangles_centroids_vy_before[i] + (mesh.domain_triangles[i].area / simulation_parameters.dt) * sumy; 
            }
            std::swap(fields.triangles_centroids_vx_before, fields.triangles_centroids_vx_after);
            std::swap(fields.triangles_centroids_vy_before, fields.triangles_centroids_vy_after);
        }
        return;
    }
}