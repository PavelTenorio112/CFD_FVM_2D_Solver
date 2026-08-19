#include"include/inviscid_burgers.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/boundary_conditions.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
#include"include/simulation_parameters.hpp"
#include"include/upwind_interpolation.hpp"
#include"include/fields_data_print.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace ib
{
    void inviscid_burguers(t::MeshData &mesh, fds::Fields &fields, sp::SimulationParameters &simulation_parameters)
    {
        int time_step_number, time_steps_number, i, j, edge_id, frame_number;
        double sumx, sumy, vxijn, vyijn, nxij, nyij, time_accumulator, image_dt;
        time_steps_number = simulation_parameters.T / simulation_parameters.dt;
        image_dt = (1.0 / simulation_parameters.frames_per_second);
        time_accumulator = 0;
        frame_number = 0;
        for(time_step_number = 0; time_step_number < time_steps_number; ++time_step_number)
        {
            time_accumulator = time_accumulator + simulation_parameters.dt;
            if(time_accumulator > image_dt)
            {
                fdp::fields_data_print(mesh, fields, frame_number);
                time_accumulator = time_accumulator - image_dt;
                ++frame_number;
            }

            bc::boundary_conditions(mesh, fields, simulation_parameters);
            ui::upwind_interpolation(mesh, fields, simulation_parameters);
            for(i=0; i<mesh.domain_triangles_number; ++i)
            {
                sumx = 0;
                sumy = 0;
                for(j=0; j<3; ++j)
                {
                    vxijn = fields.edges_vx[mesh.domain_triangles[i].edges_ids[j]];
                    vyijn = fields.edges_vy[mesh.domain_triangles[i].edges_ids[j]];
                    nxij = mesh.domain_triangles[i].normal_unitary_vectors[j][0];
                    nyij = mesh.domain_triangles[i].normal_unitary_vectors[j][1];

                    edge_id = mesh.domain_triangles[i].edges_ids[j];

                    sumx = sumx - (vxijn * nxij + vyijn * nyij) * vxijn * mesh.edges[edge_id].length; 
                    sumy = sumy - (vxijn * nxij + vyijn * nyij) * vyijn * mesh.edges[edge_id].length;
                }
                fields.triangles_centroids_vx_after[i] = fields.triangles_centroids_vx_before[i] + (simulation_parameters.dt / mesh.domain_triangles[i].area) * sumx;
                fields.triangles_centroids_vy_after[i] = fields.triangles_centroids_vy_before[i] + (simulation_parameters.dt / mesh.domain_triangles[i].area) * sumy; 
            }
            std::swap(fields.triangles_centroids_vx_before, fields.triangles_centroids_vx_after);
            std::swap(fields.triangles_centroids_vy_before, fields.triangles_centroids_vy_after);
        }
        return;
    }
}