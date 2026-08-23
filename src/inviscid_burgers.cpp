#include"include/inviscid_burgers.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/boundary_conditions.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
#include"include/simulation_parameters.hpp"
#include"include/upwind_interpolation.hpp"
#include"include/fields_data_print.hpp"
#include"include/time_step_computer.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace ib
{
    void inviscid_burguers(t::MeshData &mesh, fds::Fields &fields, sp::SimulationParameters &simulation_parameters)
    {
        int i, j, frame_number;
        double sumx, sumy, vxijn, vyijn, nxij, nyij, lij, area;
        double time_accumulator, image_dt, image_time_accumulator, simulation_dt;
        image_dt = (1.0 / simulation_parameters.frames_per_second);
        image_time_accumulator = 0.0;
        time_accumulator = 0.0;
        frame_number = 0;
        fdp::fields_data_print(mesh, fields, frame_number);
        ++frame_number;
        while(time_accumulator < simulation_parameters.total_simulation_time)
        {   
            bc::boundary_conditions(mesh, fields, simulation_parameters);
            ui::upwind_interpolation(mesh, fields, simulation_parameters);
            simulation_dt = tsc::time_step_computer(mesh, fields, simulation_parameters);
            for(i = 0; i < mesh.domain_triangles_number; ++i)
            {
                sumx = 0;
                sumy = 0;
                area = mesh.domain_triangles[i].area;
                for(j=0; j<3; ++j)
                {
                    vxijn = fields.edges_vx[mesh.domain_triangles[i].edges_ids[j]];
                    vyijn = fields.edges_vy[mesh.domain_triangles[i].edges_ids[j]];
                    nxij = mesh.domain_triangles[i].normal_unitary_vectors[j][0];
                    nyij = mesh.domain_triangles[i].normal_unitary_vectors[j][1];
                    lij = mesh.edges[mesh.domain_triangles[i].edges_ids[j]].length;

                    sumx = sumx - (vxijn * nxij + vyijn * nyij) * vxijn * lij;
                    sumy = sumy - (vxijn * nxij + vyijn * nyij) * vyijn * lij;
                }
                fields.triangles_centroids_vx_after[i] = fields.triangles_centroids_vx_before[i] + (simulation_dt / area) * sumx;
                fields.triangles_centroids_vy_after[i] = fields.triangles_centroids_vy_before[i] + (simulation_dt / area) * sumy; 
            }

            std::swap(fields.triangles_centroids_vx_before, fields.triangles_centroids_vx_after);
            std::swap(fields.triangles_centroids_vy_before, fields.triangles_centroids_vy_after);

            if(image_time_accumulator > image_dt)
            {
                fdp::fields_data_print(mesh, fields, frame_number);
                image_time_accumulator = image_time_accumulator - image_dt;
                ++frame_number;
            }

            time_accumulator = time_accumulator + simulation_dt;
            image_time_accumulator = image_time_accumulator + simulation_dt;
        }
        return;
    }
}