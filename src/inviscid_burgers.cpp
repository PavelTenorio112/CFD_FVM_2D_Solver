#include"include/inviscid_burgers.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/boundary_conditions.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
#include"include/simulation_parameters.hpp"
#include"include/upwind_interpolation.hpp"
#include "include/setup/initial_conditions.hpp"
#include"include/fields_data_print.hpp"
#include"include/time_step_computer.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <filesystem>

namespace ib
{
    void inviscid_burguers(t::MeshInfo &mesh_info, t::Nodes& nodes, t::Edges& edges, t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles, f::Fields &fields, sp::SimulationParameters &simulation_parameters)
    {
        int frame_number = 0;
        double image_dt = (1.0 / simulation_parameters.frames_per_second);
        double image_time_accumulator = 0.0;
        double time_accumulator = 0.0;
        double simulation_dt = simulation_parameters.time_step;
        
        ic::initial_conditions(mesh_info, fields);

        std::filesystem::path output_direction = "output/paraview";
        std::filesystem::create_directories(output_direction);
        fdp::fields_data_print(mesh_info, nodes, domain_triangles, fields, frame_number);
        ++frame_number;

        while(time_accumulator < simulation_parameters.total_simulation_time)
        {   
            bc::boundary_conditions(mesh_info, domain_triangles, ghost_triangles, fields, simulation_parameters);

            ui::upwind_interpolation(mesh_info, edges, domain_triangles, fields);

            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                double sum_x = 0;
                double sum_y = 0;
                double area = domain_triangles.areas[i];

                for(int j = 0; j < 3; ++j)
                {
                    int edge_ID = domain_triangles.edges_IDs[i][j];

                    auto [nxij, nyij] = domain_triangles.unit_normal_vectors[i][j];
                    auto [vxijn, vyijn] = fields.edges_midpoints_velocities[edge_ID];
                    
                    double lij = edges.lengths[edge_ID];

                    sum_x = sum_x - (vxijn * nxij + vyijn * nyij) * vxijn * lij;
                    sum_y = sum_y - (vxijn * nxij + vyijn * nyij) * vyijn * lij; 
                }

                fields.domain_triangles_centroids_velocities_after[i][0] = fields.domain_triangles_centroids_velocities_before[i][0] + (simulation_dt / area) * sum_x;
                fields.domain_triangles_centroids_velocities_after[i][1] = fields.domain_triangles_centroids_velocities_before[i][1] + (simulation_dt / area) * sum_y;
            }
            std::swap(fields.domain_triangles_centroids_velocities_before, fields.domain_triangles_centroids_velocities_after);

            if(image_time_accumulator > image_dt)
            {
                fdp::fields_data_print(mesh_info, nodes, domain_triangles, fields, frame_number);
                image_time_accumulator = image_time_accumulator - image_dt;
                ++frame_number;
            }

            time_accumulator = time_accumulator + simulation_dt;
            image_time_accumulator = image_time_accumulator + simulation_dt;
        }
        return;
    }
}