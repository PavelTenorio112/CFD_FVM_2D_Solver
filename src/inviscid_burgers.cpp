#include"include/inviscid_burgers.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/velocity_boundary_conditions.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
#include"include/simulation_parameters.hpp"
#include"include/upwind_interpolation.hpp"
#include "include/setup/initial_conditions.hpp"
#include"include/fields_data_print.hpp"
#include"include/least_squares_gradient_construction.hpp"
#include"include/setup/gradient_boundary_conditions.hpp"
#include"include/weights_based_linear_interpolation_for_gradient.hpp"
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
        double image_dt = (1.0 / simulation_parameters.frames_per_second),
               image_time_accumulator = 0.0,
               time_accumulator = 0.0,
               simulation_dt = simulation_parameters.time_step,
               kinematic_viscosity = simulation_parameters.kinematic_viscosity;
        
        ic::initial_conditions(mesh_info, fields);

        std::filesystem::path output_direction = "output/paraview";
        std::filesystem::create_directories(output_direction);
        fdp::fields_data_print(mesh_info, nodes, domain_triangles, fields, frame_number);
        ++frame_number;

        while(time_accumulator < simulation_parameters.total_simulation_time)
        {   
            vbc::velocity_boundary_conditions(mesh_info, domain_triangles, ghost_triangles, fields, simulation_parameters);
            lsgc::least_squares_gradient_construction(mesh_info, domain_triangles, fields);
            ui::upwind_interpolation(mesh_info, edges, domain_triangles, fields);
            gbc::gradient_boundary_conditions(mesh_info, domain_triangles, ghost_triangles, fields);
            wblifg::weights_based_linear_interpolation_for_gradients(mesh_info, edges, fields);

            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                double sum_x = 0.0,
                       sum_y = 0.0,
                       area = domain_triangles.areas[i];

                for(int j = 0; j < 3; ++j)
                {
                    int edge_ID = domain_triangles.edges_IDs[i][j];

                    double Dxvxij = fields.edges_midpoints_gradient[edge_ID][0][0],
                           Dyvxij = fields.edges_midpoints_gradient[edge_ID][0][1],
                           Dxvyij = fields.edges_midpoints_gradient[edge_ID][1][0],
                           Dyvyij = fields.edges_midpoints_gradient[edge_ID][1][1];

                    auto [nxij, nyij] = domain_triangles.unit_normal_vectors[i][j];
                    auto [vxijn, vyijn] = fields.edges_midpoints_velocities[edge_ID];
                    
                    double lij = edges.lengths[edge_ID];

                    sum_x = sum_x 
                          - (vxijn * nxij + vyijn * nyij) * vxijn * lij
                          + kinematic_viscosity * (Dxvxij * nxij + Dyvxij * nyij) * lij;
                    
                    sum_y = sum_y
                          - (vxijn * nxij + vyijn * nyij) * vyijn * lij
                          + kinematic_viscosity * (Dxvyij * nxij + Dyvyij * nyij) * lij;
                }

                fields.triangles_centroids_velocities_after[i][0] = fields.triangles_centroids_velocities_before[i][0] + (simulation_dt / area) * sum_x;
                fields.triangles_centroids_velocities_after[i][1] = fields.triangles_centroids_velocities_before[i][1] + (simulation_dt / area) * sum_y;
            }

            std::swap(fields.triangles_centroids_velocities_before, fields.triangles_centroids_velocities_after);

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