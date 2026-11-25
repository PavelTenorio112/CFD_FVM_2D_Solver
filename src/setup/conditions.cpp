#include "setup/conditions.hpp"
#include "geometric_preprocess/types.hpp"
#include "setup/fields.hpp"
#include "setup/simulation_parameters_reader.hpp"
#include <cmath>
#include <iostream>

namespace c
{
    void initial_conditions(t::MeshInfo &mesh_info, f::Fields &fields, t::Edges& edges, t::DomainTriangles& domain_triangles, sp::SimulationParameters& simulation_parameters)
    {
        fields.triangles_centroids_velocities_before.resize(mesh_info.triangles_number);
        fields.triangles_centroids_velocities_after.resize(mesh_info.triangles_number);
        fields.triangles_centroids_gradient.resize(mesh_info.triangles_number);

        fields.edges_midpoints_velocities.resize(mesh_info.edges_number);
        fields.edges_midpoints_gradient.resize(mesh_info.edges_number);



        //Quiescent initial condition
        for(int i = 0; i < mesh_info.domain_triangles_number; ++i) 
        {
            fields.triangles_centroids_velocities_before[i][0] = 0.0;
            fields.triangles_centroids_velocities_before[i][1] = 0.0;
        }


        
        //Computing of the soft start duration
        double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;
        for(int i = 0; i < mesh_info.boundary_edges_number; ++i)
        {
            if(edges.types[i] == 1)
            {
                x0 = domain_triangles.centroids[edges.triangles_IDs[i][0]][0];
                y0 = domain_triangles.centroids[edges.triangles_IDs[i][0]][1];
                break;
            }
        }
        for(int i = 0; i < mesh_info.boundary_edges_number; ++i)
        {
            if(edges.types[i] == 2)
            {
                x1 = domain_triangles.centroids[edges.triangles_IDs[i][0]][0];
                y1 = domain_triangles.centroids[edges.triangles_IDs[i][0]][1];
                break;
            }
        }
        simulation_parameters.t0 = std::abs(std::hypot(x1 - x0, y1 - y0) / std::hypot(simulation_parameters.vx_inlet, simulation_parameters.vy_inlet)) * 0.5;
        return;
    }
    
    void velocity_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles, f::Fields& fields, sp::SimulationParameters& simulation_parameters, double t)
    {
        for(int i = mesh_info.domain_triangles_number; i < mesh_info.triangles_number; ++i)
        {
            int type = ghost_triangles.types[i - mesh_info.domain_triangles_number];
            int inner_triangle_ID = ghost_triangles.inner_triangles_IDs[i - mesh_info.domain_triangles_number];
            double vx_inlet = simulation_parameters.vx_inlet;
            double vy_inlet = simulation_parameters.vy_inlet;
            double t0 = simulation_parameters.t0;
            double vx_inlet_t = 0.0;
            double vy_inlet_t = 0.0;
            double ghost_triangle_vx = 0.0;
            double ghost_triangle_vy = 0.0;
            auto[inner_triangle_vx, inner_triangle_vy] = fields.triangles_centroids_velocities_before[inner_triangle_ID];
            
            if(type == 1) // Inlet
            {
                if(t < t0)
                {
                    vx_inlet_t = 3 * vx_inlet * (t / t0) * (t / t0) - 2 * vx_inlet * (t / t0) * (t / t0) * (t / t0);
                    vy_inlet_t = 3 * vy_inlet * (t / t0) * (t / t0) - 2 * vy_inlet * (t / t0) * (t / t0) * (t / t0);
                }
                else
                {
                    vx_inlet_t = vx_inlet;
                    vy_inlet_t = vy_inlet;
                }
                ghost_triangle_vx = 2 * vx_inlet_t - inner_triangle_vx;
                ghost_triangle_vy = 2 * vy_inlet_t - inner_triangle_vy;
            }
            
            else if(type == 2) // Outlet
            {
                ghost_triangle_vx = inner_triangle_vx;
                ghost_triangle_vy = inner_triangle_vy;
            }

            else if(type == 3) // Walls
            {
                // no-slip condition
                ghost_triangle_vx = - inner_triangle_vx;
                ghost_triangle_vy = - inner_triangle_vy;

                /*// free- slip condition
                double nijx = 0.0, nijy = 0.0;
                for(int j = 0; j < 3; ++j)
                {
                    int neighbor_triangle_ID = domain_triangles.neighbor_triangles_IDs[inner_triangle_ID][j];
                    if(neighbor_triangle_ID == i)
                    {
                        nijx = domain_triangles.unit_normal_vectors[inner_triangle_ID][j][0];
                        nijy = domain_triangles.unit_normal_vectors[inner_triangle_ID][j][1];
                        break;
                    }
                }
                ghost_triangle_vx = inner_triangle_vx - 2 * (inner_triangle_vx * nijx + inner_triangle_vy * nijy) * nijx;
                ghost_triangle_vy = inner_triangle_vy - 2 * (inner_triangle_vx * nijx + inner_triangle_vy * nijy) * nijy;
                */
            }
            fields.triangles_centroids_velocities_before[i][0] = ghost_triangle_vx;
            fields.triangles_centroids_velocities_before[i][1] = ghost_triangle_vy;
        }
        return;
    }

    void gradient_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, [[maybe_unused]] t::GhostTriangles& ghost_triangles, f::Fields& fields)
    {
        for(int i = mesh_info.domain_triangles_number; i < mesh_info.triangles_number; ++i)
        {
            int inner_triangle_ID = ghost_triangles.inner_triangles_IDs[i - mesh_info.domain_triangles_number],
                edge_ID = ghost_triangles.boundary_edges_IDs[i - mesh_info.domain_triangles_number];
            double rijx = 0.0, rijy = 0.0;

            auto[Dxvxi, Dyvxi] = fields.triangles_centroids_gradient[inner_triangle_ID][0];
            auto[Dxvyi, Dyvyi] = fields.triangles_centroids_gradient[inner_triangle_ID][1];
            auto[vix, viy] = fields.triangles_centroids_velocities_before[inner_triangle_ID];
            auto[vkx, vky] = fields.triangles_centroids_velocities_before[i];
            

            for(int j = 0; j < 3; ++j)
            {
                if(domain_triangles.edges_IDs[inner_triangle_ID][j] == edge_ID)
                {
                    rijx = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[inner_triangle_ID][j][0];
                    rijy = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[inner_triangle_ID][j][1];
                    break;
                }
            }

            double rij = std::sqrt(rijx * rijx + rijy * rijy),
                   unit_rijx = rijx / rij,
                   unit_rijy = rijy / rij,
                   
                   Dxvxk = (unit_rijx / rij) * (vkx - vix) - Dxvxi,
                   Dyvxk = (unit_rijy / rij) * (vkx - vix) - Dyvxi,
                   Dxvyk = (unit_rijx / rij) * (vky - viy) - Dxvyi,
                   Dyvyk = (unit_rijy / rij) * (vky - viy) - Dyvyi;
                   


            fields.triangles_centroids_gradient[i][0][0] = Dxvxk;
            fields.triangles_centroids_gradient[i][0][1] = Dyvxk;
            fields.triangles_centroids_gradient[i][1][0] = Dxvyk;
            fields.triangles_centroids_gradient[i][1][1] = Dyvyk;
        }
        return;
    }
}
