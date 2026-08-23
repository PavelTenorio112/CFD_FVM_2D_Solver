#include"include/time_step_computer.hpp"
#include"include/simulation_parameters.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include<cmath>
#include<limits>
namespace tsc
{
    double time_step_computer(t::MeshData &mesh, fds::Fields &fields, sp::SimulationParameters &simulation_parameters)
    {
        int i;
        double last_minimum_dt, new_minimum_dt, l1, l2, l3, area, vx, vy, v, CFL, maximum_allowed_dt;
        CFL = simulation_parameters.CFL;
        maximum_allowed_dt = 0.02;

        last_minimum_dt = std::numeric_limits<double>::max();
        for(i = 0; i < mesh.domain_triangles_number; ++i)
        {
            l1 = mesh.edges[mesh.domain_triangles[i].edges_ids[0]].length;
            l2 = mesh.edges[mesh.domain_triangles[i].edges_ids[1]].length;
            l3 = mesh.edges[mesh.domain_triangles[i].edges_ids[2]].length;
            area = mesh.domain_triangles[i].area;
            vx = fields.triangles_centroids_vx_before[i];
            vy = fields.triangles_centroids_vy_before[i];
            v = std::sqrt(vx * vx + vy * vy);
            new_minimum_dt = (CFL * 2.0 * area)/((l1 + l2 + l3) * (v + 1e-12));

            if(new_minimum_dt < last_minimum_dt)
            {
                last_minimum_dt = new_minimum_dt;
            }
        }
        if(last_minimum_dt > maximum_allowed_dt)
        {
            last_minimum_dt = maximum_allowed_dt;
        }
        return last_minimum_dt;
    }
}