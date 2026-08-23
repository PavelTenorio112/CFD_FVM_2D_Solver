#pragma once
#include"include/simulation_parameters.hpp"
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
namespace tsc
{
    double time_step_computer(t::MeshData &mesh, fds::Fields &fields, sp::SimulationParameters &simulation_parameters);
}