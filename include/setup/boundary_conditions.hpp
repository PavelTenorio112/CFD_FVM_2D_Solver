#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
#include"include/simulation_parameters.hpp"
namespace bc
{
    void boundary_conditions(t::MeshData &mesh, fds::Fields &fields, sp::SimulationParameters &simulation_parameters);
}