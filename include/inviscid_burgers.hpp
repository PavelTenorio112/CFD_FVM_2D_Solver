#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/simulation_parameters.hpp"

namespace ib
{
    void inviscid_burguers(t::MeshData& mesh, fds::Fields& field, sp::SimulationParameters& simulation_parameters);
}