#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
namespace ui
 {
    void upwind_interpolation(const t::MeshData& mesh, fds::Fields& fields, const sp::SimulationParameters& simulation_parameters);
}