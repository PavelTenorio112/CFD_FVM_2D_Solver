#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include "include/setup/simulation_parameters_reader.hpp"

namespace c
{
    void initial_conditions(t::MeshInfo &mesh_info, f::Fields &fields);

    void velocity_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles, f::Fields& fields, sp::SimulationParameters& simulation_parameters);

    void gradient_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, [[maybe_unused]] t::GhostTriangles& ghost_triangles, f::Fields& fields);

}
