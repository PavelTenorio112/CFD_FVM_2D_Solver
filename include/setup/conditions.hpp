#pragma once
#include"geometric_preprocess/types.hpp"
#include"setup/fields.hpp"
#include "setup/simulation_parameters_reader.hpp"

namespace c
{
    void initial_conditions(t::MeshInfo &mesh_info, f::Fields &fields, t::Edges& edges, t::DomainTriangles& domain_triangles, sp::SimulationParameters& simulation_parameters);

    void velocity_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles, f::Fields& fields, sp::SimulationParameters& simulation_parameters, double time_accumulator);

    void gradient_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, [[maybe_unused]] t::GhostTriangles& ghost_triangles, f::Fields& fields);

}
