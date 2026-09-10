#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
#include"include/simulation_parameters.hpp"
namespace gbc
{
    void gradient_boundary_conditions(t::MeshInfo& mesh_info, [[maybe_unused]] t::DomainTriangles& domain_triangles, [[maybe_unused]] t::GhostTriangles& ghost_triangles, f::Fields& fields);
}