#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/simulation_parameters.hpp"

namespace ib
{
    void inviscid_burguers(t::MeshInfo &mesh_info, t::Nodes& nodes, t::Edges& edges, t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles, f::Fields &fields, sp::SimulationParameters &simulation_parameters);
}