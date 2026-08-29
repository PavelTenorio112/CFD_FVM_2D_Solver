#pragma once
#include"include/geometric_preprocess/types.hpp"
#include <filesystem>
namespace mdp
{
    void mesh_data_print(t::MeshInfo& mesh_info, t::Nodes& nodes, t::Edges& edges, t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles);
}