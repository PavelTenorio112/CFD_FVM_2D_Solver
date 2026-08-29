#pragma once
#include"include/geometric_preprocess/types.hpp"
namespace mgp
{
    void mesh_geometric_preprocess(t::MeshInfo &mesh_info, t::MeshMaps &mesh_maps, t::MeshReadingTools &mesh_reading_tools, t::Nodes &nodes,
                                   t::Edges &edges, t::DomainTriangles &domain_triangles, t::GhostTriangles &ghost_triangles);
}