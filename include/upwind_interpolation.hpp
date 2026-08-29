#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
#include"include/setup/simulation_parameters_reader.hpp"
namespace ui
 {
    void upwind_interpolation(t::MeshInfo &mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles, f::Fields &fields);
}