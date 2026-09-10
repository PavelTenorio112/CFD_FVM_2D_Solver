#pragma once
#include"include/geometric_preprocess/types.hpp"
#include"include/setup/fields.hpp"
namespace wblifg
{
    void weights_based_linear_interpolation_for_gradients(t::MeshInfo& mesh_info, t::Edges& edges, f::Fields &fields);
}