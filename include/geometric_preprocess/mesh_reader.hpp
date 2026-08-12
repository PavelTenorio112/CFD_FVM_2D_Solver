#pragma once
#include"include/geometric_preprocess/types.hpp"
#include <string>
namespace mr
{
    bool mesh_reader(const std::string& file_name, t::MeshData &mesh);
}