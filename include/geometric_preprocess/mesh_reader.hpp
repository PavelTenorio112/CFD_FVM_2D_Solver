#pragma once
#include"include/geometric_preprocess/types.hpp"
#include <string>
namespace mr
{
    bool mesh_reader(const std::string& file_name, t::MeshInfo &mesh_info, t::MeshReadingTools &mesh_reading_tools);
}