#pragma once
#include <string>
#include <types.hpp>
namespace fvm_mesh_reader
{
    bool mesh_reader_function(const std::string& file_name, fvm_types::MeshData &mesh);
}