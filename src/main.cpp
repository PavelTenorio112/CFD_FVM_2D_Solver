#include "types.hpp"
#include "mesh_reader.hpp"
#include "mesh_geometric_processing.hpp"
#include "mesh_data_printing.hpp"
#include <iostream>
#include <string>
int main()
{
    fvm_types::MeshData mesh;
    std::string file_name = "input/Mesh.msh";
    bool success = fvm_mesh_reader::mesh_reader_function(file_name, mesh);
    if(success != true)
    {
        std::cout<<"There was an error reading the mesh data."<<std::endl;
        return 1;
    }
    
    fvm_mesh_geometric_processing::mesh_geometric_processing_function(mesh);

    fvm_mesh_data_printing::mesh_data_printing_function(mesh);

    return 0;
}