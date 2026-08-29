#include "include/geometric_preprocess/types.hpp"
#include "include/geometric_preprocess/mesh_reader.hpp"
#include "include/geometric_preprocess/mesh_geometric_preprocess.hpp"
#include "include/geometric_preprocess/mesh_data_print.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/simulation_parameters_reader.hpp"
#include "include/simulation_parameters.hpp"
#include "include/inviscid_burgers.hpp"
#include <iostream>
#include <string>
#include <fstream>

int main()
{
    t::MeshInfo mesh_info;
    t::MeshMaps mesh_maps;
    t::MeshReadingTools mesh_reading_tools;
    t::Nodes nodes;
    t::Edges edges;
    t::DomainTriangles domain_triangles;
    t::GhostTriangles ghost_triangles;
    std::string file_name = "input/Mesh.msh";

    bool success = mr::mesh_reader(file_name, mesh_info, mesh_reading_tools);
    if(success != true)
    {
        return 1;
    }
    
    mgp::mesh_geometric_preprocess(mesh_info, mesh_maps, mesh_reading_tools, nodes, edges, domain_triangles, ghost_triangles);
    
    mdp::mesh_data_print(mesh_info, nodes, edges, domain_triangles, ghost_triangles);

    f::Fields fields;
    sp::SimulationParameters simulation_parameters;
    spr::simulation_parameters_reader(simulation_parameters);
    ib::inviscid_burguers(mesh_info, nodes, edges, domain_triangles, ghost_triangles, fields, simulation_parameters);

    return 0;
}