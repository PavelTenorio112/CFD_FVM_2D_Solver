#include "include/geometric_preprocess/types.hpp"
#include "include/geometric_preprocess/mesh_reader.hpp"
#include "include/geometric_preprocess/mesh_geometric_preprocess.hpp"
#include "include/geometric_preprocess/mesh_data_print.hpp"
#include "include/setup/fields.hpp"
#include "include/setup/initial_conditions.hpp"
#include "include/setup/simulation_parameters_reader.hpp"
#include "include/simulation_parameters.hpp"
#include "include/inviscid_burgers.hpp"
#include <iostream>
#include <string>
#include <fstream>

int main()
{
    t::MeshData mesh;
    std::string file_name = "input/Mesh.msh";
    bool success = mr::mesh_reader(file_name, mesh);
    if(success != true)
    {
        return 1;
    }
    mgp::mesh_geometric_preprocess(mesh);
    mdp::mesh_data_print(mesh);

    fds::Fields fields;
    sp::SimulationParameters simulation_parameters;
    spr::simulation_parameters_reader(simulation_parameters);
    ic::initial_conditions(mesh, fields);
    ib::inviscid_burguers(mesh, fields, simulation_parameters);

    return 0;
}