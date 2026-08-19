#include"include/geometric_preprocess/mesh_reader.hpp"
#include"include/geometric_preprocess/types.hpp"
#include<string>
#include<fstream>
#include<limits>
#include<iostream>
#include <filesystem>

/*
        Gmsh generates a .msh file whose entities section is divided into two sub-sections, one for nodes and another for elements. The nodes section 
    contains 4 numbers per line, the number in the first column is for the id, the following three are for the node spatial coordinates. The elements 
    section contains eight numbers per line, the first one is for the id and the second one is for the type of element (1 for edges, 2 for triangles). 
    In case we have an edge, the fourth number is for the type of boundary (1 for inlet, 2 for outlet, 3 for wall). And depending on the type of element, 
    the remaining numbers are for the ids of the used nodes. This function reads all of this information and puts them into three std::unordered_map, 
    one for nodes, one for edges and another for triangles, where the key is the id and the object is the value. It is important to remember that Gmsh
    only gives us the boundary edges, not the interior ones.
*/
namespace mr
{
    bool mesh_reader(const std::string& file_name, t::MeshData &mesh)
    {
        std::ifstream mesh_file(file_name);
        if(mesh_file.is_open())
        {
            t::Node node;
            t::Edge edge;
            t::DomainTriangle domain_triangle;
            int i, boundary_edges_and_domain_triangles_number, id, element_type, garbage_int;
            double garbage_double;
            std::string object;

            while(mesh_file >> object)
            {
                if(object == "$Nodes")
                {
                    mesh_file >> mesh.nodes_number;
                    for(i = 0; i < mesh.nodes_number; ++i)
                    {
                        mesh_file >> node.id >> node.position[0] >> node.position[1] >> garbage_double;
                        mesh.nodes_list[node.id] = node;
                    }
                }
                if(object == "$Elements")
                {
                    mesh_file >> boundary_edges_and_domain_triangles_number;
                    for(i = 0; i < boundary_edges_and_domain_triangles_number; ++i)
                    {
                        mesh_file >> id >> element_type;
                        if(element_type == 1)
                        {
                            edge.id = id;
                            mesh_file >> garbage_int >> edge.type >> garbage_int >> edge.nodes_ids[0] >> edge.nodes_ids[1];
                            mesh.boundary_edges_list[edge.id] = edge;

                        }
                        else if(element_type == 2)
                        {
                            domain_triangle.id = id;
                            mesh_file >> garbage_int >> garbage_int >> garbage_int >> domain_triangle.nodes_ids[0] >> domain_triangle.nodes_ids[1] >> domain_triangle.nodes_ids[2];
                            mesh.domain_triangles_list[domain_triangle.id] = domain_triangle;
                        }
                    }
                }
            }

            mesh.boundary_edges_number = mesh.boundary_edges_list.size();
            mesh.domain_triangles_number = mesh.domain_triangles_list.size();
            mesh_file.close();
            std::cout << "The mesh was read succesfully. " << std::endl;
            return true;
        }
        else
        {
            std::cout << "There was an error reading mesh data. " << std::endl;
            return false;
        }
    }
}