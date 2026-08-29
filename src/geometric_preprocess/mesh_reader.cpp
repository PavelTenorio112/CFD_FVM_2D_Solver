#include"include/geometric_preprocess/mesh_reader.hpp"
#include"include/geometric_preprocess/types.hpp"
#include<string>
#include<fstream>
#include<limits>
#include<iostream>
#include <filesystem>
#include<queue>

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
    bool mesh_reader(const std::string& file_name, t::MeshInfo &mesh_info, t::MeshReadingTools &mesh_reading_tools)
    {
        std::ifstream mesh_file(file_name);
        if(mesh_file.is_open())
        {
            int i, boundary_edges_and_domain_triangles_number, ID, element_type, garbage_int;
            double garbage_double;
            std::string object;
            t::NodeRead read_node;
            t::EdgeRead read_edge;
            t::DomainTriangleRead read_domain_triangle;

            mesh_info.boundary_edges_number = 0;
            mesh_info.domain_triangles_number = 0;
            while(mesh_file >> object)
            {
                if(object == "$Nodes")
                {
                    mesh_file >> mesh_info.nodes_number;
                    for(i = 0; i < mesh_info.nodes_number; ++i)
                    {
                        mesh_file >> read_node.ID >> read_node.x >> read_node.y >> garbage_double;
                        mesh_reading_tools.nodes_read_queue.push(read_node);
                    }
                }

                if(object == "$Elements")
                {
                    mesh_file >> boundary_edges_and_domain_triangles_number;
                    for(i = 0; i < boundary_edges_and_domain_triangles_number; ++i)
                    {
                        mesh_file >> ID >> element_type;
                        if(element_type == 1)
                        {
                            read_edge.ID = ID;
                            mesh_file >> garbage_int >> read_edge.type >> garbage_int >> read_edge.node_0_ID >> read_edge.node_1_ID;
                            mesh_reading_tools.boundary_edges_read_queue.push(read_edge);
                            ++mesh_info.boundary_edges_number;
                        }
                        else if(element_type == 2)
                        {
                            read_domain_triangle.ID = ID;
                            mesh_file >> garbage_int >> garbage_int >> garbage_int >> read_domain_triangle.node_0_ID >> read_domain_triangle.node_1_ID >> read_domain_triangle.node_2_ID;
                            mesh_reading_tools.domain_triangles_read_queue.push(read_domain_triangle);
                            ++mesh_info.domain_triangles_number;
                        }
                    }
                }
            }
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