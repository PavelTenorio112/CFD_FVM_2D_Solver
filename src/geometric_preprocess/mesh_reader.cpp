#include"include/geometric_preprocess/mesh_reader.hpp"
#include"types.hpp"
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
        std::ifstream file(file_name);
        if(file.is_open())
        {
            std::string line;
            t::Node node;
            t::DomainTriangle triangle;
            t::Edge edge;
            int edges_and_domain_triangles_number;
            double i, a1, a2, a3, a4, a5, a6, a7, a8;
            while(std::getline(file, line))
            {
                if(line == "$Nodes")
                {
                    file >> mesh.nodes_number;
                    for(i=1; i<= mesh.nodes_number; i++)
                    {
                        file >> a1 >> node.position[0] >> node.position[1] >> a4;
                        node.id=a1;
                        mesh.nodes_list[node.id] = node;
                    }
                    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');/* Ignora el salto de linea */
                }
                else if(line=="$Elements")
                {
                    file >> edges_and_domain_triangles_number;
                    for(i=1; i<=edges_and_domain_triangles_number; i++)
                    {
                        file>>a1>>a2;
                        if(a2 == 1)
                        {
                            file >> a3 >> a4 >> a5 >> a6 >> a7;
                            edge.id = a1;
                            edge.type = a4;
                            edge.nodes_ids[0] = std::min(a6,a7);
                            edge.nodes_ids[1] = std::max(a6,a7);
                            mesh.boundary_edges_list[edge.id] = edge;
                        }
                        else if(a2==2)
                        {
                            file >> a3 >> a4 >> a5 >>a6>>a7>>a8;
                            triangle.id=a1;
                            triangle.nodes_ids[0]=a6;
                            triangle.nodes_ids[1]=a7;
                            triangle.nodes_ids[2]=a8;
                            mesh.domain_triangles_list[triangle.id] = triangle;
                        }
                    }
                    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');/* Ignora el salto de linea*/
                }
            }
            mesh.boundary_edges_number = mesh.boundary_edges_list.size();
            mesh.domain_triangles_number = mesh.domain_triangles_list.size();
            file.close();
            std::cout << "The mesh was read succesfully. " << std::endl;
            return true;
        }
        else
        {
            return false;
        }
    }
}