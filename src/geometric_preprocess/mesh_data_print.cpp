#include"include/geometric_preprocess/mesh_data_print.hpp"
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<filesystem>





namespace mdp
{
    namespace
    {
        void nodes_data_print(t::MeshInfo& mesh_info, t::Nodes& nodes)
        {
            std::filesystem::path output_direction = "output/text";
            std::filesystem::create_directories(output_direction);

            int i;
            std::ofstream NodesData("output/text/NodesData.txt");
            if(NodesData.is_open())
            {
                NodesData << "NODES DATA\n"
                          << "-------------------------------------------------------------------------------\n"
                          << "Column headers (in order)\n"
                          << "ID\n"
                          << "X\n"
                          << "Y\n"
                          << "-------------------------------------------------------------------------------\n";
                NodesData << std::fixed << std::setprecision(8);
                for(i = 0; i < mesh_info.nodes_number; ++i)
                {
                    NodesData << i << "\t"
                              << nodes.positions[i][0] << "\t"
                              << nodes.positions[i][1] << "\n";
                }
                NodesData.close();
            }
            else
            {
                std::cout<<"There was an error printing nodes data\n";
            }
            return;
        }

        void edges_data_print(t::MeshInfo& mesh_info, t::Edges& edges)
        {
            int i;
            std::ofstream EdgesData("output/text/EdgesData.txt");
            if(EdgesData.is_open())
            {
                EdgesData << "EDGES DATA\n"
                        << "-------------------------------------------------------------------------------\n"
                        << "Column headers (in order)\n"
                        << "ID\n"
                        << "Type\n"
                        << "Node 0 ID\n"
                        << "Node 1 ID\n"
                        << "Cell 0 ID\n"
                        << "Cell 1 ID\n"
                        << "Length\n"
                        << "Middle Point (x)\n"
                        << "Middle Point(y)\n"
                        << "Weight 0\n"
                        << "Weight 1\n"
                        << "Distance from cell 0 centroid to middle point\n"
                        << "Distance from cell 1 centroid to middle point\n"
                        << "-------------------------------------------------------------------------------\n";
                EdgesData << std::fixed << std::setprecision(8);
                for(i = 0; i < mesh_info.edges_number; ++i)
                {
                    EdgesData << i << "\t"
                              << edges.types[i] << "\t"
                              << edges.nodes_IDs[i][0] << "\t"
                              << edges.nodes_IDs[i][1] << "\t"
                              << edges.triangles_IDs[i][0] << "\t"
                              << edges.triangles_IDs[i][1] << "\t"
                              << edges.lengths[i] << "\t"
                              << edges.middle_points[i][0] << "\t"
                              << edges.middle_points[i][1] << "\t"
                              << edges.geometric_weights[i][0] << "\t"
                              << edges.geometric_weights[i][1] << "\t"
                              << edges.distances_from_midpoint_to_triangles_centroids[i][0] << "\t"
                              << edges.distances_from_midpoint_to_triangles_centroids[i][1] << "\n";
                }
                EdgesData.close();
            }
            else
            {
                std::cout<<"There was an error printing edges data\n";
            }
            return;
        }

        void domain_triangles_data_print(t::MeshInfo& mesh_info, t::DomainTriangles& domain_triangles)
        {
            std::ofstream TrianglesData("output/text/TrianglesData.txt");
            if(TrianglesData.is_open())
            {
                TrianglesData << "TRIANGLES DATA\n"
                            << "-------------------------------------------------------------------------------\n"
                            << "Column headers (in order)\n"
                            << "ID\n"
                            << "Node 0 ID\n"
                            << "Node 1 ID\n"
                            << "Node 2 ID\n"
                            << "Neighbor Triangle 0 ID\n"
                            << "Neighbor Triangle 1 ID\n"
                            << "Neighbor Triangle 2 ID\n"
                            << "Area\n"
                            << "Centroid (x)\n"
                            << "Centroid (y)\n"
                            << "Normal Unitary Vector 0 (x)\n"
                            << "Normal Unitary Vector 0 (y)\n"
                            << "Normal Unitary Vector 1 (x)\n"
                            << "Normal Unitary Vector 1 (y)\n"
                            << "Normal Unitary Vector 2 (x)\n"
                            << "Normal Unitary Vector 2 (y)\n"
                            << "from self centroid to neighbor triangle 0 centroid vector (x)\n"
                            << "from self centroid to neighbor triangle 0 centroid vector (y)\n"
                            << "from self centroid to neighbor triangle 1 centroid vector (x)\n"
                            << "from self centroid to neighbor triangle 1 centroid vector (y)\n"
                            << "from self centroid to neighbor triangle 2 centroid vector (x)\n"
                            << "from self centroid to neighbor triangle 2 centroid vector (y)\n"
                            << "from self centroid to edge 0 middle point (x)\n"
                            << "from self centroid to edge 0 middle point (y)\n"
                            << "from self centroid to edge 1 middle point (x)\n"
                            << "from self centroid to edge 1 middle point (y)\n"
                            << "from self centroid to edge 2 middle point (x)\n"
                            << "from self centroid to edge 2 middle point (y)\n"
                            << "-------------------------------------------------------------------------------\n";
                TrianglesData << std::fixed << std::setprecision(8);
                for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
                {
                    TrianglesData << i << "\t";
                    for(int j = 0; j < 3; ++j)
                    {
                        TrianglesData << domain_triangles.nodes_IDs[i][j] << "\t";
                    }
                    for(int j = 0; j < 3; ++j)
                    {
                        TrianglesData << domain_triangles.neighbor_triangles_IDs[i][j] << "\t";
                    }
                    TrianglesData << domain_triangles.areas[i] << "\t"
                                  << domain_triangles.centroids[i][0] << "\t"
                                  << domain_triangles.centroids[i][1] << "\t";

                    for(int j = 0; j < 3; ++j)
                    {
                        TrianglesData << domain_triangles.unit_normal_vectors[i][j][0] << "\t"
                                      << domain_triangles.unit_normal_vectors[i][j][1] << "\t";
                    }
                    for(int j = 0; j < 3; ++j)
                    {
                        TrianglesData << domain_triangles.from_self_centroid_to_neighbor_triangles_centroids_vectors[i][j][0] << "\t"
                                      << domain_triangles.from_self_centroid_to_neighbor_triangles_centroids_vectors[i][j][1] << "\t";
                    }
                    for(int j = 0; j < 3; ++j)
                    {
                        TrianglesData << domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][j][0] << "\t"
                                      << domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][j][1] << "\t";
                    }
                }
                TrianglesData.close();
            }
            else
            {
                std::cout<<"There was an error printing triangles data\n";
            }
            return;
        }
    
        void ghost_triangles_data_print(t::MeshInfo& mesh_info, t::GhostTriangles& ghost_triangles)
        {
            std::ofstream ghost_triangles_data("output/text/GhostTrianglesData.txt");
            if(ghost_triangles_data.is_open())
            {
                ghost_triangles_data << "GHOST TRIANGLES DATA" << "\n"
                                     << "-------------------------------------------------------------------------------\n"
                                     << "Column headers (in order): " << "\n"
                                     << "Triangle ID" << "\n"
                                     << "Boundary type" << "\n"
                                     << "Inner triangle ID" << "\n"
                                     << "Centroid (x)" << "\n"
                                     << "Centroid (y)" << "\n";
                for(int i = 0; i < mesh_info.boundary_edges_number; ++i)
                {
                    ghost_triangles_data << ghost_triangles.IDs[i] << "\t"
                                         << ghost_triangles.types[i] << "\t"
                                         << ghost_triangles.inner_triangles_IDs[i] << "\t"
                                         << ghost_triangles.centroids[i][0] << "\t"
                                         << ghost_triangles.centroids[i][1] << "\n";
                }
                ghost_triangles_data.close();
            }
            else
            {
                std::cout << "There was an error printing ghost cells data. " << std::endl;
            }


        }

        void general_mesh_data_print(t::MeshInfo& mesh_info)
        {
            std::ofstream MeshData("output/text/MeshData.txt");
            if(MeshData.is_open())
            {
                MeshData << "MESH DATA" << "\n"
                         << "Nodes Number: " << mesh_info.nodes_number << "\n"
                         << "Edges Number: " << mesh_info.edges_number << "\n"
                         << "Boundary Edges Number: " << mesh_info.boundary_edges_number << "\n"
                         << "Domain Triangles Number: " << mesh_info.domain_triangles_number << "\n";
                MeshData.close();

            }
            else 
            {
                std::cout<<"There was an error printing mesh data\n";
            }
            return;
        }
    }
    void mesh_data_print(t::MeshInfo& mesh_info, t::Nodes& nodes, t::Edges& edges, t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles)
    {
        nodes_data_print(mesh_info, nodes);
        edges_data_print(mesh_info, edges);
        domain_triangles_data_print(mesh_info, domain_triangles);
        ghost_triangles_data_print(mesh_info, ghost_triangles);
        general_mesh_data_print(mesh_info);
        return;
    }
}