#include"include/geometric_preprocess/mesh_data_print.hpp"
#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>





namespace mdp
{
    namespace
    {
        void nodes_data_print(t::MeshData &mesh)
        {
            int i;
            std::ofstream NodesData("output/NodesData.txt");
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
                for(i=0; i<mesh.nodes_number; ++i)
                {
                    NodesData << mesh.nodes[i].id << "\t"
                            << mesh.nodes[i].position[0] << "\t"
                            << mesh.nodes[i].position[1] << "\n";
                }
                NodesData.close();
            }
            else
            {
                std::cout<<"There was an error printing nodes data\n";
            }
            return;
        }

        void edges_data_print(t::MeshData &mesh)
        {
            int i;
            std::ofstream EdgesData("output/EdgesData.txt");
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
                for(i=0; i<mesh.edges_number; ++i)
                {
                    EdgesData << mesh.edges[i].id << "\t"
                            << mesh.edges[i].type << "\t"
                            << mesh.edges[i].nodes_ids[0] << "\t"
                            << mesh.edges[i].nodes_ids[1] << "\t"
                            << mesh.edges[i].triangles_ids[0] << "\t"
                            << mesh.edges[i].triangles_ids[1] << "\t"
                            << mesh.edges[i].length << "\t"
                            << mesh.edges[i].middle_point[0] << "\t"
                            << mesh.edges[i].middle_point[1] << "\t"
                            << mesh.edges[i].geometric_weights[0] << "\t"
                            << mesh.edges[i].geometric_weights[1] << "\t"
                            << mesh.edges[i].distances_from_middle_point_to_centroids[0] << "\t"
                            << mesh.edges[i].distances_from_middle_point_to_centroids[1] << "\n";
                }
                EdgesData.close();
            }
            else
            {
                std::cout<<"There was an error printing edges data\n";
            }
            return;
        }

        void triangles_data_print(t::MeshData &mesh)
        {
            int i;
            std::ofstream TrianglesData("output/TrianglesData.txt");
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
                for(i=0; i<mesh.domain_triangles_number; ++i)
                {
                    TrianglesData << mesh.domain_triangles[i].id << "\t"
                                << mesh.domain_triangles[i].nodes_ids[0] << "\t"
                                << mesh.domain_triangles[i].nodes_ids[1] << "\t"
                                << mesh.domain_triangles[i].nodes_ids[2] << "\t"
                                << mesh.domain_triangles[i].neighbor_triangles_ids[0] << "\t"
                                << mesh.domain_triangles[i].neighbor_triangles_ids[1] << "\t"
                                << mesh.domain_triangles[i].neighbor_triangles_ids[2] << "\t"
                                << mesh.domain_triangles[i].area << "\t"
                                << mesh.domain_triangles[i].centroid[0] << "\t"
                                << mesh.domain_triangles[i].centroid[1] << "\t"
                                << mesh.domain_triangles[i].normal_unitary_vectors[0][0] << "\t"
                                << mesh.domain_triangles[i].normal_unitary_vectors[0][1] << "\t"
                                << mesh.domain_triangles[i].normal_unitary_vectors[1][0] << "\t"
                                << mesh.domain_triangles[i].normal_unitary_vectors[1][1] << "\t"
                                << mesh.domain_triangles[i].normal_unitary_vectors[2][0] << "\t"
                                << mesh.domain_triangles[i].normal_unitary_vectors[2][1] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[0][0] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[0][1] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[1][0] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[1][1] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[2][0] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[2][1] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[0][0] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[0][1] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[1][0] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[1][1] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[2][0] << "\t"
                                << mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[2][1] << "\n";
                }
                TrianglesData.close();
            }
            else
            {
                std::cout<<"There was an error printing triangles data\n";
            }
            return;
        }
    
        void ghost_triangles_data_print(t::MeshData mesh)
        {
            int i;
            std::ofstream ghost_triangles_data("output/GhostTrianglesData.txt");
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
                for(i=0; i<mesh.boundary_edges_number; ++i)
                {
                    ghost_triangles_data << mesh.ghost_triangles[i].id << "\t"
                                         << mesh.ghost_triangles[i].type << "\t"
                                         << mesh.ghost_triangles[i].interior_triangle_id << "\t"
                                         << mesh.ghost_triangles[i].centroid[0] << "\t"
                                         << mesh.ghost_triangles[i].centroid[1] << "\n";
                }
                ghost_triangles_data.close();
            }
            else
            {
                std::cout << "There was an error printing ghost cells data. " << std::endl;
            }


        }

        void general_mesh_data_print(t::MeshData &mesh)
        {
            std::ofstream MeshData("output/MeshData.txt");
            if(MeshData.is_open())
            {
                MeshData << "MESH DATA" << "\n"
                     << "Nodes Number: " << mesh.nodes_number << "\n"
                     << "Edges Number: " << mesh.edges_number << "\n"
                     << "Boundary Edges Number: " << mesh.boundary_edges_number << "\n"
                     << "Domain Triangles Number: " << mesh.domain_triangles_number << "\n";
                MeshData.close();

            }
            else 
            {
                std::cout<<"There was an error printing mesh data\n";
            }
            return;
        }
    }
    void mesh_data_print(t::MeshData &mesh)
    {
        nodes_data_print(mesh);
        edges_data_print(mesh);
        triangles_data_print(mesh);
        ghost_triangles_data_print(mesh);
        general_mesh_data_print(mesh);
        return;
    }
}