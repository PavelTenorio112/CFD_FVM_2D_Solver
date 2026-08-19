#include"include/geometric_preprocess/types.hpp"
#include<utility>
#include<cmath>
#include<fstream>
#include<iostream>
#include<iomanip>





namespace mgp
{
    namespace
    {
        /*
            Since our .msh file does not contain the elements and nodes ids in an ordered and consecutive way, we renumber the nodes and the elements. 
            Then, we update the nodes ids references in the triangles and boundary edges.
        */
        void arrays_generation(t::MeshData& mesh)
        {
            mesh.nodes.resize(mesh.nodes_number);
            mesh.edges.resize(mesh.boundary_edges_number);
            mesh.domain_triangles.resize(mesh.domain_triangles_number);
            int edge_node_0, edge_node_1, i;
            i = 0;
            for(auto it = mesh.nodes_list.begin(); it != mesh.nodes_list.end(); it++)
            {
                mesh.nodes[i] = it -> second;
                mesh.old_node_id_to_new_node_id_map[mesh.nodes[i].id] = i;
                mesh.nodes[i].id = i;
                ++i;
            }
            i = 0;
            for(auto it = mesh.boundary_edges_list.begin(); it != mesh.boundary_edges_list.end(); it++)
            {
                mesh.edges[i] = it -> second;
                mesh.edges[i].id = i;
                edge_node_0 = mesh.edges[i].nodes_ids[0];
                edge_node_1 = mesh.edges[i].nodes_ids[1];
                mesh.edges[i].nodes_ids[0] = std::min(mesh.old_node_id_to_new_node_id_map[edge_node_0], mesh.old_node_id_to_new_node_id_map[edge_node_1]);
                mesh.edges[i].nodes_ids[1] = std::max(mesh.old_node_id_to_new_node_id_map[edge_node_0], mesh.old_node_id_to_new_node_id_map[edge_node_1]);
                ++i;
            }
            i = 0;
            for(auto it = mesh.domain_triangles_list.begin(); it != mesh.domain_triangles_list.end(); it++)
            {
                mesh.domain_triangles[i] = it -> second;
                mesh.domain_triangles[i].id = i;
                mesh.domain_triangles[i].nodes_ids[0] = mesh.old_node_id_to_new_node_id_map[mesh.domain_triangles[i].nodes_ids[0]];
                mesh.domain_triangles[i].nodes_ids[1] = mesh.old_node_id_to_new_node_id_map[mesh.domain_triangles[i].nodes_ids[1]];
                mesh.domain_triangles[i].nodes_ids[2] = mesh.old_node_id_to_new_node_id_map[mesh.domain_triangles[i].nodes_ids[2]];
                ++i;
            }
            return;
        }
        
        /*
            We verify that the node ordering in the triangles has a counterclockwise direction. This is important to ensure a positive area and
            appropiate mathematical results.
        */
        void counterclockwise_triangles_enumeration_verifier_and_area_computing(t::MeshData& mesh)
        {
            double x1, x2, x3, y1, y2, y3, area;
            int i;
            for(i=0; i<mesh.domain_triangles_number; i++)
            {
                x1 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[0]].position[0];
                y1 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[0]].position[1];
                x2 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[1]].position[0];
                y2 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[1]].position[1];
                x3 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[2]].position[0];
                y3 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[2]].position[1];

                area = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
                if(area < 0)
                {
                    std::swap(mesh.domain_triangles[i].nodes_ids[0], mesh.domain_triangles[i].nodes_ids[1]);
                    mesh.domain_triangles[i].area = - area;
                }
                else
                {
                    mesh.domain_triangles[i].area = area;
                }
            }
            return;
        }
        
        /*
            To start obtaining the mesh topology, we create a std::unordered_map whose key is a std::pair<int, int> and the value is an int. To create
            this std::unordered_map we have to iterate over the boundary edges and populate our std::unordered_map. Then, we iterate 
            over all the triangles and extract node pairs, we verify the existence of these pairs within the std::unordered_map and if a pair 
            does not exist as a key, we insert it and assing a no occupied edge ID.
        */
        void nodes_to_edge_mapping(t::MeshData& mesh)
        {
            std::pair <int, int> nodes_pair;
            int i;
            for(i=0; i<mesh.boundary_edges_number; i++)
            {
                nodes_pair = {mesh.edges[i].nodes_ids[0], mesh.edges[i].nodes_ids[1]};
                mesh.nodes_to_edge_map[nodes_pair] = i;
            }
            mesh.edges_number = mesh.boundary_edges_number;
            for(i=0; i<mesh.domain_triangles_number; i++)
            {
                nodes_pair = {std::min(mesh.domain_triangles[i].nodes_ids[0], mesh.domain_triangles[i].nodes_ids[1]), std::max(mesh.domain_triangles[i].nodes_ids[0], mesh.domain_triangles[i].nodes_ids[1])};
                if(mesh.nodes_to_edge_map.count(nodes_pair) == 0)
                {
                    mesh.nodes_to_edge_map[nodes_pair] = mesh.edges_number;
                    mesh.edges_number++;
                }
                nodes_pair = {std::min(mesh.domain_triangles[i].nodes_ids[1], mesh.domain_triangles[i].nodes_ids[2]), std::max(mesh.domain_triangles[i].nodes_ids[1], mesh.domain_triangles[i].nodes_ids[2])};
                if(mesh.nodes_to_edge_map.count(nodes_pair) == 0)
                {
                    mesh.nodes_to_edge_map[nodes_pair] = mesh.edges_number;
                    mesh.edges_number++;
                }
                nodes_pair = {std::min(mesh.domain_triangles[i].nodes_ids[2], mesh.domain_triangles[i].nodes_ids[0]), std::max(mesh.domain_triangles[i].nodes_ids[2], mesh.domain_triangles[i].nodes_ids[0])};
                if(mesh.nodes_to_edge_map.count(nodes_pair) == 0)
                {
                    mesh.nodes_to_edge_map[nodes_pair] = mesh.edges_number;
                    mesh.edges_number++;
                }
            }
            return;
        }
        
        /*
            Using the std::unordered_map obtained previously, we add the edges IDs to each triangle.
        */
        void edges_to_triangles_assigning(t::MeshData& mesh)
        {
            int i, node_0_id, node_1_id, node_2_id;
            for(i=0; i<mesh.domain_triangles_number; i++)
            {
                node_0_id = mesh.domain_triangles[i].nodes_ids[0];
                node_1_id = mesh.domain_triangles[i].nodes_ids[1];
                node_2_id = mesh.domain_triangles[i].nodes_ids[2];

                mesh.domain_triangles[i].edges_ids[0] = mesh.nodes_to_edge_map[{std::min(node_0_id, node_1_id), std::max(node_0_id, node_1_id)}];
                mesh.domain_triangles[i].edges_ids[1] = mesh.nodes_to_edge_map[{std::min(node_1_id, node_2_id), std::max(node_1_id, node_2_id)}];
                mesh.domain_triangles[i].edges_ids[2] = mesh.nodes_to_edge_map[{std::min(node_2_id, node_0_id), std::max(node_2_id, node_0_id)}];
            }
            return;
        }
        
        /*
            Now, we create another std::unordered_map where the key is an edge ID and the value is a std::vector<int> that holds the IDs of the
            triangles that share the edge. In case an edge is used only by one domain triangle (a boundary edge), the vector will just contain one ID, 
            meaning its size will be 1.
        */
        void edge_to_triangles_mapping(t::MeshData& mesh)
        {
            int i, j, edge_id;
            for(i = 0; i < mesh.domain_triangles_number; ++i)
            {
                for(j=0; j<3; j++)
                {
                    edge_id = mesh.domain_triangles[i].edges_ids[j];
                    mesh.edge_to_triangles_map[edge_id].push_back(i);
                }
            }
            return;   
        }
        
        /*
            Using the std::unordered_map obtained previously, we assign the triangles IDs to the edges. For the boundary edges the second triangle ID 
            will be "-1".
        */
        void triangles_to_edges_assigning(t::MeshData& mesh)
        {
            int i;
            mesh.edges.resize(mesh.edges_number);
            for(i=0; i<mesh.boundary_edges_number; ++i)
            {
                mesh.edges[i].triangles_ids[0] = mesh.edge_to_triangles_map[i][0];
                mesh.edges[i].triangles_ids[1] = - 1;
            }
            for(i=mesh.boundary_edges_number; i<mesh.edges_number; ++i)
            {
                mesh.edges[i].id = i;
                mesh.edges[i].triangles_ids[0] = std::min(mesh.edge_to_triangles_map[i][0], mesh.edge_to_triangles_map[i][1]);
                mesh.edges[i].triangles_ids[1] = std::max(mesh.edge_to_triangles_map[i][0], mesh.edge_to_triangles_map[i][1]);
                mesh.edges[i].type = 0;
            }
            return;
        }
        
        /*
            For the domain edges, we assign their nodes IDs.
        */
        void nodes_to_domain_edges_assigning(t::MeshData& mesh)
        {
            int i;
            for(i=0; i<mesh.domain_triangles_number; ++i)
            {
                mesh.edges[mesh.domain_triangles[i].edges_ids[0]].nodes_ids[0] = std::min(mesh.domain_triangles[i].nodes_ids[0], mesh.domain_triangles[i].nodes_ids[1]);
                mesh.edges[mesh.domain_triangles[i].edges_ids[0]].nodes_ids[1] = std::max(mesh.domain_triangles[i].nodes_ids[0], mesh.domain_triangles[i].nodes_ids[1]);

                mesh.edges[mesh.domain_triangles[i].edges_ids[1]].nodes_ids[0] = std::min(mesh.domain_triangles[i].nodes_ids[1], mesh.domain_triangles[i].nodes_ids[2]);
                mesh.edges[mesh.domain_triangles[i].edges_ids[1]].nodes_ids[1] = std::max(mesh.domain_triangles[i].nodes_ids[1], mesh.domain_triangles[i].nodes_ids[2]);

                mesh.edges[mesh.domain_triangles[i].edges_ids[2]].nodes_ids[0] = std::min(mesh.domain_triangles[i].nodes_ids[2], mesh.domain_triangles[i].nodes_ids[0]);
                mesh.edges[mesh.domain_triangles[i].edges_ids[2]].nodes_ids[1] = std::max(mesh.domain_triangles[i].nodes_ids[2], mesh.domain_triangles[i].nodes_ids[0]);
            }
            return;
        }
        
        /*
            We iterate over all the edges to compute their middle points and their length.
        */
        void properties_computing_1(t::MeshData& mesh)
        {
            double x1, x2, y1, y2;
            for(int i=0; i<mesh.edges_number; i++)
            {
                x1 = mesh.nodes[mesh.edges[i].nodes_ids[0]].position[0];
                y1 = mesh.nodes[mesh.edges[i].nodes_ids[0]].position[1];
                x2 = mesh.nodes[mesh.edges[i].nodes_ids[1]].position[0];
                y2 = mesh.nodes[mesh.edges[i].nodes_ids[1]].position[1];
                
                mesh.edges[i].length = std::sqrt((y1 - y2) * (y1 - y2) + (x1 - x2) * (x1 - x2));
                mesh.edges[i].middle_point[0] = (x2 + x1) / 2;
                mesh.edges[i].middle_point[1] = (y2 + y1) / 2;
            }
            return;
        }
        
        /*
            Iterating over all the triangles, we compute their centroids. We also obtain the vectors that go from the centroids to the middle points and
            the normal unitary vectors.
        */
        void properties_computing_2(t::MeshData& mesh)
        {
            double x1, x2, x3, y1, y2, y3;
            int edge_id, i, j;
            for(i=0; i<mesh.domain_triangles_number; i++)
            {
                x1 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[0]].position[0];
                y1 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[0]].position[1];
                x2 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[1]].position[0];
                y2 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[1]].position[1];
                x3 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[2]].position[0];
                y3 = mesh.nodes[mesh.domain_triangles[i].nodes_ids[2]].position[1];

                mesh.domain_triangles[i].centroid[0] = (x1 + x2 + x3) / 3;
                mesh.domain_triangles[i].centroid[1] = (y1 + y2 + y3) / 3;

                mesh.domain_triangles[i].normal_unitary_vectors[0][0] = (y2 - y1) / mesh.edges[mesh.domain_triangles[i].edges_ids[0]].length;
                mesh.domain_triangles[i].normal_unitary_vectors[0][1] = (x1 - x2) / mesh.edges[mesh.domain_triangles[i].edges_ids[0]].length;

                mesh.domain_triangles[i].normal_unitary_vectors[1][0] = (y3 - y2) / mesh.edges[mesh.domain_triangles[i].edges_ids[1]].length;
                mesh.domain_triangles[i].normal_unitary_vectors[1][1] = (x2 - x3) / mesh.edges[mesh.domain_triangles[i].edges_ids[1]].length;

                mesh.domain_triangles[i].normal_unitary_vectors[2][0] = (y1 - y3) / mesh.edges[mesh.domain_triangles[i].edges_ids[2]].length;
                mesh.domain_triangles[i].normal_unitary_vectors[2][1] = (x3 - x1) / mesh.edges[mesh.domain_triangles[i].edges_ids[2]].length;

                for(j=0; j<3; j++)
                {
                    edge_id = mesh.domain_triangles[i].edges_ids[j];

                    mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[j][0] = mesh.edges[edge_id].middle_point[0] - mesh.domain_triangles[i].centroid[0];
                    mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[j][1] = mesh.edges[edge_id].middle_point[1] - mesh.domain_triangles[i].centroid[1];
                }
            }
            return;
        }
        
        /*
            We create the ghost cells, basically, we create a reflex of the boundary cells outside the domain, the only properties 
            that ghost cells have are the ID, interior cell ID, type (from the boundary edge) and centroid.
        */
        void ghost_cells_creation(t::MeshData& mesh)
        {   
            int triangle_id = mesh.domain_triangles_number;
            int interior_triangle_id, edge_id, i, j;
            double r1x, r1y, n1x, n1y, r2x, r2y, n2x, n2y, Det, DetX, DetY;
            mesh.ghost_triangles.resize(mesh.boundary_edges_number);
            for(i=0; i<mesh.boundary_edges_number; ++i)
            {
                mesh.ghost_triangles[i].id = triangle_id;
                mesh.edges[i].triangles_ids[1] = triangle_id;
                interior_triangle_id = mesh.edges[i].triangles_ids[0];
                mesh.ghost_triangles[i].interior_triangle_id = interior_triangle_id;
                mesh.ghost_triangles[i].type = mesh.edges[i].type;
                mesh.ghost_triangles[i].boundary_edge_id = mesh.edges[i].id;

                for(j=0; j<3; j++)
                {
                    edge_id = mesh.domain_triangles[interior_triangle_id].edges_ids[j];
                    if(mesh.edges[edge_id].id == mesh.edges[i].id)
                    {
                        r1x = mesh.domain_triangles[interior_triangle_id].from_self_centroid_to_edges_middle_points[j][0];
                        r1y = mesh.domain_triangles[interior_triangle_id].from_self_centroid_to_edges_middle_points[j][1];
                        n1x = mesh.domain_triangles[interior_triangle_id].normal_unitary_vectors[j][0];
                        n1y = mesh.domain_triangles[interior_triangle_id].normal_unitary_vectors[j][1];

                        n2x = - n1x;
                        n2y = - n1y;

                        Det =  n2x * n2x + n2y * n2y;
                        DetX = (r1x * n1x + r1y * n1y) * n2x - n2y * (r1x * n1y - r1y * n1x);
                        DetY = n2x * (r1x * n1y - r1y * n1x) - (r1x * n1x + r1y * n1y) * ( - n2y);

                        r2x = DetX / Det;
                        r2y = DetY / Det;

                        mesh.ghost_triangles[i].centroid[0] = mesh.edges[i].middle_point[0] - r2x;
                        mesh.ghost_triangles[i].centroid[1] = mesh.edges[i].middle_point[1] - r2y;

                        break;
                    }
                }
                triangle_id++;
            }
            mesh.ghost_triangles_number = mesh.boundary_edges_number;
            mesh.triangles_number = mesh.ghost_triangles_number + mesh.domain_triangles_number;
            return;
        }
        
        /*
            Here we iterate over all the domain triangles and compute the vectors that go from their centroids to the neighbor's ones.
        */
        void properties_computing_3(t::MeshData& mesh)
        {
            int edge_id, i, j;
            double delta_x, delta_y;
            for(i=0; i<mesh.domain_triangles_number; i++)
            {
                for(j=0; j<3; j++)
                {
                    edge_id = mesh.domain_triangles[i].edges_ids[j];
                    if(mesh.edges[edge_id].type == 0)
                    {
                        delta_x = mesh.domain_triangles[mesh.edges[edge_id].triangles_ids[1]].centroid[0] - mesh.domain_triangles[mesh.edges[edge_id].triangles_ids[0]].centroid[0];
                        delta_y = mesh.domain_triangles[mesh.edges[edge_id].triangles_ids[1]].centroid[1] - mesh.domain_triangles[mesh.edges[edge_id].triangles_ids[0]].centroid[1];
                        
                        if(mesh.edges[edge_id].triangles_ids[0] == mesh.domain_triangles[i].id)
                        {
                            mesh.domain_triangles[i].neighbor_triangles_ids[j] = mesh.edges[edge_id].triangles_ids[1];
                        }
                        else
                        {
                            mesh.domain_triangles[i].neighbor_triangles_ids[j] = mesh.edges[edge_id].triangles_ids[0];
                        }
                    }
                    else
                    {
                        delta_x = mesh.ghost_triangles[mesh.edges[edge_id].triangles_ids[1] - mesh.domain_triangles_number].centroid[0] - mesh.domain_triangles[mesh.edges[edge_id].triangles_ids[0]].centroid[0];
                        delta_y = mesh.ghost_triangles[mesh.edges[edge_id].triangles_ids[1] - mesh.domain_triangles_number].centroid[1] - mesh.domain_triangles[mesh.edges[edge_id].triangles_ids[0]].centroid[1];
                        
                        mesh.domain_triangles[i].neighbor_triangles_ids[j] = mesh.edges[edge_id].triangles_ids[1];
                    }
                    mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[j][0] = delta_x;
                    mesh.domain_triangles[i].from_self_centroid_to_neighbor_triangles_centroids_vectors[j][1] = delta_y;
                }
            }
            return;
        }
        
        /*
            Here we compute the distance from the centroids to the middle points.
        */
        void properties_computing_4(t::MeshData& mesh)
        {
            double delta_x, delta_y, distance;
            int edge_id, i, j;
            for(i=0; i<mesh.domain_triangles_number; i++)
            {
                for(j=0; j<3; j++)
                {
                    edge_id = mesh.domain_triangles[i].edges_ids[j];
                    delta_x = mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[j][0];
                    delta_y = mesh.domain_triangles[i].from_self_centroid_to_edges_middle_points[j][1];
                    distance = std::sqrt(delta_x * delta_x + delta_y * delta_y);

                    if(mesh.domain_triangles[i].id == mesh.edges[edge_id].triangles_ids[0])
                    {
                        mesh.edges[edge_id].distances_from_middle_point_to_centroids[0] = distance;
                    }
                    else
                    {
                        mesh.edges[edge_id].distances_from_middle_point_to_centroids[1] = distance;
                    }
                }
            }

            for(i=0; i<mesh.boundary_edges_number; ++i)
            {
                mesh.edges[i].distances_from_middle_point_to_centroids[1] = mesh.edges[i].distances_from_middle_point_to_centroids[0];
            }
        }
        
        /*
            The geometric weights computing is important for linear interpolating schemes.
        */
        void properties_computing_5(t::MeshData& mesh)
        {
            int i;
            double l0, l1;
            for(i=0; i<mesh.edges_number; i++)
            {
                l0 = mesh.edges[i].distances_from_middle_point_to_centroids[0];
                l1 = mesh.edges[i].distances_from_middle_point_to_centroids[1];

                mesh.edges[i].geometric_weights[0] = l1 / (l0 + l1);
                mesh.edges[i].geometric_weights[1] = l0 / (l0 + l1);
            }
            return;
        }
    }



    void mesh_geometric_preprocess(t::MeshData& mesh)
    {
        arrays_generation(mesh);
        counterclockwise_triangles_enumeration_verifier_and_area_computing(mesh);
        nodes_to_edge_mapping(mesh);
        edges_to_triangles_assigning(mesh);
        edge_to_triangles_mapping(mesh);
        triangles_to_edges_assigning(mesh);
        nodes_to_domain_edges_assigning(mesh);
        properties_computing_1(mesh);
        properties_computing_2(mesh);
        ghost_cells_creation(mesh);
        properties_computing_3(mesh);
        properties_computing_4(mesh);
        properties_computing_5(mesh);
        return;
    }
}