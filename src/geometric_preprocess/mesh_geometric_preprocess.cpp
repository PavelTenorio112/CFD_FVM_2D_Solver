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
        void memory_reserver(t::MeshInfo &mesh_info,t::Nodes &nodes, t::Edges &edges, t::DomainTriangles &domain_triangles, t::GhostTriangles &ghost_triangles)
        {

            nodes.positions.resize(mesh_info.nodes_number);

            edges.types.resize(mesh_info.boundary_edges_number);
            edges.nodes_IDs.resize(mesh_info.boundary_edges_number);
            edges.triangles_IDs.resize(mesh_info.boundary_edges_number);
            edges.lengths.resize(mesh_info.boundary_edges_number);
            edges.middle_points.resize(mesh_info.boundary_edges_number);
            edges.geometric_weights.resize(mesh_info.boundary_edges_number);
            edges.distances_from_midpoint_to_triangles_centroids.resize(mesh_info.boundary_edges_number);

            domain_triangles.nodes_IDs.resize(mesh_info.domain_triangles_number);
            domain_triangles.edges_IDs.resize(mesh_info.domain_triangles_number);
            domain_triangles.neighbor_triangles_IDs.resize(mesh_info.domain_triangles_number);
            domain_triangles.areas.resize(mesh_info.domain_triangles_number);
            domain_triangles.centroids.resize(mesh_info.domain_triangles_number);
            domain_triangles.unit_normal_vectors.resize(mesh_info.domain_triangles_number);
            domain_triangles.from_self_centroid_to_neighbor_triangles_centroids_vectors.resize(mesh_info.domain_triangles_number);
            domain_triangles.from_self_centroid_to_edges_midpoints_vectors.resize(mesh_info.domain_triangles_number);

            ghost_triangles.IDs.resize(mesh_info.boundary_edges_number);
            ghost_triangles.types.resize(mesh_info.boundary_edges_number);
            ghost_triangles.inner_triangles_IDs.resize(mesh_info.boundary_edges_number);
            ghost_triangles.boundary_edges_IDs.resize(mesh_info.boundary_edges_number);
            ghost_triangles.centroids.resize(mesh_info.boundary_edges_number);
        }

        /*
            Since our .msh file does not contain the elements and nodes ids in an ordered and consecutive way, we renumber the nodes and the elements. 
            Then, we update the nodes ids references in the triangles and boundary edges.
        */
        void arrays_generation(t::MeshMaps &mesh_maps, t::MeshReadingTools &mesh_reading_tools, t::Nodes &nodes,t::Edges &edges, t::DomainTriangles &domain_triangles)
        {
            int i;
            i = 0;
            while(!mesh_reading_tools.nodes_read_queue.empty())
            {
                double x = mesh_reading_tools.nodes_read_queue.front().x;
                double y = mesh_reading_tools.nodes_read_queue.front().y;
                int old_node_ID = mesh_reading_tools.nodes_read_queue.front().ID;

                nodes.positions[i][0] = x;
                nodes.positions[i][1] = y;
                mesh_maps.old_node_id_to_new_node_id_map[old_node_ID] = i;
                mesh_reading_tools.nodes_read_queue.pop();
                ++i;
            }

            i = 0;
            while(!mesh_reading_tools.boundary_edges_read_queue.empty())
            {
                int old_boundary_edge_node_0 = mesh_reading_tools.boundary_edges_read_queue.front().node_0_ID;
                int old_boundary_edge_node_1 = mesh_reading_tools.boundary_edges_read_queue.front().node_1_ID;
                int boundary_edge_type = mesh_reading_tools.boundary_edges_read_queue.front().type;
                int new_boundary_edge_node_0 = mesh_maps.old_node_id_to_new_node_id_map[old_boundary_edge_node_0];
                int new_boundary_edge_node_1 = mesh_maps.old_node_id_to_new_node_id_map[old_boundary_edge_node_1];

                edges.nodes_IDs[i][0] = std::min(new_boundary_edge_node_0, new_boundary_edge_node_1);
                edges.nodes_IDs[i][1] = std::max(new_boundary_edge_node_0, new_boundary_edge_node_1);
                edges.types[i] = boundary_edge_type;
                mesh_reading_tools.boundary_edges_read_queue.pop();
                ++i;
            }

            i = 0;
            while(!mesh_reading_tools.domain_triangles_read_queue.empty())
            {
                int old_domain_triangle_node_0_ID = mesh_reading_tools.domain_triangles_read_queue.front().node_0_ID;
                int old_domain_triangle_node_1_ID = mesh_reading_tools.domain_triangles_read_queue.front().node_1_ID;
                int old_domain_triangle_node_2_ID = mesh_reading_tools.domain_triangles_read_queue.front().node_2_ID;
                int new_domain_triangle_node_0_ID = mesh_maps.old_node_id_to_new_node_id_map[old_domain_triangle_node_0_ID];
                int new_domain_triangle_node_1_ID = mesh_maps.old_node_id_to_new_node_id_map[old_domain_triangle_node_1_ID];
                int new_domain_triangle_node_2_ID = mesh_maps.old_node_id_to_new_node_id_map[old_domain_triangle_node_2_ID];
                domain_triangles.nodes_IDs[i][0] = new_domain_triangle_node_0_ID;
                domain_triangles.nodes_IDs[i][1] = new_domain_triangle_node_1_ID;
                domain_triangles.nodes_IDs[i][2] = new_domain_triangle_node_2_ID;
                mesh_reading_tools.domain_triangles_read_queue.pop();
                ++i;
            }

            return;
        }
        
        /*
            We verify that the node ordering in the triangles has a counterclockwise direction. This is important to ensure a positive area and
            appropiate mathematical results.
        */
        void counterclockwise_triangles_enumeration_verifier_and_area_computing(t::MeshInfo &mesh_info, t::Nodes &nodes, t::DomainTriangles &domain_triangles)
        {
            for(int i = 0; i < mesh_info.domain_triangles_number; i++)
            {
                auto [node_0_ID, node_1_ID, node_2_ID] = domain_triangles.nodes_IDs[i];

                auto[x0, y0] = nodes.positions[node_0_ID];
                auto[x1, y1] = nodes.positions[node_1_ID];
                auto[x2, y2] = nodes.positions[node_2_ID];

                double area = 0.5 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
                if(area < 0)
                {
                    std::swap(domain_triangles.nodes_IDs[i][0], domain_triangles.nodes_IDs[i][1]);
                    domain_triangles.areas[i] = -area;
                }
                else
                {
                    domain_triangles.areas[i] = area;
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
        void nodes_to_edge_mapping(t::MeshInfo& mesh_info, t::MeshMaps &mesh_maps, t::Edges &edges, t::DomainTriangles &domain_triangles)
        {
            std::pair <int, int> nodes_pair;
            for(int i = 0; i < mesh_info.boundary_edges_number; ++i)
            {
                auto[node_0_ID, node_1_ID] = edges.nodes_IDs[i];
                nodes_pair = {node_0_ID, node_1_ID};
                mesh_maps.nodes_to_edge_map[nodes_pair] = i;
            }
            mesh_info.edges_number = mesh_info.boundary_edges_number;
            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                auto[node_0_ID, node_1_ID, node_2_ID] = domain_triangles.nodes_IDs[i];

                nodes_pair = {std::min(node_0_ID, node_1_ID), std::max(node_0_ID, node_1_ID)};
                if(mesh_maps.nodes_to_edge_map.count(nodes_pair) == 0)
                {
                    mesh_maps.nodes_to_edge_map[nodes_pair] = mesh_info.edges_number;
                    ++mesh_info.edges_number;
                }

                nodes_pair = {std::min(node_1_ID, node_2_ID), std::max(node_1_ID, node_2_ID)};
                if(mesh_maps.nodes_to_edge_map.count(nodes_pair) == 0)
                {
                    mesh_maps.nodes_to_edge_map[nodes_pair] = mesh_info.edges_number;
                    ++mesh_info.edges_number;
                }

                nodes_pair = {std::min(node_2_ID, node_0_ID), std::max(node_2_ID, node_0_ID)};
                if(mesh_maps.nodes_to_edge_map.count(nodes_pair) == 0)
                {
                    mesh_maps.nodes_to_edge_map[nodes_pair] = mesh_info.edges_number;
                    ++mesh_info.edges_number;
                }
            }
            return;
        }
        
        /*
            Using the std::unordered_map obtained previously, we add the edges IDs to each triangle.
        */
        void edges_to_triangles_assigning(t::MeshInfo& mesh_info, t::MeshMaps& mesh_maps, t::DomainTriangles& domain_triangles)
        {
            for(int i = 0; i < mesh_info.domain_triangles_number; i++)
            {
                auto [node_0_ID, node_1_ID, node_2_ID] = domain_triangles.nodes_IDs[i];

                domain_triangles.edges_IDs[i][0] = mesh_maps.nodes_to_edge_map[{std::min(node_0_ID, node_1_ID), std::max(node_0_ID, node_1_ID)}];
                domain_triangles.edges_IDs[i][1] = mesh_maps.nodes_to_edge_map[{std::min(node_1_ID, node_2_ID), std::max(node_1_ID, node_2_ID)}];
                domain_triangles.edges_IDs[i][2] = mesh_maps.nodes_to_edge_map[{std::min(node_2_ID, node_0_ID), std::max(node_2_ID, node_0_ID)}];
            }
            return;
        }
        
        /*
            Now, we create another std::unordered_map where the key is an edge ID and the value is a std::vector<int> that holds the IDs of the
            triangles that share the edge. In case an edge is used only by one domain triangle (a boundary edge), the vector will just contain one ID, 
            meaning its size will be 1.
        */
        void edge_to_triangles_mapping(t::MeshInfo& mesh_info, t::MeshMaps& mesh_maps, t::DomainTriangles& domain_triangles)
        {
            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                auto [edge_0_ID, edge_1_ID, edge_2_ID] = domain_triangles.edges_IDs[i];

                mesh_maps.edge_to_triangles_map[edge_0_ID].push_back(i);
                mesh_maps.edge_to_triangles_map[edge_1_ID].push_back(i);
                mesh_maps.edge_to_triangles_map[edge_2_ID].push_back(i);
            }
            return;   
        }
        

        void edges_vectors_resizer(t::MeshInfo& mesh_info, t::Edges& edges)
        {
            int edges_number = mesh_info.edges_number;

            edges.types.resize(edges_number);
            edges.nodes_IDs.resize(edges_number);
            edges.triangles_IDs.resize(edges_number);
            edges.lengths.resize(edges_number);
            edges.middle_points.resize(edges_number);
            edges.geometric_weights.resize(edges_number);
            edges.distances_from_midpoint_to_triangles_centroids.resize(edges_number);
        }
        
        /*
            Using the std::unordered_map obtained previously, we assign the triangles IDs to the edges. For the boundary edges the second triangle ID 
            will be "-1".
        */
        void triangles_to_edges_assigning(t::MeshInfo& mesh_info, t::MeshMaps& mesh_maps, t::Edges& edges)
        {

            for(int i = 0; i < mesh_info.boundary_edges_number; ++i)
            {
                int triangle_0_ID = mesh_maps.edge_to_triangles_map[i][0];
                edges.triangles_IDs[i][0] = triangle_0_ID;
                edges.triangles_IDs[i][1] = -1;
            }
            for(int i = mesh_info.boundary_edges_number; i < mesh_info.edges_number; ++i)
            {
                int triangle_0_ID = mesh_maps.edge_to_triangles_map[i][0];
                int triangle_1_ID = mesh_maps.edge_to_triangles_map[i][1];

                edges.triangles_IDs[i][0] = std::min(triangle_0_ID, triangle_1_ID);
                edges.triangles_IDs[i][1] = std::max(triangle_0_ID, triangle_1_ID);
                edges.types[i] = 0;
            }
            return;
        }
        
        /*
            For the domain edges, we assign their nodes IDs.
        */
        void nodes_to_domain_edges_assigning(t::MeshInfo& mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles)
        {
            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                auto [edge_0_ID, edge_1_ID, edge_2_ID] = domain_triangles.edges_IDs[i];
                auto [node_0_ID, node_1_ID, node_2_ID] = domain_triangles.nodes_IDs[i];

                edges.nodes_IDs[edge_0_ID][0] = std::min(node_0_ID, node_1_ID);
                edges.nodes_IDs[edge_0_ID][1] = std::max(node_0_ID, node_1_ID);

                edges.nodes_IDs[edge_1_ID][0] = std::min(node_1_ID, node_2_ID);
                edges.nodes_IDs[edge_1_ID][1] = std::max(node_1_ID, node_2_ID);

                edges.nodes_IDs[edge_2_ID][0] = std::min(node_2_ID, node_0_ID);
                edges.nodes_IDs[edge_2_ID][1] = std::max(node_2_ID, node_0_ID);
            }
            return;
        }
        
        /*
            We iterate over all the edges to compute their middle points and their length.
        */
        void properties_computing_1(t::MeshInfo& mesh_info, t::Nodes& nodes, t::Edges& edges)
        {
            for(int i = 0; i < mesh_info.edges_number; ++i)
            {
                auto [node_0_ID, node_1_ID] = edges.nodes_IDs[i];
                auto [x0, y0] = nodes.positions[node_0_ID];
                auto [x1, y1] = nodes.positions[node_1_ID];

                edges.lengths[i] = std::sqrt((y1 - y0) * (y1 - y0) + (x1 - x0) * (x1 - x0));

                edges.middle_points[i][0] = (x1 + x0) / 2;
                edges.middle_points[i][1] = (y1 + y0) / 2;
            }
            return;
        }
        
        /*
            Iterating over all the triangles, we compute their centroids. We also obtain the vectors that go from the centroids to the middle points and
            the normal unitary vectors.
        */
        void properties_computing_2(t::MeshInfo& mesh_info, t::Nodes& nodes, t::Edges& edges, t::DomainTriangles& domain_triangles)
        {
            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                auto [node_0_ID, node_1_ID, node_2_ID] = domain_triangles.nodes_IDs[i];
                auto [edge_0_ID, edge_1_ID, edge_2_ID] = domain_triangles.edges_IDs[i];
                auto [x0, y0] = nodes.positions[node_0_ID];
                auto [x1, y1] = nodes.positions[node_1_ID];
                auto [x2, y2] = nodes.positions[node_2_ID];

                double l0 = edges.lengths[edge_0_ID];
                double l1 = edges.lengths[edge_1_ID];
                double l2 = edges.lengths[edge_2_ID];

                domain_triangles.centroids[i][0] = (x0 + x1 + x2) / 3;
                domain_triangles.centroids[i][1] = (y0 + y1 + y2) / 3;

                domain_triangles.unit_normal_vectors[i][0][0] = (y1 - y0) / l0;
                domain_triangles.unit_normal_vectors[i][0][1] = (x0 - x1) / l0;
                domain_triangles.unit_normal_vectors[i][1][0] = (y2 - y1) / l1;
                domain_triangles.unit_normal_vectors[i][1][1] = (x1 - x2) / l1;
                domain_triangles.unit_normal_vectors[i][2][0] = (y0 - y2) / l2;
                domain_triangles.unit_normal_vectors[i][2][1] = (x2 - x0) / l2;

                domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][0][0] = edges.middle_points[edge_0_ID][0] - domain_triangles.centroids[i][0];
                domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][0][1] = edges.middle_points[edge_0_ID][1] - domain_triangles.centroids[i][1];
                domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][1][0] = edges.middle_points[edge_0_ID][0] - domain_triangles.centroids[i][0];
                domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][1][1] = edges.middle_points[edge_0_ID][1] - domain_triangles.centroids[i][1];
                domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][2][0] = edges.middle_points[edge_0_ID][0] - domain_triangles.centroids[i][0];
                domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][2][1] = edges.middle_points[edge_0_ID][1] - domain_triangles.centroids[i][1];
            }
            return;
        }
        
        /*
            We create the ghost cells, basically, we create a reflex of the boundary cells outside the domain, the only properties 
            that ghost cells have are the ID, interior cell ID, type (from the boundary edge) and centroid.
        */
        void ghost_cells_creation(t::MeshInfo& mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles)
        {   
            int triangle_ID = mesh_info.domain_triangles_number;
            double r0x = 0, r0y = 0, n0x = 0, n0y = 0, n1x = 0, n1y = 0, r1x, r1y, Det, DetX, DetY;

            for(int i = 0; i < mesh_info.boundary_edges_number; ++i)
            {
                int inner_domain_triangle_ID = edges.triangles_IDs[i][0];
                ghost_triangles.IDs[i] = triangle_ID;
                edges.triangles_IDs[i][1] = triangle_ID;
                ghost_triangles.inner_triangles_IDs[i] = edges.triangles_IDs[i][0];
                ghost_triangles.types[i] = edges.types[i];
                ghost_triangles.boundary_edges_IDs[i] = i;

                for(int j = 0; j < 3; ++j)
                {
                    if(domain_triangles.edges_IDs[inner_domain_triangle_ID][j] == i)
                    {
                        r0x = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[inner_domain_triangle_ID][j][0];
                        r0y = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[inner_domain_triangle_ID][j][1];
                        n0x = domain_triangles.unit_normal_vectors[inner_domain_triangle_ID][j][0];
                        n0y = domain_triangles.unit_normal_vectors[inner_domain_triangle_ID][j][1];
                    }
                }

                n1x = - n0x;
                n1y = - n0y;

                Det =  n1x * n1x + n1y * n1y;
                DetX = (r0x * n0x + r0y * n0y) * n1x - n1y * (r0x * n0y - r0y * n0x);
                DetY = n1x * (r0x * n0y - r0y * n0x) - (r0x * n0x + r0y * n0y) * ( - n1y);

                r1x = DetX / Det;
                r1y = DetY / Det;

                ghost_triangles.centroids[i][0] = edges.middle_points[i][0] - r1x;
                ghost_triangles.centroids[i][1] = edges.middle_points[i][0] - r1y;

                ++triangle_ID;
            }
            mesh_info.ghost_triangles_number = mesh_info.boundary_edges_number;
            mesh_info.triangles_number = mesh_info.ghost_triangles_number + mesh_info.domain_triangles_number;
            return;
        }
        
        /*
        Lets iterate over all the cells and their edges and calculate their neighbor triangles IDs.
        */
        void properties_computing_3(t::MeshInfo& mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles)
        {
            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    int edge_ID = domain_triangles.edges_IDs[i][j];
                    if(edges.triangles_IDs[edge_ID][0] == i)
                    {
                        domain_triangles.neighbor_triangles_IDs[i][j] = edges.triangles_IDs[edge_ID][1];
                    }
                    else
                    {
                        domain_triangles.neighbor_triangles_IDs[i][j] = edges.triangles_IDs[edge_ID][0];
                    }
                }
            }
            return;
        }
        
        /*
            Here we iterate over all the domain triangles and compute the vectors that go from their centroids to the neighbor's ones.
        */
        void properties_computing_4(t::MeshInfo& mesh_info, t::Edges& edges, t::DomainTriangles& domain_triangles, t::GhostTriangles& ghost_triangles)
        {
            double delta_x, delta_y;
            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                for(int j = 0; j < 3; ++j)
                {
                    int edge_ID = domain_triangles.edges_IDs[i][j];
                    int neighbor_triangle_ID = domain_triangles.neighbor_triangles_IDs[i][j];
                    if(edges.types[edge_ID] == 0)
                    {
                        delta_x = domain_triangles.centroids[neighbor_triangle_ID][0] - domain_triangles.centroids[i][0];
                        delta_y = domain_triangles.centroids[neighbor_triangle_ID][1] - domain_triangles.centroids[i][1];
                    }
                    else
                    {
                        delta_x = ghost_triangles.centroids[neighbor_triangle_ID - mesh_info.domain_triangles_number][0] - domain_triangles.centroids[i][0];
                        delta_y = ghost_triangles.centroids[neighbor_triangle_ID - mesh_info.domain_triangles_number][1] - domain_triangles.centroids[i][1];
                    }
                    domain_triangles.from_self_centroid_to_neighbor_triangles_centroids_vectors[i][j][0] = delta_x;
                    domain_triangles.from_self_centroid_to_neighbor_triangles_centroids_vectors[i][j][1] = delta_y;
                }
            }
            return;
        }
        
        /*
            Here we compute the distance from the centroids to the middle points.
        */
        void properties_computing_5(t::MeshInfo& mesh_info, t::Edges& edges , t::DomainTriangles& domain_triangles)
        {
            for(int i = 0; i < mesh_info.domain_triangles_number; ++i)
            {
                auto[edge_0_ID, edge_1_ID, edge_2_ID] = domain_triangles.edges_IDs[i];
                
                auto[delta_0_x, delta_0_y] = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][0];
                auto[delta_1_x, delta_1_y] = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][1];
                auto[delta_2_x, delta_2_y] = domain_triangles.from_self_centroid_to_edges_midpoints_vectors[i][2];



                double distance_0 = std::sqrt(delta_0_x * delta_0_x + delta_0_y * delta_0_y);
                double distance_1 = std::sqrt(delta_1_x * delta_1_x + delta_1_y * delta_1_y);
                double distance_2 = std::sqrt(delta_2_x * delta_2_x + delta_2_y * delta_2_y);

                if(edges.triangles_IDs[edge_0_ID][0] == i)
                {
                    edges.distances_from_midpoint_to_triangles_centroids[edge_0_ID][0] = distance_0;
                }
                else
                {
                    edges.distances_from_midpoint_to_triangles_centroids[edge_0_ID][1] = distance_0;
                }

                if(edges.triangles_IDs[edge_1_ID][0] == i)
                {
                    edges.distances_from_midpoint_to_triangles_centroids[edge_1_ID][0] = distance_1;
                }
                else
                {
                    edges.distances_from_midpoint_to_triangles_centroids[edge_1_ID][1] = distance_1;
                }

                if(edges.triangles_IDs[edge_2_ID][0] == i)
                {
                    edges.distances_from_midpoint_to_triangles_centroids[edge_2_ID][0] = distance_2;
                }
                else
                {
                    edges.distances_from_midpoint_to_triangles_centroids[edge_2_ID][1] = distance_2;
                }

            }
            for(int i = 0; i < mesh_info.boundary_edges_number; ++i)
            {
                edges.distances_from_midpoint_to_triangles_centroids[i][1] = edges.distances_from_midpoint_to_triangles_centroids[i][0];
            }
            return;
        }
        
        /*
            The geometric weights computing is important for linear interpolating schemes.
        */
        void properties_computing_6(t::MeshInfo& mesh_info, t::Edges& edges)
        {
            for(int i = 0; i < mesh_info.edges_number; ++i)
            {
                auto [l0, l1] = edges.distances_from_midpoint_to_triangles_centroids[i];

                double l = l0 + l1;

                edges.geometric_weights[i][0] = l1 / l;
                edges.geometric_weights[i][1] = l0 / l;
            }
            return;
        }
    }



    void mesh_geometric_preprocess(t::MeshInfo &mesh_info, t::MeshMaps &mesh_maps, t::MeshReadingTools &mesh_reading_tools, t::Nodes &nodes, t::Edges &edges, t::DomainTriangles &domain_triangles, t::GhostTriangles& ghost_triangles)
    {
        memory_reserver(mesh_info, nodes, edges, domain_triangles, ghost_triangles);
        arrays_generation(mesh_maps, mesh_reading_tools, nodes, edges, domain_triangles);
        counterclockwise_triangles_enumeration_verifier_and_area_computing(mesh_info, nodes, domain_triangles);
        nodes_to_edge_mapping(mesh_info, mesh_maps, edges, domain_triangles);
        edges_to_triangles_assigning(mesh_info, mesh_maps, domain_triangles);
        edge_to_triangles_mapping(mesh_info, mesh_maps, domain_triangles);
        edges_vectors_resizer(mesh_info, edges);
        triangles_to_edges_assigning(mesh_info, mesh_maps, edges);
        nodes_to_domain_edges_assigning(mesh_info, edges, domain_triangles);
        properties_computing_1(mesh_info, nodes, edges);
        properties_computing_2(mesh_info, nodes, edges, domain_triangles);
        ghost_cells_creation(mesh_info, edges, domain_triangles, ghost_triangles);
        properties_computing_3(mesh_info, edges, domain_triangles);
        properties_computing_4(mesh_info, edges, domain_triangles, ghost_triangles);
        properties_computing_5(mesh_info, edges, domain_triangles);
        properties_computing_6(mesh_info, edges);
        return;
    }
}