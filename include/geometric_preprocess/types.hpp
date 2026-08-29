#pragma once
#include<array>
#include<vector>
#include<unordered_map>
#include<utility>
#include<functional>
#include<queue>

namespace t
{
    struct Nodes
    {
        std::vector<std::array<double, 2>> positions;
    };
    struct Edges
    {
        std::vector <int> types;
        std::vector <double> lengths;
        std::vector <std::array<int, 2>> nodes_IDs;
        std::vector <std::array<int, 2>> triangles_IDs;
        std::vector <std::array<double, 2>> middle_points;
        std::vector <std::array<double, 2>> geometric_weights;
        std::vector <std::array<double, 2>> distances_from_midpoint_to_triangles_centroids;

    };
    struct DomainTriangles
    {
        std::vector <double> areas;
        std::vector <std::array<int, 3>> nodes_IDs;
        std::vector <std::array<int, 3>> edges_IDs;
        std::vector <std::array<int, 3>> neighbor_triangles_IDs;
        std::vector <std::array<double, 2>> centroids;
        std::vector <std::array<std::array<double, 2>, 3>> unit_normal_vectors;
        std::vector <std::array<std::array<double, 2>, 3>> from_self_centroid_to_neighbor_triangles_centroids_vectors;
        std::vector <std::array<std::array<double, 2>, 3>> from_self_centroid_to_edges_midpoints_vectors;
    };
    struct GhostTriangles
    {
        std::vector <int> IDs;
        std::vector <int> types;
        std::vector <int> inner_triangles_IDs;
        std::vector <int> boundary_edges_IDs;
        std::vector <std::array<double, 2>> centroids;
    };
    struct HashPair 
    {
        template <class T1, class T2>
        std::size_t operator () (const std::pair<T1, T2>& p) const 
        {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
    struct MeshInfo
    {
        int nodes_number;
        int edges_number;
        int boundary_edges_number;
        int domain_triangles_number;
        int ghost_triangles_number;
        int triangles_number;
    };
    struct MeshMaps
    {
        std::unordered_map <int,int> old_node_id_to_new_node_id_map;
        std::unordered_map<std::pair<int, int>, int, HashPair> nodes_to_edge_map;
        std::unordered_map<int, std::vector<int>> edge_to_triangles_map;
    };
    
    struct NodeRead
    {
        int ID;
        double x, y;
        std::array <int, 2> position;
    };
    struct EdgeRead
    {
        int ID;
        int node_0_ID, node_1_ID;
        int type;
    };
    struct DomainTriangleRead
    {
        int ID;
        int node_0_ID, node_1_ID, node_2_ID;
    };
    struct MeshReadingTools
    {
        std::queue <NodeRead> nodes_read_queue;
        std::queue <EdgeRead> boundary_edges_read_queue;
        std::queue <DomainTriangleRead> domain_triangles_read_queue;
    };
}