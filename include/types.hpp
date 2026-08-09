#pragma once
#include<array>
#include<vector>
#include<unordered_map>
#include<utility>
#include<functional> 

namespace fvm_types
{
    struct Node
    {
        int id;
        std::array <double, 2> position;
    };
    struct Edge
    {
        int id;
        int type;
        double length;
        std::array <int, 2> nodes_ids;
        std::array <int, 2> triangles_ids; 
        std::array <double, 2> middle_point;
        std::array <double, 2> distances_from_middle_point_to_centroids;
        std::array <double, 2> geometric_weights;
    };
    struct DomainTriangle
    {
        int id;
        double area;
        std::array <int, 3> nodes_ids;
        std::array <int, 3> edges_ids;
        std::array <int, 3> neighbor_triangles_ids;
        std::array <double, 2> centroid;
        std::array <std::array <double, 2>, 3> normal_unitary_vectors;
        std::array <std::array <double, 2>, 3> from_self_centroid_to_neighbor_triangles_centroids_vectors;
        std::array <std::array <double, 2>, 3> from_self_centroid_to_edges_middle_points;
    };
    struct GhostTriangle
    {
        int id;
        int type;
        int interior_triangle_id;
        std::array <double, 2> centroid;
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
    struct MeshData
    {
        int nodes_number;
        int edges_number;
        int boundary_edges_number;
        int domain_triangles_number;
        int ghost_triangles_number;

        std::vector <Node> nodes;
        std::vector <DomainTriangle> domain_triangles;
        std::vector <GhostTriangle> ghost_triangles;
        std::vector <Edge> edges;
        std::unordered_map <int,int> old_node_id_to_new_node_id_map;
        std::unordered_map <int, Node> nodes_list;
        std::unordered_map <int, DomainTriangle> domain_triangles_list;
        std::unordered_map <int, Edge> boundary_edges_list;
        std::unordered_map<std::pair<int, int>, int, HashPair> nodes_to_edge_map;
        std::unordered_map<int, std::vector<int>> edge_to_triangles_map;
    };
};