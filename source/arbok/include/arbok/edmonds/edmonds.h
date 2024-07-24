#pragma once

#include <vector>
#include <tuple>
#include <queue>
#include <limits>
#include <arbok/data_structures/dsu.h>

namespace arbok {

class Edmonds {
protected:
    struct Edge { int from, to, weight; };

    const int num_vertices;
    std::vector<Edge> edges;

    DSU dsu; // For cycle detection and component management
    std::vector<int> in_degree;  // To track incoming edges to each node

public:
    Edmonds(int n, int m) : num_vertices(n), dsu(n), in_degree(n, 0) {
        edges.reserve(m);
    }

    void create_edge(int from, int to, int weight) {
        edges.push_back({from, to, weight});
        in_degree[to]++;
    }

    long long run(int root) {
        std::vector<Edge> contracted_edges;
        std::queue<int> q;
        std::vector<bool> in_queue(num_vertices, false);

        // Enqueue nodes with no incoming edges (except the root)
        for (int i = 0; i < num_vertices; ++i) {
            if (i != root && in_degree[i] == 0) {
                q.push(i);
                in_queue[i] = true;
            }
        }

        long long answer = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop(); in_queue[u] = false;

            auto [min_edge_weight, min_edge_index] = get_min_incoming_edge(u);
            if (min_edge_index != -1) {
                answer += min_edge_weight;
                contracted_edges.push_back(edges[min_edge_index]);
            }

            // Contract u with its parent
            if (u != root) {
                int v = edges[min_edge_index].from;
                dsu.join(u, v);
                if (!in_queue[v] && in_degree[v] == 1) {
                    q.push(v);
                    in_queue[v] = true;
                }
            }
        }
        
        // Check if all edges have been contracted, indicating the arborescence exists
        if(contracted_edges.size() != num_vertices -1){
            // There's a cycle, so no arborescence exists
            answer = -1; 
        }
        return answer;
    }
   
  std::vector<int> reconstruct(int root) {
        // Track parent of each node in the arborescence
        std::vector<int> parent(num_vertices, -1);
        
        // Iterate through contracted edges to find parents
        for (const Edge& edge : contracted_edges) {
            int u = dsu.find(edge.from); // Find the representative of the source node
            int v = dsu.find(edge.to);   // Find the representative of the target node
            if (u != v) {  // If they are not the same component (avoiding self-loops in contracted graph)
                parent[v] = u; // Set 'u' as the parent of 'v' in the arborescence
            }
        }

        // Find edges in original graph corresponding to the arborescence
        std::vector<int> arborescence_edges;
        for (int i = 0; i < edges.size(); ++i) {
            int u = dsu.find(edges[i].from);
            int v = dsu.find(edges[i].to);
            if (parent[v] == u) { 
                arborescence_edges.push_back(i); // This edge is part of the arborescence
            }
        }

        return arborescence_edges;
    }


private:
    std::vector<Edge> contracted_edges; // Store the edges that form the arborescence

    // Helper function to find the minimum weight incoming edge to a node
    std::tuple<int, int> get_min_incoming_edge(int node) {
        int min_weight = std::numeric_limits<int>::max();
        int min_index = -1;
        for (int i = 0; i < edges.size(); ++i) {
            if (edges[i].to == node && edges[i].weight < min_weight) {
                min_weight = edges[i].weight;
                min_index = i;
            }
        }
        return {min_weight, min_index};
    }
};

} // end namespace arbok