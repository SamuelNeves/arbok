#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include <arbok/data_structures/dsu.h>

namespace arbok {

class Kruskal {
protected:
    struct Edge { int from, to, weight; };

    const int num_vertices;
    std::vector<Edge> edges;

    DSU dsu; // Disjoint set union to track connected components

public:
    Kruskal(int n, int m) : num_vertices(n), dsu(n) {
        edges.reserve(m);
    }

    void create_edge(int from, int to, int weight) {
        edges.push_back({from, to, weight});
    }

    // long long run(int root) {
    //     // Sort edges in non-decreasing order of their weights
    //     std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
    //         return a.weight < b.weight;
    //     });

    //     long long answer = 0;
    //     for (const Edge& edge : edges) {
    //         // Check if including this edge creates a cycle
    //         if (dsu.join(edge.from, edge.to)) {
    //             answer += edge.weight;
    //         }
    //     }
    //     return answer;
    // }
  long long run(int /*root*/) {  // Accept the 'root' argument but ignore it
        // Sort edges in non-decreasing order of their weights (unchanged)
        std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;
        });

        // Rest of the MST construction logic remains the same
        long long answer = 0;
        for (const Edge& edge : edges) {
            if (dsu.join(edge.from, edge.to)) {
                answer += edge.weight;
            }
        }
        return answer;
    }
    std::vector<int> reconstruct(int /*root*/) {
        // For Kruskal, the result is the set of edges included in the MST (edges added during `run()`)
        std::vector<int> res;
        for (int i = 0; i < edges.size(); ++i) {
            if (dsu.same(edges[i].from, edges[i].to)) { // Edge is in MST
                res.push_back(i); // Store the index of the edge in the original `edges` vector
            }
        }
        return res;
    }
};

} // end namespace arbok