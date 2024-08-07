#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// Estrutura para representar uma aresta
struct Edge {
    int u, v, weight;
};

// Função para encontrar o representante do conjunto
int findParent(int node, vector<int>& parent) {
    if (parent[node] != node) {
        parent[node] = findParent(parent[node], parent); // Compressão de caminho
    }
    return parent[node];
}

// Função para unir dois conjuntos
void unionSets(int u, int v, vector<int>& parent, vector<int>& rank) {
    int root_u = findParent(u, parent);
    int root_v = findParent(v, parent);

    if (root_u != root_v) {
        // União por rank
        if (rank[root_u] > rank[root_v]) {
            parent[root_v] = root_u;
        } else if (rank[root_u] < rank[root_v]) {
            parent[root_u] = root_v;
        } else {
            parent[root_v] = root_u;
            rank[root_u]++;
        }
    }
}

// Função que implementa o algoritmo de Kruskal
vector<Edge> kruskal(int n, vector<Edge>& edges) {
    vector<Edge> mst;
    vector<int> parent(n);
    vector<int> rank(n, 0);

    // Inicializa cada nó como seu próprio pai
    for (int i = 0; i < n; ++i) {
        parent[i] = i;
    }

    // Ordena as arestas pelo peso
    sort(edges.begin(), edges.end(), [](Edge a, Edge b) {
        return a.weight < b.weight;
    });

    // Itera sobre as arestas e adiciona à MST se não formar um ciclo
    for (const Edge& edge : edges) {
        if (findParent(edge.u, parent) != findParent(edge.v, parent)) {
            unionSets(edge.u, edge.v, parent, rank);
            mst.push_back(edge);
        }
    }

    return mst;
}

int main() {
    int n, m;
    cin >> n >> m; // Número de nós e arestas

    vector<Edge> edges(m);

    for (int i = 0; i < m; ++i) {
        cin >> edges[i].u >> edges[i].v >> edges[i].weight;
    }

    vector<Edge> mst = kruskal(n, edges);

    int totalWeight = 0;
    cout << "Arestas na MST:" << endl;
    for (const Edge& edge : mst) {
        cout << edge.u << " - " << edge.v << ": " << edge.weight << endl;
        totalWeight += edge.weight;
    }

    cout << "Peso total da MST: " << totalWeight << endl;

    return 0;
}
