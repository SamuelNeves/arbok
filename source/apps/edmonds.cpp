#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>

using namespace std;

struct Edge {
    double u;
    double v;
    double w;
    Edge(double start,  double end, double weight) : u(start), v(end), w(weight) {}
};

struct Node {
    unordered_map<double, Edge*> edgesComingIn;
    unordered_map<double, Edge*> edgesGoingOut;
    Edge* edgeInMDST = nullptr;
};

void MDST(double n);
double findACycle(double vertex, unordered_map<double, Edge*>& vertexToMinEdge, unordered_set<double>& toRoot, unordered_set<double>& currentPath, unordered_set<Edge*>& cycleEdges, unordered_set<double>& cycleNodes);
void DFS(double n, unordered_set<double>& seen);

vector<Node*> nodes;

int main() {
    double n, m;
    cin >> n >> m;

    nodes.resize(n);
    for (double i = 0; i < n; i++) nodes[i] = new Node;

    for (int i = 0; i < m; i++) {
        double u, v, w;
        cin >> u >> v >> w;

        if (v == 0 || u == v || (nodes[v]->edgesComingIn.find(u) != nodes[v]->edgesComingIn.end() && w >= nodes[v]->edgesComingIn[u]->w)) continue;

        Edge* edge = new Edge(u, v, w);
        nodes[u]->edgesGoingOut[v] = edge;
        nodes[v]->edgesComingIn[u] = edge;
    }

    unordered_set<double> seen;
    DFS(0, seen);
    if (seen.size() < n) {
        cout << "Sem solução!" << endl;
        return 0;
    }

    MDST(n);

    cout << "-------------- Algoritmo de Edmonds --------------" << endl;
    cout << "-------------- Árvore geradora mínima: -----------" << endl;
    double sum = 0;
    for (int i = 1; i < n; i++) {
        Edge* edge = nodes[i]->edgeInMDST;
        sum += edge->w;
        cout << edge->u << " " << edge->v << " " << edge->w << endl;
    }
    cout << "-------------- Soma dos pesos: " << sum << " -----------" << endl;

    cout << "-------------- Tempo de execução -----------------" << endl;

    return 0;
}

void MDST(double n) {
    unordered_map<double, Edge*> vertexToMinEdge;
    for (int i = 1; i < n; i++) {
        if (nodes[i]->edgesComingIn.size() == 0) continue;

        double minCost = 2147483647;
        Edge* minEdge;
        for (auto &edge : nodes[i]->edgesComingIn) {
            if (edge.second->w < minCost) {
                minEdge = edge.second;
                minCost = edge.second->w;
            }
        }

        vertexToMinEdge[i] = minEdge;
    }

    double cycleFound = 0;
    unordered_set<double> toRoot;
    unordered_set<double> currentPath;
    unordered_set<Edge*> cycleEdges;
    unordered_set<double> cycleNodes;
    for (auto &vertex : vertexToMinEdge) {
        currentPath.clear();
        cycleFound = findACycle(vertex.first, vertexToMinEdge, toRoot, currentPath, cycleEdges, cycleNodes);

        if (cycleFound) break;
    }

    if (!cycleFound) {
        for (auto &edge : vertexToMinEdge) nodes[edge.first]->edgeInMDST = edge.second;
        return;
    }

    n++;
    nodes.push_back(new Node); 

    unordered_map<Edge*, Edge*> changedEdges;
    unordered_set<Edge*> edgesToErase;
    for (auto &vertex : cycleNodes) {
        for (auto &edge : nodes[vertex]->edgesComingIn) {
            if (cycleNodes.find(edge.second->u) == cycleNodes.end()) {
                double updatedWeight = edge.second->w - vertexToMinEdge[vertex]->w;
                if (nodes[n - 1]->edgesComingIn.find(edge.second->u) != nodes[n - 1]->edgesComingIn.end()) {
                    if (updatedWeight >= nodes[n - 1]->edgesComingIn[edge.second->u]->w) {
                        edgesToErase.insert(edge.second);
                        continue;
                    } else {
                        edgesToErase.insert(changedEdges[nodes[n - 1]->edgesComingIn[edge.second->u]]);
                    }
                }

                changedEdges[edge.second] = new Edge(edge.second->u, edge.second->v, edge.second->w);
                nodes[edge.second->u]->edgesGoingOut.erase(edge.second->v);

                edge.second->v = n - 1;
                edge.second->w = updatedWeight;
                nodes[n - 1]->edgesComingIn[edge.second->u] = edge.second;
                nodes[edge.second->u]->edgesGoingOut[n - 1] = edge.second;
            } else {
                edgesToErase.insert(edge.second);
            }
        }
        nodes[vertex]->edgesComingIn.clear();
        for (auto &edge : nodes[vertex]->edgesGoingOut) {
            if (cycleNodes.find(edge.second->v) == cycleNodes.end()) {
                if (nodes[n - 1]->edgesGoingOut.find(edge.second->v) != nodes[n - 1]->edgesGoingOut.end()) {
                    if (edge.second->w >= nodes[n - 1]->edgesGoingOut[edge.second->v]->w) {
                        edgesToErase.insert(edge.second);
                        continue;
                    } else {
                        edgesToErase.insert(changedEdges[nodes[n - 1]->edgesGoingOut[edge.second->v]]);
                    }
                }

                changedEdges[edge.second] = new Edge(edge.second->u, edge.second->v, edge.second->w);
                nodes[edge.second->v]->edgesComingIn.erase(edge.second->u);

                edge.second->u = n - 1;
                nodes[n - 1]->edgesGoingOut[edge.second->v] = edge.second;
                nodes[edge.second->v]->edgesComingIn[n - 1] = edge.second;
            } else {
                edgesToErase.insert(edge.second);
            }
        }
        nodes[vertex]->edgesGoingOut.clear();
    }

    for (auto &edge : edgesToErase) {
        nodes[edge->u]->edgesGoingOut.erase(edge->v);
        nodes[edge->v]->edgesComingIn.erase(edge->u);
    }

    MDST(n);

    Edge* edgeGoingIntoConcatenatedNode = nodes[n - 1]->edgeInMDST;
    double node = changedEdges[edgeGoingIntoConcatenatedNode]->v;
    cycleEdges.erase(vertexToMinEdge[node]);

    for (auto &edge :  cycleEdges) {
        nodes[edge->v]->edgeInMDST = edge;
    }

    Edge* oldEdge = changedEdges[edgeGoingIntoConcatenatedNode];
    nodes[oldEdge->v]->edgeInMDST = edgeGoingIntoConcatenatedNode;

    for (auto &edge: changedEdges) {
        edge.first->u = edge.second->u;
        edge.first->v = edge.second->v;
        edge.first->w = edge.second->w;

        nodes[edge.first->v]->edgesComingIn[edge.first->u] = edge.first;
        nodes[edge.first->u]->edgesGoingOut[edge.first->v] = edge.first;
    }

    n--;
}

double findACycle(double vertex, unordered_map<double, Edge*>& vertexToMinEdge, unordered_set<double>& toRoot, unordered_set<double>& currentPath, unordered_set<Edge*>& cycleEdges, unordered_set<double>& cycleNodes) {
    currentPath.insert(vertex);
    
    Edge* edge = vertexToMinEdge[vertex];
    double parentVertex = edge->u;
    if (parentVertex == 0 || toRoot.find(parentVertex) != toRoot.end()) { 
        toRoot.insert(vertex);
        return false;
    } else if (currentPath.find(parentVertex) == currentPath.end()) { 
        double foundCycle = findACycle(parentVertex, vertexToMinEdge, toRoot, currentPath, cycleEdges, cycleNodes);
        if (foundCycle && foundCycle != parentVertex) { 
            cycleEdges.insert(edge);
            cycleNodes.insert(vertex);
            return foundCycle;
        } else {
            toRoot.insert(vertex);
            return false;
        }
    } else {
        cycleEdges.insert(edge);
        cycleNodes.insert(vertex);
        return parentVertex;
    }
}

void DFS(double n, unordered_set<double>& seen) {
    seen.insert(n);
    for (auto &edge : nodes[n]->edgesGoingOut) {
        if (seen.find(edge.first) == seen.end()) {
            DFS(edge.first, seen);
        }
    }
}