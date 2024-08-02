#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>

using namespace std;

struct Edge {
    double u;
    double v;
    double w;
    Edge(double start,  double end, double weight) {
        u = start;
        v = end;
        w = weight;
    }
};

struct Node {
    unordered_map<double, Edge*> edgesComingIn;
    unordered_map<double, Edge*> edgesGoingOut;
    Edge* edgeInMDST;
};

void MDST(double n);
double findACycle(double vertex, unordered_map<double, Edge*>& vertexToMinEdge, unordered_set<double>& toRoot, unordered_set<double>& currentPath, unordered_set<Edge*>& cycleEdges, unordered_set<double>& cycleNodes);
void DFS(double n, unordered_set<double>& seen);

vector<Node*> nodes;
// NOTE: In this program we assume that node 0 is the root node
int main() {
    double n, m;
    cin >> n;
    cin >> m;

    // Initialize empty nodes
    nodes.resize(n);

    for (double i = 0; i < n; i++) {
        nodes[i] = new Node;
    }

    // Parse edge information. Do not add edge (u,v) with weight w if:
    // 1. v = 0 since we assume 0 is our root node
    // 2. u == v since this will never be in a minimum directed spanning tree
    // 3. There is already an edge (u,v) with weight w' and w >= w'
    for (int i = 0; i < m; i++) {
        double u, v, w;
        cin >> u;
        cin >> v;
        cin >> w;

        if (v == 0 
            || u == v
            || (nodes[v]->edgesComingIn.find(u) != nodes[v]->edgesComingIn.end() 
                && w >= nodes[v]->edgesComingIn[u]->w)
        ) {
            continue;
        }

        Edge* edge = new Edge(u, v, w);
        nodes[u]->edgesGoingOut[v] = edge;
        nodes[v]->edgesComingIn[u] = edge;
    }

    // Use DFS to determine if there is a path from the root vertex to all other vertices.
    // If there is not, there is no spanning directed tree and thus no minimum directed
    // spanning tree
    unordered_set<double> seen;
    DFS(0, seen);
    if (seen.size() < n) {
        cout << "No solution!" << endl;
        return 0;
    }

    MDST(n);

    // Every vertex in the MDST (except for the root) has one and only edge going into it.
    // Sum the weights of all the edges in the MDST by traversing the singular edge going
    // into each vertex
    double sum = 0;
    for (int i = 1; i < n; i++) {
        Edge* edge = nodes[i]->edgeInMDST;
        sum += edge->w;
        cout << edge->u << " " << edge->v << " " << edge->w << endl;
    }
    cout << "Total sum: " << sum << endl;

    return 0;
}

void MDST(double n) {
    // Find the edge of minimum weight going into every vertex other than the root
    // with ties broken arbitrarily
    unordered_map<double, Edge*> vertexToMinEdge;
    for (int i = 1; i < n; i++) {
        if (nodes[i]->edgesComingIn.size() == 0) {
            continue;
        }

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

    // Determine if there is a cycle in the set of minimum weight edges for each vertex by
    // iterating over every vertex (that has a min edge in the current version of the graph) and 
    // determine whether it leads to a cycle
    double cycleFound = 0;
    // toRoot: Memoized set of vertices that lead to root node (for performance)
    unordered_set<double> toRoot;
    // currentPath: Set of vertices in current path, if we encounter one of these vertices
    // again we have found a cycle
    unordered_set<double> currentPath;
    // cycleNodes and cycleEdges: Variables to store results from function
    unordered_set<Edge*> cycleEdges;
    unordered_set<double> cycleNodes;
    for (auto &vertex : vertexToMinEdge) {
        currentPath.clear();
        cycleFound = findACycle(vertex.first, vertexToMinEdge, toRoot, currentPath, cycleEdges, cycleNodes);

        if (cycleFound) {
            break;
        }
    }

    // If there is no cycle then all the edges in vertexToMinEdge belong in the MDST
    if (!cycleFound) {
        for (auto &edge : vertexToMinEdge) {
            nodes[edge.first]->edgeInMDST = edge.second;
        }

        return;
    }

    // Create new concatenated node
    n++;
    nodes.push_back(new Node); 

    // Update the graph to have all edges with a vertex in the cycle to point to new concatenated node v_C.
    // Maintain information on what we are changing in order to reverse the changes
    unordered_map<Edge*, Edge*> changedEdges;
    unordered_set<Edge*> edgesToErase;
    for (auto &vertex : cycleNodes) {
        for (auto &edge : nodes[vertex]->edgesComingIn) {
            // Check to see if both u and v are in cycle
            if (cycleNodes.find(edge.second->u) == cycleNodes.end()) {
                // Determine updated cost of edge based on Edmonds' algorithm
                double updatedWeight = edge.second->w - vertexToMinEdge[vertex]->w;
                // Check to see if there is already an edge coming into the new vertex from
                // this vertex. If there is, update the edge if updatedWeight has a smaller weight
                // and if not, discard this edge
                if (nodes[n - 1]->edgesComingIn.find(edge.second->u) != nodes[n - 1]->edgesComingIn.end()) {
                    if (updatedWeight >= nodes[n - 1]->edgesComingIn[edge.second->u]->w) {
                        edgesToErase.insert(edge.second);
                        continue;
                    } else {
                        edgesToErase.insert(changedEdges[nodes[n - 1]->edgesComingIn[edge.second->u]]);
                    }
                }
                // Store old edge and then update edge
                changedEdges[edge.second] = new Edge(edge.second->u, edge.second->v, edge.second->w);
                nodes[edge.second->u]->edgesGoingOut.erase(edge.second->v);

                edge.second->v = n - 1;
                edge.second->w = updatedWeight;
                nodes[n - 1]->edgesComingIn[edge.second->u] = edge.second;
                nodes[edge.second->u]->edgesGoingOut[n - 1] = edge.second;
            } else {
                // If both u and v are in cycle, remove edge from graph completely
                edgesToErase.insert(edge.second);
            }
        }
        nodes[vertex]->edgesComingIn.clear();
        for (auto &edge : nodes[vertex]->edgesGoingOut) {
            // Check to see if both u and v are in cycle
            if (cycleNodes.find(edge.second->v) == cycleNodes.end()) {
                // Check to see if there is already an edge coming into the new vertex from
                // this vertex. If there is, update the edge if it has a smaller weight
                // and if not, discard this edge
                if (nodes[n - 1]->edgesGoingOut.find(edge.second->v) != nodes[n - 1]->edgesGoingOut.end()) {
                    if (edge.second->w >= nodes[n - 1]->edgesGoingOut[edge.second->v]->w) {
                        edgesToErase.insert(edge.second);
                        continue;
                    } else {
                        edgesToErase.insert(changedEdges[nodes[n - 1]->edgesGoingOut[edge.second->v]]);
                    }
                }
                // Store old edge and then update edge
                changedEdges[edge.second] = new Edge(edge.second->u, edge.second->v, edge.second->w);
                nodes[edge.second->v]->edgesComingIn.erase(edge.second->u);
                
                edge.second->u = n - 1;
                nodes[n - 1]->edgesGoingOut[edge.second->v] = edge.second;
                nodes[edge.second->v]->edgesComingIn[n - 1] = edge.second;
            } else {
                // If both u and v are in cycle, remove edge from graph completely
                edgesToErase.insert(edge.second);
            }
        }
        nodes[vertex]->edgesGoingOut.clear();
    }

    for (auto &edge : edgesToErase) {
        nodes[edge->u]->edgesGoingOut.erase(edge->v);
        nodes[edge->v]->edgesComingIn.erase(edge->u);
    }

    // Recurse on new updated graph with node v_C
    MDST(n);
    
    // Let (u, v_C) be the unique incoming edge to v_C in A'. 
    // This edge corresponds to an edge (u,v) in E with v in C. 
    // Remove the edge vertexToMinEdge(v) from C
    Edge* edgeGoingIntoConcatenatedNode = nodes[n - 1]->edgeInMDST;
    double node = changedEdges[edgeGoingIntoConcatenatedNode]->v;
    cycleEdges.erase(vertexToMinEdge[node]);

    // Add the remaining edges in the cycle to MDST
    for (auto &edge : cycleEdges) {
        nodes[edge->v]->edgeInMDST = edge;
    }

    // Insert the original edge going into concatenated node into MDST
    Edge* oldEdge = changedEdges[edgeGoingIntoConcatenatedNode];
    nodes[oldEdge->v]->edgeInMDST = edgeGoingIntoConcatenatedNode;
    // Revert all changes to edges we made
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
    if (parentVertex == 0) { 
        // Parent vertex is root vertex
        return false;
    } 

    if (toRoot.find(parentVertex) != toRoot.end()) { 
        // Parent vertex leads to root
        toRoot.insert(vertex);

        return false;
    } else if (currentPath.find(parentVertex) == currentPath.end()) { 
        // Haven't seen parent vertex before so recurse
        double foundCycle = findACycle(parentVertex, vertexToMinEdge, toRoot, currentPath, cycleEdges, cycleNodes);
        if (foundCycle && foundCycle != parentVertex) { 
            // Second check is to see if vertex does not belong in foundCycle
            cycleEdges.insert(edge);
            cycleNodes.insert(vertex);

            return foundCycle;
        } else {
            toRoot.insert(vertex);

            return false;
        }
    } else {
        // Have seen parent vertex in current path before so this is a cycle
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

