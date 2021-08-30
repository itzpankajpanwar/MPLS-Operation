#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cstring>
#define INF 1e7
using namespace std;

struct link{
    float cost, capacity;
};

struct path{
    int parent;
    float cost;
    string route;
};

struct conn_req{
    float bw_min, bw_avg, bw_max;
};

vector<vector<struct link>> topology;
vector<vector<struct conn_req>> connections;
vector<vector<struct path>> shortest_path1;
vector<vector<struct path>> shortest_path2;
string base_path = "./";
int n_nodes, n_edges, n_conn_req;

void initializeVectors(){
    topology.assign(n_nodes, vector<struct link>(n_nodes,{INF,-1}));
    shortest_path1.assign(n_nodes, vector<struct path>());
    shortest_path2.assign(n_nodes, vector<struct path>());
    connections.assign(n_nodes, vector<struct conn_req>(n_nodes, {-1,-1,-1}));
}

void readTopology(string file_name){
    string file_path = base_path + file_name;
    ifstream ifs(file_path);

    ifs >> n_nodes >> n_edges;
    initializeVectors();

    for(int i=0; i<n_edges; ++i){
        int source, dest, cost, capacity;
        ifs >> source >> dest >> cost >> capacity;
        topology[source][dest].cost = cost;
        topology[dest][source].cost = cost;
        topology[source][dest].capacity = capacity;
        topology[dest][source].capacity = capacity;
    }

    ifs.close();
}

void readConnections(string file_name){
    string file_path = base_path + file_name;
    ifstream ifs(file_path);

    ifs >> n_conn_req;
    for(int i=0; i< n_conn_req; ++i){
        int source, dest, bw_min, bw_avg, bw_max;
        ifs >> source >> dest;
        ifs >> bw_min >> bw_avg >> bw_max;
        connections[source][dest].bw_min = bw_min;
        connections[source][dest].bw_avg = bw_avg;
        connections[source][dest].bw_max = bw_max;
    }

    ifs.close();
}

void applyDijkstras(vector<vector<struct link>> graph, int s, int d, bool first){
    vector<struct path> shortest_path(n_nodes, {-1,INF,""});
    priority_queue<pair<float,int>, vector<pair<float,int>>,  greater<pair<float,int>>> pq; // creating min-heap
    pq.push(make_pair(0.0, s));
    shortest_path[s].parent = s;
    shortest_path[s].cost = 0.0;
    
    while(!pq.empty()){
        int u = pq.top().second;
        pq.pop();

        for(int v=0; v<n_nodes; ++v){
            if(shortest_path[u].cost + graph[u][v].cost < shortest_path[v].cost){
                shortest_path[v].cost = shortest_path[u].cost + graph[u][v].cost;
                shortest_path[v].parent = u;
                pq.push(make_pair(shortest_path[v].cost, v));  // decrease key
            }
        }
    }

    if(first) shortest_path1[s] = shortest_path;
    else{
       traceShortestPath(shortest_path, s, d); 
       shortest_path2[s][d].route = shortest_path[d].route;
    }
}

void removeFirstShortestPath(vector<vector<struct link>> &graph, int u, int v){
    while(v != u){
        int parent = shortest_path1[u][v].parent;
        graph[v][parent].cost = INF;
        graph[v][parent].capacity = -1;
        graph[parent][v].cost = INF;
        graph[parent][v].capacity = -1;
        v = parent;
    }
}

void traceShortestPath(vector<struct path> &shortest_path, int u, int v){
    if(shortest_path[v].parent == -1) return;

    string route = "" + to_string(v);
    int curNode = v;
    while(curNode != u){
        int parent = shortest_path[curNode].parent;
        route = to_string(parent) + " " + route;
        curNode = parent;
    }
    shortest_path[v].route = route;
}

void findShortestPaths(){
    for(int s=0; s<n_nodes; ++ s) applyDijkstras(topology, s, -1, true);

    for(int u=0; u<n_nodes; ++u)
        for(int v=0; v<n_nodes; ++v) 
            traceShortestPath(shortest_path1[u], u, v);

    for(int u=0; u<n_nodes; ++u){
        for(int v=0; v<n_nodes; ++v){
            vector<vector<struct link>> residual_graph = topology;
            removeFirstShortestPath(residual_graph, u, v);
            applyDijkstras(residual_graph, u, v, false);
            // traceShortestPath(shortest_path2, u, v);
            // cout << u << " " << v << ": " << shortest_path2[u][v].route << endl;
        }
    }
}

int main(){
    readTopology("top9.txt");
    // readConnections("NSFNET_100.txt");
    findShortestPaths();
    
    // cout << "----------------------------------" << endl;
    for(auto x:shortest_path2){
        for(auto y:x) cout << y.route << endl;
        cout << "---------------------------" << endl;
    }

    return 0;
}