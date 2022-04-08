#include <iostream>
#include <vector>
#include <stack>
#include <unordered_map>
using namespace std;

unordered_map<int, vector<int>> tarjan(vector<vector<int>>& graph); 
void dfs(int current, vector<vector<int>>& graph,
         int& id, vector<int>& ids, vector<int>& low, 
         stack<int>& Stack, vector<bool>& onStack);

unordered_map<int, vector<int>> tarjan(vector<vector<int>>& graph) {
    int size = graph.size(); 
    
    int id = 0; 
    vector<int> ids(size, -1); 
    vector<int> low(size, 0); 
    vector<bool> onStack(size, false); 
    stack<int> Stack; 
    
    for (int node = 0; node < size; node++) {
        if (ids[node] == -1) {
            dfs(node, graph, id, ids, low, Stack, onStack); 
        }
    }
    
    unordered_map<int, vector<int>> scc; 
    for (int index = 0; index < size; index++) {
        scc[low[index]].push_back(index); 
    }
    return scc; 
}

void dfs(int current, vector<vector<int>>& graph,
         int& id, vector<int>& ids, vector<int>& low, 
         stack<int>& Stack, vector<bool>& onStack) {
    Stack.push(current); 
    onStack[current] = true; 
    ids[current] = id; 
    low[current] = id; 
    id++; 
    
    for (int n : graph[current]) {
        if (ids[n] == -1) {
            dfs(n, graph, id, ids, low, Stack, onStack); 
        }
        if (onStack[n]) {
            low[current] = min(low[current], low[n]); 
        }
    }
    
    if (ids[current] == low[current]) {
        while (true) {
            int node = Stack.top(); 
            Stack.pop(); 
            onStack[node] = false;  
            if (node == current) break;
            low[node] = ids[current];
        }
    }
}

int main() {

    vector<vector<int>> graph = {
        {1, 4}, 
        {5}, 
        {1, 3, 6}, 
        {6}, 
        {0, 5}, 
        {2, 6}, 
        {7}, 
        {3}
    };

    unordered_map<int, vector<int>> scc = tarjan(graph);

    for (auto& [key, value] : scc) {
        cout << "SCC with low link " << key << ":\n";
        for (int node : value) {
            cout << node << " "; 
        }
        cout << endl; 
    }

    return 0; 
}