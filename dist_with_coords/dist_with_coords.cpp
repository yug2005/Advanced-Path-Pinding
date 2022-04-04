#include <cstdio>
#include <cassert>
#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <cmath>
#include <iostream>

using namespace std;

// See the explanations of these typedefs and constants in the starter for friend_suggestion
typedef vector<vector<vector<int>>> Adj;
typedef long long Len;
typedef vector<priority_queue<pair<Len, int>,vector<pair<Len,int>>,greater<pair<Len,int>>>> Queue;

const Len INFIN = numeric_limits<Len>::max() / 4LL;

class AStar {
    // See the descriptions of these fields in the starter for friend_suggestion
    int n_;
    Adj adj_;
    Adj cost_;
    vector<vector<Len>> distance_;
    vector<int> workset_;
    vector<bool> visited_;
    // Coordinates of the nodes
    vector<pair<Len,Len>> xy_;

public:
    AStar(int n, Adj adj, Adj cost, vector<pair<Len,Len>> xy)
        : n_(n), adj_(adj), cost_(cost), distance_(2, vector<Len>(n_, INFIN)), visited_(n), xy_(xy)
    { workset_.reserve(n); }

    // See the description of this method in the starter for friend_suggestion
    void clear() {
        for (int i = 0; i < workset_.size(); ++i) {
            int v = workset_[i];
            distance_[0][v] = distance_[1][v] = INFIN;
            visited_[v] = false;
        }
        workset_.clear();
    }

    // See the description of this method in the starter for friend_suggestion
    void visit(Queue& q, int side, int v, Len dist, int s, int t) {
        if (dist < distance_[side][v]) {
            distance_[side][v] = dist;
            Len h = getH(v, s, t); 
            q[side].push({(side == 0 ? dist + h : dist - h), v});
        }
    }

    Len getH(int v, int s, int t) {
        if (v == s || v == t) return 0; 
        Len dt = sqrt(pow(xy_[v].first - xy_[t].first, 2) + pow(xy_[v].second - xy_[t].second, 2));
        Len ds = sqrt(pow(xy_[v].first - xy_[s].first, 2) + pow(xy_[v].second - xy_[s].second, 2));
        return (dt - ds)/2; 
    }

    // Returns the distance from s to t in the graph
    Len query(int s, int t) {
        if (s == t) return 0; 
        
        clear();
        Queue q(2);
        visit(q, 0, s, 0, s, t);
        visit(q, 1, t, 0, s, t);

        pair<vector<bool>, vector<bool>> visited = {visited_, visited_};

        while (!q[0].empty() && !q[1].empty()) {
            int v = q[0].top().second;
            q[0].pop();
            if (!visited.first[v]) {
                for (int i = 0; i < adj_[0][v].size(); i++) {
                    visit(q, 0, adj_[0][v][i], distance_[0][v] + cost_[0][v][i], s, t);
                }
                if (visited.second[v]) return shortest_distance(s, t); 
                visited.first[v] = true; 
                workset_.push_back(v);
            }

            int w = q[1].top().second;
            q[1].pop();
            if (!visited.second[w]) {
                for (int i = 0; i < adj_[1][w].size(); i++) {
                    visit(q, 1, adj_[1][w][i], distance_[1][w] + cost_[1][w][i], s, t);
                }
                if (visited.first[w]) return shortest_distance(s, t);
                visited.second[w] = true; 
                workset_.push_back(w);
            }
        }
        return -1;
    }

    Len shortest_distance(int s, int t) {
        Len distance = INFIN; 
        for (int v : workset_) {
            if (distance_[0][v] + distance_[1][v] < distance) {
                distance = distance_[0][v] + distance_[1][v];
            }
        }
        return distance; 
    }
};

int main() {
    int n, m;
    scanf("%d%d", &n, &m);
    vector<pair<Len,Len>> xy(n);
    for (int i=0;i<n;++i){
        int a, b;
        scanf("%d%d", &a, &b);
        xy[i] = make_pair(a,b);
    }
    Adj adj(2, vector<vector<int>>(n));
    Adj cost(2, vector<vector<int>>(n));
    for (int i=0; i<m; ++i) {
        int u, v, c;
        scanf("%d%d%d", &u, &v, &c);
        adj[0][u-1].push_back(v-1);
        cost[0][u-1].push_back(c);
        adj[1][v-1].push_back(u-1);
        cost[1][v-1].push_back(c);
    }

    AStar astar(n, adj, cost, xy);

    int t;
    scanf("%d", &t);
    for (int i=0; i<t; ++i) {
        int u, v;
        scanf("%d%d", &u, &v);
        printf("%lld\n", astar.query(u-1, v-1));
    }
}
