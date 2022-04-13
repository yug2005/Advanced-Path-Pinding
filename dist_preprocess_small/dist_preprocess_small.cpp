#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
#include <queue>
#include <iostream>
#include <memory>
#include <cassert>

// cd dist_preprocess_small
// g++ -pipe -O2 -std=c++14 dist_preprocess_small.cpp -lm

class Graph
{
    typedef int Distance;
    typedef int Node;

    // Number of nodes
    int N;
    // Source and target
    int s, t;
    // Estimate of the distance from s to t
    int estimate = INFINITY;
    // Lists of edges outgoing from each node
    std::vector<std::vector<std::pair<int, int>>> outgoing_edges;
    // Lists of edges incoming to each node
    std::vector<std::vector<std::pair<int, int>>> incoming_edges;

    static constexpr int INFINITY = std::numeric_limits<int>::max() / 2;
    // Levels of nodes for node ordering
    std::vector<int> level;
    // Ranks of nodes - positions in the node ordering
    std::vector<int> rank;
    // Whether a node has been contracted
    std::vector<bool> contracted;

    // Distance to node v, bidistance[0][v] - from source in the forward search, bidistance[1][v] - from target
    // in the backward search.
    std::vector<std::vector<Distance>> bidistance;

    // Wrapper around STL priority_queue
    class StlHeap
    {
    public:
        using T = std::pair<Distance, Node>;
        using Queue = std::priority_queue<T, std::vector<T>, std::greater<T>>;

        StlHeap() {
            queue.reset(new Queue());
        }

        bool empty() const {
            return queue->empty();
        }

        void update(Node v, Distance d) {
            queue->push(std::make_pair(d,v));
        }

        void clear() {
            queue.reset(new Queue());
        }

        std::pair<Distance, Node> pop() {
            std::pair<Distance, Node> top = queue->top();
            queue->pop();
            return top;
        }

    private:
        std::unique_ptr<Queue> queue;
    };

    // Priority queues for forward and backward searches
    StlHeap diqueue[2];
public:
    Graph() {
        read_stdin();
        bidistance.resize(2, std::vector<int>(N, INFINITY));
    }

    int get_n() { return N;}

    std::vector<std::pair<int, int>>& get_adjacent(int v, bool forward = true) {
        if (forward) {
            return outgoing_edges[v];
        } else {
            return incoming_edges[v];
        }
    }

    void preprocess() {
        // resize the level and rank vector
        rank.resize(N);
        level.resize(N, 0);
        contracted.resize(N, false); 
        incoming_shortcuts.resize(N);
        outgoing_shortcuts.resize(N);
        // temp vector to add nodes to queue
        std::vector<std::pair<int, int>> temp;
        for (int node = 0; node < N; node++) {
            temp.push_back({0, node}); 
        }

        // Priority queue will store pairs of (importance, node) with the least important node in the head
        Queue queue(temp.begin(), temp.end());

        int nodeRank = 0; 

        // Implement the rest of the algorithm yourself
        // get the top of the queue
        while (queue.size() > 0) {
            auto top = queue.top();
            queue.pop(); 
            
            // recompute the importance of the node
            std::vector<Shortcut> nodeShortcuts; 
            top.first = do_shortcut(top.second, nodeShortcuts, level[top.second]); 
            // compare the new importance to the node in the queue
            if (!queue.empty() && top.first > queue.top().first) {
                queue.push(top);
                nodeShortcuts.clear(); 
                continue; 
            }

            // if the top node has least importance, then contract it
            // add the shortcuts
            for (Shortcut& s : nodeShortcuts) {
                // remove any outgoing edges that are the same as the shortcut 
                for (auto iterator = outgoing_edges[s.from].begin(); iterator != outgoing_edges[s.from].end();) {
                    if (iterator->first == s.to) {
                        iterator = outgoing_edges[s.from].erase(iterator);
                    }
                    else iterator++;
                }
                // remove any incoming edges that are the same as the shortcut
                for (auto iterator = incoming_edges[s.to].begin(); iterator != incoming_edges[s.to].end();) {
                    if (iterator->first == s.from) {
                        iterator = incoming_edges[s.to].erase(iterator);
                    }
                    else iterator++;
                }
                // add the shortcut to the graph
                outgoing_shortcuts[s.from].push_back(s);
                incoming_shortcuts[s.to].push_back(s);
            }
            // update the level of the neighbor nodes once v is contracted
            for (auto& e : incoming_edges[top.second]) {
                level[e.first] = std::max(level[e.first], level[top.second] + 1);
            }
            for (auto& e : outgoing_edges[top.second]) {
                level[e.first] = std::max(level[e.first], level[top.second] + 1);
            }
            // mark the node as contracted
            contracted[top.second] = true; 
            rank[top.second] = nodeRank; 
            nodeRank++; 
        }

        // remove the edges and shortcuts that do not increase in node importance
        remove_unnecessary_edges();
    }

    // Returns distance from s to t in the graph
    int query(int u, int w) {
        if (u == w) return 0; 
        s = u; 
        t = w;

        std::vector<std::vector<bool>> processed(2, std::vector<bool>(N, false));
        bidistance[0] = std::vector<int>(N, INFINITY);
        bidistance[1] = std::vector<int>(N, INFINITY);

        update(u, 0, true);
        update(w, 0, false);

        estimate = INFINITY; 

        while (!diqueue[0].empty() || !diqueue[1].empty()) {
            while (!diqueue[0].empty()) {
                auto current = diqueue[0].pop().second;
                if (processed[0][current]) break;
                if (bidistance[0][current] < estimate) {
                    // add the outgoing edges to the queue
                    for (auto& e : outgoing_edges[current]) {
                        update(e.first, bidistance[0][current] + e.second, true);
                    }
                    // add the outgoing shortcuts to the queue
                    for (auto& s : outgoing_shortcuts[current]) {
                        update(s.to, bidistance[0][current] + s.cost, true);
                    }
                }
                processed[0][current] = true; 
                if (processed[1][current] && bidistance[0][current] + bidistance[1][current] < estimate) {
                    estimate = bidistance[0][current] + bidistance[1][current];
                }
                break; 
            }
            while (!diqueue[1].empty()) {
                auto current = diqueue[1].pop().second; 
                if (processed[1][current]) break;
                if (bidistance[1][current] < estimate) {
                    // add the incoming edges to the queue
                    for (auto& e : incoming_edges[current]) {
                        update(e.first, bidistance[1][current] + e.second, false);
                    }
                    // add the incoming shortcuts to the queue
                    for (auto& s : incoming_shortcuts[current]) {
                        update(s.from, bidistance[1][current] + s.cost, false);
                    }
                }
                processed[1][current] = true;
                if (processed[0][current] && bidistance[0][current] + bidistance[1][current] < estimate) {
                    estimate = bidistance[0][current] + bidistance[1][current];
                }
                break; 
            }
        }

        if (estimate == INFINITY) {
            return -1;
        } else {
            return estimate;
        }
    }

private:
    // Try to relax the node v using distance d either in the forward or in the backward search
    void update(int v, int d, bool forward) {
        // Implement this method yourself
        const int i = forward ? 0 : 1; 
        if (d < bidistance[i][v]) {
            bidistance[i][v] = d;
            diqueue[i].update(v, d);
        }
    }

    class VertexSet
    {
    public:
        VertexSet(int n = 0) : visited(n) {}
        void resize(int n) {
            visited.resize(n);
        }
        void add(int v) {
            if (!visited[v]) {
                vertices.push_back(v);
                visited[v] = true;
            }
        }
        const std::vector<int>& get() const {
            return vertices;
        }
        const bool has(int v) {
            return visited[v];
        }
        void clear() {
            for (int v : vertices) {
                visited[v] = false;
            }
            vertices.clear();
        }

    private:
        std::vector<int> visited;
        std::vector<int> vertices;
    };
    VertexSet visited;

    // QEntry = (distance, vertex)
    typedef std::pair<int,int> QEntry;
    typedef std::priority_queue<QEntry, std::vector<QEntry>, std::greater<QEntry>> Queue;

    struct Shortcut {
        int from;
        int to;
        int cost;
    };
    // Lists of shortcuts outgoing from each node
    std::vector<std::vector<Shortcut>> outgoing_shortcuts;
    // Lists of shortcuts incoming to each node
    std::vector<std::vector<Shortcut>> incoming_shortcuts;

    // Adds all the shortcuts for the case when node v is contracted, and returns the importance of node v
    // in this case
    int do_shortcut(int v, std::vector<Shortcut>& shortcuts, int& mylevel) {
        // Implement this method yourself
        shortcuts.clear();
        // importance criteria
        int cn = 0, sc = 0; 
        int ed = sc - outgoing_edges[v].size() - incoming_edges[v].size() 
                    - outgoing_shortcuts[v].size() - incoming_shortcuts[v].size();
        // if there are no incoming or outgoing edges, then return 
        if ((incoming_edges[v].empty() && incoming_shortcuts[v].empty()) || 
            (outgoing_edges[v].empty() && outgoing_shortcuts[v].empty())) {
            // ? What do you do when there are no incoming or outgoing edges?
            return ed + cn + sc + mylevel; 
        }
        // find max outgoing edge or shortcut for stoping witness search
        int max_out = 0; 
        for (std::pair<int, int>& e : outgoing_edges[v]) {
            if (contracted[e.first]) cn++; 
            // update the max out if the current edge is longer
            if (e.second > max_out) {
                max_out = e.second;
            }
        }
        for (Shortcut& s : outgoing_shortcuts[v]) {
            if (contracted[s.to]) cn++; 
            // update the max out if the current shortcut is longer
            if (s.cost > max_out) {
                max_out = s.cost;
            }
        }

        // witness search
        // for every predecessor edge of v (u, v), add shortcuts if appropriate
        for (std::pair<int, int>& edge : incoming_edges[v]) {
            int u = edge.first;
            // if neighbor u is contracted, add to contracted neighbors
            if (contracted[u]) {
                cn++; 
                continue; 
            }
            int cost = edge.second;
            // find the max length of the witness path
            int witness_max = cost + max_out; 
            
            append_shortcuts(u, cost, v, witness_max, shortcuts);
        }

        // for every predecessor shortcut of v (s.from, v), add shortcuts if appropriate
        for (Shortcut& shortcut : incoming_shortcuts[v]) {
            int u = shortcut.from;
            // if neighbor u is contracted, add to contracted neighbors
            if (contracted[u]) {
                cn++; 
                continue; 
            }
            int cost = shortcut.cost;
            // find the max length of the witness path
            int witness_max = cost + max_out;

            append_shortcuts(u, cost, v, witness_max, shortcuts); 
        }

        sc = shortcuts.size(); 
        ed += sc; 

        // Add neighbors and shortcut cover heuristics
        return ed + cn + sc + mylevel;
    }

    void append_shortcuts(int u, int cost, int v, int witness_max, std::vector<Shortcut>& shortcuts) {
        int num_edges = outgoing_edges[v].size();
        int num_shortcuts = outgoing_shortcuts[v].size();
        // find out whether there are witness paths
        std::vector<bool> witness_found = find_witness({u, cost}, v, witness_max, 5);
        // append shortcuts if there are no witness paths 
        for (int index = 0; index < num_edges; index++) {
            auto e = outgoing_edges[v][index];
            // if neighbor is contracted, add to contracted neighbors and continue
            if (contracted[e.first] || witness_found[index]) {
                continue; 
            }
            // if there was no witness path found, then create a shortcut
            shortcuts.push_back({u, e.first, cost + e.second});
        }
        for (int index = num_edges; index < num_edges + num_shortcuts; index++) {
            auto s = outgoing_shortcuts[v][index - num_edges];
            // if neighbor is contracted, add to contracted neighbors and continue
            if (contracted[s.to] || witness_found[index]) {
                continue; 
            }
            // if there was no witness path found, then create a shortcut
            shortcuts.push_back({u, s.to, cost + s.cost});
        }
    }

    std::vector<bool> find_witness(std::pair<int, int> incoming, int v, int witness_max, int number_of_hops) {
        // perform dijkstra ignoring the node v to find the witness paths
        // if the path is shorter than the length through v, then there is a witness path
        int num_edges = outgoing_edges[v].size();
        int num_shortcuts = outgoing_shortcuts[v].size();

        std::vector<bool> witness_found(num_edges + num_shortcuts, false); 

        std::vector<int> distance(N, INFINITY);
        std::vector<int> hops(N, 0); 
        Queue queue; 
        queue.push({0, incoming.first}); 

        while (!queue.empty()) {
            // first : distance second : node
            auto top = queue.top(); 
            queue.pop();

            // check whether we are at any of the outgoing nodes
            for (int index = 0; index < num_edges; index++) {
                if (!witness_found[index]) {
                    auto e = outgoing_edges[v][index];
                    // if the path is shorter, then there is a witness path
                    if (top.second == e.first && top.first <= incoming.second + e.second){
                        witness_found[index] = true;
                    }
                }
            }
            // check whether we are at any of the outgoing shortcuts
            for (int index = num_edges; index < num_edges + num_shortcuts; index++) {
                if (!witness_found[index]) {
                    auto s = outgoing_shortcuts[v][index - num_edges];
                    // if the path is shorter, then there is a witness path
                    if (top.second == s.to && top.first <= incoming.second + s.cost) {
                        witness_found[index] = true;
                    }
                }
            }

            // if the distance is past the max length, then stop witness search
            // if the number of hops is more than max, then stop witness search
            if (top.first >= witness_max || hops[top.second] > number_of_hops) {
                break;
            }

            // add all the neighbors to the queue
            for (auto& e : outgoing_edges[top.second]) {
                // ignore node v
                if (e.first == v) {
                    continue;
                }
                // relax the outgoing edge
                int d = top.first + e.second;
                if (d < distance[e.first]) {
                    distance[e.first] = d;
                    // update the number of hops
                    hops[e.first] = hops[top.second] + 1;
                    queue.push({d, e.first});
                }
            }
        }

        distance.clear(); 
        hops.clear(); 

        return witness_found;
    }

    void remove_unnecessary_edges() { 
        // iterate through all the nodes
        for (int node = 0; node < N; node++) {
            // remove all outgoing edges that do not increase in node importance
            for (auto iterator = outgoing_edges[node].begin(); iterator != outgoing_edges[node].end();) {
                // if rank of outgoing node is less than current, remove this edge
                if (rank[iterator->first] < rank[node]) {
                    iterator = outgoing_edges[node].erase(iterator);
                } else {
                    iterator++;
                }
            }
            // remove all outgoing shortcuts that do not increase in node importance
            for (auto iterator = outgoing_shortcuts[node].begin(); iterator != outgoing_shortcuts[node].end();) {
                // if rank of outgoing node is less than current, remove this edge
                if (rank[iterator->to] < rank[node]) {
                    iterator = outgoing_shortcuts[node].erase(iterator);
                } else {
                    iterator++;
                }
            }
            // remove all incoming edges that do not increase in node importance
            for (auto iterator = incoming_edges[node].begin(); iterator != incoming_edges[node].end();) {
                // if rank of incoming node is less than current, remove this edge
                if (rank[iterator->first] < rank[node]) {
                    iterator = incoming_edges[node].erase(iterator);
                } else {
                    iterator++;
                }
            }
            // remove all incoming shortcuts that do not increase in node importance
            for (auto iterator = incoming_shortcuts[node].begin(); iterator != incoming_shortcuts[node].end();) {
                // if rank of incoming node is less than current, remove this edge
                if (rank[iterator->from] < rank[node]) {
                    iterator = incoming_shortcuts[node].erase(iterator);
                } else {
                    iterator++;
                }
            }
        }
    }

    void set_n(int n) {
        N = n;
        outgoing_edges.resize(n);
        incoming_edges.resize(n);
    }

    void add_edge_to_list(std::vector<std::pair<int,int>>& list, int w, int c) {
        for (int i = 0; i < list.size(); ++i) {
            std::pair<int, int>& p = list[i];
            if (p.first == w) {
                if (p.second > c) {
                    p.second = c;
                }
                return;
            }
        }
        list.push_back({w, c});
    }

    void add_directed_edge(int u, int v, int c) {
        add_edge_to_list(outgoing_edges[u], v, c);
        add_edge_to_list(incoming_edges[v], u, c);
    }

    void add_edge(int u, int v, int c) {
        add_directed_edge(u, v, c);
    }

    void finalize() {
        // Remove unnecessary edges
        // remove duplicates
        for (std::vector<std::pair<int, int>>& edges : outgoing_edges) {
            std::sort(edges.begin(), edges.end());
            edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        }
        for (std::vector<std::pair<int, int>>& edges : incoming_edges) {
            std::sort(edges.begin(), edges.end());
            edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        }
    }

    bool read_stdin() {
        int u,v,c,n,m;
        assert(scanf("%d %d", &n, &m) == 2);
        set_n(n);
        for (int i = 0; i < m; ++i) {
            assert(scanf("%d %d %d", &u, &v, &c) == 3);
            add_edge(u-1, v-1, c);
        }
        finalize();
        return true;
    }
};

int main() {
    Graph g;
    g.preprocess();
    std::cout << "Ready" << std::endl;

    int t;
    assert(scanf("%d", &t) == 1);
    for (int i = 0; i < t; ++i) {
        int u, v;
        assert(scanf("%d %d", &u, &v) == 2);
        printf("%d\n", g.query(u-1, v-1));
    }
}
