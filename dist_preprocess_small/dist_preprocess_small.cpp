#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
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
    int estimate = INFIN;
    // Lists of edges outgoing from each node
    std::vector<std::vector<std::pair<int, int>>> outgoing_edges;
    // Lists of edges incoming to each node
    std::vector<std::vector<std::pair<int, int>>> incoming_edges;

    static constexpr int INFIN = std::numeric_limits<int>::max() / 2;
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
        bidistance.resize(2, std::vector<int>(N, INFIN));
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
        std::vector<std::pair<int, int>> queue;
        for (int node = 0; node < N; node++) {
            queue.push_back({0, node}); 
        }
        struct min_heap {
            inline bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) {
                return a.first > b.first;
            }
        };

        // Priority queue will store pairs of (importance, node) with the least important node in the head
        // Queue queue(temp.begin(), temp.end());

        int nodeRank = 0; 

        // Implement the rest of the algorithm yourself
        // get the top of the queue
        while (queue.size() > 0) {
            auto top = queue.front(); 
            std::pop_heap(queue.begin(), queue.end(), min_heap());
            queue.pop_back(); 
            std::cout << "node : " << top.second + 1<< std::endl;
            // recompute the importance of the node
            std::vector<Shortcut> nodeShortcuts; 
            top.first = do_shortcut(top.second, nodeShortcuts, level[top.second]); 
            // compare the new importance to the node in the queue
            if (!queue.empty() && top.first > queue.front().first) {
                queue.push_back(top);
                std::push_heap(queue.begin(), queue.end(), min_heap());
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
            std::cout << "contracted : " << top.second + 1 << std::endl;
            contracted[top.second] = true; 
            rank[top.second] = nodeRank; 
            nodeRank++; 
        }

        for (int node = 0; node < N; node++) {
            std::cout << "rank of node " << node + 1 << " is " << rank[node] << std::endl;
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
        bidistance[0] = std::vector<int>(N, INFIN);
        bidistance[1] = std::vector<int>(N, INFIN);

        update(u, 0, true);
        update(w, 0, false);

        estimate = INFIN; 

        while (!diqueue[0].empty() || !diqueue[1].empty()) {
            while (!diqueue[0].empty()) {
                auto current = diqueue[0].pop().second;
                std::cout << "forward current: " << current << std::endl;
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
                std::cout << "backward current: " << current << std::endl;
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

        if (estimate == INFIN) {
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
        // int max_out = 0; 
        std::unordered_map<int, int> neighbors; 
        for (std::pair<int, int>& e : outgoing_edges[v]) {
            if (contracted[e.first]) cn++; 
            else neighbors[e.first] = e.second;
            // update the max out if the current edge is longer
            // if (e.second > max_out) {
            //     max_out = e.second;
            // }
        }
        for (Shortcut& s : outgoing_shortcuts[v]) {
            if (contracted[s.to]) cn++; 
            else neighbors[s.to] = s.cost;
            // update the max out if the current shortcut is longer
            // if (s.cost > max_out) {
            //     max_out = s.cost;
            // }
        }

        std::cout << "num neighbors: " << neighbors.size() << std::endl;

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
            // int witness_max = cost + max_out; 
            
            append_shortcuts(u, cost, v, shortcuts, neighbors);
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
            // int witness_max = cost + max_out;

            append_shortcuts(u, cost, v, shortcuts, neighbors); 
        }

        sc = shortcuts.size(); 
        ed += sc; 

        std::cout << "sc: " << sc << " cn: " << cn << " ed: " << ed << " L: " << mylevel << std::endl;

        // Add neighbors and shortcut cover heuristics
        return ed + cn + sc + mylevel;
    }

    void append_shortcuts(int u, int cost, int v, std::vector<Shortcut>& shortcuts, std::unordered_map<int, int>& neighbors) {
        int num_edges = outgoing_edges[v].size();
        int num_shortcuts = outgoing_shortcuts[v].size();

        // the nodes for which there are witnesses
        std::unordered_set<int> witness_found; 
        find_witness(u, cost, v, 3, neighbors, witness_found);
        // append shortcuts if there are no witness paths 
        for (auto& e : outgoing_edges[v]) {
            if (u != e.first && !contracted[e.first] && witness_found.find(e.first) == witness_found.end()) {
                std::cout << "adding shortcut from " << u + 1 << " to " << e.first + 1 << std::endl;
                shortcuts.push_back({u, e.first, cost + e.second});
            }
        }
        for (auto& s : outgoing_shortcuts[v]) {
            if (u != s.to && !contracted[s.to] && witness_found.find(s.to) == witness_found.end()) {
                std::cout << "adding shortcut from " << u + 1 << " to " << s.to + 1 << std::endl;
                shortcuts.push_back({u, s.to, cost + s.cost});
            }
        }
    }

    void find_witness(int u, int cost, int v, int number_of_hops, std::unordered_map<int, int>& neighbors, 
                      std::unordered_set<int>& witness_found) {
        // perform dijkstra ignoring the node v to find the witness paths
        // if the path is shorter than the length through v, then there is a witness path
        int num_edges = outgoing_edges[v].size();
        int num_shortcuts = outgoing_shortcuts[v].size();

        std::vector<bool> visited(N, false); 
        std::vector<int> distance(N, INFIN);
        std::vector<int> hops(N, 0); 
        Queue queue; 

        distance[u] = 0; 
        queue.push({0, u}); 

        while (!queue.empty()) {
            // first : distance second : node
            auto top = queue.top(); 
            auto current = top.second;
            queue.pop();

            if (visited[current]) continue; 

            // check whether we are at any of the outgoing nodes 
            if (current != u && neighbors.find(current) != neighbors.end()) {
                std::cout << "found witness path from node " << u + 1 << " to " << current + 1 << std::endl;
                // if the distance ignoring v is less than the distance to the neighbor, then there is a witness path
                if (top.first < cost + neighbors[current]) {
                    witness_found.insert(current); 
                }
                continue; 
            }

            // if the number of hops is more than max, then stop witness search
            if (hops[current] > number_of_hops) {
                break;
            }

            // add all the neighbors to the queue
            for (auto& e : outgoing_edges[current]) {
                // ignore node v
                if (e.first == v)  continue;
                // relax the outgoing edge
                int d = top.first + e.second;
                if (d < distance[e.first]) {
                    distance[e.first] = d;
                    // update the number of hops
                    hops[e.first] = hops[current] + 1;
                    queue.push({d, e.first});
                }
            }
            for (auto& s : outgoing_shortcuts[current]) {
                // ignore node v
                if (s.to == v) continue;
                // relax the outgoing shortcut
                int d = top.first + s.cost;
                if (d < distance[s.to]) {
                    distance[s.to] = d;
                    // update the number of hops
                    hops[s.to] = hops[current] + 1;
                    queue.push({d, s.to});
                }
            }

            visited[current] = true; 
        }

        visited.clear();
        distance.clear(); 
        hops.clear(); 
    }

    void remove_unnecessary_edges() { 
        // iterate through all the nodes
        for (int node = 0; node < N; node++) {
            // remove all outgoing edges that do not increase in node importance
            for (auto iterator = outgoing_edges[node].begin(); iterator != outgoing_edges[node].end();) {
                // if rank of outgoing node is less than current, remove this edge
                if (rank[iterator->first] < rank[node]) {
                    std::cout << "removing outgoing edge from " << node + 1 << " to " << iterator->first + 1 << std::endl;
                    iterator = outgoing_edges[node].erase(iterator);
                } else {
                    iterator++;
                }
            }
            // remove all outgoing shortcuts that do not increase in node importance
            for (auto iterator = outgoing_shortcuts[node].begin(); iterator != outgoing_shortcuts[node].end();) {
                // if rank of outgoing node is less than current, remove this edge
                if (rank[iterator->to] < rank[node]) {
                    std::cout << "removing outgoing shortcut from " << node + 1 << " to " << iterator->to + 1 << std::endl;
                    iterator = outgoing_shortcuts[node].erase(iterator);
                } else {
                    iterator++;
                }
            }
            // remove all incoming edges that do not increase in node importance
            for (auto iterator = incoming_edges[node].begin(); iterator != incoming_edges[node].end();) {
                // if rank of incoming node is less than current, remove this edge
                if (rank[iterator->first] < rank[node]) {
                    std::cout << "removing incoming edge from " << node + 1 << " to " << iterator->first + 1 << std::endl;
                    iterator = incoming_edges[node].erase(iterator);
                } else {
                    iterator++;
                }
            }
            // remove all incoming shortcuts that do not increase in node importance
            for (auto iterator = incoming_shortcuts[node].begin(); iterator != incoming_shortcuts[node].end();) {
                // if rank of incoming node is less than current, remove this edge
                if (rank[iterator->from] < rank[node]) {
                    std::cout << "removing incoming shortcut from " << node + 1 << " to " << iterator->from + 1 << std::endl;
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
    std::cout  << "Ready" << std::endl;

    int t;
    assert(scanf("%d", &t) == 1);
    for (int i = 0; i < t; ++i) {
        int u, v;
        assert(scanf("%d %d", &u, &v) == 2);
        printf("%d\n", g.query(u-1, v-1));
    }
}
