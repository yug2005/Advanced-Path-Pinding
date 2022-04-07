# Contraction Hierarchies

[TOC]

## Highway Hierarchies

In road networks, long distance trips usually go through long highways. 

* To get from A to B, you first merge into a highway, then into a bigger highway. 
* Then once you are close to your destination, you exit to a highway, then exit to a street, then go to B.
* Less important roads often merge into more important roads. There is a hierarchy of roads. 
* Many shortest path involve important nodes. Some nodes are also unavoidable such as nodes on a bridge. 

There are algorithms based on this idea such as “Highway Hierarchies” and “Transit Node Routing” by Sanders and Schultes.

These algorithms can be millions of times faster than Dijkstra.

Contraction hierarchies can be thousands of times faster than Dijkstra algorithm. 

### Node Ordering 

Nodes can be ordered by some “importance”.

Importance first increases, then decreases back along any shortest path.

We can also use bidirectional search. 

### Preprocessing

For contraction heirarchies, we need to preprocess the graph. 

When we are finding the shortest path, we use the preprocessed graph. 

Then reconstruct the shortest path in the initial graph. 

## Node Contraction

In the preprocessing, we eliminate nodes one by one in some order

We also add shortcuts to preserve distances. 

The output is an augmented graph with node ordering. 

We **contract** the nodes in the graph through some order. 
In the graph above, the nodes will be contracted in the order they are numbered. 

When you contract a node, you move the node to the augmented graph, and create a new shortcut edge. 
The green nodes are the contracted nodes and the blue edges are the new shortcut edges added to the graph. 

### Witness Paths

When we are performing the contraction of node $v$, we need to determine whether the path through that node is the shortest path between its neighbor nodes. Assuming that the shortest path from $u_1$ to $w_1$ is the path through the node $v$. Then we contract the node $v$, we need to add a new shortcut edge $(u_1, w_1)$. 
However, we do not always add shortcuts. Consider that the shortest path from $u_2$ to $w_2$ is actually not through the node $v$ and it is instead a different path. 
In this case, we do not need to add a new shortcut edge from the $u_2$ to $w_2$. 

* The shortest path is referred to as a **witness path**. 
* We only add shortcuts when there is no witness path $P_{uw}$ shorter than $\ell(u, v) + \ell(v, w)$ and bypassing $v$. 
* That is only add the shortcut $(u, v)$ for $\ell(u, v) + \ell(v, w)$

## Witness Search

When contracting node $v$, for any pair of edges $(u, v)$ and $(v, w)$ we want to check whether there is a witness path from $u$ to $w$ bypassing $v$ with length at most $\ell(u,v)+\ell(v,w)$, then there is no need to add a shortcut from $u$ to $w$.

**Witness search** is the search for witness paths. 

**Definition** : If there is an edge $(u,v)$, call u a **predecessor** of $v$. If there is an edge $(v,w)$, call $w$ a **successor** of $v$.

For each predecessor $u_i$ of $v$, run Dijkstra from $u_i$ ignoring $v$. This will find the witness path and help us avoid the extra shortcut nodes. 

To optimize the witness search, we can do either of the following : 

* We can stop Dijkstra when the distance from the source becomes too big.
* Limit the number of hops / traversals over edges

### Optimizing Witness Search

To stop Dijkstra at the appropriate time, we can ensure that the distance from the starting node to some node $x$ is always less than the max length of the any path from the predecessor to a successor of $v$. 

If $d(u_i, x) > \mathrm{max_{u, v}}( \ell(u, v) + \ell(v, w) )$, then there is no witness path going through $x$. 
Limiting the number of traversals through edges is another method of optimization. 

Consider only shortest paths from source with at most $k$ edges. 

* If there is no witness path found, then we add a shortcut. 
* There is a tradeoff between the preprocessing time and the augmented graph size. 
* If $k$ is small, then we won't find many witness paths and add more shortcuts. 

## Query

Once the graph is processed and the nodes are contracted, we can use bidirectional Dijkstra to find the shortest path. 

* The later the node is contracted, the more important the node is. It is higher in the contracted graph. 
* In the forward search and the backward search, the search is always going to nodes with higher importance. 

Bidirectional Dijkstra only uses **upwards** edges. 

* Don't stop the search when some node was processed by both the forward search and the backward search. 
* Stop Dijkstra only when the extracted node is already farther than the target. 

> Note you only need to store the edges in the graph that are directed upwards. 

Pseudocode for the query : 

```
estimate = infinity
Fill dist, dist^R with infinity for each node
dist[s] = 0, dist^R[t] = 0
proc = empty, proc^R = empty
while there are nodes to process:
	v = ExtractMin(dist)
	if dist[v] 
		Process(v, ...)
	if v in proc^R and dist[v] + dist^R[v] 
		estimate = dist[v] + dist^R[v]
	v^R = extractMin(dist^R)
	Repeat symmetrically for v^R
return estimate
```

## Node Ordering

The ordering of node contraction heavily influences the preprocessing and the query time. 

* We want to minimize the number of added shortcuts. 
* We want to spread the important nodes across the graph. 
* We want to minimize the number of edges in the shortest paths in the augmented graph. 

### Order by Importance

We need to introduce a measure of importance and contract the least important nodes first. 

However, when contracting the nodes, the importance can change. 

We will keep all the nodes in a priority queue by decreasing importance and the head of the queue will contain the least important node. 

* On each iteration, we will extract the least important node. However, we need to recompute its importance. 
* If it is still minimal when compared to the top of the priority queue, then we can contract the node. 
* Otherwise, we need to put it back into the priority queue with the new priority. 

### Importance Criteria

The node importance is based on the following factors : 

* Edge difference
* Number of contracted neighbors 
* Shortcut cover
* Node level

#### Edge Difference

* Want to minimize the number of edges in the augmented graph. 
* The number of added shortcuts $s(v)$, incoming degree $\mathrm{in}(v)$, outgoing degree $\mathrm{out}(v)$. 
* The **edge difference** is given by $\mathrm{ed}(v) = s(v) - \mathrm{in}(v) - \mathrm{out}(v)$.
* The number of edges increases by $\mathrm{ed}(v)$ after contracting the node $v$. 
* Contract node with small $\mathrm{ed}(v)$ first. 

#### Contracted Neighbors

* Want to spread contracted nodes across the graph. 
* Contract a node with small number of already contracted neighbors $\mathrm{cn}(v)$. 

#### Shortcut Cover

* Want to contract important nodes later. 
* Shortcut cover $\mathrm{sc}(v)$ : The number of neighbors $w$ of $v$ such that we have to add shortcut to or from $w$ after contracting $v$. 
* If shortcut cover is big, many nodes "depend" on the node $v$. 
* Contract a node with small $\mathrm{sc}(v)$. 

#### Node Level

* Node level $L(v)$ is an upper bound on the number of edges in the shortest path from $s$ to $v$ in the augmented graph. Initially, it is equal to 0 for node $v$. 
* After contracting node $v$, for neighbors $u$ of $v$, do $L(u) = \mathrm{max}(L(u), L(v) + 1)$. 
* Contract a node with small $L(v)$. 

### Importance Value of Node

We use the sum of the importance criteria above to measure the importance of a node in the graph. 
$$
\mathrm I(v) = \mathrm{ed}(v) + \mathrm{cn}(v) + \mathrm{sc}(v) + L(v)
$$
You can experiment with different weights of the four quantities to see how the preprocessing time and the query time changes. 

* Each of the four quantities are necessary for fast preprocessing and queries. 
* Find a way to compute them efficiently at any stage of the preprocessing. 
