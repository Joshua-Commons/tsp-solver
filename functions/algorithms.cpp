#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstring>
#include <set>
#include <iostream>
#include <unordered_map>
#include <queue>
#include <utility>
#include <random>

extern "C" {

    /**
     * @brief Structure to store the result of the algorithms.
     * 
     * @param cost The minimum cost of the TSP tour.
     * @param path Pointer to the array storing the order of nodes in the optimal path.
     */
    struct Result {
        double cost;
        int* path;
    };


    /**
     * @brief Solves the Traveling Salesperson Problem (TSP) using the Held-Karp algorithm.
     * 
     * The function computes the minimum cost of a TSP tour using dynamic programming and bit masking.
     * It returns both the minimum cost and the optimal path.
     * 
     * @param n The number of nodes in the TSP problem.
     * @param dist_matrix A 1D array of size n*n representing the distance matrix. The distance between
     *        node i and node j is stored at dist_matrix[i * n + j].
     * @return A structure containing the minimum cost and the optimal path.
     */
    Result held_karp(
        int n, 
        const double* dist_matrix
        ) {
        
        const double INF = std::numeric_limits<double>::infinity();
        
        // number of problem subsets (2^n)
        int num_subsets = 1 << n;

        // initialise dp table to store minimum cost to reach subset i ending at node j and parent to track previous node on the heap
        double** dp = new double*[num_subsets];
        int** parent = new int*[num_subsets];
        for (int i = 0; i < num_subsets; ++i) {
            dp[i] = new double[n];
            parent[i] = new int[n];
            for (int j = 0; j < n; ++j) {
                dp[i][j] = INF;
                parent[i][j] = -1;
            }
        }

        // base case starting at node 0
        dp[1][0] = 0;

        // fill dp table for all subsets and possible last nodes
        for (int subset = 1; subset < num_subsets; ++subset) {
            for (int last_node = 0; last_node < n; ++last_node) {
                
                // skip if last node is not in the subset
                if (!(subset & (1 << last_node))) continue;
                
                // calculate subset excluding last node
                int prev_subset = subset & ~(1 << last_node);
                for (int prev_node = 0; prev_node < n; ++prev_node) {

                    // skip if previous node not in subset
                    if (!(prev_subset & (1 << prev_node))) continue;

                    // calculate cost of new path
                    double new_cost = dp[prev_subset][prev_node] + dist_matrix[prev_node * n + last_node];
                    if (new_cost < dp[subset][last_node]) {
                        dp[subset][last_node] = new_cost;       // update dp table
                        parent[subset][last_node] = prev_node;  // track transition
                    }
                }
            }
        }

        // find minimum cost to return to start node
        double min_cost = INF;
        int end_node = -1;
        
        for (int last_node = 1; last_node < n; ++last_node) {
            double new_cost = dp[num_subsets - 1][last_node] + dist_matrix[last_node * n + 0];
            if (new_cost < min_cost) {
                min_cost = new_cost;    
                end_node = last_node;
            }
        }

        // exit if there is no valid path
        if (end_node == -1) {
            std::cerr << "Error: No valid end_node found. Check your input distance matrix.\n";
            Result result;
            result.cost = INF;
            result.path = nullptr;
            return result;
        }

        // reconstruct optimal path
        int path[n];
        int current_subset = num_subsets - 1;
        int current_node = end_node;
        for (int i = n - 1; i >= 0; --i) {
            path[i] = current_node;                                 // record current node
            int next_node = parent[current_subset][current_node];   // get previous node
            current_subset &= ~(1 << current_node);                 // remove current node from subset
            current_node = next_node;                               // move to previous node
        }

        // allocate memory for the path to return it to Python
        int* result_path = new int[n];
        memcpy(result_path, path, n * sizeof(int));

        // free allocated memory
        for (int i = 0; i < num_subsets; ++i) {
            delete[] dp[i];
            delete[] parent[i];
        }
        delete[] dp;
        delete[] parent;

        // return result
        Result result;
        result.cost = min_cost;
        result.path = result_path;
        return result;
    }


    
    /**
     * @brief Generates the minimum spanning tree using Prim's algorithm.
     * 
     * The function computes the MST for a given graph represented as an adjacency matrix.
     * It uses a priority queue to efficiently select the next edge with the minimum weight.
     * The MST is returned as a list of edges (pairs of vertices).
     * 
     * @param graph The distance matrix. 
     * @return A vector of pairs, where each pair represents an edge in the MST.
     */
    std::vector<std::pair<int, int>> compute_mst(
        const std::vector<std::vector<double>>& graph
        ) {

        // intialise variables
        int n = graph.size();
        std::vector<bool> visited(n, false);
        std::vector<int> parent(n, -1);
        std::vector<double> min_edge(n, std::numeric_limits<double>::infinity());
        std::vector<std::pair<int, int>> mst_edges;

        // push first node into the queue 
        min_edge[0] = 0.0;
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
        pq.push({0.0, 0});

        // process the queue until all nodes are visited or queue is empty
        while (!pq.empty()) {

            // get node with minimum weight
            int u = pq.top().second;
            pq.pop();

            // if in visited set skip
            if (visited[u]) continue;

            // mark node visited
            visited[u] = true;

            // if it is not the starting node, add to MST
            if (parent[u] != -1) {
                mst_edges.emplace_back(parent[u], u);
            }

            // iterate over all neighbours of current node 
            for (int v = 0; v < n; ++v) {

                // if neighbour is unvisited and weight is less than minimum add to queue 
                if (!visited[v] && graph[u][v] < min_edge[v]) {
                    min_edge[v] = graph[u][v];
                    parent[v] = u;
                    pq.push({min_edge[v], v});
                }
            }
        }

        return mst_edges;
    }


    /**
     * @brief Finds all nodes with odd degree in MST
     * 
     * @param mst_edges Pairs of MST 
     * @param n Size of MST 
     * @return Nodes with odd degree in MST
     */
    std::vector<int> find_odd_degree_nodes(
        const std::vector<std::pair<int,int>>& mst_edges, 
        int n) {
        
        // initialise vector to store degrees of each node
        std::vector<int> degree(n, 0);

        // iteratively 1 degree to each pair of nodes
        for (const auto& edge : mst_edges) {
            degree[edge.first]++;
            degree[edge.second]++;
        }

        // vector to store nodes with odd degrees 
        std::vector<int> odd_degree_nodes;

        // add nodes with odd degrees
        for (int i = 0; i < n; ++i) {
            if (degree[i] % 2 == 1) {
                odd_degree_nodes.push_back(i);
            }
        }

        return odd_degree_nodes;
    }


    /**
     * @brief Computes a greedy minimum weight perfect matching for odd-degree nodes.
     * 
     * @param odd_nodes Vector of nodes with odd degrees.
     * @param graph Distance matrix.
     * @return A vector of pairs, where each pair represents a matched edge (u, v).
     */
    std::vector<std::pair<int, int>> compute_matching(
        const std::vector<int>& odd_nodes, 
        const std::vector<std::vector<double>>& graph
        ) {

        // create set of unmatched nodes
        std::set<int> unmatched(odd_nodes.begin(), odd_nodes.end());

        // vector to store edges of matching nodes
        std::vector<std::pair<int, int>> matching;

        // continue until all nodes are matched
        while (!unmatched.empty()) {

            // remove first unmatched node from set 
            int u = *unmatched.begin();
            unmatched.erase(u);

            // initialise variables to track best match
            double min_weight = std::numeric_limits<double>::infinity();
            int best_match = -1;

            // find closest match
            for (int v : unmatched) {
                if (graph[u][v] < min_weight) {
                    min_weight = graph[u][v];
                    best_match = v;
                }
            }

            // remove best match from set and store edge
            unmatched.erase(best_match);
            matching.emplace_back(u, best_match);
        }

        return matching;
    }


    /**
     * @brief Hierholzer's algorithm to find an Eulerian circuit.
     * 
     * @param graph Distance matrix.
     * @param start_node The node to start the Eulerian circuit.
     * @return Vector of integers representing the order of nodes in the Eulerian circuit.
     */
    std::vector<int> find_eulerian_circuit(
        const std::unordered_map<int, 
        std::multiset<int>>& graph, 
        int start_node) {

        // create copy of graph to modify
        std::unordered_map<int, std::multiset<int>> graph_copy = graph;

        // initialise vectors for circuit and current exploration path
        std::vector<int> circuit;
        std::vector<int> stack = {start_node};

        // process graph until all edges have been used
        while (!stack.empty()) {

            // get current node
            int u = stack.back();

            // if there are edges from the current node, traverse one
            if (!graph_copy[u].empty()) {

                // get neighbour
                int v = *graph_copy[u].begin();

                // move across it and back again
                graph_copy[u].erase(graph_copy[u].find(v));
                graph_copy[v].erase(graph_copy[v].find(u));

                // add to stack to continue traversal
                stack.push_back(v);
            
            // final node in the circuit
            } else {
                circuit.push_back(u);
                stack.pop_back();
            }
        }

        return circuit;
    }


    /**
     * @brief Shortcut Eulerian circuit to form Hamiltonian circuit
     * 
     * @param eulerian_circuit Vector of integers representing Eulerian circuit. 
     * @param n Number of nodes in graph.
     * @return Vector of integers representing the nodes in the Hamiltonian circuit, in order.
     */
    std::vector<int> shortcut_eulerian_circuit(
        const std::vector<int>& eulerian_circuit, 
        int n
        ) {

        // initialise set of visited nodes and new circuit
        std::vector<bool> visited(n, false);
        std::vector<int> hamiltonian_circuit;

        // rather than visit node twice, skip if it appears more than once
        for (int node : eulerian_circuit) {
            if (!visited[node]) {
                visited[node] = true;
                hamiltonian_circuit.push_back(node);
            }
        }

        return hamiltonian_circuit;
    }


    /**
     * @brief Solves the Traveling Salesperson Problem (TSP) using the Christofides algorithm.
     * 
     * The function computes the minimum cost of a TSP tour by generatign a  a minimum spaning tree 
     * and altering it such that it satisfies the constraints to be a TSP solution.
     * It returns both the minimum cost and the optimal path.
     * 
     * @param n The number of nodes in the TSP problem.
     * @param dist_matrix A 1D array of size n*n representing the distance matrix. The distance between
     *        node i and node j is stored at dist_matrix[i * n + j].
     * @return Structure containing the minimum cost and the optimal path.
     */
    std::vector<int> christofides(
        int n, 
        const double* dist_matrix
        ) {
        
        // Convert flat dist_matrix to 2D vector
        std::vector<std::vector<double>> graph(n, std::vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                graph[i][j] = dist_matrix[i * n + j];
            }
        }

        // compute MST
        auto mst_edges = compute_mst(graph);

        // find odd degree nodes
        auto odd_nodes = find_odd_degree_nodes(mst_edges, n);

        // compute minimum weight perfect matching
        auto matching = compute_matching(odd_nodes, graph);

        // combine MST and matching to form Eulerian graph
        std::unordered_map<int, std::multiset<int>> eulerian_graph;
        for (const auto& edge : mst_edges) {
            eulerian_graph[edge.first].insert(edge.second);
            eulerian_graph[edge.second].insert(edge.first);
        }
        for (const auto& edge : matching) {
            eulerian_graph[edge.first].insert(edge.second);
            eulerian_graph[edge.second].insert(edge.first);
        }

        // find Eulerian circuit
        auto eulerian_circuit = find_eulerian_circuit(eulerian_graph, 0);

        // shortcut to Hamiltonian circuit
        auto hamiltonian_circuit = shortcut_eulerian_circuit(eulerian_circuit, n);

        return hamiltonian_circuit;
    }


    /**
     * @brief Improves an existing solution to a TSP by performing 3-opt swaps.
     * 
     * The func
     * 
     * @param initial_tour The initial solution to the 
     * @param dist_matrix A 1D array of size n*n representing the distance matrix. The distance between
     *        node i and node j is stored at dist_matrix[i * n + j].
     * @return Structure containing the minimum cost and the optimal path.
     */
    std::vector<int> three_opt(
        const std::vector<int>& initial_tour, 
        const std::vector<std::vector<double>>& dist_matrix,
        const int max_iter
        ) {

        // get intial tour size and initialise new working tour
        int n = initial_tour.size();
        std::vector<int> tour = initial_tour;

        // initialise iteration variables
        bool improved = true;
        int iter = 0;
        constexpr double EPS = 1e-6;

        // iterate until local (or theoretically global) minimum found
        while (improved && (iter < max_iter)) {
            improved = false;
            iter++;

            // iterate through each combination of 3 connections
            for (int i = 0; i < n - 2; ++i) {
                for (int j = i + 1; j < n - 1; ++j) {
                    for (int k = j + 1; k < n; ++k) {

                        // calculate current cost of connections
                        double current_cost = dist_matrix[tour[i]][tour[i + 1]] + 
                                              dist_matrix[tour[j]][tour[j + 1]] + 
                                              dist_matrix[tour[k]][tour[(k + 1) % n]];

                        // calculate cost of reconnection scenarios
                        double new_cost1 = dist_matrix[tour[i]][tour[j]] + 
                                           dist_matrix[tour[i + 1]][tour[k]] + 
                                           dist_matrix[tour[j + 1]][tour[(k + 1) % n]];
                        
                        double new_cost2 = dist_matrix[tour[i]][tour[k]] + 
                                           dist_matrix[tour[j + 1]][tour[i + 1]] + 
                                           dist_matrix[tour[j]][tour[(k + 1) % n]];
                        
                        double new_cost3 = dist_matrix[tour[i]][tour[j + 1]] + 
                                           dist_matrix[tour[k]][tour[i + 1]] + 
                                           dist_matrix[tour[j]][tour[(k + 1) % n]];

                        // apply best improvement
                        if (new_cost1 < current_cost - EPS) {
                            std::reverse(tour.begin() + i + 1, tour.begin() + j + 1);
                            std::reverse(tour.begin() + j + 1, tour.begin() + k + 1);
                            improved = true;
                        
                        } else if (new_cost2 < current_cost - EPS) {
                            std::reverse(tour.begin() + i + 1, tour.begin() + k + 1);
                            improved = true;
                        
                        } else if (new_cost3 < current_cost - EPS) {
                            std::reverse(tour.begin() + j + 1, tour.begin() + k + 1);
                            improved = true;
                        }

                        // exit early if improved
                        if (improved) break;
                    }
                    if (improved) break;
                }
                if (improved) break;
            }
        }

        // std::cout << "Info: Exited on iteration: " << iter << std::endl;
        return tour;
    }


    /**
     * @brief Calculates the cost of a solution to a TSP problem.
     * 
     * @param tour The solution of the TSP problem showing the order of visited nodes.
     * @param dist_matrix A 1D array of size n*n representing the distance matrix. The distance between
     *        node i and node j is stored at dist_matrix[i * n + j].
     * @return The cost of the solution.
     */
    double calculate_tour_cost(
        const std::vector<int>& tour, 
        const std::vector<std::vector<double>>& dist_matrix
        ) {
        
        // initialises cost
        double cost = 0.0;

        // iterates through each connection in tour and adds distances together
        int n = tour.size();
        for (int i = 0; i < n; ++i) {
            cost += dist_matrix[tour[i]][tour[(i + 1) % n]];
        }

        // returns cost
        return cost;
    }
        


    /**
     * @brief Solve the Traveling Salesperson Problem (TSP) using 3-opt swap until a local maximum is found.
     * 
     * @param n The number of nodes in the TSP problem.
     * @param dist_matrix A 1D array of size n*n representing the distance matrix. The distance between
     *        node i and node j is stored at dist_matrix[i * n + j].
     * @return Result A structure containing the minimum cost and the optimal path.
     */
    Result solve_tsp(int n, const double* dist_matrix) {
        
        // if size is small enough, solve with dynamic programming
        if (n <= 23) {
            return held_karp(n, dist_matrix);
        }
        
        // convert flat distance matrix to 2D vector
        std::vector<std::vector<double>> graph(n, std::vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                graph[i][j] = dist_matrix[i * n + j];
            }
        }

        // generate initial solution with Christofides
        auto christofides_tour = christofides(n, dist_matrix);
        // std::cout << "Info: Christofides cost: " << christofides_cost << std::endl;

        // refine with 3-opt 
        int max_iters = 100000;
        auto three_opt_tour = three_opt(christofides_tour, graph, max_iters);
        double three_opt_cost = calculate_tour_cost(three_opt_tour, graph);
        // std::cout << "Info: 3-Opt refined cost: " << three_opt_cost << std::endl;

        // allocate result path for python
        int* result_path = new int[n];
        std::copy(three_opt_tour.begin(), three_opt_tour.end(), result_path);

        // return result
        Result result;
        result.cost = three_opt_cost;
        result.path = result_path;
        return result;
    }

}
