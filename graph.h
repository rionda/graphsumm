#ifndef GRAPH_H
#define GRAPH_H

#include <array>
#include <cassert>

#include "util.h"

#define MAXLINE 1000

// General graph data (simple, undirected, unweighted)
int n, num_edges;                   // number of vertices and edges
ll triangles = -1;
vector<vector<vertex>> adj;         // adjacency lists; adj[i] is sorted in increasing order
vector<SparseRow> adj2;
double avg_deg, avg_dens;           // average degree/density

unordered_map<ll, int> vertices;    // consecutively numbered
vector<ll> vertex_labels;

bool verify_ids() {
    for (int i =0; i < vertex_labels.size(); i++) {
        if (vertices[vertex_labels[i]] != i) {
            return false;
        }
    }
    return true;
}

/*
 * Returns the integer representing a given vertex label.
 * If there is no integer for this label, assign one to it, in increasing
 * order.
 *
 */
vertex vertex_id(ll label) {
    auto it = vertices.find(label);
    if (it == vertices.end()) {
        int size = vertices.size();
        adj.push_back(vector<int>());
        vertex_labels.push_back(label);
        return vertices[label] = size;
    } else
        return it->second;
}

/*
 * Returns the next pair of labels in the input in_file
 *
 * Each line in the input file is either a comment if it starts with '#' or a
 * pair of integers, "a b", representing an edge between node 'a' and node 'b'.
 *
 */
bool next_edge(array<ll, 2>& e) {
    static FILE* in_file = NULL;               // input file

    // Open input file, or fall back to stdin
    if (!in_file) {
        if (!in_filename[0]) in_file = stdin;
        else in_file = fopen(in_filename, "r");
        if (in_file == NULL) {
            perror("Error opening input file");
            exit(errno);
        }
    }

    char line[MAXLINE], c;
    // Skip comments
    while ((c = fgetc(in_file)) == '#') {
        fgets(line, sizeof(line), in_file);
        if (strstr(line, "triangles:") != NULL) {
            char *start = strstr(line, ":") + 1;
            triangles = strtoll(start, NULL, 10);
            printf("Found triangles metadata: %lld\n", triangles);
        }
    }

    // Handle the case where the last line is a comment and there's no newline.
    if (feof(in_file)) {
        fclose(in_file);
        return false;
    }

    ungetc(c, in_file);

    char *a, *b;
    fgets(line, sizeof(line), in_file);
    a = strtok(line, " \t\r\n");
    b = strtok(NULL, " \t\r\n");

    if (a == NULL || b == NULL || strlen(a) == 0 || strlen(b) == 0) {
        if (!feof(in_file)) {
            fprintf(stderr, "Malformed line in input file: %s\n", line);
            exit(1);
        }
        fclose(in_file);
        return false;
    }

    e[0] = atoll(a);
    e[1] = atoll(b);
    return true;
}

void print_octave_matrix() {
    FILE* file = fopen("a.mat", "w");
    fprintf(file, "# name: A\n# type: sparse matrix\n# nnz: %i\n# rows: %i\n# columns: %i\n",
            2 * num_edges, n, n);
    for (vertex v = 0; v < n; ++v)
        for (vertex w : adj[v])
            fprintf(file, "%i %i 1\n", w + 1, v + 1);
    fclose(file);
/*
load a.mat; k = 400; [U,S,V]=svds(A,k); P=U*S*V'; norm(A-P,'fro')^2/prod(size(A))
    */
}

/*
 * Reads the graph from the input file in_file and fills in the general graph
 * data.
 *
 * Each line in the input file is either a comment if it starts with '#' or a
 * pair of integers, "a b", representing an edge between node 'a' and node 'b'.
 *
 */
void read_graph() {
    array<ll, 2> edge;
    num_edges = 0;
    while (next_edge(edge)) {
        vertex a = vertex_id(edge[0]), b = vertex_id(edge[1]);
        if (a == b) continue; // we assume no loops, although we could handle them.

        // Add neighbor to the adjacency lists of the edge end points
        // XXX We are assuming the graph is _undirected_
        adj[a].push_back(b);
        adj[b].push_back(a);
        ++num_edges;
    }
    printf("before duplicate removal: %d edges; ", num_edges);

    n = vertices.size();
    num_edges = 0;
    for (vertex i = 0; i < n; ++i) {
        sort(adj[i].begin(), adj[i].end());
        // Make undirected (remove duplicates)
        adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
        num_edges += adj[i].size();
    }
    assert(num_edges % 2 == 0);
    num_edges /= 2;
    printf("after duplicate removal: %d edges.\n", num_edges);

    // Compute a couple of simple statistics
    avg_deg = 2.0 * num_edges / n;
    avg_dens = avg_deg / n;

    print_octave_matrix();
}

/*
 * Is vertex i a neighbour of vertex j?
 *
 */
bool is_neighbour(vertex i, vertex j) {
    if (adj[i].size() > adj[j].size()) swap(i, j);
    return binary_search(adj[i].begin(), adj[i].end(), j);
}

template<class T>
T adjacency_row(int v) {
    T ret;
    for (vertex w : adj[v])
        ret[w] = 1;
    return ret;
}

/*
 * Return the full n x n adjacency matrix of the graph
 *
 * ret[v][w] == 1 iff there is an edge between v and w, 0 otherwise.
 *
 * Use example: full_adjacency_matrix<bool>() or full_adjacency_matrix<double>()
 */
template<typename T>
Matrix<T> full_adjacency_matrix() {
    Matrix<T> ret;
    ret.resize(n, vector<T>(n));
    for (vertex v = 0; v < n; ++v)
        for (vertex w : adj[v])
            ret[v][w] = 1;
    return ret;
}

/*
 * Return a random neighbor of vertex v
 *
 */
vertex random_neighbour(vertex v) {
    return adj[v][rand() % adj[v].size()];
}

/*
 * Return the list of common neighbours
 *
 */
vector<vertex> common_neighbours(vertex i, vertex j) {
    vector<vertex> ret;
    set_intersection(adj[i].begin(), adj[i].end(), adj[j].begin(), adj[j].end(), back_inserter(ret));
    return ret;
}

/*
 * Return the exact number of common neighbours
 *
 */
inline int num_common(vertex a, vertex b) {
    return common_neighbours(a, b).size();
    /*
    if (adj[a].size() > adj[b].size()) swap(a, b);

    int ret = 0;
    auto it = adj[b].begin();
    for (vertex v : adj[a]) {
        it = lower_bound(it, adj[b].end(), v);
        if (it == adj[b].end()) break;
        if (*it == v) ++ret;
    }
    return ret;
    */
}

inline int xor_neighbours(vertex i, vertex j) {
    return adj[i].size() + adj[j].size() - 2 * num_common(i, j);
}

/*
 * Compute the number of triangles in the graph
 *
 */
int triangles_number() {
    if (triangles == -1) {
        printf("Computing triangles...");
        triangles = 0;
        for (int i = 0; i < n; ++i)
            for (int j : adj[i]) {
                if (j >= i) break;

                for (int k : adj[i]) {
                    if (k >= j) break;

                    if (is_neighbour(j, k))
                        ++triangles;
                }
            }
        printf("done: %lld\n", triangles);
    } else {
        printf("Cached triangles: %lld\n", triangles);
    }

    return triangles;
}

/*
 * Compute the exact square of the adjacency matrix
 *
 */
void square_matrix() {
    double forever = (double)num_edges * n;
    double estimate = (double)num_edges * avg_deg;
    printf("m · n = %lg  rough estimate = %lg; %lf secs\n", forever, estimate, estimate / 1e6);

    adj2.resize(n);
    for (vertex i = 0; i < n; ++i)
        for (vertex k : adj[i])
            for (vertex j : adj[k])
                ++adj2[i][j];

#ifdef DEBUG
//    test_square();
//    printf("original matrix:\n"); for (vertex i = 0; i < n; ++i) { for (vertex j = 0; j < n; ++j) printf("%i ", is_neighbour(i, j)); puts(""); } printf("squared matrix:\n"); for (vertex i = 0; i < n; ++i) { for (vertex j = 0; j < n; ++j) printf("%i ", adj2[i][j]); puts(""); }
#endif
}

#endif
