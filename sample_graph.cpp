#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <set>

#include <unistd.h>

using namespace std;

// Magic numbers for sampling procedure type
#define SAMPLE_RAND 1	// uniform sampling of nodes
#define SAMPLE_FOREST 2	// forest fire

#define MAXLEN 1000 // max length for filename

char in_filename[MAXLEN];

#include "graph.h"
#include "util.h"

int sample_size = 0;
int sample_type = SAMPLE_RAND;
int sampled_edges = 0;

/*
 * Sample vertices uniformly at random and consider the induced subgraph.
 * It may be not connected and have isolated vertices, which we remove.
 *
 */
void sample_rand() {
    vector<int> sampled_vertices = sample_without_replacement(n, sample_size);
    for (int sampled_vertex : sampled_vertices) {
	for (int other_vertex : sampled_vertices) {
	    if (sampled_vertex <= other_vertex) {
		if (find(adj[sampled_vertex].begin(), adj[sampled_vertex].end(), other_vertex) != adj[sampled_vertex].end()) {
		    cerr << vertex_labels[sampled_vertex] << " " << vertex_labels[other_vertex] << endl;
		    sampled_edges++;
		}
	    }
	}
    }

}

/*
 * Sample vertices according to the forest fire approach
 *
 * XXX Assumes a connected graph
 *
 */
void sample_forest_fire() {

    // Start from a random node
    int first = rand() %n; // XXX not perfectly uninformly random
    set<int> sampled = { first };
    if (find(adj[first].begin(), adj[first].end(), first) != adj[first].end()) {
	    cerr << vertex_labels[first] << " " << vertex_labels[first] << endl;
	    sampled_edges++;
    }


    while (sampled.size() < sample_size) {
	set<int> curr_sampled;

	// Expand frontier by adding nodes we haven't touched yet
	for (int sampled_node : sampled) {
	    for (int adj_node : adj[sampled_node]) {
		if (find(sampled.begin(), sampled.end(), adj_node) == sampled.end() && find(curr_sampled.begin(), curr_sampled.end(), adj_node) == curr_sampled.end()) {
		    curr_sampled.insert(adj_node);
		}
	    }
	}

	// If at the last round we found too many new nodes to add to the
	// sample, select some of them randomly so that we have the right
	// sample size.
	if (sampled.size() + curr_sampled.size() > sample_size) {
	    vector<int> new_curr_sampled_vector(curr_sampled.begin(), curr_sampled.end());
	    random_shuffle(new_curr_sampled_vector.begin(), new_curr_sampled_vector.end());
	    new_curr_sampled_vector.resize(sample_size - sampled.size());
	    set<int> new_curr_sampled(new_curr_sampled_vector.begin(), new_curr_sampled_vector.end());
	    curr_sampled = new_curr_sampled;
	}

	// Print the edges connecting the new found nodes to the other sampled nodes
	// (both old and new)
	for (int sampled_vertex : curr_sampled) {
	    for (int other_vertex : sampled) {
		if (find(adj[sampled_vertex].begin(), adj[sampled_vertex].end(), other_vertex) != adj[sampled_vertex].end()) {
		    cerr << vertex_labels[sampled_vertex] << " " << vertex_labels[other_vertex] << endl;
		    sampled_edges++;
		}
	    }
	    for (int other_vertex : curr_sampled) {
		if (sampled_vertex <= other_vertex && find(adj[sampled_vertex].begin(), adj[sampled_vertex].end(), other_vertex) != adj[sampled_vertex].end()) {
		    cerr << vertex_labels[sampled_vertex] << " " << vertex_labels[other_vertex] << endl;
		    sampled_edges++;
		}
	    }
	}
	
	// Add the new sampled edge to the "old" list
	sampled.insert(curr_sampled.begin(), curr_sampled.end());
    }
}

void usage(const char *binary_name) {
    fprintf(stderr, "Usage: %s [-f|-u] [-s sample_size] edges_file\n", binary_name);
    fprintf(stderr, " -f : use forest-fire sampling\n");
    fprintf(stderr, " -u : use uniform node sampling\n");
    fprintf(stderr, " -s sample_size : specify sample size (default = 500)\n");
}

int main(int argc, char **argv) {
    //Parse arguments
    int opt;
    bool specified_f = false;
    bool specified_u = false;
    while ((opt = getopt(argc, argv, "fhs:u")) != -1) {
        switch (opt) {	    
	    case 'f':
		specified_f = true;
		sample_type = SAMPLE_FOREST;
		break;
	    case 'h':
		usage(argv[0]);
		return 0;
		break; // not reached
	    case 's':
		sample_size = strtol(optarg, NULL, 10);
		if (sample_size < 2) {
		    fprintf(stderr, "Sample size must be greater than 1\n");
		    return 1;
		}
		break;
	    case 'u':
		specified_u = true;
		sample_type = SAMPLE_RAND;
		break;
	    default:
		fprintf(stderr, "Wrong option.\n");
		usage(argv[0]);
		return 1;
	}
    }
    if (specified_f && specified_u) {
	fprintf(stderr, "ERROR: Only one between -f and -u can be specified. Exiting.\n");
	return 1;
    }
    if (optind == argc - 1) { 
	strncpy(in_filename, argv[optind], MAXLEN);
    } else {
	fprintf(stderr, "ERROR: wrong number of arguments passed\n");
	return 1;
    }

#ifndef DEBUG
    srand(time(NULL));
#endif

    // Read input graph
    read_graph();

    // Call the right sample procedure
    if (sample_type == SAMPLE_RAND) {
	sample_rand();
    } else if (sample_type == SAMPLE_FOREST) {
	sample_forest_fire();
    }
    printf("Sampled: %d nodes, %d edges\n", sample_size, sampled_edges);

    return 0;
}

