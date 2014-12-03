#ifndef GRASS_H
#define GRASS_H

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <map>
#include <set>
#include <vector>

#include "graph.h"
#include "summ.h"
#include "util.h"

int error_type = 0;
float sample_size_multiplier = 1.0;
double curr_reconstruct_err = 0.0; // reconstruction error
//double (*reconstruction_error)(map<pair<int,int>,double>); // function pointer
#ifdef DEBUG
double (*reconstruction_error)(); // function pointer
#endif
double (*cut_norm_error_function)(const map<pair<int,int>,double> &);

/*
 * Compute the initial densities between supernodes
 *
 * Basically, this corresponds to the adjacency matrix.
 *
 */
void build_initial_densities() {
#ifdef DEBUG
	printf("Computing initial densities...");
#endif
	// Basically build the explicit identity matrix
	for (int i = 0; i < n; i++) {
		double dens = is_neighbour(i,i) ? 1.0 : 0.0;
		pair<int,int> coordinates = make_pair(i, i);
		densities.insert(make_pair(coordinates, dens));
		for (int j = i+1; j < n; j++) {
			dens = is_neighbour(i,j) ? 1.0 : 0.0;
			coordinates = make_pair(i, j);
			densities.insert(make_pair(coordinates, dens));
			coordinates = make_pair(j, i);
			densities.insert(make_pair(coordinates, dens));
		}
	}

#ifdef DEBUG
	printf("done\n");
#endif
}

/*
 * Put each node in its own supernode and compute initial densities
 *
 */
void build_initial_supernodes() {
#ifdef DEBUG
	printf("Computing initial supernodes...");
#endif
	labels.reserve(n);
	// At the beginning, each node is in its own supernode
	for (int i = 0; i < n; ++i) {
		supernode curr;
		curr.insert(i);
		supernodes.insert(make_pair(i, curr));
		labels.push_back(i);
	}
	labels.shrink_to_fit();
#ifdef DEBUG
	assert(verify_labels());
	printf("done\n");
#endif
	// Compute initial densities
	build_initial_densities();
}

/*
 * Compute the density matrix ret representing the densities after the supernodes with
 * labels first and second have been merged.
 *
 * Return the change in reconstruction error.
 *
 */
double update_densities_after_merging(int first, int second, map<pair<int,int>,double>& ret) {
#ifdef DEBUG
	printf("Merging %d and %d\n", first, second);
	printf("Updating densities...");
#endif

	int first_size = supernodes[first].size();
	int second_size = supernodes[second].size();
	int merged_size = first_size + second_size;
	// The densities for pair of supernodes different than the two that we
	// are trying to merge do not change, so we can carry them forward.
	for (auto it : densities) {
		pair<int,int> key = it.first;
		if (key.first != first && key.first != second && key.second != first && key.second != second) {
			ret.insert(make_pair(key, it.second));
		}
	}

	// Variables to compute the change in reconstruction error
	double old_error = 0.0;
	double new_error = 0.0;

	// Insert internal density of the new supernode.
	pair<int,int> self_key = make_pair(first,first);
	pair<int,int> second_key = make_pair(second,second);
	pair<int,int> between_key = make_pair(first,second);

	double internal_edges = (densities[self_key] * pow(first_size, 2.0) + densities[second_key] * pow(second_size, 2.0)) / 2 + densities[between_key] * first_size * second_size;
	double internal_density = 2 * internal_edges / pow(merged_size, 2.0);
	ret.insert(make_pair(self_key, internal_density));

	double multip = 1;
	if (error_type == ERROR_L1 || error_type == ERROR_L2) {
		if (error_type == ERROR_L1) 
			multip = 2;
		old_error += multip * 2 * (((double) first_size) / n) * (((double) second_size) / n) * densities[between_key] * (1 - densities[between_key]);
		old_error += multip * pow(((double) first_size) / n, 2.0) * densities[self_key] * (1 - densities[self_key]);
		old_error += multip * pow(((double) second_size) / n, 2.0) * densities[second_key] * (1 - densities[second_key]);
		new_error =  multip * (pow(((double) merged_size) / n, 2.0)) * internal_density * (1 - internal_density);
	}

	// Recompute densities between the new supernode and all other supernodes
	for (auto it : supernodes) {
		if (it.first == first || it.first == second) {
			continue;
		}
		pair<int,int> cross_key = make_pair(first, it.first);
		pair<int,int> second_cross_key = make_pair(second, it.first);
		double cross_edges = densities[cross_key] * (first_size * it.second.size()) + densities[second_cross_key] * (second_size * it.second.size());
		double cross_density = cross_edges / (merged_size * it.second.size());
		ret.insert(make_pair(cross_key, cross_density));
		pair<int,int> inverted_cross_key = make_pair(it.first, first);
		ret.insert(make_pair(inverted_cross_key, cross_density));

		if (error_type == ERROR_L1 || error_type == ERROR_L2) {
			old_error += multip * 2 * (((double) first_size) / n) * (((double) it.second.size()) / n) * densities[cross_key] * (1.0 - densities[cross_key]);
			old_error += multip * 2 * (((double) second_size) / n) * (((double) it.second.size()) / n) * densities[second_cross_key] * (1.0 - densities[second_cross_key]);
			new_error += multip * 2 * (((double) merged_size) / n) * (((double) it.second.size()) / n) * cross_density * (1.0 - cross_density);
		}
	}

#ifdef DEBUG
	printf("done\n");
#endif
	if (error_type == ERROR_LC) {
		old_error = curr_reconstruct_err;
		new_error = cut_norm_error_function(ret);
	}

	return new_error - old_error;
}

/*
 * Sample a number of pairs of supernodes, to be used to select which pair to
 * merge
 *
 */
vector<array<int, 2>> sample_candidates(int sample_size) {
#ifdef DEBUG
	printf("Sampling candidates...");
#endif
#ifdef DEBUG
	printf("sample_size=%i...", sample_size);
#endif
	vector<int> supernode_keys;
	supernode_keys.reserve(supernodes.size());
	for(auto it : supernodes) {
		  supernode_keys.push_back(it.first);
	}

	set<array<int, 2>> sampled_pairs;
	int sampled = 0;
	while (sampled < sample_size) {
		int first = rand() % supernodes.size();
		int second = 0;
		do {
			second = rand() %supernodes.size();
		} while (second == first);
		array<int,2> pair = {supernode_keys[first], supernode_keys[second]};
		// Ensure order
		if (pair[0] > pair[1]) {
			int tmp = pair[1];
			pair[1] = pair[0];
			pair[0] = tmp;
		}
		// If we haven't sampled this pair yet, add it.
		if (sampled_pairs.find(pair) == sampled_pairs.end()) {
#ifdef DEBUG
			printf("(%d, %d), ", pair[0], pair[1]);
#endif
			sampled_pairs.insert(pair);
			sampled++;
		}
	}
	vector<array<int,2>> candidates(sampled_pairs.begin(), sampled_pairs.end());
#ifdef DEBUG
	printf("done\n");

#endif
	return candidates;
}

/*
 * Sample a number of pairs of supernodes, to be used to select which pair to
 * merge
 *
 * OLD CODE
 *
 */
vector<array<int, 2>> old_sample_candidates(int sample_size) {
#ifdef DEBUG
	printf("Sampling candidates...");
#endif
	vector<array<int,2>> candidates(sample_size * (sample_size -1) / 2);
#ifdef DEBUG
	printf("sample_size=%i...", sample_size);
#endif
	for (auto iter1 = supernodes.begin(); iter1 != supernodes.end(); ++iter1) {
		for (auto iter2 = supernodes.begin(); iter2 != supernodes.end(); ++iter2) {
			if (iter1->first < iter2->first) {
				array<int,2> curr;
				curr[0] = iter1->first;
				curr[1] = iter2->first;
				candidates.push_back(curr);
			}
		}
	}
	if (sample_size != candidates.size()) {
		random_shuffle(candidates.begin(), candidates.end());
		candidates.resize(sample_size);
	}
	candidates.shrink_to_fit();
#ifdef DEBUG
	printf("done\n");
#endif
	return candidates;
}

/*
 * Find and merge the "optimal" pair of supernodes
 *
 * The pair is chosen using sampling: we sample a number of pairs and merge the
 * best one.
 */
void merge_best_supernodes(int sample_size) {
#ifdef DEBUG
	printf("Merging supernodes...\n");
#endif
	// Sample candidates
	vector<array<int,2>> candidates = sample_candidates(sample_size);

	// Various variables to store the best candidate
	array<int,2> best_candidate;
	double best_new_reconstruct_err_change = 1.0; // This works because we consider the normalized error
	map<pair<int,int>,double> best_densities;

	// Evaluate candidates
	printf("%d candidates: ", sample_size);
	int iteration = 0;
	for (array<int,2> candidate : candidates) {
		if (iteration % (sample_size / 10 + 1) == 0) {
			printf("%d..", iteration);
		}
		iteration++;

		// Update the densities and compute the reconstruction error
		map<pair<int,int>,double> candidate_densities;
		double new_reconstruct_err_change =
			update_densities_after_merging(candidate[0], candidate[1], candidate_densities);
		if (new_reconstruct_err_change < best_new_reconstruct_err_change) {
			// Is this the candidate that increases the
			// reconstruction error the least?
			// If yes, then save it.
#ifdef DEBUG
			printf("\t%f new best. Old was %f\n", new_reconstruct_err_change, best_new_reconstruct_err_change);
#endif
			best_candidate = candidate;
			best_new_reconstruct_err_change = new_reconstruct_err_change;
			best_densities = candidate_densities;
		}
	}
	printf("\nMerged %d and %d\n", best_candidate[0], best_candidate[1]);

	curr_reconstruct_err += best_new_reconstruct_err_change;
	for (int node : supernodes[best_candidate[1]]) {
		labels[node] = best_candidate[0];
		supernodes[best_candidate[0]].insert(node);
	}
	supernodes.erase(best_candidate[1]);
	densities = best_densities;

#ifdef DEBUG
	assert(verify_densities(densities));
	assert(verify_labels());
	printf("done (merging supernodes)\n");
#endif
}

/*
 * Build a summary of the graph with k supernodes
 *
 * This is basically an agglomerative hierarchical clustering approach: at the
 * beginning each node is in its own supernode, then, until there are more than
 * k supernodes, two supernodes are merged in a single one. The pair of
 * supernodes to be merged is the one that makes the reconstruction error
 * increase the least.
 *
 * Since finding the actual optimal pair would take too much time, we sample a
 * number of pair of supernodes and merge the best among these.
 *
 * The parameter c controls the sample size.
 *
 */
void grass_sample_pairs(int k, float c, int my_error_type, bool approx_reconstruction_error __attribute__ ((unused))) {
	sample_size_multiplier = c;
	num_supernodes = k;
	error_type = my_error_type;

#ifdef DEBUG
	// Select the desired error function
	if (error_type == ERROR_L1) {
		reconstruction_error = &fast_l1_reconstruction_error;
	} else if (error_type == ERROR_L2) {
		reconstruction_error = &fast_l2_reconstruction_error;
	}
#endif
	if (error_type == ERROR_LC) {
		cut_norm_error_function = approx_reconstruction_error ? &cut_norm_error_approx : &cut_norm_error;
	}

	// Each node in its own supernode 
	build_initial_supernodes();

        int sample_size = 0;
	while (supernodes.size() > num_supernodes) { 
		printf("Currently: %zu supernodes\n", supernodes.size());
		//sample_size = sample_size_multiplier * (( supernodes.size() * (supernodes.size() -1)) / 2);
		sample_size = sample_size_multiplier * supernodes.size();
                if (sample_size == 0) {
                    sample_size = 1;
                }
		merge_best_supernodes(sample_size);
	}
}

/*
 * Create a summary by randomly merging supernodes
 *
 */
void random_summary(int k) {
    num_supernodes = k;

    // Each node in its own supernode 
    build_initial_supernodes();

    while (supernodes.size() > num_supernodes) { 
#ifdef DEBUG
        printf("%zu supernodes\n", supernodes.size());
#endif
	array<int,2> candidate;
	candidate[0] = candidate[1] = rand() % supernodes.size();
	do {
		candidate[1] = rand() % supernodes.size();
	} while (candidate[1] == candidate[0]);

	// Update the densities 
	map<pair<int,int>,double> new_densities;
	update_densities_after_merging(candidate[0], candidate[1], new_densities);
	densities = new_densities;
	// Merge the nodes
	for (int node : supernodes[candidate[1]]) {
		supernodes[candidate[0]].insert(node);
		// Change the labels 
		labels[node] = candidate[0];
	}
	supernodes.erase(candidate[1]);
#ifdef DEBUG
        assert(verify_labels());
#endif
    }
}

#endif

