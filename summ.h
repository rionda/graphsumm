#ifndef SUMM_H
#define SUMM_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "cutnorm.h"

#define CUT_NORM_SAMPLE_SIZE 500    // The sample size to compute the approx to the cut-norm. It's actually the square of this.
#define REC_ERR_SAMPLE_SIZE 10000   // The sample size to copmute the approx to the reconstruct error.

#define ERROR_L1 1                  // Magic numbers for error types. 1 = L1, 2 = L2, 3 = cut-norm
#define ERROR_L2 2
#define ERROR_LC 3

#define RES_MAXLINE 100000

typedef unordered_set<int> supernode;       // a set of vertex ids.

unordered_map<int, supernode> supernodes;   // the key is like a "supernode id".
map<pair<int, int>, double> densities;      // the key is a pair of supernode ids.
vector<int> labels;                         // for each node, store the id of the supernode it belongs to.

int num_supernodes = 0;                     // number of supernodes in the end
double expected_triangles = -1;
char summary_filename[maxlen];

double original_size() {
    return entropy(avg_dens) * n * n / 2 + 1 + log2(num_edges);
}

int bits(int x) {
    return x == 0 ? 1 : 1 + floor(log2(x));
}

double encoding_length(int x) {
    return 1 + bits(x) + bits(bits(x));
}

double best_encoding(double factor, double dens) {
    double edges = factor * dens;
    return factor * entropy(dens) + encoding_length(min(edges, 1 - edges)) + 1;

    if (dens == 0 || dens == 1) return 2;

    double best = infinity;
    for (double p = 1; dens * p <= edges; p *= 2) {
        double approx_dens = round(dens * p) / p;
        double used = cross_entropy(dens, approx_dens) * factor + encoding_length(approx_dens * p);
        if (used < best)
            best = used;
    }
    return best;
}

/*
 * Compute the number of bits needed to encode the summary
 *
 */
double compressed_size() {
    double bits = 0.0;
    // Same as the code above that's been commented out (for 0-1 graphs), but faster.
    for (auto first_it : supernodes) {
        for (auto second_it : supernodes) if (second_it.first <= first_it.first) {
            pair<int, int> key = make_pair(first_it.first, second_it.first);
            double dens = densities[key];
            double factor = (double) first_it.second.size() * second_it.second.size();
            if (second_it.first == first_it.first) factor /= 2;

            bits += best_encoding(factor, dens);
        }
    }

    // Bits to encode the partition. It's a function from [n] to [num_supernodes].
    bits += n * log2(num_supernodes);
    return bits;
}

/*
 * Compute the index of the summary
 *
 * The closer the index is to the density of the graph, the closer the partition
 * of the summary is to a (non-weak) regular partition.
 * The index is the expectation, over all ordered pair of vertices, of the
 * squared densities of the supernodes containing u and v.
 * In other words, it's sum_{u,v} (|U||V|/n^2) d(cl(u), cl(v))^2.
 * In the most regular partition (of a simple unweighted graph), all densities
 * are 1 or 0 (think of putting every vertex into a supernode), so the square
 * doesn't matter and you get the average density of the graph.
 * The index is used in proofs of the reg. lemma by arguing that if you take
 * two vertex sets which violate epsilon-regularity and subdivide them, the
 * index increases by poly(eps), and since the index is bounded by 1 the process
 * goes on for poly(1/\eps) steps.
 *
 */
double index() {
    double ind = 0.0;
    for (auto f_it = supernodes.begin(); f_it != supernodes.end(); ++f_it) {
        for (auto s_it = supernodes.begin(); s_it != supernodes.end(); ++s_it) {
            double dens = densities[make_pair(f_it->first, s_it->first)];
            ind += ((double) f_it->second.size() / n) * ((double) s_it->second.size() / n) * pow(dens, 2.0);
        }
    }
    return ind;
}

// works for the expected adjacency matrix
double fast_l1_reconstruction_error() {
    double ind = 0.0;
    for (auto f_it = supernodes.begin(); f_it != supernodes.end(); ++f_it) {
        for (auto s_it = supernodes.begin(); s_it != supernodes.end(); ++s_it) {
            double dens = densities[make_pair(f_it->first, s_it->first)];
            ind += 2 * ((double) f_it->second.size() / n) * ((double) s_it->second.size() / n) * dens * (1 - dens);
        }
    }
    return ind;
}

double fast_l1_reconstruction_error_with_majority() {
    double ind = 0.0;
    for (auto f_it = supernodes.begin(); f_it != supernodes.end(); ++f_it) {
        for (auto s_it = supernodes.begin(); s_it != supernodes.end(); ++s_it) {
            double dens = densities[make_pair(f_it->first, s_it->first)];
            ind += ((double) f_it->second.size() / n) * ((double) s_it->second.size() / n) * min(dens, 1 - dens);
        }
    }
    return ind;
}

// works for the expected adjacency matrix
double fast_l2_reconstruction_error() {
    double ind = 0.0;
    for (auto f_it = supernodes.begin(); f_it != supernodes.end(); ++f_it) {
        for (auto s_it = supernodes.begin(); s_it != supernodes.end(); ++s_it) {
            double dens = densities[make_pair(f_it->first, s_it->first)];
            ind += ((double) f_it->second.size() / n) * ((double) s_it->second.size() / n) * dens * (1 - dens);
        }
    }
    return ind;
}


/*
 * Compute the (inherently approximate) (normalized) cut-norm error of the
 * density matrix dens
 *
 */
double cut_norm_error(const map<pair<int,int>,double>& dens) {
#ifdef DEBUG
    printf("  computing \"exact\" cut-norm reconstruction error...");
#endif
    double ret = 0.0;

    // Build the difference matrix
    Matrix<double> m(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double a = is_neighbour(i, j) ? 1 : 0;
            auto density_it = dens.find(make_pair(labels[i],labels[j]));
#ifdef DEBUG
            assert(density_it != dens.end());
#endif
            m[i][j] = a - density_it->second;
        }
    }
    ret = cut_norm(m) / pow(n, 2.0); // compute and normalize;
#ifdef DEBUG
    printf("done\n");
#endif
    return ret;
}

/*
 * Compute the (inherently approximate) approximated (through sampling)
 * (normalized) cut-norm error of the density matrix dens
 *
 */
double cut_norm_error_approx(const map<pair<int,int>,double>& dens) {
#ifdef DEBUG
    printf("  computing approx cut-norm reconstruction error...");
#endif
    vector<int> sample = sample_without_replacement(n, CUT_NORM_SAMPLE_SIZE);
    int s = sample.size();
    double ret = 0.0;

    // Build the difference matrix with the sampled entries.
    Matrix<double> m(s, vector<double>(s));
    for (int i = 0; i < s; ++i) {
        for (int j = 0; j < s; ++j) {
            int v = sample[i], w = sample[j];
            double a = is_neighbour(v, w) ? 1 : 0;
            auto density_it = dens.find(make_pair(labels[v],labels[w]));
#ifdef DEBUG
            assert(density_it != dens.end());
#endif
            m[i][j] = a - density_it->second;
        }
    }
    //ret = fabs(cut_norm(m)) * pow((double) n / s, 2);
    ret = cut_norm(m) / pow(s, 2.0); // compute and normalize;
#ifdef DEBUG
    printf("done\n");
#endif
    return ret;
}

/*
 * Verify that the dens matrix is representing of the densities of the
 * supernodes
 *
 */
bool verify_densities(map<pair<int,int>,double> dens) {
    bool ret = true;
    for (auto dens_it = dens.begin(); dens_it != dens.end(); ++dens_it) {
        double density = dens_it->second;
        int label_first = dens_it->first.first;
        int label_second = dens_it->first.second;
        int density_denom = (long long)supernodes[label_first].size() * supernodes[label_second].size();

        /* The following is commented out because it was for the old
         * definition of densities from [LeFevreT10]
         */
        //if (label_first == label_second) {
        //  if (supernodes[label_first].size() == 1) {
        //      density_denom = 1;
        //  } else {
        //      density_denom = supernodes[label_first].size() * (supernodes[label_first].size() - 1);
        //  }

        //} else {
        //  density_denom = supernodes[label_first].size() * supernodes[label_second].size();
        //}
        int edges = 0;
        int nodes_with_label_first = 0;
        int nodes_with_label_second = 0;
        for (int i = 0; i < n; i++) {
            if (labels[i] == label_first) {
                nodes_with_label_second = 0;
                nodes_with_label_first++;
                for (int j = 0; j < n; j++) {
                    if (labels[j] == label_second) {
                        nodes_with_label_second++;
                        if (is_neighbour(i,j)) {
                            edges++;
                        }
                    }
                }
            }
        }

        double new_density = ((double) edges) / density_denom;
        if (new_density > 1.0) {
            printf("SOMETHING IS WRONG: new_density > 1.0\n");
            printf("%d %d %f %f %d %d %d\n", label_first, label_second, density, new_density, nodes_with_label_first, nodes_with_label_second, edges);
            ret = false;
            break;
        }
        if (fabs(new_density - density) > 1e-8) {
            printf("SOMETHING IS WRONG: new_density:%lf density:%lf diff=%lf\n", new_density, density, fabs(new_density -density));
            ret = false;
            break;
        }
        printf("label_first=%d nodes_with_label_first=%d supernodes[label_first].size()=%zu\n", label_first, nodes_with_label_first, supernodes[label_first].size());
        printf("label_second=%d nodes_with_label_second=%d supernodes[label_second].size()=%zu\n", label_second, nodes_with_label_second, supernodes[label_second].size());
        assert(nodes_with_label_first > 0);
        assert(nodes_with_label_second > 0);
        assert(nodes_with_label_first == supernodes[label_first].size());
        assert(nodes_with_label_second == supernodes[label_second].size());
    }
    return ret;
}

/*
 * Check that each node has label corresponding to the supernode it belongs to
 *
 */
bool verify_labels() {
    for (auto it = supernodes.begin(); it != supernodes.end(); ++it) {
        int curr_label = it->first;
        for (int node : it->second) {
            if (labels[node] != curr_label) {
                printf("MISMATCHED_LABEL: node %d label %d supernode %d\n", node, labels[node], curr_label);
                return false;
            }
        }
    }
    return true;
}

/*
 * Compute and return statistics about the degree error from the summary
 *
 */
array<double,7> expected_degree_error() {
    map<int, double> supernodes_exp_degrees;
    for (auto i : supernodes) {
        int label_i = i.first;
        supernodes_exp_degrees[label_i] = 0.0;
        for (auto j : supernodes) {
            int label_j = j.first;
	    if (label_j != label_i) {
	        supernodes_exp_degrees[label_i] += j.second.size() * densities[make_pair(label_i,label_j)];
	    } else {
		if (j.second.size() > 1) {
		    supernodes_exp_degrees[label_i] += (pow(j.second.size(), 2.0) * densities[make_pair(label_i,label_j)]) /  (j.second.size() -1);
		    }
	    }
        }
    }

    double err_max = -1000000;
    double err_min = 1000000;
    double err_sum = 0;
    double squared_err_sum = 0.0;
    double abs_err_max = -1.0;
    double abs_err_min = 1000000;
    double abs_err_sum = 0;
    for (int i = 0; i < n; i++) {
        double err = adj[i].size() - supernodes_exp_degrees[labels[i]];
        double abs_err = fabs(err);
        err_sum += err;
        abs_err_sum += abs_err;
        squared_err_sum += pow(err, 2);
        if (err > err_max) {
           err_max = err;
        } 
        
        if (err < err_min) {
            err_min = err;
        }
        if (abs_err > abs_err_max) {
            abs_err_max = abs_err;
        } 
        if (abs_err < abs_err_min) {
            abs_err_min = abs_err;
        }
    }

    array<double, 7> res = { err_sum / n, sqrt(squared_err_sum / n - pow(err_sum / n, 2)), err_max, err_min, abs_err_sum / n, abs_err_max, abs_err_min };
    return res;
}

/*
 * Compute and return statistics about the (relative) degree error from the summary
 *
 */
array<double,4> expected_rel_degree_error() {
    map<int, double> supernodes_exp_degrees;
    for (auto i : supernodes) {
        int label_i = i.first;
        supernodes_exp_degrees[label_i] = 0.0;
        for (auto j : supernodes) {
            int label_j = j.first;
	    if (label_j != label_i) {
	        supernodes_exp_degrees[label_i] += j.second.size() * densities[make_pair(label_i,label_j)];
	    } else {
	        if (j.second.size() > 1) {
		    supernodes_exp_degrees[label_i] += (pow(j.second.size(), 2.0) * densities[make_pair(label_i,label_j)]) /  (j.second.size() -1);
		}
	    }
	}
    }

    double rel_err_max = -1.0;
    double rel_err_min = n;
    double rel_err_sum = 0;
    double squared_rel_err_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double err = fabs(adj[i].size() - supernodes_exp_degrees[labels[i]]);
        double rel_err = err / adj[i].size();
        rel_err_sum += rel_err;
        squared_rel_err_sum += pow(rel_err, 2);
        if (rel_err > rel_err_max) {
            rel_err_max = rel_err;
        } 
        if (rel_err < rel_err_min) {
            rel_err_min = rel_err;
        }
    }

    array<double, 4> res = { rel_err_sum / n, sqrt(squared_rel_err_sum / n - pow(rel_err_sum / n, 2)), rel_err_max, rel_err_min };
    return res;
}

/*
 * Compute and return statistics about eigenvector centrality
 *
 */
array<double, 7> expected_eigenvector_centrality_error() {
    array<double, 7> degree_res = expected_degree_error();
    double factor = 2 * num_edges;
    array<double, 7> res;
    for (int i = 0; i < 8; i++) {
        res[i] = degree_res[i] / factor;
    }
    return res;
}

/*
 * Compute the expected number of triangles
 *
 */
double expected_triangles_number() {
    if (expected_triangles == -1.0) {
        printf("Computing expected triangles number...");
        map<pair<int,int>,double> new_densities(densities);
        for (auto it : supernodes) {
            int label_i = it.first;
            double size_i = it.second.size();
	    if (size_i > 1) {
		    new_densities[make_pair(label_i,label_i)] = densities[make_pair(label_i,label_i)] * size_i / (size_i - 1);
	    }
        }
        expected_triangles = 0.0;
        for (auto it_i = supernodes.begin(); it_i != supernodes.end(); it_i++) {
            int label_i = it_i->first;
            double size_i = it_i->second.size();
            pair<int,int> ii_pair = make_pair(label_i,label_i);
            // Handle the case when a supernode has less than 3 nodes (no
            // internal triangles possible)
            double triplets = (size_i >= 3) ? size_i * (size_i - 1) * (size_i - 2) / 6.0 : 0;
            expected_triangles += triplets * pow(new_densities[make_pair(label_i,label_i)], 3);
            auto it_j = it_i;
            it_j++;
            while (it_j != supernodes.end()) {
                int label_j = it_j->first;
                double size_j = it_j->second.size();
                double third_addend = 0.0;
                pair<int,int> jj_pair = make_pair(label_j,label_j);
                pair<int,int> ij_pair = make_pair(label_i,label_j);
                auto it_w = it_j;
                it_w++;
                while (it_w != supernodes.end()) {
                    int label_w = it_w->first;
                    int size_w = it_w->second.size();
                    third_addend += size_i * size_j * size_w * new_densities[ij_pair] * new_densities[make_pair(label_i,label_w)] * new_densities[make_pair(label_j,label_w)];
                    it_w++;
                }
                // Handle the case when a supernode has less than 2 nodes (no
                // triangles with two nodes in it possible)
                double pairs_i = (size_i >= 2) ? size_i * ((double) (size_i -1)) / 2 : 0.0;
                double pairs_j = (size_j >= 2) ? size_j * ((double) (size_j -1)) / 2 : 0.0;
                double first_addend = size_j * new_densities[ii_pair] * pairs_i;
                double second_addend = size_i * new_densities[jj_pair] * pairs_j;

                expected_triangles += pow(new_densities[ij_pair], 2) * (first_addend + second_addend + third_addend);
                it_j++;
            }
        }
        printf("done: %lf\n", expected_triangles);
    } else {
        printf("Cached expected triangles: %lf\n", expected_triangles);
    }
    return expected_triangles;
}

/*
 * Return the error for the number of triangles
 *
 */
double expected_triangles_error() {
    if (triangles == -1) {
        triangles_number();
    } else {
        printf("Triangles (cached): %lld\n", triangles);
    }
    if (expected_triangles == -1) {
        expected_triangles_number();
    } else {
        printf("Expected Triangles (cached): %lf\n", expected_triangles);
    }
    return expected_triangles - triangles;
}

/*
 * Return the error for the clustering coeffient
 *
 */
double expected_clustering_coefficient_error() {
    return 6 * expected_triangles_error() / ((ll)n * (n-1) * (n-2));
}

/*
 * Return the relative error for the clustering coeffient
 *
 */
double expected_rel_clustering_coefficient_error() {
    return expected_triangles_error() / triangles;
}

/*
 * Compute a bunch of statistics about the sizes of the
 * supernodes and the densities
 *
 */
vector<double> get_summary_statistics() {
    printf("Computing summary statistics ...");
    vector<double> stats;
    double size_max = -1;
    double size_min = n;
    double size_sum = 0;
    double squared_size_sum = 0;

    double int_dens_max = -1;
    double int_dens_min = 2;
    double int_dens_sum = 0;
    double squared_int_dens_sum = 0;

    double ext_dens_max = -1;
    double ext_dens_min = 2;
    double ext_dens_sum = 0;
    double squared_ext_dens_sum = 0;

    for (auto it_i = supernodes.begin(); it_i != supernodes.end(); it_i++) {
        int label_i = it_i->first;
        int size_i = it_i->second.size();
        size_sum += size_i;
        squared_size_sum += pow(size_i, 2);
        if (size_i > size_max) {
            size_max = size_i;
        } 
        if (size_i < size_min) {
            size_min = size_i;
        }
        double dens = densities[make_pair(label_i,label_i)];

        int_dens_sum += dens;
        squared_int_dens_sum += pow(dens, 2);
        if (dens > int_dens_max) {
            int_dens_max = dens;
        }
        if (dens < int_dens_min) {
            int_dens_min = dens;
        }
        auto it_j = it_i;
        it_j++;
        while (it_j != supernodes.end()) {
            int label_j = it_j->first;
            dens = densities[make_pair(label_i,label_j)];
            ext_dens_sum += dens;
            squared_ext_dens_sum += pow(dens, 2);
            if (dens > ext_dens_max) {
                ext_dens_max = dens;
            }
            if (dens < ext_dens_min) {
                ext_dens_min = dens;
            }
            it_j++;
        }
    }
    double dens_max = (ext_dens_max > int_dens_max) ? ext_dens_max : int_dens_max;
    double dens_min = (ext_dens_min < int_dens_min) ? ext_dens_min : int_dens_min;
    double dens_sum = ext_dens_sum + int_dens_sum;
    double squared_dens_sum = squared_ext_dens_sum + squared_int_dens_sum;

    double size_avg = size_sum / num_supernodes;
    double size_stddev = sqrt((squared_size_sum / num_supernodes) - pow(size_avg, 2));

    double int_dens_avg = int_dens_sum / num_supernodes;
    double int_dens_stddev = sqrt((squared_int_dens_sum / num_supernodes) - pow(int_dens_avg, 2));

    double ext_dens_avg = ext_dens_sum / (num_supernodes * (num_supernodes -1));
    double ext_dens_stddev = sqrt((squared_ext_dens_sum / (num_supernodes * (num_supernodes -1))) - pow(ext_dens_avg, 2));

    double dens_avg = dens_sum / pow(num_supernodes,2);
    double dens_stddev = sqrt((squared_dens_sum / pow(num_supernodes,2)) - pow(dens_avg, 2));

    stats.push_back(size_avg);
    stats.push_back(size_max);
    stats.push_back(size_min);
    stats.push_back(size_stddev);
    stats.push_back(dens_avg);
    stats.push_back(dens_max);
    stats.push_back(dens_min);
    stats.push_back(dens_stddev);
    stats.push_back(int_dens_avg);
    stats.push_back(int_dens_max);
    stats.push_back(int_dens_min);
    stats.push_back(int_dens_stddev);
    stats.push_back(ext_dens_avg);
    stats.push_back(ext_dens_max);
    stats.push_back(ext_dens_min);
    stats.push_back(ext_dens_stddev);

    printf("done\n");
    return stats;
}

/*
 * Compute adjacency error
 * 
 * This is the l1 reconstruction error, but we want statistics about it.
 * 
 * This is the old function that took O(n^2);
 *
 */
array<double,4> old_expected_adjacency_error() {
    map<pair<int,int>,double> new_densities(densities);
    for (auto it : supernodes) {
        int label_i = it.first;
        int size_i = it.second.size();
	if (size_i > 1) {
            new_densities[make_pair(label_i,label_i)] = new_densities[make_pair(label_i,label_i)] * size_i / (size_i - 1);
	}
    }
    array<double,4> ret;
    double sum = 0.0;
    double squared_sum = 0.0;
    double max = -1;
    double min = n;
    for (int i = 0; i < n; ++i) {
        int label_i = labels[i]; 
        for (int j = i; j < n; ++j) {
            int label_j = labels[j] ; 
            double edge = is_neighbour(i, j) ? 1.0 : 0.0;
	    if (new_densities[make_pair(label_i, label_j)] > 0) {
		double err = fabs(edge - new_densities[make_pair(label_i, label_j)]);
		sum += err;
		squared_sum += pow(err, 2.0);
		if (err > max) {
		    max = err;
		}
		if (err < min) {
		    min = err;
		}
            }
        }
    }
    ret[0] = sum / pow(n,2);
    ret[1] = sqrt((squared_sum / pow(n,2)) - pow(ret[0], 2));
    ret[2] = max;
    ret[3] = min;
    return ret;
}

/*
 * Compute adjacency error
 * 
 * This is the l1 reconstruction error, but we want statistics about it.
 * 
 * This is the new function that takes O(k^2);
 *
 */
array<double,4> expected_adjacency_error() {
    map<pair<int,int>,double> new_densities(densities);
    for (auto it : supernodes) {
        int label_i = it.first;
        int size_i = it.second.size();
	if (size_i > 1) {
            new_densities[make_pair(label_i,label_i)] = new_densities[make_pair(label_i,label_i)] * size_i / (size_i - 1);
	}
    }
    array<double,4> ret;
    double sum = 0.0;
    double squared_sum = 0.0;
    double max = -1;
    double min = n;
    for (auto it_i = supernodes.begin(); it_i != supernodes.end(); ++it_i) {
        int label_i = it_i->first;
        double err_if_not_edge = new_densities[make_pair(label_i, label_i)];
        double err_if_edge = 1.0 - err_if_not_edge;
        double max_err = ( err_if_edge > 0.5) ? err_if_edge : err_if_not_edge;
        double min_err = ( err_if_edge < 0.5) ? err_if_edge : err_if_not_edge;
	int total_possible_edges = supernodes[label_i].size() * (supernodes[label_i].size() -1) / 2;
	//int total_possible_edges = ((supernodes[label_i].size() *
	//			(supernodes[label_i].size()-1)) / 2) +
	//	supernodes[label_i].size();
	// There is no way to compute the number of internal edges only from the
	// internal density and the supernode size.
	// (commented out as we ignore self loops when building the graph)
	/*
	 * int self_loops = 0;
	for (int node_id : supernodes[label_i]) {
		if (binary_search(adj[node_id].begin(), adj[node_id].end(), node_id)) {
			self_loops++;
		}
	}
	int present_edges = self_loops + ((err_if_not_edge * pow(supernodes[label_i].size(),2.0) - self_loops) / 2);
	*/
	int present_edges = total_possible_edges * new_densities[make_pair(label_i,label_i)];
	if (present_edges > 0) {
   	    if (max_err > max) {
		max = max_err;
            }
    	    if (min_err < min) {
		min = min_err;
            }
	}
	sum += err_if_not_edge * (total_possible_edges - present_edges) +
		err_if_edge * present_edges;
	squared_sum += pow(err_if_not_edge, 2.0) * (total_possible_edges -
			present_edges) + pow(err_if_edge, 2.0) * present_edges;

        auto it_j = it_i;
        it_j++;
        while (it_j != supernodes.end()) {
	    int label_j = it_j->first;
	    err_if_not_edge = densities[make_pair(label_i, label_j)];
	    err_if_edge = 1.0 - err_if_not_edge;
	    max_err = ( err_if_edge > 0.5) ? err_if_edge : err_if_not_edge;
            min_err = ( err_if_edge < 0.5) ? err_if_edge : err_if_not_edge;
	    total_possible_edges = supernodes[label_i].size() * supernodes[label_j].size();
	    present_edges = err_if_not_edge * total_possible_edges;
	    if (present_edges > 0) {
   	        if (max_err > max) {
	    	max = max_err;
                }
    	        if (min_err < min) {
	    	min = min_err;
                }
	    }
	    sum += err_if_not_edge * (total_possible_edges - present_edges) +
		    err_if_edge * present_edges;
	    squared_sum += pow(err_if_not_edge, 2.0) * (total_possible_edges -
			    present_edges) + pow(err_if_edge, 2.0) *
		    present_edges;
	    it_j++;
	}
    }
    ret[0] = sum / pow(n,2);
    ret[1] = sqrt((squared_sum / pow(n,2)) - pow(ret[0], 2));
    ret[2] = max;
    ret[3] = min;
    return ret;
}

/*
 * Run various queries on the summary and compute errors
 *
 */
array<double,16> get_query_errors() {
    printf("Getting query errors...\n");
    printf("\tAdjacency error...");
    array<double,4> adj_err = expected_adjacency_error();
    printf("\tdone\n");
    printf("\tRelative degree error...");
    array<double,4> rel_deg_err = expected_rel_degree_error();
    printf("done\n");
    printf("\tDegree error...");
    array<double,7> deg_err = expected_degree_error();
    printf("done\n");
    printf("\tRelative clustering coefficient error...");
    double clust_err = expected_rel_clustering_coefficient_error();
    printf("done\n");
    array<double,16> ret = {adj_err[0], adj_err[1], adj_err[2], adj_err[3], deg_err[0], deg_err[1], deg_err[2], deg_err[3], deg_err[4], deg_err[5], deg_err[6], rel_deg_err[0], rel_deg_err[1], rel_deg_err[2], rel_deg_err[3], clust_err};
    printf("done\n");
    return ret;
}

/*
 * Compute and print to stdout a number of statistics about the summary
 *
 */
vector<double> analyze_summary(bool approx_reconstruction_error, bool run_queries) {
    vector<double> ret;

    // Compute final density matrix
    vector<vector<double>> density;
    density.resize(num_supernodes);
    for (int i = 0; i < num_supernodes; ++i) {
        density[i].resize(num_supernodes);
    }
    // Since the supernodes id may not be consecutive, we need to find them
    // and give them consecutive indexes.
    unordered_map<int,int> id_to_index;
    int *index_to_id = new int[num_supernodes];
    int curr_index = 0;
    for (auto it = densities.begin(); it != densities.end(); ++it){
        pair<int,int> key = it->first;
        auto first_index_it = id_to_index.find(key.first);
        int first_index = 0;
        if (first_index_it  == id_to_index.end()) {
            first_index = curr_index;
            id_to_index.insert(make_pair(key.first, curr_index));
            index_to_id[curr_index] = key.first;
            curr_index++;
        } else {
            first_index = first_index_it->second;
        }

        auto second_index_it = id_to_index.find(key.second);
        int second_index = 0;
        if (id_to_index.find(key.second)  == id_to_index.end()) {
            second_index = curr_index;
            id_to_index.insert(make_pair(key.second, curr_index));
            index_to_id[curr_index] = key.second;
            curr_index++;
        } else {
            second_index = second_index_it->second;
        }
        density[first_index][second_index] = it->second;
    }

    // Print the partition.
    // For each supernode, print the list of vertices to stderr
    for (int i = 0; i < num_supernodes; ++i) {
        auto curr_supernode_it = supernodes.find(index_to_id[i]);
        for (int id : curr_supernode_it->second) {
            cerr << vertex_labels[id] << " "; // We print the original label for the vertex
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n"); // Double newline to separate from the next part
    delete[] index_to_id;

    // Print the final density matrix on stderr
    for (int i = 0; i < num_supernodes; ++i) {
        for (int j = 0; j < num_supernodes; ++j) {
            fprintf(stderr, "%lf ", density[i][j]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n"); // Double newline to separate from the next part

    // Print a 'rounded' version of the final density matrix
    if (num_supernodes < 100) {
        printf("majority matrix: \n");
        for (int i = 0; i < num_supernodes; ++i) {
            for (int j = 0; j < num_supernodes; ++j)
                printf("%i", density[i][j] >= 0.5);
            printf("\n");
        }
    }

    // Compute stats and print them
    double ind = index(), max_ind = avg_dens;
    printf("\nk = %i\nindex = %lf (%lf) max_index = %lf\n", num_supernodes, ind, ind / max_ind, max_ind);
    double orig_size = original_size(),
       comp_size = compressed_size(),
       ratio = comp_size / orig_size;
    printf("original size = %lf bits; compressed size = %lf bits (%lf)\n", orig_size, comp_size, ratio);

    // L1_max_err is the error that you get (for a 0-1 graph) if you use a summary with k = 1.
    // XXX is it obvious that it is the max?
    double L1_max_err = 2 * avg_dens * (1 - avg_dens);
    double L1_err = fast_l1_reconstruction_error();
    printf("L1 reconstruction error with expected adjacency matrix = %lf (%lf)\n", L1_err, L1_err / L1_max_err);

    double L1_maj_err = fast_l1_reconstruction_error_with_majority();
    printf("L1 reconstruction error with majority matrix = %lf (%lf)\n", L1_maj_err, L1_maj_err / min(avg_dens, 1 - avg_dens));

    // L2_max_err is the error that you get (for a 0-1 graph) if you use a summary with k = 1.
    // XXX is it obvious that it is the max?
    double L2_max_err = avg_dens * (1 - avg_dens);
    double L2_err = fast_l2_reconstruction_error();
    printf("Squared L2 reconstruction error = %lf (%lf)\n", L2_err, L2_err / L2_max_err);

#ifdef DEBUG
    assert(fabs(L1_err - 2 * L2_err) < 1e-8);
    assert(fabs(ind - (avg_dens - L2_err)) < 1e-8);
#endif

    // Free some memory as the SDP solver may need it
    //adj2.clear(); adj2.shrink_to_fit();//debug

    double cut_norm_err = approx_reconstruction_error ? cut_norm_error_approx(densities) : cut_norm_error(densities);
    cut_norm_err = fabs(cut_norm_err);
    printf("cut norm error = %lf (%lf)\n", cut_norm_err, cut_norm_err / avg_dens);


    // Add stats to the return list, so we can print them to stderr
    ret.push_back(ind);
    ret.push_back(max_ind);
    ret.push_back(comp_size);
    ret.push_back(orig_size);
    ret.push_back(ratio);
    ret.push_back(L1_err);
    ret.push_back(L1_max_err);
    ret.push_back(L2_err);
    ret.push_back(L2_max_err);
    ret.push_back(cut_norm_err);

    for (double stat : get_summary_statistics()) {
        ret.push_back(stat);
    }

    if (run_queries) {
	printf("Running queries...\n");
        double triangle_number_err = expected_triangles_error();
        printf("Triangle error = %lf, (%lld, %lf, %lf)\n", triangle_number_err, triangles, expected_triangles, triangle_number_err / triangles);
        //array<double, 7> degree_error_stats = expected_degree_error();
        //printf("Average degree errors: %lf %lf\n", degree_error_stats[0], degree_error_stats[4]);
        
        for (double stat: get_query_errors()) {
            ret.push_back(stat);
        }
    } else {
	printf("Not running queries.\n");
	for (int i = 0; i < 16; ++i) {
		ret.push_back(-1);
	}
    }

    return(ret);
}

/* 
 * Read a summary by parsing own output.
 *
 * Assumes read_graph() has been called somewhere else.
 *
 */
void read_summary() {
    int blank_lines = 0;
    int curr_supernode_id = 0;
    FILE* res_file = NULL;               // summ stderr result file

    res_file = fopen(summary_filename, "r");
    if (res_file == NULL) {
        perror("Error opening input file");
        exit(errno);
    }

    labels.assign(n,-1);

    while (blank_lines < 2) {
        char *line = (char *) malloc(RES_MAXLINE);
        fgets(line, RES_MAXLINE, res_file);
        if (strlen(line) == 1) {
            if (blank_lines == 1) {
                //printf("%d\n", curr_supernode_id);
                break;
            } else {
#ifdef DEBUG
                for (int i = 0; i < n; i++)
                    assert(labels[i] != -1);
                assert(verify_ids());
                assert(verify_labels());
                assert(vertex_labels.size() == n);
#endif
                blank_lines = 1;
                curr_supernode_id = 0; // reset counter
                continue;
            }
        }
        if (blank_lines == 0) { // Reading the summary partition
            supernode curr;
            char *token;
            for (; (token = strsep(&line, " \n")) != NULL;) {
               if (*token != '\0') {
                   vertex id = vertex_id(atoll(token));
#ifdef DEBUG
                   assert(labels[id] == -1);
#endif
                   labels[id] = curr_supernode_id;
                   curr.insert(id);
               }
            }
            supernodes.insert(make_pair(curr_supernode_id, curr));
        } else { // blank_lines == 1. Reading the densities matrix
            int j = 0;
            char *token;
            for (; (token = strsep(&line, " \n")) != NULL;) {
               if (*token != '\0') {
                    double dens = strtod(token, NULL);
                    //printf("%d %d %f\n", curr_supernode_id, j, dens);
                    pair<int,int> coordinates = make_pair(curr_supernode_id, j);
                    densities.insert(make_pair(coordinates, dens));
                    j++;
               }
            }
        }
        curr_supernode_id++;
        free(line);
    } // END WHILE
#ifdef DEBUG
    assert(verify_labels());
    assert(verify_ids());
    // Commented out because it requires too much precision
    //assert(verify_densities(densities));
#endif
    char *line = (char *) malloc(RES_MAXLINE);
    // Go to last line
    while ((fgets(line, RES_MAXLINE, res_file)) != NULL) {
        if (feof(res_file)) {
            break;
        }
    }
    if (line != NULL) {
        char *linecopy = line;
        char *token;
        int count = 0;
        for (; (token = strsep(&line, ", \n")) != NULL && count < 34;) {
            if (*token != '\0') {
                count++;
            }
        }
        token = strsep(&line, ", \n"); // skip empty token.
        double expected_triangles_error = strtod(token, NULL);
        expected_triangles = triangles + expected_triangles_error;
        num_supernodes = curr_supernode_id;
        free(linecopy);
    } else {
        fclose(res_file);
        perror("Error reading summary file");
        exit(errno);
    }
    fclose(res_file);
}

void print_summary_as_full_matrix() {
    for (int i = 0; i < n; ++i) {
	int label_i = labels[i];
        for (int j = 0; j < n; ++j) {
	    int label_j = labels[j];
	    printf("%f ", densities[make_pair(label_i,label_j)]);
	}
	printf("\n");
    }
}

#endif

