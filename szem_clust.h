#ifndef SZEMCLUST_H
#define SZEMCLUST_H

#include <chrono>
#include <cmath>
#include <set>
#include <vector>
using namespace std;

#include "graph.h"
#include "summ.h"
#include "util.h"

#define MIN_IMPROVEMENT 1e-3        // Minimum relative improvement from one iteration to the next in greedy k-means / k-median. If lower than this, stop.

#define MINI_BATCH_SIZE 1000        // Size for mini-batch clustering algorithm.

int p_exponent;
int error_type;
int dimension;                      // number of coordinates used in column sampling
int max_iters;

#include "sparse.h"

void prepare_matrix() {
    matrix.resize(n);
    dimension = n;
    modp.resize(n);
    if (error_type != ERROR_LC) {
        for (vertex v = 0; v < n; ++v) {
            for (vertex w : adj[v])
                matrix[v].push_back(Entry(w, 1));
            modp[v] = adj[v].size();
        }
    } else {
        // Divide the squared adjacency matrix by n.
        for (vertex v = 0; v < n; ++v) {
            for (auto p : adj2[v]) {
                double val = (double)p.second / n;
                matrix[v].push_back(Entry(p.first, val));
                modp[v] += fabs(val);
            }
            adj2[v].clear();
        }
        adj2.clear(); adj2.shrink_to_fit();
    }
}

// The following are only meaningful as long as the cluster centers are one of the points
// (i.e., in initialization and local improvements)
vector<vertex> facilities;
vector<int> fac_index;

struct FacilityHeap;
vector<FacilityHeap> facilityHeaps;

Point find_center_l1(const supernode& sup) {
    vector<vector<double>> m(dimension);
    for (vertex v : sup)
        for (Entry e : matrix[v]) 
            m[e.vert].push_back(e.val);

    Point center(dimension);
    for (int i = 0; i < dimension; ++i) if (!m[i].empty())
        center[i] = median(m[i]);
    return center;
}

Point find_center_l2(const supernode& sup) {
    Point center(dimension);
    for (vertex v : sup)
        for (Entry e : matrix[v]) 
            center[e.vert] += e.val;

    for (int i = 0; i < n; ++i)
        center[i] /= sup.size();
    return center;
}

Point find_center(const supernode& sup) {
    return error_type != ERROR_L1 ? find_center_l2(sup) : find_center_l1(sup);
}

void assign_distance_measure(int p) {
    assert(p == 1 || p == 2);
    p_exponent = p;
}

void compute_supernodes() {
    // Create supernodes from labels
    supernodes.clear();
    for (vertex v = 0; v < n; ++v)
        supernodes[labels[v]].insert(v);
}

bool compute_nonempty_supernodes() {
    bool empty, change = false;
    do {
        compute_supernodes();
        empty = false;
        for (int i = 0; i < num_supernodes; ++i)
            if (supernodes[i].empty()) {
                change = empty = true;
                printf("  empty cluster!\n");
                labels[rand() % n] = i;
            }
    } while (empty);
    return change;
}

/*
 * Initialize centers for clusters.
 * Uses kmeans++ and guarantees an (16+8 ln k)-competitive clustering.
 *
 */
void initialize_facilities() {
    vector<vertex> candidates = sample_without_replacement(n, n);
    labels.resize(n, -1);
    fac_index.resize(n, -1);

    FenwickTree tree(n);
    double inf = 4 * n + 1;          // infinity, but don't want to use too large values because we will use subtractions
    for (int i = 0; i < n; ++i)
        tree.update(i, inf);

    for (int i = 0; i < num_supernodes; ++i) {
        if (i % 10 == 0)
            printf(" ...%i ", i);

        int bestv = choose_from_distribution(tree);
//        printf("bestv = %i d = %lf %lf %i\n", bestv, tree.value(bestv), n * tree.value(bestv), p_exponent);
        labels[bestv] = i;
        facilities.push_back(bestv);
        fac_index[bestv] = facilities.size() - 1;
        tree.set_value(bestv, 0);

        static double factor = 1.0 / pow(2, p_exponent);
        Point p = get_point(bestv);
        vector<double> bound(num_supernodes);
        for (int j = 0; j < i; ++j)
            bound[j] = factor * quick_distance_p(facilities[j], p, modp[bestv]);

        for (vertex v = 0; v < n; ++v) {
            // The probability is proportional to quick_distance^p. For 0-1 matrices this is the same as just quick_distance_l1,
            // so we avoid taking roots.
            double prev = tree.value(v);
            if (prev == 0) continue;

            if (i == 0 || prev > bound[labels[v]]) {
                double d = quick_distance_p(v, p, modp[bestv]);
                if (d < prev) {
                    tree.update(v, d - prev);
                    labels[v] = i;
                }
            }
        }
    }

    double cost = tree.partial_sum(n - 1);
    printf("\n    cost = %.11lf final_estimate = %lf (%lf) avg_cost = %.11lf\n", cost, cost / n / n, cost / (avg_deg * (n - avg_deg)) / (p_exponent == 1 ? 2 : 1), cost / n);

    compute_nonempty_supernodes();
}

/*
void accelerated_lp_cluster_old() {
    puts("allocating centers");
    vector<Point> centers;              // point associated with each cluster center
    for (vertex v: facilities)
        centers.push_back(get_point(v));

    vector<vector<int>> center_sup(num_supernodes);
    vector<double> modc(num_supernodes);

    for (int i = 0; i < num_supernodes; ++i) {
        center_sup[i].clear();
        modc[i] = 0;
        for (int j = 0; j < dimension; ++j)
            if (centers[i][j] != 0) {
                center_sup[i].push_back(j);
                modc[i] += centers[i][j] * centers[i][j];
            }
    }

    // Compute intra-center distances
    puts("intra-center distances");
    Matrix<double> center_dist(num_supernodes, vector<double>(num_supernodes));
    for (int i = 0; i < num_supernodes; ++i)
        for (int j = i; j < num_supernodes; ++j)
            center_dist[i][j] = center_dist[j][i] = 
                quick_distance(facilities[i], centers[j], modc[j]);
                //quick_distance(facilities[i], facilities[j]);

    puts("bounds");
    // Assign each point to its closest facility
    vector<double> upper(n, infinity);
    Matrix<double> lower(n, vector<double>(num_supernodes, 0));

    labels.resize(n, -1);
    for (int i = 0; i < num_supernodes; ++i)
        for (vertex v = 0; v < n; ++v)
            if (labels[v] < 0 || upper[v] > center_dist[i][labels[v]] / 2) {
                double d = quick_distance(facilities[i], v);
                if (d < upper[v]) {
                    upper[v] = d;
                    labels[v] = i;
                }
                lower[v][i] = d;
            } else {
                lower[v][i] = upper[v];     // why is this not in the paper?
            }

    // Ensure non-empty clusters even when some distances are zero
    for (int i = 0; i < num_supernodes; ++i)
        labels[facilities[i]] = i;

    double prev_cost = infinity;
    for (int it = 0; it < max_iters; ++it) {
        printf("  Iteration: %d\n", it);
        printf("  average center support = %lf\n", accumulate(modc.begin(), modc.end(), 0.0) / num_supernodes);

        // Compute intracenter distances and half the distance between each center and another center (sc[])
        if (it != 0) {
            for (int i = 0; i < num_supernodes; ++i)
                for (int j = i; j < num_supernodes; ++j) {
                    double d = 0;
                    for (int k : center_sup[i])
                        if (centers[j][k])
                            d += centers[i][k] * centers[j][k];
                    d = sqrt(max(0, modc[i] + modc[j] - 2 * d));
                    center_dist[i][j] = center_dist[j][i] = d;
                }
        }

        vector<double> sc(num_supernodes, infinity);
        for (int i = 0; i < num_supernodes; ++i) {
            for (int j = 0; j < num_supernodes; ++j)
                if (i != j && center_dist[i][j] < sc[i])
                    sc[i] = center_dist[i][j];
            sc[i] /= 2;
        }

//for (int i = 0; i < num_supernodes; ++i) { int nz = 0; for (int j = 0; j < n; ++j) if (fabs(centers[i][j]) > 1e-4) ++nz; printf("%i: %i\n", i, nz); }//debug

        for (vertex x = 0; x < n; ++x) {
            if (x % 1000 == 0)
                printf("  %i", x);
            if (upper[x] > sc[labels[x]]) {
                for (int c = 0; c < num_supernodes; ++c)
                    if (c != labels[x] && upper[x] > lower[x][c] && upper[x] > center_dist[labels[x]][c] / 2) {
                        double d = quick_distance(x, centers[c], modc[c]);
                        lower[x][c] = d;
                        if (lower[x][c] < upper[x]) {
                            labels[x] = c;
                            upper[x] = lower[x][c];
                        }
                    }
            }
        }
        puts("");

        // Find new centers and reassign empty clusters
        bool empty = false;
        compute_supernodes();
        for (int c = 0; c < num_supernodes; ++c)
            if (supernodes[c].empty()) {
                empty = true;
                printf("  empty cluster!\n");

                int v = rand() % n;
                for (int x = 0; x < n; ++x)
                    lower[x][c] = max(lower[x][labels[v]] - upper[v], 0);
                labels[v] = c;
                upper[v] = 0;
            }
        if (empty) continue;

        for (int c = 0; c < num_supernodes; ++c) {
//printf("%i (%i)..", c, (int)supernodes[c].size());//debug
            Point mc = find_center(supernodes[c]);

            double d_mc = pref_distance(centers[c], mc);
            for (vertex x = 0; x < n; ++x)
                lower[x][c] = max(lower[x][c] - d_mc, 0);
            centers[c] = mc;
        }

        for (int i = 0; i < num_supernodes; ++i) {
            center_sup[i].clear();
            modc[i] = 0;
            for (int j = 0; j < dimension; ++j)
                if (centers[i][j] != 0) {
                    center_sup[i].push_back(j);
                    modc[i] += centers[i][j] * centers[i][j];
                }
        }

        // Note: unlike Elkan, we compute exact upper values here. We can afford this, as it's n distance computations.
        // This also simplified some stuff in step 3.
        double cost = 0;
        for (vertex x = 0; x < n; ++x) {
//printf("%i->c[%i]\n", x, labels[x]); // debug
            upper[x] = quick_distance(x, centers[labels[x]], modc[labels[x]]);//pref_distance(get_point(x), centers[labels[x]]);
            lower[x][labels[x]] = upper[x];
            cost += pow(upper[x], p_exponent);
        }

        printf("    cost = %.11lf final_estimate = %lf (%lf) avg_cost = %.11lf factor = %lf\n",
                cost, cost / n / n, cost / (avg_deg * (n - avg_deg)) / (p_exponent == 1 ? 2 : 1), cost / n, cost / prev_cost);

        if (cost >= (1.0 - MIN_IMPROVEMENT) * prev_cost)
            break;
        prev_cost = cost;
    }

    compute_nonempty_supernodes();
//    reducer.clear();//debug
#ifdef DEBUG
    verify_labels();
#endif
}
*/

void accelerated_lp_cluster() {
    vector<double> modc(num_supernodes);

    for (int i = 0; i < num_supernodes; ++i) {
        modc[i] = 0;
        for (Entry e : matrix[facilities[i]]) {
            //int j = e.vert; 
	    double val = e.val;
            modc[i] += val * val;
        }
    }

    // Compute intra-center distances
    puts("intra-center distances");
    Matrix<double> center_dist(num_supernodes, vector<double>(num_supernodes));
    for (int i = 0; i < num_supernodes; ++i) {
        Point centeri = get_point(facilities[i]);
        for (int j = i; j < num_supernodes; ++j)
            center_dist[i][j] = center_dist[j][i] = 
                quick_distance(facilities[j], centeri, modc[i]);
    }

    puts("bounds");
    // Assign each point to its closest facility
    vector<double> upper(n, infinity);
    Matrix<float> lower(n, vector<float>(num_supernodes, 0));
    puts("initializing bounds");

    labels.resize(n, -1);
    for (int i = 0; i < num_supernodes; ++i)
        for (vertex v = 0; v < n; ++v)
            if (labels[v] < 0 || upper[v] > center_dist[i][labels[v]] / 2) {
                double d = quick_distance(facilities[i], v);
                if (d < upper[v]) {
                    upper[v] = d;
                    labels[v] = i;
                }
                lower[v][i] = d;
            } else {
                lower[v][i] = upper[v];     // why is this not in the paper?
            }

    // Ensure non-empty clusters even when some distances are zero
    for (int i = 0; i < num_supernodes; ++i)
        labels[facilities[i]] = i;

    puts("allocating centers");
    vector<SparsePoint> centers2;
    for (vertex v: facilities)
        centers2.push_back(get_sparse_point(v));

    double prev_cost = infinity;
    for (int it = 0; it < max_iters; ++it) {
        printf("  Iteration: %d\n", it);
//        printf("  average center support = %lf\n", accumulate(modc.begin(), modc.end(), 0.0) / num_supernodes);

        puts("    computing intracenter distances...");
        // Compute intracenter distances and half the distance between each center and another center (sc[])
        if (it != 0) {
            for (int i = 0; i < num_supernodes; ++i) {
                for (int j = i; j < num_supernodes; ++j) {
                    center_dist[i][j] = center_dist[j][i] = euclidean_dist(centers2[i], centers2[j]);
                }
            }
        }

        vector<double> sc(num_supernodes, infinity);
        for (int i = 0; i < num_supernodes; ++i) {
            for (int j = 0; j < num_supernodes; ++j)
                if (i != j && center_dist[i][j] < sc[i])
                    sc[i] = center_dist[i][j];
            sc[i] /= 2;
        }

//for (int i = 0; i < num_supernodes; ++i) { int nz = 0; for (int j = 0; j < n; ++j) if (fabs(centers[i][j]) > 1e-4) ++nz; printf("%i: %i\n", i, nz); }//debug

        puts("    assigning points to clusters...");
        for (vertex x = 0; x < n; ++x) {
            SparsePoint px;
            if (x % 5000 == 0)
                printf("      %i", x);

            if (upper[x] > sc[labels[x]]) {
                for (int c = 0; c < num_supernodes; ++c)
                    if (c != labels[x] && upper[x] > lower[x][c] && upper[x] > center_dist[labels[x]][c] / 2) {
                        if (px.vec.empty()) px = get_sparse_point(x);
                        double d = euclidean_dist(px, centers2[c]);
                        lower[x][c] = d;
                        if (lower[x][c] < upper[x]) {
                            labels[x] = c;
                            upper[x] = lower[x][c];
                        }
                    }
            }
        }
        puts("");

        // Reassign empty clusters
        puts("    assigning empty clusters...");
        bool empty = false;
        compute_supernodes();
        for (int c = 0; c < num_supernodes; ++c)
            if (supernodes[c].empty()) {
                empty = true;
                printf("  empty cluster!\n");

                int v = rand() % n;
                for (int x = 0; x < n; ++x)
                    lower[x][c] = max(lower[x][labels[v]] - upper[v], 0);
                labels[v] = c;
                upper[v] = 0;
            }
        if (empty) continue;

        // Find new centers
        puts("    finding new centers...");
        for (int c = 0; c < num_supernodes; ++c) {
            if (c % 100 == 0) printf("      %i", c);
//printf("%i (%i)..", c, (int)supernodes[c].size());//debug
            SparsePoint mc = make_sparse(find_center(supernodes[c]));

            double d_mc = euclidean_dist(centers2[c], mc);
//            double d_mc = pref_distance(centers[c], mc);
            for (vertex x = 0; x < n; ++x)
                lower[x][c] = max(lower[x][c] - d_mc, 0);
            centers2[c] = mc;
        }
        puts("");

        // Note: unlike Elkan, we compute exact upper values here. We can afford this, as it's n distance computations.
        // This also simplified some stuff in step 3.
        puts("    recomputing distances...");
        double cost = 0;
        for (int c = 0; c < num_supernodes; ++c) {
            if (c % 100 == 0) printf("      %i", c);

            Point p = centers2[c].full_point(n);
            double m = centers2[c].mod22();

            for (vertex x : supernodes[c]) {
                upper[x] = quick_distance(x, p, m);
                lower[x][labels[x]] = upper[x];
                cost += pow(upper[x], p_exponent);
            }
        }

        /*
        for (vertex x = 0; x < n; ++x) {
            if (x % 1000 == 0)
                printf("      %i [%i %lf]", x, matrix[x].size(), centers2[labels[x]].vec.size());
//printf("%i->c[%i]\n", x, labels[x]); // debug
            upper[x] = euclidean_dist(get_sparse_point(x), centers2[labels[x]]);
//            upper[x] = quick_distance(x, centers[labels[x]], modc[labels[x]]);//pref_distance(get_point(x), centers[labels[x]]);
            lower[x][labels[x]] = upper[x];
            cost += pow(upper[x], p_exponent);
        }
        */
        puts("");

        printf("    cost = %.11lf final_estimate = %lf (%lf) avg_cost = %.11lf factor = %lf\n",
                cost, cost / n / n, cost / (avg_deg * (n - avg_deg)) / (p_exponent == 1 ? 2 : 1), cost / n, cost / prev_cost);

        if (cost >= (1.0 - MIN_IMPROVEMENT) * prev_cost)
            break;
        prev_cost = cost;
    }

    compute_nonempty_supernodes();
//    reducer.clear();//debug
#ifdef DEBUG
    verify_labels();
#endif
}


// min heap: dist[i / 2] <= dist[i], i = 1..N
struct FacilityHeap {
    vertex src;
    int N;

    typedef pair<int, double> pid;
    vector<pid> heap;                // heap[1], ..., heap[k] = facility indices
    vector<int> pos;                 // heap[pos[i]] = i for i=0..num_supernodes-1

    FacilityHeap(int max_size = num_supernodes) : src(-1), N(0), heap(max_size + 1), pos(max_size) {}

    int closest_label() { return heap[1].first; }

    int closest() { return facilities[heap[1].first]; }

    double dist(int index) { return quick_distance_p(src, facilities[index]); }

    double lowest_cost() { return heap[1].second; }

    double second_lowest_cost() {
        if (N >= 3)
            return min(heap[2].second, heap[3].second);
        else if (N >= 2)
            return heap[2].second;
        else
            return infinity;
    }

    // Cost if we replace facility[i] with f
    double new_cost(int i, vertex f) {
        if (heap[1].first != i)
            return min(lowest_cost(), quick_distance_p(src, f));
        return min(second_lowest_cost(), quick_distance_p(src, f));
    }

    void sift_up(int i, int x, double d) {
        while (i > 1 && heap[i / 2].second > d) {
            heap[i] = heap[i / 2];
            pos[heap[i].first] = i;
            i /= 2;
        }
        heap[i] = pid(x, d);
        pos[x] = i;
    }

    void sift_down(int i, int x, double d) {
        while (2 * i <= N) {
            i *= 2;
            if (i < N && heap[i + 1].second <= heap[i].second) ++i;

            if (heap[i].second < d) {
                heap[i / 2] = heap[i];
                pos[heap[i / 2].first] = i / 2;
            } else {
                i /= 2;
                break;
            }
        }
        heap[i] = pid(x, d);
        pos[x] = i;
    }

    void put_into_place(int x) {
        int index = pos[x];
        double d = dist(x);
        if (index > 1 && heap[index / 2].second > d)
            sift_up(index, x, d);
        else
            sift_down(index, x, d);
    }

    void heapify() {
        for (int i = N / 2; i >= 1; --i)
            sift_down(i, heap[i].first, heap[i].second);
    }

    bool check_heap() {
        if (!(heap.size() >= N + 1)) return false;
        for (int i = 1; i <= N; ++i)
            if (!(heap[i].second == dist(heap[i].first))) {
                printf("v = %i i = %i %lf %lf\n", src, i, dist(heap[i].first), heap[i].second);
                return false;
            }
        for (int i = N; i > 1; --i)
            if (!(heap[i / 2].second <= heap[i].second))
                return false;

        if (!(pos.size() == N)) return false;
        for (int i = N; i > 0; --i)
            if (!(pos[heap[i].first] = i)) return false;

        return true;
    }

    void clear() {
        vector<pid>().swap(heap);
        vector<int>().swap(pos);
    }
};

// Cost of assigment to facilities[], using precomputed FacilityHeaps
double assignment_cost() {
    double ret = 0;
    for (vertex v = 0; v < n; ++v)
        ret += facilityHeaps[v].lowest_cost();
    return ret;
}

double assignment_cost_after_swap(int i, vertex f) {
    double ret = 0;
    for (vertex v = 0; v < n; ++v)
        ret += facilityHeaps[v].new_cost(i, f);
    return ret;
}

void open_all_facilities() {
    facilityHeaps.resize(n);
    printf("  Opening facilities... ");
    for (vertex v = 0; v < n; ++v) {
        facilityHeaps[v].src = v;
        facilityHeaps[v].N = num_supernodes;
        for (int i = 0; i < num_supernodes; ++i) {
            facilityHeaps[v].heap[i + 1] = make_pair(i, facilityHeaps[v].dist(i));
            facilityHeaps[v].pos[i] = i + 1;
        }
        facilityHeaps[v].heapify();
    }
    printf("done!\n");
}

// Swap facility f with vertex v
bool try_improvement(int index, vertex new_f, double& cost, int& improvements, int& tried, int &p1, int p2) {
    if (fac_index[new_f] >= 0) return false;    // it makes no sense to swap two facilities
    ++tried;
    int prev = facilities[index];

    double next_cost = assignment_cost_after_swap(index, new_f);
    if (next_cost < cost) {
        cost = next_cost;
        facilities[index] = new_f;
        fac_index[new_f] = index;
        for (vertex v = 0; v < n; ++v)
            facilityHeaps[v].put_into_place(index);
        ++improvements;
        printf("    #%i/%i: cost = %lf   avg_cost = %lf (swapped f[%i] = %i with v = %i; pair %i/%i) %i\n",
                improvements, tried, cost, cost / n, index, prev, new_f, p1, p2, max_att_impr);
        return true;
    } else {
        return false;
    }
}

bool try_center_improvements(double& cost, int& improvements, int& tried) {
    printf("  center improvements...\n");

    // Find supernodes
    supernodes.clear();
    for (vertex v = 0; v < n; ++v) {
        int f = facilityHeaps[v].closest();
        supernodes[f].insert(v);
    }

    bool improved = false;
    for (int index = 0; tried < max_att_impr && index < facilities.size(); ++index) {
        // Find point within each supernode closest to its center
        int f = facilities[index];
        Point center = find_center(supernodes[f]);

        int best = 0;
        double mind = infinity;
        for (vertex v : supernodes[f]) {
            double d = pref_distance_to_p(center, get_point(v));
            if (d < mind) {
                best = v;
                mind = d;
            }
        }

        if (try_improvement(index, best, cost, improvements, tried, index, facilities.size()))
            improved = true;
    }
    return improved;
}

bool try_other_improvements(double& cost, int& improvements, int& tried) {
    printf("  other improvements...\n");

    // Find a random permutation of pairs (facility, vertex)
    vector<pair<int, vertex>> pairs;
    for (int index = 0; index < facilities.size(); ++index)
        for (vertex v = 0; v < n; ++v) {
//            if ((*dist[v].begin()).second != facilities[index])    // this would eliminate center improvements
                pairs.push_back(make_pair(index, v));
        }
    random_shuffle(pairs.begin(), pairs.end());

    bool improved = false;
    for (int pi = 0; tried < max_att_impr && pi < pairs.size(); ++pi) {
        auto p = pairs[pi];
        int index = p.first;
        vertex v = p.second;

        if (try_improvement(index, v, cost, improvements, tried, pi, pairs.size()))
            improved = true;
    }
    return improved;
}

// Approximation algorithm for k-medians and k-means
// In practice we limit the number of improvements for faster running time.
void local_improvements() {
    open_all_facilities();
    double cost = assignment_cost();
    printf("  Cost of assignment = %lf   avg_cost = %lf\n", cost, cost / n);

    int improvements = 0, tried = 0;
    for (;;) {
        bool improved = false;
        // Note that we try center improvements while possible, emulating k-means
        if (try_center_improvements(cost, improvements, tried) ||
            try_other_improvements(cost, improvements, tried)) {
            improved = true;
            if (tried >= max_att_impr) break;
        }
        if (!improved) {
            puts("Done! Guaranteed 5-approximation to clustering problem!\n");
            break;
        }
    }

    for (vertex v = 0; v < n; ++v)
        labels[v] = facilityHeaps[v].closest_label();
    compute_nonempty_supernodes();

    facilityHeaps.clear();
    facilityHeaps.shrink_to_fit();
}

// Take a partition and compute the associated densities. Assumes consecutively numbered supernodes.
void compute_densities() {
    for (vertex v = 0; v < n; ++v)
        for (vertex w : adj[v])
            ++densities[make_pair(labels[v], labels[w])];

    for (int i = 0; i < num_supernodes; ++i) {
        for (int j = 0; j < num_supernodes; ++j) {
            int a = supernodes[i].size(), b = supernodes[j].size();
            double t = (double)a * b;
            if (t != 0)
                densities[make_pair(i,j)] /= t;
//            printf("%i,%i: %lf (%i, %i)\n", i, j, densities[make_pair(i, j)], a, b);
        }
    }
#ifdef DEBUG
    assert(verify_densities(densities));
#endif
}

void random_summary() {
    labels.resize(n);
    for (int i = 0; i < n; ++i)
        labels[i] = rand() % num_supernodes;
    for (int i = 0; i < n; ++i)
        supernodes[labels[i]].insert(i);
}

// Create the partition using k-means or k-medians
//
// The appropriate algorithm and function is chosen according to the parameters
//
// Return a vector containing at position 0 the time needed to perform
// clustering+density computation, and at position 1 the time needed to square
// the adjacency matrix (or 0).
//
//
vector<double> summarize_by_clustering(int k, int _error_type, __attribute__((unused)) bool use_mini_batch, bool use_approx_alg, bool do_use_random_summary, __attribute__((unused)) int JL_dims, int iters) {
    max_iters = iters;
//srand(0);//debug
    vector<double> ret(2);
    std::chrono::steady_clock::time_point start_clustering, end_clustering, start_squaring, end_squaring, start_jl, end_jl;

    num_supernodes = min(k, n);
    error_type = _error_type;
    if (error_type == ERROR_L2)
        assign_distance_measure(2);
    else
        assign_distance_measure(1);

    if (error_type == ERROR_LC) {
        // Square the matrix. We divide it by n in compute_points()
        printf("Squaring matrix...\n");
        start_squaring = std::chrono::steady_clock::now();
        if (error_type == ERROR_LC)
            square_matrix(); //debug
        else
            printf("  Skipping!\n");
        end_squaring = std::chrono::steady_clock::now();
        // Compute elapsed time
        auto elapsed_squaring = std::chrono::duration_cast<std::chrono::milliseconds>(end_squaring - start_squaring);
        ret[1] = elapsed_squaring.count() / 1000.0; // Divide to get seconds
        printf("Matrix Squaring Elapsed time: %f secs\n\n", ret[1]);

        start_clustering = chrono::steady_clock::now();


        /*
        // Project to lower dimensions
        if (JL_dims >= 0) {
            printf("Reducing dimensions (to %d)...\n", JL_dims);
            start_jl = chrono::steady_clock::now();
            if (error_type == ERROR_L1)
                JL_dims = n;

            if (error_type == ERROR_LC) {
                puts("no");
                exit(1);//debug
    //            compute_points_adj2(JL_dims);
            }
            else {
    //            compute_points_adj(JL_dims);//debug
            }
            end_jl = chrono::steady_clock::now();
            auto elapsed_jl = std::chrono::duration_cast<std::chrono::milliseconds>(end_jl - start_jl);
            printf("JL Elapsed time: %f secs\n\n", elapsed_jl.count() / 1000.0);
        } else
            printf("  Skipping dimensionality reduction!\n\n");
            */
//        quick_distance = quick_distance_general;            // override
//        quick_distance_to_p = quick_distance_to_p_general;  // override
    }

    start_clustering = chrono::steady_clock::now();

    if (do_use_random_summary) {
        puts("Using random summary!"); // useless, just for testing
        random_summary();
    } else {
        prepare_matrix();
        puts("Clustering...");

        printf("k-means++ initialization...\n");
        auto start_init_cl = chrono::steady_clock::now();
        initialize_facilities();
        auto end_init_cl = chrono::steady_clock::now();
        auto elapsed_init_cl = std::chrono::duration_cast<std::chrono::milliseconds>(end_init_cl - start_init_cl);
        printf("Clustering Initialization Elapsed time: %lf secs\n\n", elapsed_init_cl.count() / 1000.0);

        if (use_approx_alg) {
            printf("Local improvements...\n");
            auto start_local_impr = chrono::steady_clock::now();
            local_improvements();
            auto end_local_impr = chrono::steady_clock::now();
            auto elapsed_local_impr = std::chrono::duration_cast<std::chrono::milliseconds>(end_local_impr - start_local_impr);
            printf("Local Improvements Elapsed time: %lf secs\n\n", elapsed_local_impr.count() / 1000.0);
        }

        printf("k-means / k-median...\n");
        if (max_iters > 0)
            accelerated_lp_cluster();
    }

    // Compute densities
    printf("Computing densities...\n");
    compute_densities();
    end_clustering = std::chrono::steady_clock::now();

    // Compute elapsed time
    auto elapsed_clustering = std::chrono::duration_cast<std::chrono::milliseconds>(end_clustering - start_clustering);
    ret[0] = elapsed_clustering.count() / 1000.0; // Divide to get seconds
    printf("Clustering Elapsed time: %f secs\n\n", ret[0]);

    return ret;
}

#endif
