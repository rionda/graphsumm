#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>
#include <random>
#include <numeric>

#include "fenwick.h"
using namespace std;

typedef long long ll;
const double infinity = 1e150; // the square must not overflow

typedef int vertex;

template<typename T>
using Matrix = vector<vector<T>>;

typedef vector<double> Point;
typedef unordered_map<vertex, int> SparseRow;

/*
 * Generate random numbers
 *
 */
double rand01() { return (double)rand() / RAND_MAX; }
double randab(double a, double b) { return a + rand01() * (b - a); }

/* 
 * Remove the element at index 'index' from a vector. May reorder elements.
 *
 */
template<typename T>
void remove_element(vector<T>& v, int index) {
    swap(v[index], v[v.size() - 1]);
    v.pop_back();
}

/* 
 * Compute the squared Frobenius norm of a matrix (i.e., it's l2 norm)
 *
 */
template<typename T>
double frobenius_squared(const Matrix<T>& m) {
    double sum = 0;
    for (int i = 0; i < m.size(); ++i)
        for (double a : m[i])
            sum += a * a;
    return sum;
}

/*
 * Print a Point to stdout (component by component)
 *
 */
void print(const Point& v) {
    for (int i = 0; i < v.size(); ++i) printf("%g ", v[i]);
    printf("\n");
}

/*
 * Print a vector of integers to stdout (component by component)
 *
 */
void print(const vector<int>& v) {
    for (int i = 0; i < v.size(); ++i) printf("%i ", v[i]);
    printf("\n");
}

/*
 * Print a vector of points to stdout (component by component)
 *
 */
void print(const vector<Point>& v) { for (Point p : v) print(p); }

vector<int> random_permutation(int n) {
    vector<int> v;
    for (int i = 0; i < n; ++i)
        v.push_back(i);
    random_shuffle(v.begin(), v.end());
    return v;
}

/*
 * Sample k numbers from [0,n) without replacement
 *
 */
vector<int> sample_without_replacement(int n, int k) {
    k = min(k, n);
    vector<int> ret;
    for (int i = 0; i < n; ++i) ret.push_back(i);
    random_shuffle(ret.begin(), ret.end());
    ret.resize(k);
    ret.shrink_to_fit();
    return ret;
}

/*
 * Sample k elements from a vector without replacement
 *
 */
template<typename T>
vector<T> sample_without_replacement(const vector<T>& S, int k) {
    vector<T> ret;
    vector<int> indices = sample_without_replacement(S.size(), k);
    for (int i : indices)
        ret.push_back(S[i]);
    return ret;
}

/*
 * Convert a Sparse Row to a vector of n dimensions
 *
 */
Point convert(int n, const SparseRow& r) {
    vector<double> ret(n);
    for (auto p : r)
        ret[p.first] = p.second;
    return ret;
}

/*
 * Compute the dot product of two Points
 *
 */
double dot_product(const Point& v, const Point& w) {
#ifdef DEBUG
    assert(v.size() == w.size());
#endif
    return inner_product(v.begin(), v.end(), w.begin(), 0.0);
}

/*
 * Compute the dot product between a SparseRow and a Point
 *
 *
 */
double dot_product(const SparseRow& s, const Point& w) {
    double ret = 0;
    if (s.size() > w.size()) {
        for (int i = 0; i < w.size(); ++i) {
            auto p = s.find(i);
            if (p != s.end())
                ret += p->second * w[i];
        }
    } else {
        for (auto p : s)
            ret += p.second * w[p.first];
    }
    return ret;
}

/*
 * Compute the dot product between a vector of int and a point
 *
 */
double dot_product(const vector<int> &v, const Point&w) {
#ifdef DEBUG
    assert(v.size() == w.size());
#endif
    double ret = 0.0;
    for (int i = 0; i < v.size(); i++) {
        ret += v[i] * w[i];
    }
    return ret;
}

/*
 * Compute the squared l2 distance between two points
 * Not a distance (satisfies doubled triangle inequality instead).
 *
 */
double l2_2_distance(const Point& v, const Point& w) {
    double ret = 0;
    for (int i = 0; i < v.size(); ++i) {
        double d = v[i] - w[i];
        ret += d * d;
    }
    return ret;
}

double l2_distance(const Point& v, const Point& w) {
    return sqrt(l2_2_distance(v, w));
}

/*
 * Compute the L1 distance between two points
 *
 */
double l1_distance(const Point& v, const Point& w) {
    double ret = 0;
    for (int i = 0; i < v.size(); ++i) {
        double d = v[i] - w[i];
        ret += fabs(d);
    }
    return ret;
}

/*
 * Normalize a vector in place
 *
 */
void normalize(vector<double>& v) {
    double length = sqrt(dot_product(v, v));
    for (int i = 0; i < v.size(); ++i)
        v[i] /= length;
}

/*
 * Gram-Schmidt orthogonalization
 *
 */
void orthonormalize(vector<Point>& v) {
    int d = v.size(), n = v[0].size();
    for (int i = 0; i < d; ++i) {
        normalize(v[i]);
        for (int j = i + 1; j < d; ++j) {
            double proj = dot_product(v[i], v[j]);
            for (int k = 0; k < n; ++k)
                v[j][k] -= proj * v[i][k];
        }
    }
}

/*
 * Create a random vector of n dimensions with entries distributed according to
 * a standard normal distribution
 */
Point gaussian_vector(int dimensions) {
    static default_random_engine generator;
    static normal_distribution<double> distribution(0, 1);

    Point ret;
    for (int j = 0; j < dimensions; ++j)
        ret.push_back(distribution(generator));
    return ret;
}

/*
 * Return a basis of a random d-dimensional subspace
 *
 */
vector<Point> random_basis(int dimensions, int d) {
    // Pick d random vectors with Gaussian coordinates
    vector<Point> v;
    for (int i = 0; i < d; ++i)
        v.push_back(gaussian_vector(dimensions));
    return v;
}

/*
 * Project a point on a space using a basis
 * This is actually matrix multiplication.
 */
template<typename T>
Point projection(const T& v, const vector<Point>& basis) {
    Point ret;
    for (int i = 0; i < basis.size(); ++i)
        ret.push_back(dot_product(v, basis[i]));
    return ret;
}

/*
 * Johnson-Lindenstrauss transform
 *
 */
struct JohnsonLindenstrauss {
    vector<Point> basis;

    // Go from orig_dims dimensions to d
    JohnsonLindenstrauss(int orig_dims, int d, bool database_friendly = true) {
        if (d < orig_dims) {
            if (database_friendly) {
                // 'Database-friendly random projections' by Achlioptas
                basis.resize(d, Point(orig_dims));
                for (int i = 0; i < d; ++i)
                    for (int j = 0; j < orig_dims; ++j)
                        basis[i][j] = ((rand() & 1) ? 1 : -1) / sqrt(d);
            } else {
                // Create an orthonormal random base using standard gaussians for the components
                basis = random_basis(orig_dims, d);
                orthonormalize(basis);

                // Scale the basis appropriately
                double scale = sqrt((double) orig_dims / d);
                for (Point &p : basis)
                    for (int i = 0; i < p.size(); ++i)
                        p[i] *= scale;
            }
        } else {
            basis.resize(orig_dims, Point(orig_dims));
            for (int i = 0; i < orig_dims; ++i) basis[i][i] = 1;
        }
    }

    template<typename T>
    Point transform(const T& p) {
        return projection(p, basis);
    }
};

struct DimReducer {
    int dim;
    bool do_reduction;
    double scale;
    vector<vector<bool>> rand_choices;

    DimReducer() {}

    DimReducer(int o, int d) { prepare(o, d); }

    void clear() { rand_choices.clear(); rand_choices.shrink_to_fit(); }

    void prepare(int orig_dims, int d) {
        do_reduction = (d < orig_dims);
        if (do_reduction) {
            // 'Database-friendly random projections' by Achlioptas
            dim = d;
            rand_choices.resize(d, vector<bool>(orig_dims));
            for (int i = 0; i < d; ++i)
                for (int j = 0; j < orig_dims; ++j)
                    rand_choices[i][j] = rand() & 1;
            scale = 1.0 / sqrt(d);
        } else {
            dim = orig_dims;
        }
    }

    Point transform(const SparseRow& p) {
        Point ret(dim);
        if (do_reduction) {
            for (int i = 0; i < dim; ++i)
                for (auto z : p) {
                    double f = rand_choices[i][z.first] ? scale : -scale;
                    ret[i] += f * z.second;
                }
        } else {
            for (auto z : p) ret[z.first] = 1;
        }
        return ret;
    }

    Point transform(const Point& p) {
        if (!do_reduction) return p;
        Point ret(dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < p.size(); ++j) {
                double f = rand_choices[i][j] ? scale : -scale;
                ret[i] += f * p[j];
            }
        return ret;
    }
};

/*
 * Reduce the dimensionality of a vector of points of dimension orig_dims to dimension d
 * using Johnson-Lindenstrauss transform
 *
 */
vector<Point> johnson_lindenstrauss(const vector<Point>& points, int orig_dims, int d) {
    vector<Point> ret;
    DimReducer reducer(orig_dims, d);
    for (int i = 0; i < points.size(); ++i)
        ret.push_back(reducer.transform(points[i]));
    return ret;
}

/*
 * Test function for the Johnson-Lindenstrauss transform
 *
 */
void test_jl() {
    const double eps = 0.2;
    const int m = 300;
    const int d = ceil(4 * log(m) / (pow(eps, 2) / 2 - pow(eps, 3) / 3));
    const int orig_dims = 2 * d;

    vector<Point> v(m);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < orig_dims; ++j)
            v[i].push_back(rand() % 100 - 50);

    double dist[m][m], dist2[m][m];

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
            dist[i][j] = l2_distance(v[i], v[j]);

    printf("eps = %lf m = %i orig_dims = %i d = %i\n", eps, m, orig_dims, d);
    for (;;) {
        vector<Point> v_jl = johnson_lindenstrauss(v, orig_dims, d);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j)
                dist2[i][j] = l2_distance(v_jl[i], v_jl[j]);

        double avg_err = 0, max_err = 0;
        int total = 0, worst_i = -1, worst_j = -1;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j) if (i != j) {
                double error = fabs(dist2[i][j] / dist[i][j] - 1);
                if (error > max_err) {
                    max_err = error;
                    worst_i = i;
                    worst_j = j;
                }
                max_err = max(max_err, error);
                avg_err += error;
#ifdef DEBUG                
                printf("%lf %lf\n", dist[i][j], dist2[i][j]);
#endif                
                ++total;
            }
        avg_err /= total;
        printf("%i %i  %lf %lf\n", worst_i, worst_j, dist[worst_i][worst_j], dist2[worst_i][worst_j]);
        printf("eps = %lf avg_err = %lf max_err = %lf\n", eps, avg_err, max_err);
        if (max_err < eps) break;
    }
    puts("passed!\n");
}

/*
 * Compute the entropy of a number x in (0,1)
 *
 */
double entropy(double x) {
    if (x == 0 || x == 1) return 0;
    return -(x * log2(x) + (1 - x) * log2(1 - x));
}

double cross_entropy(double p, double q) {
    if (q == 0 || q == 1) {
        if (p == 1) return 0;
        else return infinity;
    }
    return -(p * log2(q) + (1 - p) * log2(1 - q));
}

// prob doesn't need to be normalized.
int choose_from_distribution(const vector<double>& prob) {
    int n = prob.size();
    vector<double> sums(n + 1);
    for (int i = 0; i < n; ++i)
        sums[i + 1] = sums[i] + prob[i];
    double chosen = rand01() * sums.back();
    int index = upper_bound(sums.begin(), sums.end(), chosen) - sums.begin() - 1;
    return index;
}

int choose_from_distribution(const FenwickTree& tree) {
    return tree.upper_bound(rand01() * tree.partial_sum(tree.n - 1));
}

// Compute the median in linear time
template<typename T>
T median(vector<T> v) {
    int m = v.size() / 2;
    nth_element(v.begin(), v.begin() + m, v.end());
    return v[m];
//    return (v.size() & 1) ? v[m] : ((double)v[m] + *min_element(&v[0], &v[m])) / 2;
}

double jaccard(vector<int>& v, vector<int>& w) {
    if (v.empty() && w.empty()) return 0;

    vector<int> q;
    set_intersection(v.begin(), v.end(), w.begin(), w.end(), back_inserter(q));
    return (double)q.size() / (v.size() + w.size() - q.size());
}

#endif
