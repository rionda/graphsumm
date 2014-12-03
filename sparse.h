struct Entry {
    vertex vert;
    double val;
    Entry(vertex a = 0, double b = 0) : vert(a), val(b) {}

    bool operator<(const Entry& e) const {
        if (vert != e.vert) return vert < e.vert;
        return val < e.val;
    }
};

vector<vector<Entry>> matrix;
vector<double> modp;

inline Point get_point(vertex v) {
    Point p(dimension);
    for (Entry e : matrix[v])
        p[e.vert] = e.val;
    return p;
}

/*
inline double sparse_sum1(const Point& p) {
    double ret = 0;
    for (Entry e : p)
        ret += fabs(e.val);
    return ret;
}

inline double sparse_sum2(const Point& p) {
    double ret = 0;
    for (Entry e : p)
        ret += e.val * e.val;
    return ret;
}

inline double sparse_sump(const Point& p) {
    return p_exponent == 2 ? sparse_sum2(p) : sparse_sum1(p);
}
*/

inline double sparse_dot_product(int v, const Point& p) {
    double ret = 0;
    for (Entry e : matrix[v]) 
        if (p[e.vert] != 0)
            ret += e.val * p[e.vert];
    return ret;
}

inline double sparse_dot_product(int v, int w) {
    auto a = matrix[v].begin(), b = matrix[w].begin();
    double ret = 0;
    while (a < matrix[v].end() && b < matrix[w].end()) {
        if (a->vert < b->vert) ++a;
        else if (b->vert < a->vert) ++b;
        else {
            ret += a->val * b->val;
            ++a; ++b;
        }
    }
    return ret;
}

inline double quick_l2_2_distance(vertex v, const Point& p, double modpp) {
    return max(0, modp[v] + modpp - 2 * sparse_dot_product(v, p));
}

inline double quick_l2_2_distance(vertex v, vertex w) {
    return max(0, modp[v] + modp[w] - 2 * sparse_dot_product(v, w));
}



inline double sparse_l1_correction(int v, const Point& p) {
    double ret = 0;
    for (Entry e : matrix[v]) 
        if (p[e.vert] != 0)
            ret += fabs(e.val - p[e.vert]) - (fabs(e.val) + fabs(p[e.vert]));
    return ret;
}

inline double sparse_l1_correction(int v, int w) {
    auto a = matrix[v].begin(), b = matrix[w].begin();
    double ret = 0;
    while (a < matrix[v].end() && b < matrix[w].end()) {
        if (a->vert < b->vert) ++a;
        else if (b->vert < a->vert) ++b;
        else {
            ret += fabs(a->val - b->val) - (fabs(a->val) + fabs(b->val));
            ++a; ++b;
        }
    }
    return ret;
}


inline double quick_l1_distance(vertex v, const Point& p, double modpp) {
    return modp[v] + modpp + sparse_l1_correction(v, p);
}

inline double quick_l1_distance(vertex v, vertex w) {
    return modp[v] + modp[w] + sparse_l1_correction(v, w);
}








inline double quick_distance_p(vertex v, const Point& p, double modpp) {
    return p_exponent == 2 ? quick_l2_2_distance(v, p, modpp) : quick_l1_distance(v, p, modpp);
}

inline double quick_distance(vertex v, const Point& p, double modpp) {
    return p_exponent == 2 ? sqrt(quick_l2_2_distance(v, p, modpp)) : quick_l1_distance(v, p, modpp);
}

inline double quick_distance_p(vertex v, vertex w) {
    return p_exponent == 2 ? quick_l2_2_distance(v, w) : quick_l1_distance(v, w);
}

inline double quick_distance(vertex v, vertex w) {
    return p_exponent == 2 ? sqrt(quick_l2_2_distance(v, w)) : quick_l1_distance(v, w);
}

inline double pref_distance_to_p(const Point& v, const Point& w) {
    return p_exponent == 2 ? l2_2_distance(v, w) : l1_distance(v, w);
}

inline double pref_distance(const Point& v, const Point& w) {
    return p_exponent == 2 ? l2_distance(v, w) : l1_distance(v, w);
}



struct SparsePoint {
    vector<Entry> vec;
    mutable double _mod22 = -1;

    double operator[](int i) {
        auto it = lower_bound(vec.begin(), vec.end(), Entry(i, -infinity));
        if (it != vec.end() && it->vert == i)
            return it->val;
        else
            return 0;
    }

    Point full_point(int x) const {
        Point ret(x);
        for (Entry e : vec)
            ret[e.vert] = e.val;
        return ret;
    }

    double mod22() const {
        if (_mod22 < 0) {
            _mod22 = 0;
            for (Entry a : vec)
                _mod22 += a.val * a.val;
        }
        return _mod22;
    }
};

SparsePoint make_sparse(const Point& p) {
    SparsePoint ret;
    for (int i = 0; i < p.size(); ++i)
        if (p[i] != 0)
            ret.vec.push_back(Entry(i, p[i]));
    return ret;
}

inline int find_next(const vector<Entry>& v, int start, int x) {
    int end = v.size();
    if (v.back() < x) return end;

    while (end - start != 1) {
        int m = (start + end) / 2;
        if (v[m].vert <= x)
            start = m;
        else
            end = m;
    }
    return start;
}

double dot_product(const SparsePoint& p, const SparsePoint& q) {
    /*
    const vector<Entry>& v = p.vec, w = q.vec;
    if (v.empty() || w.empty()) return 0;

    double ret = 0;
    for (int i = 0, j = 0; i < v.size() && j < w.size();) {
        if (v[i].vert == w[j].vert) {
            ret += v[i].val * w[j].val;
            ++i;
            ++j;
        } else if (v[i].vert < w[j].vert) {
            i = find_next(v, i, w[j].vert);
        } else {
            j = find_next(w, j, v[i].vert);
        }
    }
    return ret;
    */


    double ret = 0;
    auto a = p.vec.begin(), b = q.vec.begin();
    while (a != p.vec.end() && b != q.vec.end()) {
        if (a->vert < b->vert) ++a;
        else if (a->vert > b->vert) ++b;
        else {
            ret += a->val * b->val;
            ++a; ++b;
        }
    }
    return ret;
}

inline double dist22(const SparsePoint& p, const SparsePoint& q) {
    double d = dot_product(p, q);
    d = max(0, p.mod22() + q.mod22() - 2 * d);
    return d;
}

inline double euclidean_dist(const SparsePoint& p, const SparsePoint& q) {
    return sqrt(dist22(p, q));
}


inline SparsePoint get_sparse_point(vertex v) {
    SparsePoint p;
    p.vec = matrix[v];
    return p;
}
