#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
using namespace std;

// Data structure to keep partial sums in O(log n) per operation
struct FenwickTree {
    int n;
    vector<double> tree;

    FenwickTree(int _n) : n(_n), tree(n + 1) {}

    // idx in [-1, n]
    double partial_sum(int idx) const {
        ++idx;
        double sum = 0;
        while (idx > 0){
            sum += tree[idx];
            idx -= (idx & -idx);
        }
        return sum;
    }

    double value(int idx) const {
        ++idx;
        double sum = tree[idx];
        if (idx > 0) {
            int z = idx - (idx & -idx);
            idx--;
            while (idx != z){
                sum -= tree[idx];
                // substract tree frequency between y and "the same path"
                idx -= (idx & -idx);
            }
        }
        return sum;
    }

    void update(int idx, double delta) {
        ++idx;
        while (idx <= n){
            tree[idx] += delta;
            idx += (idx & -idx);
        }
    }

    void set_value(int idx, double new_val) {
        update(idx, new_val - value(idx));
    }

    double upper_bound(double val) const {
        int a = -1, b = n - 1;
        while (b - a != 1) {
            int m = (a + b) / 2;
            if (partial_sum(m) > val)
                b = m;
            else
                a = m;
        }
        return b;
    }

    void check() const {
        double acum = 0;
        assert(upper_bound(value(0) - 1) == 0);
        for (int i = 0; i < n; ++i) {
            acum += value(i);
//            printf("%lf (%lf) ", value(i), partial_sum(i));
            assert(fabs(acum - partial_sum(i)) < 1e-6);

            double m = (partial_sum(i) + partial_sum(i - 1)) / 2;
            if (value(i) != 0)
                assert(upper_bound(m) == i);
        }
        puts("\nChecked!\n");
    }
};
