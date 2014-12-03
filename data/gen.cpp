// g++-4.8 -std=c++11 -Wall -Wextra -Wno-sign-compare -Wshadow -Werror -pedantic -o gen gen.cpp 
#include <cstdio>
#include <ctime>
#include <vector>
using namespace std;

#include "../util.h"

const int maxn = 2000;
long n, num_edges;
bool matrix[maxn][maxn];

void create_graph(long N, long k) {
    n = N;
    vector<vector<double>> density;
    density.reserve(k);
    for (int i = 0; i < k; ++i) {
	 vector<double> tmp(k, 0.0);
	 density[i] = tmp;
    }
    for (int i = 0; i < k; ++i) {
        for (int j = i; j < k; ++j) {
            density[i][j] = rand01();
            density[j][i] = density[i][j];
        }
    }

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            double d = density[i % k][j % k];
            matrix[i][j] = (rand01() <= d);

            if (i == j) matrix[i][j] = false;
            matrix[j][i] = matrix[i][j];
        }

    vector<int> perm(n, 0);
    for (int i = 0; i < n; ++i) perm[i] = i;
    random_shuffle(&perm[0], &perm[n]);

    double index = 0;
    for (int i = 0; i < k; ++i)
        for (int j = i + 1; j < k; ++j)
            index += pow(density[i][j] / k, 2);
    num_edges = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (matrix[i][j]) {
                printf("%i %i\n", perm[i], perm[j]);
                ++num_edges;
            }
    fprintf(stderr, "index = %lf\n", index);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s nodes clusters\n", argv[0]);
	return 1;
    }
    srand(time(NULL));
    create_graph(strtol(argv[1], NULL, 10), strtol(argv[2], NULL, 10));
    return 0;
}

