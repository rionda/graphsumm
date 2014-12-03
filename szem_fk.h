#ifndef SZEMFK_H
#define SZEMFK_H
#include "cutnorm.h"
void szemeredi_frieze_kannan() {
    Matrix residual = full_adjacency_matrix();
    for (vertex v = 0; v < n; ++v)
        for (vertex w = 0; w < n; ++w)
            residual[v][w] -= 2.0 * num_edges / (n * n);
//    auto& m = residual; n = 2; m.resize(n, vector<double>(n)); m[0][0] = 0.5; m[0][1] = m[1][0] = -1.5; m[1][1] = 3.5;

    for (int t = 0; t < 100; ++t) {
        vector<int> S;
        printf("\n***frob = %lf\n", frobenius_squared(residual));

        double r = densest_subgraph(residual, S), dens1 = r / S.size();
        printf("***densest subgraph cut = %lf; size = %i; dens1 = %lf; dec_frob = %lf\n", r, (int)S.size(), dens1, dens1 * dens1);

        double norm = cut_norm_single(residual, S);
        printf("***cut norm single >= |%lf|; size = %i;", norm, (int)S.size());

        double dens = cut_value(S, S, residual) / S.size() / S.size();
        printf(" dens1 = %lf dens2 = %lf dec_frob = %lf\n", dens * S.size(), dens, pow(dens * S.size(), 2));
        for (int v : S)
            for (int w : S)
                residual[v][w] -= dens;
    }
}

#endif

