void graph_partition(int L) {
    int K = centers.size();
    ll pow2K = 1ll << K;
    vector<double> outside(pow2K), cost(pow2K);
    for (ll S = 0; S < pow2K; ++S) {
        int size = 0;
        for (int i = 0; i < K; ++i) if (S & (1ll << i)) size += blocks[i].size();
        outside[S] = pow((double)n / L - size, 2);

        for (int i = 0; i < K; ++i) if (S & (1ll << i))
            for (int j = 0; j < K; ++j) if (!(S & (1ll << j)))
                outside[S] += density[i][j] * blocks[i].size() * blocks[j].size();
    }

    vector<vector<ll>> best(L + 1);
    for (int l = 1; l <= L; ++l) best[l].resize(pow2K);

    cost = outside;
    for (int l = 2; l <= L; ++l)
        for (ll S = pow2K - 1; S >= 0; --S) {
            cost[S] = infinity;
            ll T = S;
            for (;;) {
                T = (T - 1) & S;
                if (cost[T] + outside[T] < cost[S]) {
                    cost[S] = cost[T] + outside[T];
                    best[l][S] = T;
                }
                if (T == 0) break;
            }
        }

    printf("total cost = %lf\n", cost[pow2K - 1]);
    ll S = pow2K - 1;
    for (int l = L; l > 1; --l) {
        ll T = best[l][S];
        int size = 0;
        for (int i = 0; i < K; ++i) if (T & (1ll << i)) size += blocks[i].size();
        printf("%lli (%i, %lf):", T, size, cost[T]);

        for (int i = 0; i < K; ++i)
            if (T & (1 << i))
                printf("%i ", i);
        printf("\n");
        S = S & ~T;
    }
}
