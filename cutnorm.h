#ifndef CUTNORM_H
#define CUTNORM_H

// #define DEBUG_CN

#include <cmath>
#include <vector>

#include <sdpa_call.h>
#include <f2c.h>
#include <fcntl.h>

#include "util.h"

extern "C" void dpotrf_(char *uplo,int *jb,double *A,int *lda,int *info);
extern "C" void spstrf_(char *uplo,int *jb,double *A,int *lda,int* piv, int* rank, double* tol, double* work, int *info);

void printMatrix(double* ele, int dim, char* printFormat, FILE* fpout)
{
  fprintf(fpout,"[\n");
  for (int i=0; i<dim; ++i) {
    fprintf(fpout,"[ ");
    for (int j=0; j<dim-1; ++j) {
      fprintf(fpout,printFormat,ele[i+dim*j]);
      fprintf(fpout," ");
    }
    fprintf(fpout,printFormat,ele[i+dim*(dim-1)]);
    fprintf(fpout,"]; \n");
  }
  fprintf(fpout,"]; \n");
}

double cut_value(const vector<int>& S, const vector<int>& T, const Matrix<double>& m) {
    double sum = 0;
    for (vertex v : S)
        for (vertex w : T)
            sum += m[v][w];
    return sum;
}

vector<int> round_gaussian(const Matrix<double>& m) {
    vector<int> ret;

    Point p = gaussian_vector(m[0].size());
    for (int i = 0; i < m.size(); ++i)
        ret.push_back(dot_product(m[i], p) >= 0 ? 1 : -1);
    return ret;
}

double infinity_to_1_norm(const vector<vector<double>>& m, vector<int>& sol) {
#ifdef DEBUG_CN
    print(m);
#endif
    int m_size = m.size();

    SDPA Problem1;
    Problem1.setDisplay(0);
    Problem1.setParameterType(SDPA::PARAMETER_DEFAULT);

    Problem1.inputConstraintNumber(2 * m_size);            // number of matrices (excluding constant one)
    Problem1.inputBlockNumber(1);                     // number of blocks
    Problem1.inputBlockSize(1, 2 * m_size);                // size of first block
    Problem1.inputBlockType(1, SDPA::SDP);

    Problem1.initializeUpperTriangleSpace();
    for (int i = 0; i < 2 * m_size; ++i)
        Problem1.inputCVec(i + 1, 1);

    for (int i = 0, t = 1; i < m_size; ++i)
        for (int j = 0; j < m_size; ++j, ++t)
            Problem1.inputElement(0, 1, i + 1, m_size + j + 1, m[i][j]);     // matrix, block, coordinates, value

    for (int i = 0; i < 2 * m_size; ++i)
        Problem1.inputElement(i + 1, 1, i + 1, i + 1, 1);

    Problem1.initializeUpperTriangle();
    Problem1.initializeSolve();

#ifdef DEBUG_CN
    // if necessary, dump input data and initial point
    Problem1.writeInputSparse((char*)"tmp.dat-s",(char*)"%+8.3e");
    Problem1.writeInitSparse((char*)"tmp.ini-s",(char*)"%+8.3e");
#endif

    Problem1.solve();

#ifdef DEBUG_CN
    double value = Problem1.getPrimalObj();
    printf("  value = %lf\n", value);
#endif

    double* solY = Problem1.getResultYMat(1);
//    printMatrix(solY, 2 * m_size, (char*)"%+8.3e", stdout);

    // Cholesky factorization
    int info, mdim = 2 * m_size;
    char UPLO = 'L';
//    double tol = 0; int piv, rank; double* work = (double*)malloc(2 * m_size * m_size * sizeof(double));
//    spstrf_(&UPLO, &m_size, solY, &m_size, &piv, &rank, &tol, work, &info);

    for (int i = 0; i < mdim; ++i) solY[i * mdim] += 1e-14;
    dpotrf_(&UPLO, &mdim, solY, &mdim, &info); // who the fuck named these functions?
    if (info) printf("Cholesky failed!!! error %i\n", info);
    for (int i = 0; i < mdim; ++i)
        for (int j = 0; j < i; ++j) solY[j * m_size + i] = 0;
//    printMatrix(solY, 2 * m_size, (char*)"%+8.3e", stdout);

    vector<Point> v(mdim, Point(mdim));
    for (int i = 0; i < mdim; ++i)
        for (int j = 0; j < mdim; ++j)
            v[i][j] = solY[i * mdim + j];

    double rounded_value = 0;
    for (int t = 0; t < 10; ++t) {
        auto cur_sol = round_gaussian(v);
        double val = 0;
        for (int i = 0; i < m_size; ++i)
            for (int j = 0; j < m_size; ++j)
                val += cur_sol[i] * m[i][j] * cur_sol[j + m_size];
        if (t == 0 || fabs(val) > fabs(rounded_value)) {
            sol = cur_sol;
            rounded_value = val;
        }
    }

#ifdef DEBUG_CN
    printf("  sol = "); for (int i = 0; i < mdim; ++i) printf("%i ", sol[i]); puts("");
    printf("  rounded value = %lf\n", rounded_value);
#endif
    return rounded_value;
}

vector<int> positive_neighbours(const vector<int>& S, int sgn, const Matrix<double>& m) {
    vector<int> T;
    for (vertex w = 0; w < m.size(); ++w) {
        double sum = 0;
        for (vertex v : S)
            sum += m[v][w];
        if (sum * sgn > 0)
            T.push_back(w);
    }
    return T;
}

vector<int> positive_neighbours(int v, int sgn, const Matrix<double>& m) {
    vector<int> S;
    S.push_back(v);
    return positive_neighbours(S, sgn, m);
}

// beware: assumes symmetric matrix (the second call to positive_neighbours should actually use the transpose)
void extend_cut(const Matrix<double>& m, vector<int>& S, vector<int>& T) {
    double val = cut_value(S, T, m);
#ifdef DEBUG_CN
    printf("\n  extending cut with value %lf\n", val); print(S); print(T); printf("\n");
#endif
    int sgn = val > 0 ? 1 : -1;
    T = positive_neighbours(S, sgn, m);
    S = positive_neighbours(T, sgn, m);
#ifdef DEBUG_CN
    val = cut_value(S, T, m);
    printf("\n  new cut with value %lf\n", val); print(S); print(T); printf("\n");
#endif
}

void cut_norm_frieze_kannan(const Matrix<double>& m, vector<int>& S) {
    S.clear();
    vertex v = rand() % m.size();
    int sgn = (rand() & 1) ? 1 : -1;
    S = positive_neighbours(v, sgn, m);
}

void densest_subgraph_aux(const Matrix<double>& m, vector<int>& bestS) {
    int m_size = m.size();
    vector<double> deg(m_size), neg_dev(m_size);
    for (int i = 0; i < m_size; ++i)
        for (int j = 0; j < m_size; ++j)
            deg[i] += m[i][j];

    vector<int> S;
    for (int i = 0; i < m_size; ++i)
        S.push_back(i);
    double bestval = -infinity, sum = cut_value(S, S, m);
    for (int it = 0; it < m_size; ++it) {
        // Compute density and update best value
        double val = sum / S.size();
        if (val > bestval) {
            bestval = val;
            bestS = S;
        }
#ifdef DEBUG_CN
        printf("  it = %i val = %lf bestval = %lf\n", it, val, bestval);
#endif

        // Pick the vertex with the smallest degree
        int v = min_element(&deg[0], &deg[m_size]) - &deg[0];

        // Remove v from S and update degrees
        deg[v] = infinity;
        int idx = lower_bound(S.begin(), S.end(), v) - S.begin();
        remove_element(S, idx);
        for (int i : S) {
            deg[i] -= m[i][v];
            sum -= m[i][v] + m[v][i];
        }
    }
}

double densest_subgraph(Matrix<double> m, vector<int>& S) {
    int m_size = m.size();
    vector<int> T;

    densest_subgraph_aux(m, S);
    for (int i = 0; i < m_size; ++i) 
    for (int j = 0; j < m_size; ++j)
            m[i][j] = -m[i][j];
    densest_subgraph_aux(m, T);

    double retS = cut_value(S, S, m),
           retT = -cut_value(T, T, m);
    if (fabs(retS) < fabs(retT)) {
        S = T;
        swap(retS, retT);
    }
    return retS;
}

void cut_norm_grothendieck(const Matrix<double>& m, vector<int>& S, vector<int>& T) {
    int m_size = m.size();

    Matrix<double> m2 = m;
    m2.resize(m_size + 1);
    for (int i = 0; i < m_size; ++i) m2[i].push_back(0);
    m2[m_size].resize(m_size + 1);
    for (int i = 0; i < m_size; ++i)
        for (int j = 0; j < m_size; ++j) {
            m2[i][m_size] -= m[i][j];
            m2[m_size][j] -= m[i][j];
        }

    vector<int> sol;
    infinity_to_1_norm(m2, sol);
    assert(sol.size() == 2 * (m_size + 1));

    if (sol[m_size] > 0) for (int i = 0; i <= m_size; ++i) sol[i] *= -1;
    if (sol[2 * m_size + 1] > 0) for (int i = m_size + 1; i <= 2 * m_size + 1; ++i) sol[i] *= -1;

    S.clear();
    for (int i = 0; i < m_size; ++i) if (sol[i] > 0) S.push_back(i);
    T.clear();
    for (int i = 0; i < m_size; ++i) if (sol[i + m_size + 1] > 0) T.push_back(i);

#ifdef DEBUG_CN
    printf("  new val = %lf\n", cut_value(S, T, m));
#endif
}

double cut_norm(const Matrix<double>& m, vector<int>& bestS, vector<int>& bestT) {
    double bestval = 0;
    for (int t = 0; t < 4; ++t) {
        vector<int> S, T;
        if (t == 0) cut_norm_frieze_kannan(m, S);
        else if (t == 1) densest_subgraph(m, S);
        else cut_norm_grothendieck(m, S, T);
        extend_cut(m, S, T);
#ifdef DEBUG_CN
        print(S); print(T);
#endif
        double val = cut_value(S, T, m);
        if (t == 0 || fabs(val) > fabs(bestval)) {
            bestval = val;
            bestS = S;
            bestT = T;
        }
#ifdef DEBUG_CN
        printf("  best so far = %lf\n", bestval);
#endif
    }

    return bestval; // We don't return the fabs() on purpose. Some other functions expect signed value.
}

double from_ST_to_S(vector<int>& S, const vector<int>& T, const Matrix<double> &m) {
    vector<vector<int>> possib(4);

    set_difference(S.begin(), S.end(), T.begin(), T.end(), back_inserter(possib[0]));
    set_difference(T.begin(), T.end(), S.begin(), S.end(), back_inserter(possib[1]));
    set_intersection(S.begin(), S.end(), T.begin(), T.end(), back_inserter(possib[2]));
    set_union(S.begin(), S.end(), T.begin(), T.end(), back_inserter(possib[3]));

    int best = 0;
    double bestw = -1, absbestw = -1;
    for (int i = 0; i < 4; ++i) {
        double w = cut_value(possib[i], possib[i], m);
#ifdef DEBUG_CN
        printf("  %i %i %lf\n", i, (int)possib[i].size(), w);
#endif
        if (fabs(w) > absbestw || (fabs(w) == absbestw && possib[i].size() < possib[best].size())) {
            bestw = w;
            absbestw = fabs(bestw);
            best = i;
        }
    }
#ifdef DEBUG_CN
    printf("  chosen %i\n", best);
#endif

    S = possib[best];
    return bestw;
}

double cut_norm(const Matrix<double>& m) {
    vector<int> S, T;
    double val = cut_norm(m, S, T);
    int sgn = val > 0 ? 1 : -1;
    T = positive_neighbours(S, sgn, m);
    S = positive_neighbours(T, sgn, m);
    return cut_value(S, T, m);
}

double cut_norm_single(const Matrix<double>& m, vector<int>& S) {
    vector<int> T;
    cut_norm(m, S, T);
    return from_ST_to_S(S, T, m);
}

double cut_norm_single(const Matrix<double>& m) {
    vector<int> S;
    return cut_norm_single(m, S);
}

#endif
