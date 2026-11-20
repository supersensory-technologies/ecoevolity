// GPU spike: phycoeval's per-branch matrix-exponential step (cf_expmvCOMPLEX on
// QMatrix) in FP64 on the GPU, against ecoevolity's own CPU code as the oracle.
//
// For each allele count N, P pattern vectors go through exp(Q t) x with the
// 12-pole Caratheodory-Fejer approximation, exactly as MatrixExponentiator does:
// one structured complex solve (Q^T + offset I) x = v per pattern and pole, then
// w = Re(sum_i c_i x_i). Kernel A: one thread per (pattern, pole) solve.
// Kernel B: one thread per (pattern, element) combines the 12 poles in the
// CPU's order.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "matrix.hpp"  // ecoevolity: QMatrix, MatrixExponentiator (the oracle)

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    std::fprintf(stderr, "CUDA %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); std::exit(2); } } while (0)

struct C { double re, im; };

// The same operation order as ecoevolity's Complex.
__device__ __forceinline__ C cdiv(C n, C d) {
    double f = d.re * d.re + d.im * d.im;
    return C{(n.re * d.re + n.im * d.im) / f, (n.im * d.re - n.re * d.im) / f};
}
__device__ __forceinline__ C cmul(C m1, C m2) {
    return C{m2.re * m1.re - m2.im * m1.im, m2.re * m1.im + m2.im * m1.re};
}

constexpr int MAXN = 64;
__constant__ double c_zr[12], c_zi[12], c_cr[12], c_ci[12];

// (M_n + offset I) x = y for one tridiagonal central block (general u, v > 0 case).
__device__ void central_block(const C* y, C K, int n, double u, double v, C* x, C* d, C* e) {
    d[0] = K;
    e[0] = y[0];
    for (int r = 1; r <= n; ++r) {
        C m = cdiv(C{r * u, 0.0}, d[r - 1]);
        d[r].re = K.re + r * (v - u) - m.re * (n - r + 1.0) * v;
        d[r].im = K.im + -m.im * (n - r + 1.0) * v;
        C t = cmul(m, e[r - 1]);
        e[r] = C{y[r].re - t.re, y[r].im - t.im};
    }
    x[n] = cdiv(e[n], d[n]);
    for (int xr = n; xr >= 1; --xr) {
        int r = xr - 1;
        C t{e[r].re - (double)(n - r) * v * x[r + 1].re, e[r].im - (double)(n - r) * v * x[r + 1].im};
        x[r] = cdiv(t, d[r]);
    }
}

// One (pattern, pole): solve (Q^T + offset I) x = vc, with vc = w / t.
__global__ void solve_kernel(const double* w, int L, int N, double u, double v, double coal,
                             double t, int P, C* xout) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= P * 12) return;
    int p = id / 12, pole = id % 12;
    const double* wp = w + (size_t)p * L;
    C* x = xout + (size_t)id * L;
    C offset{-c_zr[pole] / t, -c_zi[pole] / t};
    C xn[MAXN + 1], yn[MAXN + 1], d[MAXN + 1], e[MAXN + 1];
    x[0] = C{0, 0};
    int xptr = L - 1 - N, yptr = L - 1 - N;
    for (int i = 0; i <= N; ++i) yn[i] = C{wp[yptr + i] / t, 0.0};
    C K{-coal * (double)N * (N - 1.0) / 2.0 - N * v + offset.re, offset.im};
    central_block(yn, K, N, u, v, xn, d, e);
    for (int i = 0; i <= N; ++i) x[xptr + i] = xn[i];
    for (int n = N - 1; n >= 1; --n) {
        xptr -= (n + 1);
        for (int r = 0; r <= n; ++r) {
            double f1 = coal * (double)(n - r) * (n + 1.0) / 2.0, f2 = coal * (double)r * (n + 1.0) / 2.0;
            yn[r] = C{f1 * xn[r].re + f2 * xn[r + 1].re, f1 * xn[r].im + f2 * xn[r + 1].im};
        }
        yptr -= (n + 1);
        for (int r = 0; r <= n; ++r) {
            double yr = (yptr + r) >= 1 ? wp[yptr + r] / t : 0.0;  // vc[0] is 0
            yn[r] = C{yr - yn[r].re, 0.0 - yn[r].im};
        }
        C Kn{-coal * (double)n * (n - 1.0) / 2.0 - n * v + offset.re, offset.im};
        central_block(yn, Kn, n, u, v, xn, d, e);
        for (int i = 0; i <= n; ++i) x[xptr + i] = xn[i];
    }
}

// w[j] = Re(sum_i c_i x_i[j]) in the CPU's pole order.
__global__ void combine_kernel(const C* xs, int L, int P, double* out) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= P * L) return;
    int p = id / L, j = id % L;
    if (j == 0) { out[id] = 0.0; return; }
    C acc{0, 0};
    for (int i = 0; i < 12; ++i) {
        C xi = xs[((size_t)p * 12 + i) * L + j];
        acc.re += xi.re * c_cr[i] - xi.im * c_ci[i];
        acc.im += xi.re * c_ci[i] + xi.im * c_cr[i];
    }
    out[id] = acc.re;
}

__global__ void empty_kernel() {}

int main(int argc, char** argv) {
    const int P = argc > 1 ? std::atoi(argv[1]) : 477;
    const int reps = 30;
    // Fabian's run: base frequency 0.5 (u = v = 1), theta = 2 * ploidy * N_e * mu
    // with N_e near 0.0005, branch lengths around the root height (~0.06).
    const double u = 1.0, v = 1.0, theta = 2 * 2 * 0.0005 * 1.0, t = 0.03;
    const double coal = 2.0 / theta;

    // The CF poles and residues, copied from ecoevolity's MatrixExponentiator
    // (consts.cpp is generated from matrix.hpp).
    extern const double ECO_ZR[12], ECO_ZI[12], ECO_CR[12], ECO_CI[12];
    CK(cudaMemcpyToSymbol(c_zr, ECO_ZR, sizeof(double) * 12));
    CK(cudaMemcpyToSymbol(c_zi, ECO_ZI, sizeof(double) * 12));
    CK(cudaMemcpyToSymbol(c_cr, ECO_CR, sizeof(double) * 12));
    CK(cudaMemcpyToSymbol(c_ci, ECO_CI, sizeof(double) * 12));

    // Launch latency: one empty kernel, synchronised.
    empty_kernel<<<1, 1>>>(); CK(cudaDeviceSynchronize());
    auto l0 = std::chrono::steady_clock::now();
    for (int i = 0; i < 1000; ++i) { empty_kernel<<<1, 1>>>(); CK(cudaDeviceSynchronize()); }
    double launch_us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - l0).count() / 1000;
    std::printf("launch+sync round trip: %.1f us\n", launch_us);
    std::printf("N\tlen\tcpu_ms\tgpu_kernel_ms\tgpu_e2e_ms\tmax_rel_err\tmax_abs_err\n");

    MatrixExponentiator me;
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    for (int N = 1; N <= 45; ++N) {
        QMatrix q(N, u, v, theta);
        const int L = (int)q.get_nrows() + 1;  // xcol has a leading 0
        std::vector<double> w((size_t)P * L, 0.0);
        for (int p = 0; p < P; ++p) {
            double s = 0;
            for (int j = 1; j < L; ++j) { w[(size_t)p * L + j] = unif(rng); s += w[(size_t)p * L + j]; }
            for (int j = 1; j < L; ++j) w[(size_t)p * L + j] /= s;
        }
        // CPU oracle: ecoevolity's own code, one thread.
        std::vector<double> cpu((size_t)P * L);
        auto c0 = std::chrono::steady_clock::now();
        int cpu_reps = N < 20 ? 5 : 2;
        for (int k = 0; k < cpu_reps; ++k)
            for (int p = 0; p < P; ++p) {
                std::vector<double> xcol(w.begin() + (size_t)p * L, w.begin() + (size_t)(p + 1) * L);
                QMatrix qp(N, u, v, theta);
                std::vector<double> y = me.cf_expmvCOMPLEX(t, qp, xcol);
                std::copy(y.begin(), y.end(), cpu.begin() + (size_t)p * L);
            }
        double cpu_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - c0).count() / cpu_reps;

        double *dw, *dout; C* dx;
        CK(cudaMalloc(&dw, sizeof(double) * P * L));
        CK(cudaMalloc(&dout, sizeof(double) * P * L));
        CK(cudaMalloc(&dx, sizeof(C) * (size_t)P * 12 * L));
        std::vector<double> gpu((size_t)P * L);
        cudaEvent_t e0, e1; CK(cudaEventCreate(&e0)); CK(cudaEventCreate(&e1));
        int tb = 128;
        auto run = [&](bool timed_kernel) {
            CK(cudaMemcpy(dw, w.data(), sizeof(double) * P * L, cudaMemcpyHostToDevice));
            if (timed_kernel) CK(cudaEventRecord(e0));
            solve_kernel<<<(P * 12 + tb - 1) / tb, tb>>>(dw, L, N, u, v, coal, t, P, dx);
            combine_kernel<<<(P * L + 255) / 256, 256>>>(dx, L, P, dout);
            if (timed_kernel) CK(cudaEventRecord(e1));
            CK(cudaMemcpy(gpu.data(), dout, sizeof(double) * P * L, cudaMemcpyDeviceToHost));
        };
        run(false);  // warm
        float kms_total = 0;
        auto g0 = std::chrono::steady_clock::now();
        for (int k = 0; k < reps; ++k) {
            run(true);
            float kms; CK(cudaEventElapsedTime(&kms, e0, e1)); kms_total += kms;
        }
        double e2e_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - g0).count() / reps;
        CK(cudaGetLastError());

        double max_rel = 0, max_abs = 0;
        for (int p = 0; p < P; ++p)
            for (int j = 1; j < L; ++j) {
                double a = cpu[(size_t)p * L + j], b = gpu[(size_t)p * L + j];
                double ad = std::fabs(a - b);
                max_abs = std::max(max_abs, ad);
                if (std::fabs(a) > 1e-12) max_rel = std::max(max_rel, ad / std::fabs(a));
            }
        std::printf("%d\t%d\t%.3f\t%.3f\t%.3f\t%.2e\t%.2e\n", N, L - 1, cpu_ms, kms_total / reps, e2e_ms, max_rel, max_abs);
        cudaFree(dw); cudaFree(dout); cudaFree(dx);
        cudaEventDestroy(e0); cudaEventDestroy(e1);
    }
    return 0;
}
