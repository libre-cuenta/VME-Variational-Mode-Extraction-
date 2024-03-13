#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <array>
#include <string_view>

using namespace std;

std::vector<std::complex<double>> fftshift(std::vector<std::complex<double>> f_hat) {
    int T = f_hat.size();
    std::vector<std::complex<double>> f_hat_shifted(T);
    int ltemp = T / 2;
    for (int i = 0; i < ltemp; i++) {
        f_hat_shifted[i] = f_hat[i + ltemp];
        f_hat_shifted[i + ltemp] = f_hat[i];
    }
    return f_hat_shifted;
}

template<typename T> std::vector<std::complex<double>> fft(T f) {
    int N = f.size();
    std::vector<std::complex<double>> f_hat(N);
    for (int k = 0; k < N; k++) {
        std::complex<double> sum = 0;
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<double>(0, -2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

template<typename T> std::vector<std::complex<double>> ifft(T f) {
    int N = f.size();
    std::vector<std::complex<double>> f_hat(N);
    for (int k = 0; k < N; k++) {
        std::complex<double> sum = 0;
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<double>(0, 2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
    return result;
}

std::vector<double> arange(double start, double end, double step) {
    std::vector<double> result;
    for (double i = start; i < end; i += step) {
        result.push_back(i);
    }
    return result;
}

std::vector<double> dot_product(std::vector<double> a, std::vector<double> b) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] * b[i];
    }
    return result;
}

std::vector<double> abs(std::vector<double> a) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = std::abs(a[i]);
    }
    return result;
}

double sum(std::vector<double> a) {
    double result = 0;
    for (int i = 0; i < a.size(); i++) {
        result += a[i];
    }
    return result;
}

std::vector<double> sort(std::vector<double> a) {
    std::sort(a.begin(), a.end());
    return a;
}

std::vector<double> exp(std::vector<double> a) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = std::exp(a[i]);
    }
    return result;
}

std::vector<double> random(int num) {
    std::vector<double> result(num);
    for (int i = 0; i < num; i++) {
        result[i] = (double)rand() / RAND_MAX;
    }
    return result;
}


vector<vector<complex<double>>> VME(vector<double> signal, double Alpha, double omega_int, double fs, double tau, double tol) {
    int save_T = signal.size();
    int T = save_T;
    vector<double> f_mir(2*T);
    for (int i = 0; i < T/2; i++) {
        f_mir[i] = signal[i];
    }
    for (int i = T/2; i < 3*T/2; i++) {
        f_mir[i] = signal[T - i - 1];
    }
    for (int i = 3*T/2; i < 2*T; i++) {
        f_mir[i] = signal[i - T/2];
    }
    vector<double> f(T);
    for (int i = 0; i < T; i++) {
        f[i] = f_mir[T - i - 1];
    }
    T = f.size();
    vector<double> t(T);

    for (int i = 0; i < T; i++) {
        t[i] = i / (double)T;
    }

    double eps = numeric_limits<double>::epsilon();
    double udiff = tol + eps;
    vector<double> omega_axis(T);
    for (int i = 0; i < T; i++) {
        omega_axis[i] = i - 0.5 - 1 / (double)T;
    }

    vector<complex<double>> f_hat(T);
    for (int i = 0; i < T; i++) {
        f_hat[i] = complex<double>(0, 0);
    }
    f_hat = fftshift(fft<vector<double>>(f));

    vector<complex<double>> f_hat_onesided = f_hat;
    for (int i = 0; i < T/2; i++) {
        f_hat_onesided[i] = complex<double>(0, 0);
    }

    int N = 500;

    vector<complex<double>> omega_d(N);
    for (int i = 0; i < N; i++) {
        omega_d[i] = complex<double>(0, 0);
    }
    omega_d[0] = complex<double>(omega_int / fs, 0);

    vector<vector<complex<double>>> lamb(N, vector<complex<double>>(T));
    vector<vector<complex<double>>> u_hat_d(N, vector<complex<double>>(T));
    int n = 0;

    auto count_uhat = [&]() {
        vector<complex<double>> result(T);
        for (int i = 0; i < T; i++) {
            result[i] = (f_hat_onesided[i] + (u_hat_d[n][i] * pow(Alpha, 2) * pow(omega_axis[i] - omega_d[n], 4)) + (lamb[n][i] / 2.0)) / (((pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n], 4)) * ((2 * pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n], 4)));
        }
        return result;
    };

    auto count_omegad = [&]() {
        complex<double> A = 0;
        complex<double> B = 0;
        for (int i = T/2; i < T; i++) {
            A += omega_axis[i] * pow(abs(u_hat_d[n+1][i]), 2);
            B += pow(abs(u_hat_d[n+1][i]), 2);
        }
        return A / B;
    };

    auto count_lamb = [&]() {
        vector<complex<double>> result(T);
        for (int i = 0; i < T; i++) {
            result[i] = lamb[n][i] + (tau * (f_hat_onesided[i] - (u_hat_d[n+1][i] + (((pow(Alpha, 2) * pow(omega_axis[i] - omega_d[n+1], 4)) * (f_hat_onesided[i] - u_hat_d[n+1][i])) / (2 * pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n+1], 4)))));
        }
        return result;
    };

    auto count_udiff = [&](double udiff) {
        complex<double> result = udiff;
        for (int i = 0; i < T; i++) {
            result += ((u_hat_d[n][i] - u_hat_d[n-1][i]) * conj(u_hat_d[n][i] - u_hat_d[n-1][i]) *= (1 / T));
        }
        return result;
    };

    while (udiff > tol && n < N-1) {
        u_hat_d[n+1] = count_uhat();
        omega_d[n+1] = count_omegad();
        lamb[n+1] = count_lamb();
        n = n + 1;
        udiff = numeric_limits<double>::epsilon();
        complex<double> udiff_buf = count_udiff(udiff);
        udiff = abs(udiff_buf);
    }

    N = min(N, n) - 1;
    vector<complex<double>> omega(N);
    for (int i = 0; i < N; i++) {
        omega[i] = omega_d[i];
    }
    vector<complex<double>> u_hatd(T);
    for (int i = 0; i < T; i++) {
        u_hatd[i] = complex<double>(0, 0);
    }

    for (int i = T/2; i < T; i++) {
        u_hatd[i] = u_hat_d[N][i];
    }
    for (int i = T/2; i > 0; i--) {
        u_hatd[T - i] = conj(u_hat_d[N][i]);
    }
    u_hatd[0] = conj(u_hatd[T-1]);

    vector<complex<double>> u_d(T);
    for (int i = 0; i < T; i++) {
        u_d[i] = 0;
    }
    u_d = fftshift(ifft<vector<complex<double>>>(u_hatd));

    vector<complex<double>> u_d_final(T/2);
    for (int i = T/4; i < 3*T/4; i++) {
        u_d_final[i - T/4] = u_d[i];
    }
    vector<complex<double>> u_hatd_final = fftshift(fft<vector<complex<double>>>(u_d));

    vector<vector<complex<double>>> output = { u_d_final , u_hatd_final , omega };
    return output;
}
