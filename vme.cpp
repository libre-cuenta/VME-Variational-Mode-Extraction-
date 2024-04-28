#include "HeadVME.hpp"

template<typename SignalDT>
VMEDecomposition VME(vector<SignalDT> signal, double Alpha, double omega_int, double fs, double tau, double tol) {
    F f_hat_onesided<SignalDT>(signal);
    int T = f_hat_onesided.getT();
    Epsilon eps(tol);
    OmegaAxis omega_axis(T);
    NumbIter N;
    
    Omega omega_d(N(), omega_int, fs);
    Lamb lamb(N(), T);
    Ud u_hat_d(N(), T);

    int n = 0;
    while (eps.getUdiff() > tol && n < N() - 1)
    {
        //u_hat_d[n+1] = count_uhat();
        u_hat_d.update_u_hat(T, n, Alpha, f_hat_onesided.getF(), omega_axis, omega_d, lamb);
        //omega_d[n+1] = count_omegad();
        omega_d.update_omega(T, n, u_hat_d, omega_axis);
        if (!u_hat_d.isFinite() || !omega_d.isFinite())
        {
            break;
        }

        //lamb[n+1] = count_lamb();
        lamb.update_lamb(T, n, tau, Alpha, f_hat_onesided.getF(), u_hat_d, omega_axis, omega_d);
        n+=1;

        eps.update_udiff(udiff, T, n, u_hat_d);
    }

    N.SetMinToN(n);
    omega_d.minOmega(N());
    u_hat_d.ud(T, N());

    return new VMEDecomposition(u_hat_d.getU_hatD(), u_hat_d.getU_D(), omega_d());
}
