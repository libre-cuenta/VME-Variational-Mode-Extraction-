from vme import VME

def VME_K(K, f, Alpha, omega_int, fs, tau, tol):
    f1 = f[::]
    u = list()
    u_hat = list()
    omega = list()

    for i in range(K):
        u_s, u_hat_s, omega_s = VME(f1, Alpha, omega_int, fs, tau, tol)
        f1 = f1 - u_s
        
        u.append(u_s)
        u_hat.append(u_hat_s)
        omega.append(omega_s)
    
    res = f1
    return u, res, u_hat, omega
