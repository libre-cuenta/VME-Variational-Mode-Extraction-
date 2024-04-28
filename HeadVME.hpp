#include "header.hpp"
//F<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb

template<typename SignalDT>
class F;           //���������������� ��������� ���

class Ud;          //������������ (������������� �������)
class OmegaAxis;   //���������� ��������� �������
class NumbIter;    //����������� ���������� ��������
class Epsilon;     //������� ��������� �������� (������������ �������� ������������� �������� ������(����� ��������) � ������� ��������)
class Omega;       //�����-�������� ��� ������������
class Lamb;        //������-��������
class VMEDecomposition; //����� ������������


template<typename SignalDT>
class F
{
private:
    vector<complex<double>> f_hat_onesided;
    int T;
public:
    //������������ � ����������� �������, ��������� �������� ������� ������� (fs) � ����� ������� �
    F(vector<SignalDT> signal);

    vector<complex<double>> operator()() const;
    complex<double> operator[](int i) const;
    F<SignalDT>& operator=(vector<complex<double>> value);
    vector<complex<double>> getF() const;
    size_t size() const;
    int getT();
};

class Epsilon
{
private:
    complex<double> udiff;
    complex<double> count_udiff(complex<double> udiff, int T, int n, Ud& u_hat_d);
public:
    //double eps = numeric_limits<double>::epsilon();
    //double udiff = tol + eps;
    //������������� ������� ���������
    Epsilon(double tol);

    complex<double> getUdiff() const;
    void setUdiff(complex<double> value);
    //���������� ������� ��������� ����� ��������� ��������
    void update_udiff(complex<double> udiff, int T, int n, Ud& u_hat_d);

};
class OmegaAxis
{
private:
    vector<double> f;
public:
    //��������� ���������� ������� ������
    OmegaAxis(int T);

    vector<double> operator()() const;
    double operator[](int i) const;
    OmegaAxis& operator=(vector<double> value);
    size_t size() const;
};
class NumbIter
{
private:
    //int N = 500;
    int N = 500;
public:
    NumbIter();

    int operator()() const;
    NumbIter& operator += (const int& counter);
    void SetMinToN(int n);

};
class Lamb
{
private:
    vector<vector<complex<double>>> lamb;
    vector<complex<double>> count_lamb(int T, int n, double tau, double Alpha, vector<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d);
public:
    //������������� ������
    Lamb(int N = 500, int T);
    //���������� ������
    void update_lamb(int T, int n, double tau, double Alpha, vector<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d);
    vector<vector<complex<double>>> operator()() const;
    vector<complex<double>>& operator[](int i);
    size_t size_1() const;
    size_t size_2() const;
};
class Omega
{
private:
    vector<complex<double>> omega_d;//����������� ������� ����� ������ ��������
    bool check_finite = true;
    complex<double> count_omegad(int T, int n, Ud& u_hat_d, OmegaAxis& omega_axis);
public:
    //������������� �������� ����� ��� ������ ����
    Omega(int N = 500, double omega_int, double fs);
    //�������� ������������
    bool isFinite();
    vector<complex<double>> operator()() const;
    complex<double>& operator[](int i);
    //���������� �������� �����
    void update_omega(int T, int n, Ud& u_hat_d, OmegaAxis& omega_axis);
    //��������� ���������� �������� �����
    void minOmega(int n);
};
class Ud
{
private:
    vector<vector<complex<double>>> u_hat_d; //�������� ����� ����� ������ ��������
    bool check_finite = true;
    vector<complex<double>> u_hatd_final; //���� �� ��������� �������
    vector<complex<double>> u_d_final; //���� �� ��������� �������
    vector<complex<double>> count_uhat(int T, int n, double Alpha, vector<complex<double>>& f_hat_onesided, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb);
public:
    Ud(int N = 500, int T);
    //���������� ����
    void update_u_hat(int T, int n, double Alpha, vector<complex<double>>& f_hat_onesided, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb);
    //�������� ������������
    bool isFinite();
    vector<vector<complex<double>>> operator()() const;
    vector<complex<double>>& operator[](int i);
    size_t size_1() const;
    size_t size_2() const;
    //��������� �������� ������������
    void ud(int T, int N);

    vector<complex<double>> getU_D() const;
    vector<complex<double>> getU_hatD() const;
};
class VMEDecomposition
{
private:
    vector<complex<double>> u_hatd; //���� �� ��������� �������
    vector<complex<double>> u_d; //���� �� ��������� �������
    vector<complex<double>> omega_d; //����������� �������
public:
    VMEDecomposition(vector<complex<double>> u_hatd, vector<complex<double>> u_d, vector<complex<double>> omega_d);
    vector<complex<double>> GetUhatD() const;
    vector<complex<double>> GetUD() const;
    vector<complex<double>> GetOmegaD() const;
};

template<typename SignalDT> 
VMEDecomposition VME(vector<SignalDT> signal, double Alpha, double omega_int, double fs, double tau, double tol);