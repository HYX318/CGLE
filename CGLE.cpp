#include "CGLE.h"
#include <cmath>
#include <cassert>

#define PI 3.141592653

// 构造函数
CGLE::CGLE(int space_steps, int time_segments, int segment_steps) 
    : n(space_steps), M(time_segments), N(segment_steps) {
    init_parameters();
}

// 初始化参数
void CGLE::init_parameters() {
    omega0 = 1.6512;
    l = 2 * PI;
    beta = PI / 15;
    alpha = 1.0;

    // 计算派生参数
    double mu1 = -0.003;
    double mu2 = -0.005;
    double epc = 0.92363;
    double tauc = 1.59602;
    epsilon = epc + mu1;
    tau = tauc + mu2;

    // 网格参数
    hh = l * PI / (n - 1);  // 空间步长
    k = tau / N;            // 时间步长
    m = N + M * N + 1;      // 总时间步数
    p = N + 1;
    r1 = k / (pow(hh, 2));  // 扩散系数
}

// 初始化状态变量
void CGLE::initialize_state() {
    u.resize(n, Vec1D(m, 0.01));
    v.resize(n, Vec1D(m, 0.01));
}

// 运行CGLE模拟
void CGLE::run_simulation() {
    for (int j = p - 1; j < m - 1; ++j) {
        // 内部点计算
        for (int i = 1; i < n - 1; ++i) {
            assert(j - N >= 0 && "Time index out of bounds");
            
            // 计算u的下一时间步
            u[i][j+1] = u[i][j] + k * (
                r1 * (u[i+1][j] - 2 * u[i][j] + u[i-1][j]) +
                (1 - epsilon * cos(beta)) * u[i][j] -
                (omega0 + epsilon * sin(beta)) * v[i][j] +
                epsilon * cos(beta) * u[i][j-N] + 
                epsilon * sin(beta) * v[i][j-N] -
                (pow(u[i][j], 2) + pow(v[i][j], 2)) * u[i][j] +
                alpha * pow(v[i][j], 3) + 
                alpha * pow(u[i][j], 2) * v[i][j]
            );
            
            // 计算v的下一时间步
            v[i][j+1] = v[i][j] + k * (
                r1 * (v[i+1][j] - 2 * v[i][j] + v[i-1][j]) +
                (omega0 + epsilon * sin(beta)) * u[i][j] +
                (1 - epsilon * cos(beta)) * v[i][j] -
                epsilon * sin(beta) * u[i][j-N] + 
                epsilon * cos(beta) * v[i][j-N] -
                (pow(u[i][j], 2) + pow(v[i][j], 2)) * v[i][j] -
                alpha * pow(u[i][j], 3) + 
                alpha * pow(v[i][j], 2) * u[i][j]
            );
        }
        
        // 边界条件
        u[0][j+1] = u[1][j+1];
        u[n-1][j+1] = u[n-2][j+1];
        v[0][j+1] = v[1][j+1];
        v[n-1][j+1] = v[n-2][j+1];
    }
}

// 提取计算结果
std::tuple<Vec2D, Vec2D, Vec1D, Vec1D> CGLE::get_results() const {
    int size = M * N + 1;
    Vec2D u1(size, Vec1D(n));
    Vec2D v1(size, Vec1D(n));
    Vec1D re(size), im(size);
    
    // 提取数据
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < size; ++j) {
            u1[j][i] = u[i][j + N];
            v1[j][i] = v[i][j + N];
        }
    }
    
    // 提取实部和虚部
    for (int j = 0; j < size; ++j) {
        re[j] = u1[j][1];
        im[j] = v1[j][1];
    }
    
    return {u1, v1, re, im};
}

// 生成Poincare截面数据
std::tuple<Vec1D, Vec1D> CGLE::get_poincare_section() const {
    int t_length = tau * M / k + 1;
    int inte = 0;
    Vec1D uu, vv, uud;
    
    // 提取数据段
    for (int i = 2 * N + inte; i < t_length; ++i) {
        uu.push_back(u[1][i]);
        vv.push_back(v[1][i]);
    }
    
    for (int i = N + inte; i < t_length - N; ++i) {
        uud.push_back(u[1][i]);
    }
    
    // 计算截面
    Vec1D Poincarex, Poincarey;
    const double threshold = 0.0;
    for (int nn = 0; nn < (int)uud.size() - 1; ++nn) {
        if ((uud[nn] - threshold) * (uud[nn + 1] - threshold) < 0) {
            Poincarex.push_back(uu[nn]);
            Poincarey.push_back(vv[nn]);
        }
    }
    
    return {Poincarex, Poincarey};
}