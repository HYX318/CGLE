#ifndef CGLE_H
#define CGLE_H

#include <vector>
#include <tuple>

// 类型别名定义
using Vec1D = std::vector<double>;
using Vec2D = std::vector<Vec1D>;

class CGLE {
private:
    // 模型参数
    double omega0, l, beta, alpha;
    double epsilon, tau, k, hh, r1;
    int M, n, N, m, p;
    Vec2D u, v;  // 状态变量

    // 初始化参数
    void init_parameters();

public:
    // 构造函数
    CGLE(int space_steps = 20, int time_segments = 500, int segment_steps = 200);
    
    // 初始化状态
    void initialize_state();
    
    // 运行模拟
    void run_simulation();
    
    // 获取结果
    std::tuple<Vec2D, Vec2D, Vec1D, Vec1D> get_results() const;
    
    // 生成Poincare截面数据
    std::tuple<Vec1D, Vec1D> get_poincare_section() const;
};

#endif // CGLE_H