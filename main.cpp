#include "CGLE.h"
#include "matplotlibcpp.h"
#include <iostream>

namespace plt = matplotlibcpp;

// 可视化结果
void plot_results(const Vec1D& re, const Vec1D& im,
                 const Vec1D& poincare_x, const Vec1D& poincare_y) {
    // 绘制实部vs虚部
    plt::figure();
    plt::plot(re, im);
    plt::title("Phase Portrait (Real vs Imaginary)");
    plt::xlabel("Real Part");
    plt::ylabel("Imaginary Part");
    plt::pause(1);
    
    // 绘制Poincare截面
    // 绘制Poincare截面
    if (!poincare_x.empty() && !poincare_y.empty()) {
        plt::figure();
        int half = poincare_x.size() / 2;
        
        // 第一段散点（红色）
        Vec1D x1(poincare_x.begin(), poincare_x.begin() + half);
        Vec1D y1(poincare_y.begin(), poincare_y.begin() + half);
        plt::scatter(x1, y1, 5.0, {{"color", "red"}});
        
        // 第二段散点（蓝色）
        Vec1D x2(poincare_x.begin() + half, poincare_x.end());
        Vec1D y2(poincare_y.begin() + half, poincare_y.end());
        plt::scatter(x2, y2, 5.0, {{"color", "blue"}});
        
        plt::title("Poincare Section");
        plt::xlabel("u");
        plt::ylabel("v");
        plt::pause(1);
    }    
    plt::show();
}

int main() {
    try {
        // 初始化模型（空间步数、时间段数、每段步数）
        CGLE model(20, 80, 200);
        
        // 初始化状态并运行模拟
        model.initialize_state();
        model.run_simulation();
        std::cout << "Simulation completed successfully" << std::endl;
        
        // 获取结果并可视化
        auto [u1, v1, re, im] = model.get_results();
        auto [poincare_x, poincare_y] = model.get_poincare_section();
        plot_results(re, im, poincare_x, poincare_y);
    } 
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}