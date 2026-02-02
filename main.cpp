#include "CGLE.h"
#include "matplotlibcpp.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>

namespace plt = matplotlibcpp;

// 可视化结果
void plot_results(const Vec1D& re, const Vec1D& im,
                 const Vec1D& poincare_x, const Vec1D& poincare_y) {
    plt::figure();
    plt::plot(re, im);
    plt::title("Phase Portrait (Real vs Imaginary)");
    plt::xlabel("Real Part");
    plt::ylabel("Imaginary Part");
    plt::pause(1);
    
    if (!poincare_x.empty() && !poincare_y.empty()) {
        plt::figure();
        int half = poincare_x.size() / 2;
        Vec1D x1(poincare_x.begin(), poincare_x.begin() + half);
        Vec1D y1(poincare_y.begin(), poincare_y.begin() + half);
        plt::scatter(x1, y1, 5.0, {{"color", "red"}});
        
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

// 解析命名命令行参数
bool parse_command_line(int argc, char* argv[], int& space_steps, int& time_segments, int& segment_steps) {
    // 默认值（与原代码一致）
    space_steps = 20;
    time_segments = 50000;
    segment_steps = 200;

    // 遍历参数（跳过argv[0]：程序名）
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        // 处理--time（时间段数M）
        if (arg == "--time" && i + 1 < argc) {
            try {
                time_segments = std::stoi(argv[++i]);
                if (time_segments <= 0) {
                    std::cerr << "错误：--time 参数必须大于0！" << std::endl;
                    return false;
                }
            } catch (const std::exception& e) {
                std::cerr << "错误：--time 后必须跟正整数！" << std::endl;
                return false;
            }
        }
        // 处理--space（空间步数）
        else if (arg == "--space" && i + 1 < argc) {
            try {
                space_steps = std::stoi(argv[++i]);
                if (space_steps <= 0) {
                    std::cerr << "错误：--space 参数必须大于0！" << std::endl;
                    return false;
                }
            } catch (const std::exception& e) {
                std::cerr << "错误：--space 后必须跟正整数！" << std::endl;
                return false;
            }
        }
        // 处理--step（每段步数N）
        else if (arg == "--step" && i + 1 < argc) {
            try {
                segment_steps = std::stoi(argv[++i]);
                if (segment_steps <= 0) {
                    std::cerr << "错误：--step 参数必须大于0！" << std::endl;
                    return false;
                }
            } catch (const std::exception& e) {
                std::cerr << "错误：--step 后必须跟正整数！" << std::endl;
                return false;
            }
        }
        // 未知参数
        else {
            std::cerr << "未知参数：" << arg << std::endl;
            return false;
        }
    }

    return true;
}

// 打印用法提示
void print_usage(const char* prog_name) {
    std::cerr << "用法：" << prog_name << " [可选参数]" << std::endl;
    std::cerr << "可选参数：" << std::endl;
    std::cerr << "  --space <数值>   空间步数（默认：20）" << std::endl;
    std::cerr << "  --time  <数值>   时间段数M（默认：50000）" << std::endl;
    std::cerr << "  --step  <数值>   每段步数N（默认：200）" << std::endl;
    std::cerr << "示例：" << std::endl;
    std::cerr << "  " << prog_name << " --time 50000 --space 20 --step 200" << std::endl;
}

int main(int argc, char* argv[]) {
    try {
        int space_steps, time_segments, segment_steps;

        // 解析参数（失败则打印用法并退出）
        if (!parse_command_line(argc, argv, space_steps, time_segments, segment_steps)) {
            print_usage(argv[0]);
            return 1;
        }

        // 初始化模型（参数从命令行传入）
        CGLE model(space_steps, time_segments, segment_steps);
        
        // 运行模拟
        model.initialize_state();
        model.run_simulation();
        std::cout << "模拟完成！参数：" << std::endl;
        std::cout << "  空间步数：" << space_steps << std::endl;
        std::cout << "  时间段数(M)：" << time_segments << std::endl;
        std::cout << "  每段步数(N)：" << segment_steps << std::endl;
        
        // 可视化结果
        auto [u1, v1, re, im] = model.get_results();
        auto [poincare_x, poincare_y] = model.get_poincare_section();
        plot_results(re, im, poincare_x, poincare_y);
    } 
    catch (const std::exception& e) {
        std::cerr << "运行错误：" << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}