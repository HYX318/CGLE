cmake_minimum_required(VERSION 3.15)
project(CGLE)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找Python3（包含NumPy）
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
if(NOT Python3_FOUND)
    message(FATAL_ERROR "Python3 with NumPy not found. Please install Python3 and NumPy.")
endif()

# 包含matplotlibcpp头文件路径（根据实际路径调整）
include_directories(
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
    "./matplotlib-cpp/"  # 替换为你的matplotlibcpp.h所在目录
)

# 添加可执行文件
add_executable(CGLE main.cpp CGLE.cpp)

# 链接Python库
target_link_libraries(CGLE 
    PRIVATE 
    Python3::Python 
    Python3::NumPy
)

# 可选：输出配置信息，用于调试路径问题
message(STATUS "Python3 Include Dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "NumPy Include Dir: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "Python3 Libraries: ${Python3_LIBRARIES}")





