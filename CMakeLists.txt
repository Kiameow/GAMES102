cmake_minimum_required(VERSION 3.21)

project(GAMES102 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Qt 5.15 查找策略（方便协作者开箱即用）：
#   1. 优先使用 -DCMAKE_PREFIX_PATH=... 或环境变量 QT_PREFIX 指定的路径；
#   2. 都没有时，回退到常见默认安装位置（如 C:/Qt/5.15.2/msvc2019_64）。
# ---------------------------------------------------------------------------
find_package(Qt5 5.15 QUIET COMPONENTS Widgets)
if(NOT Qt5_FOUND)
    if(NOT CMAKE_PREFIX_PATH)
        foreach(candidate IN ITEMS
                "C:/Qt/5.15.2/msvc2019_64"
                "C:/Qt/5.15.2/msvc2019")
            if(EXISTS "${candidate}/lib/cmake/Qt5Widgets")
                set(CMAKE_PREFIX_PATH "${candidate}")
                break()
            endif()
        endforeach()
    endif()
    find_package(Qt5 5.15 REQUIRED COMPONENTS Widgets)
endif()

# ---------------------------------------------------------------------------
# 主程序
# ---------------------------------------------------------------------------
add_executable(GAMES102 WIN32
    src/main.cpp
    src/MainWindow.h
    src/MainWindow.cpp
    src/PlotWidget.h
    src/PlotWidget.cpp
    src/math/CurveMath.h
    src/math/CurveMath.cpp
)

target_include_directories(GAMES102 PRIVATE src third_party/eigen)

target_link_libraries(GAMES102 PRIVATE Qt5::Widgets)

if(WIN32)
    # WIN32 子系统下需要 qtmain 提供 WinMain 入口
    target_link_libraries(GAMES102 PRIVATE Qt5::WinMain)
endif()

if(MSVC)
    # 源码中的中文字符串按 UTF-8 编译，避免在中文系统上乱码
    target_compile_options(GAMES102 PRIVATE /utf-8)
endif()

# ---------------------------------------------------------------------------
# 数学库自测（纯 C++，无 Qt 依赖，验证曲线算法数值正确性）
# 构建: cmake --build build\msvc2022 --config Debug --target curve_math_test
# 运行: build\msvc2022\Debug\curve_math_test.exe
# ---------------------------------------------------------------------------
add_executable(curve_math_test
    tests/curve_math_test.cpp
    src/math/CurveMath.h
    src/math/CurveMath.cpp
)
target_include_directories(curve_math_test PRIVATE src third_party/eigen)
target_compile_features(curve_math_test PRIVATE cxx_std_17)
if(MSVC)
    target_compile_options(curve_math_test PRIVATE /utf-8)
endif()
