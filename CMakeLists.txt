# CMakeList.txt: 顶层 CMake 项目文件，在此处执行全局配置
# 并包含子项目。
#
cmake_minimum_required (VERSION 3.8)

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

set(CMAKE_CXX_STANDARD 17)
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Building in Debug mode.")
    set(VTK_DIR C:/source/VTK/lib/cmake/vtk-9.3)
    set(ITK_DIR C:/source/ITK/lib/cmake/ITK-5.1)

elseif(CMAKE_BUILD_TYPE MATCHES Release)
    message(STATUS "Building in Release mode.")
    set(VTK_DIR C:/source/VTK/lib/cmake/vtk-9.3)
    set(ITK_DIR C:/source/ITK/lib/cmake/ITK-5.1)
else()
    message(STATUS "Building in an unknown mode: ${CMAKE_BUILD_TYPE}")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories(${CMAKE_SOURCE_DIR}/src/head)


# 使用系统的 Eigen3
set(ITK_USE_SYSTEM_EIGEN ON)

# 使用标准的 find_package 来查找 Eigen3
set(Eigen3_DIR "C:/source/eigen-3.4.0/build")
find_package(Eigen3 REQUIRED)

# 包含 Eigen 的头文件
include_directories(${EIGEN3_INCLUDE_DIR})

# 如果需要手动包含 UseEigen3.cmake，保留这一行
include("C:/source/eigen-3.4.0/cmake/UseEigen3.cmake")

# 设置 ITK 和 VTK 的查找路径
set(CMAKE_PREFIX_PATH "C:/source/ITK" "C:/source/VTK" ${CMAKE_PREFIX_PATH})

project ("liver_cut")

# 包含子项目。
add_subdirectory ("src")
