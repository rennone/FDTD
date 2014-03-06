#要求するcmakeのバージョン
cmake_minimum_required(VERSION 2.8)

file(GLOB srcFile ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

#プロジェクト名
project(fdtd)

#GLUT, OpenGLのフォルダを見つける?

find_package(GLUT)
find_package(OpenGL)

set(CMAKE_C_FLAGS "-Wall -std=c99")

#Debug Mode
set(CMAKE_C_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_BUILD_TYPE debug)

include_directories(/opt/local/include)
include_directories(~/glew-1.10.0/include)
link_directories(~/glew-1.10.0/lib)
add_executable(fdtd ${srcFile})

target_link_libraries(fdtd /System/Library/Frameworks/GLUT.framework ${OPENGL_LIBRARY} m)
target_link_libraries(fdtd GLEW.a)

