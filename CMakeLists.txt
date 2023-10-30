cmake_minimum_required(VERSION 3.28.0)

project(GameOfLife)

set(GAME_NAME "game")
add_executable(${GAME_NAME})

include_directories("${CMAKE_SOURCE_DIR}/dependencies/include")
target_compile_definitions(${GAME_NAME} PRIVATE GLEW_STATIC)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(
        ${GAME_NAME} PRIVATE 
        -msse4.1
        -pedantic -Wall -Wextra -Werror 
        -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function
        -Wno-unused-but-set-variable -Wno-unused-label -Wno-unused-private-field
    )
endif()
target_compile_features(${GAME_NAME} PUBLIC cxx_std_17)

file(GLOB_RECURSE GAME_SOURCES "src/*.cpp")
target_sources(${GAME_NAME} PRIVATE ${GAME_SOURCES})

target_link_libraries(${GAME_NAME} "${CMAKE_SOURCE_DIR}/dependencies/libs/GLFW/glfw3.lib")
target_link_libraries(${GAME_NAME} "${CMAKE_SOURCE_DIR}/dependencies/libs/GLEW/glew32s.lib")

target_link_libraries(${GAME_NAME} opengl32.dll gdi32.dll user32.dll kernel32.dll)
