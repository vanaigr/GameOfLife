cmake_minimum_required(VERSION 3.28.0)

project(GameOfLife)

set(GAME_NAME "game")
add_executable(${GAME_NAME})

include_directories("${CMAKE_SOURCE_DIR}/dependencies/include")
target_compile_definitions(${GAME_NAME} PRIVATE GLEW_STATIC)
target_compile_options(${GAME_NAME} PUBLIC -msse4.1)

file(GLOB_RECURSE GAME_SOURCES "src/*.cpp")
target_sources(${GAME_NAME} PRIVATE ${GAME_SOURCES})

target_link_libraries(${GAME_NAME} "${CMAKE_SOURCE_DIR}/dependencies/libs/GLFW/glfw3.lib")
target_link_libraries(${GAME_NAME} "${CMAKE_SOURCE_DIR}/dependencies/libs/GLEW/glew32s.lib")

target_link_libraries(${GAME_NAME} opengl32.dll gdi32.dll user32.dll kernel32.dll)
