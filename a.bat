@echo off
cmake -B build -G Ninja && cmake --build build && .\build\game