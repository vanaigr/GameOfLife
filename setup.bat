@echo off

md dependencies\libs\GLFW
md dependencies\libs\GLEW
md dependencies\include\GLFW
md dependencies\include\GLEW
md glslang\

echo setting up glfw
curl -sSL -o glfw-3.3.2.bin.WIN64.zip https://github.com/glfw/glfw/releases/download/3.3.2/glfw-3.3.2.bin.WIN64.zip
tar -xf glfw-3.3.2.bin.WIN64.zip
move glfw-3.3.2.bin.WIN64\lib-vc2019\glfw3.lib dependencies\libs\GLFW
move glfw-3.3.2.bin.WIN64\include\GLFW\glfw3.h dependencies\include\GLFW
del glfw-3.3.2.bin.WIN64.zip
rmdir /s /q glfw-3.3.2.bin.WIN64

echo settng up glew
curl -sSL -o glew-2.2.0-win32.zip https://github.com/nigels-com/glew/releases/download/glew-2.2.0/glew-2.2.0-win32.zip
tar -xf glew-2.2.0-win32.zip
move glew-2.2.0\lib\Release\x64\glew32s.lib dependencies\libs\GLEW
move glew-2.2.0\include\GL\glew.h dependencies\include\GLEW
del glew-2.2.0-win32.zip
rmdir /s /q glew-2.2.0

echo settng up glslangValidator
curl -sSL -o glslang-master-windows-x64-Release.zip https://github.com/KhronosGroup/glslang/releases/download/master-tot/glslang-master-windows-x64-Release.zip
tar -xf glslang-master-windows-x64-Release.zip -C glslang
move glslang\bin\glslangValidator.exe .\dependencies
del glslang-master-windows-x64-Release.zip
rmdir /s /q glslang
