# Game of Life

This is a high-performance implementation of [Conway's Game of Life](https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life), written in C++ with OpenGL. The simulation is optimized using SIMD instructions and multithreading.

![image](https://github.com/user-attachments/assets/8159fe09-9bab-4a1d-ab74-d888776e0fcc)

## Controls

- `Space` - pause/unpause the game
- `Enter` - clear the board
- `=` - increase brush size
- `-` - decrease brush size
- `Left Mouse Button` - draw cells
- `Right Mouse Button` - erase cells
- `Scroll Wheel` - zoom in/out
- `Middle Mouse Button` - pan
- `Esc` - exit the game

## Running the game

Your CPU must support SSE4.1 to run the game.

The project assumes the compiler is clang.

### Windows

Download GLEW, GLFW, and GLSL validator:

```batch
setup
```

Build the game:

```batch
mkdir build
cmake -B build -G Ninja
cmake --build build
```

Run the game:

```batch
.\build\game
```

### Linux

Install development versions of GLEW and GLFW 3.

Build the game:

```batch
mkdir build
cmake -B build -G Ninja
cmake --build build
```

Run the game:

```batch
./build/game
```
