#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <math.h> 

#include "ShaderUtil.h"
#include "misc.h"
#include "Vector.h"

//delete
void printFrameRate();

const int windowWidth = 1920, windowHeight = 1080;

const uint32_t gridWidth = 200, gridHeight = 200;
const size_t gridSize = gridWidth * gridHeight;
const float size = std::min((float)windowHeight / (float)gridHeight, (float)windowWidth / (float)gridWidth); //cell size in pixels

enum class GridCell : unsigned char
{
	DEAD,
	ALIVE,
	WALL
};


GridCell g1[gridSize] = { GridCell::DEAD };
GridCell g2[gridSize] = { GridCell::DEAD };

GridCell (*grid)[gridSize] = &g1;
GridCell(*gridBuffer)[gridSize] = &g2;

const uint32_t packedGridSize = (gridWidth * gridHeight) / (32 / 2) + 1;
uint32_t packedGrid[packedGridSize] = { 0 };

bool gridUpdate = true;

vec2 offset(0, 0);
//float tx = 0, ty = 0;
float desiredScale = 2.0, currentScale = desiredScale;

float lensDistortion = -0.1;


enum class BrushMode : bool {
	CELL,
	WALL
};

int32_t brushSize = 0;
BrushMode brushMode = BrushMode::CELL;


bool pan = false;
enum class PaintMode : unsigned char {
	NONE,
	PAINT,
	DELETE
};
PaintMode paintMode = PaintMode::NONE;


vec2 mousePos(0, 0), pmousePos(0, 0);

std::chrono::steady_clock::time_point lastUpdateTime, curTime;
uint32_t gridUpdatesPerSecond = 20;

float r1 = 1; //w key not pressed
float r2 = 0; //normalized mouseX

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS && key == GLFW_KEY_TAB) { //debug info
		std::cout
			<< "r1=" << r1 << std::endl
			<< "r2=" << r2 << std::endl;
	}
	if (key == GLFW_KEY_ESCAPE) {
		printFrameRate();
		exit(0);
	}
	if (key == GLFW_KEY_W && action == GLFW_PRESS) r1 = 0;
	else if(key == GLFW_KEY_W && action == GLFW_RELEASE) r1 = 1;
	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_ENTER) {
			std::fill(std::begin(*grid), std::end(*grid), GridCell::DEAD);
		}
		else if (key == GLFW_KEY_1) {
			brushMode = BrushMode::CELL;
		}
		else if (key == GLFW_KEY_2) {
			brushMode = BrushMode::WALL;
		}
		/*else if (key == 51) {
			mode = 2;
		}*/
		else if (key == GLFW_KEY_SPACE) { //sace
			gridUpdate = !gridUpdate;
		}

		return;
	}
	else if ((key == GLFW_KEY_EQUAL) && brushSize < 10) {
		brushSize++;
	}
	else if (key == GLFW_KEY_MINUS && brushSize >= 1) {
		brushSize--;
	}
}

static void cursor_position_callback(GLFWwindow* window, double mousex, double mousey)
{
	r2 = mousex / windowWidth;
	mousePos = vec2(mousex, mousey);
	//std::cout << 'm' << mouseX << ',' << mouseX << std::endl;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	//std::cout << 'b' << (action == GLFW_PRESS) << std::endl;

	if (action == GLFW_PRESS) {
		if (button == GLFW_MOUSE_BUTTON_MIDDLE) pan = true;
		else if (button == GLFW_MOUSE_BUTTON_LEFT) paintMode = PaintMode::PAINT;
		else if (button == GLFW_MOUSE_BUTTON_RIGHT) paintMode = PaintMode::DELETE;
	}
	else if (action == GLFW_RELEASE) {
		if (button == GLFW_MOUSE_BUTTON_MIDDLE) pan = false;
		else if (button == GLFW_MOUSE_BUTTON_LEFT) paintMode = PaintMode::NONE;
		else if (button == GLFW_MOUSE_BUTTON_RIGHT) paintMode = PaintMode::NONE;
	}
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	const float scaleWheelFac = 0.04f;

	desiredScale += desiredScale * yoffset * scaleWheelFac;

	//std::cout << 's' << desiredScale << std::endl;
}

/*void updateGrid() {
	const GridCell(&g)[] = (GridCell(&)[]) grid;
	const size_t gw = gridWidth;
	const auto isCell = [g](size_t index) -> bool { return g[index] == GridCell::ALIVE; };
	const auto updatedCell = [](const GridCell cell, const unsigned char aliveNeighbours) -> GridCell { 
		if (cell == GridCell::DEAD && aliveNeighbours == 3) {
			return GridCell::ALIVE;
		}
		else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
			return GridCell::DEAD;
		}
	};
	const size_t topRow = gw * (gridHeight - 1);
	const size_t lastElement = gw - 1;
	unsigned char r1[gridWidth];
	unsigned char r2[gridWidth];
	unsigned char (*topRowNeighbours)[gridWidth] = &r1;
	unsigned char (*curRowNeighbours)[gridWidth] = &r2;

	for (size_t i = 1; i < gridWidth-1; i++) {
		(*topRowNeighbours)[i] = isCell(topRow + i - 1) + isCell(topRow + i) + isCell(topRow + i + 1);
	}
	(*topRowNeighbours)[0] = isCell(topRow + lastElement) + isCell(topRow + 0) + isCell(topRow + 1);
	(*topRowNeighbours)[lastElement] = isCell(topRow + lastElement - 1) + isCell(topRow + lastElement) + isCell(topRow + 0);

	for (size_t row = 0; row < gridHeight-1; row++) {
		const size_t rowIndex = gw * row;
		for (size_t i = 1; i < gridWidth - 1; i++) {
			const size_t index = i + rowIndex;
			if (g[index] != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCell(index + gw - 1) + isCell(index + gw) + isCell(index + gw + 1);
				const unsigned char aliveNeighbours = (*topRowNeighbours)[i] + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;

				(*gridBuffer)[index] = updatedCell(g[index], aliveNeighbours);
			}

			(*curRowNeighbours)[i] = isCell(index - 1) + isCell(index) + isCell(index + 1);
		}

		{ //first element
			const unsigned char bottomRowNeighbours = isCell((rowIndex + gw) + lastElement) + isCell((rowIndex + gw) + 0) + isCell((rowIndex + gw) + 1);
			const unsigned char aliveNeighbours = (*topRowNeighbours)[0] + isCell(rowIndex + lastElement) + isCell(rowIndex + 1) + bottomRowNeighbours;

			(*gridBuffer)[rowIndex] = updatedCell(g[rowIndex], aliveNeighbours);

			(*curRowNeighbours)[0] = isCell(rowIndex + lastElement) + isCell(rowIndex + 0) + isCell(rowIndex + 1);
		}
		{ //last element
			const unsigned char bottomRowNeighbours = isCell((rowIndex + gw) + lastElement - 1) + isCell((rowIndex + gw) + lastElement) + isCell((rowIndex + gw) + 1);
			const unsigned char aliveNeighbours = (*topRowNeighbours)[rowIndex + lastElement] + isCell(rowIndex + lastElement - 1) + isCell(rowIndex + 0) + bottomRowNeighbours;

			(*gridBuffer)[rowIndex] = updatedCell(g[rowIndex], aliveNeighbours);

			(*curRowNeighbours)[0] = isCell(rowIndex + lastElement - 1) + isCell(rowIndex + lastElement) + isCell(rowIndex + 0);
		}

		auto tmp = topRowNeighbours;
		topRowNeighbours = curRowNeighbours;
		curRowNeighbours = tmp;
	}

	{ //bottom row
		for (size_t i = 1; i < gridWidth - 1; i++) {
			const size_t index = i + topRow;
			const unsigned char bottomRowNeighbours = isCell(i - 1) + isCell(i) + isCell(i + 1);
			const unsigned char aliveNeighbours = (*topRowNeighbours)[i] + isCell(index - 1) + isCell(index + 1) + bottomRowNeighbours;

			(*gridBuffer)[index] = updatedCell(g[index], aliveNeighbours);
		}

		{ //first element
			const unsigned char bottomRowNeighbours = isCell(0 + lastElement) + isCell(0 + 0) + isCell(0 + 1);
			const unsigned char aliveNeighbours = (*topRowNeighbours)[0] + isCell(topRow + lastElement) + isCell(topRow + 1) + bottomRowNeighbours;

			(*gridBuffer)[topRow + 0] = updatedCell(g[topRow + 0], aliveNeighbours);
		}
		{ //last element
			const unsigned char bottomRowNeighbours = isCell(0 + lastElement-1) + isCell(0 + lastElement) + isCell(0 + 0);
			const unsigned char aliveNeighbours = (*topRowNeighbours)[lastElement] + isCell(topRow + lastElement - 1) + isCell(topRow + 0) + bottomRowNeighbours;

			(*gridBuffer)[topRow + lastElement] = updatedCell(g[topRow + lastElement], aliveNeighbours);
		}
	}

	auto tmp = grid;
	grid = gridBuffer;
	gridBuffer = tmp;
}*/

GridCell updatedCell(const uint32_t index) {
	GridCell cell = (*grid)[index];

	if (cell != GridCell::WALL) {
		uint32_t aliveNeighbours = 0;

		for (int yo = -1; yo <= 1; yo++) {
			for (int xo = -1; xo <= 1; xo++) {
				if (xo != 0 || yo != 0) {
					const auto offset = vec2i{ xo, yo };
					int x = ((index + xo) + gridWidth) % gridWidth,
						y = (((index / gridWidth) + yo) + gridHeight) % gridHeight;
					int offsetedIndex = x + gridWidth * y;
					if ((*grid)[offsetedIndex] == GridCell::ALIVE)
						aliveNeighbours++;
				}
			}
		}

		if (cell == GridCell::DEAD && aliveNeighbours == 3) {
			return GridCell::ALIVE;
		}
		else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
			return GridCell::DEAD;
		}
	}

	return cell;
}

void updateGrid() {
	for (size_t i = 0; i < gridSize; i++) {
		(*gridBuffer)[i] = updatedCell(i);
	}

	auto tmp = grid;
	grid = gridBuffer;
	gridBuffer = tmp;
}

void updateGrid_() {
	const auto &g = *grid;
	const auto gw = gridWidth;
	const auto gh = gridHeight;
	const auto updatedCell = [](const GridCell &cell, const unsigned char aliveNeighbours) -> GridCell {
		if (cell == GridCell::DEAD && aliveNeighbours == 3) {
			return GridCell::ALIVE;
		}
		else if (cell == GridCell::ALIVE && (aliveNeighbours < 2 || aliveNeighbours > 3)) {
			return GridCell::DEAD;
		}
		else return cell;
	};
	const auto isCell = [&g](size_t index) -> bool { return g[index] == GridCell::ALIVE; };
	const auto asIndex = [gw, gh](const int32_t row, const int32_t i) -> size_t { return size_t(mod(row, gh) * int32_t(gw) + mod(i, gw)); };
	const auto isCellCoord = [&asIndex, &isCell](const int32_t row, const int32_t i) -> bool { return isCell(asIndex(row, i)); };

	unsigned char r1[gw];
	unsigned char r2[gw];
	unsigned char(*topRowNeighbours)[gridWidth] = &r1;
	unsigned char(*curRowNeighbours)[gridWidth] = &r2;

	const int topRow = (gh - 1) * gw;
	for (size_t i = 0; i < gw; i++) {
		(*topRowNeighbours)[i] = isCellCoord(topRow, i - 1) + isCellCoord(topRow, i) + isCellCoord(topRow, i + 1);
	}

	for (int32_t row = 0; row < gh; row++) {
		for (int32_t i = 0; i < gw; i++) {
			const size_t index = asIndex(row, i);
			if (g[index] != GridCell::WALL) {
				const unsigned char bottomRowNeighbours = isCellCoord(row+1, i-1) + isCellCoord(row + 1, i) + isCellCoord(row + 1, i + 1);
				//const unsigned char topRowNeighbours = isCellCoord(row - 1, i - 1) + isCellCoord(row - 1, i) + isCellCoord(row - 1, i + 1);
				const unsigned char aliveNeighbours = (*topRowNeighbours)[i] + isCellCoord(row, i - 1) + isCellCoord(row, i + 1) + bottomRowNeighbours;

				(*gridBuffer)[index] = updatedCell(g[index], aliveNeighbours);

				//std::cout << (int)aliveNeighbours << " at " << row << ',' << i << std::endl;
			}

			(*curRowNeighbours)[i] = isCellCoord(row, i - 1) + isCellCoord(row, i) + isCellCoord(row, i + 1);
		}

		//auto tmp2 = topRowNeighbours;
		//topRowNeighbours = curRowNeighbours;
		//curRowNeighbours = tmp2;
	}

	auto tmp = grid;
	grid = gridBuffer;
	gridBuffer = tmp;
}

vec2 lensDistortio(vec2 coord, float intensity) {
	float x = coord.x, y = coord.y;
	float w2 = windowWidth / 2., h2 = windowHeight / 2;
	float xc = x - w2, yc = y - h2;
	float dist = sqrt(xc * xc + yc * yc);
	float maxDist = sqrt(windowWidth * float(windowWidth) + windowHeight * float(windowHeight)) / 2.;
	float distortion = dist / maxDist;
	float newX = x - distortion * intensity * (xc);
	float newY = y - distortion * intensity * (yc);
	return vec2(newX, newY);
}

vec2 applyLensDistortion(vec2 coord) {
	return lensDistortio(coord, lensDistortion);
}

vec2 distortedScreenToGlobal(vec2 coord) {
	return vec2(((coord.x - windowWidth / 2.) / currentScale + windowWidth / 2. - offset.x) / size, ((coord.y - windowHeight / 2.) / currentScale + windowHeight / 2. - offset.y) / size);
}


vec2 screenToGlobal(vec2 coord) {
	return distortedScreenToGlobal(applyLensDistortion(coord));
}

vec2 mouseToGlobal() {
	return screenToGlobal(mousePos);
}


vec2 globalToScreen(vec2 coord) {
	//x = (sx - width / 2f)/scale + width / 2f - tx;
	//(x + tx - width/2f) * scale + width/2f = sx;

	return vec2((coord.x * size + offset.x - windowWidth / 2.) * currentScale + windowWidth / 2., (coord.y * size + offset.y - windowHeight / 2.) * currentScale + windowHeight / 2.);
}

vec2i globalAsCell(vec2 coord) {
	float cx = coord.x;
	float cy = coord.y;
	//int cellX = (int)(fmod(cx,gridWidth) + gridWidth) % gridWidth;
	//int cellY = (int)(fmod(cy, gridHeight) + gridHeight) % gridHeight;
	int cellX = int(modf(cx, gridWidth));
	int cellY = int(modf(cy, gridHeight));
	return vec2i{ cellX, cellY };
}


void updateState() {
	curTime = std::chrono::steady_clock::now();

	vec2 global = mouseToGlobal();
	vec2i cell = globalAsCell(global);
	const int x = cell.x, y = cell.y;

	if (paintMode == PaintMode::PAINT) {
		for (int32_t yo = -brushSize; yo <= brushSize; yo++) {
			for (int32_t xo = -brushSize; xo <= brushSize; xo++) {
				const size_t index = ((x + xo + gridWidth) % gridWidth) + (gridWidth * ((y + yo + gridHeight) % gridHeight));
				if (brushMode == BrushMode::WALL) { (*grid)[index] = GridCell::WALL; }
				else if (brushMode == BrushMode::CELL && (*grid)[index] != GridCell::WALL) { (*grid)[index] = GridCell::ALIVE; }
			}
		}
	}
	else if (paintMode == PaintMode::DELETE) {
		for (int32_t xo = -brushSize; xo <= brushSize; xo++) {
			for (int32_t yo = -brushSize; yo <= brushSize; yo++) {
				(*grid)[((x + xo + gridWidth) % gridWidth) + (gridWidth * ((y + yo + gridHeight) % gridHeight))] = GridCell::DEAD;
			}
		}
	}

	auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(curTime - lastUpdateTime).count();
	if (gridUpdate && (elapsedTime >= 1000.0 / gridUpdatesPerSecond)) {
		lastUpdateTime = curTime;
		updateGrid();
	}

	currentScale = lerp(currentScale, desiredScale, 0.3);
	if (pan) {
		const vec2  dmouse = applyLensDistortion(mousePos);
		const vec2 dpmouse = applyLensDistortion(pmousePos);
		const vec2 diff = dmouse - dpmouse;
		offset += diff / currentScale;
		//tx += (mouseX - pmouseX) / currentScale;
		//ty += (mouseY - pmouseY) / currentScale;
	}
	pmousePos = mousePos;
}

//delete + delete function declaration at the start of the file
auto startTime = std::chrono::steady_clock::now();
uint64_t frameCount = 0;

void printFrameRate() {
	auto curTime = std::chrono::steady_clock::now();

	auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(curTime - startTime).count();
	std::cout << "frame ms " << (elapsedTime / frameCount) << std::endl;
}

int main(void)
{
	std::srand(75489385988);
	for (size_t i = 0; i < gridSize; i++) {
		/*auto row = i / gridWidth;
		auto column = mod(i, gridWidth);
		if ((row == 3) || (row == 4) || (row == 5))
			if ((column == 3) || (column == 4) || (column == 5))
				(*grid)[i] = GridCell::ALIVE;*/
		if((std::rand() / (double) 0x7fff) > 0.34)(*grid)[i] = GridCell::ALIVE;
	}

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	//glfwWindowHint(GLFW_SAMPLES, 121);
	//glEnable(GL_MULTISAMPLE);

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(windowWidth, windowHeight, "Game ofLife", glfwGetPrimaryMonitor(), NULL);
	//window = glfwCreateWindow(windowWidth, windowHeight, "Game ofLife", NULL, NULL);

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	GLenum err = glewInit();
	if (err != GLEW_OK)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		glfwTerminate();
		return -1;
	}

	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));


	// TODO: Create and compile shaders here (vertex and fragment shaders)
	// and finally draw something with modern OpenGL!
	ShaderUtil shaderUtil;
	unsigned int progId;
	shaderUtil.Load("shaders/vs.shader", "shaders/fs.shader", progId);

	GLuint ssbo;
	glGenBuffers(1, &ssbo);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(uint32_t) * packedGridSize, packedGrid, GL_DYNAMIC_DRAW);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	shaderUtil.Use();

	//uniforms set
	glUniform1i(glGetUniformLocation(progId, "width"), GLint(windowWidth));
	glUniform1i(glGetUniformLocation(progId, "height"), GLint(windowHeight));

	glUniform1i(glGetUniformLocation(progId, "gridWidth"), GLint(gridWidth));
	glUniform1i(glGetUniformLocation(progId, "gridHeight"), GLint(gridHeight));

	glUniform1f(glGetUniformLocation(progId, "size"), GLfloat(size));
	glUniform1f(glGetUniformLocation(progId, "lensDistortion"), GLfloat(lensDistortion));

	GLint curScaleP = glGetUniformLocation(progId, "currentScale");
	GLint txP = glGetUniformLocation(progId, "tx");
	GLint tyP = glGetUniformLocation(progId, "ty");


	GLint r1P = glGetUniformLocation(progId, "r1");
	GLint r2P = glGetUniformLocation(progId, "r2");
	//GLint fieldP = glGetUniformLocation(progId, "field");


	//callbacks
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, scroll_callback);


	lastUpdateTime = curTime = std::chrono::steady_clock::now();

	/* Loop until the user closes the window */
	auto tim = std::chrono::steady_clock::now();
	
	while (!glfwWindowShouldClose(window))
	{
		glUniform1f(curScaleP, GLfloat(currentScale));
		glUniform1f(txP, GLfloat(offset.x));
		glUniform1f(tyP, GLfloat(offset.y));

		glUniform1f(r1P, r1);
		glUniform1f(r2P, r2);

		for (size_t i = 0; i < gridSize; i++) {
			unsigned int isWall = (*grid)[i] == GridCell::WALL;
			unsigned int isCell = (*grid)[i] == GridCell::ALIVE;
			unsigned int mask = (((isWall << 1) & 0b10) | isCell) & 0b11;

			unsigned int arrIndex = i / (32 / 2);
			unsigned int arrShift = (i % (32 / 2)) * 2;
			packedGrid[arrIndex] = (packedGrid[arrIndex] & ~(3 << (arrShift))) | (mask << (arrShift));
		}

		GLuint ssbo;
		glGenBuffers(1, &ssbo);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
		glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(uint32_t) * packedGridSize, packedGrid, GL_DYNAMIC_DRAW);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
		

		while ((err = glGetError()) != GL_NO_ERROR)
		{
			std::cout << err << std::endl;
		}

		//glClear(GL_COLOR_BUFFER_BIT);

		glDrawArrays(GL_POINTS, 0, 1);

		glfwSwapBuffers(window);

		glfwPollEvents();

		updateState();

		frameCount++;
	}

	shaderUtil.Delete();

	glfwTerminate();
	return 0;
}