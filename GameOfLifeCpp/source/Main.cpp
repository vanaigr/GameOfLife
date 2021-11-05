#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <math.h> 

#include "Misc.h"
#include "Vector.h"

#include "Grid.h"

#include <thread>

#include"Timer.h"
#include"AutoTimer.h"

#include"PerlinNoise.h"
#include"MedianCounter.h"

#include"ShaderLoader.h"

#include "data.h"

const uint32_t counterSampleSize = 200;
UMedianCounter 
	set{ counterSampleSize }, 
	bufferSet{ counterSampleSize }, 
	draw{ counterSampleSize }, 
	postProcessing{ counterSampleSize },
	swap{ counterSampleSize }, 
	update{ counterSampleSize }, 
	fieldUpdateWait{ counterSampleSize };

UMedianCounter microsecPerFrame{ counterSampleSize };

//#define FULLSCREEN

#ifdef FULLSCREEN
const uint32_t windowWidth = 1920, windowHeight = 1080;
#else
const uint32_t windowWidth = 800, windowHeight = 800;
#endif // FULLSCREEN

const vec2<double> windowSize{ windowWidth, windowHeight };

const uint32_t gridWidth = 2'000, gridHeight = 2'000;
const uint32_t gridSize = gridWidth * gridHeight;
const uint32_t numberOfTasks = 1;
std::unique_ptr<Field> grid;

const double cellSize_px = std::min((float)windowHeight / (float)gridHeight, (float)windowWidth / (float)gridWidth); //cell size in pixels

bool gridUpdate = true;


std::chrono::steady_clock::time_point lastScreenUpdateTime;
const uint32_t screenUpdatesPerSecond = 40;
const double screenUpdateTimeMs = 1000.0 / screenUpdatesPerSecond;

const double chasingSpeed = .5;
double deltaSizeChange = 0;
double desiredSize = 4.0, currentSize = desiredSize;

vec2<double> deltaPosChange(0, 0);
vec2<double> desiredPosition(0, 0), currentPosition = desiredPosition;


float lensDistortion = 0.17;

int32_t brushSize = 0;

bool pan = false;
enum class PaintMode : unsigned char {
	NONE,
	PAINT,
	DELETE
};
PaintMode paintMode = PaintMode::NONE;


vec2<double> mousePos(0, 0), pmousePos(0, 0);

std::chrono::steady_clock::time_point lastGridUpdateTime;
uint32_t gridUpdatesPerSecond = 1;

std::chrono::steady_clock::time_point curTime;
float r1 = 1; //w key not pressed
float r2 = 0; //normalized mouseX


GLuint frameBuffer, frameBufferTexture;

void printMouseCellInfo();

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) noexcept {
	
	if (action == GLFW_PRESS && key == GLFW_KEY_TAB) { //debug info
		const auto printC = [](const std::string label, const UMedianCounter& counter) {
			std::cout << label << '=' << counter.median()  << ", max=" << counter.max() << std::endl;
		};
		std::cout << "r1(w key not pressed)" << '=' << r1 << std::endl;
		std::cout << "r2(normalized mouse x)" << '=' << r2 << std::endl;
		printC("set", set);
		printC("draw", draw);
		printC("post processing", postProcessing);
		printC("swap", swap);
		printC("update", update);
		printC("field wait", fieldUpdateWait);

		const auto mpf = microsecPerFrame.median();
		const auto maxfps = microsecPerFrame.max();
		std::cout << "fps " << float(1'000'000 / (mpf)) << " (" << (mpf / 1'000) << "ms" << ", maximum: " << (maxfps / 1'000) << "ms)" << std::endl;
	}
	if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS) {
		printMouseCellInfo();
	}
	if (key == GLFW_KEY_ESCAPE) {
		exit(0);
	}
	if (key == GLFW_KEY_W && action == GLFW_PRESS) r1 = 0;
	else if(key == GLFW_KEY_W && action == GLFW_RELEASE) r1 = 1;
	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_ENTER) {
			grid->fill(FieldCell::DEAD);
		}
		else if (key == GLFW_KEY_SPACE) { //sace
			gridUpdate = !gridUpdate;
		}
	}

	if (key == GLFW_KEY_GRAVE_ACCENT) {
		if (action == GLFW_PRESS) {
			pan = true;
		}
		else if (action == GLFW_RELEASE) {
			pan = false;
		}
	}

	if ((key == GLFW_KEY_EQUAL) && brushSize < 30) {
		brushSize++;
	}
	else if (key == GLFW_KEY_MINUS && brushSize >= 1) {
		brushSize--;
	}
}

vec2<double> lensDistortio(vec2<double> coord, double intensity) {
	auto x = coord.x, y = coord.y;
	auto w2 = windowWidth / 2.0, h2 = windowHeight / 2.0;
	auto xc = x - w2, yc = y - h2;
	auto dist = sqrt(xc * xc + yc * yc);
	auto maxDist = sqrt(windowWidth * float(windowWidth) + windowHeight * float(windowHeight)) / 2.;
	auto distortion = dist / maxDist;
	auto newX = x - distortion * intensity * (xc);
	auto newY = y - distortion * intensity * (yc);
	return vec2<double>(newX, newY);
}

vec2<double> applyLensDistortion(vec2<double> coord) {
	return lensDistortio(coord, lensDistortion);
}

vec2<double> distortedScreenToGlobal(vec2<double> coord) {
	return (((coord - (windowSize / 2.0)) * currentSize) + (windowSize / 2.0) + currentPosition) / cellSize_px;
}


vec2<double> screenToGlobal(vec2<double> coord) {
	return distortedScreenToGlobal(applyLensDistortion(coord));
}

vec2<double> mouseToGlobal() {
	return screenToGlobal(mousePos);
}

vec2<double> globalToScreen(vec2<double> coord) {
	return ((coord * cellSize_px) - currentPosition - (windowSize / 2.0)) / currentSize + (windowSize / 2.0);
}

vec2i globalAsCell(vec2<double> coord) {
	float cx = coord.x;
	float cy = coord.y;
	int cellX = int(misc::modf(cx, gridWidth));
	int cellY = int(misc::modf(cy, gridHeight));
	return vec2i{ cellX, cellY };
}

void printMouseCellInfo() {
	const auto mouseCellCoord = globalAsCell(mouseToGlobal());
	const auto mouseCellIndex = grid->coordAsIndex(mouseCellCoord);
	const auto mouseCell = grid->cellAtCoord(mouseCellCoord);

	printf(
		"mouse is at (%d; %d), index=%d (%d mod 16), cell:%s (%d)" "\n",
		mouseCellCoord.x, mouseCellCoord.y, mouseCellIndex, mouseCellIndex % 16, fieldCell::asString(mouseCell), mouseCell
	);

	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++)
			std::cout << (fieldCell::isAlive(grid->cellAtCoord(mouseCellCoord + vec2i(j, i))) ? '1' : '0') << ' ';
		std::cout << std::endl;
	}
}

static void cursor_position_callback(GLFWwindow* window, double mousex, double mousey) noexcept {
	r2 = mousex / windowWidth;
	mousePos = vec2<double>(mousex, mousey);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) noexcept {
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

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) noexcept {
	const float scaleWheelFac = 0.04f;

	desiredSize -= desiredSize * yoffset * scaleWheelFac;

	desiredPosition = (currentPosition - (windowSize / 2.0 * currentSize) + (applyLensDistortion(mousePos) * currentSize)) - (applyLensDistortion(mousePos) * desiredSize) + (windowSize / 2.0 * desiredSize);
}

void window_size_callback(GLFWwindow* window, int width, int height) noexcept {
	/*int oldWidth, oldHeight;
	glfwGetWindowSize(window, &width, &height);

	glBindTexture(GL_TEXTURE_2D, frameBufferTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glBindTexture(GL_TEXTURE_2D, 0);*/
}

void updateState() {
	curTime = std::chrono::steady_clock::now();

	vec2<double> global = mouseToGlobal();
	vec2i cell = globalAsCell(global);

	const auto gridUpdateElapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(curTime - lastGridUpdateTime).count();
	if (gridUpdate && (gridUpdateElapsedTime >= 1000.0 / gridUpdatesPerSecond)) {
		lastGridUpdateTime = curTime;
		Timer<> t{};
		grid->updateGeneration();
		fieldUpdateWait.add(t.elapsedTime());
	}

	if (paintMode != PaintMode::NONE) {
		for (int32_t yo = -brushSize; yo <= brushSize; yo++) {
			for (int32_t xo = -brushSize; xo <= brushSize; xo++) {
				const vec2i offset{ xo, yo };
				const auto coord = cell + offset;
				if (paintMode == PaintMode::PAINT) {
					grid->setCellAtCoord(cell, FieldCell::ALIVE);
				}
				else if (paintMode == PaintMode::DELETE) {
					grid->setCellAtCoord(coord + offset, FieldCell::DEAD);
				}
			}
		}
	}

	const auto screenUpdateElapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(curTime - lastScreenUpdateTime).count();
	const uint32_t updateCount = (screenUpdateElapsedTime / screenUpdateTimeMs);
	if (updateCount > 0) {
		if (pan) {
			const vec2<double>  dmouse = applyLensDistortion(mousePos);
			const vec2<double> dpmouse = applyLensDistortion(pmousePos);
			const vec2<double> diff = dmouse - dpmouse;
			desiredPosition -= diff * currentSize;
		}

		lastScreenUpdateTime = curTime; //not quite correct

		const auto prevSize = currentSize;
		const auto prevPos = currentPosition;
		for (uint32_t i = 0; i < updateCount; i++) {
			currentSize = misc::lerp<double>(currentSize, desiredSize, chasingSpeed); //linear lerp relative to world, not viewport
			currentPosition = misc::vec2lerp(currentPosition, desiredPosition, chasingSpeed);
		}
		deltaSizeChange = misc::lerp<double>(deltaSizeChange, (currentSize - prevSize) / currentSize, chasingSpeed);
		deltaPosChange = misc::vec2lerp(deltaPosChange, (currentPosition - prevPos), chasingSpeed);

		pmousePos = mousePos;
	}
}

//delete + delete function declaration at the start of the file
auto startTime = std::chrono::steady_clock::now();

struct BufferData {
	GLFWwindow* offscreen_context = nullptr;
	GLuint bufferP;
	bool& isOffset;
	uint32_t offset;

	uint32_t currentOffset() {
		return offset * isOffset;
	}
};

int main(void)
{
	GLFWwindow* window;

	if (!glfwInit())
		return -1;

	//glfwWindowHint(GLFW_SAMPLES, 121);
	//glEnable(GL_MULTISAMPLE);

#ifdef FULLSCREEN
	window = glfwCreateWindow(windowWidth, windowHeight, "Game ofLife", glfwGetPrimaryMonitor(), NULL);
#else
	window = glfwCreateWindow(windowWidth, windowHeight, "Game ofLife", NULL, NULL);
#endif // !FULLSCREEN

	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);

    //glfwSwapInterval(0);

	GLenum err = glewInit();
	if (err != GLEW_OK)
	{
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		glfwTerminate();
		return -1;
	}

	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	//callbacks
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetScrollCallback(window, scroll_callback); 
	glfwSetWindowSizeCallback(window, window_size_callback);

	GLuint programId = glCreateProgram();
	ShaderLoader sl{};
	sl.addScreenSizeTriangleStripVertexShader();
	sl.addShaderFromProjectFilePath("shaders/fs.shader", GL_FRAGMENT_SHADER, "Main shader");

	sl.attachShaders(programId);

	glLinkProgram(programId);
	glValidateProgram(programId);

	sl.deleteShaders();

	glUseProgram(programId);

	siv::PerlinNoise noise{};
	std::srand(75489385988);

	GLuint ssbo;
	glGenBuffers(1, &ssbo);
	GLenum status;
	if ((status = glCheckFramebufferStatus(GL_FRAMEBUFFER)) != GL_FRAMEBUFFER_COMPLETE) {
		fprintf(stderr, "glCheckFramebufferStatus: error %u", status);
		return -1;
	}

	grid = std::unique_ptr<Field>{ new Field(gridWidth, gridHeight, numberOfTasks, ssbo, window, false) };

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glBufferData(GL_SHADER_STORAGE_BUFFER, misc::roundUpIntTo(grid->size_bytes() * 2, 4), NULL, GL_DYNAMIC_DRAW);
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, grid->size_bytes(), grid->rawData());
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	
	for (size_t i = 0; i < gridSize; i++) {
		const double freq = 0.05;
		const auto coord = grid->indexAsCoord(i);
		//const auto a = misc::map<double>(noise.accumulatedOctaveNoise2D_0_1(coord.x * freq, coord.y * freq, 2), 0, 1, -0.4, 0.9);
		if ((std::rand() / (double)0x7fff) /*+ a*/ > 0.6) grid->setCellAtIndex(i, FieldCell::ALIVE);
	}

	grid->startAllGridTasks();


	//uniforms set
	glUniform1i(glGetUniformLocation(programId, "width"), windowWidth);
	glUniform1i(glGetUniformLocation(programId, "height"), windowHeight);

	glUniform1i(glGetUniformLocation(programId, "gridWidth"), gridWidth);
	glUniform1i(glGetUniformLocation(programId, "gridHeight"), gridHeight);
	glUniform1ui(glGetUniformLocation(programId, "gridWidth_actual"), grid->width_actual());

	glUniform1f(glGetUniformLocation(programId, "cellSize_px"), cellSize_px);

	GLint sizeP = glGetUniformLocation(programId, "size");
	GLint posP = glGetUniformLocation(programId, "pos");

	GLint mousePosP = glGetUniformLocation(programId, "mousePos");

	glUniform1f(glGetUniformLocation(programId, "lensDistortion"), lensDistortion);

	GLint is2ndBufferP = glGetUniformLocation(programId, "is2ndBuffer");
	glUniform1ui(glGetUniformLocation(programId, "bufferOffset_bytes"), grid->size_bytes());

	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &frameBufferTexture);
	glBindTexture(GL_TEXTURE_2D, frameBufferTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, windowWidth, windowHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glBindTexture(GL_TEXTURE_2D, 0);

	glGenFramebuffers(1, &frameBuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, frameBufferTexture, 0);
	{
		GLenum status;
		if ((status = glCheckFramebufferStatus(GL_FRAMEBUFFER)) != GL_FRAMEBUFFER_COMPLETE) {
			fprintf(stderr, "glCheckFramebufferStatus: error %u", status);
			return 0;
		}
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	GLuint postProcessingProg = glCreateProgram();
	ShaderLoader ppsl{};
	ppsl.addScreenSizeTriangleStripVertexShader();
	ppsl.addShaderFromProjectFilePath("shaders/pp_chromAbb.shader", GL_FRAGMENT_SHADER, "Postprocessing chromatic abberation shader");
	//ppsl.addShaderFromProjectFilePath("shaders/pp_vignette.shader", GL_FRAGMENT_SHADER, "Postprocessing vignette shader");

	ppsl.attachShaders(postProcessingProg);
	glLinkProgram(postProcessingProg);
	glValidateProgram(postProcessingProg);
	ppsl.deleteShaders();

	glUseProgram(postProcessingProg);

	GLint frameBufferP = glGetUniformLocation(postProcessingProg, "frameBuffer");
	GLint textureSizeP = glGetUniformLocation(postProcessingProg, "textureSize");

	GLint ppDeltaScaleChangeP = glGetUniformLocation(postProcessingProg, "deltaScaleChange");
	GLint deltaOffsetChangeP = glGetUniformLocation(postProcessingProg, "deltaOffsetChange");

	GLint ppSizeP = glGetUniformLocation(postProcessingProg, "size");

	glUniform1f(glGetUniformLocation(postProcessingProg, "lensDistortion"), lensDistortion);


	setProgramId(programId);
	setSSBOHandle(ssbo);

	GLint r1P = glGetUniformLocation(programId, "r1");
	GLint r2P = glGetUniformLocation(programId, "r2");

	lastGridUpdateTime = lastScreenUpdateTime = curTime = std::chrono::steady_clock::now();

	auto tim = std::chrono::steady_clock::now();

	while (!glfwWindowShouldClose(window))
	{
		glUseProgram(programId);

		Timer<> frame{};
		{

			Timer<> t{};
			glUniform1d(sizeP, currentSize);

			glUniform2d(posP, currentPosition.x, currentPosition.y); 

			glUniform2f(mousePosP, mousePos.x, mousePos.y);

			glUniform1f(r1P, r1);
			glUniform1f(r2P, r2);

			glUniform1ui(is2ndBufferP, !grid->isFieldBufferWriteOffset());

			set.add(t.elapsedTime());
		}

		while ((err = glGetError()) != GL_NO_ERROR)
		{
			fprintf(stderr, "Error %u\n", err);
		}

		//glClear(GL_COLOR_BUFFER_BIT);
		{
			Timer<> t{};
			glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
			glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			draw.add(t.elapsedTime());
		}

		{
			Timer<> t{};
			glUseProgram(postProcessingProg);
			glBindTexture(GL_TEXTURE_2D, frameBufferTexture);

			glUniform1i(frameBufferP, 0); 
			glUniform2f(textureSizeP, windowWidth, windowHeight);

			glUniform1f(ppDeltaScaleChangeP, deltaSizeChange);
			glUniform2f(deltaOffsetChangeP, deltaPosChange.x, deltaPosChange.y);

			glUniform1d(ppSizeP, currentSize);

			glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
			glBindTexture(GL_TEXTURE_2D, 0);
			postProcessing.add(t.elapsedTime());
		}

		{
			Timer<> t{};
			glfwSwapBuffers(window);
			swap.add(t.elapsedTime());
		}

		{
			Timer<> t{};
			glfwPollEvents();

			updateState();
			update.add(t.elapsedTime());
		}

		//std::cout << frame.elapsedTime() << ::std::endl;
		microsecPerFrame.add(frame.elapsedTime());
	}

	glDeleteTextures(1, &frameBufferTexture);
	glDeleteFramebuffers(1, &frameBuffer);

	glfwTerminate();
	return 0;
}