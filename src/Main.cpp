#include "glew.h"
#include <GLFW/glfw3.h>
#include <cmath>
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

#include<atomic>
#include<mutex>

#include<vector>

#include<type_traits>

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

struct WindowSize {
    vec2i windowSize;
    vec2d windowSizeD;
    vec2d winSizeH;
    double normalizeCoeff;
    double lensDistortionCoeff;
};
static WindowSize winSize;

struct Space {
    vec2d vpPosDesired, vpPos; //center of the viewport's position in global space
    double vpSize; //total vertical size of viewport in global space
};
static Space space;

static double zoomLevel;
static vec2d dVpPosDesired, dVpPos;
static double dVpSizeDesired, dVpSize;

vec2d desiredZoomPoint, zoomPoint; //TODO: init
bool isZoomChanged;

const uint32_t gridWidth = 50000, gridHeight = 50000;
const uint32_t gridSize = gridWidth * gridHeight;
const uint32_t numberOfTasks = 10;
std::unique_ptr<Field> grid;

static bool gridUpdate = true;

float lensDistortion = 0.07;

int32_t brushSize = 0;

bool pan = false;
enum class PaintMode : unsigned char {
    NONE,
    PAINT,
    DELETE
};
PaintMode paintMode = PaintMode::NONE;


vec2d mousePos(0), pmousePos(0, 0);

std::chrono::steady_clock::time_point lastGridUpdateTime;
uint32_t gridUpdatesPerSecond = 5;

std::chrono::steady_clock::time_point curTime;
float r1 = 1; //w key not pressed
float r2 = 0; //normalized mouseX


GLuint frameBuffer, frameBufferTexture;


std::chrono::steady_clock::time_point lastScreenUpdateTime;
const uint32_t screenUpdatesPerSecond = 40;
const double screenUpdateTimeMs = 1000.0 / screenUpdatesPerSecond;

const double chasingSpeed = 0.4; //.4;

static GLuint packedGrid1 = 0, packedGrid2 = 0;
static size_t field_size_bytes;
static std::atomic_bool isBufferSecond{ true };
static std::mutex gpuBufferLock{ };


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
        grid->halt();
        exit(0);
    }
    if (key == GLFW_KEY_W && action == GLFW_PRESS) r1 = 0;
    if (key == GLFW_KEY_Q && action == GLFW_PRESS) isBufferSecond = !isBufferSecond;
    else if(key == GLFW_KEY_W && action == GLFW_RELEASE) r1 = 1;
    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_ENTER) {
            grid->fill(fieldCell::cellDead);
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

    if ((key == GLFW_KEY_EQUAL) && brushSize < 75) {
        brushSize++;
    }
    else if (key == GLFW_KEY_MINUS && brushSize >= 1) {
        brushSize--;
    }
}


static void updateWinSize(int const width, int const height) {
    auto const winSizeD = vec2d{ double(width), double(height) };

    auto w = WindowSize{};
    w.windowSize = vec2i{ width, height };
    w.windowSizeD = winSizeD;
    w.winSizeH = winSizeD * 0.5;
    w.normalizeCoeff = 1.0 / winSizeD.y;
    w.lensDistortionCoeff = 1.0/w.winSizeH.length() * lensDistortion;

    winSize = w;
}

static double getVpSizeDesired() {
    return std::pow(2, zoomLevel);
}
static double vpSizeToZoomLevel(double const size) {
    return std::log2(size);
}

static void updateSpace() {
    auto const size = std::min(gridWidth, gridHeight);
    auto const pos = vec2d{ gridWidth * 0.5, gridHeight * 0.5 };
    zoomLevel = vpSizeToZoomLevel(size);
    space = Space{ pos, pos, size };
}


static vec2d applyLensDistortion(vec2d const screenCoord) {
    auto const centerCoord = screenCoord - winSize.winSizeH;
    return screenCoord - centerCoord * centerCoord.length() * winSize.lensDistortionCoeff;
}

static vec2d mouseToGlobal() {
    return (applyLensDistortion(mousePos) * winSize.normalizeCoeff - 0.5)
        * space.vpSize + space.vpPos;
}

vec2i globalAsCell(vec2d coord) {
    return vec2{
        int(misc::modf(coord.x, gridWidth)),
        int(misc::modf(coord.y, gridHeight))
    };
}

void printMouseCellInfo() {
    const auto mouseCellCoord = globalAsCell(mouseToGlobal());
    const auto mouseCellIndex = grid->coordAsIndex(mouseCellCoord);
    const auto mouseCell = grid->cellAtCoord(mouseCellCoord);

    printf(
        "mouse is at (%d; %d), index=%lld (%lld mod 16), cell:%s (%d)" "\n",
        mouseCellCoord.x, mouseCellCoord.y, mouseCellIndex, mouseCellIndex % 16, fieldCell::asString(mouseCell), int(mouseCell)
    );

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++)
            std::cout << (grid->cellAtCoord(mouseCellCoord + vec2i(j, i)) ? '1' : '0') << ' ';
        std::cout << std::endl;
    }
}

static void cursor_position_callback(GLFWwindow* window, double mousex, double mousey) noexcept {
    r2 = mousex / winSize.windowSize.x;
    mousePos = vec2d(mousex, mousey);
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
    const float scaleWheelFac = 0.1f;

    zoomLevel -= yoffset * scaleWheelFac;
    isZoomChanged = true;
}

void window_size_callback(GLFWwindow* window, int width, int height) noexcept {
    /*int oldWidth, oldHeight;
    glfwGetWindowSize(window, &width, &height);

    glBindTexture(GL_TEXTURE_2D, frameBufferTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glBindTexture(GL_TEXTURE_2D, 0);*/
}

class GLFieldOutput final : public FieldOutput {
    class GLBufferedFieldOutput final : public FieldOutput {
        static void writeModification(FieldModification fm) {
            glBufferSubData(GL_SHADER_STORAGE_BUFFER,
                fm.startIndex_int * sizeOfBatch
                , fm.size_int * sizeOfBatch
                , fm.data);
        }

        class InnerFieldOutput final : public FieldOutput {
        public:
            InnerFieldOutput() = default;

            virtual void write(FieldModification fm) override {
                writeModification(fm);
            }

            virtual std::unique_ptr<FieldOutput> tryBatched() const override
            {
                return std::unique_ptr<InnerFieldOutput>(new InnerFieldOutput());
            }

            virtual std::unique_ptr<FieldOutput> batched() const override
            {
                return std::unique_ptr<InnerFieldOutput>(new InnerFieldOutput());
            }

            ~InnerFieldOutput() noexcept = default;
        };
    private:
        GLFieldOutput const& parent;
        std::unique_lock<std::mutex> lock;

        static constexpr auto sizeOfBatch = sizeof std::remove_pointer<decltype(FieldModification::data)>::type();/*
            to convert from batches to bytes
        */
    public:
        GLBufferedFieldOutput(
            GLFieldOutput const& parent_,
            std::unique_lock<std::mutex> _lock
        ) : parent{ parent_ }, lock(std::move(_lock)) {
            if (!glfwGetCurrentContext()) glfwMakeContextCurrent(parent.context);

            GLuint fieldBufferHandle;// { !isBufferSecond ^ parent.isBuffer ? packedGrid1 : packedGrid2 };
            if (parent.isBuffer) {
                if (isBufferSecond) fieldBufferHandle = packedGrid2;
                else fieldBufferHandle = packedGrid1;
            }
            else
                if (isBufferSecond) fieldBufferHandle = packedGrid1;
                else fieldBufferHandle = packedGrid2;

            glBindBuffer(GL_SHADER_STORAGE_BUFFER, fieldBufferHandle);
        }

        GLBufferedFieldOutput(GLBufferedFieldOutput const&) = delete;
        GLBufferedFieldOutput& operator=(GLBufferedFieldOutput const&) = delete;

        virtual void write(FieldModification fm) override {
            writeModification(fm);
        }

        // TODO: unsafe. What if it outlives us?
        virtual std::unique_ptr<FieldOutput> tryBatched() const override
        {
            return std::unique_ptr<InnerFieldOutput>(new InnerFieldOutput());
        }
        virtual std::unique_ptr<FieldOutput> batched() const override
        {
            return std::unique_ptr<InnerFieldOutput>(new InnerFieldOutput());
        }

        ~GLBufferedFieldOutput() noexcept {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
            glFinish();
        }
    };
private:
    GLFWwindow* context;
    bool isBuffer;
public:
    GLFieldOutput(bool isBuffer_, GLFWwindow* window_origin) : isBuffer{ isBuffer_ } {
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        auto* context_ = glfwCreateWindow(2, 2, "", NULL, window_origin);
        if (!context_) {
            ::std::cout << "window for output(current) is null";
            exit(-1);
        }
        context = context_;
    }

    GLFieldOutput(GLFieldOutput const&) = delete;
    GLFieldOutput& operator=(GLFieldOutput const&) = delete;

    std::unique_ptr<FieldOutput> tryBatched() const override {
        std::unique_lock<std::mutex> lock{ gpuBufferLock, std::try_to_lock_t{} };
        if(!lock.owns_lock()) return std::unique_ptr<FieldOutput>{};

        return std::unique_ptr<FieldOutput>(new GLBufferedFieldOutput{ *this, std::move(lock) });
    }

    std::unique_ptr<FieldOutput> batched() const override {
        std::unique_lock<std::mutex> lock{ gpuBufferLock };
        return std::unique_ptr<FieldOutput>(new GLBufferedFieldOutput{ *this, std::move(lock) });
    }

    ~GLFieldOutput() override {
        glfwDestroyWindow(context);
    }
};




void updateState() {
    curTime = std::chrono::steady_clock::now();

    vec2d global = mouseToGlobal();
    vec2i cell = globalAsCell(global);

    const auto gridUpdateElapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(curTime - lastGridUpdateTime).count();
    if(
        gridUpdate
        && (gridUpdateElapsedTime >= 1000.0 / gridUpdatesPerSecond)
        && grid->tryFinishGeneration()
    ) {
        lastGridUpdateTime = curTime;
        Timer<> t{};
        fieldUpdateWait.add(t.elapsedTime());
        isBufferSecond = !isBufferSecond;
        grid->startNewGeneration();
    }

    if (paintMode != PaintMode::NONE) {
        std::vector<Cell> cells{};
        auto const side = brushSize+1 - -brushSize;
        assert(side >= 0);
        auto const size = size_t(side * side);
        cells.reserve(size);

        for (int32_t yo = -brushSize; yo <= brushSize; yo++) {
            for (int32_t xo = -brushSize; xo <= brushSize; xo++) {
                const vec2i offset{ xo, yo };
                const auto coord = cell + offset;

                cells.push_back(
                    Cell{
                        paintMode == PaintMode::PAINT ? fieldCell::cellAlive : fieldCell::cellDead,
                        grid->coordAsIndex(cell + offset)
                    }
                );
                /*if (paintMode == PaintMode::PAINT) {
                    grid->setCellAtCoord(cell + offset, FieldCell::ALIVE);
                }
                else if (paintMode == PaintMode::DELETE) {
                    grid->setCellAtCoord(coord + offset, FieldCell::DEAD);
                }*/
            }
        }
        assert(cells.size() == size);

        grid->setCells(&cells[0], size);
    }

    auto const vpSizeDesired = getVpSizeDesired();

    const auto screenUpdateElapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(curTime - lastScreenUpdateTime).count();
    lastScreenUpdateTime = curTime;
    if (pan) {
        const vec2d  dmouse = applyLensDistortion(mousePos);
        const vec2d dpmouse = applyLensDistortion(pmousePos);
        auto const diff = (dmouse - dpmouse) * winSize.normalizeCoeff;
        space.vpPosDesired -= diff * space.vpSize;
        desiredZoomPoint -= (dmouse - dpmouse);
        dVpPosDesired -= (dmouse - dpmouse);
    }
    if (isZoomChanged) {
        auto const pivot = mouseToGlobal();
        space.vpPosDesired = pivot + (space.vpPos - pivot) / space.vpSize * vpSizeDesired;

        // reset current zoom position to avoid bleed from previous zoom position which can be far away.
        auto diff = zoomPoint - desiredZoomPoint;
        auto resetCurrent = diff.dot(diff) < 1;

        desiredZoomPoint = mousePos;
        if(resetCurrent) zoomPoint = desiredZoomPoint;

        isZoomChanged = false;
    }

    auto const coeff = 1 - std::pow(
        1 - chasingSpeed,
        screenUpdateElapsedTime / screenUpdateTimeMs
    );

    const auto prevSize = space.vpSize;
    const auto prevPos = space.vpPos;

    space.vpSize = misc::lerp(space.vpSize, vpSizeDesired, coeff);
    space.vpPos = misc::vec2lerp(space.vpPos, space.vpPosDesired, coeff);
    zoomPoint = misc::vec2lerp(zoomPoint, desiredZoomPoint, coeff);

    dVpSizeDesired += (space.vpSize - prevSize) / space.vpSize;

    dVpSize = misc::lerp(dVpSize, dVpSizeDesired, coeff);
    dVpPos = misc::vec2lerp(dVpPos, dVpPosDesired, coeff);

    dVpSizeDesired = misc::lerp(dVpSizeDesired, 0.0, coeff);
    dVpPosDesired = misc::vec2lerp(dVpPosDesired, vec2d(0.0), coeff);

    pmousePos = mousePos;
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

void APIENTRY glDebugOutput(GLenum source, GLenum type, unsigned int id, GLenum severity,
                            GLsizei length, const char *message, const void *userParam) {
    if(id == 131169 || id == 131185 || id == 131218 || id == 131204) return;
    std::cout << message << '\n';
}

int main() {
    if(!glfwInit()) return -1;

    GLFWwindow* window;
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);

#define FULLSCREEN
#ifdef FULLSCREEN
    window = glfwCreateWindow(1920, 1080, "Game ofLife", glfwGetPrimaryMonitor(), NULL);
#else
    window = glfwCreateWindow(800, 800, "Game ofLife", NULL, NULL);
#endif

    if(!window) {
        glfwTerminate();
        return 1;
    }

    int width, height;
    glfwGetWindowSize(window, &width, &height);
    updateWinSize(width, height);
    updateSpace();

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    GLenum err = glewInit();
    if(err != GLEW_OK) {
        std::cerr << "Error:\n" << glewGetErrorString(err) << '\n';
        glfwTerminate();
        return 2;
    }

    int flags;
    glGetIntegerv(GL_CONTEXT_FLAGS, &flags);
    if (flags & GL_CONTEXT_FLAG_DEBUG_BIT) {
        glEnable(GL_DEBUG_OUTPUT);
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
        glDebugMessageCallback(glDebugOutput, nullptr);
        glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE);
    }

    //callbacks
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);

    GLuint mainProg = glCreateProgram();
    ShaderLoader sl{};
    sl.addScreenSizeTriangleStripVertexShader();
    sl.addShaderFromProjectFilePath("shaders/fs.shader", GL_FRAGMENT_SHADER, "Main shader");

    sl.attachShaders(mainProg);

    glLinkProgram(mainProg);
    glValidateProgram(mainProg);

    sl.deleteShaders();

    glUseProgram(mainProg);

    siv::PerlinNoise noise{};
    std::srand(2474941956);

    glGenBuffers(1, &packedGrid1);
    glGenBuffers(1, &packedGrid2);

    if(GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        status != GL_FRAMEBUFFER_COMPLETE
    ) {
        std::cerr << "framebuffer error: " << status << '\n';
        return 2;
    }

    auto const current_outputs = [window]() -> std::unique_ptr<FieldOutput> {
        return std::unique_ptr<FieldOutput>( new GLFieldOutput{ false, window } );
    };
    auto const buffer_outputs = [window]() -> std::unique_ptr<FieldOutput> {
        return std::unique_ptr<FieldOutput>(new GLFieldOutput{ true, window });
    };

    const auto currrrr = current_outputs();
    const auto bufffff = buffer_outputs();

    grid = std::unique_ptr<Field>{ new Field(
        gridWidth, gridHeight, numberOfTasks,
        current_outputs, buffer_outputs
    ) };

    field_size_bytes = grid->size_bytes();

    //{
    //    AutoTimer<> t{ "set" };
    //    for (size_t i = 0; i < gridSize; i++) {
    //        const double freq = 0.05;
    //        const auto coord = grid->indexAsCoord(i);
    //        //const auto a = misc::map<double>(noise.accumulatedOctaveNoise2D_0_1(coord.x * freq, coord.y * freq, 2), 0, 1, -0.4, 0.9);
    //        if ((std::rand() / (double)0x7fff) /*+ a*/ > 0.6) grid->setCellAtIndex(i, FieldCell::ALIVE);
    //    }
    //}

    {
        auto* const rawData = grid->rawData();
        for (size_t i = 0; i < field_size_bytes / 4; i++) {
            uint32_t cells{ 0 };
            for (unsigned j{ 0 }; j < 2; ++j) {
                cells = (cells << 16) | std::rand() % 32767u;
                static_assert(32767u == 0x7fff, "");
            }
            rawData[i] = cells;
        }
    }

    gridUpdate = false;
    grid->startCurGeneration();

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, packedGrid1);
    glBufferData(GL_SHADER_STORAGE_BUFFER, misc::roundUpIntTo(field_size_bytes, 4), NULL, GL_DYNAMIC_DRAW);
    //glBufferData(GL_SHADER_STORAGE_BUFFER, misc::roundUpIntTo(field_size_bytes * 2, 4), NULL, GL_DYNAMIC_DRAW);
    //if(isWritingSecondBuffer) glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, field_size_bytes, grid->rawData());
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, packedGrid1);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, packedGrid2);
    glBufferData(GL_SHADER_STORAGE_BUFFER, misc::roundUpIntTo(field_size_bytes, 4), NULL, GL_DYNAMIC_DRAW);
    //if (!isWritingSecondBuffer) glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, field_size_bytes, grid->rawData());
    //glBufferSubData(GL_SHADER_STORAGE_BUFFER, !isWritingSecondBuffer * field_size_bytes, field_size_bytes, grid->rawData());
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, packedGrid2);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    currrrr->write({ 0, misc::intDivCeil(field_size_bytes, 4), grid->rawData() });

    //uniforms set
    glUniform2i(glGetUniformLocation(mainProg, "winSize"), winSize.windowSize.x, winSize.windowSize.y);

    glUniform1i(glGetUniformLocation(mainProg, "gridWidth"), gridWidth);
    glUniform1i(glGetUniformLocation(mainProg, "gridHeight"), gridHeight);
    glUniform1ui(glGetUniformLocation(mainProg, "gridWidth_actual"), grid->width_actual());

    GLint vpSizeP = glGetUniformLocation(mainProg, "vpSize");
    GLint vpPosP = glGetUniformLocation(mainProg, "vpPos");

    GLint mousePosP = glGetUniformLocation(mainProg, "mousePos");
    GLint zoomPointP = glGetUniformLocation(mainProg, "zoomPoint");

    glUniform1f(glGetUniformLocation(mainProg, "lensDistortion"), lensDistortion);

    GLint is2ndBufferP = glGetUniformLocation(mainProg, "is2ndBuffer");
    glUniform1ui(glGetUniformLocation(mainProg, "bufferOffset_bytes"), grid->size_bytes());

    GLint mDeltaScaleChangeP = glGetUniformLocation(mainProg, "deltaScaleChange");


    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &frameBufferTexture);
    glBindTexture(GL_TEXTURE_2D, frameBufferTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, winSize.windowSize.x, winSize.windowSize.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
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

    GLint ppDeltaSizeChangeP = glGetUniformLocation(postProcessingProg, "deltaSizeChange");
    GLint deltaOffsetChangeP = glGetUniformLocation(postProcessingProg, "deltaOffsetChange");

    glUniform1f(glGetUniformLocation(postProcessingProg, "lensDistortion"), lensDistortion);
    GLint ppZoomPointP = glGetUniformLocation(postProcessingProg, "zoomPoint");

    lastGridUpdateTime = lastScreenUpdateTime = curTime = std::chrono::steady_clock::now();

    auto tim = std::chrono::steady_clock::now();

    while (!glfwWindowShouldClose(window))
    {
        glUseProgram(mainProg);

        Timer<> frame{};
        {

            Timer<> t{};
            glUniform2f(vpPosP, space.vpPos.x, space.vpPos.y);
            glUniform1f(vpSizeP, space.vpSize);

            glUniform2f(mousePosP, mousePos.x, mousePos.y);
            glUniform2f(zoomPointP, zoomPoint.x, zoomPoint.y);

            glUniform1ui(is2ndBufferP, !isBufferSecond);

            glUniform1f(mDeltaScaleChangeP, dVpSize);

            glFinish();

            set.add(t.elapsedTime());
        }

        while ((err = glGetError()) != GL_NO_ERROR)
        {
            fprintf(stderr, "Error %u\n", err);
        }

        //glClear(GL_COLOR_BUFFER_BIT);
        {
            //std::unique_lock<std::mutex> lock{ gpuBufferLock };
            {
                Timer<> t{};

                glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
                glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
                glBindFramebuffer(GL_FRAMEBUFFER, 0);
                draw.add(t.elapsedTime());
            }

            {
                auto const windowSize = winSize.windowSizeD;
                Timer<> t{};
                glUseProgram(postProcessingProg);
                glBindTexture(GL_TEXTURE_2D, frameBufferTexture);

                glUniform1i(frameBufferP, 0);
                glUniform2f(textureSizeP, windowSize.x, windowSize.y);

                glUniform1f(ppDeltaSizeChangeP, dVpSize);
                glUniform2f(deltaOffsetChangeP, dVpPos.x / windowSize.x, -dVpPos.y / windowSize.y);
                glUniform2f(ppZoomPointP, zoomPoint.x / windowSize.x, 1 - zoomPoint.y / windowSize.y);

                glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
                glBindTexture(GL_TEXTURE_2D, 0);
                postProcessing.add(t.elapsedTime());
            }

            {
                Timer<> t{};
                glfwSwapBuffers(window);
                swap.add(t.elapsedTime());
            }
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
