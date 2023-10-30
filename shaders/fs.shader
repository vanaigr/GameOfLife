#version 430

#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif

uniform double size;
uniform dvec2 pos;

uniform float cellSize_px;


uniform float deltaScaleChange;
//uniform vec2 deltaOffsetChange;

uniform int width;
uniform int height;

uniform float lensDistortion;

uniform int gridWidth;
uniform int gridHeight;
uniform uint gridWidth_actual;

uniform float r1;
uniform float r2;

uniform vec2 mousePos, zoomPoint;

vec2 windowSize() {
    return vec2(width, height);
}

uniform uint is2ndBuffer;
uniform uint bufferOffset_bytes;

layout(std430, binding = 1) buffer Grid1
{
    uint grid[200];
} packedGrid1;

layout(std430, binding = 2) buffer Grid2
{
    uint grid[200];
} packedGrid2;

uint cellAt(const uint index) {
    const uint row = index / gridWidth;
    const uint col = index % gridWidth;

    const uint index_actual = row * gridWidth_actual + col;
    const uint index_actual_buffer = index_actual;// +is2ndBuffer * bufferOffset_bytes * 8;

    const uint arrIndex = index_actual_buffer / 32;
    const uint arrShift = (index_actual_buffer % 32);
    if(is2ndBuffer != 0)
        return ((packedGrid2.grid[arrIndex]) >> arrShift) & 1;
    return ((packedGrid1.grid[arrIndex]) >> arrShift) & 1;
}

layout(origin_upper_left) in vec4 gl_FragCoord;
out vec4 color;

vec2 distortedScreenToGlobal(const vec2 coord) {
    const vec2 wh2 = vec2(width, height)*0.5;
    return vec2(((coord - wh2) * size + wh2 + pos) / cellSize_px);
}

ivec2 globalAsCell(const vec2 coord) {
    const float cx = coord.x;
    const float cy = coord.y;
    const int cellX = int(mod(cx, gridWidth));
    const int cellY = int(mod(cy, gridHeight));
    return ivec2(cellX, cellY);
}

vec2 applyLensDistortion(const vec2 coord, const float intensity) {
    const float x = coord.x, y = coord.y;
    const float w2 = width / 2.0, h2 = height / 2.0;
    const float xc = x - w2, yc = y - h2;
    const float dist_sq = (xc * xc + yc * yc);
    const float maxDist_sq = (w2 * w2 + h2 * h2);
    const float distortion = sqrt(dist_sq / maxDist_sq);
    const float newX = x - distortion * intensity * xc;
    const float newY = y - distortion * intensity * yc;
    return vec2(newX, newY);
}

const vec3 wallColor = vec3(30.0 / 255.0, 240.0 / 255.0, 20.0 / 255.0);
const vec3 cellColor = vec3(30.0 / 255.0, 30.0 / 255.0, 30.0 / 255.0);
const vec3 bkgColor = vec3(225.0 / 255.0, 225.0 / 255.0, 225.0 / 255.0);

vec3 colorForCoords(const vec2 screenCoords) {
    const vec2 global = distortedScreenToGlobal(screenCoords);
    const ivec2 cellInt = globalAsCell(global);
    const int cellX = cellInt.x;
    const int cellY = cellInt.y;

    const uint index = cellX + gridWidth * cellY;

    float edgeMask = 0;
    if (cellX == 0 || cellX == gridWidth - 1 || cellY == 0 || cellY == gridHeight - 1) edgeMask = .5;

    const uint mask = cellAt(index);

    const bool isCell = mask == 1;
    const bool isNothing = mask == 0;


    const float dx = mod(global.x, 1.0);
    const float dy = mod(global.y, 1.0);
    const float padding = 1 / 15.0;
    const bool isPadding = (dx < padding || dx > 1 - padding) || (dy < padding || dy > 1 - padding);
    

    const vec3 cell = float(isCell) * float(!isPadding) * cellColor;
    const vec3 bkg = float(isNothing) * bkgColor;
    const vec3 cellPadding = float(isCell) * float(isPadding) * bkgColor;
    vec3 col = cell + bkg + cellPadding;
    col = mix(col, vec3(.5, .5, .5), edgeMask);
    col = mix(col, vec3(.5, .5, .5), float(isPadding));
    return col;
}

vec3 colorForCoordsChromaticAbberation(const vec2 coord, const float caIntens) {
    const vec2 newCoord = applyLensDistortion(coord, caIntens);
    return colorForCoords(newCoord);
}

vec3 col(const vec2 coord) {
    const vec2 windowCenter = vec2(width, height)*0.5;
    const vec2 distortedCoord =  applyLensDistortion(
        applyLensDistortion(coord, lensDistortion) - zoomPoint + windowCenter,
        -deltaScaleChange / 3
    ) + zoomPoint - windowCenter;

    return colorForCoords(distortedCoord);
}

float rand(const vec2 co) {
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 sampleN(const vec2 coord, const uint n) {
    const vec2 pixelCoord = floor(coord);
    const float fn = float(n);

    vec3 result = vec3(0, 0, 0);
    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++) {
            const vec2 offset = vec2(rand(coord.xy), rand(coord.xy + pixelCoord)) / fn;
            const vec2 coord = pixelCoord + vec2(i / fn, j / fn) + offset;

            const vec3 sampl = col(coord);
            result += sampl;
        }
    }

    return result / (fn * fn);
}

void main(void) {
    const vec2 coord = gl_FragCoord.xy;
    color = vec4(sampleN(coord, 3), 1);
}
