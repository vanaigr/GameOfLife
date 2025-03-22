#version 430

uniform float deltaScaleChange;
//uniform vec2 deltaOffsetChange;

uniform vec2 vpPos;
uniform float vpSize;
uniform ivec2 winSize;

uniform float lensDistortion;

uniform int gridWidth;
uniform int gridHeight;
uniform uint gridWidth_actual;

uniform vec2 mousePos, zoomPoint;

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
    return (coord / winSize.y - 0.5) * vpSize + vpPos;
}

ivec2 globalAsCell(const vec2 coord) {
    return ivec2(
        int(mod(coord.x, gridWidth )),
        int(mod(coord.y, gridHeight))
    );
}

vec2 applyLensDistortion(const vec2 screenCoord, const float intensity) {
    const vec2 size2 = vec2(winSize) * 0.5;
    const vec2 centerCoord = screenCoord - size2;
    const float size = dot(centerCoord, centerCoord);
    const float maxSize = dot(size2, size2);
    const float coeff = sqrt(size / maxSize);
    return screenCoord - centerCoord * coeff * intensity;
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
    //if (cellX == 0 || cellX == gridWidth - 1 || cellY == 0 || cellY == gridHeight - 1) edgeMask = .5;

    const uint mask = cellAt(index);

    const bool isCell = mask == 1;
    const bool isNothing = mask == 0;


    const float dx = mod(global.x, 1.0);
    const float dy = mod(global.y, 1.0);
    const float padding = 1 / 15.0;
    const bool isPadding = false; //(dx < padding || dx > 1 - padding) || (dy < padding || dy > 1 - padding);


    const vec3 cell = float(isCell) * float(!isPadding) * cellColor;
    const vec3 bkg = float(isNothing) * bkgColor;
    const vec3 cellPadding = float(isCell) * float(isPadding) * bkgColor;
    vec3 col = cell + bkg + cellPadding;
    col = mix(col, vec3(.5, .5, .5), edgeMask);
    col = mix(col, vec3(.5, .5, .5), float(isPadding));
    return col;
}

vec3 col(const vec2 coord) {
    const vec2 windowCenter = vec2(winSize) * 0.5;
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
