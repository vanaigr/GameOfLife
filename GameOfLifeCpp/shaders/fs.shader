#version 430

#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif

uniform float currentScale;
uniform float size;
uniform float tx;
uniform float ty;

uniform float deltaScaleChange;
//uniform vec2 deltaOffsetChange;

uniform int width;
uniform int height;

uniform float lensDistortion;

uniform int gridWidth;
uniform int gridHeight;

uniform float r1;
uniform float r2;

uniform vec2 mousePos;

vec2 windowSize() {
    return vec2(width, height);
}

uniform uint is2ndBuffer;
uniform uint bufferOffset_bytes;

layout(std430, binding = 1) buffer Grid
{
    uint grid[];
} packedGrid;

uint cellAt(uint index) {
    uint bufIndex = index + is2ndBuffer * bufferOffset_bytes * 8;
    uint arrIndex = bufIndex / 32;
    uint arrShift = (bufIndex % 32);
    return ((packedGrid.grid[arrIndex]) >> arrShift) & 1;
}

layout(origin_upper_left) in vec4 gl_FragCoord;
out vec4 color;

vec2 distortedScreenToGlobal(vec2 coord) {
    return vec2(((coord.x - width / 2f) / currentScale + width / 2f - tx) / size, ((coord.y - height / 2f) / currentScale + height / 2f - ty) / size);
}

ivec2 globalAsCell(vec2 coord) {
    float cx = coord.x;
    float cy = coord.y;
    int cellX = int(mod(cx, gridWidth));
    int cellY = int(mod(cy, gridHeight));
    return ivec2(cellX, cellY);
}

vec2 applyLensDistortion(vec2 coord, float intensity) {
    float x = coord.x, y = coord.y;
    float w2 = width / 2f, h2 = height / 2f;
    float xc = x - w2, yc = y - h2;
    float dist = sqrt(xc * xc + yc * yc);
    float maxDist = sqrt(w2 * w2 + h2 * h2);
    float distortion = dist / maxDist;
    float newX = x - distortion * intensity * (xc);
    float newY = y - distortion * intensity * (yc);
    return vec2(newX, newY);
}

const vec4 wallColor = vec4(30.0 / 255.0, 240.0 / 255.0, 20.0 / 255.0, 1.0);
const vec4 cellColor = vec4(30.0 / 255.0, 30.0 / 255.0, 30.0 / 255.0, 1.0);
const vec4 bkgColor = vec4(225.0 / 255.0, 225.0 / 255.0, 225.0 / 255.0, 1.0);

vec4 colorForCoords(vec2 screenCoords) {
    vec2 global = distortedScreenToGlobal(screenCoords);
    ivec2 cellInt = globalAsCell(global);
    int cellX = cellInt.x;
    int cellY = cellInt.y;

    uint index = cellX + gridWidth * cellY;

    float edgeMask = 0;
    if (cellX == 0 || cellX == gridWidth - 1 || cellY == 0 || cellY == gridHeight - 1) edgeMask = .5;

    //uint arrIndex = index / (32 / 2);
    //uint arrShift = (index % (32 / 2)) * 2;

    //uint mask = ((packedGrid.grid[arrIndex]) >> arrShift) & 3;

    uint mask = cellAt(index);

    bool isCell = mask == 1;
    bool isNothing = mask == 0;

    bool isPadding;// = (dx < pv || dx > padding-pv) || (dy < pv || dy > padding-pv)
    {
        float dx = mod(global.x, 1.0);
        float dy = mod(global.y, 1.0);
        float padding = 1 / 15.0;
        isPadding = (dx < padding || dx > 1 - padding) || (dy < padding || dy > 1 - padding);
    }

    vec4 cell = float(isCell) * float(!isPadding) * cellColor;
    vec4 bkg = float(isNothing) * bkgColor;
    vec4 cellPadding = float(isCell) * float(isPadding) * bkgColor;
    vec4 col = cell + bkg + cellPadding;
    col = mix(col, vec4(.5, .5, .5, 1), edgeMask);
    col = mix(col, vec4(.5, .5, .5, 1), isPadding);
    return col;
}

vec4 colorForCoordsChromaticAbberation(vec2 coord, float caIntens) {
    vec2 newCoord = applyLensDistortion(coord, caIntens);
    return colorForCoords(newCoord);
}

vec4 col(vec2 coord) {
    //float v = vignette(coord, 0.97, 1, 1);

    vec2 distortedCoord = applyLensDistortion(coord, lensDistortion + (deltaScaleChange / 5));

    ////colorForCoords
    //float red = colorForCoordsChromaticAbberation(distortedCoord - deltaOffsetChange / 32.0 * currentScale, -0.002 + deltaScaleChange / 10.0).r;
    //float green = colorForCoordsChromaticAbberation(distortedCoord, 0).g;
    //float blue = colorForCoordsChromaticAbberation(distortedCoord + deltaOffsetChange / 32.0 * currentScale, 0.002 - deltaScaleChange / 10.0).b;

    //return vec4(red, green, blue, 1.0) * (1.0 - v);//colorForCoords(distortedCoord);

    return colorForCoords(distortedCoord);
}

vec4 sample2(vec2 coord) {
    vec4 c1 = col(coord + vec2(0.25, 0.25));
    vec4 c2 = col(coord + vec2(0.25, -0.25));
    vec4 c3 = col(coord + vec2(-0.25, 0.25));
    vec4 c4 = col(coord + vec2(-0.25, -0.25));

    return c1 * .25 + c2 * .25 + c3 * .25 + c4 * .25;
}

vec4 sampleN(vec2 coord, uint n) {
    int intn = int(n);
    int intn2 = intn / 2;
    float fn = float(n);
    vec4 result = vec4(.0, .0, .0, .0);
    for (int i = -intn2; i <= intn2; ++i) {
        for (int j = -intn2; j <= intn2; ++j) {
            vec4 sampl = col(coord + vec2(1.0 / (fn + 2) * i, 1.0 / (fn + 2) * j));
            result += sampl;
        }
    }
    return result / (fn * fn);
}

void main(void) {
    vec2 coord = gl_FragCoord.xy;

    color = sampleN(coord, 3);
}