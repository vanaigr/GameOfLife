#version 430

#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif

uniform sampler2D frameBuffer;
uniform vec2 textureSize;

in vec4 gl_FragCoord;
out vec4 color;

uniform float deltaScaleChange;
uniform vec2 deltaOffsetChange;
uniform float currentScale;

vec4 colorForCoords(vec2 coords) {
    return texture2D(frameBuffer, coords.xy);
}

vec2 applyLensDistortion(vec2 coord, float intensity) {
    float x = coord.x, y = coord.y;
    float w2 = 1 / 2.0, h2 = 1 / 2.0;
    float xc = x - w2, yc = y - h2;
    float dist = sqrt(xc * xc + yc * yc);
    float maxDist = sqrt(w2 * w2 + h2 * h2);
    float distortion = dist / maxDist;
    float newX = x - distortion * intensity * (xc);
    float newY = y - distortion * intensity * (yc);
    return vec2(newX, newY);
}


vec4 colorForCoordsChromaticAbberation(vec2 coord, float caIntens) {
    vec2 newCoord = applyLensDistortion(coord, caIntens);
    return colorForCoords(newCoord);
}

void main(void) {
    vec2 coord = gl_FragCoord.xy / textureSize;

    float red = colorForCoordsChromaticAbberation(coord - deltaOffsetChange.x / textureSize.x / 32.0 * currentScale, -0.002 + deltaScaleChange / 10.0).r;
    float green = colorForCoordsChromaticAbberation(coord, 0).g;
    float blue = colorForCoordsChromaticAbberation(coord + deltaOffsetChange.y / textureSize.y / 32.0 * currentScale, 0.002 - deltaScaleChange / 10.0).b;

    color = vec4(red, green, blue, 1.0);
}