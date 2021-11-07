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
uniform double size;

uniform float lensDistortion;

vec4 colorForCoords(vec2 coord) {
    return texture2D(frameBuffer, coord / textureSize);
}

vec2 applyLensDistortion(vec2 coord, float intensity) {
    float x = coord.x, y = coord.y;
    float w2 = textureSize.x / 2.0, h2 = textureSize.y / 2.0;
    float xc = x - w2, yc = y - h2;
    float dist = sqrt(xc * xc + yc * yc);
    float maxDist = sqrt(w2 * w2 + h2 * h2);
    float distortion = dist / maxDist;
    float newX = x - distortion * intensity * (xc);
    float newY = y - distortion * intensity * (yc);
    return vec2(newX, newY);
}


vec4 colorForCoordsChromaticAbberation(dvec2 coord, double caIntens) {
    vec2 newCoord = applyLensDistortion(vec2(coord), float(caIntens));
    return colorForCoords(newCoord);
}

vec4 col(vec2 coord) {
    vec2 distortedCoord = coord;// applyLensDistortion(coord, lensDistortion + (deltaScaleChange / 5.0));

    float red = colorForCoordsChromaticAbberation(distortedCoord, 0).r;
    float green = colorForCoordsChromaticAbberation(distortedCoord - deltaOffsetChange / 13.0 / size, -0.004 - deltaScaleChange / 4.0 / size).g;
    float blue = colorForCoordsChromaticAbberation(distortedCoord - deltaOffsetChange / 13.0 / size, -0.004 - deltaScaleChange / 4.0 / size).b;

    return vec4(red, green, blue, 1.0);
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