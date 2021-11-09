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

vec3 col(vec2 coord) {
    vec2 distortedCoord = coord;// applyLensDistortion(coord, lensDistortion + (deltaScaleChange / 5.0));

    float red = colorForCoordsChromaticAbberation(distortedCoord, 0).r;
    float green = colorForCoordsChromaticAbberation(distortedCoord + deltaOffsetChange / 13.0 / size, -0.004 + deltaScaleChange / 6.0 / size).g;
    float blue = colorForCoordsChromaticAbberation(distortedCoord + deltaOffsetChange / 13.0 / size, -0.004 + deltaScaleChange / 6.0 / size).b;

    return vec3(red, green, blue);
}

float rand(vec2 co) {
    return fract(sin(dot(co + vec2(11111.0, 22222.0), vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 sampleN(vec2 coord, uint n) {
    vec2 pixelCoord = floor(coord);
    float fn = float(n);

    vec3 result = vec3(0, 0, 0);
    for (uint i = 1; i <= n; i++) {
        for (uint j = 1; j <= n; j++) {
            vec2 coord = pixelCoord + vec2(i / fn, j / fn);
            vec2 offset = vec2(rand(coord.xy), rand(coord.yx)) / fn;
            coord -= offset;

            vec3 sampl = col(coord);
            result += sampl;
        }
    }

    return result / (fn * fn);
}

void main(void) {
    vec2 coord = gl_FragCoord.xy;

    color = vec4(sampleN(coord, 2), 1.0);
}