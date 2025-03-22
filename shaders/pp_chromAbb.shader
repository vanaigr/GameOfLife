#version 430

#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif

uniform sampler2D frameBuffer;
uniform vec2 textureSize;

layout(origin_upper_left) in vec4 gl_FragCoord;
out vec4 color;

uniform float deltaSizeChange;
uniform vec2 deltaOffsetChange;
uniform double size;

uniform float lensDistortion;

uniform vec2 zoomPoint;

vec3 colorForCoords(const vec2 coord) {
    return texture(frameBuffer, vec2(coord.x, 1 - coord.y)).xyz;
}

vec2 applyLensDistortion(const vec2 coord, const float intensity) {
    const float x = coord.x, y = coord.y;
    const float w2 = 1 / 2.0, h2 = 1 / 2.0;
    const float xc = x - w2, yc = y - h2;
    const float dist = sqrt(xc * xc + yc * yc);
    const float maxDist = sqrt(w2 * w2 + h2 * h2);
    const float distortion = dist / maxDist;
    const float newX = x - distortion * intensity * (xc);
    const float newY = y - distortion * intensity * (yc);
    return vec2(newX, newY);
}


vec3 colorForCoordsChromaticAbberation(const dvec2 coord, const double caIntens, const vec2 deltaPosChange, const float scaleIntensity) {
    vec2 newCoord = applyLensDistortion(vec2(coord), float(caIntens));
    newCoord = applyLensDistortion(newCoord - applyLensDistortion(zoomPoint, float(caIntens)) + vec2(1, 1) / 2.0 + deltaPosChange, scaleIntensity) + applyLensDistortion(zoomPoint, float(caIntens)) - vec2(1, 1) / 2.0;
    return colorForCoords(newCoord);
}

vec3 col(const vec2 coord) {
    const vec2 distortedCoord = coord / textureSize;// applyLensDistortion(coord  / textureSize, lensDistortion + (deltaScaleChange / 5.0));
    //const vec3 empty = colorForCoordsChromaticAbberation(distortedCoord, 0);
    const float red = colorForCoordsChromaticAbberation(distortedCoord, 0.003, -deltaOffsetChange / 140, -deltaSizeChange / 60.0).r;
    const float green = colorForCoordsChromaticAbberation(distortedCoord, -0.003, deltaOffsetChange / 140, deltaSizeChange / 60.0).g;
    const float blue = colorForCoordsChromaticAbberation(distortedCoord, -0.003, deltaOffsetChange / 140, deltaSizeChange / 60.0).b;

    return vec3(red, green, blue);
    //return empty;
}

float rand(const vec2 co) {
    return fract(sin(dot(co + vec2(11111.0, 22222.0), vec2(12.9898, 78.233))) * 43758.5453);
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

float vignette(vec2 coord, float startDistance, float endDistance, float intensity) {
    float x = coord.x, y = coord.y;
    float w2 = textureSize.x / 2.0, h2 = textureSize.y / 2.0;
    float xc = x - w2, yc = y - h2;
    float dist = sqrt(xc * xc + yc * yc);
    float maxDist = sqrt(w2 * w2 + h2 * h2);
    float distortion = dist / maxDist;

    return smoothstep(startDistance, endDistance, distortion) * intensity;
}

void main(void) {
    const vec2 coord = gl_FragCoord.xy;

    color = vec4(sampleN(coord, 3), 1.0);

    float v = vignette(coord, .6, 1, .15);

    color = mix(color, vec4(0, 0, 0, 1), v);
}
