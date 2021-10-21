#version 430

#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif

uniform sampler2D frameBuffer;
uniform vec2 textureSize;

in vec4 gl_FragCoord;
out vec4 color;

float vignette(vec2 coord, float startDistance, float endDistance, float intensity) {
    float x = coord.x, y = coord.y;
    float w2 = 1.0 / 2.0, h2 = 1.0 / 2.0;
    float xc = x - w2, yc = y - h2;
    float dist = sqrt(xc * xc + yc * yc);
    float maxDist = sqrt(w2 * w2 + h2 * h2);
    float distortion = dist / maxDist;

    return smoothstep(startDistance, endDistance, distortion) * intensity;
}

float distanceToEdges_px(float coord, float size) {
    float size2 = size / 2.0;
    return size2 - abs((coord - size2));
}

/*float squareVignette(vec2 coord, float startDistance, float endDistance, float startGrad, float endGrad, float intensity) {
    float w2 = width / 2f, h2 = height / 2f;
    float whMin2 = min(width, height) / 2.0;
    float start = whMin2 * (1 - startDistance);
    float end = whMin2 * (1 - endDistance);
    float diff = end - start;
    float x = coord.x, y = coord.y;
    float dx = distanceToEdges_px(x, width);
    float dy = distanceToEdges_px(y, height);
    float xc = w2 - dx;
    float yc = h2 - dy;
    float xd = max(xc - (w2 - start), 0);
    float yd = max(yc - (h2 - start), 0);
    float d = sqrt(xd * xd + yd * yd);
    float maxD = sqrt(diff * diff + diff * diff);
    float distortion = d / maxD;

    return smoothstep(startGrad, endGrad, distortion) * intensity;
}*/

float squareVignette(vec2 coord, float radius, float padding, float startGrad, float endGrad) {
    float width = textureSize.x, height = textureSize.y;
    float w2 = width / 2f, h2 = height / 2f;
    float whMin2 = min(width, height) / 2.0;
    float x = coord.x, y = coord.y;

    float p_px = padding * whMin2;
    float r_px = radius * (whMin2 - p_px);

    float dx = distanceToEdges_px(x, width) - p_px;
    float dy = distanceToEdges_px(y, height) - p_px;

    float xd = max(r_px - dx, 0);
    float yd = max(r_px - dy, 0);

    float d = sqrt(xd * xd + yd * yd);
    float maxD = r_px;

    float effect = d / maxD;
    return smoothstep(startGrad, endGrad, effect);
}

float edgeVignette(vec2 coord, float startDistance, float endDistance, float intensity) {
    float width = textureSize.x, height = textureSize.y;

    float whMin = min(width, height);
    float whMin2 = whMin / 2.0;
    float distToStartHalf_px = (1 - startDistance) * whMin2;
    float distToEndHalf_px = (1 - endDistance) * whMin2;
    float distToXBorderHalf_px = distanceToEdges_px(coord.x, width);
    float distToYBorderHalf_px = distanceToEdges_px(coord.y, height);
    float xV = smoothstep(distToStartHalf_px, distToEndHalf_px, distToXBorderHalf_px);
    float yV = smoothstep(distToStartHalf_px, distToEndHalf_px, distToYBorderHalf_px);
    return max(xV, yV) * intensity;
}

void main() {
    vec2 coord = gl_FragCoord.xy;

    float v = squareVignette(coord, 0.06, 0.003, .9, 1.1);

    color = texture2D(frameBuffer, coord.xy / textureSize.xy) * (1.0 - v);
}