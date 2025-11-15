#version 330 core
out vec4 FragColor; // we'll write to red channel
uniform float intensity; // per-particle multiplier
uniform float sigma; // controls gaussian shape

void main(){
    // gl_PointCoord ranges [0,1] across the point
    vec2 uv = gl_PointCoord - vec2(0.5);
    float d2 = dot(uv, uv);
    // gaussian
    float g = exp(-d2 / (2.0 * sigma * sigma));
    // write density into red channel
    FragColor = vec4(g * intensity, 0.0, 0.0, 1.0);
}
