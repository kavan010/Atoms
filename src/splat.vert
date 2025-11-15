#version 330 core
layout(location=0) in vec3 inPosition; // world position
uniform mat4 view;
uniform mat4 projection;
uniform float pointScreenSize; // size in pixels
void main(){
    vec4 clip = projection * view * vec4(inPosition, 1.0);
    gl_Position = clip;
    // convert desired pixel size to gl_PointSize using clip.w
    gl_PointSize = pointScreenSize * (1.0 / clip.w); // approximate; tune as needed
}
