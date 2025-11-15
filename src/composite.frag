#version 330 core
out vec4 FragColor;
in vec2 TexCoord;
uniform sampler2D uDensity;
uniform float exposure; // controls overall brightness

vec3 palette(float v) {
    // simple mapping: blue -> cyan -> yellow -> red
    return mix(vec3(0.0,0.0,1.0), vec3(1.0,0.0,0.0), smoothstep(0.0, 1.0, v));
}

void main(){
    float d = texture(uDensity, TexCoord).r;
    float mapped = 1.0 - exp(-d * exposure); // tone-like mapping
    vec3 col = palette(mapped);
    float alpha = clamp(d * 5.0, 0.0, 1.0); // you control
    FragColor = vec4(col, alpha);
}
