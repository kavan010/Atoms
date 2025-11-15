#version 330 core
out vec4 FragColor;
in vec2 TexCoord;
uniform sampler2D uTex;
uniform vec2 uDir; // (1/width,0) or (0,1/height)

void main(){
    float kernel[5] = float[](0.204164, 0.304005, 0.093913, 0.010381, 0.000558); // example weights (symmetric)
    vec2 uv = TexCoord;
    vec4 sum = texture(uTex, uv) * kernel[0];
    for(int i=1;i<5;i++){
        sum += texture(uTex, uv + uDir * float(i)) * kernel[i];
        sum += texture(uTex, uv - uDir * float(i)) * kernel[i];
    }
    FragColor = sum;
}
