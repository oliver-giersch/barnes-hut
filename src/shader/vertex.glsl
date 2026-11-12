#version 450 core

layout(location=0) in vec3 pos;
uniform mat4 view;
uniform float r;
out vec3 vert_color;

void main() {
    gl_Position = view * vec4(pos, 1.0);
    float d = 2.0 * r;
    vert_color.r = (pos.x + r) / d;
    vert_color.g = (pos.y + r) / d;
    vert_color.b = (pos.z + r) / d;
}

