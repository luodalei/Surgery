#version 330 core
layout (location = 0) in vec3 aPos;

out vec3 FragPos;

uniform mat4 modelView;
uniform mat4 projection;

void main()
{
    FragPos = vec3(modelView * vec4(aPos, 1.0));
    
    gl_Position = projection * vec4(FragPos, 1.0);
}