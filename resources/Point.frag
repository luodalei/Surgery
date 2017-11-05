#version 120

// interpolated vertex position in eye coordinate space (from vertex shader)
//varying vec4 vPosition;
varying vec4 vColor;
//in vec4 vColor;

void main (void)
{
	gl_FragColor = vColor;
}