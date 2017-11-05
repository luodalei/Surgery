#version 120
#extension GL_EXT_geometry_shader4 : enable

//layout(points) in;
//layout(points, max_vertices = 1) out;
varying vec4 vColor[];
varying vec4 vColorNew;

void main(void)
{
	gl_Position = gl_PositionIn[0];
	vColorNew = vColor[0];
	EmitVertex();
	EndPrimitive();
}