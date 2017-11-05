#version 120
#extension GL_EXT_geometry_shader4 : enable

//layout(points) in;
//layout(triangle_strip, max_vertices = 200) out;

uniform float radius;
uniform int lats;
uniform int longs;

varying vec4 vColor[];
varying vec4 vColorNew;

const float PI = 3.141592653;

vec3 GenPoint(float _u, float _v, float _radius)
{
	return vec3(cos(_u)*sin(_v)*_radius, cos(_v)*radius, sin(_u)*sin(_v)*_radius);
}

void BuildSphere(vec4 _center, float _radius, int _lats, int _longs)
{
	float startU = 0; //u is longitude
	float startV = 0; //v is latitude
	float endU = 2 * PI;
	float endV = PI;

	float stepU = (endU - startU) / _longs;
	float stepV = (endV - startV) / _lats;

	for(int i = 0; i < _longs; i++)
	{
		for(int j = 0; j < _lats; j++)
		{
			float u = i * stepU + startU;
			float v = j * stepV + startV;
			float un = (i + 1 == _longs) ? endU : (i + 1) * stepU + startU;
			float vn = (j + 1 == _lats) ? endV : (j + 1) * stepV + startV;

			vec3 p0 = GenPoint(u, v, _radius);
			vec3 p1 = GenPoint(u, vn, _radius);
			vec3 p2 = GenPoint(un, v, _radius);
			vec3 p3 = GenPoint(un, vn, _radius);

			gl_Position = vec4(_center.xyz + p0, _center.w);
			vColorNew = vColor[0];
			EmitVertex();
			gl_Position = vec4(_center.xyz + p1, _center.w);
			EmitVertex();
			gl_Position = vec4(_center.xyz + p2, _center.w);
			EmitVertex();
			gl_Position = vec4(_center.xyz + p3, _center.w);
			EmitVertex();
		}
	}
}

void main(void)
{
	//BuildSphere(gl_in[0].gl_Position, radius, lats, longs);
	BuildSphere(gl_PositionIn[0], radius, lats, longs);
	EndPrimitive();
}