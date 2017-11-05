#version 330 core

layout(points) in;
layout(line_strip, max_vertices = 256) out;

uniform float radius;
uniform int lats;
uniform int longs;

in vec3 FragPos[]; 

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

			gl_Position = _center + vec4(p0, 0.0);
			EmitVertex();
			gl_Position = _center + vec4(p1, 0.0);
			EmitVertex();
			gl_Position = _center + vec4(p2, 0.0);
			EmitVertex();
			gl_Position = _center + vec4(p3, 0.0);
			EmitVertex();
		}
	}
}

void BuildSphere2(vec4 _center, int _lats, int _longs)
{
	for(int i = 0; i < _lats; i++)
	{
		float lat0 = PI * (-0.5 + (i - 1) / _lats);
		float z0 = sin(lat0);
		float zr0 = cos(lat0);

		float lat1 = PI * (-0.5 + i / _lats);
		float z1 = sin(lat1);
		float zr1 = cos(lat1);

		for(int j = 0; j < _longs; j++)
		{
			float lng = 2 * PI * (j - 1) / _longs;
			float x = cos(lng);
			float y = sin(lng);

			vec3 pos1 = vec3(x * zr0, y * zr0, z0);
			gl_Position = vec4(_center.xyz + pos1, _center.w);
			EmitVertex();

			vec3 pos2 = vec3(x * zr1, y * zr1, z1);
			gl_Position = vec4(_center.xyz + pos2, _center.w);
			EmitVertex();
		}
	}
}

void main(void)
{
	BuildSphere(gl_in[0].gl_Position, radius, lats, longs);
	//BuildSphere(vec4(FragPos[0], 1.0), radius, lats, longs);
	//BuildSphere2(gl_in[0].gl_Position, lats, longs);
	EndPrimitive();
}