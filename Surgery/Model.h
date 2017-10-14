#pragma once
#include "Utilities.h"

class Model
{
public:
	Model();
	Model(char* filename, GLenum usage = GL_STATIC_DRAW);
	~Model();

	void LoadFromFile(char* filename);
	void SetUpGLBuffer(GLenum usage = GL_STATIC_DRAW);
	void DeleteGLBuffer();

public:
	std::vector<Vert> vertices;
	std::vector<int> indices;

	GLuint VAO;
	GLuint VBO;
};

//void SetUpGL_Model(Model& model, GLuint& _VAO, GLuint& _VBO, GLenum usage = GL_STATIC_DRAW);
//void DrawGL_Model(Model& model, const GLuint _VAO, GLenum mode, bool useIndex = true);