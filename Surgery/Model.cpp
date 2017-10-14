#include "Model.h"



Model::Model()
{
}

Model::Model(char* filename, GLenum usage)
{
	LoadFromFile(filename);
	SetUpGLBuffer(usage);
}

Model::~Model()
{
	DeleteGLBuffer();
}

void Model::LoadFromFile(char* filename)
{
	readObj(filename, this->vertices, this->indices);
}

void Model::SetUpGLBuffer(GLenum usage)
{
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, this->vertices.size() * sizeof(Vert), this->vertices.data(), usage);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(VAO);
	{
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		{
			//position
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (GLvoid*)0);

			//normal
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (GLvoid*)offsetof(Vert, normal));

			//uv
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vert), (GLvoid*)offsetof(Vert, uv));
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
	glBindVertexArray(0);
}

void Model::DeleteGLBuffer()
{
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);
}

/*void SetUpGL_Model(Model& model, GLuint& _VAO, GLuint& _VBO, GLenum usage)
{
	glGenVertexArrays(1, &_VAO);
	glGenBuffers(1, &_VBO);

	glBindBuffer(GL_ARRAY_BUFFER, _VBO);
	glBufferData(GL_ARRAY_BUFFER, model.vertices.size() * sizeof(Vert), model.vertices.data(), usage);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(_VAO);
	{
		glBindBuffer(GL_ARRAY_BUFFER, _VBO);
		{
			//position
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (GLvoid*)0);

			//normal
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (GLvoid*)offsetof(Vert, normal));

			//uv
			glEnableVertexAttribArray(2);
			glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vert), (GLvoid*)offsetof(Vert, uv));
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
	glBindVertexArray(0);
}

void DrawGL_Model(Model& model, const GLuint _VAO, GLenum mode, bool useIndex)
{
	glBindVertexArray(_VAO);
	if (useIndex)
	{
		glDrawElements(mode, model.indices.size(), GL_UNSIGNED_INT, model.indices.data());
	}
	else
	{
		glDrawArrays(mode, 0, model.vertices.size());
	}
	glBindVertexArray(0);
}*/