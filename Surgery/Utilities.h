#pragma once
#include <glm\glm.hpp>
#include <glad\glad.h>
#include <string>
#include <vector>
#include <iostream>

//SOIL2
#include <SOIL2.h>

//assimp
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

using namespace std;

//texture
typedef struct _Texture
{
	GLuint id;
	string type;
	int width, height;
} Texture;

//vertex
typedef struct _Vert
{
	glm::vec3 pos;
	glm::vec3 normal;
	glm::vec2 uv;
} Vert;

//obj SOA
typedef struct _ObjInfo
{
	std::vector<glm::vec3> positions;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec2> uvs;
	std::vector<int> vIndices;
	std::vector<int> nIndices;
	std::vector<int> tIndices;
} ObjInfo;

//light
typedef struct _Light
{
	glm::vec3 lightPos = glm::vec3(0.0f);
	glm::vec3 lightAmbient = glm::vec3(0.2f);
	glm::vec3 lightDiffuse = glm::vec3(0.5f);
	glm::vec3 lightSpecular = glm::vec3(1.0f);
	float lightConstant = 1.0f;
	float lightLinear = 0.0f;
	float lightQuadratic = 0.0f;
} Light;


unsigned TextureFromFile(const char* filename, int &width, int &height, bool isRGBA = true);
Texture LoadTexture(string filename, string typeName, bool isRGBA = true);

void readObj(char *filename, ObjInfo& objInfo);
void readObj(char *filename, std::vector<Vert>& vertices, std::vector<int>& indices);
void processNode(aiNode* node, const aiScene* scene, std::vector<Vert>& vertices, std::vector<int>& indices);
void processMesh(aiMesh* mesh, const aiScene* scene, std::vector<Vert>& vertices, std::vector<int>& indices);

void readASC(char *filename, std::vector<glm::dvec3> &pointPos, int &pointCount);

void MatToArray(glm::mat4 mat, double* dst);
void ArrayToMat(double* src, glm::mat4& mat);