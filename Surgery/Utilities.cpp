#include "Utilities.h"
#include <fstream>

/*unsigned TextureFromFile(const char* filename, int &width, int &height, bool isRGBA)
{
	GLuint textureID;
	glGenTextures(1, &textureID);


	unsigned char* image = SOIL_load_image(filename, &width, &height, 0, SOIL_LOAD_RGBA);
	if (image == 0)
	{
		std::cout << "Failed to load texture" << std::endl;
	}
	else
	{
		//assign texture to id
		glBindTexture(GL_TEXTURE_2D, textureID);
		glTexImage2D(GL_TEXTURE_2D, 0, isRGBA ? GL_RGBA : GL_RGB, width, height, 0, isRGBA ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, image);
		glGenerateMipmap(GL_TEXTURE_2D);

		// Parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glBindTexture(GL_TEXTURE_2D, 0);
	}

	SOIL_free_image_data(image);

	return textureID;
}

Texture LoadTexture(string filename, string typeName, bool isRGBA)
{
	Texture tex;
	int imgWidth, imgHeight;
	tex.id = TextureFromFile(filename.c_str(), imgWidth, imgHeight, isRGBA);
	tex.type = typeName;
	tex.width = imgWidth;
	tex.height = imgHeight;

	return tex;
}*/

void readObj(char *filename, ObjInfo& objInfo)
{
	FILE *file;

	file = fopen(filename, "r");
	if (!file)
	{
		fprintf(stderr, "readSurface() failed: can't open data file \"%s\".\n",
			filename);
		getchar();
		exit(1);
	}

	char buf[128];
	float vx, vy, vz;
	while (fscanf(file, "%s", buf) != EOF)
	{
		switch (buf[0])
		{
		case '#':
			fgets(buf, sizeof(buf), file);
			break;
		case 'v':
			switch (buf[1])
			{
			case '\0':
				fscanf(file, "%f %f %f", &vx, &vy, &vz);
				objInfo.positions.push_back(glm::vec3(vx, vy, vz));
				break;
			case 'n':
				fscanf(file, "%f %f %f", &vx, &vy, &vz);
				objInfo.normals.push_back(glm::vec3(vx, vy, vz));
				break;
			case 't':
				fscanf(file, "%f %f", &vx, &vy);
				objInfo.uvs.push_back(glm::vec2(vx, vy));
				break;
			default:
				break;
			}
			break;
		case 'f':
			switch (buf[1])
			{
			case '\0':
				int v1, t1, n1;
				int v2, t2, n2;
				int v3, t3, n3;
				int v4, t4, n4;
				fscanf(file, "%d/%d/%d", &v1, &t1, &n1);
				fscanf(file, "%d/%d/%d", &v2, &t2, &n2);
				fscanf(file, "%d/%d/%d", &v3, &t3, &n3);
				objInfo.vIndices.push_back(v1); objInfo.tIndices.push_back(t1); objInfo.nIndices.push_back(n1);
				objInfo.vIndices.push_back(v2); objInfo.tIndices.push_back(t2); objInfo.nIndices.push_back(n2);
				objInfo.vIndices.push_back(v3); objInfo.tIndices.push_back(t3); objInfo.nIndices.push_back(n3);

				if (fscanf(file, "%d/%d/%d", &v4, &t4, &n4) == 3)
				{
					objInfo.vIndices.push_back(v2); objInfo.tIndices.push_back(t2); objInfo.nIndices.push_back(n2);
					objInfo.vIndices.push_back(v3); objInfo.tIndices.push_back(t3); objInfo.nIndices.push_back(n3);
					objInfo.vIndices.push_back(v4); objInfo.tIndices.push_back(t4); objInfo.nIndices.push_back(n4);
				}
				break;
			default:
				break;
			}
			break;
		default:
			fgets(buf, sizeof(buf), file);
			break;
		}
	}
}

void readObj(char *filename, std::vector<Vert>& vertices, std::vector<int>& indices)
{
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(string(filename), aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_CalcTangentSpace);
	if (!scene || scene->mFlags == AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) // if is Not Zero
	{
		cout << "ERROR::ASSIMP:: " << importer.GetErrorString() << endl;
		return;
	}

	// Process ASSIMP's root node recursively
	processNode(scene->mRootNode, scene, vertices, indices);
}

void processNode(aiNode* node, const aiScene* scene, std::vector<Vert>& vertices, std::vector<int>& indices)
{
	// Process each mesh located at the current node
	for (GLuint i = 0; i < node->mNumMeshes; i++)
	{
		// The node object only contains indices to index the actual objects in the scene. 
		// The scene contains all the data, node is just to keep stuff organized (like relations between nodes).
		aiMesh* mesh = scene->mMeshes[node->mMeshes[i]];
		processMesh(mesh, scene, vertices, indices);
	}
	// After we've processed all of the meshes (if any) we then recursively process each of the children nodes
	for (GLuint i = 0; i < node->mNumChildren; i++)
	{
		processNode(node->mChildren[i], scene, vertices, indices);
	}
}

void processMesh(aiMesh* mesh, const aiScene* scene, std::vector<Vert>& vertices, std::vector<int>& indices)
{
	for (GLuint i = 0; i < mesh->mNumVertices; i++)
	{
		Vert vertex;
		glm::vec3 vector;

		//position
		vector.x = mesh->mVertices[i].x;
		vector.y = mesh->mVertices[i].y;
		vector.z = mesh->mVertices[i].z;
		vertex.pos = vector;

		// Normals
		if (mesh->mNormals)
		{
			vector.x = mesh->mNormals[i].x;
			vector.y = mesh->mNormals[i].y;
			vector.z = mesh->mNormals[i].z;
			vertex.normal = vector;
		}
		else
		{
			vertex.normal = glm::vec3(0.0f, 0.0f, 0.0f);
		}

		// Texture Coordinates
		if (mesh->mTextureCoords[0]) // Does the mesh contain texture coordinates?
		{
			glm::vec2 vec;
			// A vertex can contain up to 8 different texture coordinates. We thus make the assumption that we won't 
			// use models where a vertex can have multiple texture coordinates so we always take the first set (0).
			vec.x = mesh->mTextureCoords[0][i].x;
			vec.y = mesh->mTextureCoords[0][i].y;
			vertex.uv = vec;
		}
		else
		{
			vertex.uv = glm::vec2(0.0f, 0.0f);
		}
		vertices.push_back(vertex);
	}

	for (GLuint i = 0; i < mesh->mNumFaces; i++)
	{
		aiFace face = mesh->mFaces[i];
		for (GLuint j = 0; j < face.mNumIndices; j++)
			indices.push_back(face.mIndices[j]);
	}
}

bool _AllisNum(std::string str)
{
	for (int i = 0; i < str.size(); i++)
	{
		int tmp = (int)str[i];
		if (tmp >= 48 && tmp <= 57)
		{
			continue;
		}
		else
		{
			return false;
		}
	}

	return true;
}


bool _ReadASC(char *filename, std::vector<glm::dvec3> &_vertices, int &pointCount)
{
	/*std::ifstream dataFile;

	dataFile.open(filename);

	if (!dataFile)
	{
		cout << "文件路径有误，读取失败，按任意键退出。" << endl;
		getchar();
	}

	dataFile >> pointCount;

	//读入坐标
	pointPos.clear();

	while (dataFile.peek() != EOF)
	{
		int pointNum = 0;
		glm::dvec3 tmpPos;
		float t1, t2, t3, t4, t5, t6;

		dataFile >> pointNum >> tmpPos.x >> tmpPos.y >> tmpPos.z >> t1 >> t2 >> t3 >> t4 >> t5 >> t6;
		pointPos.push_back(tmpPos);
	}


	dataFile.close();*/

	std::ifstream dataFile;

	dataFile.open(filename);

	if (!dataFile)
	{
		std::cout << "文件路径有误，读取失败，按任意键退出。" << std::endl;
		getchar();
		return false;
	}

	std::string str;
	_vertices.clear();
	unsigned index = 0;

	while (dataFile.peek() != EOF)
	{
		dataFile >> str;
		if (str == "particles")
		{
			dataFile >> pointCount;
		}
		else if (_AllisNum(str) && index < pointCount)
		{
			double posX, posY, posZ;
			dataFile >> posX >> posY >> posZ;
			_vertices.push_back(glm::vec3(posX, posY, posZ));
		}
		std::getline(dataFile, str);
	}

	dataFile.close();

	return true;
}

void MatToArray(glm::mat4 mat, double* dst)
{
	for (unsigned i = 0; i < 4; i++)
	{
		for (unsigned j = 0; j < 4; j++)
		{
			dst[i * 4 + j] = mat[i][j];
		}
	}
}

void ArrayToMat(double* src, glm::mat4& mat)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = src[i * 4 + j];
		}
	}
}