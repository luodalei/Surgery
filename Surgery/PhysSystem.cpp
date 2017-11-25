#include "PhysSystem.h"
#include <fstream>


PhysSystem::PhysSystem(std::string filename)
{
	vertices = cVertexArray::create(false, false, true, false, false, false);

	pointCloud = new cMultiPoint();
	stretches = new cMultiSegment();

	this->addChild(pointCloud);
	//this->addChild(stretches);

	LoadFromASC(filename);

	isShowPointCloud = true;
}


PhysSystem::~PhysSystem()
{
	
}

bool AllStringisNum(std::string str)
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

bool ReadFromASC(const char *filename, cVertexArrayPtr _vertices, int &pointCount)
{
	std::ifstream dataFile;

	dataFile.open(filename);

	if (!dataFile)
	{
		std::cout << "文件路径有误，读取失败，按任意键退出。" << std::endl;
		getchar();
		return false;
	}

	std::string str;
	_vertices->clear();
	unsigned index = 0;

	while (dataFile.peek() != EOF)
	{
		dataFile >> str;
		if (str == "particles")
		{
			dataFile >> pointCount;
		}
		else if (AllStringisNum(str) && index < pointCount)
		{
			double posX, posY, posZ;
			dataFile >> posX >> posY >> posZ;
			//pointPos[index++] = cVector3d(posX, posY, posZ);

			unsigned index = _vertices->newVertex();

			_vertices->setLocalPos(index, cVector3d(posX, posY, posZ));
			_vertices->setNormal(index, cVector3d(1, 0, 0));
			_vertices->setTexCoord(index, cVector3d(0, 0, 0));
			_vertices->setColor(index, cColorf(1.0, 1.0, 1.0));
		}
		std::getline(dataFile, str);
	}

	dataFile.close();

	return true;
}

bool PhysSystem::LoadFromASC(std::string filename)
{
	int pointCount = 0;
	if (!ReadFromASC(filename.c_str(), vertices, pointCount))
		return false;

	pointCloud->m_vertices = vertices;
	pointCloud->m_points->m_vertices = vertices;

	stretches->m_vertices = vertices;
	stretches->m_segments->m_vertices = vertices;

	for (unsigned i = 0; i < vertices->getNumElements(); i++)
	{
		pointCloud->newPoint(i);


	}

	return true;
}
