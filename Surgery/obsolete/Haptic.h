#pragma once

#include <HL/hl.h>
#include <HDU/hduMatrix.h>
#include <HDU/hduError.h>

#include <HLU/hlu.h>

#include "Model.h"

#define CURSOR_SIZE_PIXELS 20

class Haptic
{
public:
	Haptic();
	~Haptic();

	static void InitHL();
	static void ExitHandler(void);
	static void DrawSceneHaptics(Model* model, const GLuint _VAO, GLenum mode, bool useIndex = true);
	static void UpdateWorkspace(glm::mat4 _modelView, glm::mat4 _projection, int _viewport[4]);
	static glm::mat4 GetProxyTransform();

public:
	// Haptic device and rendering context handles. 
	static HHD ghHD;
	static HHLRC ghHLRC;

	// Shape id for shape we will render haptically. 
	static HLuint gShapeId;

	static double gCursorScale;
};

