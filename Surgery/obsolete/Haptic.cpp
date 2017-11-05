#include "Haptic.h"

HHD Haptic::ghHD = HD_INVALID_HANDLE;
HHLRC Haptic::ghHLRC = 0;
HLuint Haptic::gShapeId;
double Haptic::gCursorScale;

Haptic::Haptic()
{
}


Haptic::~Haptic()
{
}

void Haptic::InitHL()
{
	HDErrorInfo error;

	ghHD = hdInitDevice(HD_DEFAULT_DEVICE);
	if (HD_DEVICE_ERROR(error = hdGetError()))
	{
		hduPrintError(stderr, &error, "Failed to initialize haptic device");
		fprintf(stderr, "Press any key to exit");
		getchar();
	}

	ghHLRC = hlCreateContext(ghHD);
	hlMakeCurrent(ghHLRC);

	// Enable optimization of the viewing parameters when rendering
	// geometry for OpenHaptics.
	//hlEnable(HL_HAPTIC_CAMERA_VIEW);//这一句打开，会莫名出现一个小窗口！

	// Generate id for the shape.
	gShapeId = hlGenShapes(1);

	hlTouchableFace(HL_FRONT);
}

void Haptic::ExitHandler(void)
{
	// Deallocate the sphere shape id we reserved in initHL.
	hlDeleteShapes(gShapeId, 1);

	// Free up the haptic rendering context.
	hlMakeCurrent(NULL);
	if (ghHLRC != NULL)
	{
		hlDeleteContext(ghHLRC);
	}

	// Free up the haptic device.
	if (ghHD != HD_INVALID_HANDLE)
	{
		hdDisableDevice(ghHD);
	}
}

void Haptic::DrawSceneHaptics(Model* model, const GLuint _VAO, GLenum mode, bool useIndex)
{
	// Start haptic frame.  (Must do this before rendering any haptic shapes.)
	hlBeginFrame();

	// Set material properties for the shapes to be drawn.
	hlMaterialf(HL_FRONT_AND_BACK, HL_STIFFNESS, 0.7f);
	hlMaterialf(HL_FRONT_AND_BACK, HL_DAMPING, 0.1f);
	hlMaterialf(HL_FRONT_AND_BACK, HL_STATIC_FRICTION, 0.2f);
	hlMaterialf(HL_FRONT_AND_BACK, HL_DYNAMIC_FRICTION, 0.3f);

	// Start a new haptic shape.  Use the feedback buffer to capture OpenGL 
	// geometry for haptic rendering.
	hlBeginShape(HL_SHAPE_FEEDBACK_BUFFER, gShapeId);

	// Use OpenGL commands to create geometry.
	//DrawGL_Model(*model, _VAO, mode, useIndex);

	// End the shape.
	hlEndShape();

	// End the haptic frame.
	hlEndFrame();
}

void Haptic::UpdateWorkspace(glm::mat4 _modelView, glm::mat4 _projection, int _viewport[4])
{
	double modelView[16];
	double projection[16];

	MatToArray(_modelView, modelView);
	MatToArray(_projection, projection);

	hlMatrixMode(HL_TOUCHWORKSPACE);
	hlLoadIdentity();

	// Fit haptic workspace to view volume.
	hluFitWorkspace(projection);

	// Compute cursor scale.
	gCursorScale = hluScreenToModelScale(modelView, projection, _viewport);
	gCursorScale *= CURSOR_SIZE_PIXELS;
}

glm::mat4 Haptic::GetProxyTransform()
{
	double proxyTransform[16];
	hlGetDoublev(HL_PROXY_TRANSFORM, proxyTransform);

	glm::mat4 transformMat;
	ArrayToMat(proxyTransform, transformMat);

	double proxyPosition[3];
	hlGetDoublev(HL_PROXY_POSITION, proxyPosition);


	return transformMat;
}