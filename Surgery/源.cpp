#if defined(NANOGUI_GLAD)
#if defined(NANOGUI_SHARED) && !defined(GLAD_GLAPI_EXPORT)
#define GLAD_GLAPI_EXPORT
#endif



#include "Program.h"
#include "Window.h"

#include "Model.h"
#include "Haptic.h"

#include <glm/gtc/matrix_transform.hpp>
#include <time.h>

//#include <glad/glad.h>
#else
#if defined(__APPLE__)
#define GLFW_INCLUDE_GLCOREARB
#else
#define GL_GLEXT_PROTOTYPES
#endif
#endif

//#include <GLFW/glfw3.h>

#include <nanogui/nanogui.h>

#include <iostream>


using namespace std;


//nanogui
nanogui::Screen *screen = nullptr;
nanogui::Screen *screen_light = nullptr;
nanogui::ref<nanogui::Window> nanoguiWindow;
nanogui::ref<nanogui::Window> nanoguiWindow_light;


//window
int16_t IMAGE_WIDTH = 1320;
int16_t IMAGE_HEIGHT = 800;
Window *g_pWindow;


// Mouse Movement 
int oldMouseX, oldMouseY;
int currentMouseX, currentMouseY;
bool isLeftDown = false, isRightDown = false, isMiddleDown = false;
float xTranslate = 0.0f, yTranslate = 0.0f, zTranslate = 0.0f;
float xRotation = 0.0f, yRotation = 0.0f, zRotation = 0.0f;
float zoom = 1.0f;

//Modern GL
struct GLInfo
{
	/*static GLuint headVAO;
	static GLuint headVBO;
	static GLuint headEBO;

	static GLuint cursorVAO;
	static GLuint cursorVBO;*/

	static GLuint pointCloudVAO;
	static GLuint pointCloudVBO;
};

/*GLuint GLInfo::headVBO = -1;
GLuint GLInfo::headVAO = -1;
GLuint GLInfo::headEBO = -1;

GLuint GLInfo::cursorVAO = -1;
GLuint GLInfo::cursorVBO = -1;*/

GLuint GLInfo::pointCloudVAO = -1;
GLuint GLInfo::pointCloudVBO = -1;

//shader
Program *g_pHeadProgram;
Program *g_pPointCloudProgram;

//model material
typedef struct _ModelMat
{
	nanogui::Color modelSurfaceColor = nanogui::Color(1.0f, 0.84f, 0.0f, 1.0f);
	glm::vec3 materialSpecular = glm::vec3(0.1f);
	float materialShininess = 32.0f;
} ModelMat;
ModelMat headMat;

//point cloud
std::vector<glm::dvec3> pointCloud;
int pointCount;

//head
Model model;
Model cursor;
bool isDrawModel = true;

//light
Light light;


//Function prototypes
void drawPointCloud();
void setUpPointCloud();
void deleteMesh();

// callback
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double x, double y);
void windowsize_callback(GLFWwindow* window, int width, int height); 
void mouseButton_callback(GLFWwindow* window, int button, int action, int mods);

int main(int argc, char** argv)
{
	std::srand(time(NULL));

	//create window
	g_pWindow = new Window("Fur Simulation", IMAGE_WIDTH, IMAGE_HEIGHT);
	g_pWindow->SetKeyCallback(key_callback);
	g_pWindow->SetCursorPosCallback(mouse_callback);
	g_pWindow->SetWindowSizeCallback(windowsize_callback);
	g_pWindow->SetMouseButtonCallback(mouseButton_callback);
	g_pWindow->SetCharCallback([](GLFWwindow *, unsigned int codepoint) {
		screen->charCallbackEvent(codepoint);
		screen_light->charCallbackEvent(codepoint);
	});
	g_pWindow->SetDropCallback([](GLFWwindow *, int count, const char **filenames) {
		screen->dropCallbackEvent(count, filenames);
		screen_light->dropCallbackEvent(count, filenames);
	});
	g_pWindow->SetScrollCallback([](GLFWwindow *, double x, double y) {
		screen->scrollCallbackEvent(x, y);
		screen_light->scrollCallbackEvent(x, y);
	});
	g_pWindow->SetFrameBufferSizeCallback([](GLFWwindow *, int width, int height) {
		screen->resizeCallbackEvent(width, height);
		screen_light->resizeCallbackEvent(width, height);
	});

	//nanogui
	{
		// Create a nanogui screen and pass the glfw pointer to initialize
		bool enabled = true;
		screen = new nanogui::Screen();
		screen->initialize(g_pWindow->GetGLFWWindow(), true);
		nanogui::FormHelper *gui = new nanogui::FormHelper(screen);
		nanoguiWindow = gui->addWindow(Eigen::Vector2i(10, 10), "Dashboard");

		gui->addGroup("Model Config");
		gui->addVariable("Draw Model", isDrawModel);
		gui->addVariable("Model Color", headMat.modelSurfaceColor);
		gui->addVariable("Material Specular X", headMat.materialSpecular.x)->setSpinnable(true);
		gui->addVariable("Material Specular Y", headMat.materialSpecular.y)->setSpinnable(true);
		gui->addVariable("Material Specular Z", headMat.materialSpecular.z)->setSpinnable(true);
		gui->addVariable("Material Shininess", headMat.materialShininess)->setSpinnable(true);

		screen->setVisible(true);
		screen->performLayout();
	}
	{
		screen_light = new nanogui::Screen();
		screen_light->initialize(g_pWindow->GetGLFWWindow(), true);
		nanogui::FormHelper *gui = new nanogui::FormHelper(screen_light);
		nanoguiWindow_light = gui->addWindow(Eigen::Vector2i(IMAGE_WIDTH - 210, 10), "Light");
		
		gui->addGroup("Light Config");
		gui->addVariable("Light Position X", light.lightPos.x)->setSpinnable(true);
		gui->addVariable("Light Position Y", light.lightPos.y)->setSpinnable(true);
		gui->addVariable("Light Position Z", light.lightPos.z)->setSpinnable(true);
		gui->addVariable("Light Ambient R", light.lightAmbient.r)->setSpinnable(true);
		gui->addVariable("Light Ambient G", light.lightAmbient.g)->setSpinnable(true);
		gui->addVariable("Light Ambient B", light.lightAmbient.b)->setSpinnable(true);
		gui->addVariable("Light Diffuse R", light.lightDiffuse.r)->setSpinnable(true);
		gui->addVariable("Light Diffuse G", light.lightDiffuse.g)->setSpinnable(true);
		gui->addVariable("Light Diffuse B", light.lightDiffuse.b)->setSpinnable(true);
		gui->addVariable("Light Specular R", light.lightSpecular.r)->setSpinnable(true);
		gui->addVariable("Light Specular G", light.lightSpecular.g)->setSpinnable(true);
		gui->addVariable("Light Specular B", light.lightSpecular.b)->setSpinnable(true);
		gui->addVariable("Constant Damping", light.lightConstant)->setSpinnable(true);
		gui->addVariable("Linear Damping", light.lightLinear)->setSpinnable(true);
		gui->addVariable("Quadratic Damping", light.lightQuadratic)->setSpinnable(true);

		screen_light->setVisible(true);
		screen_light->performLayout();
	}

	//point cloud
	readASC("../resources/PointClouds/bunny.asc", pointCloud, pointCount);

	//models
	model.LoadFromFile("../resources/bunny.obj");
	model.SetUpGLBuffer(GL_STATIC_DRAW);
	cursor.LoadFromFile("../resources/Cone.obj");
	cursor.SetUpGLBuffer(GL_STATIC_DRAW);


	//load shader & init
	std::vector<ShaderYH> shaders;
	shaders.push_back(ShaderYH::ShaderYHFromFile("../resources/HeadShader.vs", GL_VERTEX_SHADER));
	shaders.push_back(ShaderYH::ShaderYHFromFile("../resources/HeadShader.frag", GL_FRAGMENT_SHADER));
	g_pHeadProgram = new Program(shaders);

	shaders.clear();
	shaders.push_back(ShaderYH::ShaderYHFromFile("../resources/pointCloudShader.vs", GL_VERTEX_SHADER));
	shaders.push_back(ShaderYH::ShaderYHFromFile("../resources/pointCloudShader.frag", GL_FRAGMENT_SHADER));
	g_pPointCloudProgram = new Program(shaders);

	//set up point cloud
	setUpPointCloud();


	//haptic
	Haptic::InitHL();

	// Game loop
	while (!g_pWindow->ShouldClose())
	{
		// Check if any events have been activated (key pressed, mouse moved etc.) and call corresponding response functions
		glfwPollEvents();
	
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);

		int pixelWidth, pixelHeight;
		glfwGetFramebufferSize(g_pWindow->GetGLFWWindow(), &pixelWidth, &pixelHeight);
		glViewport(0, 0, pixelWidth, pixelHeight);

		glm::mat4 projMat = glm::perspective(45.0f, (float)IMAGE_WIDTH / IMAGE_HEIGHT, 0.1f, 1000.0f);
		glm::mat4 modelView = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -10.0f));
		modelView = glm::translate(modelView, glm::vec3(xTranslate, yTranslate, 0.0f));
		
		modelView = glm::rotate(modelView, yRotation, glm::vec3(1, 0, 0));
		modelView = glm::rotate(modelView, xRotation, glm::vec3(0, 1, 0));
		modelView = glm::scale(modelView, glm::vec3(zoom, zoom, zoom));

		//model
		if (isDrawModel)
		{
			g_pHeadProgram->Use();
			{
				g_pHeadProgram->SetUniform("projection", projMat);
				g_pHeadProgram->SetUniform("modelView", modelView);

				g_pHeadProgram->SetUniform("lightPosition", light.lightPos);
				g_pHeadProgram->SetUniform("lightAmbient", light.lightAmbient);
				g_pHeadProgram->SetUniform("lightDiffuse", light.lightDiffuse);
				g_pHeadProgram->SetUniform("lightSpecular", light.lightSpecular);

				g_pHeadProgram->SetUniform("lightConstant", light.lightConstant);
				g_pHeadProgram->SetUniform("lightLinear", light.lightLinear);
				g_pHeadProgram->SetUniform("lightQuadratic", light.lightQuadratic);

				g_pHeadProgram->SetUniform("surfaceColor", glm::vec3(headMat.modelSurfaceColor.r(), headMat.modelSurfaceColor.g(), headMat.modelSurfaceColor.b()));
				g_pHeadProgram->SetUniform("materialSpecular", headMat.materialSpecular);
				g_pHeadProgram->SetUniform("materialShininess", headMat.materialShininess);

				glBindVertexArray(model.VAO);
				glDrawElements(GL_TRIANGLES, model.indices.size(), GL_UNSIGNED_INT, model.indices.data());
				glBindVertexArray(0);

				//haptic
				int viewport[4] = {0, 0, pixelWidth, pixelHeight};
				double modelViewArray[16];
				MatToArray(modelView, modelViewArray);
				Haptic::UpdateWorkspace(modelView, projMat, viewport);
				hlMatrixMode(HL_MODELVIEW);
				hlLoadMatrixd(modelViewArray);

				hlBeginFrame();
				{
					hlMaterialf(HL_FRONT_AND_BACK, HL_STIFFNESS, 0.7f);
					hlMaterialf(HL_FRONT_AND_BACK, HL_DAMPING, 0.1f);
					hlMaterialf(HL_FRONT_AND_BACK, HL_STATIC_FRICTION, 0.2f);
					hlMaterialf(HL_FRONT_AND_BACK, HL_DYNAMIC_FRICTION, 0.3f);

					hlBeginShape(HL_SHAPE_FEEDBACK_BUFFER, Haptic::gShapeId);
					{
						glBindVertexArray(model.VAO);
						glDrawElements(GL_TRIANGLES, model.indices.size(), GL_UNSIGNED_INT, model.indices.data());
						glBindVertexArray(0);
					}
					hlEndShape();
				}
				hlEndFrame();

				//cursor
				double proxyTransform[16];
				hlGetDoublev(HL_PROXY_TRANSFORM, proxyTransform);
				glm::mat4 proxyTransMat;
				ArrayToMat(proxyTransform, proxyTransMat);
				proxyTransMat = glm::scale(proxyTransMat, glm::vec3(Haptic::gCursorScale));
				g_pHeadProgram->SetUniform("modelView", modelView * proxyTransMat);
				glBindVertexArray(cursor.VAO);
				glDrawElements(GL_TRIANGLES, cursor.indices.size(), GL_UNSIGNED_INT, cursor.indices.data());
				glBindVertexArray(0);
			}
			g_pHeadProgram->End();
		}
	
		//point cloud
		/*g_pPointCloudProgram->Use();
		{
			g_pPointCloudProgram->SetUniform("projection", projMat);
			g_pPointCloudProgram->SetUniform("modelView", modelView);

			g_pPointCloudProgram->SetUniform("surfaceColor", glm::vec3(headMat.modelSurfaceColor.r(), headMat.modelSurfaceColor.g(), headMat.modelSurfaceColor.b()));

			glBindVertexArray(GLInfo::pointCloudVAO);
			glDrawArrays(GL_POINTS, 0, pointCount);
			glBindVertexArray(0);
		}
		g_pPointCloudProgram->End();*/

		//haptic render
		/*int viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		Haptic::UpdateWorkspace(modelView, projMat, viewport);
		Haptic::DrawSceneHaptics(&model, GLInfo::headVAO, GL_TRIANGLES, true);*/

		//cursor
		/*g_pHeadProgram->Use();
		{
			glm::mat4 cursorModelView = Haptic::GetProxyTransform();//¾ØÕóÓÐÎÊÌâ£¡
			cursorModelView = glm::scale(cursorModelView, glm::vec3(Haptic::gCursorScale, Haptic::gCursorScale, Haptic::gCursorScale));

			g_pHeadProgram->SetUniform("projection", projMat);
			g_pHeadProgram->SetUniform("modelView", modelView * cursorModelView);

			g_pHeadProgram->SetUniform("lightPosition", light.lightPos);
			g_pHeadProgram->SetUniform("lightAmbient", light.lightAmbient);
			g_pHeadProgram->SetUniform("lightDiffuse", light.lightDiffuse);
			g_pHeadProgram->SetUniform("lightSpecular", light.lightSpecular);

			g_pHeadProgram->SetUniform("lightConstant", light.lightConstant);
			g_pHeadProgram->SetUniform("lightLinear", light.lightLinear);
			g_pHeadProgram->SetUniform("lightQuadratic", light.lightQuadratic);

			g_pHeadProgram->SetUniform("surfaceColor", glm::vec3(headMat.modelSurfaceColor.r(), headMat.modelSurfaceColor.g(), headMat.modelSurfaceColor.b()));
			g_pHeadProgram->SetUniform("materialSpecular", headMat.materialSpecular);
			g_pHeadProgram->SetUniform("materialShininess", headMat.materialShininess);

			DrawGL_Model(cursor, GLInfo::cursorVAO, GL_TRIANGLES, true);
		}
		g_pHeadProgram->End();*/


		// Draw nanogui
		screen->drawContents();
		screen->drawWidgets();
		screen_light->drawContents();
		screen_light->drawWidgets();
	
		g_pWindow->SwapBuffers();
	}
	
	// Terminate GLFW, clearing any resources allocated by GLFW.
	glfwTerminate();

	delete g_pWindow;
	delete g_pHeadProgram;
	delete g_pPointCloudProgram;

	//delete haptic
	Haptic::ExitHandler();

	return 0;
}


void drawPointCloud()
{
	glBindVertexArray(GLInfo::pointCloudVAO);
	glDrawArrays(GL_POINTS, 0, pointCount);
	glBindVertexArray(0);
}

void setUpPointCloud()
{
	glGenVertexArrays(1, &GLInfo::pointCloudVAO);
	glGenBuffers(1, &GLInfo::pointCloudVBO);

	glBindBuffer(GL_ARRAY_BUFFER, GLInfo::pointCloudVBO);
	glBufferData(GL_ARRAY_BUFFER, pointCloud.size() * sizeof(glm::dvec3), pointCloud.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(GLInfo::pointCloudVAO);
	{
		glBindBuffer(GL_ARRAY_BUFFER, GLInfo::pointCloudVBO);
		{
			//pos
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(glm::dvec3), (GLvoid*)0);
		}
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
	glBindVertexArray(0);
}

void deleteMesh()
{
	glDeleteVertexArrays(1, &GLInfo::pointCloudVAO);
	glDeleteBuffers(1, &GLInfo::pointCloudVBO);
}



void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
	screen->keyCallbackEvent(key, scancode, action, mode);
	screen_light->keyCallbackEvent(key, scancode, action, mode);

	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case GLFW_KEY_UP:
			
			break;
		case GLFW_KEY_DOWN:
			
			break;
		case GLFW_KEY_LEFT:
			
			break;
		case GLFW_KEY_RIGHT:
			
			break;
		}
	}
}

void mouse_callback(GLFWwindow* window, double x, double y)
{
	screen->cursorPosCallbackEvent(x, y);
	screen_light->cursorPosCallbackEvent(x, y);

	currentMouseX = x;
	currentMouseY = y;

	if (isLeftDown)
	{
		xRotation += (currentMouseX - oldMouseX) / 125.0f;
		yRotation += (currentMouseY - oldMouseY) / 125.0f;
	}

	if (isMiddleDown)
	{
		xTranslate += (currentMouseX - oldMouseX) / 125.0f;
		yTranslate -= (currentMouseY - oldMouseY) / 125.0f;
	}

	if (isRightDown)
	{
		zoom -= (currentMouseY - oldMouseY) / 125.0f;
	}

	oldMouseX = currentMouseX;
	oldMouseY = currentMouseY;
}

void windowsize_callback(GLFWwindow* window, int width, int height)
{
	IMAGE_WIDTH = width;
	IMAGE_HEIGHT = height;

	nanoguiWindow->setPosition(Eigen::Vector2i(10, 10));
	nanoguiWindow_light->setPosition(Eigen::Vector2i(IMAGE_WIDTH - 210, 10));
}

void mouseButton_callback(GLFWwindow* window, int button, int action, int mods)
{
	screen->mouseButtonCallbackEvent(button, action, mods);
	screen_light->mouseButtonCallbackEvent(button, action, mods);

	if (action == GLFW_RELEASE)
	{
		switch (button)
		{
		case GLFW_MOUSE_BUTTON_LEFT:
			isLeftDown = false;
			break;
		case GLFW_MOUSE_BUTTON_RIGHT:
			isRightDown = false;
			break;
		case GLFW_MOUSE_BUTTON_MIDDLE:
			isMiddleDown = false;
			break;
		}
	}

	if (action == GLFW_PRESS)
	{
		switch (button)
		{
		case GLFW_MOUSE_BUTTON_LEFT:
			isLeftDown = true;
			break;
		case GLFW_MOUSE_BUTTON_RIGHT:
			isRightDown = true;
			break;
		case GLFW_MOUSE_BUTTON_MIDDLE:
			isMiddleDown = true;
			break;
		}
	}
}