#include "chai3d.h"

//#include "GLFW\glfw3.h"
#include "WindowYH.h"
#include "PhysModel.h"



using namespace chai3d;
using namespace std;

// stereo Mode
/*
C_STEREO_DISABLED:            Stereo is disabled
C_STEREO_ACTIVE:              Active stereo for OpenGL NVDIA QUADRO cards
C_STEREO_PASSIVE_LEFT_RIGHT:  Passive stereo where L/R images are rendered next to each other
C_STEREO_PASSIVE_TOP_BOTTOM:  Passive stereo where L/R images are rendered above each other
*/
cStereoMode stereoMode = C_STEREO_DISABLED;


// mirrored display
bool mirroredDisplay = false;


//------------------------------------------------------------------------------
// DECLARED VARIABLES
//------------------------------------------------------------------------------

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cSpotLight *light;

#define FINALLY_USE_HAPTIC
#ifdef FINALLY_USE_HAPTIC
// a haptic device handler
cHapticDeviceHandler* handler;

// a pointer to the current haptic device
cGenericHapticDevicePtr hapticDevice;

// a virtual tool representing the haptic device in the scene
cToolCursor* tool;

// define the radius of the tool (sphere)
double toolRadius = 0.02;

// a frequency counter to measure the simulation haptic rate
cFrequencyCounter freqCounterHaptics;

// haptic thread
cThread* hapticsThread;
#endif
// a few shape primitives that compose our scene
cMultiMesh* g_pBunny;

// a font for rendering text
cFontPtr font;

// a label to display the rate [Hz] at which the simulation is running
cLabel* labelRates;

// a small scope to display the interaction force signal
cScope* scope;

// a flag to indicate if the haptic simulation currently running
bool simulationRunning = false;

// a flag to indicate if the haptic simulation has terminated
bool simulationFinished = true;

// a frequency counter to measure the simulation graphic rate
cFrequencyCounter freqCounterGraphics;

// a handle to window display context
WindowYH* g_pWindow = NULL;

//particle物理模型
PhysModel* g_pPhysModel;

// state machine 
enum STATE 
{
	STATE_IDLE,
	STATE_MODIFY_OBJECT,
	STATE_MOVE_CAMERA
};
STATE state = STATE_IDLE;

//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

// callback when the window display is resized
void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height);

// callback when an error GLFW occurs
void errorCallback(int error, const char* a_description);

// callback when a key is pressed
void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);

// this function renders the scene
void updateGraphics(void);

// this function contains the main haptics simulation loop
void updateHaptics(void);

// this function closes the application
void close(void);

//yyh添加

//phys deformation
void PhysDeformation(PhysModel* physModel);


//==============================================================================
/*
DEMO:   04-shapes.cpp

In this example we illustrate how to create some simple scene composed
of a few shape primitives that include the sphere, cylinder and line.

In the haptic threads we compute the interaction between the tool (cursor)
and the different object composing the scene.

A widget is also placed on the front plane to display information relative
to the position of the cylinder.
*/
//==============================================================================

//nanogui
#ifdef FINALLY_USE_NANOGUI
#include <nanogui\nanogui.h>
nanogui::Screen *screen = nullptr;
nanogui::ref<nanogui::Window> nanoguiWindow;

typedef struct _NanoGUIInfo
{
	//只有同时满足静态常量整型才能在类内初始化
	//Graphics Model
	static bool isDrawSurface;
	static nanogui::Color materialAmbient;
	static nanogui::Color materialDiffuse;
	static nanogui::Color materialSpecular;
	static unsigned materialShininess;

	static float transparencyLevel;

	//Phys Moddel
	static bool isDrawParticles;
	static nanogui::Color particleColor;
} NanoGUIInfo;
bool NanoGUIInfo::isDrawSurface = true;
nanogui::Color NanoGUIInfo::materialAmbient(0.3f, 0.3f, 0.3f, 1.0f);
nanogui::Color NanoGUIInfo::materialDiffuse(0.7f, 0.7f, 0.7f, 1.0f);
nanogui::Color NanoGUIInfo::materialSpecular(0.1f, 0.1f, 0.1f, 1.0f);
unsigned NanoGUIInfo::materialShininess = 32;
float NanoGUIInfo::transparencyLevel = 0.7f;

bool NanoGUIInfo::isDrawParticles = true;
nanogui::Color NanoGUIInfo::particleColor(1.0f, 1.0f, 1.0f, 1.0f);
#endif
// convert to resource path
#define RESOURCE_PATH(p)    (char*)((string(p)).c_str())

int main(int argc, char* argv[])
{
	//--------------------------------------------------------------------------
	// INITIALIZATION
	//--------------------------------------------------------------------------

	cout << endl;
	cout << "-----------------------------------" << endl;
	cout << "CHAI3D" << endl;
	cout << "Demo: 04-shapes" << endl;
	cout << "Copyright 2003-2016" << endl;
	cout << "-----------------------------------" << endl << endl << endl;
	cout << "Keyboard Options:" << endl << endl;
	cout << "[m] - Enable/Disable vertical mirroring" << endl;
	cout << "[q] - Exit application" << endl;
	cout << endl << endl;


	//--------------------------------------------------------------------------
	// OPEN GL - WINDOW DISPLAY
	//--------------------------------------------------------------------------

	g_pWindow = new WindowYH("Surgery", 1300, 768);

	glfwSetErrorCallback(errorCallback);
	g_pWindow->SetKeyCallback([](GLFWwindow *window, int key, int scancode, int action, int mods) {
		keyCallback(window, key, scancode, action, mods);
#ifdef FINALLY_USE_NANOGUI
		screen->keyCallbackEvent(key, scancode, action, mods);
#endif
	});
#ifdef FINALLY_USE_NANOGUI
	g_pWindow->SetCursorPosCallback([](GLFWwindow *window, double x, double y) {
		screen->cursorPosCallbackEvent(x, y);
	});
	g_pWindow->SetMouseButtonCallback([](GLFWwindow *window, int button, int action, int modifiers) {
		screen->mouseButtonCallbackEvent(button, action, modifiers);
	});
	g_pWindow->SetCharCallback([](GLFWwindow *window, unsigned int codepoint) {
		screen->charCallbackEvent(codepoint);
	});
	g_pWindow->SetDropCallback([](GLFWwindow *window, int count, const char** filenames) {
		screen->dropCallbackEvent(count, filenames);
	});
	g_pWindow->SetScrollCallback([](GLFWwindow *window, double x, double y) {
		screen->scrollCallbackEvent(x, y);
	});
	g_pWindow->SetFrameBufferSizeCallback([](GLFWwindow *window, int width, int height) {
		screen->resizeCallbackEvent(width, height);
	});
#endif
	g_pWindow->SetWindowSizeCallback([](GLFWwindow *window, int width, int height) {
		windowSizeCallback(window, width, height);
#ifdef FINALLY_USE_NANOGUI
		nanoguiWindow->setPosition(Eigen::Vector2i(10, 10));
#endif
	});

#ifdef FINALLY_USE_NANOGUI
	//nanogui
	{
		// Create a nanogui screen and pass the glfw pointer to initialize
		screen = new nanogui::Screen();
		screen->initialize(g_pWindow->GetGLFWWindow(), true);

		// Create nanogui gui
		bool enabled = true;
		nanogui::FormHelper *gui = new nanogui::FormHelper(screen);
		nanoguiWindow = gui->addWindow(Eigen::Vector2i(10, 10), "DashBoard");
		gui->addGroup("Geaphics");
		gui->addVariable("Draw Surface", NanoGUIInfo::isDrawSurface);
		gui->addVariable("Material Ambient", NanoGUIInfo::materialAmbient);
		gui->addVariable("Material Diffuse", NanoGUIInfo::materialDiffuse);
		gui->addVariable("Material Specular", NanoGUIInfo::materialSpecular);
		gui->addVariable("Material Shiness", NanoGUIInfo::materialShininess)->setSpinnable(true);
		gui->addVariable("Transparency Level", NanoGUIInfo::transparencyLevel)->setSpinnable(true);

		gui->addGroup("Physics");
		gui->addVariable("Draw Particles", NanoGUIInfo::isDrawParticles);
		gui->addVariable("Particle Color", NanoGUIInfo::particleColor);
		
		screen->setVisible(true);
		screen->performLayout();
	}
#endif
	//--------------------------------------------------------------------------
	// WORLD - CAMERA - LIGHTING
	//--------------------------------------------------------------------------

	// create a new world.
	world = new cWorld();

	// set the background color of the environment
	world->setBackgroundColor(cColorf(0.2f, 0.2f, 0.2f));

	// create a camera and insert it into the virtual world
	camera = new cCamera(world);
	world->addChild(camera);

	// position and orient the camera
	camera->set(cVector3d(3.0, 0.0, 0.0),    // camera position (eye)
		cVector3d(0.0, 0.0, 0.0),    // lookat position (target)
		cVector3d(0.0, 0.0, 1.0));   // direction of the (up) vector

									 // set the near and far clipping planes of the camera
	camera->setClippingPlanes(0.01, 10.0);

	// set stereo mode
	camera->setStereoMode(stereoMode);

	// set stereo eye separation and focal length (applies only if stereo is enabled)
	camera->setStereoEyeSeparation(0.02);
	camera->setStereoFocalLength(3.0);

	// set vertical mirrored display mode
	camera->setMirrorVertical(mirroredDisplay);

	//enable multi-pass rendering to handle transparent objects 模型可以透明
	//camera->setUseMultipassTransparency(true);//暂时注释掉，方便调试physmodel着色器

	//聚光灯
	light = new cSpotLight(world);

	// add light to world
	camera->addChild(light);

	// enable light source
	light->setEnabled(true);

	// position the light source
	light->setLocalPos(3.5, 2.0, 0.0);

	// define the direction of the light beam
	light->setDir(-3.5, -2.0, 0.0);

	// set light cone half angle
	light->setCutOffAngleDeg(50);

	//--------------------------------------------------------------------------
	// HAPTIC DEVICES / TOOLS
	//--------------------------------------------------------------------------
#ifdef FINALLY_USE_HAPTIC
	handler = new cHapticDeviceHandler();

	// get access to the first available haptic device found
	handler->getDevice(hapticDevice, 0);

	// retrieve information about the current haptic device
	cHapticDeviceInfo hapticDeviceInfo = hapticDevice->getSpecifications();

	// create a tool (cursor) and insert into the world
	tool = new cToolCursor(world);
	world->addChild(tool);

	// connect the haptic device to the virtual tool
	tool->setHapticDevice(hapticDevice);

	// map the physical workspace of the haptic device to a larger virtual workspace.
	tool->setWorkspaceRadius(1.0);

	// define a radius for the virtual tool (sphere)
	tool->setRadius(toolRadius);

	// create a white cursor
	tool->m_hapticPoint->m_sphereProxy->m_material->setWhiteAliceBlue();

	// oriente tool with camera
	tool->setLocalRot(camera->getLocalRot());

	// haptic forces are enabled only if small forces are first sent to the device;
	// this mode avoids the force spike that occurs when the application starts when 
	// the tool is located inside an object for instance. 
	tool->setWaitForSmallForce(true);

	// start the haptic tool
	tool->start();
#endif

	//--------------------------------------------------------------------------
	// CREATING SHAPES
	//--------------------------------------------------------------------------
#ifdef FINALLY_USE_HAPTIC
	// read the scale factor between the physical workspace of the haptic
	// device and the virtual workspace defined for the tool
	double workspaceScaleFactor = tool->getWorkspaceScaleFactor();

	// get properties of haptic device
	double maxStiffness = hapticDeviceInfo.m_maxLinearStiffness / workspaceScaleFactor;
	double maxLinearForce = cMin(hapticDeviceInfo.m_maxLinearForce, 7.0);

#endif
	////////////////////////////////////////////////////////////////////////////
	// SHAPE - Soft bunny
	////////////////////////////////////////////////////////////////////////////
	g_pBunny = new cMultiMesh();
	//world->addChild(g_pBunny);//暂时注释掉，方便调试physmodel着色器
	bool fileLoad = g_pBunny->loadFromFile("../resources/SoftBunny-1858.obj");
	if (!fileLoad)
	{
		cout << "Error - 3D Model failed to load correctly." << endl;
	}

	g_pBunny->scale(0.3);

	// rotate the object 90 degrees
	g_pBunny->rotateAboutGlobalAxisDeg(cVector3d(1, 0, 0), 90);
	g_pBunny->rotateAboutGlobalAxisDeg(cVector3d(0, 0, 1), 90);

	// enable haptic rendering on both sides of triangles
	cMaterial mat;
#ifdef FINALLY_USE_NANOGUI
	mat.m_ambient = cColorf(NanoGUIInfo::materialAmbient.r(), NanoGUIInfo::materialAmbient.g(), NanoGUIInfo::materialAmbient.b());
	mat.m_diffuse = cColorf(NanoGUIInfo::materialDiffuse.r(), NanoGUIInfo::materialDiffuse.g(), NanoGUIInfo::materialDiffuse.b());
	mat.m_specular = cColorf(NanoGUIInfo::materialSpecular.r(), NanoGUIInfo::materialSpecular.g(), NanoGUIInfo::materialSpecular.b());
	mat.setShininess(NanoGUIInfo::materialShininess);
#else
	mat.m_ambient = cColorf(0.3, 0.3, 0.3);
	mat.m_diffuse = cColorf(0.7, 0.7, 0.7);
	mat.m_specular = cColorf(0.1, 0.1, 0.1);
	mat.setShininess(32);
#endif
	mat.setHapticTriangleSides(true, false);
	g_pBunny->setMaterial(mat);

	//enable face-culling
	g_pBunny->setUseCulling(true);
#ifdef FINALLY_USE_NANOGUI
	g_pBunny->setTransparencyLevel(NanoGUIInfo::transparencyLevel);
#else
	g_pBunny->setTransparencyLevel(0.7);
#endif

#ifdef FINALLY_USE_HAPTIC
	// create collision detector
	g_pBunny->createAABBCollisionDetector(toolRadius);//必须要有！

	// set stiffness property
	g_pBunny->setStiffness(maxStiffness, true);

	// define some haptic friction properties
	g_pBunny->setFriction(0.1, 0.2, true);
#endif

#ifdef FINALLY_USE_NANOGUI
	g_pBunny->setEnabled(NanoGUIInfo::isDrawSurface);
#else
	g_pBunny->setEnabled(true);
#endif

	////////////////////////////////////////////////////////////////////////////
	// PhysModel
	////////////////////////////////////////////////////////////////////////////
	g_pPhysModel = new PhysModel();
	//读入asc文件	
	g_pPhysModel->LoadFromASC("../resources/PointClouds/SoftBunny-1858-43250p.asc");
	world->addChild(g_pPhysModel);

	g_pPhysModel->scale(0.3);

	// rotate the object 90 degrees
	g_pPhysModel->rotateAboutGlobalAxisDeg(cVector3d(1, 0, 0), 90);
	g_pPhysModel->rotateAboutGlobalAxisDeg(cVector3d(0, 0, 1), 90);

	g_pPhysModel->SetPointSize(2.0);
#ifdef FINALLY_USE_NANOGUI
	g_pPhysModel->SetPointColor(cColorf(NanoGUIInfo::particleColor.r(), NanoGUIInfo::particleColor.g(), NanoGUIInfo::particleColor.b(), NanoGUIInfo::particleColor.w()));

	g_pPhysModel->setEnabled(NanoGUIInfo::isDrawParticles);
#else
	g_pPhysModel->SetPointColor(cColorf(1.0, 1.0, 1.0));
	//g_pPhysModel->DrawAsSphere(true); //把点画成球
	g_pPhysModel->setEnabled(true);
#endif

	////////////////////////////////////////////////////////////////////////////
	// Shader - Program
	////////////////////////////////////////////////////////////////////////////
	// create program shader
	/*cShaderProgramPtr shaderProgram = cShaderProgram::create(C_SHADER_FONG_VERT, C_SHADER_FONG_FRAG);
	shaderProgram->setUniformi("uShadowMap", C_TU_SHADOWMAP);
	g_pBunny->setShaderProgram(shaderProgram);*/

	cShaderPtr vertexShader = cShader::create(C_VERTEX_SHADER);
	vertexShader->loadSourceFile(RESOURCE_PATH("../resources/Point.vert"));
	cShaderPtr geometryShader = cShader::create(C_GEOMETRY_SHADER);
	geometryShader->loadSourceFile(RESOURCE_PATH("../resources/PointTry.geom"));
	cShaderPtr fragmentShader = cShader::create(C_FRAGMENT_SHADER);
	fragmentShader->loadSourceFile(RESOURCE_PATH("../resources/Point.frag"));
	
	cShaderProgramPtr programShader = cShaderProgram::create();
	programShader->attachShader(vertexShader);
	//programShader->attachShader(geometryShader);
	programShader->attachShader(fragmentShader);

	programShader->setGeometryInputType(GL_POINTS);
	programShader->setGeometryOutputType(GL_POINTS);
	programShader->setGeometryVerticesOut(1);

	programShader->linkProgram();

	g_pPhysModel->setShaderProgram(programShader);

	
	//--------------------------------------------------------------------------
	// WIDGETS
	//--------------------------------------------------------------------------

	// create a font
	font = NEW_CFONTCALIBRI20();

	// create a label to display the haptic and graphic rate of the simulation
	labelRates = new cLabel(font);
	camera->m_frontLayer->addChild(labelRates);

	// create a scope to plot haptic device position data
	scope = new cScope();
	camera->m_frontLayer->addChild(scope);
	scope->setSize(400, 100);
	scope->setRange(0.0, 5.0);
	scope->setSignalEnabled(true, false, false, false);
	scope->setShowPanel(false);
	scope->m_colorSignal0.setRedCrimson();


	//--------------------------------------------------------------------------
	// START SIMULATION
	//--------------------------------------------------------------------------
#ifdef FINALLY_USE_HAPTIC
	// create a thread which starts the main haptics rendering loop
	hapticsThread = new cThread();
	hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

	// setup callback when application exits
	atexit(close);
#endif

	//--------------------------------------------------------------------------
	// MAIN GRAPHIC LOOP
	//--------------------------------------------------------------------------
	windowSizeCallback(g_pWindow->GetGLFWWindow(), g_pWindow->GetWidth(), g_pWindow->GetHeight());

	// main graphic loop
	while (!g_pWindow->ShouldClose())
	{
		// process events
		glfwPollEvents();

		//phys defomation

		// render graphics
		updateGraphics();

#ifdef FINALLY_USE_NANOGUI
		//nanogui
		{
			// Draw nanogui
			screen->drawContents();
			screen->drawWidgets();
		}
#endif
		// swap buffers
		g_pWindow->SwapBuffers();

		// signal frequency counter
		freqCounterGraphics.signal(1);
	}

	// close window
	glfwDestroyWindow(/*window*/g_pWindow->GetGLFWWindow());

	// terminate GLFW library
	glfwTerminate();

	// exit
	return (0);
}

//------------------------------------------------------------------------------

void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height)
{
	// update window size
	g_pWindow->SetWidthAndHeight(a_width, a_height);

	// update position of scope
	scope->setLocalPos((0.5 * (a_width - scope->getWidth())), 120);
}

//------------------------------------------------------------------------------

void errorCallback(int a_error, const char* a_description)
{
	cout << "Error: " << a_description << endl;
}

//------------------------------------------------------------------------------

void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
	// filter calls that only include a key press
	if ((a_action != GLFW_PRESS) && (a_action != GLFW_REPEAT))
	{
		return;
	}

	// option - exit
	else if ((a_key == GLFW_KEY_ESCAPE) || (a_key == GLFW_KEY_Q))
	{
		glfwSetWindowShouldClose(a_window, GLFW_TRUE);
	}
	// option - toggle vertical mirroring
	else if (a_key == GLFW_KEY_M)
	{
		mirroredDisplay = !mirroredDisplay;
		camera->setMirrorVertical(mirroredDisplay);
	}
}

//------------------------------------------------------------------------------

void close(void)
{
	// stop the simulation
	simulationRunning = false;

	// wait for graphics and haptics loops to terminate
	while (!simulationFinished) { cSleepMs(100); }

#ifdef FINALLY_USE_HAPTIC
	// close haptic device
	tool->stop();

	// delete resources
	delete hapticsThread;
	delete handler;
#endif

	delete world;
}

//------------------------------------------------------------------------------

//phys deformation
void PhysDeformation(PhysModel* physModel)
{
	//copy vertices pos & vel to GPU

	//GPU process
	 //update nerghboring
	 
	 //pbd

	 //correct position

	//copy vertices pos & vel back to CPU
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void updateGraphics(void)
{
	/////////////////////////////////////////////////////////////////////
	// UPDATE WIDGETS
	/////////////////////////////////////////////////////////////////////
#ifdef FINALLY_USE_HAPTIC
	// update haptic and graphic rate data
	labelRates->setText(cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " + cStr(freqCounterHaptics.getFrequency(), 0) + " Hz");
#endif
	// update position of label
	labelRates->setLocalPos((int)(0.5 * (g_pWindow->GetWidth() - labelRates->getWidth())), 15);

#ifdef FINALLY_USE_HAPTIC
	// update value of scope
	scope->setSignalValues(tool->getDeviceGlobalForce().length());
#endif

	/////////////////////////////////////////////////////////////////////
	// UPDATE MODEL
	/////////////////////////////////////////////////////////////////////

	// update object normals
	if (state == STATE_MODIFY_OBJECT)
	{
		g_pBunny->computeAllNormals();
	}

#ifdef FINALLY_USE_NANOGUI
	//更新模型材质
	cMaterial mat;
	mat.m_ambient = cColorf(NanoGUIInfo::materialAmbient.r(), NanoGUIInfo::materialAmbient.g(), NanoGUIInfo::materialAmbient.b());
	mat.m_diffuse = cColorf(NanoGUIInfo::materialDiffuse.r(), NanoGUIInfo::materialDiffuse.g(), NanoGUIInfo::materialDiffuse.b());
	mat.m_specular = cColorf(NanoGUIInfo::materialSpecular.r(), NanoGUIInfo::materialSpecular.g(), NanoGUIInfo::materialSpecular.b());
	mat.setShininess(NanoGUIInfo::materialShininess);
	g_pBunny->setMaterial(mat);

	g_pBunny->setTransparencyLevel(NanoGUIInfo::transparencyLevel);
	g_pBunny->setEnabled(NanoGUIInfo::isDrawSurface);

	/////////////////////////////////////////////////////////////////////
	// UPDATE Phys MODEL
	/////////////////////////////////////////////////////////////////////
	g_pPhysModel->SetPointColor(cColorf(NanoGUIInfo::particleColor.r(), NanoGUIInfo::particleColor.g(), NanoGUIInfo::particleColor.b(), NanoGUIInfo::particleColor.w()));
	g_pPhysModel->setEnabled(NanoGUIInfo::isDrawParticles);
#endif


	/////////////////////////////////////////////////////////////////////
	// RENDER SCENE
	/////////////////////////////////////////////////////////////////////

	// update shadow maps (if any)
	world->updateShadowMaps(false, mirroredDisplay);

	// render world
	camera->renderView(g_pWindow->GetWidth(), g_pWindow->GetHeight());

	// wait until all GL commands are completed
	glFinish();

	// check for any OpenGL errors
	GLenum err;
	err = glGetError();
	if (err != GL_NO_ERROR) cout << "Error:  %s\n" << gluErrorString(err);
}

//------------------------------------------------------------------------------
#ifdef FINALLY_USE_HAPTIC
void updateHaptics(void)
{
	// precision clock
	cPrecisionClock clock;
	clock.reset();

	// simulation in now running
	simulationRunning = true;
	simulationFinished = false;

	// current tool position
	cVector3d toolGlobalPos;        // global world coordinates
	cVector3d toolLocalPos;         // local coordinates

	// previous tool position
	cVector3d prevToolGlobalPos;    // global world coordinates
	cVector3d prevToolLocalPos;     // local coordinates

	state = STATE_IDLE;

	// main haptic simulation loop
	while (simulationRunning)
	{
		/////////////////////////////////////////////////////////////////////
		// SIMULATION TIME    
		/////////////////////////////////////////////////////////////////////

		// stop the simulation clock
		clock.stop();

		// read the time increment in seconds
		double timeInterval = clock.getCurrentTimeSeconds();

		// restart the simulation clock
		clock.reset();
		clock.start();


		/////////////////////////////////////////////////////////////////////
		// HAPTIC FORCE COMPUTATION
		/////////////////////////////////////////////////////////////////////

		// compute global reference frames for each object
		world->computeGlobalPositions(true);

		// update position and orientation of tool
		tool->updateFromDevice();

		// compute interaction forces
		tool->computeInteractionForces();

		// read user switch
		bool userSwitch0 = tool->getUserSwitch(0);//控制摄像机移动
		bool userSwitch1 = tool->getUserSwitch(1);//让物体变形

		// update tool position
		toolGlobalPos = tool->getDeviceGlobalPos();
		toolLocalPos = tool->getDeviceLocalPos();

		if (state == STATE_MOVE_CAMERA && !userSwitch0)
		{
			state = STATE_IDLE;

			// enable haptic interaction with map
			g_pBunny->setHapticEnabled(true, true);
		}
		else if (((state == STATE_MODIFY_OBJECT) && (!userSwitch1)) || ((state == STATE_MODIFY_OBJECT) && (!tool->isInContact(g_pBunny))))
		{
			state = STATE_IDLE;

			// enable haptic interaction with map
			g_pBunny->setHapticEnabled(true, true);

			// disable forces, 模型可能发生改变，先把力反馈关掉，重新建立AABB碰撞检测
			tool->setForcesOFF();

			// update bounding box (can take a little time)
			g_pBunny->createAABBCollisionDetector(toolRadius);

			// enable forces again
			tool->setForcesON();
		}
		else if (state == STATE_IDLE && userSwitch0)
		{
			state = STATE_MOVE_CAMERA;

			// disable haptic interaction with object
			g_pBunny->setHapticEnabled(false, true);

		}
		else if (state == STATE_IDLE && userSwitch1)
		{
			// start deforming object
			if (tool->isInContact(g_pBunny))
			{
				state = STATE_MODIFY_OBJECT;

				// disable haptic interaction with map
				g_pBunny->setHapticEnabled(false, true);
			}
		}
		else if (state == STATE_MOVE_CAMERA)
		{
			// compute tool offset
			cVector3d offset = toolLocalPos - prevToolLocalPos;

			// compute new coordinates for camera in spherical coordinates
			double radius = camera->getSphericalRadius() - 2 * offset.x();
			double azimuthDeg = camera->getSphericalAzimuthDeg() - 40 * offset.y();
			double polarDeg = camera->getSphericalPolarDeg() + 40 * offset.z();

			// update coordinates
			camera->setSphericalDeg(radius, polarDeg, azimuthDeg);

			// oriente tool with camera
			tool->setLocalRot(camera->getLocalRot());
		}
		else if (state == STATE_MODIFY_OBJECT)
		{
			//deforming object

		}
		

		// store tool position
		prevToolLocalPos = toolLocalPos;
		prevToolGlobalPos = toolGlobalPos;

		// send forces to haptic device
		tool->applyToDevice();

		// signal frequency counter
		freqCounterHaptics.signal(1);
	}

	// exit haptics thread
	simulationFinished = true;
}
#endif
//------------------------------------------------------------------------------
