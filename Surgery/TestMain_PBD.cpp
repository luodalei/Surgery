#include "chai3d.h"

//#include "GLFW\glfw3.h"
#include "WindowYH.h"
#include "PhysModel.h"
#include "PBDSystem.h"


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
typedef struct _YHInt2 {
	int x, y;

	_YHInt2(int _x, int _y) :
		x(_x),
		y(_y)
	{
	}

	_YHInt2()
	{
	}
} YHInt2;

std::vector<cShapeSphere*> g_pSpheres;
const uint g_numSpheres = 25;
const uint g_numEdges = 72;
const uint g_numSphere_X = 5;
const uint g_numSphere_Y = 5;
const uint g_numSphere_Z = 3;

cVector3d g_InitPoses[g_numSpheres] = {
	cVector3d(-0.5, 0.5, 0), cVector3d(-0.25, 0.5, 0), cVector3d(0, 0.5, 0), cVector3d(0.25, 0.5, 0), cVector3d(0.5, 0.5, 0),
	cVector3d(-0.5, 0.25, 0), cVector3d(-0.25, 0.25, 0), cVector3d(0, 0.25, 0), cVector3d(0.25, 0.25, 0), cVector3d(0.5, 0.25, 0),
	cVector3d(-0.5, 0, 0), cVector3d(-0.25, 0, 0), cVector3d(0.0, 0.0, 0), cVector3d(0.25, 0, 0), cVector3d(0.5, 0, 0),
	cVector3d(-0.5, -0.25, 0), cVector3d(-0.25, -0.25, 0), cVector3d(0, -0.25, 0), cVector3d(0.25, -0.25, 0), cVector3d(0.5, -0.25, 0),
	cVector3d(-0.5, -0.5, 0), cVector3d(-0.25, -0.5, 0), cVector3d(0, -0.5, 0), cVector3d(0.25, -0.5, 0), cVector3d(0.5, -0.5, 0),

	/*cVector3d(-0.5, 0.5, -0.25), cVector3d(-0.25, 0.5, -0.25), cVector3d(0, 0.5, -0.25), cVector3d(0.25, 0.5, -0.25), cVector3d(0.5, 0.5, -0.25),
	cVector3d(-0.5, 0.25, -0.25), cVector3d(-0.25, 0.25, -0.25), cVector3d(0, 0.25, -0.25), cVector3d(0.25, 0.25, -0.25), cVector3d(0.5, 0.25, -0.25),
	cVector3d(-0.5, 0, -0.25), cVector3d(-0.25, 0, -0.25), cVector3d(0.0, 0.0, -0.25), cVector3d(0.25, 0, -0.25), cVector3d(0.5, 0, -0.25),
	cVector3d(-0.5, -0.25, -0.25), cVector3d(-0.25, -0.25, -0.25), cVector3d(0, -0.25, -0.25), cVector3d(0.25, -0.25, -0.25), cVector3d(0.5, -0.25, -0.25),
	cVector3d(-0.5, -0.5, -0.25), cVector3d(-0.25, -0.5, -0.25), cVector3d(0, -0.5, -0.25), cVector3d(0.25, -0.5, -0.25), cVector3d(0.5, -0.5, -0.25)*/
};

YHInt2 g_Edges[g_numEdges] = {
	//行
	YHInt2(0, 1), YHInt2(1, 2), YHInt2(2, 3), YHInt2(3, 4),
	YHInt2(5, 6), YHInt2(6, 7), YHInt2(7, 8), YHInt2(8, 9),
	YHInt2(10, 11), YHInt2(11, 12), YHInt2(12, 13), YHInt2(13, 14),
	YHInt2(15, 16), YHInt2(16, 17), YHInt2(17, 18), YHInt2(18, 19),
	YHInt2(20, 21), YHInt2(21, 22), YHInt2(22, 23), YHInt2(23, 24),

	//列
	YHInt2(0, 5), YHInt2(5, 10), YHInt2(10, 15), YHInt2(15, 20),
	YHInt2(1, 6), YHInt2(6, 11), YHInt2(11, 16), YHInt2(16, 21),
	YHInt2(2, 7), YHInt2(7, 12), YHInt2(12, 17), YHInt2(17, 22),
	YHInt2(3, 8), YHInt2(8, 13), YHInt2(13, 18), YHInt2(18, 23),
	YHInt2(4, 9), YHInt2(9, 14), YHInt2(14, 19), YHInt2(19, 24),

	//右斜
	YHInt2(0, 6), YHInt2(1, 7), YHInt2(2, 8), YHInt2(3, 9),
	YHInt2(5, 11), YHInt2(6, 12), YHInt2(7, 13), YHInt2(8, 14),
	YHInt2(10, 16), YHInt2(11, 17), YHInt2(12, 18), YHInt2(13, 19),
	YHInt2(15, 21), YHInt2(16, 22), YHInt2(17, 23), YHInt2(18, 24),

	//左斜
	YHInt2(1, 5), YHInt2(2, 6), YHInt2(3, 7), YHInt2(4, 8),
	YHInt2(6, 10), YHInt2(7, 11), YHInt2(8, 12), YHInt2(9, 13),
	YHInt2(11, 15), YHInt2(12, 16), YHInt2(13, 17), YHInt2(14, 18),
	YHInt2(16, 20), YHInt2(17, 21), YHInt2(18, 22), YHInt2(19, 23),
};

double g_InitInvMasses[g_numSpheres] = {
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

std::vector<cVector3d> g_SpherePoses;
std::vector<cVector3d> g_SphereOldPoses;
std::vector<cVector3d> g_SphereVels;
std::vector<double> g_SphereInvMasses;
double g_dt = 0.01;

typedef struct _YHDistanceConstraint {

public:
	_YHDistanceConstraint() 
	{

	}

	bool InitConstraint(cVector3d* originalPos, const unsigned i1, const unsigned i2)
	{
		vIdx[0] = i1;
		vIdx[1] = i2;

		const cVector3d &x1_0 = originalPos[i1];
		const cVector3d &x2_0 = originalPos[i2];

		restLength = (x2_0 - x1_0).length();

		return true;
	}

	bool SolveConstraint(cVector3d* pos, double* invMasses)
	{
		const unsigned i1 = vIdx[0];
		const unsigned i2 = vIdx[1];

		cVector3d &x1 = pos[i1];
		cVector3d &x2 = pos[i2];

/*		std::cout << "i1, i2: " << i1 << " " << i2 << std::endl;
		std::cout << "x1: " << x1.x() << " " << x1.y() << " " << x1.z() << std::endl;
		std::cout << "x2: " << x2.x() << " " << x2.y() << " " << x2.z() << std::endl;
*/
		const double invMass1 = invMasses[i1];
		const double invMass2 = invMasses[i2];

		const double wSum = invMass1 + invMass2;
		if (wSum == 0.0)
			return false;

//		std::cout << "wSum: " << wSum << std::endl;

		cVector3d n = x1 - x2;
		double d = n.length();
		n.normalize();

//		std::cout << "d: " << d << std::endl;

		cVector3d corr, corr1, corr2;
		if (d < restLength)
			corr = compressionStiffness * n * (d - restLength) / wSum;
		else
			corr = stretchStiffness * n * (d - restLength) / wSum;

//		std::cout << "corr: " << corr << std::endl;

		//这里跟position-based-dynamics-master项目中负号相反
		corr1 = -invMass1 * corr;
		corr2 = invMass2 * corr;

		if (invMass1 != 0.0)
			x1 += corr1;
		if (invMass2 != 0.0)
			x2 += corr2;

/*		std::cout << "new-x1: " << x1.x() << " " << x1.y() << " " << x1.z() << std::endl;
		std::cout << "new-x2: " << x2.x() << " " << x2.y() << " " << x2.z() << std::endl << std::endl;
*/
		return true;
	}

public:
	//indices of linked verts
	unsigned vIdx[2];

	//distance constraint
	double restLength;
	
	static double compressionStiffness;
	static double stretchStiffness;

} YHDistanceConstraint;

double YHDistanceConstraint::compressionStiffness = 0.8;
double YHDistanceConstraint::stretchStiffness = 0.8;

std::vector<YHDistanceConstraint*> g_DistanceConstraints;

//physModel接收tool的力，按下F开启/关闭
bool applyForceToPhysModel = false;

// current tool position
cVector3d toolGlobalPos;        // global world coordinates
cVector3d toolLocalPos;         // local coordinates

								// previous tool position
cVector3d prevToolGlobalPos;    // global world coordinates
cVector3d prevToolLocalPos;     // local coordinates


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
//compute vel & pos
void SemiImplicitEuler(const double dt, const double invMass, const cVector3d force, cVector3d& position, cVector3d& velocity);

//velocity update first order
void VelocityUpdateFirstOrder(const double dt, const double invMass, const cVector3d &position, const cVector3d &oldPosition, cVector3d &velocity);

//phys deformation
void PhysDeformation(cVector3d toolPos, cVector3d toolForce, double dst, bool applyForceOnPhysModel);


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
	camera->set(cVector3d(10.0, 0.0, 0.0),    // camera position (eye)
		cVector3d(0.0, 0.0, 0.0),    // lookat position (target)
		cVector3d(0.0, 0.0, 1.0));   // direction of the (up) vector

									 // set the near and far clipping planes of the camera
	camera->setClippingPlanes(0.01, 100.0);

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
	// Sphere Model
	////////////////////////////////////////////////////////////////////////////
	g_pSpheres.resize(g_numSpheres);
	g_SpherePoses.clear();
	g_SphereOldPoses.clear();

	cMaterial matSphere;
	matSphere.setWhiteAliceBlue();
	for (unsigned i = 0; i < g_numSpheres; i++)
	{
		g_pSpheres[i] = new cShapeSphere(0.04);
		world->addChild(g_pSpheres[i]);
		g_pSpheres[i]->setMaterial(matSphere);
		g_pSpheres[i]->setShowEnabled(true);

		//初始化位置
		g_SpherePoses.push_back(g_InitPoses[i]);
		g_SphereOldPoses.push_back(g_InitPoses[i]);
	}

	//初始化速度
	g_SphereVels.resize(g_numSpheres, cVector3d(0, 0, 0));

	//初始化质量
	for (unsigned i = 0; i < g_numSpheres; i++)
	{
		g_SphereInvMasses.push_back(g_InitInvMasses[i]);
	}

	////////////////////////////////////////////////////////////////////////////
	// Constraints
	////////////////////////////////////////////////////////////////////////////
	cMultiSegment* segments = new cMultiSegment();
	world->addChild(segments);

	for (unsigned i = 0; i < g_numSpheres; i++)
	{
		segments->newVertex(g_InitPoses[i]);
	}
	for (unsigned i = 0; i < g_numEdges; i++)
	{
		segments->newSegment(g_Edges[i].x, g_Edges[i].y);
	}

	segments->setLineColor(cColorf(0.0, 1.0, 0.0));
	segments->setUseDisplayList(true);
	segments->setEnabled(true);

	//初始化约束
	for (unsigned i = 0; i < g_numEdges; i++)
	{
		const unsigned i1 = g_Edges[i].x;
		const unsigned i2 = g_Edges[i].y;

		YHDistanceConstraint *c = new YHDistanceConstraint();
		bool res = c->InitConstraint(g_InitPoses, i1, i2);
		if (res)
		{
			g_DistanceConstraints.push_back(c);
		}
	}

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

		// render graphics
		updateGraphics(); //first time register buffer to cuda

		//phys deformation
		cVector3d toolForce(0, -10, 0);
		//std::cout << toolGlobalPos.x() << " " << toolGlobalPos.y() << " " << toolGlobalPos.z() << std::endl;
		PhysDeformation(toolGlobalPos, toolForce, 0.1, applyForceToPhysModel);

		//更新圆球位置
		for (unsigned i = 0; i < g_numSpheres; i++)
		{
			g_pSpheres[i]->setLocalPos(g_SpherePoses[i]);
		}

		//更新约束位置
		for (unsigned i = 0; i < g_numSpheres; i++)
		{
			segments->m_vertices->setLocalPos(i, g_SpherePoses[i]);
			segments->markForUpdate(false);
		}
		

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
	glfwDestroyWindow(g_pWindow->GetGLFWWindow());

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
	//yyh option - apply force on physModel
	else if (a_key == GLFW_KEY_F)
	{
		applyForceToPhysModel = !applyForceToPhysModel;
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

//semi implicit
void SemiImplicitEuler(const double dt, const double invMass, const cVector3d force, cVector3d& position, cVector3d& velocity)
{
	if (invMass != 0.0)
	{
		velocity += force * invMass * dt;
		position += velocity * dt;
	}
}

void VelocityUpdateFirstOrder(const double dt, const double invMass, const cVector3d &position, const cVector3d &oldPosition, cVector3d &velocity)
{
	if (invMass != 0.0)
	{
		velocity = (1.0 / dt) * (position - oldPosition);
	}
}

//phys deformation
void PhysDeformation(cVector3d toolPos, cVector3d toolForce, double dst, bool applyForceOnPhysModel)
{
	/*
	update neighboring

	update distance constraint (add/delete)
	update volume constraint (add/delete)
	...

	do pbd
	1. init velocity
	2. predict pos
	3. use constraints correct pos
	4. update pos & velocity
	*/

	cMaterial matSphereRed;
	matSphereRed.setRedDark();
	cMaterial matSphereWhite;
	matSphereWhite.setWhiteAliceBlue();

	cVector3d force;
	if (applyForceOnPhysModel)
		force = toolForce;
	else
		force = cVector3d(0, 0, 0);

	//计算速度/预测位置
	for (unsigned i = 0; i < g_numSpheres; i++)
	{
		g_SphereOldPoses[i] = g_SpherePoses[i];

		cVector3d& pos = g_SpherePoses[i];
		cVector3d& vel = g_SphereVels[i];
		double invMass = g_SphereInvMasses[i];
		cVector3d forceOnParticle(0, 0, 0);

		if (cDistanceSq(toolPos, pos) <= dst * dst)
		{
			if (applyForceOnPhysModel)
				std::cout << "开启受力&接触到球体" << std::endl;
			forceOnParticle = force;
			g_pSpheres[i]->setMaterial(matSphereRed);
		}
		else
		{
			forceOnParticle = cVector3d(0, 0, 0);
			g_pSpheres[i]->setMaterial(matSphereWhite);
		}
		SemiImplicitEuler(g_dt, invMass, forceOnParticle, pos, vel); 
	}

	//执行约束
	for (unsigned i = 0; i < g_DistanceConstraints.size(); i++)
	{
		g_DistanceConstraints[i]->SolveConstraint(g_SpherePoses.data(), g_SphereInvMasses.data());
	}
	
	////更新速度
	for (unsigned i = 0; i < g_numSpheres; i++)
	{
		cVector3d& pos = g_SpherePoses[i];
		cVector3d& oldPos = g_SphereOldPoses[i];
		cVector3d& vel = g_SphereVels[i];
		double invMass = g_SphereInvMasses[i];
		VelocityUpdateFirstOrder(g_dt, invMass, pos, oldPos, vel);
	}
}

//phys phase change
void PhysPhaseChange(PhysModel* physModel)
{
	/*	compute heat absorbed
	if(particle is SOLID)
	{
	1. heat from air

	2. heat from neighboring
	}

	3. heat from heat source

	compute phase change
	1. compute delT from absorbed heat
	if(particle is SOLID)
	{
	if(T reaches threshold T)
	Gen n SOLID2 particle in SOLID neighboring
	SOLID delete;
	else
	SOLID T += delT
	}
	if(particle is SOLID2)
	{
	if(T reaches melting point)
	{
	compute percentage of liquid r = Qhole / (Lf * density * volume)
	compute mass of free liquid
	SOLID2 mass -= mass of free liquid

	if(SOLID2 mass <=0)
	SOLID2 delete;

	while(mass of free liquid >= mass of a liquid particle)
	{
	Gen a liquid particle in SOLID2 neighboring
	mass of free liquid -= liquid particle mass
	}
	}
	else
	SOLID2 T += delT
	}
	if(particle is LIQUID)
	{
	if(T reaches boiling point)
	{
	compute percentage of gas r = Qhole / (Lf * density * volume)
	compute mass of free gas
	LIQUID mass -= mass of free gas

	if(LIQUID mass <=0)
	LIQUID delete;
	}
	else
	LIQUID T += delT
	}

	*/
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

		// update tool position
		toolGlobalPos = tool->getDeviceGlobalPos();
		toolLocalPos = tool->getDeviceLocalPos();

		// compute transformation from world to tool (haptic device)
		cTransform world_T_tool = tool->getDeviceGlobalTransform();


		if (state == STATE_MOVE_CAMERA && !userSwitch0)
		{
			state = STATE_IDLE;
		}
		else if (state == STATE_IDLE && userSwitch0)
		{
			state = STATE_MOVE_CAMERA;

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