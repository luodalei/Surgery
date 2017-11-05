#include "WindowYH.h"



WindowYH::WindowYH(string name, int width, int height) :
screenWidth(width),
screenHeight(height),
swapInterval(1)
{
	// Init GLFW
	if (!glfwInit())
	{
		cout << "failed initialization" << endl;
	}
	glfwSetTime(0);

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_STEREO, GL_FALSE);
	/*glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	glfwWindowHint(GLFW_SAMPLES, 0);
	glfwWindowHint(GLFW_RED_BITS, 8);
	glfwWindowHint(GLFW_GREEN_BITS, 8);
	glfwWindowHint(GLFW_BLUE_BITS, 8);
	glfwWindowHint(GLFW_ALPHA_BITS, 8);
	glfwWindowHint(GLFW_STENCIL_BITS, 8);
	glfwWindowHint(GLFW_DEPTH_BITS, 24);
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
*/
	m_pWindow = glfwCreateWindow(screenWidth, screenHeight, name.c_str(), nullptr, nullptr); // Windowed
	if (m_pWindow == nullptr)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
	}
	glfwMakeContextCurrent(m_pWindow);
	glfwSwapInterval(swapInterval);

	

#if defined(NANOGUI_GLAD)
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
		throw std::runtime_error("Could not initialize GLAD!");
	glGetError(); // pull and ignore unhandled errors like GL_INVALID_ENUM
#endif
}


WindowYH::~WindowYH()
{
}

bool WindowYH::ShouldClose()
{
	return glfwWindowShouldClose(m_pWindow);
}

void WindowYH::SetKeyCallback(GLFWkeyfun keyFunc)
{
	glfwSetKeyCallback(m_pWindow, keyFunc);
}

void WindowYH::SetCursorPosCallback(GLFWcursorposfun mouseFunc)
{
	glfwSetCursorPosCallback(m_pWindow, mouseFunc);
}

void WindowYH::SetScrollCallback(GLFWscrollfun scroolFunc)
{
	glfwSetScrollCallback(m_pWindow, scroolFunc);
}

void WindowYH::SetWindowSizeCallback(GLFWwindowsizefun windowSizeFun)
{
	glfwSetWindowSizeCallback(m_pWindow, windowSizeFun);
}

void WindowYH::SetMouseButtonCallback(GLFWmousebuttonfun cbfun)
{
	glfwSetMouseButtonCallback(m_pWindow, cbfun);
}

void WindowYH::SetCharCallback(GLFWcharfun cbfun)
{
	glfwSetCharCallback(m_pWindow, cbfun);
}

void WindowYH::SetDropCallback(GLFWdropfun cbfun)
{
	glfwSetDropCallback(m_pWindow, cbfun);
}

void WindowYH::SetFrameBufferSizeCallback(GLFWframebuffersizefun cbfun)
{
	glfwSetFramebufferSizeCallback(m_pWindow, cbfun);
}

GLuint WindowYH::GetWidth()
{
	return screenWidth;
}

GLuint WindowYH::GetHeight()
{
	return screenHeight;
}

GLFWwindow* WindowYH::GetGLFWWindow()
{
	return m_pWindow;
}

void WindowYH::SetWidthAndHeight(GLuint width, GLuint height)
{
	screenWidth = width;
	screenHeight = height;
}

void WindowYH::SwapBuffers()
{
	glfwSwapBuffers(m_pWindow);
}