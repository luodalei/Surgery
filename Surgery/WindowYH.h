#pragma once

#include <GLAD\glad.h>
#include <GLFW\glfw3.h>
#include <GLFW\glfw3native.h>

//STL
#include <string>
#include <iostream>

using namespace std;


class WindowYH
{
public:
public:
	WindowYH(string name, int width, int height);
	~WindowYH();
	bool ShouldClose();
	void SetKeyCallback(GLFWkeyfun keyFunc);
	void SetCursorPosCallback(GLFWcursorposfun mouseFunc);
	void SetScrollCallback(GLFWscrollfun scroolFunc);
	void SetWindowSizeCallback(GLFWwindowsizefun windowSizeFun);
	void SetMouseButtonCallback(GLFWmousebuttonfun cbfun);
	void SetCharCallback(GLFWcharfun cbfun);
	void SetDropCallback(GLFWdropfun cbfun);
	void SetFrameBufferSizeCallback(GLFWframebuffersizefun cbfun);

	GLuint GetWidth();
	GLuint GetHeight();
	GLFWwindow* GetGLFWWindow();
	void SetWidthAndHeight(GLuint width, GLuint height);

	void SwapBuffers();

private:
	GLuint screenWidth;
	GLuint screenHeight;
	int swapInterval;
	GLFWwindow* m_pWindow;
};

