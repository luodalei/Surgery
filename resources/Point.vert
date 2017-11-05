#version 120

//attribute vec3 aPosition;
//in vec3 aPosition;             //new version
//attribute vec3 aNormal;
//attribute vec3 aTexCoord;
attribute vec4 aColor;
//in vec4 aColor;                //new version
//attribute vec3 aTangent;
//attribute vec3 aBitangent;

// vertex position and normal in eye coordinate space
//varying vec4 vPosition;
varying vec4 vColor;
//out vec4 vColor;               //new version

//----------------------------------------------------------------------
// Main vertex shader code.
//----------------------------------------------------------------------

void main(void)
{
    // pass along a transformed vertex position, normal, and texture
    //vPosition = gl_ModelViewMatrix * gl_Vertex;
    vColor = aColor;
    
    // fixed function vertex transform
    gl_Position = ftransform();
}