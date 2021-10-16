#pragma once

#include <GL/glew.h>
#include <stdint.h>

void setProgramId(GLuint programId);
GLuint programId();

void setBufferWriteOffset(uint32_t offset);
void setSSBOHandle(GLuint handle);

uint32_t bufferWriteOffset();
GLuint ssboHandle();

void swapBuffers();