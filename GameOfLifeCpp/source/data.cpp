#include "data.h"

GLuint programId_;
GLuint ssboHandle_;
uint32_t bufferOffset_;
bool isOffset = false;

void setProgramId(GLuint programId) {
	programId_ = programId;
}

GLuint programId() {
	return programId_;
}

void setBufferWriteOffset(uint32_t offset) {
	bufferOffset_ = offset;
}

uint32_t bufferWriteOffset() {
	return bufferOffset_ * isOffset;
}

void swapBuffers() {
	isOffset = !isOffset;
	GLint is2ndBufferP = glGetUniformLocation(programId_, "is2ndBuffer");
	glUniform1i(is2ndBufferP, isOffset);
}

void setSSBOHandle(GLuint handle) {
	ssboHandle_ = handle;
}

GLuint ssboHandle() {
	return ssboHandle_;
}