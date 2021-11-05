#include "data.h"

GLuint programId_;
GLuint ssboHandle_;
uint32_t bufferOffset_;

void setProgramId(GLuint programId) {
	programId_ = programId;
}

GLuint programId() {
	return programId_;
}

void setSSBOHandle(GLuint handle) {
	ssboHandle_ = handle;
}

GLuint ssboHandle() {
	return ssboHandle_;
}