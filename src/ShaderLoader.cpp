#include "ShaderLoader.h"
#include<GLEW/glew.h>
#include<vector>
#include<string>
#include<iostream>
#include<fstream>


struct ShaderLoader::Impl {
    std::vector<unsigned int> shaders{};
};

ShaderLoader::ShaderLoader() : pimpl{ new Impl{} } {}
ShaderLoader::~ShaderLoader() = default;

void ShaderLoader::addShaderFromCode(const std::string& shaderCode, const unsigned int shaderType, const std::string& name) {
    unsigned int shader_id = glCreateShader(shaderType);

    const char* c_source = shaderCode.c_str();
    glShaderSource(shader_id, 1, &c_source, nullptr);
    glCompileShader(shader_id);

    GLint result;
    glGetShaderiv(shader_id, GL_COMPILE_STATUS, &result);

    if (result == GL_FALSE)
    {
        int length;
        glGetShaderiv(shader_id, GL_INFO_LOG_LENGTH, &length);

        std::unique_ptr<GLchar[]> strInfoLog{ new GLchar[length + 1] };
        glGetShaderInfoLog(shader_id, length, &length, strInfoLog.get());

        fprintf(stderr, "Compilation error in shader \"%s\": %s\n", name.c_str(), strInfoLog.get());
    }

    pimpl->shaders.push_back(shader_id);
}

void ShaderLoader::addShaderFromProjectFilePath(const std::string& fileName, const unsigned int shaderType, const std::string& name) {
    std::ifstream file(fileName);
    if(file.fail()) fprintf(stderr, "Shader file does not exist: %s\n", fileName.c_str());
    const std::string shaderCode((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    this->addShaderFromCode(shaderCode, shaderType, name);
}

void ShaderLoader::attachShaders(const unsigned int programId) {
    for (auto it = pimpl->shaders.begin(); it != pimpl->shaders.end(); it++) {
        glAttachShader(programId, *it);
    }
}
void ShaderLoader::deleteShaders() {
    for (auto it = pimpl->shaders.begin(); it != pimpl->shaders.end(); it++) {
        glDeleteShader(*it);
    }
    pimpl->shaders.clear();
}
