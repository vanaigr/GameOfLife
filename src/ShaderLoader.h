#pragma once

#include<string>
#include<memory>
#include<GL/glew.h>

class ShaderLoader
{
private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
public:
    ShaderLoader();
    ~ShaderLoader();

public:
    void addShaderFromCode(const std::string& shaderCode, const unsigned int shaderType, const std::string& name = "");
    void addShaderFromProjectFilePath(const std::string& shaderCode, const unsigned int shaderType, const std::string& name = "");

    void addScreenSizeTriangleStripVertexShader(const std::string& name = "");
    void addScreenCoordFragmentShader(const std::string& name = "");

    void attachShaders(const unsigned int programId);
    void deleteShaders();

    void linkProgram(const unsigned int programId) {
        glLinkProgram(programId);

        GLint linkStatus;
        glGetProgramiv(programId, GL_LINK_STATUS, &linkStatus);
        if (linkStatus == GL_FALSE) {
            GLint logLength;
            glGetProgramiv(programId, GL_INFO_LOG_LENGTH, &logLength);

            char *logMessage = (char *)malloc(logLength);
            glGetProgramInfoLog(programId, logLength, NULL, logMessage);
            fprintf(stderr, "Shader Program Linking Error:\n%s\n", logMessage);
            free(logMessage);
        }
    }
};


inline void ShaderLoader::addScreenSizeTriangleStripVertexShader(const std::string& name) {
    this->addShaderFromCode(
        "\n#version 430"
        "\nprecision mediump float;"
        "\nvoid main(void){"
        "\ngl_Position = vec4("
        "\n    2 * (gl_VertexID / 2) - 1,"
        "\n    2 * (gl_VertexID % 2) - 1,"
        "\n    0.0,"
        "\n    1.0);"
        "\n}"
        ,
        GL_VERTEX_SHADER,
        name
    );
}
inline void ShaderLoader::addScreenCoordFragmentShader(const std::string& name) {
    this->addShaderFromCode(
        "\n#version 430"
        "\nprecision mediump float;"
        "\nlayout(origin_upper_left) in vec4 gl_FragCoord;"
        "\nout vec4 color;"
        "\nvoid main(void){"
        "\ncolor = gl_FragCoord;"
        "\n}"
        ,
        GL_FRAGMENT_SHADER,
        name
    );
}
