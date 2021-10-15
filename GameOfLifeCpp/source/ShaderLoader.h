#pragma once

#include<string>
#include<memory>

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
	void addShaderFromProjectFileName(const std::string& shaderCode, const unsigned int shaderType, const std::string& name = "");

	void attachShaders(const unsigned int programId);
	void deleteShaders();
};

