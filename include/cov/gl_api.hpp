#pragma once

#ifdef _WIN32
#include <Windows.h>
#endif
#include <GL/gl.h>

#ifndef APIENTRY
#define APIENTRY
#endif

namespace cov::gl {

bool load();

extern GLuint (APIENTRY* CreateShader)(GLenum);
extern void (APIENTRY* ShaderSource)(GLuint, GLsizei, const char* const*, const GLint*);
extern void (APIENTRY* CompileShader)(GLuint);
extern void (APIENTRY* GetShaderiv)(GLuint, GLenum, GLint*);
extern void (APIENTRY* GetShaderInfoLog)(GLuint, GLsizei, GLsizei*, char*);
extern GLuint (APIENTRY* CreateProgram)();
extern void (APIENTRY* AttachShader)(GLuint, GLuint);
extern void (APIENTRY* LinkProgram)(GLuint);
extern void (APIENTRY* GetProgramiv)(GLuint, GLenum, GLint*);
extern void (APIENTRY* GetProgramInfoLog)(GLuint, GLsizei, GLsizei*, char*);
extern void (APIENTRY* UseProgram)(GLuint);
extern void (APIENTRY* DeleteShader)(GLuint);
extern void (APIENTRY* DeleteProgram)(GLuint);
extern GLint (APIENTRY* GetUniformLocation)(GLuint, const char*);
extern void (APIENTRY* Uniform1i)(GLint, GLint);
extern void (APIENTRY* Uniform1f)(GLint, GLfloat);
extern void (APIENTRY* Uniform2f)(GLint, GLfloat, GLfloat);
extern void (APIENTRY* Uniform3f)(GLint, GLfloat, GLfloat, GLfloat);
extern void (APIENTRY* ActiveTexture)(GLenum);
extern void (APIENTRY* TexImage3D)(GLenum, GLint, GLint, GLsizei, GLsizei, GLsizei,
                                   GLint, GLenum, GLenum, const void*);

} // namespace cov::gl
