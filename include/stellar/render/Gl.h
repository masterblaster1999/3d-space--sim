#pragma once

// Small OpenGL function loader using SDL_GL_GetProcAddress.
// This is intentionally minimal (not a full GLAD replacement).

#include <SDL.h>
#include <SDL_opengl.h>
// SDL provides a portable OpenGL extension header that defines the PFNGL*PROC
// function pointer typedefs on platforms where system OpenGL headers may be
// incomplete (notably Windows' legacy gl.h).
#include <SDL_opengl_glext.h>

// Some header combinations (notably on Windows) can still miss PFN typedefs for
// core OpenGL 1.1 functions (e.g. glGenTextures). The missing typedefs cause
// MSVC to fall back to the old C "default-int" assumption, which then produces
// a huge cascade of confusing syntax errors.
//
// SDL_opengl_glext.h typically provides these, but we add small, guarded
// fallbacks so the project builds reliably across toolchains.

#ifndef PFNGLGENTEXTURESPROC
typedef void (APIENTRYP PFNGLGENTEXTURESPROC)(GLsizei n, GLuint* textures);
#endif

#ifndef PFNGLBINDTEXTUREPROC
typedef void (APIENTRYP PFNGLBINDTEXTUREPROC)(GLenum target, GLuint texture);
#endif

#ifndef PFNGLTEXIMAGE2DPROC
typedef void (APIENTRYP PFNGLTEXIMAGE2DPROC)(GLenum target, GLint level, GLint internalformat,
                                            GLsizei width, GLsizei height, GLint border,
                                            GLenum format, GLenum type, const void* pixels);
#endif

#ifndef PFNGLTEXPARAMETERIPROC
typedef void (APIENTRYP PFNGLTEXPARAMETERIPROC)(GLenum target, GLenum pname, GLint param);
#endif

#ifndef PFNGLDELETETEXTURESPROC
typedef void (APIENTRYP PFNGLDELETETEXTURESPROC)(GLsizei n, const GLuint* textures);
#endif

namespace stellar::render::gl {

bool load();

// Debug helpers
const char* glVersionString();

extern PFNGLCREATESHADERPROC CreateShader;
extern PFNGLSHADERSOURCEPROC ShaderSource;
extern PFNGLCOMPILESHADERPROC CompileShader;
extern PFNGLGETSHADERIVPROC GetShaderiv;
extern PFNGLGETSHADERINFOLOGPROC GetShaderInfoLog;
extern PFNGLDELETESHADERPROC DeleteShader;

extern PFNGLCREATEPROGRAMPROC CreateProgram;
extern PFNGLATTACHSHADERPROC AttachShader;
extern PFNGLLINKPROGRAMPROC LinkProgram;
extern PFNGLGETPROGRAMIVPROC GetProgramiv;
extern PFNGLGETPROGRAMINFOLOGPROC GetProgramInfoLog;
extern PFNGLUSEPROGRAMPROC UseProgram;
extern PFNGLDELETEPROGRAMPROC DeleteProgram;

extern PFNGLGETUNIFORMLOCATIONPROC GetUniformLocation;
extern PFNGLUNIFORM1IPROC Uniform1i;
extern PFNGLUNIFORM1FPROC Uniform1f;
extern PFNGLUNIFORM3FPROC Uniform3f;
extern PFNGLUNIFORMMATRIX4FVPROC UniformMatrix4fv;

extern PFNGLGENVERTEXARRAYSPROC GenVertexArrays;
extern PFNGLBINDVERTEXARRAYPROC BindVertexArray;
extern PFNGLDELETEVERTEXARRAYSPROC DeleteVertexArrays;

extern PFNGLGENBUFFERSPROC GenBuffers;
extern PFNGLBINDBUFFERPROC BindBuffer;
extern PFNGLBUFFERDATAPROC BufferData;
extern PFNGLBUFFERSUBDATAPROC BufferSubData;
extern PFNGLDELETEBUFFERSPROC DeleteBuffers;

extern PFNGLENABLEVERTEXATTRIBARRAYPROC EnableVertexAttribArray;
extern PFNGLVERTEXATTRIBPOINTERPROC VertexAttribPointer;
extern PFNGLVERTEXATTRIBDIVISORPROC VertexAttribDivisor;

extern PFNGLDRAWARRAYSINSTANCEDPROC DrawArraysInstanced;
extern PFNGLDRAWELEMENTSINSTANCEDPROC DrawElementsInstanced;

extern PFNGLACTIVETEXTUREPROC ActiveTexture;
extern PFNGLGENTEXTURESPROC GenTextures;
extern PFNGLBINDTEXTUREPROC BindTexture;
extern PFNGLTEXIMAGE2DPROC TexImage2D;
extern PFNGLTEXPARAMETERIPROC TexParameteri;
extern PFNGLGENERATEMIPMAPPROC GenerateMipmap;
extern PFNGLDELETETEXTURESPROC DeleteTextures;

} // namespace stellar::render::gl
