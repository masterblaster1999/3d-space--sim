#include "stellar/render/Gl.h"

#include "stellar/core/Log.h"

namespace stellar::render::gl {

PFNGLCREATESHADERPROC CreateShader = nullptr;
PFNGLSHADERSOURCEPROC ShaderSource = nullptr;
PFNGLCOMPILESHADERPROC CompileShader = nullptr;
PFNGLGETSHADERIVPROC GetShaderiv = nullptr;
PFNGLGETSHADERINFOLOGPROC GetShaderInfoLog = nullptr;
PFNGLDELETESHADERPROC DeleteShader = nullptr;

PFNGLCREATEPROGRAMPROC CreateProgram = nullptr;
PFNGLATTACHSHADERPROC AttachShader = nullptr;
PFNGLLINKPROGRAMPROC LinkProgram = nullptr;
PFNGLGETPROGRAMIVPROC GetProgramiv = nullptr;
PFNGLGETPROGRAMINFOLOGPROC GetProgramInfoLog = nullptr;
PFNGLUSEPROGRAMPROC UseProgram = nullptr;
PFNGLDELETEPROGRAMPROC DeleteProgram = nullptr;

PFNGLGETUNIFORMLOCATIONPROC GetUniformLocation = nullptr;
PFNGLUNIFORM1IPROC Uniform1i = nullptr;
PFNGLUNIFORM1FPROC Uniform1f = nullptr;
PFNGLUNIFORM3FPROC Uniform3f = nullptr;
PFNGLUNIFORMMATRIX4FVPROC UniformMatrix4fv = nullptr;

PFNGLGENVERTEXARRAYSPROC GenVertexArrays = nullptr;
PFNGLBINDVERTEXARRAYPROC BindVertexArray = nullptr;
PFNGLDELETEVERTEXARRAYSPROC DeleteVertexArrays = nullptr;

PFNGLGENBUFFERSPROC GenBuffers = nullptr;
PFNGLBINDBUFFERPROC BindBuffer = nullptr;
PFNGLBUFFERDATAPROC BufferData = nullptr;
PFNGLBUFFERSUBDATAPROC BufferSubData = nullptr;
PFNGLDELETEBUFFERSPROC DeleteBuffers = nullptr;

PFNGLENABLEVERTEXATTRIBARRAYPROC EnableVertexAttribArray = nullptr;
PFNGLVERTEXATTRIBPOINTERPROC VertexAttribPointer = nullptr;
PFNGLVERTEXATTRIBDIVISORPROC VertexAttribDivisor = nullptr;

PFNGLDRAWARRAYSINSTANCEDPROC DrawArraysInstanced = nullptr;
PFNGLDRAWELEMENTSINSTANCEDPROC DrawElementsInstanced = nullptr;

PFNGLACTIVETEXTUREPROC ActiveTexture = nullptr;
PFNGLGENTEXTURESPROC GenTextures = nullptr;
PFNGLBINDTEXTUREPROC BindTexture = nullptr;
PFNGLTEXIMAGE2DPROC TexImage2D = nullptr;
PFNGLTEXPARAMETERIPROC TexParameteri = nullptr;
PFNGLGENERATEMIPMAPPROC GenerateMipmap = nullptr;
PFNGLDELETETEXTURESPROC DeleteTextures = nullptr;

template <class T>
static T loadProc(const char* name) {
  return reinterpret_cast<T>(SDL_GL_GetProcAddress(name));
}

bool load() {
  CreateShader = loadProc<PFNGLCREATESHADERPROC>("glCreateShader");
  ShaderSource = loadProc<PFNGLSHADERSOURCEPROC>("glShaderSource");
  CompileShader = loadProc<PFNGLCOMPILESHADERPROC>("glCompileShader");
  GetShaderiv = loadProc<PFNGLGETSHADERIVPROC>("glGetShaderiv");
  GetShaderInfoLog = loadProc<PFNGLGETSHADERINFOLOGPROC>("glGetShaderInfoLog");
  DeleteShader = loadProc<PFNGLDELETESHADERPROC>("glDeleteShader");

  CreateProgram = loadProc<PFNGLCREATEPROGRAMPROC>("glCreateProgram");
  AttachShader = loadProc<PFNGLATTACHSHADERPROC>("glAttachShader");
  LinkProgram = loadProc<PFNGLLINKPROGRAMPROC>("glLinkProgram");
  GetProgramiv = loadProc<PFNGLGETPROGRAMIVPROC>("glGetProgramiv");
  GetProgramInfoLog = loadProc<PFNGLGETPROGRAMINFOLOGPROC>("glGetProgramInfoLog");
  UseProgram = loadProc<PFNGLUSEPROGRAMPROC>("glUseProgram");
  DeleteProgram = loadProc<PFNGLDELETEPROGRAMPROC>("glDeleteProgram");

  GetUniformLocation = loadProc<PFNGLGETUNIFORMLOCATIONPROC>("glGetUniformLocation");
  Uniform1i = loadProc<PFNGLUNIFORM1IPROC>("glUniform1i");
  Uniform1f = loadProc<PFNGLUNIFORM1FPROC>("glUniform1f");
  Uniform3f = loadProc<PFNGLUNIFORM3FPROC>("glUniform3f");
  UniformMatrix4fv = loadProc<PFNGLUNIFORMMATRIX4FVPROC>("glUniformMatrix4fv");

  GenVertexArrays = loadProc<PFNGLGENVERTEXARRAYSPROC>("glGenVertexArrays");
  BindVertexArray = loadProc<PFNGLBINDVERTEXARRAYPROC>("glBindVertexArray");
  DeleteVertexArrays = loadProc<PFNGLDELETEVERTEXARRAYSPROC>("glDeleteVertexArrays");

  GenBuffers = loadProc<PFNGLGENBUFFERSPROC>("glGenBuffers");
  BindBuffer = loadProc<PFNGLBINDBUFFERPROC>("glBindBuffer");
  BufferData = loadProc<PFNGLBUFFERDATAPROC>("glBufferData");
  BufferSubData = loadProc<PFNGLBUFFERSUBDATAPROC>("glBufferSubData");
  DeleteBuffers = loadProc<PFNGLDELETEBUFFERSPROC>("glDeleteBuffers");

  EnableVertexAttribArray = loadProc<PFNGLENABLEVERTEXATTRIBARRAYPROC>("glEnableVertexAttribArray");
  VertexAttribPointer = loadProc<PFNGLVERTEXATTRIBPOINTERPROC>("glVertexAttribPointer");
  VertexAttribDivisor = loadProc<PFNGLVERTEXATTRIBDIVISORPROC>("glVertexAttribDivisor");

  DrawArraysInstanced = loadProc<PFNGLDRAWARRAYSINSTANCEDPROC>("glDrawArraysInstanced");
  DrawElementsInstanced = loadProc<PFNGLDRAWELEMENTSINSTANCEDPROC>("glDrawElementsInstanced");

  ActiveTexture = loadProc<PFNGLACTIVETEXTUREPROC>("glActiveTexture");
  if (!ActiveTexture) {
    // Some drivers expose this via ARB extension entry points.
    ActiveTexture = loadProc<PFNGLACTIVETEXTUREPROC>("glActiveTextureARB");
  }
  GenTextures = loadProc<PFNGLGENTEXTURESPROC>("glGenTextures");
  BindTexture = loadProc<PFNGLBINDTEXTUREPROC>("glBindTexture");
  TexImage2D = loadProc<PFNGLTEXIMAGE2DPROC>("glTexImage2D");
  TexParameteri = loadProc<PFNGLTEXPARAMETERIPROC>("glTexParameteri");
  GenerateMipmap = loadProc<PFNGLGENERATEMIPMAPPROC>("glGenerateMipmap");
  if (!GenerateMipmap) {
    // Some older drivers expose mipmap generation as an extension entry point.
    GenerateMipmap = loadProc<PFNGLGENERATEMIPMAPPROC>("glGenerateMipmapEXT");
  }
  DeleteTextures = loadProc<PFNGLDELETETEXTURESPROC>("glDeleteTextures");

  // Fallback for entry points that may be exported directly by the OpenGL library
  // rather than returned by SDL_GL_GetProcAddress (common on Windows for GL 1.1/1.3 funcs).
  if (!GenTextures)   GenTextures   = reinterpret_cast<PFNGLGENTEXTURESPROC>(&::glGenTextures);
  if (!BindTexture)   BindTexture   = reinterpret_cast<PFNGLBINDTEXTUREPROC>(&::glBindTexture);
  if (!TexImage2D)    TexImage2D    = reinterpret_cast<PFNGLTEXIMAGE2DPROC>(&::glTexImage2D);
  if (!TexParameteri) TexParameteri = reinterpret_cast<PFNGLTEXPARAMETERIPROC>(&::glTexParameteri);
  if (!DeleteTextures)DeleteTextures= reinterpret_cast<PFNGLDELETETEXTURESPROC>(&::glDeleteTextures);

  // Note: Do NOT fall back to &::glActiveTexture or &::glGenerateMipmap on Windows.
  // The system OpenGL import library (opengl32.dll) only guarantees prototypes/exports for
  // OpenGL 1.1, so referencing these symbols can fail to compile even if the runtime supports
  // them. These must be loaded via SDL_GL_GetProcAddress.

  const bool ok =
      CreateShader && ShaderSource && CompileShader && GetShaderiv && GetShaderInfoLog && DeleteShader &&
      CreateProgram && AttachShader && LinkProgram && GetProgramiv && GetProgramInfoLog && UseProgram && DeleteProgram &&
      GetUniformLocation && Uniform1i && Uniform1f && Uniform3f && UniformMatrix4fv &&
      GenVertexArrays && BindVertexArray && DeleteVertexArrays &&
      GenBuffers && BindBuffer && BufferData && BufferSubData && DeleteBuffers &&
      EnableVertexAttribArray && VertexAttribPointer && VertexAttribDivisor &&
      DrawArraysInstanced && DrawElementsInstanced &&
      ActiveTexture && GenTextures && BindTexture && TexImage2D && TexParameteri && GenerateMipmap && DeleteTextures;

  if (!ok) {
    stellar::core::log(stellar::core::LogLevel::Error, "OpenGL loader: missing required functions (context too old?)");
  }

  return ok;
}

const char* glVersionString() {
  const auto* s = glGetString(GL_VERSION);
  return s ? reinterpret_cast<const char*>(s) : "(null)";
}

} // namespace stellar::render::gl
