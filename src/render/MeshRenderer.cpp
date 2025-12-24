#include "stellar/render/MeshRenderer.h"

#include "stellar/render/Gl.h"

#include <cstring>

namespace stellar::render {

MeshRenderer::~MeshRenderer() {
  if (instanceVbo_) {
    gl::DeleteBuffers(1, &instanceVbo_);
    instanceVbo_ = 0;
  }
}

static const char* kVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec2 aUv;

layout(location=3) in vec3 iPos;
layout(location=4) in vec3 iScale;
layout(location=5) in vec4 iQuat; // (x,y,z,w)
layout(location=6) in vec3 iColor;

uniform mat4 uView;
uniform mat4 uProj;

out vec2 vUv;
out vec3 vColor;
out vec3 vNrm;

vec3 quatRotate(vec4 q, vec3 v) {
  // q assumed normalized
  vec3 t = 2.0 * cross(q.xyz, v);
  return v + q.w * t + cross(q.xyz, t);
}

void main() {
  vec3 local = aPos * iScale;
  vec3 pos = quatRotate(iQuat, local) + iPos;
  gl_Position = uProj * uView * vec4(pos, 1.0);
  vUv = aUv;
  vColor = iColor;
  vNrm = quatRotate(iQuat, aNrm);
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
in vec2 vUv;
in vec3 vColor;
in vec3 vNrm;

uniform sampler2D uTex;

out vec4 FragColor;

void main() {
  vec3 n = normalize(vNrm);
  vec3 l = normalize(vec3(0.4, 0.8, 0.2));
  float diff = max(dot(n,l), 0.0);
  vec3 tex = texture(uTex, vUv).rgb;
  vec3 col = tex * vColor * (0.35 + 0.65*diff);
  FragColor = vec4(col, 1.0);
}
)GLSL";

bool MeshRenderer::init(std::string* outError) {
  if (!shader_.build(kVS, kFS, outError)) return false;

  gl::GenBuffers(1, &instanceVbo_);

  // default identity matrices
  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void MeshRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void MeshRenderer::drawInstances(const std::vector<InstanceData>& instances) {
  if (!mesh_ || instances.empty()) return;

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform1i("uTex", 0);

  if (tex_) tex_->bind(0);

  // Bind mesh VAO and configure instance attributes
  mesh_->bind();

  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(instances.size() * sizeof(InstanceData)),
                 instances.data(),
                 GL_DYNAMIC_DRAW);

  // location 3: vec3 position
  gl::EnableVertexAttribArray(3);
  gl::VertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, px));
  gl::VertexAttribDivisor(3, 1);

  // location 4: vec3 scale
  gl::EnableVertexAttribArray(4);
  gl::VertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, sx));
  gl::VertexAttribDivisor(4, 1);

  // location 5: vec4 quaternion
  gl::EnableVertexAttribArray(5);
  gl::VertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, qx));
  gl::VertexAttribDivisor(5, 1);

  // location 6: vec3 color
  gl::EnableVertexAttribArray(6);
  gl::VertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, cr));
  gl::VertexAttribDivisor(6, 1);

  mesh_->drawInstanced(static_cast<std::uint32_t>(instances.size()));
}

} // namespace stellar::render
