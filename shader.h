
class ShaderBase{
 public:
  virtual void shadeVertex(Vec4f vertex) = 0;
  virtual void shadeFragment(Vec3i uvd) = 0; 

};
