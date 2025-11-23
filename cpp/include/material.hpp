#pragma once
typedef unsigned int SAINT;
typedef double SADB;
class conductivity
{
public:
  conductivity()
  {
    k_x = -1;
    k_y = -1;
    k_z = -1;
  }
  SADB k_x;
  SADB k_y;
  SADB k_z;
};
class material
{
public:
  material()
  {
    density = -1.0;
  };
  SADB density;
  conductivity k;
};
