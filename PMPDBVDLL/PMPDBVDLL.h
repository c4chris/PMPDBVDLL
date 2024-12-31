#pragma once
extern "C" {
#ifdef PMPDBVDLL_EXPORTS
#define PMPDBVDLL_API __declspec(dllexport)
#else
#define PMPDBVDLL_API __declspec(dllimport)
#endif

//
//  Use this file to import your target's public headers that you would like to expose to Swift.
//

#ifdef _MSC_VER
#define _Nonnull
#define _Nullable
#endif

typedef struct Space_Point Point3D;
struct Space_Point
{
  double x, y, z;
};

#if defined(_MSC_VER) || defined(__ANDROID__)
typedef bool Boolean;
#else
typedef _Bool Boolean;
#endif
typedef double DBL;

typedef struct FSpace_Point FPoint3D;
struct FSpace_Point
{
  float x, y, z;
};

extern const char * _Nullable bundleRsrcDir;
extern const char * _Nullable downloadDir;
extern const char * _Nullable tempDir;

struct SwiftCallbacks
{
  void (*_Nonnull GL_atm)(FPoint3D A, short layer, unsigned short group, unsigned short atom, char chain, long color);
  void (*_Nonnull GL_bnd)(FPoint3D A, FPoint3D B, short layer, unsigned short group, unsigned short atomA, unsigned short atomB, char chain, long color);
  void (*_Nonnull GL_Atom_Color)(const char * _Nonnull name, double red, double green, double blue, Boolean isMetallic);
  void (*_Nonnull GL_All_Invisible)(short layer);
  void (*_Nonnull GL_ForceRebuild)(short layer);
  void (*_Nonnull GL_New_Shape)(void);
  void (*_Nonnull GL_Add_Vertex)(Point3D V, Point3D N, unsigned int color);
  void (*_Nonnull GL_Draw_Shape)(short layer, char type);
};
typedef struct SwiftCallbacks SwiftCallbacks;

PMPDBVDLL_API void SwiftCBSetup(const SwiftCallbacks *_Nonnull callbacks);
PMPDBVDLL_API void GL_init_Atom_Colors(void);
PMPDBVDLL_API void GLRender(void);
PMPDBVDLL_API short InitGlobals(void);
PMPDBVDLL_API void doPDBinput(const char * _Nullable, const char * _Nullable);

extern int (*_Nullable printfDelegate)(const char* fmt, va_list va);
}
