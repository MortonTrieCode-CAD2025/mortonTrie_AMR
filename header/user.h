#ifndef USER_H_INCLUDED
#define USER_H_INCLUDED

#include "Constants.h"
#include "Grid_Class.h"
#if (C_MAP_TYPE == 3)
#include "mortonTrieDDF.hpp"
#endif

// template<typename T_Phy_DDF>
#if (C_MAP_TYPE == 3)
// In C_MAP_TYPE=3 (mortonTrie) mode, f[C_Q] and fcol[C_Q] each share a single
// underlying MortonTrieDDF (one trie skeleton + bitmap, fat leaves holding
// `double[BRANCH][C_Q]`). The 27 MortonTrieDDFSlot views expose the same
// `D_map_define<double>`-compatible interface (at / operator[] / insert /
// exists / count) so existing call sites like `df.f[i_q].at(code)` keep
// working without changes.
struct _s_DDF
{
    MortonTrieDDF f_storage;
    MortonTrieDDF fcol_storage;
    MortonTrieDDFSlot f[C_Q];
    MortonTrieDDFSlot fcol[C_Q];

    _s_DDF()
    {
        for (int q = 0; q < C_Q; ++q)
        {
            f[q].bind(&f_storage, static_cast<uint8_t>(q));
            fcol[q].bind(&fcol_storage, static_cast<uint8_t>(q));
        }
    }
};
#else
struct _s_DDF
{
    /**
     * @todo array ddf
     */
    // D_map_define<D_Phy_DDF> f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18;
    D_map_define<D_Phy_DDF> f[C_Q];

    // D_map_define<D_Phy_DDF> fcol0, fcol1, fcol2, fcol3, fcol4, fcol5, fcol6, fcol7, fcol8, fcol9, fcol10, fcol11, fcol12, fcol13, fcol14, fcol15, fcol16, fcol17, fcol18;
    D_map_define<D_Phy_DDF> fcol[C_Q];
};
#endif

// Type alias for one direction of f[C_Q] / fcol[C_Q]. Use this in code that
// previously held `D_map_define<D_Phy_DDF>*` pointers into _s_DDF::f, so the
// type stays correct under all C_MAP_TYPE backends.
#if (C_MAP_TYPE == 3)
using D_DDFSlot = MortonTrieDDFSlot;
#else
using D_DDFSlot = D_map_define<D_Phy_DDF>;
#endif

/**
 * @brief Any variable described with {x, y, z}
 * @date 2023/6/1
 */
template<typename D_T>
struct Coord
{
    D_T x;
    D_T y;
    D_T z;

    Coord() = default;

    Coord(const Coord& c) : x(c.x), y(c.y), z(c.z) {};
    
    Coord(const D_T& construct_x, const D_T& construct_y, const D_T& construct_z)
        : x(construct_x), y(construct_y), z(construct_z)
    {}

    Coord& operator=(const Coord& right)
    {
        x = right.x;
        y = right.y;
        z = right.z;
        return *this;
    }

    Coord& operator-=(const Coord& right)
    {
        x -= right.x;
        y -= right.y;
        z -= right.z;
        return *this;
    }

    Coord& operator+=(const Coord& right)
    {
        x += right.x;
        y += right.y;
        z += right.z;
        return *this;
    }

    Coord& operator/(const D_real& right)
    {
        x = x/right;
        y = y/right;
        z = z/right;
        return *this;
    }

    Coord<D_T> norm()
    {
        D_T magnitude = sqrtf64(x*x + y*y + z*z);
        return Coord<D_T>(x/magnitude, y/magnitude, z/magnitude);
    }
};

template<typename D_T> inline
D_T dot_product(const Coord<D_T>& vec1, const Coord<D_T>& vec2)
{
    return (vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
};

// Cross product
//          | i  j  k|
//  U x V = |u1 u2 u3| = (u2*v3-u3*v2)i + (u3*v1-u1*v3)j + (u1*v2-u2*v1)k
//          |v1 v2 v3|
//
template<typename D_T> inline
Coord<D_T> cross_product(const Coord<D_T>& d1, const Coord<D_T>& d2)
{
    return Coord<D_T>(d1.y * d2.z - d1.z * d2.y, 
                      d1.z * d2.x - d1.x * d2.z, 
                      d1.x * d2.y - d1.y * d2.x);
};

template<typename D_T> inline
Coord<D_T> operator+(const Coord<D_T>& left, const Coord<D_T>& right)
{
    return Coord<D_T>(left.x + right.x,
                      left.y + right.y,
                      left.z + right.z);
};

template<typename D_T> inline
Coord<D_T> operator-(const Coord<D_T>& left, const Coord<D_T>& right)
{
    return Coord<D_T>(left.x - right.x,
                      left.y - right.y,
                      left.z - right.z);
};

using D_vec = Coord<D_real>;
using D_uvw = Coord<D_Phy_Velocity> ;

// template<typename T_Phy_DDF, typename T_Phy_V, typename T_Phy_Rho>
struct User
{
    _s_DDF df;
    // D_uvw *velocity;
    D_map_define<D_uvw> velocity;
    // D_uvw *v_average;
    D_map_define<D_uvw> v_average;
    // D_uvw *v_old;
    D_map_define<D_uvw> v_old;
    // D_Phy_Rho *density;
    D_map_define<D_Phy_Rho> density;
    // D_Phy_Rho *rho_average;
    D_map_define<D_Phy_Rho> rho_average;
};

#endif
