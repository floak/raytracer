#pragma once

#ifndef __UTILS_H__
#define __UTILS_H__

//#include <bits/stdc++.h>
#include <omp.h>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <functional>
#include <random>
#include "vec3.h"


#define stds std::

using namespace std;
using std::shared_ptr;
using std::make_shared;
////一个读写图片的库stb图像解码库，
////image loading/decoding from file/memory: JPG, PNG, TGA, BMP, PSD, GIF, HDR, PIC
//#ifndef STB_IMAGE_IMPLEMENTATION
//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"
//#endif // STB_IMAGE_IMPLEMENTATION   


//typedef是定义了一种类型的新别名，不是简单的字符串替换，所以它比宏来得稳健(可移植性
//typedef double ld;
//const ld PI = acos(-1); //反余弦
const double pi = 3.14159265358979f;
const double eps = 1e-6;    //
//const ld INF = 1 << 20;
const double infinity = std::numeric_limits<double>::infinity();

//const ld min_p[3] = { 100, 100, 100 };
//const ld max_p[3] = { 10000, 10000, 10000 };

//enum:枚举量，
enum Refl_t { DIFF, SPEC, REFR }; //反射项reflection type

inline double degrees_to_radians(double degrees) { return degrees * pi / 180;}

inline double ffmin(double a, double b) { return a <= b ? a : b; }
inline double ffmax(double a, double b) { return a >= b ? a : b; }

static std::uniform_real_distribution<double> distribution;
static std::mt19937 generator;

//inline double random_double() {
//	static std::uniform_real_distribution<double> distribution;
//    static std::mt19937 generator;
//    static std::function<double()> rand_generator = std::bind(distribution, generator);
//    return rand_generator();
//}

inline double clamp(double x, double min, double max) { if (x < min) return min; if (x > max) return max;
    return x;
}

// Common Headers

//#include "ray.h"
//#include "vec3.h"

#endif // __UTILS_H__

