#pragma once

#ifndef __RAY_H__
#define __RAY_H__ 

#include "utils.h"
#include "vec3.h"


class ray {
public:
    ray() {}
    ray(const vec3& o, const vec3& d,double time=0.0) : orig(o), dir(d) ,tm(time){}

    vec3 o() const { return orig; }
    vec3 d() const { return dir; }
    double time() const { return tm; }

    vec3 at(double t) const { return orig + t * dir; }
public:
    vec3 orig;
    vec3 dir;
    double tm;
};
#endif // __RAY_H__

