#pragma once
#include "utils.h"
#include "vec3.h"
#include "ray.h"

class aabb {
	vec3 _min; vec3 _max;
public:
	aabb() {}
	aabb(const vec3& a, const vec3& b) { _min = a; _max = b; }
	vec3 min() const { return _min; }
	vec3 max() const { return _max; }


	virtual bool hit(const ray& r, double t_min, double t_max) const;
	/*bool hit(const ray& r, double tmin, double tmax) const {
		for (int i = 0; i < 3; i++) {
			auto t0 = ffmin((_min[i] - r.o()[i]) / r.d()[i], (_max[i] - r.o()[i]) / r.d()[i]);
			auto t1 = ffmax((_min[i] - r.o()[i]) / r.d()[i], (_max[i] - r.o()[i]) / r.d()[i]);
			tmin = ffmax(t0, tmin);  tmax = ffmin(t1, tmax);
			if (tmax <= tmin) return false;
		}
		return true;
	}*/
};

inline bool aabb::hit(const ray& r, double tmin, double tmax) const {
	for (int a = 0; a < 3; a++) {
		auto invD = 1.0f / r.d()[a];
		auto t0 = (min()[a] - r.o()[a]) * invD;
		auto t1 = (max()[a] - r.o()[a]) * invD;
		if (invD < 0.0f)
			std::swap(t0, t1);
		tmin = t0 > tmin ? t0 : tmin;
		tmax = t1 < tmax ? t1 : tmax;
		if (tmax <= tmin)
			return false;
	}
	return true;
}

aabb surrounding_box(aabb box0, aabb box1) {
	vec3 small(ffmin(box0.min().x(), box1.min().x()),
		ffmin(box0.min().y(), box1.min().y()),
		ffmin(box0.min().z(), box1.min().z()));
	vec3 big(ffmax(box0.max().x(), box1.max().x()),
		ffmax(box0.max().y(), box1.max().y()),
		ffmax(box0.max().z(), box1.max().z()));
	return aabb(small, big);
}

