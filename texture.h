#pragma once
#include "utils.h"
#include "vec3.h"
#include "ray.h"


class texture {
public:
	virtual vec3 value(double u, double v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public:
	constant_texture(){}
	constant_texture(vec3 c): color(c){}

	constant_texture(double red,double green,double blue):constant_texture(vec3(red,green,blue)){}
	virtual vec3 value(double u, double v, const vec3& p) const { return color; }
private:
	vec3 color;
};

class checker_texture :public texture {
public:
	shared_ptr<texture> odd;
	shared_ptr<texture> even;
	checker_texture(){}
	checker_texture(shared_ptr<texture> t0,shared_ptr<texture> t1):even(t0),odd(t1){}
	checker_texture(vec3 c1, vec3 c2)
		: even(make_shared<constant_texture>(c1)), odd(make_shared<constant_texture>(c2)) {}

	virtual vec3 value(double u, double v, const vec3& p) const {
		auto sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());
		if (sines < 0)
			return odd->value(u, v, p);
		else
			return even->value(u, v, p);

	}
};