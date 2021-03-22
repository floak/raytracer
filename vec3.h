#pragma once

#ifndef __VEC3_H__
#define __VEC3_H__
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "utils.h"
#include <cmath>

//using namespace std;

class vec3 {
public:
	//ld x, y, z;
	//Vec3(ld x_ = 0, ld y_ = 0, ld z_ = 0) : x(x_), y(y_), z(z_) {} //??

	//Vec3 operator - () const { return Vec3(-x, -y, -z); }
	//Vec3 operator + (const Vec3& a) const { return Vec3(x + a.x, y + a.y, z + a.z); }
	//Vec3 operator - (const Vec3& a) const { return Vec3(x - a.x, y - a.y, z - a.z); }
	//Vec3 operator*(double a) const { return Vec3(x * a, y * a, z * a); }
	//Vec3 mult(const Vec3& a) const { return Vec3(x * a.x, y * a.y, z * a.z); }
	//Vec3& norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
	//double dot(const Vec3& a) const { return x * a.x + y * a.y + z * a.z; }          //内积
	//Vec3 operator%(Vec3& a) const { return Vec3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x) }  //叉积
	vec3() : e{ 0,0,0 } {}
	vec3(double e0, double e1, double e2) :e{ e0,e1,e2 } {}
	//vec3(double x,double y,double z):e{x,y,z}{}
	double x() const { return e[0]; } double y() const { return e[1]; } double z() const { return e[2]; }
	//vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
	double operator[](int i) const { return e[i]; }
	double& operator[](int i) { return e[i]; }
	vec3 operator - () const { return vec3(-e[0], -e[1], -e[2]); }
	vec3 operator += (const vec3& a) { vec3(e[0] += a.e[0], e[1] += a.e[1], e[2] += a.e[2]); return *this; }
	vec3 operator*=(const double a) { vec3(e[0] *= a, e[1] *= a, e[2] *= a); return *this; }
	vec3 operator/=(const double a) { return *this *= 1 / a; }
	double length() const {return sqrt(e[0]*e[0]+ e[1] * e[1]+ e[2] * e[2]);}

	double length_squared() const {return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];}

	//inline static vec3 random() {
	//	return vec3(random_double(), random_double(), random_double());
	//}

	bool near_zero() const {
		// Return true if the vector is close to zero in all dimensions.
		const auto s = 1e-8;
		return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
	}

	void write_color(std::ostream& out,int samples_per_pixel) {
		auto scale = 1.0 / samples_per_pixel;  auto r = sqrt(scale * e[0]);auto g = sqrt(scale * e[1]);auto b = sqrt(scale * e[2]);
		
		out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
			<< static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
			<< static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
	}
public:
	double e[3];
	/*vec3 operator + (ld p) const { return vec3(x + p, y + p, z + p); }
	vec3 operator - (ld p) const { return vec3(x - p, y - p, z - p); }
	vec3 operator * (ld p) const { return vec3(x * p, y * p, z * p); }
	vec3 operator / (ld p) const { return vec3(x / p, y / p, z / p); }

	bool operator*/

};

inline std::ostream& operator<<(std::ostream& out, const vec3& v) {return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];}
inline vec3 operator+(const vec3& u, const vec3& v) {return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);}
inline vec3 operator-(const vec3& u, const vec3& v) {return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);}
inline vec3 operator*(const vec3& u, const vec3& v) {return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);}
inline vec3 operator*(double t, const vec3& v) {return vec3(t * v.e[0], t * v.e[1], t * v.e[2]);}
inline vec3 operator*(const vec3& v, double t) {return t * v;}
inline vec3 operator/(vec3 v, double t) {return (1 / t) * v;}

inline double dot(const vec3 &u, const vec3 &v) { return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2]; }
inline vec3 cross(const vec3& u, const vec3& v) { return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1], u.e[2] * v.e[0] - u.e[0] * v.e[2], u.e[0] * v.e[1] - u.e[1] * v.e[0]); }
inline vec3 unit_vector(vec3 v) {return v / v.length();}
//
//vec3 random_in_unit_sphere() {
//	while (true) {
//		auto p = vec3::random(-1, 1);
//		if (p.length_squared() >= 1) continue;
//		return p;
//	}
//}
vec3 reflect(const vec3& v, const vec3& n) { return v - 2 * dot(v, n) * n; }

const double random(){
	static std::mt19937 mt;
	static std::uniform_real_distribution<double> rtrand;
	return rtrand(mt);
}

const double random_double(double min, double max) { return min + (max - min) * random(); }

const int random_int(int min, int max) {
	// Returns a random integer in [min,max].
	return static_cast<int>(random_double(min, max + 1));
}

const vec3 random_unit_sphere() {
	vec3 p;
	do {
		p = 2.0 * vec3(random(), random(), random()) - vec3(1, 1, 1);
	} while (dot(p, p) >= 1.0);
	return p;
}

const vec3 random_unit_vector() {
	auto a = random() + 2.0 * pi;
	auto z = 2.0 * random() - 1.0;
	auto r = sqrt(1 - z * z);
	return vec3(r * cos(a), r * sin(a), z);
}

vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
	auto cos_theta = dot(-uv, n);
	vec3 r_out_parallel = etai_over_etat * (uv + cos_theta * n);
	vec3 r_out_perp = -sqrt(1.0 - r_out_parallel.length_squared()) * n;
	return r_out_parallel + r_out_perp;
}

//单位圆盘随机点函数
const vec3 random_unit_disk() {
	vec3 p; do { p = 2.0 * vec3(random(), random(), 0) - vec3(1, 1, 0); } while (dot(p, p) >= 1.0);
	return p;
}
#endif // !__VEC3_H__