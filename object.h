#pragma once
#ifndef SPHERE_H
#define SPHERE_H
#include <iostream>
#include "vec3.h"
//#include "ray.h"
#include "hittable.h"
#include "aabb.h"
class sphere : public hittable {
public:
	sphere() {}
	sphere(vec3 cen, double r,shared_ptr<material> m) : center(cen), radius(r), mat_ptr(m) {};

	virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
	virtual bool bounding_box(double time0, double time1, aabb& output_box) const;
public:
	vec3 center;
	double radius;
	shared_ptr<material> mat_ptr;
private:
	static void get_sphere_uv(const vec3& p, double& u, double& v) {
		// p: a given point on the sphere of radius one, centered at the origin.
		// u: returned value [0,1] of angle around the Y axis from X=-1.
		// v: returned value [0,1] of angle from Y=-1 to Y=+1.
		//     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
		//     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
		//     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

		auto theta = acos(-p.y());
		auto phi = atan2(-p.z(), p.x()) + pi;

		u = phi / (2 * pi);
		v = theta / pi;
	}
};

//void get_sphere_uv(const vec3& p, double& u, double& v) {
//	auto phi = atan2(p.z(), p.x());
//	auto theta = asin(p.y());
//	u = 1 - (phi + pi) / (2 * pi);
//	v = (theta + pi / 2) / pi;
//}

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
	//p(t)=a+tb
	//b*b*t^2+2t*(a-c)*b+(a-c)*(a-c)-r^2=0
	vec3 oc = r.o() - center; //a-c
	auto a = r.d().length_squared(); //bb
	auto hb = dot(oc, r.d());
	auto c = oc.length_squared() - radius * radius;
	auto disceiminant = hb * hb - a * c;
	//get_sphere_uv((rec.p - center) / radius, rec.u, rec.v);
	if (disceiminant > 0) {
		auto temp = (-hb - sqrt(disceiminant)) / a;
		if (temp< t_max && temp > t_min) {
			rec.t = temp; rec.p = r.at(rec.t); rec.normal = (rec.p - center) / radius;
			vec3 outward_normal = (rec.p - center) / radius;
			rec.set_face_normal(r, outward_normal); rec.mat_ptr = mat_ptr;
			get_sphere_uv((rec.p - center) / radius, rec.u, rec.v);
			return true;
		}
		temp = (-hb + sqrt(disceiminant)) / a;
		if (temp< t_max && temp > t_min) {
			rec.t = temp; rec.p = r.at(rec.t); rec.normal = (rec.p - center) / radius;
			vec3 outward_normal = (rec.p - center) / radius;
			rec.set_face_normal(r, outward_normal); 
			get_sphere_uv((rec.p - center) / radius, rec.u, rec.v);
			rec.mat_ptr = mat_ptr;
			return true;
		}
	}
	return false;
}

bool sphere::bounding_box(double t0, double t1, aabb& output_box) const {
	output_box = aabb(
		center - vec3(radius, radius, radius),
		center + vec3(radius, radius, radius));
	return true;
}

class moving_sphere : public hittable {
public:
	moving_sphere() {}
	moving_sphere(
		vec3 cen0, vec3 cen1, double t0, double t1, double r, shared_ptr<material> m)
		: center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(m)
	{};

	virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
	virtual bool bounding_box(double _time0, double _time1, aabb& output_box) const;
	vec3 center(double time) const;

public:
	vec3 center0, center1;
	double time0, time1;
	double radius;
	shared_ptr<material> mat_ptr;
};

vec3 moving_sphere::center(double time) const {
	return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
}

bool moving_sphere::hit(
	const ray& r, double t_min, double t_max, hit_record& rec) const {
	vec3 oc = r.o() - center(r.time());
	auto a = r.d().length_squared();
	auto half_b = dot(oc, r.d());
	auto c = oc.length_squared() - radius * radius;

	auto discriminant = half_b * half_b - a * c;

	if (discriminant > 0) {
		auto root = sqrt(discriminant);

		auto temp = (-half_b - root) / a;
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.at(rec.t);
			vec3 outward_normal = (rec.p - center(r.time())) / radius;
			rec.set_face_normal(r, outward_normal);
			rec.mat_ptr = mat_ptr;
			return true;
		}

		temp = (-half_b + root) / a;
		if (temp < t_max && temp > t_min) {
			rec.t = temp;
			rec.p = r.at(rec.t);
			vec3 outward_normal = (rec.p - center(r.time())) / radius;
			rec.set_face_normal(r, outward_normal);
			rec.mat_ptr = mat_ptr;
			return true;
		}
	}
	return false;
}

bool moving_sphere::bounding_box(double t0, double t1, aabb& output_box) const {
	aabb box0(
		center(t0) - vec3(radius, radius, radius),
		center(t0) + vec3(radius, radius, radius));
	aabb box1(
		center(t1) - vec3(radius, radius, radius),
		center(t1) + vec3(radius, radius, radius));
	output_box = surrounding_box(box0, box1);
	return true;
}
#endif