#include <iostream>
#include "vec3.h"
#include "ray.h"
#include "object.h"
#include "utils.h"
#include "hittable_list.h"
#include "camera.h"
#include <random>
//#include <math.h>
//#include <cmath>
#include "material.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "texture.h"
#include "texture_img.h"
#include "rectangle.h"
#include "bvh.h"
#include "aabb.h"
//#include "box.h"
#define stds std::

using namespace std;


vec3 ray_color(const ray& r, const vec3& background, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return vec3(0, 0, 0);

    if (!world.hit(r, 0.001, infinity, rec)) return background; 
    ray scattered;
    vec3 attenuation;
    vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered)) return emitted;
            
    return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
 
}

//    vec3 unit_direction = unit_vector(r.d());
//    auto t = 0.5 * (unit_direction.y() + 1.0);
//    return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
//}

hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambertian>(vec3(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(vec3(0, -1000, 0), 1000, ground_material));

    auto material2 = make_shared<lambertian>(vec3(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(vec3(-4, 1, 0), 1.0, material2));


    return world;
}


hittable_list two_spheres() {
    hittable_list objects;

    auto checker = make_shared<checker_texture>(vec3(0.2, 0.3, 0.1), vec3(0.9, 0.9, 0.9));

    objects.add(make_shared<sphere>(vec3(0, -10, 0), 10, make_shared<lambertian>(checker)));
    objects.add(make_shared<sphere>(vec3(0, 10, 0), 10, make_shared<lambertian>(checker)));

    return objects;
}

hittable_list earth() {
    int nx, ny, nn;
    unsigned char* texture_data = stbi_load("earthmap.jpg", &nx, &ny, &nn, 0);

    auto earth_surface =
        make_shared<lambertian>(make_shared<image_texture>(texture_data, nx, ny));
    auto globe = make_shared<sphere>(vec3(0, 0, 0), 2, earth_surface);

    return hittable_list(globe);
}

hittable_list simple_light() {
    hittable_list objects;

    auto pertext = make_shared<constant_texture>(vec3(0.73, 0.73, 0.73));
    objects.add(make_shared<sphere>(vec3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(vec3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    auto difflight = make_shared<diffuse_light>(vec3(4, 4, 4));
    objects.add(make_shared<sphere>(vec3(0, 7, 0), 2, difflight));
    //objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));

    return objects;
}

//hittable_list cornell_box() {
//    hittable_list objects;
//
//    auto red = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.65, 0.05, 0.05)));
//    auto white = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.73, 0.73, 0.73)));
//    auto green = make_shared<lambertian>(make_shared<constant_texture>(vec3(0.12, 0.45, 0.15)));
//    auto light = make_shared<diffuse_light>(make_shared<constant_texture>(vec3(15, 15, 15)));
//
//    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
//    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
//    objects.add(make_shared<xz_rect>(213, 343, 227, 332, 554, light));
//    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
//    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
//    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
//    //objects.add(make_shared<box>(vec3(130, 0, 65), vec3(295, 165, 230), white));
//    //objects.add(make_shared<box>(vec3(265, 0, 295), vec3(430, 330, 460), white));
//    return objects;
//}


int main() {
    //// Image
    ////const char* filename = "n.ppm";
    const vec3 background(0,0,0);
    const auto aspect_ratio = 2.0;
    const int image_width = 400;
    const int image_height = 200;
    const int samples_per_pixel = 100;
    const int max_depth = 250;
    //// Render
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    auto world = simple_light();

    vec3 lookfrom(26,3,-10);
    vec3 lookat(0,2,0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.0;
    auto vfov = 20.0;
    camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus,0.0,1.0);
    for (int j = image_height - 1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining:" << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            /*vec3 color(double(i) / image_width, double(j) / image_height, 0.2);
            color.write_color(std::cout);*/
            vec3 color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random()) / image_width; auto v = (j + random()) / image_height;
                ray r = cam.get_ray(u, v);
                color += ray_color(r,background,world, max_depth);
            }
            color.write_color(std::cout, samples_per_pixel);
        }
    }
    ////std::cout << ir << ' ' << ig << ' ' << ib << '\n';
    //return 0;
    std::cerr << "\nDone.\n";
    //return basel(argc, argv);

}

