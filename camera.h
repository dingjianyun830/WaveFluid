#ifndef CAMERAH
#define CAMERAH

#include "ray.h"
#define M_PI 3.1415926

vec3 random_in_unit_disk()
{
	vec3 p;
	do
	{
		p = 2.0*vec3(randd(), randd(), 0) - vec3(1, 1, 0);
	} while (dot(p, p) >= 1.0);
	return p;
}

class camera
{
public:
	vec3 origin;
	vec3 lower_left_corner;
	vec3 horizontal;
	vec3 vertical;
	vec3 u, v, w;
	double lens_radius;
	double time0, time1;

	camera()
	{
		lower_left_corner = vec3(-2.0, -1.0, -1.0);
		horizontal = vec3(4.0, 0.0, 0.0);
		vertical = vec3(0.0, 2.0, 0.0);
		origin = vec3(0.0, 0.0, 0.0);
	}

	camera(double vfov, double aspect)
	{// vfov is top to bottom in degrees
		double theta = vfov * M_PI/180;
		double half_height = tan(theta / 2);
		double half_width = aspect*half_height;
		lower_left_corner = vec3(-half_width, -half_height, -1.0);
		horizontal = vec3(2 * half_width, 0.0, 0.0);
		vertical = vec3(0.0, 2 * half_height, 0.0);
		origin = vec3(0.0, 0.0, 0.0);
	}

	camera(vec3 lookfrom, vec3 lookat, vec3 vup, double vfov, double aspect)
	{// vfov is top to bottom in degrees
		double theta = vfov*M_PI / 180;
		double half_height = tan(theta / 2);
		double half_width = aspect * half_height;
		origin = lookfrom;
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);
		lower_left_corner = vec3(-half_width, -half_height, -1.0);
		lower_left_corner = origin - half_width*u - half_height*v - w;
		horizontal = 2 * half_width*u;
		vertical = 2 * half_height*v;
	}

	camera(vec3 lookfrom, vec3 lookat, vec3 vup, double vfov, double aspect, double aperture, double focus_dist, double t0, double t1)
	{// vfov is top to bottom in degrees
		time0 = t0;
		time1 = t1;
		lens_radius = aperture / 2;
		double theta = vfov*M_PI / 180;
		double half_height = tan(theta / 2);
		double half_width = aspect * half_height;
		origin = lookfrom;
		w = unit_vector(lookfrom - lookat);
		u = unit_vector(cross(vup, w));
		v = cross(w, u);
		lower_left_corner = origin - half_width*focus_dist*u - half_height*focus_dist*v - focus_dist*w;
		horizontal = 2 * half_width*focus_dist*u;
		vertical = 2 * half_height*focus_dist*v;
	}

	ray get_ray(double s, double t)
	{
		vec3 rd = lens_radius*random_in_unit_disk();
		vec3 offset = rd.x()*u + rd.y()*v;
		double time = time0 + randd()*(time1 - time0);
		return ray(origin + offset, lower_left_corner + s*horizontal + t*vertical - origin - offset,time);
	}
};
#endif // !CAMERAH
