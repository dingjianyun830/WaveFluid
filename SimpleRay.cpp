#include "vec3.h"
#include "vec4.h"
#include "camera.h"
#include "ray.h"
#include "perlin.h" 
#include "Wave.h"
#include <opencv2\opencv.hpp>
#include <vector>
#include <string>
#include <time.h>   
//time_t timer;
//struct tm curr_time;

using namespace std;
// define some parameters
// set the size of picture
int ImgWidth = 256;
int ImgHeight = 256;
vector<vec3> refractedRay;
vector<ray> OldRay;
vector<ray> NewRay;

vector<vec3> getNormal(double normx[GRIDW][GRIDH], double normy[GRIDW][GRIDH], double normz[GRIDW][GRIDH])
{
	vector<vec3> n;
	n.assign(GRIDH*GRIDW, vec3(0, 0, 0));
	for (int i = 0; i < GRIDH; i++)
	{
		for (int j = 0; j < GRIDW; j++)
		{
			n[i*GRIDW + j] = vec3(-normx[i][j], -normy[i][j], -normz[i][j]);
		}
	}
	return n;
}

vector<double> getHeight(double avgHeight[GRIDW][GRIDH])
{
	vector<double> h;
	h.assign(GRIDH*GRIDW, 0.0);
	for (int i = 0; i < GRIDH; i++)
	{
		for (int j = 0; j < GRIDW; j++)
		{
			h[i*GRIDW + j] = -avgHeight[i][j];
		}
	}
	return h;
}

vec4 MergeColor(const vec4 & col1, const vec4 & col2, double f)
{
	vec4 Mcolor;
	vec4 base = col1 / 255;
	vec4 blend = col2 / 255;
	Mcolor[3] = base[3] + (1 - base[3])*blend[3]*f;
	Mcolor[0] = base[0]*(1 - base[3]) + (blend[0]*blend[3] *f + base[0]*(1 - blend[3] *f))*base[3];
	Mcolor[1] = base[1]*(1 - base[3]) + (blend[1]*blend[3] *f + base[1]*(1 - blend[3] *f))*base[3];
	Mcolor[2] = base[2]*(1 - base[3]) + (blend[2]*blend[3] *f + base[2]*(1 - blend[3] *f))*base[3];
	return Mcolor;
}

bool refract(const vec3& v, const vec3& n, double ni_over_nt, vec3& refracted)
{
	vec3 uv = unit_vector(v);
	double dt = dot(uv, n);
	double discriminant = 1.0 - ni_over_nt*(1 - dt*dt);
	if (discriminant > 0)
	{
		refracted = ni_over_nt*(uv - dt*n) - sqrt(discriminant)*n;
		return true;
	}
	else
	{
		return false;
	}
}

void Dyecolor(cv::Mat tImage, cv::Mat fImage, int time, vector<double> gridHeight)
{
	// create the Mat to save the result
	double f = 0.0;
	cv::Mat rImage(ImgWidth, ImgHeight, CV_8UC4);
	for (int i = 0; i < ImgHeight; i++)
	{
		for (int j = 0; j < ImgWidth; j++)
		{
			ray r = NewRay[i*ImgHeight + j];
			double gh = gridHeight[i*ImgHeight + j];
			//double gh = 0.5;
			vec3 p = r.point_at_parameter(gh);
			cv::Vec4b bgraN = rImage.at<cv::Vec4b>(i, j);
			int u = int(127.99 * (p[0] + 1.0));
			int v = int(127.99 * (p[1] + 1.0));
			int uu = int(127.99 * (r.origin()[0] + 1.0));
			int vv = int(127.99 * (r.origin()[1] + 1.0));
			if (u < 0 || v < 0 || uu < 0 || vv < 0)
			{
				rImage.at<cv::Vec4b>(i, j)[0] = 0;
				rImage.at<cv::Vec4b>(i, j)[1] = 0;
				rImage.at<cv::Vec4b>(i, j)[2] = 0;
				rImage.at<cv::Vec4b>(i, j)[3] = 0;
			}
			else if (u > 255 || v > 255 || uu > 255 || vv > 255)
			{
				rImage.at<cv::Vec4b>(i, j)[0] = 0;
				rImage.at<cv::Vec4b>(i, j)[1] = 0;
				rImage.at<cv::Vec4b>(i, j)[2] = 0;
				rImage.at<cv::Vec4b>(i, j)[3] = 0;
			}
			else
			{
				vec4 Tcolor = vec4(tImage.at<cv::Vec4b>(u, v)[0], 
									tImage.at<cv::Vec4b>(u, v)[1],
									tImage.at<cv::Vec4b>(u, v)[2],
									tImage.at<cv::Vec4b>(u, v)[3]);
				vec4 Fcolor = vec4(fImage.at<cv::Vec4b>(uu, vv)[0],
									fImage.at<cv::Vec4b>(uu, vv)[1],
									fImage.at<cv::Vec4b>(uu, vv)[2],
									fImage.at<cv::Vec4b>(uu, vv)[3]);
				vec4 Mcolor = MergeColor(Tcolor, Fcolor, f);
				rImage.at<cv::Vec4b>(i, j)[0] = int(255.99*Mcolor[0]);
				rImage.at<cv::Vec4b>(i, j)[1] = int(255.99*Mcolor[1]);
				rImage.at<cv::Vec4b>(i, j)[2] = int(255.99*Mcolor[2]);
				rImage.at<cv::Vec4b>(i, j)[3] = int(255.99*Mcolor[3]);
			}
		}
	}
	string fileName = "testImage" + to_string(time) + ".png";
	cv::imwrite(fileName, rImage);
}

void ComputeRayTracing(vector<vec3> gridNormal, vector<double> gridHeight)
{
	// set the 3D location of camera and plate
	vec3 lower_left_corner(-1.0, -1.0, -2.0);
	vec3 horizontal(2.0, 0.0, 0.0);
	vec3 vertical(0.0, 2.0, 0.0);
	vec3 origin(0.0, 0.0, 0.0);

	refractedRay.assign(ImgWidth*ImgWidth,vec3(0,0,0));
	bool m_bRefract = true;
	for (int j = 0; j < ImgHeight; j++)
	{
		for (int i = 0; i < ImgWidth; i++)
		{
			double u = double(i) / double(ImgWidth);
			double v = double(j) / double(ImgHeight);

			ray r(origin, lower_left_corner + u*horizontal + v*vertical, 0);
			OldRay.push_back(r);
			double ni_over_nt = 1 / 1.33;
			vec3 normal = gridNormal[j*ImgHeight + i];
			double gh = gridHeight[j*ImgHeight + i];
			//vec3 normal = vec3(0,0,-0.5);
			//double gh = 0.5;
			if (m_bRefract)
			{
				refract(r.point_at_parameter(2 - gh), normal, ni_over_nt, refractedRay[j*ImgHeight + i]);
				ray r_out(r.point_at_parameter(2 - gh), refractedRay[j*ImgHeight + i], 0.0);
				NewRay.push_back(r_out);
			}
			else
			{
				NewRay.push_back(r);
			}
		}
	}
}

void WriteResult(vector<vec3> gridNormal, vector<double> gridHeight, int time)
{
	string Output = "result" + to_string(time) + ".txt";
	FILE *fp1 = fopen(Output.c_str(), "w+");
	for (int i = 0; i < gridNormal.size(); i++)
	{
		vec3 n = gridNormal[i];
		double h = gridHeight[i];
		fprintf(fp1, "%f %f %f %f\n", n[0], n[1], n[2], h);
	}
	fclose(fp1);
}

double myTime()
{
	double start = static_cast<double>(cv::getTickCount());
	double time = ((double)cv::getTickCount() - start) / cv::getTickFrequency();
	return time;
}

int main()
{
	// open the texture image
	cv::Mat img1 = cv::imread("bottom.jpg");
	cv::Mat img2 = cv::imread("floor.jpg");
	cv::Mat texture1(ImgWidth, ImgHeight, CV_8UC4);
	cv::Mat texture2(ImgWidth, ImgHeight, CV_8UC4);
	cv::cvtColor(img1, texture1, CV_BGR2BGRA);
	cv::cvtColor(img2, texture2, CV_BGR2BGRA);
	cv::imwrite("t1.png", texture1);
	cv::imwrite("t2.png", texture2);
	//cout<<texture1.channels();
	if (img1.empty() || img2.empty()) // Check for invalid input
	{
		std::cout << "Could not open or find the image" << std::endl;
		return -1;
	}

	// set a water floor
	double t, dt_total, t_old;

	// Initialize simulation
	init_vertices();
	init_grid();
	adjust_grid();

	// Initialize timer

	t_old = myTime() - 10;
	int count = 0;
	while (count<200)
	{
		t = myTime();
		dt_total = t - t_old;
		t_old = t;

		// Safety - iterate if dt_total is too large
		while (dt_total > 0.f)
		{
			// Select iteration time step
			dt = dt_total > MAX_DELTA_T ? MAX_DELTA_T : dt_total;
			dt_total -= dt;

			// Calculate wave propagation
			calc_grid();
		}

		// Compute height of each vertex
		adjust_grid();

		vector<vec3> gridNormal = getNormal(normx, normy, normz);
		vector<double> gridHeight = getHeight(avgHeight);
		if (count % 10 == 0)
		{
			ComputeRayTracing(gridNormal, gridHeight);
			Dyecolor(texture1, texture2, count, gridHeight);
			WriteResult(gridNormal, gridHeight, count);
		}
		count++;
	}

	return 1;
}