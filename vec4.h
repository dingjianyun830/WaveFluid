#pragma once
#include<iostream>
#include<stdlib.h>
#include<math.h>

class vec4
{
public:
	double e[4];
	vec4() {}
	vec4(double e0, double e1, double e2, double e3)
	{
		e[0] = e0;
		e[1] = e1;
		e[2] = e2;
		e[3] = e3;
	}

	inline double x() const { return e[0]; }
	inline double y() const { return e[1]; }
	inline double z() const { return e[2]; }
	inline double w() const { return e[3]; }
	inline double r() const { return e[0]; }
	inline double g() const { return e[1]; }
	inline double b() const { return e[2]; }
	inline double a() const { return e[3]; }

	inline const vec4& operator+() const { return *this; }
	inline vec4 operator-() const { return vec4(-e[0], -e[1], -e[2], -e[3]); }
	inline double operator[] (int i) const { return e[i]; }
	inline double& operator[] (int i) { return e[i]; }

	inline vec4& operator+=(const vec4 &v2);
	inline vec4& operator-=(const vec4 &v2);
	inline vec4& operator*=(const vec4 &v2);
	inline vec4& operator/=(const vec4 &v2);
	inline vec4& operator*=(const double t);
	inline vec4& operator/=(const double t);

	inline double length() const
	{
		return sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2] + e[3]*e[3]);
	}

	inline double squared_length() const
	{
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2] + e[3] * e[3];
	}

	inline void make_unit_vector();
};

inline std::istream& operator >> (std::istream &is, vec4 &t)
{
	is >> t.e[0] >> t.e[1] >> t.e[2] >> t.e[3];
	return is;
}

inline std::ostream& operator >> (std::ostream &os, const vec4 &t)
{
	os << t.e[0] << t.e[1] << t.e[2] << t.e[3];
	return os;
}

inline void vec4::make_unit_vector()
{
	double k = 1.0 / sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
	e[0] *= k;
	e[1] *= k;
	e[2] *= k;
	e[3] *= k;
}

inline vec4 operator+(const vec4 &v1, const vec4 &v2)
{
	return vec4(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2], v1.e[3] + v2.e[3]);
}

inline vec4 operator-(const vec4 &v1, const vec4 &v2)
{
	return vec4(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2], v1.e[3] - v2.e[3]);
}

inline vec4 operator*(const vec4 &v1, const vec4 &v2)
{
	return vec4(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2], v1.e[3] * v2.e[3]);
}

inline vec4 operator/(const vec4 &v1, const vec4 &v2)
{
	return vec4(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2], v1.e[3] * v2.e[3]);
}

inline vec4 operator*(double t, const vec4 &v)
{
	return vec4(t*v.e[0], t*v.e[1], t*v.e[2], t*v.e[3]);
}

inline vec4 operator/(vec4 v, double t)
{
	return vec4(v.e[0] / t, v.e[1] / t, v.e[2] / t, v.e[3] / t);
}

inline double dot(const vec4 &v1, const vec4 &v2)
{
	return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2] + v1.e[3] * v2.e[3];
}


inline vec4& vec4::operator+=(const vec4 &v)
{
	e[0] += v.e[0];
	e[1] += v.e[1];
	e[2] += v.e[2];
	e[3] += v.e[3];
	return *this;
}

inline vec4& vec4::operator-=(const vec4 &v)
{
	e[0] -= v.e[0];
	e[1] -= v.e[1];
	e[2] -= v.e[2];
	e[3] -= v.e[3];
	return *this;
}

inline vec4& vec4::operator*=(const vec4 &v)
{
	e[0] *= v.e[0];
	e[1] *= v.e[1];
	e[2] *= v.e[2];
	e[3] *= v.e[3];
	return *this;
}

inline vec4& vec4::operator/=(const vec4 &v)
{
	e[0] /= v.e[0];
	e[1] /= v.e[1];
	e[2] /= v.e[2];
	e[3] /= v.e[3];
	return *this;
}

inline vec4& vec4::operator*=(const double t)
{
	e[0] *= t;
	e[1] *= t;
	e[2] *= t;
	e[3] *= t;
	return *this;
}

inline vec4& vec4::operator/=(const double t)
{
	e[0] /= t;
	e[1] /= t;
	e[2] /= t;
	e[3] /= t;
	return *this;
}

inline vec4 unit_vector(vec4 v)
{
	return v / v.length();
}