#pragma once

#include<stdint.h>
#include<cmath>
#include<iostream>

template<typename C>
struct vec2 {
	C x, y;

	constexpr vec2() = default;
	constexpr vec2(const C value) : x(value), y(value) {};
	constexpr vec2(const C x_, const C y_) : x(x_), y(y_) {};

	constexpr vec2<C> operator+(const vec2<C>& o) const {
		return vec2(x + o.x, y + o.y);
	}
	constexpr vec2<C> operator-(const vec2<C>& o) const {
		return vec2(x - o.x, y - o.y);
	}
	constexpr vec2<C> operator/(const vec2<C>& o) const {
		return vec2(x / o.x, y / o.y);
	}
	constexpr vec2<C> operator/(const C o) const {
		return vec2(x / o, y / o);
	}
	constexpr vec2<C> operator*(const vec2<C>& o) const {
		return vec2(x * o.x, y * o.y);
	}
	constexpr vec2<C> operator*(const C o) const {
		return vec2(x * o, y * o);
	}
	vec2<C>& operator+=(const vec2<C>& o) {
		x += o.x;
		y += o.y;
		return *this;
	}

	vec2<C>& operator-=(const vec2<C>& o) {
		x -= o.x;
		y -= o.y;
		return *this;
	}

	constexpr C dot(vec2<C> o) const {
		return x * o.x + y * o.y;
	}

	constexpr C lengthSuqare() const {
		return this->dot(*this);
	}

	constexpr C length() const {
		return std::sqrt(this->lengthSuqare());
	}

    friend std::ostream &operator<<(std::ostream &o, vec2<C> const &it) {
        return o << '(' << it.x << ", " << it.y << ')';
    }
};

using vec2i = vec2<int>;
using vec2d = vec2<double>;
