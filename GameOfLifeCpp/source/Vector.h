#pragma once

#include<stdint.h>

template<typename C>
struct vec2
{
	C x, y;

public:
	constexpr vec2() = default;
	constexpr vec2(const C value) : x(value), y(value) {};
	constexpr vec2(const C x_, const C y_) : x(x_), y(y_) {};

	inline constexpr vec2<C> operator+(const vec2<C>& o) const {
		return vec2(x + o.x, y + o.y);
	}
	inline constexpr vec2<C> operator-(const vec2<C>& o) const {
		return vec2(x - o.x, y - o.y);
	}
	inline constexpr vec2<C> operator/(const vec2<C>& o) const {
		return vec2(x / o.x, y / o.y);
	}
	inline constexpr vec2<C> operator/(const C o) const {
		return vec2(x / o, y / o);
	}
	inline constexpr vec2<C> operator*(const vec2<C>& o) const {
		return vec2(x * o.x, y * o.y);
	}
	inline constexpr vec2<C> operator*(const C o) const {
		return vec2(x * o, y * o);
	}
	inline vec2<C>& operator+=(const vec2<C>& o) {
		x += o.x;
		y += o.y;
		return *this;
	}

	inline vec2<C>& operator-=(const vec2<C>& o) {
		x -= o.x;
		y -= o.y;
		return *this;
	}

	inline constexpr C dot(vec2<C> o) const {
		return x * o.x + y * o.y;
	}

	inline constexpr C lengthSuqare() const {
		return this->dot(*this);
	}

	inline constexpr C length() const {
		return sqrt(this->lengthSuqare());
	}
};

struct vec2i
{
	int32_t x, y;

public:
	inline constexpr vec2i(const int32_t x_, const int32_t y_) : x(x_), y(y_) {};

	inline constexpr vec2i operator+(const vec2i o) const {
		return vec2i(x + o.x, y + o.y);
	}
	inline constexpr vec2i operator-(const vec2i o) const {
		return vec2i(x - o.x, y - o.y);
	}
};

struct vec2ui
{
	uint32_t x, y;

public:
	inline constexpr vec2ui(const uint32_t x_, const uint32_t y_) : x(x_), y(y_) {};

	inline constexpr vec2ui operator+(const vec2ui o) const {
		return vec2ui(x + o.x, y + o.y);
	}
	inline constexpr vec2ui operator-(const vec2ui o) const {
		return vec2ui(x - o.x, y - o.y);
	}
};