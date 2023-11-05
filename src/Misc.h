#pragma once
#include <type_traits>
#include <stdint.h>
#include "Vector.h"
#include <cassert>


#define print_msg(msg) { std::cout << (msg) << std::endl; }

#include<iostream>
#include<mutex>
static std::mutex m{};

inline void print(const char* msg, uint32_t arg) {
    std::lock_guard<std::mutex> g{ m };
    std::cout << msg << arg << std::endl;
}

namespace misc {
    inline constexpr vec2i indexToCoord(int32_t const index, int32_t const width) {
        auto const div = index / width;
        auto const mod = index % width;
        if(mod < 0) return { div-1, width + mod };
        else return { div, mod };
    }

    template<class T>
    inline constexpr T lerp(T a, T b, T f) noexcept {
        return a + f * (b - a);
    }

    template<class T>
    inline constexpr T unlerp(T a, T b, T f) noexcept {
        return (f - a) / (b - a);
    }

    template<class C>
    inline constexpr vec2<C> vec2lerp(const vec2<C> a, const vec2<C> b, const C f) noexcept {
        return vec2<C>(misc::lerp(a.x, b.x, f), misc::lerp(a.y, b.y, f));
    }

    template <typename E>
    constexpr inline typename std::underlying_type<E>::type to_underlying(const E e) noexcept {
        return static_cast<typename std::underlying_type<E>::type>(e);
    }

    inline float modf(const float x, const float y) noexcept {
        return x - y * floor(x / y);
    }

    inline constexpr int32_t mod(const int32_t x, const int32_t y) noexcept {
        int32_t mod = x % y;
        // if the signs are different and modulo not zero, adjust result
        if ((x ^ y) < 0 && mod != 0) {
            mod += y;
        }
        return mod;
    }

    inline constexpr uint32_t umod(const uint32_t x, const uint32_t y) noexcept {
        return x % y;
    }

    template<class Type>
    inline constexpr const Type& max(const Type& t1, const Type& t2) {
        if (t1 >= t2) {
            return t1;
        }
        return t2;
    }

    template<class Type>
    inline constexpr const Type& min(const Type& t1, const Type& t2) {
        if (t1 <= t2) {
            return t1;
        }
        return t2;
    }

    inline constexpr uint32_t roundUpIntTo(uint32_t number, uint32_t round) {
        const auto remainder = (number - 1) % round;
        const auto result = number + (round - remainder - 1);
        assert(result % round == 0);
        return result;
    }

    inline constexpr int32_t roundDownIntTo(int32_t number, int32_t round) {
        const auto remainder = misc::mod(number, round);
        const auto result = number - remainder;
        assert((result % round == 0) || (result <= number));
        return result;
    }

    inline constexpr uint32_t intDivCeil(uint32_t number, uint32_t round) {
        return roundUpIntTo(number, round) / round;
    }

    inline constexpr int32_t intDivFloor(int32_t number, int32_t round) {
        return roundDownIntTo(number, round) / round;
    }

    template<class T>
    inline constexpr T clamp(const T value, const T b1, const T b2) {
        if (b1 > b2) return min(b1, max(value, b2));
        else return min(b2, max(value, b1));
    }

    template<class T>
    inline constexpr T map(const T value, const T va, const T vb, const T a, const T b) {
        return lerp(a, b, unlerp(va, vb, value));
    }
}
