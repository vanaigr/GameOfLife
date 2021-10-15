#pragma once
#include <type_traits>
#include <stdint.h>
#include "Vector.h"
#include <cassert>

namespace misc {
    template<class T>
    constexpr T lerp(T a, T b, T f) noexcept {
        return a + f * (b - a);
    }

    template<class T>
    inline constexpr T unlerp(T a, T b, T f) noexcept {
        return (f - a) / (b - a);
    }

    extern constexpr vec2 vec2lerp(vec2 a, vec2 b, float f) noexcept;

    template <typename E>
    constexpr inline typename std::underlying_type<E>::type to_underlying(const E e) noexcept {
        return static_cast<typename std::underlying_type<E>::type>(e);
    }

    extern float modf(const float x, const float y) noexcept;

    extern constexpr int32_t mod(const int32_t x, const int32_t y) noexcept;

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