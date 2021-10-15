#pragma once

constexpr float lerp(float a, float b, float f) noexcept
{
    return a + f * (b - a);
}

template <typename E>
constexpr typename std::underlying_type<E>::type to_underlying(E e) noexcept {
    return static_cast<typename std::underlying_type<E>::type>(e);
}

 float modf(const float x, const float y) noexcept {
    return x - y * floor(x / y);
}

constexpr int32_t mod(const int32_t x, const int32_t y) noexcept {
    int32_t mod = x % y;
    // if the signs are different and modulo not zero, adjust result
    if ((x^y) < 0 && mod != 0) {
        mod += y;
    }
    return mod;
}