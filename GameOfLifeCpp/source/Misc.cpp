#include "Misc.h"
#include "Vector.h"


    constexpr vec2 misc::vec2lerp(vec2 a, vec2 b, float f) noexcept
    {
        return vec2(misc::lerp(a.x, b.x, f), misc::lerp(a.y, b.y, f));
    }


    float misc::modf(const float x, const float y) noexcept {
        return x - y * floor(x / y);
    }

    constexpr int32_t misc::mod(const int32_t x, const int32_t y) noexcept {
        int32_t mod = x % y;
        // if the signs are different and modulo not zero, adjust result
        if ((x ^ y) < 0 && mod != 0) {
            mod += y;
        }
        return mod;
    }
