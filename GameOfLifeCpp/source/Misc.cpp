#include "Misc.h"
#include "Vector.h"


    constexpr vec2 misc::vec2lerp(vec2 a, vec2 b, float f) noexcept
    {
        return vec2(misc::lerp(a.x, b.x, f), misc::lerp(a.y, b.y, f));
    }


    float misc::modf(const float x, const float y) noexcept {
        return x - static_cast<double>(y) * floor(x / static_cast<double>(y));
    }
