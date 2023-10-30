#include "Misc.h"
#include "Vector.h"





    float misc::modf(const float x, const float y) noexcept {
        return x - static_cast<double>(y) * floor(x / static_cast<double>(y));
    }
