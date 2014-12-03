#include "debug.h"

void
generate_debug_colors(std::vector<math::Vec4f> & colors) {
    for (float s = 1.0f; s > 0.0f; s -= 0.4) {
        for (float v = 1.0f; v > 0.0f; v -= 0.3) {
            for (float h = 0.0f; h < 360.0f; h += 30.0f) {
                float c = v * s;
                float x = c * (1.0f - fabs(fmod(h / 60.0f, 2.0f) - 1.0f));
                float m = v - c;

                math::Vec4f color;
                if (0 <= h && h < 60)
                    color = math::Vec4f(c, x, 0.0f, 1.0f);
                if (60 <= h && h < 120)
                    color = math::Vec4f(x, c, 0.0f, 1.0f);
                if (120 <= h && h < 180)
                    color = math::Vec4f(0.0f, c, x, 1.0f);
                if (180 <= h && h < 240)
                    color = math::Vec4f(0.0f, x, c, 1.0f);
                if (240 <= h && h < 300)
                    color = math::Vec4f(x, 0.0f, c, 1.0f);
                if (300 <= h && h < 360)
                    color = math::Vec4f(c, 0.0f, x, 1.0f);

                color = color + math::Vec4f(m, m, m, 0.0f);
                colors.push_back(color);
            }
        }
    }
}
