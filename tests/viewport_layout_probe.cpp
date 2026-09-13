#include "cov/viewer_layout.hpp"
#include <iomanip>
#include <iostream>

int main() {
    std::cout << std::setprecision(17);
    int id = 0, fw = 0, fh = 0;
    float w = 0, h = 0, scale = 0, fov = 0;
    while (std::cin >> id >> w >> h >> fw >> fh >> scale >> fov) {
        const auto layout = cov::viewer_layout(w, h, fw, fh, scale);
        const auto& s = layout.scene;
        const auto& c = layout.controls;
        const auto& v = layout.framebuffer;
        const auto p = cov::scene_projection(v.width, v.height, fov);
        std::cout << "{\"id\":" << id << ",\"controls\":[" << c.x << ',' << c.y << ',' << c.width << ',' << c.height
            << "],\"scene\":[" << s.x << ',' << s.y << ',' << s.width << ',' << s.height
            << "],\"viewport\":[" << v.x << ',' << v.y << ',' << v.width << ',' << v.height
            << "],\"projection\":[" << p.aspect << ',' << p.tan_half_vertical_fov
            << "],\"hit\":[" << s.contains(s.x+s.width*.5f,s.y+s.height*.5f) << ','
            << s.contains(s.x+s.width,s.y+s.height*.5f) << ','
            << s.contains(c.x+c.width*.5f,c.y+c.height*.5f) << "]}\n";
    }
    return std::cin.eof() ? 0 : 2;
}
