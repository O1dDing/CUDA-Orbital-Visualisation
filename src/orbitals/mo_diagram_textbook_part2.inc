}

struct Canvas {
    int width = 1;
    int height = 1;
    std::vector<std::uint8_t> rgba;

    Canvas(const int w, const int h)
        : width(std::max(1, w)), height(std::max(1, h)),
          rgba(static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u, 255u) {}

    void pixel(const int x, const int y, const std::uint8_t r,
               const std::uint8_t g, const std::uint8_t b) {
        if (x < 0 || y < 0 || x >= width || y >= height) return;
        const std::size_t i = (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                               static_cast<std::size_t>(x)) * 4u;
        rgba[i] = r; rgba[i + 1] = g; rgba[i + 2] = b; rgba[i + 3] = 255;
    }

    void line(int x0, int y0, const int x1, const int y1,
              const std::uint8_t r, const std::uint8_t g,
              const std::uint8_t b, const int thickness = 1) {
        const int dx = std::abs(x1 - x0);
        const int sx = x0 < x1 ? 1 : -1;
        const int dy = -std::abs(y1 - y0);
        const int sy = y0 < y1 ? 1 : -1;
        int error = dx + dy;
        for (;;) {
            const int radius = std::max(0, thickness / 2);
            for (int oy = -radius; oy <= radius; ++oy) {
                for (int ox = -radius; ox <= radius; ++ox) pixel(x0 + ox, y0 + oy, r, g, b);
            }
            if (x0 == x1 && y0 == y1) break;
            const int e2 = 2 * error;
            if (e2 >= dy) { error += dy; x0 += sx; }
            if (e2 <= dx) { error += dx; y0 += sy; }
        }
    }
};

std::array<std::uint8_t, 7> glyph_rows(char c) {
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    switch (c) {
        case '0': return {14,17,19,21,25,17,14}; case '1': return {4,12,4,4,4,4,14};
        case '2': return {14,17,1,2,4,8,31}; case '3': return {30,1,1,14,1,1,30};
        case '4': return {2,6,10,18,31,2,2}; case '5': return {31,16,16,30,1,1,30};
        case '6': return {14,16,16,30,17,17,14}; case '7': return {31,1,2,4,8,8,8};
        case '8': return {14,17,17,14,17,17,14}; case '9': return {14,17,17,15,1,1,14};
        case 'A': return {14,17,17,31,17,17,17}; case 'B': return {30,17,17,30,17,17,30};
        case 'C': return {14,17,16,16,16,17,14}; case 'D': return {30,17,17,17,17,17,30};
        case 'E': return {31,16,16,30,16,16,31}; case 'F': return {31,16,16,30,16,16,16};
        case 'G': return {14,17,16,23,17,17,15}; case 'H': return {17,17,17,31,17,17,17};
        case 'I': return {14,4,4,4,4,4,14}; case 'J': return {7,2,2,2,2,18,12};
        case 'K': return {17,18,20,24,20,18,17}; case 'L': return {16,16,16,16,16,16,31};
        case 'M': return {17,27,21,21,17,17,17}; case 'N': return {17,25,21,19,17,17,17};
        case 'O': return {14,17,17,17,17,17,14}; case 'P': return {30,17,17,30,16,16,16};
        case 'Q': return {14,17,17,17,21,18,13}; case 'R': return {30,17,17,30,20,18,17};
        case 'S': return {15,16,16,14,1,1,30}; case 'T': return {31,4,4,4,4,4,4};
        case 'U': return {17,17,17,17,17,17,14}; case 'V': return {17,17,17,17,17,10,4};
        case 'W': return {17,17,17,21,21,21,10}; case 'X': return {17,17,10,4,10,17,17};
        case 'Y': return {17,17,10,4,4,4,4}; case 'Z': return {31,1,2,4,8,16,31};
        case '-': return {0,0,0,31,0,0,0}; case '.': return {0,0,0,0,0,12,12};
        case '/': return {1,2,2,4,8,8,16}; case '+': return {0,4,4,31,4,4,0};
        case ':': return {0,12,12,0,12,12,0}; case '*': return {0,21,14,31,14,21,0};
        case '?': return {14,17,1,2,4,0,4}; default: return {0,0,0,0,0,0,0};
    }
}

void draw_text(Canvas& canvas, int x, const int y, const std::string& text,
               const std::uint8_t r = 30, const std::uint8_t g = 35,
               const std::uint8_t b = 43, const int scale = 2) {
    for (const char c : text) {
        if (c == ' ') { x += 4 * scale; continue; }
        const auto rows = glyph_rows(c);
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if ((rows[static_cast<std::size_t>(row)] & (1u << (4 - col))) == 0) continue;
                for (int yy = 0; yy < scale; ++yy) {
                    for (int xx = 0; xx < scale; ++xx) {
                        canvas.pixel(x + col * scale + xx, y + row * scale + yy, r, g, b);
                    }
                }
            }
        }
        x += 6 * scale;
    }
}

void draw_arrow(Canvas& canvas, const int x, const int baseline, const bool up) {
    const int tip = up ? baseline - 24 : baseline + 24;
    canvas.line(x, baseline, x, tip, 20, 20, 20, 2);
    const int dy = up ? 6 : -6;
    canvas.line(x, tip, x - 5, tip + dy, 20, 20, 20, 2);
    canvas.line(x, tip, x + 5, tip + dy, 20, 20, 20, 2);
}

std::uint32_t crc32(const std::uint8_t* data, const std::size_t size) {
    std::uint32_t crc = 0xffffffffu;
    for (std::size_t i = 0; i < size; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit) {
            const std::uint32_t mask = 0u - (crc & 1u);
            crc = (crc >> 1u) ^ (0xedb88320u & mask);
        }
    }
    return ~crc;
}

std::uint32_t adler32(const std::uint8_t* data, const std::size_t size) {
    std::uint32_t a = 1u, b = 0u;
    for (std::size_t i = 0; i < size; ++i) {
        a = (a + data[i]) % 65521u;
        b = (b + a) % 65521u;
    }
    return (b << 16u) | a;
}

void push_u32(std::vector<std::uint8_t>& out, const std::uint32_t value) {
    out.push_back(static_cast<std::uint8_t>(value >> 24u));
    out.push_back(static_cast<std::uint8_t>(value >> 16u));
    out.push_back(static_cast<std::uint8_t>(value >> 8u));
    out.push_back(static_cast<std::uint8_t>(value));
}

void append_chunk(std::vector<std::uint8_t>& out, const char type[4],
                  const std::vector<std::uint8_t>& payload) {
    push_u32(out, static_cast<std::uint32_t>(payload.size()));
    const std::size_t begin = out.size();
    out.insert(out.end(), type, type + 4);
    out.insert(out.end(), payload.begin(), payload.end());
    push_u32(out, crc32(out.data() + begin, 4 + payload.size()));
}

bool save_png(const Canvas& canvas, const std::filesystem::path& path, std::string* error) {
    std::vector<std::uint8_t> raw;
    for (int y = 0; y < canvas.height; ++y) {
        raw.push_back(0);
        const std::size_t begin = static_cast<std::size_t>(y) *
                                  static_cast<std::size_t>(canvas.width) * 4u;
        raw.insert(raw.end(), canvas.rgba.begin() + static_cast<std::ptrdiff_t>(begin),
                   canvas.rgba.begin() + static_cast<std::ptrdiff_t>(begin +
                       static_cast<std::size_t>(canvas.width) * 4u));
    }

    std::vector<std::uint8_t> zlib{0x78, 0x01};
    std::size_t position = 0;
    while (position < raw.size()) {
        const std::size_t block = std::min<std::size_t>(65535u, raw.size() - position);
        const bool final = position + block == raw.size();
        zlib.push_back(final ? 0x01 : 0x00);
        const std::uint16_t len = static_cast<std::uint16_t>(block);
        const std::uint16_t nlen = static_cast<std::uint16_t>(~len);
        zlib.push_back(static_cast<std::uint8_t>(len));
        zlib.push_back(static_cast<std::uint8_t>(len >> 8u));
        zlib.push_back(static_cast<std::uint8_t>(nlen));
        zlib.push_back(static_cast<std::uint8_t>(nlen >> 8u));
        zlib.insert(zlib.end(), raw.begin() + static_cast<std::ptrdiff_t>(position),
                    raw.begin() + static_cast<std::ptrdiff_t>(position + block));
        position += block;
    }
    push_u32(zlib, adler32(raw.data(), raw.size()));

    std::vector<std::uint8_t> png{137,80,78,71,13,10,26,10};
    std::vector<std::uint8_t> ihdr;
    push_u32(ihdr, static_cast<std::uint32_t>(canvas.width));
    push_u32(ihdr, static_cast<std::uint32_t>(canvas.height));
    ihdr.insert(ihdr.end(), {8, 6, 0, 0, 0});
    append_chunk(png, "IHDR", ihdr);
    append_chunk(png, "IDAT", zlib);
    append_chunk(png, "IEND", {});

    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open PNG output path";
        return false;
    }
    out.write(reinterpret_cast<const char*>(png.data()), static_cast<std::streamsize>(png.size()));
    return static_cast<bool>(out);
}

void svg_arrow(std::ostream& out, const double x, const double y, const bool up) {
    const double tip = y + (up ? -24.0 : 24.0);
    const double dy = up ? 6.0 : -6.0;
    out << "<line x1='" << x << "' y1='" << y << "' x2='" << x << "' y2='" << tip
        << "' stroke='#111827' stroke-width='2'/>\n";
    out << "<path d='M " << x - 5 << ' ' << tip + dy << " L " << x << ' ' << tip
        << " L " << x + 5 << ' ' << tip + dy
        << "' fill='none' stroke='#111827' stroke-width='2'/>\n";
}

void svg_energy_axis(std::ostream& out, const int height) {
    out << "<line x1='62' y1='" << height - 70
        << "' x2='62' y2='115' stroke='#111827' stroke-width='2'/>\n";
    out << "<path d='M56 125 L62 113 L68 125' fill='none' stroke='#111827' stroke-width='2'/>\n";
    out << "<text x='38' y='" << height / 2 << "' transform='rotate(-90 38 " << height / 2
        << ")' font-family='sans-serif' font-size='15' fill='#374151'>Energy</text>\n";
}

} // namespace

