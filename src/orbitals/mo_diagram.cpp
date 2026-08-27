#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string_view>
#include <vector>

namespace cov {
namespace {

std::string json_escape(std::string_view value) {
    std::ostringstream out;
    for (const unsigned char c : value) {
        switch (c) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<int>(c) << std::dec << std::setfill(' ');
                } else {
                    out << static_cast<char>(c);
                }
        }
    }
    return out.str();
}

std::string csv_escape(const std::string& value) {
    if (value.find_first_of(",\"\n\r") == std::string::npos) return value;
    std::string escaped = "\"";
    for (const char c : value) {
        if (c == '"') escaped += "\"\"";
        else escaped += c;
    }
    escaped += '"';
    return escaped;
}

const char* spin_name(const Spin spin) noexcept {
    return spin == Spin::Beta ? "Beta" : "Alpha";
}

const char* region_name(const OrbitalRegion region) noexcept {
    switch (region) {
        case OrbitalRegion::Core: return "core";
        case OrbitalRegion::Valence: return "valence";
        default: return "virtual";
    }
}

const char* classification_name(const DiagramClassification classification) noexcept {
    switch (classification) {
        case DiagramClassification::SymmetryGrouped: return "symmetry-grouped";
        case DiagramClassification::SalcUnavailable: return "salc-unavailable";
        default: return "simple";
    }
}

std::vector<std::size_t> choose_level_indices(const std::vector<OrbitalMetadata>& metadata,
                                              const std::size_t selected,
                                              const std::size_t neighbourhood) {
    std::vector<std::size_t> visible;
    visible.reserve(metadata.size());
    for (std::size_t i = 0; i < metadata.size(); ++i) {
        if (metadata[i].visible) visible.push_back(i);
    }
    if (visible.empty()) return visible;
    if (neighbourhood == 0 || visible.size() <= neighbourhood * 2 + 1) return visible;

    auto it = std::lower_bound(visible.begin(), visible.end(), selected);
    std::size_t pos = it == visible.end() ? visible.size() - 1
                                          : static_cast<std::size_t>(it - visible.begin());
    if (it != visible.end() && *it != selected && pos > 0) {
        const std::size_t left = visible[pos - 1];
        const std::size_t right = visible[pos];
        if (selected - left <= right - selected) --pos;
    }

    const std::size_t begin = pos > neighbourhood ? pos - neighbourhood : 0;
    const std::size_t end = std::min(visible.size(), pos + neighbourhood + 1);
    return std::vector<std::size_t>(visible.begin() + static_cast<std::ptrdiff_t>(begin),
                                    visible.begin() + static_cast<std::ptrdiff_t>(end));
}

struct DiagramLayout {
    double e_min = -1.0;
    double e_max = 1.0;
    std::vector<std::string> columns;
    std::map<std::string, int> column_index;
};

DiagramLayout make_layout(const MODiagramData& data) {
    DiagramLayout layout;
    if (!data.levels.empty()) {
        layout.e_min = data.levels.front().metadata.energy_hartree;
        layout.e_max = layout.e_min;
        for (const auto& level : data.levels) {
            layout.e_min = std::min(layout.e_min, level.metadata.energy_hartree);
            layout.e_max = std::max(layout.e_max, level.metadata.energy_hartree);
        }
        if (std::abs(layout.e_max - layout.e_min) < 1.0e-9) {
            layout.e_min -= 0.5;
            layout.e_max += 0.5;
        } else {
            const double pad = (layout.e_max - layout.e_min) * 0.08;
            layout.e_min -= pad;
            layout.e_max += pad;
        }
    }

    if (data.plan.classification == DiagramClassification::SymmetryGrouped) {
        std::set<std::string> distinct;
        for (const auto& level : data.levels) {
            distinct.insert(level.metadata.symmetry.empty() ? "?" : level.metadata.symmetry);
        }
        layout.columns.assign(distinct.begin(), distinct.end());
    }
    if (layout.columns.empty()) layout.columns.push_back("MO");
    for (std::size_t i = 0; i < layout.columns.size(); ++i) {
        layout.column_index[layout.columns[i]] = static_cast<int>(i);
    }
    return layout;
}

double level_y(const double energy_ha,
               const DiagramLayout& layout,
               const int height,
               const int top = 90,
               const int bottom = 70) {
    const double t = (energy_ha - layout.e_min) /
                     std::max(1.0e-12, layout.e_max - layout.e_min);
    return static_cast<double>(height - bottom) -
           t * static_cast<double>(height - top - bottom);
}

double column_x(const MODiagramLevel& level,
                const DiagramLayout& layout,
                const int width,
                const int left = 120,
                const int right = 120) {
    std::string key = "MO";
    if (layout.columns.size() > 1 || layout.columns.front() != "MO") {
        key = level.metadata.symmetry.empty() ? "?" : level.metadata.symmetry;
    }
    const auto it = layout.column_index.find(key);
    const int index = it == layout.column_index.end() ? 0 : it->second;
    const int count = static_cast<int>(layout.columns.size());
    const double span = static_cast<double>(width - left - right);
    return static_cast<double>(left) +
           (static_cast<double>(index) + 0.5) * span / static_cast<double>(std::max(1, count));
}

std::string xml_escape(std::string_view value) {
    std::string out;
    out.reserve(value.size());
    for (const char c : value) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out += c; break;
        }
    }
    return out;
}

struct Canvas {
    int width = 0;
    int height = 0;
    std::vector<std::uint8_t> rgba;

    Canvas(const int w, const int h)
        : width(std::max(1, w)), height(std::max(1, h)),
          rgba(static_cast<std::size_t>(std::max(1, w)) *
               static_cast<std::size_t>(std::max(1, h)) * 4u, 255u) {
        clear(14, 20, 29, 255);
    }

    void clear(const std::uint8_t r, const std::uint8_t g,
               const std::uint8_t b, const std::uint8_t a) {
        for (std::size_t i = 0; i < rgba.size(); i += 4) {
            rgba[i] = r; rgba[i + 1] = g; rgba[i + 2] = b; rgba[i + 3] = a;
        }
    }

    void pixel(const int x, const int y,
               const std::uint8_t r, const std::uint8_t g,
               const std::uint8_t b, const std::uint8_t a = 255) {
        if (x < 0 || y < 0 || x >= width || y >= height) return;
        const std::size_t i = (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                               static_cast<std::size_t>(x)) * 4u;
        rgba[i] = r; rgba[i + 1] = g; rgba[i + 2] = b; rgba[i + 3] = a;
    }

    void line(int x0, int y0, const int x1, const int y1,
              const std::uint8_t r, const std::uint8_t g,
              const std::uint8_t b, const int thickness = 1) {
        const int dx = std::abs(x1 - x0);
        const int sx = x0 < x1 ? 1 : -1;
        const int dy = -std::abs(y1 - y0);
        const int sy = y0 < y1 ? 1 : -1;
        int err = dx + dy;
        for (;;) {
            const int radius = std::max(0, thickness / 2);
            for (int oy = -radius; oy <= radius; ++oy) {
                for (int ox = -radius; ox <= radius; ++ox) pixel(x0 + ox, y0 + oy, r, g, b);
            }
            if (x0 == x1 && y0 == y1) break;
            const int e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
    }
};

// Compact 5x7 bitmap font. Lower-case letters intentionally reuse upper-case
// shapes: the PNG is a readable overview, while exact UTF-8/symmetry metadata
// remains in SVG/JSON/CSV.
std::array<std::uint8_t, 7> glyph_rows(char c) {
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    switch (c) {
        case '0': return {14,17,19,21,25,17,14};
        case '1': return {4,12,4,4,4,4,14};
        case '2': return {14,17,1,2,4,8,31};
        case '3': return {30,1,1,14,1,1,30};
        case '4': return {2,6,10,18,31,2,2};
        case '5': return {31,16,16,30,1,1,30};
        case '6': return {14,16,16,30,17,17,14};
        case '7': return {31,1,2,4,8,8,8};
        case '8': return {14,17,17,14,17,17,14};
        case '9': return {14,17,17,15,1,1,14};
        case 'A': return {14,17,17,31,17,17,17};
        case 'B': return {30,17,17,30,17,17,30};
        case 'C': return {14,17,16,16,16,17,14};
        case 'D': return {30,17,17,17,17,17,30};
        case 'E': return {31,16,16,30,16,16,31};
        case 'F': return {31,16,16,30,16,16,16};
        case 'G': return {14,17,16,23,17,17,15};
        case 'H': return {17,17,17,31,17,17,17};
        case 'I': return {14,4,4,4,4,4,14};
        case 'J': return {7,2,2,2,2,18,12};
        case 'K': return {17,18,20,24,20,18,17};
        case 'L': return {16,16,16,16,16,16,31};
        case 'M': return {17,27,21,21,17,17,17};
        case 'N': return {17,25,21,19,17,17,17};
        case 'O': return {14,17,17,17,17,17,14};
        case 'P': return {30,17,17,30,16,16,16};
        case 'Q': return {14,17,17,17,21,18,13};
        case 'R': return {30,17,17,30,20,18,17};
        case 'S': return {15,16,16,14,1,1,30};
        case 'T': return {31,4,4,4,4,4,4};
        case 'U': return {17,17,17,17,17,17,14};
        case 'V': return {17,17,17,17,17,10,4};
        case 'W': return {17,17,17,21,21,21,10};
        case 'X': return {17,17,10,4,10,17,17};
        case 'Y': return {17,17,10,4,4,4,4};
        case 'Z': return {31,1,2,4,8,16,31};
        case '-': return {0,0,0,31,0,0,0};
        case '.': return {0,0,0,0,0,12,12};
        case '/': return {1,2,2,4,8,8,16};
        case '+': return {0,4,4,31,4,4,0};
        case ':': return {0,12,12,0,12,12,0};
        case '?': return {14,17,1,2,4,0,4};
        case '_': return {0,0,0,0,0,0,31};
        default: return {0,0,0,0,0,0,0};
    }
}

void draw_text(Canvas& canvas, int x, const int y, const std::string& text,
               const std::uint8_t r = 210, const std::uint8_t g = 220,
               const std::uint8_t b = 234, const int scale = 2) {
    for (const char c : text) {
        if (c == ' ') {
            x += 4 * scale;
            continue;
        }
        const auto rows = glyph_rows(c);
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if ((rows[static_cast<std::size_t>(row)] & (1u << (4 - col))) == 0) continue;
                for (int yy = 0; yy < scale; ++yy) {
                    for (int xx = 0; xx < scale; ++xx) {
                        canvas.pixel(x + col * scale + xx,
                                     y + row * scale + yy, r, g, b);
                    }
                }
            }
        }
        x += 6 * scale;
    }
}

void draw_electron(Canvas& canvas, const int x, const int y,
                   const bool up, const std::uint8_t r,
                   const std::uint8_t g, const std::uint8_t b) {
    const int tip = up ? y - 17 : y + 17;
    canvas.line(x, y, x, tip, r, g, b, 2);
    const int dy = up ? 5 : -5;
    canvas.line(x, tip, x - 4, tip + dy, r, g, b, 2);
    canvas.line(x, tip, x + 4, tip + dy, r, g, b, 2);
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
    constexpr std::uint32_t mod = 65521u;
    std::uint32_t a = 1u;
    std::uint32_t b = 0u;
    for (std::size_t i = 0; i < size; ++i) {
        a = (a + data[i]) % mod;
        b = (b + a) % mod;
    }
    return (b << 16u) | a;
}

void push_u32_be(std::vector<std::uint8_t>& out, const std::uint32_t value) {
    out.push_back(static_cast<std::uint8_t>((value >> 24u) & 0xffu));
    out.push_back(static_cast<std::uint8_t>((value >> 16u) & 0xffu));
    out.push_back(static_cast<std::uint8_t>((value >> 8u) & 0xffu));
    out.push_back(static_cast<std::uint8_t>(value & 0xffu));
}

void append_png_chunk(std::vector<std::uint8_t>& png,
                      const char type[4],
                      const std::vector<std::uint8_t>& payload) {
    push_u32_be(png, static_cast<std::uint32_t>(payload.size()));
    const std::size_t crc_begin = png.size();
    png.insert(png.end(), type, type + 4);
    png.insert(png.end(), payload.begin(), payload.end());
    push_u32_be(png, crc32(png.data() + crc_begin, 4 + payload.size()));
}

bool write_png_rgba(const Canvas& canvas,
                    const std::filesystem::path& path,
                    std::string* error) {
    std::vector<std::uint8_t> raw;
    raw.reserve(static_cast<std::size_t>(canvas.height) *
                (1u + static_cast<std::size_t>(canvas.width) * 4u));
    for (int y = 0; y < canvas.height; ++y) {
        raw.push_back(0); // PNG filter: None
        const std::size_t offset = static_cast<std::size_t>(y) *
                                   static_cast<std::size_t>(canvas.width) * 4u;
        raw.insert(raw.end(), canvas.rgba.begin() + static_cast<std::ptrdiff_t>(offset),
                   canvas.rgba.begin() + static_cast<std::ptrdiff_t>(offset +
                       static_cast<std::size_t>(canvas.width) * 4u));
    }

    std::vector<std::uint8_t> zlib;
    zlib.push_back(0x78);
    zlib.push_back(0x01); // valid zlib header, fastest/no compression
    std::size_t pos = 0;
    while (pos < raw.size()) {
        const std::size_t block = std::min<std::size_t>(65535u, raw.size() - pos);
        const bool final = pos + block == raw.size();
        zlib.push_back(final ? 0x01 : 0x00); // stored DEFLATE block
        const std::uint16_t len = static_cast<std::uint16_t>(block);
        const std::uint16_t nlen = static_cast<std::uint16_t>(~len);
        zlib.push_back(static_cast<std::uint8_t>(len & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>((len >> 8u) & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>(nlen & 0xffu));
        zlib.push_back(static_cast<std::uint8_t>((nlen >> 8u) & 0xffu));
        zlib.insert(zlib.end(), raw.begin() + static_cast<std::ptrdiff_t>(pos),
                    raw.begin() + static_cast<std::ptrdiff_t>(pos + block));
        pos += block;
    }
    push_u32_be(zlib, adler32(raw.data(), raw.size()));

    std::vector<std::uint8_t> png = {137,80,78,71,13,10,26,10};
    std::vector<std::uint8_t> ihdr;
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.width));
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.height));
    ihdr.push_back(8); // bit depth
    ihdr.push_back(6); // RGBA
    ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0);
    append_png_chunk(png, "IHDR", ihdr);
    append_png_chunk(png, "IDAT", zlib);
    append_png_chunk(png, "IEND", {});

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        if (error) *error = "Unable to open PNG output path";
        return false;
    }
    file.write(reinterpret_cast<const char*>(png.data()),
               static_cast<std::streamsize>(png.size()));
    if (!file) {
        if (error) *error = "Failed while writing PNG";
        return false;
    }
    return true;
}

} // namespace

MODiagramData build_mo_diagram_data(const Wavefunction& wavefunction,
                                    const MODiagramOptions& options) {
    MODiagramData data;
    data.plan = choose_diagram_plan(wavefunction);
    data.frontier = find_frontier_orbitals(wavefunction.orbitals,
                                           options.filter.occupation_threshold);
    data.metadata = build_orbital_metadata(wavefunction,
                                           options.selected_index,
                                           options.degeneracy,
                                           options.filter);

    const auto chosen = choose_level_indices(data.metadata,
                                             options.selected_index,
                                             options.neighbourhood);
    data.levels.reserve(chosen.size());
    for (const std::size_t index : chosen) {
        MODiagramLevel level;
        level.metadata = data.metadata[index];
        level.electrons = electron_glyphs_for_orbital(
            wavefunction.orbitals[index], data.frontier.separate_spin_sets);
        level.homo = data.frontier.homo && *data.frontier.homo == index;
        level.lumo = data.frontier.lumo && *data.frontier.lumo == index;
        data.levels.push_back(std::move(level));
    }
    return data;
}

bool write_mo_diagram_svg(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open SVG output path";
        return false;
    }

    const int width = std::max(640, options.width);
    const int height = std::max(480, options.height);
    const DiagramLayout layout = make_layout(data);

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << ' ' << height << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"#0e141d\"/>\n";
    out << "<text x=\"32\" y=\"36\" fill=\"#e5edf8\" font-family=\"sans-serif\" font-size=\"22\">"
        << "MO diagram</text>\n";
    out << "<text x=\"32\" y=\"60\" fill=\"#8fa1b8\" font-family=\"sans-serif\" font-size=\"13\">"
        << xml_escape(classification_name(data.plan.classification)) << " · "
        << xml_escape(data.plan.machine_reason) << "</text>\n";

    if (data.plan.classification == DiagramClassification::SalcUnavailable) {
        out << "<text x=\"32\" y=\"82\" fill=\"#e7b45b\" font-family=\"sans-serif\" font-size=\"13\">"
            << "SALC not confidently available; showing reliable orbital ordering.</text>\n";
    }

    const double column_span = static_cast<double>(width - 240) /
                               static_cast<double>(std::max<std::size_t>(1, layout.columns.size()));
    for (std::size_t i = 0; i < layout.columns.size(); ++i) {
        const double x = 120.0 + (static_cast<double>(i) + 0.5) * column_span;
        out << "<text x=\"" << x << "\" y=\"82\" text-anchor=\"middle\" fill=\"#8fa1b8\" "
            << "font-family=\"sans-serif\" font-size=\"13\">"
            << xml_escape(layout.columns[i]) << "</text>\n";
    }

    for (const auto& level : data.levels) {
        const double y = level_y(level.metadata.energy_hartree, layout, height);
        double x = column_x(level, layout, width);
        if (level.metadata.degeneracy_size > 1) {
            const std::string& label = level.metadata.display_label;
            const auto dash = label.find('-');
            if (dash != std::string::npos && dash + 1 < label.size()) {
                const char suffix = label[dash + 1];
                if (suffix >= 'a' && suffix <= 'z') {
                    x += (static_cast<int>(suffix - 'a') -
                          (static_cast<int>(level.metadata.degeneracy_size) - 1) * 0.5) * 14.0;
                }
            }
        }

        const char* colour = level.metadata.selected ? "#69a6ff" :
                             level.homo ? "#5fd6a2" :
                             level.lumo ? "#e7b45b" :
                             level.metadata.region == OrbitalRegion::Virtual ? "#8693a5" : "#c8d5e6";
        const double x0 = x - 38.0;
        const double x1 = x + 38.0;
        out << "<line x1=\"" << x0 << "\" y1=\"" << y << "\" x2=\"" << x1
            << "\" y2=\"" << y << "\" stroke=\"" << colour << "\" stroke-width=\""
            << (level.metadata.selected ? 4 : 2) << "\"/>\n";

        if (level.electrons.alpha > 0) {
            out << "<text x=\"" << (x - 12) << "\" y=\"" << (y - 5)
                << "\" fill=\"#eaf2ff\" font-family=\"sans-serif\" font-size=\"20\">↑</text>\n";
        }
        if (level.electrons.beta > 0) {
            out << "<text x=\"" << (x + 4) << "\" y=\"" << (y - 5)
                << "\" fill=\"#eaf2ff\" font-family=\"sans-serif\" font-size=\"20\">↓</text>\n";
        }

        const std::string energy = format_energy(level.metadata.energy_hartree,
                                                 options.energy_unit,
                                                 options.energy_unit == EnergyUnit::JoulePerMol ? 1 : 5);
        out << "<text x=\"" << (x + 48) << "\" y=\"" << (y + 4)
            << "\" fill=\"#d4dfed\" font-family=\"sans-serif\" font-size=\"12\">"
            << xml_escape(level.metadata.display_label) << " · " << xml_escape(energy)
            << "</text>\n";
    }

    out << "</svg>\n";
    if (!out) {
        if (error) *error = "Failed while writing SVG";
        return false;
    }
    return true;
}

bool write_mo_diagram_png(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    const int width = std::max(640, options.width);
    const int height = std::max(480, options.height);
    Canvas canvas(width, height);
    const DiagramLayout layout = make_layout(data);

    draw_text(canvas, 28, 24, "MO DIAGRAM", 230, 238, 249, 3);
    draw_text(canvas, 28, 52,
              std::string("MODE: ") + classification_name(data.plan.classification),
              142, 161, 184, 2);
    if (data.plan.classification == DiagramClassification::SalcUnavailable) {
        draw_text(canvas, 28, 72, "SALC NOT CONFIDENTLY AVAILABLE", 231, 180, 91, 2);
    }

    for (const auto& level : data.levels) {
        const int y = static_cast<int>(std::lround(level_y(level.metadata.energy_hartree,
                                                           layout, height)));
        int x = static_cast<int>(std::lround(column_x(level, layout, width)));
        if (level.metadata.degeneracy_size > 1) {
            const auto dash = level.metadata.display_label.find('-');
            if (dash != std::string::npos && dash + 1 < level.metadata.display_label.size()) {
                const char suffix = level.metadata.display_label[dash + 1];
                if (suffix >= 'a' && suffix <= 'z') {
                    x += static_cast<int>(std::lround(
                        (static_cast<int>(suffix - 'a') -
                         (static_cast<int>(level.metadata.degeneracy_size) - 1) * 0.5) * 14.0));
                }
            }
        }

        std::array<std::uint8_t, 3> colour{200, 213, 230};
        if (level.metadata.region == OrbitalRegion::Virtual) colour = {134, 147, 165};
        if (level.homo) colour = {95, 214, 162};
        if (level.lumo) colour = {231, 180, 91};
        if (level.metadata.selected) colour = {105, 166, 255};
        canvas.line(x - 38, y, x + 38, y, colour[0], colour[1], colour[2],
                    level.metadata.selected ? 4 : 2);

        if (level.electrons.alpha > 0) draw_electron(canvas, x - 9, y - 2, true, 236, 243, 252);
        if (level.electrons.beta > 0) draw_electron(canvas, x + 9, y - 2, false, 236, 243, 252);

        draw_text(canvas, x + 48, y - 7, level.metadata.display_label,
                  214, 225, 239, 2);
    }

    return write_png_rgba(canvas, path, error);
}

bool write_mo_diagram_json(const MODiagramData& data,
                           const MODiagramOptions& options,
                           const std::filesystem::path& path,
                           std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open JSON output path";
        return false;
    }

    out << "{\n";
    out << "  \"schema\": \"cov.mo-diagram.v1\",\n";
    out << "  \"classification\": \"" << classification_name(data.plan.classification) << "\",\n";
    out << "  \"strict_salc_available\": " << (data.plan.strict_salc_available ? "true" : "false") << ",\n";
    out << "  \"classification_reason\": \"" << json_escape(data.plan.machine_reason) << "\",\n";
    out << "  \"energy_unit\": \"" << energy_unit_symbol(options.energy_unit) << "\",\n";
    out << "  \"degeneracy_tolerance_hartree\": " << std::setprecision(12)
        << options.degeneracy.tolerance_hartree << ",\n";
    out << "  \"filter\": {\n";
    out << "    \"virtual_window_hartree\": " << options.filter.virtual_window_hartree << ",\n";
    out << "    \"core_energy_cutoff_hartree\": " << options.filter.core_energy_cutoff_hartree << "\n";
    out << "  },\n";
    out << "  \"orbitals\": [\n";

    const auto& metadata = data.metadata;
    for (std::size_t i = 0; i < metadata.size(); ++i) {
        const auto& item = metadata[i];
        out << "    {\"index\": " << item.orbital_index
            << ", \"raw_mo\": " << item.raw_mo_number
            << ", \"label\": \"" << json_escape(item.display_label) << "\""
            << ", \"energy_hartree\": " << std::setprecision(15) << item.energy_hartree
            << ", \"energy_display\": " << convert_hartree(item.energy_hartree, options.energy_unit)
            << ", \"occupation\": " << item.occupation
            << ", \"spin\": \"" << spin_name(item.spin) << "\""
            << ", \"symmetry\": \"" << json_escape(item.symmetry) << "\""
            << ", \"degeneracy_size\": " << item.degeneracy_size
            << ", \"region\": \"" << region_name(item.region) << "\""
            << ", \"visible\": " << (item.visible ? "true" : "false")
            << ", \"selected\": " << (item.selected ? "true" : "false") << '}';
        if (i + 1 != metadata.size()) out << ',';
        out << '\n';
    }
    out << "  ]\n}\n";
    if (!out) {
        if (error) *error = "Failed while writing JSON";
        return false;
    }
    return true;
}

bool write_mo_diagram_csv(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open CSV output path";
        return false;
    }

    out << "internal_index,raw_mo,display_label,energy_hartree,energy_display,energy_unit,"
           "occupation,spin,symmetry,degeneracy_size,region,visible,selected\n";
    for (const auto& item : data.metadata) {
        out << item.orbital_index << ',' << item.raw_mo_number << ','
            << csv_escape(item.display_label) << ',' << std::setprecision(15)
            << item.energy_hartree << ','
            << convert_hartree(item.energy_hartree, options.energy_unit) << ','
            << energy_unit_symbol(options.energy_unit) << ',' << item.occupation << ','
            << spin_name(item.spin) << ',' << csv_escape(item.symmetry) << ','
            << item.degeneracy_size << ',' << region_name(item.region) << ','
            << (item.visible ? "true" : "false") << ','
            << (item.selected ? "true" : "false") << '\n';
    }
    if (!out) {
        if (error) *error = "Failed while writing CSV";
        return false;
    }
    return true;
}

MODiagramExportResult export_mo_diagram_bundle(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::filesystem::path& base_path) {
    MODiagramExportResult result;
    const MODiagramData data = build_mo_diagram_data(wavefunction, options);

    std::filesystem::path base = base_path;
    if (base.has_extension()) base.replace_extension();
    std::string error;

    result.svg = write_mo_diagram_svg(data, options, base.string() + ".mo.svg", &error);
    if (!result.svg && result.error.empty()) result.error = error;
    error.clear();
    result.png = write_mo_diagram_png(data, options, base.string() + ".mo.png", &error);
    if (!result.png && result.error.empty()) result.error = error;
    error.clear();
    result.json = write_mo_diagram_json(data, options, base.string() + ".mo.json", &error);
    if (!result.json && result.error.empty()) result.error = error;
    error.clear();
    result.csv = write_mo_diagram_csv(data, options, base.string() + ".mo.csv", &error);
    if (!result.csv && result.error.empty()) result.error = error;

    return result;
}

} // namespace cov
