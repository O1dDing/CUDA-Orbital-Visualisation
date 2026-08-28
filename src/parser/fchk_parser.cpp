#include "cov/fchk_parser.hpp"

#include "cov/density.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cov {
namespace {

std::string trim(std::string value) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
    return value;
}

double parse_real(std::string token) {
    for (char& c : token) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    std::size_t used = 0;
    const double value = std::stod(token, &used);
    if (used != token.size()) {
        throw std::runtime_error("Invalid FCHK real value: " + token);
    }
    return value;
}

long long parse_integer(const std::string& token) {
    std::size_t used = 0;
    const long long value = std::stoll(token, &used);
    if (used != token.size()) {
        throw std::runtime_error("Invalid FCHK integer value: " + token);
    }
    return value;
}

struct Header {
    std::string label;
    char type = '\0';
    bool array = false;
    std::size_t count = 0;
    std::string scalar;
};

bool parse_header(const std::string& line, Header& header) {
    if (line.size() < 41u) return false;

    header = {};
    header.label = trim(line.substr(0, 40));
    if (header.label.empty()) return false;

    const std::string tail = line.substr(40);
    std::istringstream stream(tail);
    stream >> header.type;
    if (!stream || (header.type != 'I' && header.type != 'R' &&
                    header.type != 'C' && header.type != 'L')) {
        return false;
    }

    const auto npos = tail.find("N=");
    if (npos != std::string::npos) {
        header.array = true;
        std::istringstream count_stream(tail.substr(npos + 2u));
        unsigned long long count = 0;
        count_stream >> count;
        if (!count_stream || count > std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("Invalid FCHK array length for field: " + header.label);
        }
        header.count = static_cast<std::size_t>(count);
    } else {
        std::string scalar;
        std::getline(stream, scalar);
        header.scalar = trim(scalar);
    }
    return true;
}

struct FchkRecords {
    std::unordered_map<std::string, long long> integers;
    std::unordered_map<std::string, double> reals;
    std::unordered_map<std::string, std::string> strings;
    std::unordered_map<std::string, std::vector<long long>> integer_arrays;
    std::unordered_map<std::string, std::vector<double>> real_arrays;
};

FchkRecords read_records(std::ifstream& input) {
    FchkRecords records;
    std::string line;
    while (std::getline(input, line)) {
        Header header;
        if (!parse_header(line, header)) continue;

        if (!header.array) {
            if (header.type == 'I') {
                records.integers[header.label] = parse_integer(header.scalar);
            } else if (header.type == 'R') {
                records.reals[header.label] = parse_real(header.scalar);
            } else {
                records.strings[header.label] = header.scalar;
            }
            continue;
        }

        // Numeric arrays are the authoritative wavefunction payload we need.
        // Character/logical arrays are intentionally ignored here; their data
        // lines naturally fail parse_header() and are skipped by the outer loop.
        if (header.type != 'I' && header.type != 'R') continue;

        std::vector<std::string> tokens;
        tokens.reserve(header.count);
        while (tokens.size() < header.count && std::getline(input, line)) {
            std::istringstream values(line);
            std::string token;
            while (values >> token) {
                tokens.push_back(token);
                if (tokens.size() == header.count) break;
            }
        }
        if (tokens.size() != header.count) {
            throw std::runtime_error("Unexpected end of FCHK array: " + header.label);
        }

        if (header.type == 'I') {
            auto& values = records.integer_arrays[header.label];
            values.reserve(header.count);
            for (const auto& token : tokens) values.push_back(parse_integer(token));
        } else {
            auto& values = records.real_arrays[header.label];
            values.reserve(header.count);
            for (const auto& token : tokens) values.push_back(parse_real(token));
        }
    }
    return records;
}

const std::vector<long long>& require_int_array(const FchkRecords& records,
                                                const std::string& label) {
    const auto it = records.integer_arrays.find(label);
    if (it == records.integer_arrays.end()) {
        throw std::runtime_error("Required FCHK field is missing: " + label);
    }
    return it->second;
}

const std::vector<double>& require_real_array(const FchkRecords& records,
                                              const std::string& label) {
    const auto it = records.real_arrays.find(label);
    if (it == records.real_arrays.end()) {
        throw std::runtime_error("Required FCHK field is missing: " + label);
    }
    return it->second;
}

long long require_integer(const FchkRecords& records, const std::string& label) {
    const auto it = records.integers.find(label);
    if (it == records.integers.end()) {
        throw std::runtime_error("Required FCHK field is missing: " + label);
    }
    return it->second;
}

const std::vector<double>* find_real_array(const FchkRecords& records,
                                           const std::string& label) {
    const auto it = records.real_arrays.find(label);
    return it == records.real_arrays.end() ? nullptr : &it->second;
}

std::string atomic_symbol(const int z) {
    static constexpr std::array<const char*, 119> symbols{{
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    }};
    if (z <= 0 || z >= static_cast<int>(symbols.size())) {
        return "X" + std::to_string(z);
    }
    return symbols[static_cast<std::size_t>(z)];
}

std::uint32_t checked_u32(const long long value, const std::string& label) {
    if (value < 0 || static_cast<unsigned long long>(value) >
                         std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("FCHK integer is out of range: " + label);
    }
    return static_cast<std::uint32_t>(value);
}

void append_shell(Wavefunction& wf,
                  const std::uint32_t atom_index,
                  const std::uint8_t l,
                  const bool pure,
                  const std::vector<double>& exponents,
                  const std::vector<double>& coefficients,
                  const std::size_t begin,
                  const std::size_t count) {
    if (begin + count > exponents.size() || begin + count > coefficients.size()) {
        throw std::runtime_error("FCHK primitive/contraction array is shorter than shell metadata");
    }

    Shell shell;
    shell.atom_index = atom_index;
    shell.primitive_offset = static_cast<std::uint32_t>(wf.primitives.size());
    shell.primitive_count = static_cast<std::uint32_t>(count);
    shell.basis_offset = wf.basis_count;
    shell.angular_momentum = l;
    shell.pure = pure ? 1u : 0u;

    for (std::size_t p = 0; p < count; ++p) {
        Primitive primitive;
        primitive.exponent = static_cast<float>(exponents[begin + p]);
        primitive.coefficient = static_cast<float>(coefficients[begin + p]);
        wf.primitives.push_back(primitive);
    }

    wf.basis_count += shell_basis_count(shell);
    wf.shells.push_back(shell);
}

void set_global_purity_flags(Wavefunction& wf) {
    auto all_pure = [&](const std::uint8_t l) {
        bool seen = false;
        bool pure = true;
        for (const Shell& shell : wf.shells) {
            if (shell.angular_momentum != l) continue;
            seen = true;
            pure = pure && shell.pure != 0;
        }
        return seen && pure;
    };
    wf.pure_d = all_pure(2);
    wf.pure_f = all_pure(3);
    wf.pure_g = all_pure(4);
}

void append_orbital_set(Wavefunction& wf,
                        const std::vector<double>& energies,
                        const std::vector<double>& coefficients,
                        const Spin spin,
                        const std::uint32_t occupied_count,
                        const bool restricted) {
    const std::size_t basis = wf.basis_count;
    if (basis == 0) throw std::runtime_error("FCHK orbital data has zero basis functions");
    if (coefficients.size() != energies.size() * basis) {
        throw std::runtime_error(
            "FCHK MO coefficient dimension does not equal orbital_count * basis_count");
    }

    wf.orbitals.reserve(wf.orbitals.size() + energies.size());
    for (std::size_t i = 0; i < energies.size(); ++i) {
        MolecularOrbital mo;
        mo.energy_hartree = energies[i];
        mo.spin = spin;
        if (restricted) {
            const float alpha = i < wf.alpha_electrons ? 1.0f : 0.0f;
            const float beta = i < wf.beta_electrons ? 1.0f : 0.0f;
            mo.occupation = alpha + beta;
        } else {
            mo.occupation = i < occupied_count ? 1.0f : 0.0f;
        }
        mo.coefficients.resize(basis);
        const std::size_t offset = i * basis;
        for (std::size_t j = 0; j < basis; ++j) {
            mo.coefficients[j] = static_cast<float>(coefficients[offset + j]);
        }
        wf.orbitals.push_back(std::move(mo));
    }
}

void validate_density(std::vector<double>& density,
                      const std::size_t basis,
                      const std::string& label) {
    if (density.empty()) return;
    const std::size_t expected = basis * (basis + 1u) / 2u;
    if (density.size() != expected) {
        throw std::runtime_error(label + " dimension does not match packed AO basis dimension");
    }
}

} // namespace

Wavefunction parse_fchk(const std::filesystem::path& path,
                        const FchkParseOptions& options) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Unable to open FCHK file: " + path.string());
    }

    Wavefunction wf;
    wf.source = WavefunctionSource::Fchk;
    std::getline(input, wf.source_title);
    std::getline(input, wf.source_route);
    wf.source_title = trim(wf.source_title);
    wf.source_route = trim(wf.source_route);

    const FchkRecords records = read_records(input);

    const auto atom_count_ll = require_integer(records, "Number of atoms");
    if (atom_count_ll <= 0) throw std::runtime_error("FCHK Number of atoms must be positive");
    const std::size_t atom_count = static_cast<std::size_t>(atom_count_ll);
    if (atom_count > options.max_atoms) {
        throw std::runtime_error("FCHK atom count exceeds configured maximum of " +
                                 std::to_string(options.max_atoms));
    }

    const auto& atomic_numbers = require_int_array(records, "Atomic numbers");
    const auto& coordinates = require_real_array(records, "Current cartesian coordinates");
    if (atomic_numbers.size() != atom_count || coordinates.size() != atom_count * 3u) {
        throw std::runtime_error("FCHK atom-coordinate dimensions are inconsistent");
    }

    wf.atoms.reserve(atom_count);
    for (std::size_t i = 0; i < atom_count; ++i) {
        Atom atom;
        atom.atomic_number = static_cast<int>(atomic_numbers[i]);
        atom.symbol = atomic_symbol(atom.atomic_number);
        atom.x = coordinates[3u * i + 0u];
        atom.y = coordinates[3u * i + 1u];
        atom.z = coordinates[3u * i + 2u];
        wf.atoms.push_back(std::move(atom));
    }

    wf.alpha_electrons = checked_u32(require_integer(records, "Number of alpha electrons"),
                                     "Number of alpha electrons");
    wf.beta_electrons = checked_u32(require_integer(records, "Number of beta electrons"),
                                    "Number of beta electrons");

    const auto& shell_types = require_int_array(records, "Shell types");
    const auto& primitive_counts = require_int_array(records, "Number of primitives per shell");
    const auto& shell_to_atom = require_int_array(records, "Shell to atom map");
    const auto& exponents = require_real_array(records, "Primitive exponents");
    const auto& contractions = require_real_array(records, "Contraction coefficients");
    const auto* sp_contractions = find_real_array(records, "P(S=P) Contraction coefficients");

    if (shell_types.size() != primitive_counts.size() ||
        shell_types.size() != shell_to_atom.size()) {
        throw std::runtime_error("FCHK shell metadata arrays have inconsistent lengths");
    }
    if (exponents.size() != contractions.size()) {
        throw std::runtime_error("FCHK primitive exponent/contraction arrays have inconsistent lengths");
    }

    std::size_t primitive_cursor = 0;
    for (std::size_t s = 0; s < shell_types.size(); ++s) {
        const long long count_ll = primitive_counts[s];
        if (count_ll <= 0) throw std::runtime_error("FCHK shell has no primitives");
        const std::size_t count = static_cast<std::size_t>(count_ll);
        if (primitive_cursor + count > exponents.size()) {
            throw std::runtime_error("FCHK shell primitive counts exceed primitive arrays");
        }

        const long long atom_map = shell_to_atom[s];
        if (atom_map <= 0 || static_cast<std::size_t>(atom_map) > wf.atoms.size()) {
            throw std::runtime_error("FCHK shell-to-atom map contains an invalid atom index");
        }
        const std::uint32_t atom_index = static_cast<std::uint32_t>(atom_map - 1);
        const int shell_type = static_cast<int>(shell_types[s]);

        if (shell_type == -1) {
            if (!sp_contractions || sp_contractions->size() != contractions.size()) {
                throw std::runtime_error(
                    "FCHK SP shell requires P(S=P) Contraction coefficients");
            }
            append_shell(wf, atom_index, 0u, false,
                         exponents, contractions, primitive_cursor, count);
            append_shell(wf, atom_index, 1u, false,
                         exponents, *sp_contractions, primitive_cursor, count);
        } else {
            const int l = shell_type < 0 ? -shell_type : shell_type;
            if (l < 0 || l > 4) {
                throw std::runtime_error(
                    "FCHK shell angular momentum above g is not supported yet: type=" +
                    std::to_string(shell_type));
            }
            const bool pure = shell_type <= -2;
            append_shell(wf, atom_index, static_cast<std::uint8_t>(l), pure,
                         exponents, contractions, primitive_cursor, count);
        }
        primitive_cursor += count;
    }
    if (primitive_cursor != exponents.size()) {
        throw std::runtime_error("FCHK primitive arrays contain unused entries");
    }

    const auto expected_basis = checked_u32(require_integer(records, "Number of basis functions"),
                                            "Number of basis functions");
    if (wf.basis_count != expected_basis) {
        throw std::runtime_error(
            "FCHK shell-derived basis count " + std::to_string(wf.basis_count) +
            " does not match Number of basis functions " + std::to_string(expected_basis));
    }
    set_global_purity_flags(wf);

    const auto* alpha_energies = find_real_array(records, "Alpha Orbital Energies");
    const auto* alpha_coefficients = find_real_array(records, "Alpha MO coefficients");
    const auto* beta_energies = find_real_array(records, "Beta Orbital Energies");
    const auto* beta_coefficients = find_real_array(records, "Beta MO coefficients");

    if ((alpha_energies == nullptr) != (alpha_coefficients == nullptr)) {
        throw std::runtime_error("FCHK alpha orbital energies/coefficients are incomplete");
    }
    if ((beta_energies == nullptr) != (beta_coefficients == nullptr)) {
        throw std::runtime_error("FCHK beta orbital energies/coefficients are incomplete");
    }

    if (alpha_energies && alpha_coefficients) {
        const bool unrestricted = beta_energies && beta_coefficients;
        append_orbital_set(wf, *alpha_energies, *alpha_coefficients,
                           Spin::Alpha, wf.alpha_electrons, !unrestricted);
        if (unrestricted) {
            append_orbital_set(wf, *beta_energies, *beta_coefficients,
                               Spin::Beta, wf.beta_electrons, false);
        }
    }

    if (options.require_orbitals && wf.orbitals.empty()) {
        throw std::runtime_error("No molecular orbitals found in FCHK file");
    }

    if (options.keep_density) {
        if (const auto* total = find_real_array(records, "Total SCF Density")) {
            wf.total_density_packed = *total;
            validate_density(wf.total_density_packed, wf.basis_count, "Total SCF Density");
            wf.total_density_provenance = DataProvenance::Producer;
        } else if (options.reconstruct_density_if_missing && !wf.orbitals.empty()) {
            wf.total_density_packed = reconstruct_total_density_packed(wf);
            wf.total_density_provenance = DataProvenance::Derived;
        }

        if (const auto* spin = find_real_array(records, "Spin SCF Density")) {
            wf.spin_density_packed = *spin;
            validate_density(wf.spin_density_packed, wf.basis_count, "Spin SCF Density");
            wf.spin_density_provenance = DataProvenance::Producer;
        }
    }

    return wf;
}

} // namespace cov
