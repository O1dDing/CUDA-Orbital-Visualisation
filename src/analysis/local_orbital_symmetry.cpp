#include "cov/local_orbital_symmetry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace cov {
namespace {

using Mat3 = std::array<double, 9>;
using Vec3 = std::array<double, 3>;
using ComponentWeights = std::array<double, 5>;

constexpr double kWeightEpsilon = 1.0e-14;
constexpr double kMinimumCopyConfidence = 0.55;

double determinant(const Mat3& matrix) noexcept {
    return matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7])
         - matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6])
         + matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
}

bool valid_rotation(const Mat3& rotation) noexcept {
    for (const double value : rotation) {
        if (!std::isfinite(value)) return false;
    }
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            double value = 0.0;
            for (std::size_t k = 0; k < 3; ++k) {
                value += rotation[3 * k + row] * rotation[3 * k + column];
            }
            const double expected = row == column ? 1.0 : 0.0;
            maximum_error = std::max(maximum_error, std::abs(value - expected));
        }
    }
    return maximum_error <= 1.0e-5 &&
           std::abs(std::abs(determinant(rotation)) - 1.0) <= 1.0e-5;
}

Vec3 transpose_multiply(const Mat3& matrix, const Vec3& vector) noexcept {
    return {
        matrix[0] * vector[0] + matrix[3] * vector[1] + matrix[6] * vector[2],
        matrix[1] * vector[0] + matrix[4] * vector[1] + matrix[7] * vector[2],
        matrix[2] * vector[0] + matrix[5] * vector[1] + matrix[8] * vector[2],
    };
}

Mat3 multiply(const Mat3& left, const Mat3& right) noexcept {
    Mat3 result{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            for (std::size_t k = 0; k < 3; ++k) {
                result[3 * row + column] +=
                    left[3 * row + k] * right[3 * k + column];
            }
        }
    }
    return result;
}

Mat3 transpose(const Mat3& matrix) noexcept {
    return {
        matrix[0], matrix[3], matrix[6],
        matrix[1], matrix[4], matrix[7],
        matrix[2], matrix[5], matrix[8],
    };
}

double frobenius_dot(const Mat3& left, const Mat3& right) noexcept {
    double result = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) result += left[i] * right[i];
    return result;
}

constexpr double kInvSqrt2 = 0.70710678118654752440;
constexpr double kInvSqrt6 = 0.40824829046386301637;

// Orthonormal tensor representatives of the five real d functions.  Their
// order is dz2, dxz, dyz, dx2-y2, dxy.  Gaussian pure-shell indices 1 and 2
// carry minus signs; those signs are applied when reading coefficients.
constexpr std::array<Mat3, 5> d_tensor_basis{{
    Mat3{-kInvSqrt6, 0.0, 0.0,
          0.0, -kInvSqrt6, 0.0,
          0.0, 0.0, 2.0 * kInvSqrt6},
    Mat3{0.0, 0.0, kInvSqrt2,
          0.0, 0.0, 0.0,
          kInvSqrt2, 0.0, 0.0},
    Mat3{0.0, 0.0, 0.0,
          0.0, 0.0, kInvSqrt2,
          0.0, kInvSqrt2, 0.0},
    Mat3{kInvSqrt2, 0.0, 0.0,
          0.0, -kInvSqrt2, 0.0,
          0.0, 0.0, 0.0},
    Mat3{0.0, kInvSqrt2, 0.0,
          kInvSqrt2, 0.0, 0.0,
          0.0, 0.0, 0.0},
}};

Mat3 pure_d_tensor(std::span<const float> coefficients) noexcept {
    Mat3 tensor{};
    constexpr std::array<double, 5> signs{1.0, -1.0, -1.0, 1.0, 1.0};
    for (std::size_t component = 0; component < 5; ++component) {
        const double coefficient = signs[component] * coefficients[component];
        for (std::size_t entry = 0; entry < tensor.size(); ++entry) {
            tensor[entry] += coefficient * d_tensor_basis[component][entry];
        }
    }
    return tensor;
}

Mat3 cartesian_d_tensor(std::span<const float> coefficients) noexcept {
    // Internal Cartesian order is xx, yy, zz, xy, xz, yz.  Diagonal basis
    // functions are x^2/sqrt(3), etc.; cross terms are xy, xz and yz.
    constexpr double kInvSqrt3 = 0.57735026918962576451;
    Mat3 tensor{
        coefficients[0] * kInvSqrt3, 0.5 * coefficients[3],
            0.5 * coefficients[4],
        0.5 * coefficients[3], coefficients[1] * kInvSqrt3,
            0.5 * coefficients[5],
        0.5 * coefficients[4], 0.5 * coefficients[5],
            coefficients[2] * kInvSqrt3,
    };
    const double trace_third = (tensor[0] + tensor[4] + tensor[8]) / 3.0;
    tensor[0] -= trace_third;
    tensor[4] -= trace_third;
    tensor[8] -= trace_third;
    return tensor;
}

Mat3 to_reference_tensor(const Mat3& tensor_input,
                         const Mat3& rotation_reference_to_input) noexcept {
    return multiply(
        multiply(transpose(rotation_reference_to_input), tensor_input),
        rotation_reference_to_input);
}

void accumulate_d_weights(ComponentWeights& weights,
                          const Mat3& tensor_reference) noexcept {
    for (std::size_t component = 0; component < d_tensor_basis.size(); ++component) {
        const double amplitude =
            frobenius_dot(tensor_reference, d_tensor_basis[component]);
        weights[component] += amplitude * amplitude;
    }
}

std::optional<std::size_t> basis_component(
    const MetalAOShell shell,
    std::string_view token) noexcept {
    if (shell == MetalAOShell::S) return token == "s" ? std::optional<std::size_t>{0}
                                                       : std::nullopt;
    if (shell == MetalAOShell::P) {
        if (token == "p_x") return 0;
        if (token == "p_y") return 1;
        if (token == "p_z") return 2;
        return std::nullopt;
    }
    if (token == "d_z2") return 0;
    if (token == "d_xz") return 1;
    if (token == "d_yz") return 2;
    if (token == "d_x2-y2") return 3;
    if (token == "d_xy") return 4;
    return std::nullopt;
}

double score_basis_functions(const ComponentWeights& weights,
                             const MetalAOShell shell,
                             std::string_view functions) noexcept {
    double score = 0.0;
    while (!functions.empty()) {
        const std::size_t separator = functions.find(',');
        std::string_view token = functions.substr(0, separator);
        while (!token.empty() && token.front() == ' ') token.remove_prefix(1);
        while (!token.empty() && token.back() == ' ') token.remove_suffix(1);
        if (const auto component = basis_component(shell, token)) {
            score += weights[*component];
        }
        if (separator == std::string_view::npos) break;
        functions.remove_prefix(separator + 1);
    }
    return score;
}

struct ShellEvidence {
    MetalAOShell shell = MetalAOShell::S;
    ComponentWeights components{};
    double weight = 0.0;
    double fraction = 0.0;
};

double evidence_floor(const MetalAOShell shell) noexcept {
    return shell == MetalAOShell::D ? 0.02 : 0.08;
}

std::optional<LocalIrrepAssignment> classify_copy(
    const ShellEvidence& evidence,
    std::string_view point_group) {
    const auto decomposition = decompose_metal_ao_shell(point_group, evidence.shell);
    if (!decomposition || decomposition->empty()) return std::nullopt;

    struct ScoredCopy {
        const IrrepCopy* copy = nullptr;
        double score = 0.0;
    };
    struct ScoredLabel {
        std::string_view label;
        double score = 0.0;
        std::vector<ScoredCopy> copies;
    };
    std::vector<ScoredLabel> labels;
    double candidate_sum = 0.0;
    for (const auto& copy : *decomposition) {
        const double score = score_basis_functions(
            evidence.components, evidence.shell, copy.basis_functions);
        candidate_sum += score;
        auto label=std::find_if(labels.begin(),labels.end(),[&](const auto& item) {
            return item.label==copy.label;
        });
        if (label==labels.end()) {
            labels.push_back(ScoredLabel{copy.label,score,{{&copy,score}}});
        } else {
            label->score+=score;
            label->copies.push_back({&copy,score});
        }
    }
    if (labels.empty() || candidate_sum <= kWeightEpsilon) return std::nullopt;
    const auto best_label=std::max_element(
        labels.begin(),labels.end(),[](const auto& left,const auto& right) {
            return left.score<right.score;
        });
    const double confidence = best_label->score / candidate_sum;
    if (!std::isfinite(confidence) || confidence < kMinimumCopyConfidence) {
        return std::nullopt;
    }

    const auto best_copy=std::max_element(
        best_label->copies.begin(),best_label->copies.end(),
        [](const auto& left,const auto& right) {return left.score<right.score;});
    const double copy_confidence=best_label->score>kWeightEpsilon
        ?best_copy->score/best_label->score:0.0;
    const bool copy_determined=std::isfinite(copy_confidence) &&
        copy_confidence>=kMinimumCopyConfidence;
    return LocalIrrepAssignment{
        std::string(best_label->label),evidence.shell,
        copy_determined?best_copy->copy->copy_index:std::uint8_t{0},
        confidence,copy_determined
            ?std::string(best_copy->copy->basis_functions):std::string{}};
}

} // namespace

std::optional<LocalIrrepAssignment> classify_local_irrep_by_dimension(
    const std::string_view point_group,
    const MetalAOShell shell,
    const std::size_t dimension) {
    if (dimension==0u) return std::nullopt;
    const auto decomposition=decompose_metal_ao_shell(point_group,shell);
    if (!decomposition) return std::nullopt;
    const IrrepCopy* match=nullptr;
    std::string_view label;
    bool repeated_copy=false;
    for (const auto& copy:*decomposition) {
        if (copy.dimension!=dimension) continue;
        if (match==nullptr) {
            match=&copy;
            label=copy.label;
        } else if (copy.label!=label) {
            return std::nullopt;
        } else {
            repeated_copy=true;
        }
    }
    if (match==nullptr) return std::nullopt;
    return LocalIrrepAssignment{
        std::string(match->label),shell,
        repeated_copy?std::uint8_t{0}:match->copy_index,1.0,
        repeated_copy?std::string{}:std::string(match->basis_functions)};
}

std::optional<LocalIrrepAssignment> classify_local_metal_irrep(
    const Wavefunction& wavefunction,
    std::span<const std::size_t> orbital_indices,
    const std::size_t metal_atom,
    const std::string_view point_group,
    const std::array<double, 9>& rotation_reference_to_input) {
    if (metal_atom >= wavefunction.atoms.size() || orbital_indices.empty() ||
        !find_point_group(point_group) ||
        !valid_rotation(rotation_reference_to_input)) {
        return std::nullopt;
    }
    for (const std::size_t orbital_index : orbital_indices) {
        if (orbital_index >= wavefunction.orbitals.size()) return std::nullopt;
    }

    std::array<ShellEvidence, 3> evidence{{
        {MetalAOShell::S, {}, 0.0, 0.0},
        {MetalAOShell::P, {}, 0.0, 0.0},
        {MetalAOShell::D, {}, 0.0, 0.0},
    }};

    for (const Shell& shell : wavefunction.shells) {
        if (shell.atom_index != metal_atom || shell.angular_momentum > 2) continue;
        const std::size_t count = shell_basis_count(shell);
        const std::size_t offset = shell.basis_offset;
        for (const std::size_t orbital_index : orbital_indices) {
            const auto& coefficients = wavefunction.orbitals[orbital_index].coefficients;
            if (offset > coefficients.size() || count > coefficients.size() - offset) {
                return std::nullopt;
            }
            const std::span<const float> local(coefficients.data() + offset, count);
            if (shell.angular_momentum == 0) {
                const double amplitude = local[0];
                evidence[0].components[0] += amplitude * amplitude;
            } else if (shell.angular_momentum == 1) {
                Vec3 input{};
                if (shell.pure) {
                    // Gaussian real-p order is pz, -px, -py.
                    input = {-local[1], -local[2], local[0]};
                } else {
                    input = {local[0], local[1], local[2]};
                }
                const Vec3 reference =
                    transpose_multiply(rotation_reference_to_input, input);
                for (std::size_t component = 0; component < 3; ++component) {
                    evidence[1].components[component] +=
                        reference[component] * reference[component];
                }
            } else {
                const Mat3 input = shell.pure ? pure_d_tensor(local)
                                              : cartesian_d_tensor(local);
                accumulate_d_weights(
                    evidence[2].components,
                    to_reference_tensor(input, rotation_reference_to_input));
            }
        }
    }

    double total_weight = 0.0;
    for (auto& shell_evidence : evidence) {
        shell_evidence.weight = 0.0;
        for (const double component : shell_evidence.components) {
            shell_evidence.weight += component;
        }
        total_weight += shell_evidence.weight;
    }
    if (total_weight <= kWeightEpsilon || !std::isfinite(total_weight)) {
        return std::nullopt;
    }
    for (auto& shell_evidence : evidence) {
        shell_evidence.fraction = shell_evidence.weight / total_weight;
    }
    std::stable_sort(
        evidence.begin(), evidence.end(),
        [](const ShellEvidence& left, const ShellEvidence& right) {
            return left.fraction > right.fraction;
        });

    // Prefer the dominant central-metal shell.  If its internal copy is
    // genuinely ambiguous, a lower-weight shell may still provide a clean
    // assignment, but only above the conservative shell-specific floor.
    for (const auto& shell_evidence : evidence) {
        if (shell_evidence.fraction + kWeightEpsilon <
            evidence_floor(shell_evidence.shell)) {
            continue;
        }
        if (auto assignment = classify_copy(shell_evidence, point_group)) {
            return assignment;
        }
    }
    return std::nullopt;
}

} // namespace cov
