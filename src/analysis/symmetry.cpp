#include "cov/symmetry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace cov {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<double, 9>;

struct RotationAxis {
    Vec3 axis{};
    int order = 1;
    double error = 0.0;
};

struct ImproperAxis {
    Vec3 axis{};
    int order = 1;
    double error = 0.0;
};

struct ReflectionPlane {
    Vec3 normal{};
    double error = 0.0;
};

Vec3 add(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 sub(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 scale(const Vec3& a, const double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    };
}

double norm2(const Vec3& a) {
    return dot(a, a);
}

double norm(const Vec3& a) {
    return std::sqrt(norm2(a));
}

bool normalize(Vec3& a) {
    const double n = norm(a);
    if (!(n > 1.0e-14) || !std::isfinite(n)) return false;
    a = scale(a, 1.0 / n);
    return true;
}

void canonicalize_axis(Vec3& a) {
    if (!normalize(a)) return;
    for (double value : a) {
        if (std::abs(value) <= 1.0e-12) continue;
        if (value < 0.0) a = scale(a, -1.0);
        break;
    }
}

Vec3 mat_vec(const Mat3& m, const Vec3& v) {
    return {
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    };
}

Mat3 mat_mul(const Mat3& a, const Mat3& b) {
    Mat3 out{};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            double value = 0.0;
            for (int k = 0; k < 3; ++k) value += a[3 * r + k] * b[3 * k + c];
            out[3 * r + c] = value;
        }
    }
    return out;
}

Mat3 identity_matrix() {
    return {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    };
}

Mat3 inversion_matrix() {
    return {
        -1.0, 0.0, 0.0,
         0.0,-1.0, 0.0,
         0.0, 0.0,-1.0,
    };
}

Mat3 rotation_matrix(Vec3 axis, const double angle) {
    normalize(axis);
    const double x = axis[0];
    const double y = axis[1];
    const double z = axis[2];
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double t = 1.0 - c;
    return {
        t*x*x + c,     t*x*y - s*z,   t*x*z + s*y,
        t*x*y + s*z,   t*y*y + c,     t*y*z - s*x,
        t*x*z - s*y,   t*y*z + s*x,   t*z*z + c,
    };
}

Mat3 reflection_matrix(Vec3 normal) {
    normalize(normal);
    const double x = normal[0];
    const double y = normal[1];
    const double z = normal[2];
    return {
        1.0 - 2.0*x*x, -2.0*x*y,       -2.0*x*z,
        -2.0*y*x,       1.0 - 2.0*y*y, -2.0*y*z,
        -2.0*z*x,       -2.0*z*y,       1.0 - 2.0*z*z,
    };
}

Mat3 improper_matrix(const Vec3& axis, const int order) {
    const Mat3 r = rotation_matrix(axis, 2.0 * kPi / static_cast<double>(order));
    const Mat3 h = reflection_matrix(axis);
    return mat_mul(h, r);
}

bool same_axis(const Vec3& a, const Vec3& b, const double cosine_tolerance = 1.0e-7) {
    return std::abs(dot(a, b)) >= 1.0 - cosine_tolerance;
}

void add_axis(std::vector<Vec3>& axes,
              Vec3 axis,
              const std::size_t maximum) {
    if (axes.size() >= maximum) return;
    if (!normalize(axis)) return;
    canonicalize_axis(axis);
    for (const auto& existing : axes) {
        if (same_axis(axis, existing, 2.0e-7)) return;
    }
    axes.push_back(axis);
}

Vec3 molecular_centre(const Wavefunction& wf) {
    Vec3 centre{0.0, 0.0, 0.0};
    double weight_sum = 0.0;
    for (const auto& atom : wf.atoms) {
        const double weight = static_cast<double>(std::max(1, atom.atomic_number));
        centre[0] += weight * atom.x;
        centre[1] += weight * atom.y;
        centre[2] += weight * atom.z;
        weight_sum += weight;
    }
    if (weight_sum > 0.0) centre = scale(centre, 1.0 / weight_sum);
    return centre;
}

std::vector<Vec3> centred_positions(const Wavefunction& wf, const Vec3& centre) {
    std::vector<Vec3> positions;
    positions.reserve(wf.atoms.size());
    for (const auto& atom : wf.atoms) {
        positions.push_back({atom.x - centre[0], atom.y - centre[1], atom.z - centre[2]});
    }
    return positions;
}

bool validate_operation(const Wavefunction& wf,
                        const std::vector<Vec3>& positions,
                        const Mat3& matrix,
                        const double tolerance,
                        std::vector<std::size_t>& permutation,
                        double& max_error) {
    const std::size_t n = wf.atoms.size();
    permutation.assign(n, n);
    std::vector<bool> used(n, false);
    max_error = 0.0;

    const double tolerance2 = tolerance * tolerance;
    for (std::size_t i = 0; i < n; ++i) {
        const Vec3 transformed = mat_vec(matrix, positions[i]);
        std::size_t best = n;
        double best_d2 = std::numeric_limits<double>::infinity();
        for (std::size_t j = 0; j < n; ++j) {
            if (used[j] || wf.atoms[j].atomic_number != wf.atoms[i].atomic_number) continue;
            const Vec3 delta = sub(transformed, positions[j]);
            const double d2 = norm2(delta);
            if (d2 < best_d2) {
                best_d2 = d2;
                best = j;
            }
        }
        if (best == n || best_d2 > tolerance2 || !std::isfinite(best_d2)) return false;
        used[best] = true;
        permutation[i] = best;
        max_error = std::max(max_error, std::sqrt(best_d2));
    }
    return true;
}

void append_operation(MolecularSymmetry& result,
                      const SymmetryOperationKind kind,
                      const int order,
                      const int power,
                      const Vec3& axis,
                      const Mat3& matrix,
                      std::vector<std::size_t> permutation,
                      const double error) {
    SymmetryOperation op;
    op.kind = kind;
    op.order = order;
    op.power = power;
    op.axis_or_normal = axis;
    op.matrix = matrix;
    op.atom_permutation = std::move(permutation);
    op.max_mapping_error_bohr = error;
    result.operations.push_back(std::move(op));
    result.max_mapping_error_bohr = std::max(result.max_mapping_error_bohr, error);
}

bool is_linear_geometry(const std::vector<Vec3>& positions,
                        const double tolerance,
                        Vec3& axis) {
    std::size_t farthest = positions.size();
    double farthest_norm2 = 0.0;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const double r2 = norm2(positions[i]);
        if (r2 > farthest_norm2) {
            farthest_norm2 = r2;
            farthest = i;
        }
    }
    if (farthest == positions.size() || farthest_norm2 <= tolerance * tolerance) return true;
    axis = positions[farthest];
    normalize(axis);
    for (const auto& p : positions) {
        if (norm(cross(p, axis)) > tolerance) return false;
    }
    canonicalize_axis(axis);
    return true;
}

std::vector<Vec3> candidate_axes(const std::vector<Vec3>& positions,
                                 const SymmetryOptions& options) {
    std::vector<Vec3> axes;
    axes.reserve(std::min<std::size_t>(options.maximum_candidate_axes, 256));

    add_axis(axes, {1.0, 0.0, 0.0}, options.maximum_candidate_axes);
    add_axis(axes, {0.0, 1.0, 0.0}, options.maximum_candidate_axes);
    add_axis(axes, {0.0, 0.0, 1.0}, options.maximum_candidate_axes);

    for (const auto& p : positions) add_axis(axes, p, options.maximum_candidate_axes);

    for (std::size_t i = 0; i < positions.size() && axes.size() < options.maximum_candidate_axes; ++i) {
        for (std::size_t j = i + 1; j < positions.size() && axes.size() < options.maximum_candidate_axes; ++j) {
            add_axis(axes, cross(positions[i], positions[j]), options.maximum_candidate_axes);
            add_axis(axes, add(positions[i], positions[j]), options.maximum_candidate_axes);
            add_axis(axes, sub(positions[i], positions[j]), options.maximum_candidate_axes);
        }
    }
    return axes;
}

int count_axes_with_order_at_least(const std::vector<RotationAxis>& axes, const int order) {
    int count = 0;
    for (const auto& axis : axes) if (axis.order >= order) ++count;
    return count;
}

int maximum_rotation_order(const std::vector<RotationAxis>& axes) {
    int order = 1;
    for (const auto& axis : axes) order = std::max(order, axis.order);
    return order;
}

int maximum_improper_order(const std::vector<ImproperAxis>& axes) {
    int order = 1;
    for (const auto& axis : axes) order = std::max(order, axis.order);
    return order;
}

const RotationAxis* principal_rotation_axis(const std::vector<RotationAxis>& axes,
                                            const int order) {
    for (const auto& axis : axes) if (axis.order == order) return &axis;
    return nullptr;
}

bool has_plane_relation(const std::vector<ReflectionPlane>& planes,
                        const Vec3& axis,
                        const bool normal_parallel_axis) {
    for (const auto& plane : planes) {
        const double alignment = std::abs(dot(plane.normal, axis));
        if (normal_parallel_axis) {
            if (alignment > 1.0 - 2.0e-5) return true;
        } else {
            if (alignment < 2.0e-5) return true;
        }
    }
    return false;
}

int perpendicular_c2_count(const std::vector<RotationAxis>& axes,
                           const Vec3& principal) {
    int count = 0;
    for (const auto& axis : axes) {
        if (axis.order < 2) continue;
        if (std::abs(dot(axis.axis, principal)) < 2.0e-5) ++count;
    }
    return count;
}

} // namespace

MolecularSymmetry analyse_molecular_symmetry(const Wavefunction& wavefunction,
                                             const SymmetryOptions& options) {
    MolecularSymmetry result;
    if (wavefunction.atoms.empty()) return result;

    result.centre_bohr = molecular_centre(wavefunction);
    const auto positions = centred_positions(wavefunction, result.centre_bohr);
    for (const auto& p : positions) {
        result.molecular_radius_bohr = std::max(result.molecular_radius_bohr, norm(p));
    }
    result.tolerance_bohr = std::max(options.absolute_tolerance_bohr,
                                     options.relative_tolerance * std::max(1.0, result.molecular_radius_bohr));

    std::vector<std::size_t> permutation;
    double error = 0.0;
    if (validate_operation(wavefunction, positions, identity_matrix(), result.tolerance_bohr,
                           permutation, error)) {
        append_operation(result, SymmetryOperationKind::Identity, 1, 0,
                         {0.0, 0.0, 1.0}, identity_matrix(), std::move(permutation), error);
    }

    if (wavefunction.atoms.size() == 1u) {
        result.linear = true;
        result.point_group = "Kh";
        return result;
    }

    Vec3 linear_axis{0.0, 0.0, 1.0};
    result.linear = is_linear_geometry(positions, result.tolerance_bohr, linear_axis);

    permutation.clear();
    error = 0.0;
    if (validate_operation(wavefunction, positions, inversion_matrix(), result.tolerance_bohr,
                           permutation, error)) {
        result.has_inversion = true;
        append_operation(result, SymmetryOperationKind::Inversion, 2, 1,
                         {0.0, 0.0, 0.0}, inversion_matrix(), std::move(permutation), error);
    }

    if (result.linear) {
        result.point_group = result.has_inversion ? "Dinfh" : "Cinfv";
        return result;
    }

    const auto axes = candidate_axes(positions, options);
    std::vector<RotationAxis> rotation_axes;
    std::vector<ImproperAxis> improper_axes;
    std::vector<ReflectionPlane> planes;

    for (const Vec3& axis : axes) {
        int best_order = 1;
        double best_error = 0.0;
        std::vector<std::size_t> best_permutation;
        Mat3 best_matrix = identity_matrix();

        for (int order = std::max(2, options.maximum_rotation_order); order >= 2; --order) {
            const Mat3 matrix = rotation_matrix(axis, 2.0 * kPi / static_cast<double>(order));
            permutation.clear();
            error = 0.0;
            if (validate_operation(wavefunction, positions, matrix, result.tolerance_bohr,
                                   permutation, error)) {
                best_order = order;
                best_error = error;
                best_permutation = permutation;
                best_matrix = matrix;
                break;
            }
        }
        if (best_order >= 2) {
            bool duplicate = false;
            for (const auto& existing : rotation_axes) {
                if (same_axis(axis, existing.axis)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                rotation_axes.push_back({axis, best_order, best_error});
                append_operation(result, SymmetryOperationKind::ProperRotation,
                                 best_order, 1, axis, best_matrix,
                                 std::move(best_permutation), best_error);
            }
        }

        int best_improper = 1;
        double best_improper_error = 0.0;
        std::vector<std::size_t> best_improper_permutation;
        Mat3 best_improper_matrix = identity_matrix();
        for (int order = std::max(3, options.maximum_rotation_order); order >= 3; --order) {
            const Mat3 matrix = improper_matrix(axis, order);
            permutation.clear();
            error = 0.0;
            if (validate_operation(wavefunction, positions, matrix, result.tolerance_bohr,
                                   permutation, error)) {
                best_improper = order;
                best_improper_error = error;
                best_improper_permutation = permutation;
                best_improper_matrix = matrix;
                break;
            }
        }
        if (best_improper >= 3) {
            bool duplicate = false;
            for (const auto& existing : improper_axes) {
                if (same_axis(axis, existing.axis)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                improper_axes.push_back({axis, best_improper, best_improper_error});
                append_operation(result, SymmetryOperationKind::ImproperRotation,
                                 best_improper, 1, axis, best_improper_matrix,
                                 std::move(best_improper_permutation), best_improper_error);
            }
        }

        const Mat3 mirror = reflection_matrix(axis);
        permutation.clear();
        error = 0.0;
        if (validate_operation(wavefunction, positions, mirror, result.tolerance_bohr,
                               permutation, error)) {
            bool duplicate = false;
            for (const auto& existing : planes) {
                if (same_axis(axis, existing.normal)) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                planes.push_back({axis, error});
                append_operation(result, SymmetryOperationKind::Reflection,
                                 2, 1, axis, mirror, std::move(permutation), error);
            }
        }
    }

    const int max_order = maximum_rotation_order(rotation_axes);
    const int max_improper = maximum_improper_order(improper_axes);
    const int c5_axes = count_axes_with_order_at_least(rotation_axes, 5);
    const int c4_axes = count_axes_with_order_at_least(rotation_axes, 4);
    const int c3_axes = count_axes_with_order_at_least(rotation_axes, 3);
    const int c2_axes = count_axes_with_order_at_least(rotation_axes, 2);

    if (c5_axes >= 6 && c3_axes >= 10) {
        result.point_group = result.has_inversion ? "Ih" : "I";
    } else if (c4_axes >= 3) {
        result.point_group = result.has_inversion ? "Oh" : "O";
    } else if (c3_axes >= 4 && c2_axes >= 3) {
        if (result.has_inversion) result.point_group = "Th";
        else result.point_group = planes.size() >= 6u ? "Td" : "T";
    } else if (max_order >= 2) {
        const RotationAxis* principal = principal_rotation_axis(rotation_axes, max_order);
        if (principal) {
            const bool horizontal = has_plane_relation(planes, principal->axis, true);
            const bool vertical = has_plane_relation(planes, principal->axis, false);
            const int perpendicular_c2 = perpendicular_c2_count(rotation_axes, principal->axis);
            const bool dihedral = perpendicular_c2 >= max_order;
            if (dihedral) {
                if (horizontal) result.point_group = "D" + std::to_string(max_order) + "h";
                else if (vertical) result.point_group = "D" + std::to_string(max_order) + "d";
                else result.point_group = "D" + std::to_string(max_order);
            } else if (horizontal) {
                result.point_group = "C" + std::to_string(max_order) + "h";
            } else if (vertical) {
                result.point_group = "C" + std::to_string(max_order) + "v";
            } else if (max_improper > max_order && max_improper > 2) {
                result.point_group = "S" + std::to_string(max_improper);
            } else {
                result.point_group = "C" + std::to_string(max_order);
            }
        }
    }

    if (result.point_group.empty()) {
        if (result.has_inversion) result.point_group = "Ci";
        else if (!planes.empty()) result.point_group = "Cs";
        else result.point_group = "C1";
    }

    return result;
}

void derive_point_group_from_geometry(Wavefunction& wavefunction,
                                      const SymmetryOptions& options) {
    if (wavefunction.point_group_provenance == DataProvenance::Producer) return;
    const MolecularSymmetry symmetry = analyse_molecular_symmetry(wavefunction, options);
    if (!symmetry.available()) return;
    wavefunction.point_group_detected = symmetry.point_group;
    if (wavefunction.point_group_used.empty()) wavefunction.point_group_used = symmetry.point_group;
    wavefunction.point_group_provenance = DataProvenance::Derived;
}

} // namespace cov
