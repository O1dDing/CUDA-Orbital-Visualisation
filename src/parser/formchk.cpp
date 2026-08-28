#include "cov/formchk.hpp"

#include "cov/fchk_overlap.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace cov {
namespace {

std::string quoted_shell_argument(const std::string& value, const char* label) {
    if (value.find('"') != std::string::npos ||
        value.find('\r') != std::string::npos ||
        value.find('\n') != std::string::npos) {
        throw std::runtime_error(std::string(label) +
                                 " contains characters unsafe for formchk invocation");
    }
    return '"' + value + '"';
}

struct TemporaryFileGuard {
    std::filesystem::path path;
    ~TemporaryFileGuard() {
        std::error_code ec;
        if (!path.empty()) std::filesystem::remove(path, ec);
    }
};

} // namespace

Wavefunction parse_gaussian_chk_via_formchk(const std::filesystem::path& chk_path,
                                            const FchkParseOptions& options) {
    std::error_code ec;
    if (!std::filesystem::exists(chk_path, ec) || ec) {
        throw std::runtime_error("Gaussian CHK file does not exist: " + chk_path.string());
    }

    std::string executable;
    if (const char* configured = std::getenv("COV_FORMCHK"); configured && *configured) {
        executable = configured;
    } else {
#ifdef _WIN32
        executable = "formchk.exe";
#else
        executable = "formchk";
#endif
    }

    const auto stamp = std::chrono::high_resolution_clock::now()
                           .time_since_epoch().count();
    TemporaryFileGuard output{
        std::filesystem::temp_directory_path() /
        ("cov_formchk_" + std::to_string(stamp) + ".fchk")
    };

    const std::string command =
        quoted_shell_argument(executable, "formchk executable") + " " +
        quoted_shell_argument(chk_path.string(), "CHK path") + " " +
        quoted_shell_argument(output.path.string(), "temporary FCHK path");

    const int code = std::system(command.c_str());
    if (code != 0) {
        throw std::runtime_error(
            "Gaussian formchk failed with exit code " + std::to_string(code) +
            ". Install Gaussian formchk or set COV_FORMCHK to its executable path.");
    }

    if (!std::filesystem::exists(output.path, ec) || ec ||
        std::filesystem::file_size(output.path, ec) == 0 || ec) {
        throw std::runtime_error(
            "Gaussian formchk reported success but did not produce a usable FCHK file");
    }

    Wavefunction wf=parse_fchk(output.path, options);
    (void)enrich_fchk_overlap_from_file(wf, output.path);
    return wf;
}

} // namespace cov
