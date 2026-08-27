#pragma once

#include <filesystem>
#include <string>

namespace cov {

struct FileDialogResult {
    bool supported = true;
    bool cancelled = false;
    std::filesystem::path path;
    std::string error;

    [[nodiscard]] bool selected() const noexcept {
        return supported && !cancelled && error.empty() && !path.empty();
    }
};

[[nodiscard]] FileDialogResult open_molden_file_dialog();

} // namespace cov
