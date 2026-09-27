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

[[nodiscard]] FileDialogResult open_wavefunction_file_dialog();

// Compatibility alias for the PR #2 viewer call-site. FCHK is now included in
// the dialog and the function may be removed once all callers use the generic name.
[[nodiscard]] inline FileDialogResult open_molden_file_dialog() {
    return open_wavefunction_file_dialog();
}

} // namespace cov
