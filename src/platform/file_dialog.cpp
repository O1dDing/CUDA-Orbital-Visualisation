#include "cov/file_dialog.hpp"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#include <commdlg.h>

#include <array>
#include <sstream>

namespace cov {

FileDialogResult open_wavefunction_file_dialog() {
    FileDialogResult result;

    std::array<wchar_t, 32768> buffer{};
    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = nullptr;
    ofn.lpstrFile = buffer.data();
    ofn.nMaxFile = static_cast<DWORD>(buffer.size());
    ofn.lpstrFilter =
        L"Gaussian FCHK (*.fchk;*.fch)\0"
        L"*.fchk;*.fch\0"
        L"Molden wavefunction (*.molden;*.molden.input;*.molden.inp)\0"
        L"*.molden;*.molden.input;*.molden.inp\0"
        L"Supported wavefunctions (*.fchk;*.fch;*.molden;*.molden.input;*.molden.inp)\0"
        L"*.fchk;*.fch;*.molden;*.molden.input;*.molden.inp\0"
        L"All files (*.*)\0"
        L"*.*\0\0";
    ofn.nFilterIndex = 1;
    ofn.lpstrTitle = L"Open wavefunction";
    ofn.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST |
                OFN_NOCHANGEDIR | OFN_DONTADDTORECENT;

    if (GetOpenFileNameW(&ofn) != FALSE) {
        result.path = std::filesystem::path(buffer.data());
        return result;
    }

    const DWORD error = CommDlgExtendedError();
    if (error == 0) {
        result.cancelled = true;
        return result;
    }

    std::ostringstream message;
    message << "Windows Open File dialog failed (CommDlgExtendedError=" << error << ')';
    result.error = message.str();
    return result;
}

} // namespace cov

#else

namespace cov {

FileDialogResult open_wavefunction_file_dialog() {
    FileDialogResult result;
    result.supported = false;
    result.error = "Native Open File dialog is currently implemented for Windows only";
    return result;
}

} // namespace cov

#endif
