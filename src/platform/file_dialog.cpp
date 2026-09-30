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

FileDialogResult open_wavefunction_file_dialog(ui::Language language, bool nbo_input) {
    FileDialogResult result;

    std::array<wchar_t, 32768> buffer{};
    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = nullptr;
    ofn.lpstrFile = buffer.data();
    ofn.nMaxFile = static_cast<DWORD>(buffer.size());
    const auto local = [&](const wchar_t* en, const wchar_t* zh,
                           const wchar_t* ja, const wchar_t* fr) {
        switch(language) {
            case ui::Language::ChineseSimplified: return zh;
            case ui::Language::Japanese: return ja;
            case ui::Language::French: return fr;
            default: return en;
        }
    };
    std::wstring filters;
    const auto add_filter = [&](const wchar_t* label, const wchar_t* pattern) {
        filters += label; filters.push_back(L'\0');
        filters += pattern; filters.push_back(L'\0');
    };
    if(nbo_input) {
        add_filter(local(L"NBO data", L"NBO 数据", L"NBO データ", L"Données NBO"),
                   L"*.covnbopkg;*.log;*.out;*.nbo;*.47");
    } else {
        add_filter(local(L"Calculation files", L"计算文件", L"計算ファイル", L"Fichiers de calcul"),
                   L"*.fchk;*.fch;*.chk;*.molden;*.molden.input;*.molden.inp;*.covnbopkg");
    }
    add_filter(L"Gaussian FCHK / FCH", L"*.fchk;*.fch");
    add_filter(L"Gaussian CHK (formchk)", L"*.chk");
    add_filter(L"Molden", L"*.molden;*.molden.input;*.molden.inp");
    add_filter(local(L"All files", L"所有文件", L"すべてのファイル", L"Tous les fichiers"), L"*.*");
    filters.push_back(L'\0');
    ofn.lpstrFilter = filters.c_str();
    ofn.nFilterIndex = 1;
    ofn.lpstrTitle = nbo_input
        ? local(L"Open NBO data", L"打开 NBO 数据", L"NBO データを開く", L"Ouvrir des données NBO")
        : local(L"Open calculation", L"打开计算文件", L"計算ファイルを開く", L"Ouvrir un calcul");
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
    message << error;
    const char* detail = language == ui::Language::ChineseSimplified ? "无法打开文件窗口。错误码：" :
        language == ui::Language::Japanese ? "ファイル選択を開けません。エラー：" :
        language == ui::Language::French ? "La fenêtre de sélection ne s’ouvre pas. Code : " :
        "The file picker could not open. Error: ";
    result.error = detail + message.str();
    return result;
}

} // namespace cov

#else

namespace cov {

FileDialogResult open_wavefunction_file_dialog(ui::Language, bool) {
    FileDialogResult result;
    result.supported = false;
    result.error = "Native Open File dialog is currently implemented for Windows only";
    return result;
}

} // namespace cov

#endif
