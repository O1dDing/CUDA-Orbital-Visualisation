#include "cov/formchk.hpp"

#include "cov/fchk_overlap.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

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

#ifdef _WIN32
    // Invoke the converter directly: no cmd.exe, no console window, and no
    // shell expansion of user-selected paths. Quote using Windows argv rules.
    const auto quote = [](const std::wstring& value) {
        std::wstring result=L"\"";
        std::size_t slashes=0;
        for(wchar_t ch:value) {
            if(ch==L'\\') { ++slashes; continue; }
            result.append(ch==L'"'?slashes*2+1:slashes,L'\\');
            result+=ch;slashes=0;
        }
        result.append(slashes*2,L'\\');result+=L'"';return result;
    };
    std::wstring command=quote(std::filesystem::path(executable).wstring())+L" "+
        quote(chk_path.wstring())+L" "+quote(output.path.wstring());
    STARTUPINFOW startup{};startup.cb=sizeof(startup);
    PROCESS_INFORMATION process{};
    if(!CreateProcessW(nullptr,command.data(),nullptr,nullptr,FALSE,
                       CREATE_NO_WINDOW,nullptr,nullptr,&startup,&process))
        throw std::runtime_error("Cannot start Gaussian formchk (Windows error "+
            std::to_string(GetLastError())+"). Set COV_FORMCHK to its executable path.");
    CloseHandle(process.hThread);
    const DWORD wait=WaitForSingleObject(process.hProcess,INFINITE);
    DWORD exit_code=1;
    const bool obtained=wait==WAIT_OBJECT_0 && GetExitCodeProcess(process.hProcess,&exit_code);
    CloseHandle(process.hProcess);
    if(!obtained)throw std::runtime_error("Cannot obtain Gaussian formchk exit status");
    const auto code=exit_code;
#else
    const std::string command =
        quoted_shell_argument(executable, "formchk executable") + " " +
        quoted_shell_argument(chk_path.string(), "CHK path") + " " +
        quoted_shell_argument(output.path.string(), "temporary FCHK path");

    const int code = std::system(command.c_str());
#endif
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
