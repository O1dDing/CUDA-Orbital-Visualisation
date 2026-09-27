#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdlib>

namespace cov {

// A batch supervisor assigns this budget once per process. Do not infer a
// fresh entitlement from the host's hardware count inside each analysis.
inline std::size_t analysis_thread_budget() noexcept {
    const char* value=std::getenv("COV_CPU_THREADS");
    if (!value || *value<'0' || *value>'9') return 1;
    char* end=nullptr;
    const auto parsed=std::strtoul(value,&end,10);
    if (end==value || *end!='\0' || parsed==0) return 1;
    return std::clamp<std::size_t>(parsed,1,12);
}

} // namespace cov
