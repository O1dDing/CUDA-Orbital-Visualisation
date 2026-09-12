#include "cov/mo_diagram.hpp"

namespace cov {

const char* annotation_source_name_legacy(AnnotationSource source) noexcept;

const char* annotation_source_name(const AnnotationSource source) noexcept {
    if (source==AnnotationSource::Derived) return "derived";
    return annotation_source_name_legacy(source);
}

} // namespace cov
