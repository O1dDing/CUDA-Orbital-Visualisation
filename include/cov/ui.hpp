#pragma once

namespace cov::ui {

enum class Language {
    English = 0,
    ChineseSimplified,
    Japanese,
    French,
    Count,
};

enum class Text {
    AppTitle = 0,
    Tagline,
    LanguageLabel,
    FileSection,
    MoldenPath,
    Load,
    IdleHint,
    WavefunctionSection,
    Atoms,
    Shells,
    BasisFunctions,
    Orbitals,
    ShellConvention,
    OrbitalSection,
    MoldenMO,
    InternalIndex,
    Energy,
    Occupation,
    Spin,
    Symmetry,
    RenderingSection,
    Isovalue,
    Grid,
    RecomputeGrid,
    ResetCamera,
    PerformanceSection,
    CUDADevice,
    LastKernel,
    GPUResident,
    InteractionHint,
    IsovalueHint,
    ExperimentalNote,
    Ready,
    Parsing,
    Loaded,
    GridUpdated,
    Error,
    NoneValue,
    Alpha,
    Beta,
    FontStatus,
    Count,
};

enum class Tone {
    Neutral,
    Accent,
    Success,
    Danger,
};

[[nodiscard]] const char* tr(Text key, Language language) noexcept;
[[nodiscard]] const char* language_name(Language language) noexcept;

void apply_theme(float scale = 1.0f);
bool configure_fonts(float pixel_size = 17.0f);
[[nodiscard]] const char* font_status() noexcept;

void section_title(const char* label);
void begin_card(const char* id, float height);
void end_card();
void status_badge(const char* label, Tone tone);

} // namespace cov::ui
