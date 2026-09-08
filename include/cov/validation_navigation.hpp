#pragma once
#include <imgui.h>

struct ImGuiWindow;
namespace cov::validation {
struct NavigationTarget {
    ImVec2 lo, hi;
    ImGuiWindow* window = nullptr;
    ImVec2 clip_lo, clip_hi;
};
enum class NavigationKind { Ready, Move, Wheel, Unreachable };
struct NavigationStep {
    NavigationKind kind = NavigationKind::Unreachable;
    ImVec2 mouse{}, wheel{};
    ImGuiWindow* scrolling_window = nullptr;
};
// Observe current ImGui geometry; propose only ordinary pointer/wheel input.
// A pointer move occupies its own frame, so ImGui's input trickling cannot
// leave our wheel queued for the native driver's next ClearEventsQueue().
[[nodiscard]] NavigationStep plan_navigation(const NavigationTarget& target);
[[nodiscard]] bool navigation_target_visible(const NavigationTarget& target);
}
