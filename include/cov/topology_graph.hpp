#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace cov {
// A graph witness, not a bond detector or a complete chemical classification.
// Each closed simple cycle contains hub and they share no other vertex.
struct CycleUnionWitness {
    std::uint32_t hub=0;
    std::array<std::vector<std::uint32_t>,2> cycles;
    bool additional_connection_without_hub=false;
};
std::optional<CycleUnionWitness> find_cycle_union_at_hub(
    std::uint32_t vertex_count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,
    std::uint32_t hub);
bool verify_cycle_union(std::uint32_t vertex_count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,
    const CycleUnionWitness& witness);

// A stronger, explicitly scoped association with two electronic channels.
// Each ring has at least two atoms exclusive to its channel, has no atom of
// the other channel except hub, and may include intervening non-channel atoms.
// Failure to associate is not proof that the molecule has no ring union.
struct ChannelCycleSearchResult {
    std::optional<CycleUnionWitness> witness;
    bool unscoped_cycle_union_available=false;
    std::array<std::size_t,2> exclusive_channel_atoms{};
    std::uint64_t qualifying_first_cycles_examined=0;
};
ChannelCycleSearchResult find_channel_cycle_union_at_hub(
    std::uint32_t vertex_count,
    const std::vector<std::pair<std::uint32_t,std::uint32_t>>& edges,
    std::uint32_t hub,
    const std::array<std::vector<std::uint32_t>,2>& channels);
}
