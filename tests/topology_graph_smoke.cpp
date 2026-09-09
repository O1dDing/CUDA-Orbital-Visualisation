#include "cov/topology_graph.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>

int main() {
    using Edge=std::pair<std::uint32_t,std::uint32_t>;
    const std::vector<Edge> rings{{0,1},{1,2},{2,0},{0,3},{3,4},{4,0}};
    const std::array<std::vector<std::uint32_t>,2> channels{{{1,2},{3,4}}};
    auto check=[&](std::uint32_t n,const std::vector<Edge>& edges,
                   const std::array<std::vector<std::uint32_t>,2>& assigned,
                   bool expected,bool additional=false,std::uint32_t hub=0) {
        const auto result=cov::find_channel_cycle_union_at_hub(n,edges,hub,assigned);
        if(result.witness.has_value()!=expected)return false;
        if(!expected)return true;
        const auto& witness=*result.witness;
        if(!cov::verify_cycle_union(n,edges,witness) || witness.hub!=hub ||
           witness.additional_connection_without_hub!=additional)return false;
        for(std::size_t ring=0;ring<2;++ring) {
            std::size_t matching=0;
            for(const auto atom:witness.cycles[ring]) {
                if(atom==hub)continue;
                if(std::find(assigned[1-ring].begin(),assigned[1-ring].end(),atom)!=assigned[1-ring].end())return false;
                matching+=std::find(assigned[ring].begin(),assigned[ring].end(),atom)!=assigned[ring].end();
            }
            if(matching<2)return false;
        }
        return true;
    };
    if(!check(5,rings,channels,true))return 1;
    auto extra=rings;extra.push_back({2,3});if(!check(5,extra,channels,true,true))return 2;
    extra=rings;extra.push_back({2,5});extra.push_back({5,3});if(!check(6,extra,channels,true,true))return 3;
    if(!check(5,{{0,1},{0,2},{0,3},{0,4}},channels,false))return 4;
    if(!check(5,rings,{{{1,2,3,4},{1,2,3,4}}},false))return 5;
    if(!check(5,rings,{{{1,3},{2,4}}},false))return 6;
    const std::vector<Edge> extended{{0,1},{1,5},{5,2},{2,0},{0,3},{3,6},{6,4},{4,0}};
    if(!check(7,extended,channels,true))return 7;
    if(!check(7,extended,{{{0,1,2},{0,3,4}}},true))return 8;
    extra=rings;extra.push_back({5,6});extra.push_back({6,7});extra.push_back({7,8});
    if(!check(9,extra,{{{5,6},{7,8}}},false))return 9;
    const auto ambiguous=cov::find_channel_cycle_union_at_hub(5,rings,0,{{{1,2,3,4},{1,2,3,4}}});
    if(!ambiguous.unscoped_cycle_union_available || ambiguous.witness)return 10;
    std::vector<std::uint32_t> order{0,1,2,3,4};
    do {
        std::vector<Edge> edges;
        for(const auto& [a,b]:rings)edges.push_back({order[b],order[a]});
        const std::array<std::vector<std::uint32_t>,2> reordered{{{order[1],order[2]},{order[3],order[4]}}};
        if(!check(5,edges,reordered,true,false,order[0]))return 11;
    } while(std::next_permutation(order.begin(),order.end()));
    try {cov::find_channel_cycle_union_at_hub(5,rings,0,{{{1,9},{3,4}}});return 12;}
    catch(const std::invalid_argument&) {}
    std::cout<<"Explicit structural ring paths and channel scope controls passed\n";
}
