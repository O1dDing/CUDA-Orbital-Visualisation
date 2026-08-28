namespace {
void draw_selected_chemistry(const Wavefunction& wavefunction,
                             std::size_t selected_index,
                             const OrbitalUIState& state,
                             Language language,
                             float ui_scale);
}

void draw_energy_diagram(const Wavefunction& wavefunction,
                         const std::size_t selected_index,
                         OrbitalUIState& state,
                         const Language language,
                         const float ui_scale,
                         OrbitalUIActions& actions) {
    if (wavefunction.orbitals.empty()) {
        ImGui::TextDisabled("%s",tr(Text::NoOrbitals,language));
        return;
    }

    MODiagramOptions options;
    options.energy_unit=state.energy_unit;
    options.energy_axis_mode=state.energy_axis_mode;
    options.degeneracy=state.degeneracy;
    options.filter=state.filter;
    options.selected_index=selected_index;
    options.neighbourhood=static_cast<std::size_t>(
        std::max(3,state.diagram_neighbourhood));
    options.nonlinear_minimum_gap_weight=0.055;
    const MODiagramData data=build_mo_diagram_data(wavefunction,options);

    ImGui::TextDisabled("%s",tr(Text::EnergyDiagram,language));
    ImGui::TextDisabled("%s",data.selection.summary.c_str());
    ImGui::TextDisabled("%s",tr(Text::EnergyScale,language));
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::LinearEnergyScale,language),
            state.energy_axis_mode==EnergyAxisMode::Linear)) {
        state.energy_axis_mode=EnergyAxisMode::Linear;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::NonlinearFocus,language),
            state.energy_axis_mode==EnergyAxisMode::NonlinearFocus)) {
        state.energy_axis_mode=EnergyAxisMode::NonlinearFocus;
    }
    ImGui::SetNextItemWidth(155.0f*ui_scale);
    ImGui::SliderInt("##diagram_neighbourhood",
        &state.diagram_neighbourhood,3,32,"%d",
        ImGuiSliderFlags_AlwaysClamp);
    ImGui::SameLine();
    ImGui::TextDisabled("%s",tr(Text::AroundSelected,language));

    const float height=430.0f*ui_scale;
    ImGui::InvisibleButton("##energy_diagram_canvas",
        ImVec2(-1.0f,height),ImGuiButtonFlags_MouseButtonLeft);
    const ImVec2 p0=ImGui::GetItemRectMin();
    const ImVec2 p1=ImGui::GetItemRectMax();
    ImDrawList* draw=ImGui::GetWindowDrawList();
    draw->AddRectFilled(p0,p1,IM_COL32(12,18,27,235),7.0f*ui_scale);
    draw->AddRect(p0,p1,IM_COL32(43,58,77,220),7.0f*ui_scale);
    if (!data.levels.empty()) {
        const float top=p0.y+34.0f*ui_scale;
        const float bottom=p1.y-24.0f*ui_scale;
        const float left=p0.x+62.0f*ui_scale;
        const float right=p1.x-30.0f*ui_scale;
        const float lane_span=right-left;
        const float axis_x=left-28.0f*ui_scale;
        draw->AddLine(ImVec2(axis_x,bottom),ImVec2(axis_x,top),
            IM_COL32(129,148,171,255),1.4f*ui_scale);
        draw->AddTriangleFilled(
            ImVec2(axis_x,top-5.0f*ui_scale),
            ImVec2(axis_x-4.0f*ui_scale,top+3.0f*ui_scale),
            ImVec2(axis_x+4.0f*ui_scale,top+3.0f*ui_scale),
            IM_COL32(129,148,171,255));
        draw->AddText(ImVec2(p0.x+7.0f*ui_scale,p0.y+7.0f*ui_scale),
            IM_COL32(129,148,171,255),
            state.energy_axis_mode==EnergyAxisMode::Linear
                ?tr(Text::LinearEnergyScale,language)
                :tr(Text::NonlinearEnergyScale,language));

        for (const double fraction:{0.0,0.25,0.5,0.75,1.0}) {
            const double energy=energy_at_fraction(
                fraction,data.energy_transform);
            const float y=map_energy_y(
                energy,data.energy_transform,top,bottom);
            draw->AddLine(ImVec2(axis_x-4.0f*ui_scale,y),
                ImVec2(axis_x+4.0f*ui_scale,y),
                IM_COL32(100,118,140,255),1.0f);
            const std::string tick=format_energy(
                energy,state.energy_unit,3);
            draw->AddText(ImVec2(axis_x+6.0f*ui_scale,
                y-ImGui::GetTextLineHeight()*0.5f),
                IM_COL32(100,118,140,255),tick.c_str());
        }

        std::map<long long,std::vector<std::size_t>> coincident;
        for (std::size_t i=0;i<data.levels.size();++i) {
            coincident[static_cast<long long>(std::llround(
                data.levels[i].layout_energy_hartree*1.0e12))].push_back(i);
        }
        std::vector<float> x(data.levels.size(),left+0.5f*lane_span);
        for (const auto& [key,ids]:coincident) {
            (void)key;
            const float spacing=std::min(
                105.0f*ui_scale,
                lane_span/std::max(2.0f,static_cast<float>(ids.size())));
            for (std::size_t j=0;j<ids.size();++j) {
                x[ids[j]]+=(static_cast<float>(j)-
                    (static_cast<float>(ids.size())-1.0f)*0.5f)*spacing;
            }
        }

        std::vector<DiagramPoint> points;
        points.reserve(data.levels.size());
        for (std::size_t i=0;i<data.levels.size();++i) {
            const auto& level=data.levels[i];
            const float y=map_energy_y(level.layout_energy_hartree,
                data.energy_transform,top,bottom);
            const float half=std::min(35.0f*ui_scale,lane_span*0.18f);
            ImU32 colour=IM_COL32(196,210,227,255);
            if (level.metadata.region==OrbitalRegion::Virtual) {
                colour=IM_COL32(126,143,165,255);
            }
            if (level.annotation.delocalised_pi.available) {
                colour=ImGui::ColorConvertFloat4ToU32(kPiColour);
            } else if (level.annotation.family=="sigma") {
                colour=ImGui::ColorConvertFloat4ToU32(kSigmaColour);
            } else if (level.annotation.family=="delta") {
                colour=ImGui::ColorConvertFloat4ToU32(kDeltaColour);
            } else if (level.annotation.family=="phi") {
                colour=ImGui::ColorConvertFloat4ToU32(kPhiColour);
            }
            if (level.homo) colour=IM_COL32(79,210,157,255);
            if (level.lumo) colour=IM_COL32(228,176,82,255);
            if (level.metadata.selected) colour=IM_COL32(92,151,255,255);
            draw->AddLine(ImVec2(x[i]-half,y),ImVec2(x[i]+half,y),
                colour,level.metadata.selected?3.2f:1.8f);
            if (level.annotation.delocalised_pi.available) {
                draw->AddLine(ImVec2(x[i]-half,y+4.0f*ui_scale),
                    ImVec2(x[i]+half,y+4.0f*ui_scale),
                    ImGui::ColorConvertFloat4ToU32(kPiColour),
                    1.4f*ui_scale);
            }
            if (level.electrons.alpha>0) {
                draw_arrow(draw,ImVec2(x[i]-7.0f*ui_scale,y-2.0f),
                    true,IM_COL32(234,242,252,255));
            }
            if (level.electrons.beta>0) {
                draw_arrow(draw,ImVec2(x[i]+7.0f*ui_scale,y-2.0f),
                    false,IM_COL32(234,242,252,255));
            }
            if (level.metadata.selected) {
                draw->AddText(ImVec2(x[i]+half+4.0f*ui_scale,
                    y-ImGui::GetTextLineHeight()*0.5f),
                    colour,level.metadata.display_label.c_str());
            }
            points.push_back({&level,ImVec2(x[i]-half,y),
                              ImVec2(x[i]+half,y),y});
        }

        for (const auto& group:pi_diagram_groups(data)) {
            if (!group.chemistry || group.levels.size()<3u) continue;
            float y_min=bottom;
            float y_max=top;
            float x_max=left;
            for (const auto index:group.levels) {
                if (index>=points.size()) continue;
                y_min=std::min(y_min,points[index].y);
                y_max=std::max(y_max,points[index].y);
                x_max=std::max(x_max,points[index].right.x);
            }
            const float bracket_x=std::min(
                right-54.0f*ui_scale,x_max+18.0f*ui_scale);
            const ImU32 colour=ImGui::ColorConvertFloat4ToU32(kPiColour);
            draw->AddLine(ImVec2(bracket_x,y_min-7.0f*ui_scale),
                ImVec2(bracket_x,y_max+7.0f*ui_scale),
                colour,2.0f*ui_scale);
            for (const auto index:group.levels) {
                if (index>=points.size()) continue;
                draw->AddLine(ImVec2(points[index].right.x,
                                     points[index].y),
                    ImVec2(bracket_x,points[index].y),
                    colour,1.2f*ui_scale);
            }
            const std::string label=pi_descriptor(*group.chemistry);
            draw->AddText(ImVec2(bracket_x+5.0f*ui_scale,
                0.5f*(y_min+y_max)-ImGui::GetTextLineHeight()*0.5f),
                colour,label.c_str());
        }

        if (ImGui::IsItemHovered()) {
            const ImVec2 mouse=ImGui::GetIO().MousePos;
            const DiagramPoint* nearest=nullptr;
            float best=1.0e9f;
            for (const auto& point:points) {
                if (mouse.x<point.left.x-10.0f*ui_scale ||
                    mouse.x>point.right.x+10.0f*ui_scale) continue;
                const float distance=std::abs(mouse.y-point.y);
                if (distance<best && distance<=8.0f*ui_scale) {
                    best=distance;
                    nearest=&point;
                }
            }
            if (nearest && nearest->level) {
                draw_level_tooltip(data,*nearest->level,state,language);
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    actions.select_orbital=
                        nearest->level->metadata.orbital_index;
                }
            }
        }
    }

    if (ImGui::Button(tr(Text::ExportBundle,language),
                      ImVec2(-1.0f,0.0f))) {
        actions.export_diagram=true;
    }
    draw_selected_chemistry(
        wavefunction,selected_index,state,language,ui_scale);
}

} // namespace cov::ui
