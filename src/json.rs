use crate::model::{AppState, Chart, SingularAnalysis, SingularityType};

pub fn state_to_json(state: &AppState) -> String {
    let analysis = state.singular_locus();
    state_to_json_with_analysis(state, &analysis, true)
}

pub fn state_to_json_without_analysis(state: &AppState) -> String {
    state_to_json_with_analysis(
        state,
        &SingularAnalysis {
            components: Vec::new(),
            intersections: Vec::new(),
        },
        false,
    )
}

fn state_to_json_with_analysis(
    state: &AppState,
    analysis: &SingularAnalysis,
    include_singular_reports: bool,
) -> String {
    let mut out = String::new();
    out.push_str("{\n");
    push_field(&mut out, 1, "schema", "\"blowup-script/v1\"", true);
    push_field(&mut out, 1, "base_field", "\"Q\"", true);
    push_field(
        &mut out,
        1,
        "ambient",
        &format!(
            "{{\"kind\":\"projective\",\"hypersurface_dimension\":{},\"projective_space_dimension\":{},\"coordinates\":{}}}",
            state.projective_dim,
            state.projective_dim + 1,
            json_string_array(&state.projective_vars)
        ),
        true,
    );
    push_field(
        &mut out,
        1,
        "initial_hypersurface",
        &format!(
            "{{\"polynomial\":\"{}\",\"normalized_polynomial\":\"{}\",\"degree\":{}}}",
            escape_json(&state.initial_polynomial_text),
            escape_json(&state.initial_polynomial.to_string()),
            state.degree
        ),
        true,
    );

    out.push_str("  \"charts\": [\n");
    for (index, chart) in state.charts.iter().enumerate() {
        out.push_str(&chart_to_json(chart, 2));
        if index + 1 != state.charts.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ],\n");

    out.push_str("  \"current_singular_components\": [\n");
    for (component_index, component) in analysis.components.iter().enumerate() {
        out.push_str("    {\n");
        push_field(
            &mut out,
            3,
            "index",
            &(component_index + 1).to_string(),
            true,
        );
        push_field(
            &mut out,
            3,
            "dimension",
            &component.dimension.to_string(),
            true,
        );
        push_field(
            &mut out,
            3,
            "generic_multiplicity",
            &component.multiplicity.to_string(),
            true,
        );
        push_field(
            &mut out,
            3,
            "generic_singularity_type",
            &optional_singularity_type_to_json(component.generic_singularity_type.as_ref()),
            true,
        );
        push_field(
            &mut out,
            3,
            "generic_singularity_chart",
            &component
                .generic_singularity_chart
                .map(|chart_id| chart_id.to_string())
                .unwrap_or_else(|| "null".to_string()),
            true,
        );
        out.push_str("      \"generic_singularity_point\": ");
        if let Some(assignments) = &component.generic_singularity_point {
            out.push('{');
            for (assignment_index, (name, value)) in assignments.iter().enumerate() {
                if assignment_index > 0 {
                    out.push(',');
                }
                out.push_str(&format!("\"{}\":\"{}\"", escape_json(name), value));
            }
            out.push_str("},\n");
        } else {
            out.push_str("null,\n");
        }
        out.push_str("      \"charts\": [\n");
        for (chart_index, chart) in component.charts.iter().enumerate() {
            out.push_str("        {\n");
            push_field(&mut out, 5, "chart", &chart.chart_id.to_string(), true);
            push_field(
                &mut out,
                5,
                "chart_label",
                &format!("\"{}\"", escape_json(&chart.chart_label)),
                true,
            );
            push_field(
                &mut out,
                5,
                "variables",
                &json_string_array(&chart.variables),
                true,
            );
            push_field(&mut out, 5, "dimension", &chart.dimension.to_string(), true);
            push_field(
                &mut out,
                5,
                "generic_multiplicity",
                &chart.multiplicity.to_string(),
                true,
            );
            push_field(
                &mut out,
                5,
                "singularity_type",
                &optional_singularity_type_to_json(chart.singularity_type.as_ref()),
                true,
            );
            push_field(
                &mut out,
                5,
                "affine_linear_center",
                if chart.affine_linear_center {
                    "true"
                } else {
                    "false"
                },
                true,
            );
            out.push_str("          \"prime_ideal_generators\": [");
            for (generator_index, generator) in chart.ideal_generators.iter().enumerate() {
                if generator_index > 0 {
                    out.push(',');
                }
                out.push_str(&format!(
                    "{{\"name\":\"f{}\",\"polynomial\":\"{}\"}}",
                    generator_index,
                    escape_json(&generator.to_string())
                ));
            }
            out.push_str("],\n");
            out.push_str("          \"coordinate_center\": ");
            if let Some(assignments) = &chart.coordinate_assignments {
                out.push('{');
                for (assignment_index, (name, value)) in assignments.iter().enumerate() {
                    if assignment_index > 0 {
                        out.push(',');
                    }
                    out.push_str(&format!("\"{}\":\"{}\"", escape_json(name), value));
                }
                out.push_str("}\n");
            } else {
                out.push_str("null\n");
            }
            out.push_str("        }");
            if chart_index + 1 != component.charts.len() {
                out.push(',');
            }
            out.push('\n');
        }
        out.push_str("      ]\n");
        out.push_str("    }");
        if component_index + 1 != analysis.components.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ],\n");

    out.push_str("  \"current_singular_intersections\": [\n");
    for (intersection_index, intersection) in analysis.intersections.iter().enumerate() {
        out.push_str("    {\n");
        push_field(
            &mut out,
            3,
            "index",
            &(intersection_index + 1).to_string(),
            true,
        );
        push_field(
            &mut out,
            3,
            "components",
            &json_usize_array(&intersection.component_indices),
            true,
        );
        push_field(
            &mut out,
            3,
            "chart",
            &intersection.chart_id.to_string(),
            true,
        );
        push_field(
            &mut out,
            3,
            "chart_label",
            &format!("\"{}\"", escape_json(&intersection.chart_label)),
            true,
        );
        push_field(
            &mut out,
            3,
            "variables",
            &json_string_array(&intersection.variables),
            true,
        );
        out.push_str("      \"ideal_generators\": [");
        for (generator_index, generator) in intersection.ideal_generators.iter().enumerate() {
            if generator_index > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"name\":\"f{}\",\"polynomial\":\"{}\"}}",
                generator_index,
                escape_json(&generator.to_string())
            ));
        }
        out.push_str("],\n");
        out.push_str("      \"coordinates\": ");
        if let Some(assignments) = &intersection.coordinate_assignments {
            out.push('{');
            for (assignment_index, (name, value)) in assignments.iter().enumerate() {
                if assignment_index > 0 {
                    out.push(',');
                }
                out.push_str(&format!("\"{}\":\"{}\"", escape_json(name), value));
            }
            out.push_str("},\n");
        } else {
            out.push_str("null,\n");
        }
        push_field(
            &mut out,
            3,
            "singularity_type",
            &optional_singularity_type_to_json(intersection.singularity_type.as_ref()),
            false,
        );
        out.push_str("    }");
        if intersection_index + 1 != analysis.intersections.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ],\n");

    out.push_str("  \"current_singular_locus\": [\n");
    let reports = if include_singular_reports {
        state.singular_reports()
    } else {
        Vec::new()
    };
    for (index, report) in reports.iter().enumerate() {
        out.push_str("    {\n");
        push_field(&mut out, 3, "chart", &report.chart_id.to_string(), true);
        push_field(
            &mut out,
            3,
            "chart_label",
            &format!("\"{}\"", escape_json(&report.chart_label)),
            true,
        );
        out.push_str("      \"equations\": [");
        for (equation_index, equation) in report.equations.iter().enumerate() {
            if equation_index > 0 {
                out.push(',');
            }
            let name = if equation_index == 0 {
                "f".to_string()
            } else {
                format!("d/d{}", report.variables[equation_index - 1])
            };
            out.push_str(&format!(
                "{{\"name\":\"{}\",\"polynomial\":\"{}\"}}",
                escape_json(&name),
                escape_json(&equation.to_string())
            ));
        }
        out.push_str("],\n");
        out.push_str("      \"rational_point_samples\": [");
        for (point_index, point) in report.point_options.iter().enumerate() {
            if point_index > 0 {
                out.push(',');
            }
            out.push_str("{\"coordinates\":{");
            for (coord_index, (name, value)) in
                report.variables.iter().zip(point.values.iter()).enumerate()
            {
                if coord_index > 0 {
                    out.push(',');
                }
                out.push_str(&format!("\"{}\":\"{}\"", escape_json(name), value));
            }
            out.push_str("}}");
        }
        out.push_str("],\n");
        push_field(
            &mut out,
            3,
            "sample_limit_hit",
            if report.sample_limit_hit {
                "true"
            } else {
                "false"
            },
            false,
        );
        out.push_str("    }");
        if index + 1 != reports.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ],\n");

    out.push_str("  \"blowups\": [\n");
    for (index, blowup) in state.blowups.iter().enumerate() {
        out.push_str("    {\n");
        push_field(&mut out, 3, "stage", &blowup.stage.to_string(), true);
        push_field(
            &mut out,
            3,
            "input_chart",
            &blowup.input_chart.to_string(),
            true,
        );
        push_field(
            &mut out,
            3,
            "multiplicity",
            &blowup.multiplicity.to_string(),
            true,
        );
        out.push_str("      \"center\": {\"type\":\"affine_coordinate_center\",\"assignments\":{");
        for (assignment_index, (name, value)) in blowup.center.assignments.iter().enumerate() {
            if assignment_index > 0 {
                out.push(',');
            }
            out.push_str(&format!("\"{}\":\"{}\"", escape_json(name), value));
        }
        out.push_str("}},\n");
        push_field(
            &mut out,
            3,
            "output_charts",
            &json_usize_array(&blowup.output_charts),
            false,
        );
        out.push_str("    }");
        if index + 1 != state.blowups.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ]\n");
    out.push_str("}\n");
    out
}

fn chart_to_json(chart: &Chart, indent: usize) -> String {
    let pad = "  ".repeat(indent);
    let mut out = String::new();
    out.push_str(&format!("{pad}{{\n"));
    push_field(&mut out, indent + 1, "id", &chart.id.to_string(), true);
    if let Some(parent) = chart.parent {
        push_field(&mut out, indent + 1, "parent", &parent.to_string(), true);
    } else {
        push_field(&mut out, indent + 1, "parent", "null", true);
    }
    push_field(
        &mut out,
        indent + 1,
        "label",
        &format!("\"{}\"", escape_json(&chart.label)),
        true,
    );
    push_field(
        &mut out,
        indent + 1,
        "active",
        if chart.active { "true" } else { "false" },
        true,
    );
    push_field(
        &mut out,
        indent + 1,
        "variables",
        &json_string_array(&chart.variables),
        true,
    );
    push_field(
        &mut out,
        indent + 1,
        "polynomial",
        &format!("\"{}\"", escape_json(&chart.polynomial.to_string())),
        true,
    );
    push_field(
        &mut out,
        indent + 1,
        "substitutions",
        &json_string_array(&chart.substitutions),
        false,
    );
    out.push_str(&format!("{pad}}}"));
    out
}

fn optional_singularity_type_to_json(singularity_type: Option<&SingularityType>) -> String {
    singularity_type
        .map(singularity_type_to_json)
        .unwrap_or_else(|| "null".to_string())
}

fn singularity_type_to_json(singularity_type: &SingularityType) -> String {
    format!(
        "{{\"label\":\"{}\",\"multiplicity\":{},\"is_singular\":{},\"embedding_dimension\":{},\"tangent_cone_degree\":{},\"tangent_cone\":\"{}\",\"quadratic_rank\":{}}}",
        escape_json(&singularity_type.label),
        singularity_type.multiplicity,
        if singularity_type.is_singular {
            "true"
        } else {
            "false"
        },
        singularity_type.embedding_dimension,
        singularity_type
            .tangent_cone_degree
            .map(|degree| degree.to_string())
            .unwrap_or_else(|| "null".to_string()),
        escape_json(&singularity_type.tangent_cone.to_string()),
        singularity_type
            .quadratic_rank
            .map(|rank| rank.to_string())
            .unwrap_or_else(|| "null".to_string())
    )
}

fn push_field(out: &mut String, indent: usize, name: &str, value: &str, trailing_comma: bool) {
    out.push_str(&"  ".repeat(indent));
    out.push_str(&format!("\"{name}\": {value}"));
    if trailing_comma {
        out.push(',');
    }
    out.push('\n');
}

fn json_string_array(values: &[String]) -> String {
    let items = values
        .iter()
        .map(|value| format!("\"{}\"", escape_json(value)))
        .collect::<Vec<_>>();
    format!("[{}]", items.join(","))
}

fn json_usize_array(values: &[usize]) -> String {
    let items = values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>();
    format!("[{}]", items.join(","))
}

pub fn escape_json(value: &str) -> String {
    let mut escaped = String::new();
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            ch if ch.is_control() => escaped.push_str(&format!("\\u{:04x}", ch as u32)),
            ch => escaped.push(ch),
        }
    }
    escaped
}
