use crate::json::{state_to_json, state_to_json_without_analysis};
use crate::model::{
    AppState, AutoResolutionCenter, AutoResolutionOptions, AutoResolutionProgress,
    AutoResolutionResult, AutoResolutionStatus, Chart, SingularComponent, SingularComponentChart,
    SingularIntersection, SingularityType,
};
use crate::poly::Poly;
use crate::rational::Rational;
use std::collections::BTreeMap;
use std::io::{self, IsTerminal, Read, Write};
use std::path::Path;
use std::process::{Command, Stdio};

const COMPONENT_POLY_DETAIL_CHART_LIMIT: usize = 3;
const POLYNOMIAL_DISPLAY_LIMIT: usize = 100;

struct LineEditor {
    history: Vec<String>,
}

struct RawTerminalMode {
    saved_state: String,
}

impl LineEditor {
    fn new() -> Self {
        Self {
            history: Vec::new(),
        }
    }

    fn read_line(&mut self, prompt: &str, use_history: bool) -> Result<Option<String>, String> {
        if !use_history || !io::stdin().is_terminal() || !io::stdout().is_terminal() {
            return read_line_fallback(prompt);
        }

        let _raw_mode = RawTerminalMode::enter()?;
        print!("{prompt}");
        io::stdout()
            .flush()
            .map_err(|err| format!("failed to flush stdout: {err}"))?;

        let mut input = io::stdin();
        let mut buffer = Vec::<char>::new();
        let mut cursor = 0usize;
        let mut history_index = None::<usize>;
        let mut draft = Vec::<char>::new();

        loop {
            let Some(byte) = read_byte(&mut input)? else {
                println!();
                return Ok(None);
            };

            match byte {
                b'\r' | b'\n' => {
                    println!();
                    let line = buffer.iter().collect::<String>();
                    if use_history && !line.trim().is_empty() && self.history.last() != Some(&line)
                    {
                        self.history.push(line.clone());
                    }
                    return Ok(Some(line));
                }
                3 => {
                    println!("^C");
                    return Ok(Some(String::new()));
                }
                4 if buffer.is_empty() => {
                    println!();
                    return Ok(None);
                }
                8 | 127 => {
                    if cursor > 0 {
                        cursor -= 1;
                        buffer.remove(cursor);
                        render_line(prompt, &buffer, cursor)?;
                    }
                }
                27 => {
                    handle_escape_sequence(
                        &mut input,
                        prompt,
                        &mut buffer,
                        &mut cursor,
                        &mut history_index,
                        &mut draft,
                        &self.history,
                        use_history,
                    )?;
                }
                byte if byte.is_ascii_graphic() || byte == b' ' => {
                    history_index = None;
                    if cursor == buffer.len() {
                        buffer.push(byte as char);
                        cursor += 1;
                        print!("{}", byte as char);
                        io::stdout()
                            .flush()
                            .map_err(|err| format!("failed to flush stdout: {err}"))?;
                    } else {
                        buffer.insert(cursor, byte as char);
                        cursor += 1;
                        render_line(prompt, &buffer, cursor)?;
                    }
                }
                _ => {}
            }
        }
    }
}

impl RawTerminalMode {
    fn enter() -> Result<Self, String> {
        let output = Command::new("stty")
            .arg("-g")
            .stdin(Stdio::inherit())
            .output()
            .map_err(|err| format!("failed to query terminal mode: {err}"))?;
        if !output.status.success() {
            return Err("failed to query terminal mode with stty".to_string());
        }
        let saved_state = String::from_utf8_lossy(&output.stdout).trim().to_string();
        let status = Command::new("stty")
            .args(["raw", "-echo", "min", "1", "time", "0"])
            .stdin(Stdio::inherit())
            .status()
            .map_err(|err| format!("failed to enter raw terminal mode: {err}"))?;
        if !status.success() {
            return Err("failed to enter raw terminal mode with stty".to_string());
        }
        Ok(Self { saved_state })
    }
}

impl Drop for RawTerminalMode {
    fn drop(&mut self) {
        let _ = Command::new("stty")
            .arg(&self.saved_state)
            .stdin(Stdio::inherit())
            .status();
    }
}

fn read_line_fallback(prompt: &str) -> Result<Option<String>, String> {
    print!("{prompt}");
    io::stdout()
        .flush()
        .map_err(|err| format!("failed to flush stdout: {err}"))?;
    let mut line = String::new();
    let bytes = io::stdin()
        .read_line(&mut line)
        .map_err(|err| format!("failed to read input: {err}"))?;
    if bytes == 0 {
        Ok(None)
    } else {
        Ok(Some(line.trim_end_matches(['\r', '\n']).to_string()))
    }
}

fn read_byte(input: &mut io::Stdin) -> Result<Option<u8>, String> {
    let mut byte = [0u8; 1];
    match input
        .read(&mut byte)
        .map_err(|err| format!("failed to read input: {err}"))?
    {
        0 => Ok(None),
        _ => Ok(Some(byte[0])),
    }
}

fn handle_escape_sequence(
    input: &mut io::Stdin,
    prompt: &str,
    buffer: &mut Vec<char>,
    cursor: &mut usize,
    history_index: &mut Option<usize>,
    draft: &mut Vec<char>,
    history: &[String],
    use_history: bool,
) -> Result<(), String> {
    let Some(first) = read_byte(input)? else {
        return Ok(());
    };
    if first != b'[' {
        return Ok(());
    }
    let Some(second) = read_byte(input)? else {
        return Ok(());
    };

    match second {
        b'A' if use_history && !history.is_empty() => {
            if history_index.is_none() {
                *draft = buffer.clone();
                *history_index = Some(history.len() - 1);
            } else if let Some(index) = history_index.as_mut() {
                *index = index.saturating_sub(1);
            }
            if let Some(index) = *history_index {
                *buffer = history[index].chars().collect();
                *cursor = buffer.len();
                render_line(prompt, buffer, *cursor)?;
            }
        }
        b'B' if use_history => {
            if let Some(index) = history_index.as_mut() {
                if *index + 1 < history.len() {
                    *index += 1;
                    *buffer = history[*index].chars().collect();
                } else {
                    *history_index = None;
                    *buffer = draft.clone();
                }
                *cursor = buffer.len();
                render_line(prompt, buffer, *cursor)?;
            }
        }
        b'C' => {
            if *cursor < buffer.len() {
                *cursor += 1;
                render_line(prompt, buffer, *cursor)?;
            }
        }
        b'D' => {
            if *cursor > 0 {
                *cursor -= 1;
                render_line(prompt, buffer, *cursor)?;
            }
        }
        b'H' => {
            *cursor = 0;
            render_line(prompt, buffer, *cursor)?;
        }
        b'F' => {
            *cursor = buffer.len();
            render_line(prompt, buffer, *cursor)?;
        }
        b'3' => {
            if read_byte(input)? == Some(b'~') && *cursor < buffer.len() {
                buffer.remove(*cursor);
                render_line(prompt, buffer, *cursor)?;
            }
        }
        _ => {}
    }

    Ok(())
}

fn render_line(prompt: &str, buffer: &[char], cursor: usize) -> Result<(), String> {
    let rendered = buffer.iter().collect::<String>();
    print!("\r{prompt}{rendered}\x1b[K");
    let right = buffer.len().saturating_sub(cursor);
    if right > 0 {
        print!("\x1b[{right}D");
    }
    io::stdout()
        .flush()
        .map_err(|err| format!("failed to flush stdout: {err}"))
}

pub fn run() -> Result<(), String> {
    println!("blowup-script: projective hypersurface blowups over Q");
    println!("Type 'help' after initialization for commands.\n");

    let mut editor = LineEditor::new();
    let Some(mut state) = initialize_state(&mut editor)? else {
        return Ok(());
    };
    let mut undo_stack = Vec::<AppState>::new();
    print_singular_reports(&state);

    loop {
        let Some(line) = editor.read_line("blowup> ", true)? else {
            break;
        };
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }

        match handle_command(&mut state, &line, &mut undo_stack) {
            Ok(CommandFlow::Continue) => {}
            Ok(CommandFlow::Quit) => break,
            Err(err) => println!("error: {err}"),
        }
    }

    Ok(())
}

enum CommandFlow {
    Continue,
    Quit,
}

fn initialize_state(editor: &mut LineEditor) -> Result<Option<AppState>, String> {
    loop {
        let dimension = prompt(
            editor,
            "dimension (1 curve in P^2, 2 surface in P^3, 3 threefold in P^4): ",
        )?;
        if is_quit_command(&dimension) {
            return Ok(None);
        }
        let projective_dim = match dimension.trim().parse::<usize>() {
            Ok(1..=3) => dimension.trim().parse::<usize>().unwrap(),
            _ => {
                println!("Enter 1, 2, or 3.");
                continue;
            }
        };

        let polynomial = prompt(editor, "homogeneous polynomial: ")?;
        if is_quit_command(&polynomial) {
            return Ok(None);
        }
        match AppState::new(projective_dim, &polynomial) {
            Ok(state) => return Ok(Some(state)),
            Err(err) => {
                println!("error: {err}");
                println!("Try again.\n");
            }
        }
    }
}

fn prompt(editor: &mut LineEditor, label: &str) -> Result<String, String> {
    editor
        .read_line(label, false)?
        .map(|line| line.trim().to_string())
        .ok_or_else(|| "input ended before initialization completed".to_string())
}

fn is_quit_command(line: &str) -> bool {
    matches!(line.trim(), "quit" | "exit")
}

fn handle_command(
    state: &mut AppState,
    line: &str,
    undo_stack: &mut Vec<AppState>,
) -> Result<CommandFlow, String> {
    let parts = line.split_whitespace().collect::<Vec<_>>();
    match parts.as_slice() {
        ["help"] => {
            print_help();
            Ok(CommandFlow::Continue)
        }
        ["charts"] => {
            print_charts(state);
            Ok(CommandFlow::Continue)
        }
        ["singular"] | ["sing"] => {
            print_singular_reports(state);
            Ok(CommandFlow::Continue)
        }
        ["intersections"] | ["intersect"] => {
            print_intersections(state);
            Ok(CommandFlow::Continue)
        }
        ["analyze", rest @ ..] => {
            handle_analyze(state, rest)?;
            Ok(CommandFlow::Continue)
        }
        ["resolve", rest @ ..] | ["auto-resolve", rest @ ..] => {
            let options = parse_auto_resolution_options(rest)?;
            let previous = state.clone();
            println!(
                "Resolving crepant double-cover branch centers (max {} {})",
                options.max_steps,
                pluralize("step", options.max_steps)
            );
            let result =
                state.resolve_crepant_with_progress(options, print_auto_resolution_progress)?;
            let changed = !result.steps.is_empty();
            print_auto_resolution_result(&result);
            if changed {
                undo_stack.push(previous);
            }
            Ok(CommandFlow::Continue)
        }
        ["undo"] => {
            if let Some(previous) = undo_stack.pop() {
                *state = previous;
                println!("undid previous state-changing command");
                print_singular_reports(state);
            } else {
                println!("nothing to undo");
            }
            Ok(CommandFlow::Continue)
        }
        ["set-bound", value] => {
            let bound = value
                .parse::<i32>()
                .map_err(|_| format!("invalid search bound '{value}'"))?;
            if bound < 1 {
                return Err("search bound must be at least 1".to_string());
            }
            let previous = state.clone();
            state.search_bound = bound;
            undo_stack.push(previous);
            println!("rational point search bound set to {bound}");
            Ok(CommandFlow::Continue)
        }
        ["blowup", rest @ ..] => {
            let previous = state.clone();
            handle_blowup(state, rest)?;
            undo_stack.push(previous);
            print_singular_reports(state);
            Ok(CommandFlow::Continue)
        }
        ["save", path] => {
            let json = state_to_json(state);
            std::fs::write(Path::new(path), json)
                .map_err(|err| format!("failed to save '{path}': {err}"))?;
            println!("saved {path}");
            Ok(CommandFlow::Continue)
        }
        ["save-raw", path] => {
            let json = state_to_json_without_analysis(state);
            std::fs::write(Path::new(path), json)
                .map_err(|err| format!("failed to save '{path}': {err}"))?;
            println!("saved {path} without a singular-component audit");
            Ok(CommandFlow::Continue)
        }
        ["quit"] | ["exit"] => Ok(CommandFlow::Quit),
        [unknown, ..] => Err(format!("unknown command '{unknown}'. Type 'help'.")),
        [] => Ok(CommandFlow::Continue),
    }
}

fn handle_blowup(state: &mut AppState, parts: &[&str]) -> Result<(), String> {
    let force = parts.contains(&"--force");
    let filtered = parts
        .iter()
        .copied()
        .filter(|part| *part != "--force")
        .collect::<Vec<_>>();

    match filtered.as_slice() {
        ["component", index] => {
            let index = index
                .parse::<usize>()
                .map_err(|_| format!("invalid component index '{index}'"))?;
            state.blowup_singular_component(index, force)?;
            println!("blew up singular component {index}");
            Ok(())
        }
        ["option", index] => {
            let index = index
                .parse::<usize>()
                .map_err(|_| format!("invalid option index '{index}'"))?;
            state.blowup_option(index, force)?;
            println!("blew up singular point option {index}");
            Ok(())
        }
        ["intersection", index] => {
            let index = index
                .parse::<usize>()
                .map_err(|_| format!("invalid intersection index '{index}'"))?;
            state.blowup_intersection(index, force)?;
            println!("blew up singular component intersection {index}");
            Ok(())
        }
        ["point", chart_id, values @ ..] => {
            let chart_id = parse_chart_id(chart_id)?;
            let chart = state
                .chart(chart_id)
                .ok_or_else(|| format!("chart {chart_id} not found"))?
                .clone();
            let assignments = parse_point_assignments(&chart, values)?;
            state.blowup_coordinate_center(chart_id, assignments, force)?;
            println!("blew up point on chart {chart_id}");
            Ok(())
        }
        ["center", chart_id, assignments @ ..] => {
            let chart_id = parse_chart_id(chart_id)?;
            if assignments.is_empty() {
                return Err("center requires assignments such as x1=0 x2=0".to_string());
            }
            let assignments = parse_assignments(assignments)?;
            state.blowup_coordinate_center(chart_id, assignments, force)?;
            println!("blew up coordinate center on chart {chart_id}");
            Ok(())
        }
        _ => Err(
            "usage: blowup component N | blowup intersection N | blowup option N | blowup point CHART values... | blowup center CHART var=value ..."
                .to_string(),
        ),
    }
}

fn parse_auto_resolution_options(parts: &[&str]) -> Result<AutoResolutionOptions, String> {
    let mut options = AutoResolutionOptions::default();
    let mut index = 0usize;

    if parts.get(index) == Some(&"crepant") {
        index += 1;
    }

    while index < parts.len() {
        match parts[index] {
            "--max-steps" => {
                index += 1;
                let value = parts
                    .get(index)
                    .ok_or_else(|| "--max-steps requires a value".to_string())?;
                options.max_steps = parse_positive_usize(value, "max steps")?;
            }
            value if value.starts_with("--max-steps=") => {
                let value = value
                    .split_once('=')
                    .map(|(_, value)| value)
                    .unwrap_or_default();
                options.max_steps = parse_positive_usize(value, "max steps")?;
            }
            unknown => {
                return Err(format!(
                    "unknown resolve option '{unknown}'. Usage: resolve [crepant] [--max-steps N]"
                ));
            }
        }
        index += 1;
    }

    Ok(options)
}

fn parse_positive_usize(input: &str, label: &str) -> Result<usize, String> {
    let value = input
        .parse::<usize>()
        .map_err(|_| format!("invalid {label} '{input}'"))?;
    if value == 0 {
        return Err(format!("{label} must be at least 1"));
    }
    Ok(value)
}

fn handle_analyze(state: &AppState, parts: &[&str]) -> Result<(), String> {
    match parts {
        ["point", chart_id, values @ ..] => {
            let chart_id = parse_chart_id(chart_id)?;
            let chart = state
                .chart(chart_id)
                .ok_or_else(|| format!("chart {chart_id} not found"))?;
            let assignments = parse_point_assignments(chart, values)?;
            let singularity_type = state.analyze_point(chart_id, assignments.clone())?;
            println!("point analysis on chart {chart_id}");
            println!("  point: {}", format_assignments(&assignments));
            println!("  type: {}", format_singularity_type(&singularity_type));
            Ok(())
        }
        ["intersection", index] => {
            let index = index
                .parse::<usize>()
                .map_err(|_| format!("invalid intersection index '{index}'"))?;
            if index == 0 {
                return Err("intersection numbering starts at 1".to_string());
            }
            let intersections = state.singular_locus().intersections;
            let mut intersection = intersections
                .get(index - 1)
                .cloned()
                .ok_or_else(|| format!("no singular component intersection {index}"))?;
            if let Some(assignments) = &intersection.coordinate_assignments {
                if assignments.len() == intersection.variables.len() {
                    intersection.singularity_type =
                        Some(state.analyze_point(intersection.chart_id, assignments.clone())?);
                }
            }
            print_intersection(index, &intersection);
            Ok(())
        }
        _ => Err(
            "usage: analyze point CHART values... | analyze point CHART var=value ... | analyze intersection N"
                .to_string(),
        ),
    }
}

fn parse_chart_id(input: &str) -> Result<usize, String> {
    input
        .parse::<usize>()
        .map_err(|_| format!("invalid chart id '{input}'"))
}

fn parse_assignments(parts: &[&str]) -> Result<BTreeMap<String, Rational>, String> {
    let mut assignments = BTreeMap::new();
    for part in parts {
        let (name, value) = part
            .split_once('=')
            .ok_or_else(|| format!("expected assignment var=value, got '{part}'"))?;
        if name.trim().is_empty() {
            return Err(format!("empty variable name in '{part}'"));
        }
        assignments.insert(name.trim().to_string(), Rational::parse(value.trim())?);
    }
    Ok(assignments)
}

fn parse_point_assignments(
    chart: &Chart,
    parts: &[&str],
) -> Result<BTreeMap<String, Rational>, String> {
    if parts.is_empty() {
        return Err(format!(
            "chart {} expects {} point coordinates in order: {}",
            chart.id,
            chart.variables.len(),
            chart.variables.join(", ")
        ));
    }

    let mut assignments = if parts.iter().any(|part| part.contains('=')) {
        if !parts.iter().all(|part| part.contains('=')) {
            return Err(
                "use either positional coordinates or named coordinates, not a mix".to_string(),
            );
        }
        parse_assignments(parts)?
    } else {
        if parts.len() != chart.variables.len() {
            return Err(format!(
                "chart {} expects {} point coordinates in order: {}",
                chart.id,
                chart.variables.len(),
                chart.variables.join(", ")
            ));
        }
        chart
            .variables
            .iter()
            .cloned()
            .zip(parts.iter())
            .map(|(name, value)| Rational::parse(value).map(|parsed| (name, parsed)))
            .collect::<Result<BTreeMap<_, _>, _>>()?
    };

    for name in assignments.keys() {
        if !chart.variables.contains(name) {
            return Err(format!("unknown variable '{name}' on chart {}", chart.id));
        }
    }
    let missing = chart
        .variables
        .iter()
        .filter(|name| !assignments.contains_key(*name))
        .cloned()
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        return Err(format!(
            "point on chart {} is missing coordinates for: {}",
            chart.id,
            missing.join(", ")
        ));
    }

    assignments.retain(|name, _| chart.variables.contains(name));
    Ok(assignments)
}

fn print_help() {
    println!(
        "\
commands:
  singular                 recompute singular-locus component decomposition
  intersections            show rational intersection points of singular components
  charts                   show active affine charts and equations
  analyze point C a b ...  analyze local singularity type at a point
  analyze point C v=a ...  analyze local singularity type with named coordinates
  analyze intersection N   analyze a displayed component intersection point
  resolve                  automatically perform crepant branch blowups
  resolve crepant          same as resolve; accepts --max-steps N
  blowup component N       blow up singular component N on every listed chart
  blowup intersection N    blow up displayed singular component intersection point
  blowup option N          blow up detected singular rational point N
  blowup point C a b ...   blow up point on chart C in chart variable order
  blowup point C v=a ...   blow up point on chart C with named coordinates
  blowup center C v=a ...  blow up coordinate center v=a, ...
  set-bound N              set rational point search height for singular samples
  undo                     restore the state before the previous state-changing command
  save FILE                save standardized JSON result
  save-raw FILE            checkpoint JSON without component decomposition
  quit, exit               exit

Add --force to a blowup command to skip singular-locus containment validation.
In an interactive terminal, use Up/Down for command history and Left/Right to edit."
    );
}

fn print_auto_resolution_result(result: &AutoResolutionResult) {
    println!("\nCrepant resolution");
    if result.steps.is_empty() {
        println!("  no automatic blowups performed");
    } else {
        for (index, step) in result.steps.iter().enumerate() {
            println!("  Step {}", index + 1);
            println!(
                "    center: {}",
                format_auto_resolution_center(&step.center)
            );
            println!(
                "    input charts: {}",
                format_usize_list(&step.input_chart_ids)
            );
            println!(
                "    dimension {}, ambient codimension {}, multiplicity {}",
                step.center_dimension, step.ambient_codimension, step.multiplicity
            );
            println!(
                "    blowup stages: {}",
                format_usize_list(&step.blowup_stages)
            );
            println!(
                "    output charts: {}",
                format_usize_list(&step.output_chart_ids)
            );
        }
    }

    match result.status {
        AutoResolutionStatus::Resolved => {
            println!("  status: resolved; no singular components remain")
        }
        AutoResolutionStatus::NoCrepantCenter => println!(
            "  status: stopped; {} singular component(s) remain and no supported crepant center was found",
            result.remaining_components
        ),
        AutoResolutionStatus::StepLimitReached => println!(
            "  status: stopped at step limit; {} singular component(s) and {} intersection(s) remain",
            result.remaining_components, result.remaining_intersections
        ),
    }
    println!();
}

fn print_auto_resolution_progress(progress: AutoResolutionProgress) {
    match progress {
        AutoResolutionProgress::CheckingSingularLocus {
            completed_steps,
            max_steps,
        } => {
            if completed_steps == 0 {
                println!("  checking singular locus before step 1/{max_steps}...");
            } else {
                println!(
                    "  checking singular locus after {} {}...",
                    completed_steps,
                    pluralize("step", completed_steps)
                );
            }
        }
        AutoResolutionProgress::SingularLocusFound {
            components,
            intersections,
            ..
        } => println!(
            "  found {} singular {} and {} {}",
            components,
            pluralize("component", components),
            intersections,
            pluralize("intersection", intersections)
        ),
        AutoResolutionProgress::BlowingUp {
            step,
            max_steps,
            center,
            center_dimension,
            ambient_codimension,
            multiplicity,
            input_chart_ids,
        } => {
            println!(
                "  step {step}/{max_steps}: blowing up {}",
                format_auto_resolution_center(&center)
            );
            println!(
                "    input charts: {}; dimension {}, ambient codimension {}, multiplicity {}",
                format_usize_list(&input_chart_ids),
                center_dimension,
                ambient_codimension,
                multiplicity
            );
        }
        AutoResolutionProgress::StepFinished {
            step,
            blowup_stages,
            output_chart_ids,
        } => {
            println!(
                "    finished step {step}; blowup stages: {}; output charts: {}",
                format_usize_list(&blowup_stages),
                format_usize_list(&output_chart_ids)
            );
        }
    }
    let _ = io::stdout().flush();
}

fn format_auto_resolution_center(center: &AutoResolutionCenter) -> String {
    match center {
        AutoResolutionCenter::Component { index } => format!("component {index}"),
        AutoResolutionCenter::Intersection { index } => format!("intersection {index}"),
    }
}

fn pluralize(noun: &str, count: usize) -> String {
    if count == 1 {
        noun.to_string()
    } else {
        format!("{noun}s")
    }
}

fn format_usize_list(values: &[usize]) -> String {
    if values.is_empty() {
        return "(none)".to_string();
    }
    values
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ")
}

fn print_charts(state: &AppState) {
    println!("active charts:");
    for chart in state.active_charts() {
        println!(
            "  chart {} [{}], vars ({})",
            chart.id,
            chart.label,
            chart.variables.join(", ")
        );
        println!("    f = {}", chart.polynomial);
    }
}

fn print_singular_reports(state: &AppState) {
    println!("\nSingular locus");
    let analysis = state.singular_locus();
    let components = &analysis.components;
    if components.is_empty() {
        println!("  no singular components found");
        println!();
        return;
    }

    println!("  components: {}", components.len());
    if !analysis.intersections.is_empty() {
        println!("  intersections: {}", analysis.intersections.len());
    }
    println!();

    for (index, component) in components.iter().enumerate() {
        print_component(index + 1, component);
    }
    if !analysis.intersections.is_empty() {
        println!("Intersections");
        print_intersections_from_list(&analysis.intersections);
    }
    println!("Next: use `blowup component N` or `blowup intersection N`.");
    println!();
}

fn print_component(index: usize, component: &SingularComponent) {
    println!("Component {index}");
    println!("  dimension: {}", component.dimension);
    println!("  generic multiplicity: {}", component.multiplicity);
    if let Some(singularity_type) = &component.generic_singularity_type {
        if let Some(chart_id) = component.generic_singularity_chart {
            println!("  generic type on chart {chart_id}:");
        } else {
            println!("  generic type:");
        }
        print_singularity_type_details("    ", singularity_type);
        if let Some(point) = &component.generic_singularity_point {
            println!("    sample: {}", format_assignments(point));
        }
    }
    if component.charts.len() > COMPONENT_POLY_DETAIL_CHART_LIMIT {
        let chart_ids = component
            .charts
            .iter()
            .map(|chart| chart.chart_id.to_string())
            .collect::<Vec<_>>()
            .join(", ");
        println!(
            "  charts: {chart_ids} ({} charts; ideals omitted)",
            component.charts.len()
        );
        println!();
        return;
    }
    println!("  charts:");
    for chart in &component.charts {
        println!("    chart {} [{}]", chart.chart_id, chart.chart_label);
        println!("      variables: {}", chart.variables.join(", "));
        println!("      ideal: {}", format_component_polys(chart));
        if chart.dimension != component.dimension || chart.multiplicity != component.multiplicity {
            println!(
                "      local dimension: {}, local multiplicity: {}",
                chart.dimension, chart.multiplicity
            );
        }
        if let Some(assignments) = &chart.coordinate_assignments {
            if !assignments.is_empty() {
                println!("      center: {}", format_assignments(assignments));
            }
        } else if chart.affine_linear_center {
            let rendered = chart
                .ideal_generators
                .iter()
                .enumerate()
                .map(|(index, poly)| format!("f{index}=0 ({poly})"))
                .collect::<Vec<_>>()
                .join(", ");
            println!("      affine-linear center: {rendered}");
        }
        if let Some(singularity_type) = &chart.singularity_type {
            println!("      local type:");
            print_singularity_type_details("        ", singularity_type);
        }
    }
    println!();
}

fn format_component_polys(chart: &SingularComponentChart) -> String {
    if chart.ideal_generators.is_empty() {
        return "(0)".to_string();
    }

    chart
        .ideal_generators
        .iter()
        .enumerate()
        .map(|(index, poly)| format!("f{index} = {poly}"))
        .collect::<Vec<_>>()
        .join("; ")
}

fn print_intersections(state: &AppState) {
    println!("\nSingular component intersections");
    let intersections = state.singular_locus().intersections;
    if intersections.is_empty() {
        println!("  no rational component intersection points found");
        println!();
        return;
    }

    print_intersections_from_list(&intersections);
    println!();
}

fn print_intersections_from_list(intersections: &[SingularIntersection]) {
    for (index, intersection) in intersections.iter().enumerate() {
        print_intersection(index + 1, intersection);
    }
}

fn print_intersection(index: usize, intersection: &SingularIntersection) {
    println!("  Intersection {index}");
    println!(
        "    components: {}",
        intersection
            .component_indices
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ")
    );
    println!(
        "    chart: {} [{}]",
        intersection.chart_id, intersection.chart_label
    );
    if let Some(assignments) = &intersection.coordinate_assignments {
        println!("    point: {}", format_assignments(assignments));
    } else {
        println!(
            "    isolated ideal: {}",
            intersection
                .ideal_generators
                .iter()
                .enumerate()
                .map(|(generator_index, generator)| {
                    format!("f{generator_index}=0 ({generator})")
                })
                .collect::<Vec<_>>()
                .join("; ")
        );
    }
    println!("    variables: {}", intersection.variables.join(", "));
    if let Some(singularity_type) = &intersection.singularity_type {
        println!("    local type:");
        print_singularity_type_details("      ", singularity_type);
    }
}

fn format_assignments(assignments: &BTreeMap<String, Rational>) -> String {
    assignments
        .iter()
        .map(|(name, value)| format!("{name}={value}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn format_singularity_type(singularity_type: &SingularityType) -> String {
    let mut parts = vec![
        singularity_type.label.clone(),
        format!("multiplicity {}", singularity_type.multiplicity),
    ];
    if let Some(degree) = singularity_type.tangent_cone_degree {
        parts.push(format!(
            "tangent cone degree {degree}: {}",
            format_polynomial_summary(&singularity_type.tangent_cone)
        ));
    }
    if let Some(rank) = singularity_type.quadratic_rank {
        parts.push(format!("quadratic rank {rank}"));
    }
    parts.join("; ")
}

fn print_singularity_type_details(indent: &str, singularity_type: &SingularityType) {
    println!("{indent}label: {}", singularity_type.label);
    println!("{indent}multiplicity: {}", singularity_type.multiplicity);
    if let Some(degree) = singularity_type.tangent_cone_degree {
        println!(
            "{indent}tangent cone: degree {degree}, {}",
            format_polynomial_summary(&singularity_type.tangent_cone)
        );
    }
    if let Some(rank) = singularity_type.quadratic_rank {
        println!("{indent}quadratic rank: {rank}");
    }
}

fn format_polynomial_summary(poly: &Poly) -> String {
    let rendered = poly.to_string();
    if rendered.len() <= POLYNOMIAL_DISPLAY_LIMIT {
        rendered
    } else {
        format!("polynomial omitted ({} terms)", poly.terms.len())
    }
}
