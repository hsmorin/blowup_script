use crate::grobner::{
    decompose_ideal, generic_multiplicity, groebner_key, ideal_contains_all, ideal_dimension,
    reduced_groebner_basis,
};
use crate::parser::parse_polynomial;
use crate::poly::Poly;
use crate::rational::Rational;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::process::{Command, Stdio};

#[derive(Clone, Debug)]
pub struct AppState {
    pub projective_dim: usize,
    pub projective_vars: Vec<String>,
    pub initial_polynomial_text: String,
    pub initial_polynomial: Poly,
    pub degree: usize,
    pub charts: Vec<Chart>,
    pub blowups: Vec<BlowupRecord>,
    pub next_chart_id: usize,
    pub search_bound: i32,
    cached_singular_analysis: Option<SingularAnalysis>,
}

#[derive(Clone, Debug)]
pub struct Chart {
    pub id: usize,
    pub parent: Option<usize>,
    pub label: String,
    pub variables: Vec<String>,
    pub polynomial: Poly,
    pub active: bool,
    pub substitutions: Vec<String>,
    projective_coordinates: Vec<Poly>,
    blowup_coordinates: Vec<BlowupCoordinateChart>,
}

#[derive(Clone, Debug)]
struct BlowupCoordinateChart {
    stage: usize,
    coordinates: Vec<Poly>,
}

#[derive(Clone, Debug)]
pub struct BlowupRecord {
    pub stage: usize,
    pub input_chart: usize,
    pub center: CoordinateCenter,
    pub multiplicity: usize,
    pub output_charts: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct CoordinateCenter {
    pub assignments: BTreeMap<String, Rational>,
}

#[derive(Clone, Debug)]
pub struct SingularityReport {
    pub chart_id: usize,
    pub chart_label: String,
    pub variables: Vec<String>,
    pub equations: Vec<Poly>,
    pub point_options: Vec<SingularityPoint>,
    pub sample_limit_hit: bool,
}

#[derive(Clone, Debug)]
pub struct SingularComponent {
    pub dimension: usize,
    pub multiplicity: usize,
    pub generic_singularity_type: Option<SingularityType>,
    pub generic_singularity_chart: Option<usize>,
    pub generic_singularity_point: Option<BTreeMap<String, Rational>>,
    pub charts: Vec<SingularComponentChart>,
}

#[derive(Clone, Debug)]
pub struct SingularComponentChart {
    pub chart_id: usize,
    pub chart_label: String,
    pub variables: Vec<String>,
    pub ideal_generators: Vec<Poly>,
    pub dimension: usize,
    pub multiplicity: usize,
    pub coordinate_assignments: Option<BTreeMap<String, Rational>>,
    pub affine_linear_center: bool,
    pub singularity_type: Option<SingularityType>,
}

#[derive(Clone, Debug)]
pub struct SingularAnalysis {
    pub components: Vec<SingularComponent>,
    pub intersections: Vec<SingularIntersection>,
}

#[derive(Clone, Debug)]
pub struct SingularIntersection {
    pub chart_id: usize,
    pub chart_label: String,
    pub variables: Vec<String>,
    pub component_indices: Vec<usize>,
    pub ideal_generators: Vec<Poly>,
    pub coordinate_assignments: Option<BTreeMap<String, Rational>>,
    pub singularity_type: Option<SingularityType>,
}

#[derive(Clone, Debug)]
pub struct SingularityType {
    pub label: String,
    pub multiplicity: usize,
    pub tangent_cone: Poly,
    pub tangent_cone_degree: Option<usize>,
    pub quadratic_rank: Option<usize>,
    pub embedding_dimension: usize,
    pub is_singular: bool,
}

#[derive(Clone, Debug)]
pub struct SingularityPoint {
    pub chart_id: usize,
    pub values: Vec<Rational>,
}

#[derive(Clone, Debug)]
pub struct AutoResolutionOptions {
    pub max_steps: usize,
}

#[derive(Clone, Debug)]
pub struct AutoResolutionResult {
    pub status: AutoResolutionStatus,
    pub steps: Vec<AutoResolutionStep>,
    pub remaining_components: usize,
    pub remaining_intersections: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AutoResolutionStatus {
    Resolved,
    NoCrepantCenter,
    StepLimitReached,
}

#[derive(Clone, Debug)]
pub struct AutoResolutionStep {
    pub center: AutoResolutionCenter,
    pub center_dimension: usize,
    pub ambient_codimension: usize,
    pub multiplicity: usize,
    pub input_chart_ids: Vec<usize>,
    pub blowup_stages: Vec<usize>,
    pub output_chart_ids: Vec<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AutoResolutionCenter {
    Component { index: usize },
    Intersection { index: usize },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AutoResolutionProgress {
    CheckingSingularLocus {
        completed_steps: usize,
        max_steps: usize,
    },
    SingularLocusFound {
        completed_steps: usize,
        components: usize,
        intersections: usize,
    },
    BlowingUp {
        step: usize,
        max_steps: usize,
        center: AutoResolutionCenter,
        center_dimension: usize,
        ambient_codimension: usize,
        multiplicity: usize,
        input_chart_ids: Vec<usize>,
    },
    StepFinished {
        step: usize,
        blowup_stages: Vec<usize>,
        output_chart_ids: Vec<usize>,
    },
}

#[derive(Clone, Debug)]
enum ComponentCenter {
    Coordinate {
        chart_id: usize,
        assignments: BTreeMap<String, Rational>,
        multiplicity: usize,
    },
    AffineLinear {
        chart_id: usize,
        generators: Vec<Poly>,
        multiplicity: usize,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum BlowupTransform {
    Strict,
    DoubleCoverBranch,
}

impl BlowupTransform {
    fn exceptional_division_power(self, multiplicity: usize) -> usize {
        match self {
            Self::Strict => multiplicity,
            Self::DoubleCoverBranch => double_cover_branch_division_power(multiplicity),
        }
    }
}

#[derive(Clone, Debug)]
struct ProjectiveLinearComponent {
    generators: Vec<Poly>,
    dimension: usize,
}

#[derive(Clone, Debug)]
struct AutoResolutionCandidate {
    center: AutoResolutionCenter,
    center_dimension: usize,
    ambient_codimension: usize,
    multiplicity: usize,
    input_chart_ids: Vec<usize>,
    score: (usize, usize, usize, usize),
}

#[derive(Clone, Copy, Debug)]
struct SingularComputationOptions {
    include_singularity_types: bool,
    positive_dimensional_only: bool,
}

impl SingularComputationOptions {
    const LOCUS_ONLY: Self = Self {
        include_singularity_types: false,
        positive_dimensional_only: false,
    };

    const WITH_TYPES: Self = Self {
        include_singularity_types: true,
        positive_dimensional_only: false,
    };

    const CREPANT_CURVES_ONLY: Self = Self {
        include_singularity_types: false,
        positive_dimensional_only: true,
    };
}

impl SingularityPoint {
    pub fn assignments(&self, variables: &[String]) -> BTreeMap<String, Rational> {
        variables
            .iter()
            .cloned()
            .zip(self.values.iter().cloned())
            .collect()
    }
}

impl Default for AutoResolutionOptions {
    fn default() -> Self {
        Self { max_steps: 32 }
    }
}

impl AppState {
    pub fn new(projective_dim: usize, polynomial_text: &str) -> Result<Self, String> {
        if !(1..=3).contains(&projective_dim) {
            return Err("dimension must be 1, 2, or 3".to_string());
        }

        let projective_vars = (0..projective_dim + 2)
            .map(|index| format!("x{index}"))
            .collect::<Vec<_>>();
        let initial_polynomial = parse_polynomial(polynomial_text, &projective_vars)?;
        if initial_polynomial.is_zero() {
            return Err("the hypersurface polynomial must be nonzero".to_string());
        }

        let degree = initial_polynomial
            .homogeneous_degree()
            .ok_or_else(|| "projective hypersurface polynomial must be homogeneous".to_string())?;

        let mut charts = Vec::new();
        for chart_index in 0..projective_vars.len() {
            let affine_vars = projective_vars
                .iter()
                .enumerate()
                .filter_map(|(index, name)| {
                    if index == chart_index {
                        None
                    } else {
                        Some(name.clone())
                    }
                })
                .collect::<Vec<_>>();

            let replacements = projective_vars
                .iter()
                .enumerate()
                .map(|(index, name)| {
                    if index == chart_index {
                        Poly::constant(&affine_vars, Rational::one())
                    } else {
                        let affine_index = affine_vars
                            .iter()
                            .position(|affine_name| affine_name == name)
                            .expect("affine variable present");
                        Poly::var(&affine_vars, affine_index)
                    }
                })
                .collect::<Vec<_>>();

            let polynomial = initial_polynomial.substitute(&replacements);
            charts.push(Chart {
                id: chart_index,
                parent: None,
                label: format!("{} != 0", projective_vars[chart_index]),
                variables: affine_vars,
                polynomial,
                active: true,
                substitutions: vec![format!("{} = 1", projective_vars[chart_index])],
                projective_coordinates: replacements,
                blowup_coordinates: Vec::new(),
            });
        }

        Ok(Self {
            projective_dim,
            projective_vars,
            initial_polynomial_text: polynomial_text.trim().to_string(),
            initial_polynomial,
            degree,
            next_chart_id: charts.len(),
            charts,
            blowups: Vec::new(),
            search_bound: 2,
            cached_singular_analysis: None,
        })
    }

    pub fn active_charts(&self) -> impl Iterator<Item = &Chart> {
        self.charts.iter().filter(|chart| chart.active)
    }

    pub fn singular_reports(&self) -> Vec<SingularityReport> {
        self.active_charts()
            .map(|chart| {
                let mut equations = Vec::with_capacity(chart.variables.len() + 1);
                equations.push(chart.polynomial.clone());
                equations.extend(chart.polynomial.partials());
                let (point_options, sample_limit_hit) =
                    enumerate_singular_points(chart, &equations, self.search_bound, 50);

                SingularityReport {
                    chart_id: chart.id,
                    chart_label: chart.label.clone(),
                    variables: chart.variables.clone(),
                    equations,
                    point_options,
                    sample_limit_hit,
                }
            })
            .collect()
    }

    #[cfg(test)]
    fn singular_analysis(&self) -> SingularAnalysis {
        self.compute_singular_analysis(SingularComputationOptions::WITH_TYPES)
    }

    pub fn singular_locus(&self) -> SingularAnalysis {
        if let Some(analysis) = &self.cached_singular_analysis {
            return analysis.clone();
        }
        self.compute_singular_analysis(SingularComputationOptions::WITH_TYPES)
    }

    #[cfg(test)]
    fn singular_components(&self) -> Vec<SingularComponent> {
        self.compute_singular_components(SingularComputationOptions::WITH_TYPES)
    }

    fn compute_singular_analysis(&self, options: SingularComputationOptions) -> SingularAnalysis {
        let components = self.compute_singular_components(options);
        let intersections = self.compute_singular_component_intersections(&components, options);
        SingularAnalysis {
            components,
            intersections,
        }
    }

    fn compute_singular_components(
        &self,
        options: SingularComputationOptions,
    ) -> Vec<SingularComponent> {
        if let Some(components) = self.projective_linear_singular_components(options) {
            return components;
        }

        let mut grouped = Vec::<SingularComponent>::new();
        let mut key_to_index = BTreeMap::<String, usize>::new();

        let local_components_by_chart = compute_local_singular_components_in_parallel(
            self.active_charts().collect(),
            options.positive_dimensional_only,
        );
        for (chart_id, local_components) in local_components_by_chart {
            let chart = self
                .chart(chart_id)
                .expect("parallel singular computation returned an existing chart");
            for ideal_generators in local_components {
                if ideal_generators.iter().any(Poly::is_nonzero_constant) {
                    continue;
                }

                let dimension =
                    ideal_dimension(&ideal_generators, chart.variables.len()).unwrap_or_default();
                let linear_change = linear_center_change(&ideal_generators, &chart.variables);
                let multiplicity = generic_multiplicity_on_component(
                    chart,
                    &ideal_generators,
                    linear_change.as_ref(),
                );
                let coordinate_assignments =
                    coordinate_assignments_from_basis(&ideal_generators, &chart.variables);
                let affine_linear_center = linear_change.is_some();
                let singularity_type = if options.include_singularity_types {
                    coordinate_assignments
                        .as_ref()
                        .filter(|assignments| assignments.len() == chart.variables.len())
                        .map(|assignments| analyze_chart_point_singularity(chart, assignments))
                } else {
                    None
                };
                let key = self.component_key(
                    chart,
                    &ideal_generators,
                    dimension,
                    &coordinate_assignments,
                );
                let chart_component = SingularComponentChart {
                    chart_id: chart.id,
                    chart_label: chart.label.clone(),
                    variables: chart.variables.clone(),
                    ideal_generators,
                    dimension,
                    multiplicity,
                    coordinate_assignments,
                    affine_linear_center,
                    singularity_type,
                };

                if let Some(index) = key_to_index.get(&key).copied() {
                    if !grouped[index]
                        .charts
                        .iter()
                        .any(|entry| entry.chart_id == chart_component.chart_id)
                    {
                        grouped[index].dimension = grouped[index].dimension.max(dimension);
                        grouped[index].multiplicity = grouped[index].multiplicity.max(multiplicity);
                        grouped[index].charts.push(chart_component);
                    }
                } else {
                    key_to_index.insert(key, grouped.len());
                    grouped.push(SingularComponent {
                        dimension,
                        multiplicity,
                        generic_singularity_type: None,
                        generic_singularity_chart: None,
                        generic_singularity_point: None,
                        charts: vec![chart_component],
                    });
                }
            }
        }

        if options.include_singularity_types {
            for component in &mut grouped {
                if component.dimension == 1 {
                    if let Some((chart_id, point, singularity_type)) =
                        generic_line_singularity_type(component, &self.charts)
                    {
                        component.generic_singularity_type = Some(singularity_type);
                        component.generic_singularity_chart = Some(chart_id);
                        component.generic_singularity_point = Some(point);
                    }
                }
            }
        }

        grouped.sort_by(|left, right| {
            right
                .dimension
                .cmp(&left.dimension)
                .then_with(|| left.charts[0].chart_id.cmp(&right.charts[0].chart_id))
                .then_with(|| component_display_key(left).cmp(&component_display_key(right)))
        });
        grouped
    }

    #[cfg(test)]
    fn singular_component_intersections(&self) -> Vec<SingularIntersection> {
        let components = self.compute_singular_components(SingularComputationOptions::WITH_TYPES);
        self.compute_singular_component_intersections(
            &components,
            SingularComputationOptions::WITH_TYPES,
        )
    }

    fn compute_singular_component_intersections(
        &self,
        components: &[SingularComponent],
        options: SingularComputationOptions,
    ) -> Vec<SingularIntersection> {
        let mut by_key = BTreeMap::<String, SingularIntersection>::new();

        for chart in self.active_charts() {
            for left_index in 0..components.len() {
                let Some(left_chart) = components[left_index]
                    .charts
                    .iter()
                    .find(|component_chart| component_chart.chart_id == chart.id)
                else {
                    continue;
                };

                for right_index in left_index + 1..components.len() {
                    let Some(right_chart) = components[right_index]
                        .charts
                        .iter()
                        .find(|component_chart| component_chart.chart_id == chart.id)
                    else {
                        continue;
                    };

                    let bases = component_intersection_bases(left_chart, right_chart, chart);
                    for basis in bases {
                        if basis.iter().any(Poly::is_nonzero_constant) {
                            continue;
                        }
                        if ideal_dimension(&basis, chart.variables.len()) != Some(0) {
                            continue;
                        }

                        let coordinate_assignments =
                            coordinate_assignments_from_basis(&basis, &chart.variables);
                        let singularity_type = if options.include_singularity_types {
                            coordinate_assignments
                                .as_ref()
                                .filter(|assignments| assignments.len() == chart.variables.len())
                                .map(|assignments| {
                                    analyze_chart_point_singularity(chart, assignments)
                                })
                        } else {
                            None
                        };
                        let key = if let Some(assignments) = &coordinate_assignments {
                            self.global_point_key(chart, assignments)
                                .unwrap_or_else(|| {
                                    format!(
                                        "chart:{}:point:{}",
                                        chart.id,
                                        assignment_key(assignments)
                                    )
                                })
                        } else {
                            format!("chart:{}:ideal:{}", chart.id, groebner_key(&basis))
                        };

                        let entry = by_key.entry(key).or_insert_with(|| SingularIntersection {
                            chart_id: chart.id,
                            chart_label: chart.label.clone(),
                            variables: chart.variables.clone(),
                            component_indices: Vec::new(),
                            ideal_generators: basis.clone(),
                            coordinate_assignments,
                            singularity_type,
                        });
                        for component_index in [left_index + 1, right_index + 1] {
                            if !entry.component_indices.contains(&component_index) {
                                entry.component_indices.push(component_index);
                            }
                        }
                        entry.component_indices.sort();
                    }
                }
            }
        }

        let mut intersections = by_key.into_values().collect::<Vec<_>>();
        intersections.sort_by(|left, right| {
            left.chart_id
                .cmp(&right.chart_id)
                .then_with(|| left.component_indices.cmp(&right.component_indices))
                .then_with(|| {
                    assignment_key_option(&left.coordinate_assignments)
                        .cmp(&assignment_key_option(&right.coordinate_assignments))
                })
                .then_with(|| {
                    groebner_key(&left.ideal_generators).cmp(&groebner_key(&right.ideal_generators))
                })
        });
        intersections
    }

    pub fn singular_point_options(&self) -> Vec<SingularityPoint> {
        self.singular_reports()
            .into_iter()
            .flat_map(|report| report.point_options)
            .collect()
    }

    pub fn analyze_point(
        &self,
        chart_id: usize,
        assignments: BTreeMap<String, Rational>,
    ) -> Result<SingularityType, String> {
        let chart = self
            .chart(chart_id)
            .ok_or_else(|| format!("chart {chart_id} not found"))?;
        require_full_point_assignment(chart, &assignments)?;
        Ok(analyze_chart_point_singularity(chart, &assignments))
    }

    pub fn resolve_crepant_with_progress<F>(
        &mut self,
        options: AutoResolutionOptions,
        mut on_progress: F,
    ) -> Result<AutoResolutionResult, String>
    where
        F: FnMut(AutoResolutionProgress),
    {
        if !(1..=2).contains(&self.projective_dim) {
            return Err(
                "automatic crepant resolution currently supports curves and surfaces only"
                    .to_string(),
            );
        }

        let mut steps = Vec::new();
        loop {
            let completed_steps = steps.len();
            if completed_steps >= options.max_steps
                && std::env::var("BLOWUP_SKIP_STEP_LIMIT_AUDIT")
                    .is_ok_and(|value| matches!(value.as_str(), "1" | "true" | "yes"))
            {
                return Ok(AutoResolutionResult {
                    status: AutoResolutionStatus::StepLimitReached,
                    steps,
                    remaining_components: 0,
                    remaining_intersections: 0,
                });
            }
            on_progress(AutoResolutionProgress::CheckingSingularLocus {
                completed_steps,
                max_steps: options.max_steps,
            });
            let curve_components =
                self.compute_singular_components(SingularComputationOptions::CREPANT_CURVES_ONLY);
            let curve_candidate = curve_components
                .iter()
                .enumerate()
                .filter(|(_, component)| component.dimension > 0)
                .filter_map(|(index, component)| {
                    self.crepant_component_candidate(index + 1, component)
                })
                .min_by_key(|candidate| candidate.score);
            let curve_phase_active = curve_components
                .iter()
                .any(|component| component.dimension > 0);
            let components = if curve_phase_active {
                curve_components
            } else {
                self.compute_singular_components(SingularComputationOptions::LOCUS_ONLY)
            };
            let intersections = if curve_phase_active {
                Vec::new()
            } else {
                self.compute_singular_component_intersections(
                    &components,
                    SingularComputationOptions::LOCUS_ONLY,
                )
            };
            let analysis = SingularAnalysis {
                components,
                intersections,
            };
            on_progress(AutoResolutionProgress::SingularLocusFound {
                completed_steps,
                components: analysis.components.len(),
                intersections: analysis.intersections.len(),
            });
            if analysis.components.is_empty() {
                self.cache_typed_singular_analysis(analysis);
                return Ok(AutoResolutionResult {
                    status: AutoResolutionStatus::Resolved,
                    steps,
                    remaining_components: 0,
                    remaining_intersections: 0,
                });
            }

            let Some(candidate) = curve_candidate
                .clone()
                .or_else(|| self.choose_crepant_center(&analysis))
            else {
                self.cache_typed_singular_analysis(analysis.clone());
                return Ok(AutoResolutionResult {
                    status: AutoResolutionStatus::NoCrepantCenter,
                    steps,
                    remaining_components: analysis.components.len(),
                    remaining_intersections: analysis.intersections.len(),
                });
            };

            if steps.len() >= options.max_steps {
                return Ok(AutoResolutionResult {
                    status: AutoResolutionStatus::StepLimitReached,
                    steps,
                    remaining_components: analysis.components.len(),
                    remaining_intersections: analysis.intersections.len(),
                });
            }

            let step_number = steps.len() + 1;
            on_progress(AutoResolutionProgress::BlowingUp {
                step: step_number,
                max_steps: options.max_steps,
                center: candidate.center.clone(),
                center_dimension: candidate.center_dimension,
                ambient_codimension: candidate.ambient_codimension,
                multiplicity: candidate.multiplicity,
                input_chart_ids: candidate.input_chart_ids.clone(),
            });
            let previous_blowup_count = self.blowups.len();
            match candidate.center {
                AutoResolutionCenter::Component { index } => {
                    self.blowup_singular_component_with_transform(
                        index,
                        false,
                        BlowupTransform::DoubleCoverBranch,
                    )?;
                }
                AutoResolutionCenter::Intersection { index } => {
                    self.blowup_intersection_with_transform(
                        index,
                        false,
                        BlowupTransform::DoubleCoverBranch,
                    )?;
                }
            }

            let step = candidate.into_step(&self.blowups[previous_blowup_count..]);
            on_progress(AutoResolutionProgress::StepFinished {
                step: step_number,
                blowup_stages: step.blowup_stages.clone(),
                output_chart_ids: step.output_chart_ids.clone(),
            });
            steps.push(step);
        }
    }

    fn cache_typed_singular_analysis(&mut self, mut analysis: SingularAnalysis) {
        for component in &mut analysis.components {
            for component_chart in &mut component.charts {
                component_chart.singularity_type = component_chart
                    .coordinate_assignments
                    .as_ref()
                    .filter(|assignments| assignments.len() == component_chart.variables.len())
                    .and_then(|assignments| {
                        self.chart(component_chart.chart_id)
                            .map(|chart| analyze_chart_point_singularity(chart, assignments))
                    });
            }
        }
        for intersection in &mut analysis.intersections {
            intersection.singularity_type = intersection
                .coordinate_assignments
                .as_ref()
                .filter(|assignments| assignments.len() == intersection.variables.len())
                .and_then(|assignments| {
                    self.chart(intersection.chart_id)
                        .map(|chart| analyze_chart_point_singularity(chart, assignments))
                });
        }
        self.cached_singular_analysis = Some(analysis);
    }

    pub fn blowup_intersection(&mut self, display_index: usize, force: bool) -> Result<(), String> {
        self.blowup_intersection_with_transform(display_index, force, BlowupTransform::Strict)
    }

    fn blowup_intersection_with_transform(
        &mut self,
        display_index: usize,
        force: bool,
        transform: BlowupTransform,
    ) -> Result<(), String> {
        if display_index == 0 {
            return Err("intersection numbering starts at 1".to_string());
        }

        let intersections = self
            .compute_singular_analysis(SingularComputationOptions::LOCUS_ONLY)
            .intersections;
        let intersection = intersections
            .get(display_index - 1)
            .ok_or_else(|| format!("no singular component intersection {display_index}"))?;
        let assignments = intersection.coordinate_assignments.clone().ok_or_else(|| {
            format!(
                "intersection {display_index} is zero-dimensional but not a rational coordinate point"
            )
        })?;
        self.blowup_coordinate_center_with_transform(
            intersection.chart_id,
            assignments,
            force,
            transform,
        )
    }

    pub fn blowup_singular_component(
        &mut self,
        display_index: usize,
        force: bool,
    ) -> Result<(), String> {
        self.blowup_singular_component_with_transform(display_index, force, BlowupTransform::Strict)
    }

    fn blowup_singular_component_with_transform(
        &mut self,
        display_index: usize,
        force: bool,
        transform: BlowupTransform,
    ) -> Result<(), String> {
        if display_index == 0 {
            return Err("singular component numbering starts at 1".to_string());
        }

        let computation_options = match transform {
            BlowupTransform::Strict => SingularComputationOptions::LOCUS_ONLY,
            BlowupTransform::DoubleCoverBranch => SingularComputationOptions::CREPANT_CURVES_ONLY,
        };
        let components = self.compute_singular_components(computation_options);
        let component = components
            .get(display_index - 1)
            .ok_or_else(|| format!("no singular component {display_index}"))?;

        let mut centers = Vec::new();
        for chart_component in &component.charts {
            if let Some(assignments) = chart_component.coordinate_assignments.clone() {
                if assignments.len() < 2 {
                    return Err(format!(
                        "component {display_index} on chart {} has fewer than two coordinate equations",
                        chart_component.chart_id
                    ));
                }
                centers.push(ComponentCenter::Coordinate {
                    chart_id: chart_component.chart_id,
                    assignments,
                    multiplicity: chart_component.multiplicity,
                });
            } else if chart_component.affine_linear_center {
                centers.push(ComponentCenter::AffineLinear {
                    chart_id: chart_component.chart_id,
                    generators: chart_component.ideal_generators.clone(),
                    multiplicity: chart_component.multiplicity,
                });
            } else {
                return Err(format!(
                    "component {display_index} on chart {} is not a coordinate or affine-linear center",
                    chart_component.chart_id
                ));
            }
        }

        if centers.is_empty() {
            return Err(format!(
                "singular component {display_index} has no chart centers"
            ));
        }

        let mut next_state = self.clone();
        let stage = next_state.next_blowup_stage();
        for center in centers {
            match center {
                ComponentCenter::Coordinate {
                    chart_id,
                    assignments,
                    multiplicity,
                } => next_state.blowup_coordinate_center_with_transform_at_stage(
                    chart_id,
                    assignments,
                    force,
                    transform,
                    stage,
                    Some(multiplicity),
                )?,
                ComponentCenter::AffineLinear {
                    chart_id,
                    generators,
                    multiplicity,
                } => next_state.blowup_affine_linear_center_at_stage(
                    chart_id,
                    &generators,
                    force,
                    transform,
                    stage,
                    multiplicity,
                )?,
            }
        }
        *self = next_state;
        Ok(())
    }

    pub fn blowup_option(&mut self, option_index: usize, force: bool) -> Result<(), String> {
        let options = self.singular_point_options();
        let option = options
            .get(option_index)
            .ok_or_else(|| format!("no singular point option {option_index}"))?;
        let chart = self
            .chart(option.chart_id)
            .ok_or_else(|| format!("chart {} no longer exists", option.chart_id))?;
        let assignments = option.assignments(&chart.variables);
        self.blowup_coordinate_center(option.chart_id, assignments, force)
    }

    pub fn blowup_coordinate_center(
        &mut self,
        chart_id: usize,
        assignments: BTreeMap<String, Rational>,
        force: bool,
    ) -> Result<(), String> {
        self.blowup_coordinate_center_with_transform(
            chart_id,
            assignments,
            force,
            BlowupTransform::Strict,
        )
    }

    fn blowup_coordinate_center_with_transform(
        &mut self,
        chart_id: usize,
        assignments: BTreeMap<String, Rational>,
        force: bool,
        transform: BlowupTransform,
    ) -> Result<(), String> {
        let stage = self.next_blowup_stage();
        self.blowup_coordinate_center_with_transform_at_stage(
            chart_id,
            assignments,
            force,
            transform,
            stage,
            None,
        )
    }

    fn blowup_coordinate_center_with_transform_at_stage(
        &mut self,
        chart_id: usize,
        assignments: BTreeMap<String, Rational>,
        force: bool,
        transform: BlowupTransform,
        stage: usize,
        known_multiplicity: Option<usize>,
    ) -> Result<(), String> {
        let chart_position = self
            .charts
            .iter()
            .position(|chart| chart.id == chart_id && chart.active)
            .ok_or_else(|| format!("active chart {chart_id} not found"))?;
        let chart = self.charts[chart_position].clone();

        if assignments.len() < 2 {
            return Err(
                "a coordinate blowup center must assign at least two variables".to_string(),
            );
        }

        let mut indexed_assignments = BTreeMap::new();
        for (name, value) in &assignments {
            let index = chart
                .polynomial
                .variable_index(name)
                .ok_or_else(|| format!("unknown variable '{name}' on chart {chart_id}"))?;
            indexed_assignments.insert(index, value.clone());
        }

        self.blowup_prepared_coordinate_center(
            chart_position,
            chart,
            assignments,
            indexed_assignments,
            force,
            Vec::new(),
            transform,
            stage,
            known_multiplicity,
        )
    }

    fn blowup_affine_linear_center_at_stage(
        &mut self,
        chart_id: usize,
        generators: &[Poly],
        force: bool,
        transform: BlowupTransform,
        stage: usize,
        multiplicity: usize,
    ) -> Result<(), String> {
        let chart_position = self
            .charts
            .iter()
            .position(|chart| chart.id == chart_id && chart.active)
            .ok_or_else(|| format!("active chart {chart_id} not found"))?;
        let chart = self.charts[chart_position].clone();
        let change = linear_center_change(generators, &chart.variables).ok_or_else(|| {
            format!("component on chart {chart_id} is not an affine-linear center")
        })?;
        if !force {
            let equations = singularity_equations(&chart);
            let rows = generators
                .iter()
                .map(affine_linear_row)
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| format!("center on chart {chart_id} is not affine-linear"))?;
            if !affine_linear_subspace_annihilates(&equations, &rows) {
                return Err(format!(
                    "center is not contained in chart {chart_id}'s singular locus; use --force to override"
                ));
            }
        }

        let center_indices = change
            .indexed_assignments
            .keys()
            .copied()
            .collect::<Vec<_>>();
        let mut output_charts = Vec::new();
        let mut new_charts = Vec::new();
        for exceptional_index in center_indices.iter().copied() {
            let blowup =
                blowup_substitutions(&chart, &change.indexed_assignments, exceptional_index);
            let composed_replacements = change
                .replacements
                .iter()
                .map(|replacement| replacement.substitute(&blowup.polys))
                .collect::<Vec<_>>();
            let transformed = chart.polynomial.substitute(&composed_replacements);
            let exceptional_power = transform.exceptional_division_power(multiplicity);
            let transformed_polynomial =
                transformed.divide_by_var_power(exceptional_index, exceptional_power)?;
            let projective_coordinates =
                substitute_polys(&chart.projective_coordinates, &composed_replacements);
            let mut blowup_coordinates =
                substitute_blowup_coordinates(&chart.blowup_coordinates, &composed_replacements);
            blowup_coordinates.push(BlowupCoordinateChart {
                stage,
                coordinates: exceptional_projective_coordinates(
                    &chart.variables,
                    &center_indices,
                    exceptional_index,
                ),
            });
            let new_id = self.next_chart_id;
            self.next_chart_id += 1;
            output_charts.push(new_id);
            let mut display = change.display.clone();
            display.extend(blowup.display);
            new_charts.push(Chart {
                id: new_id,
                parent: Some(chart.id),
                label: format!(
                    "blowup {} of chart {}, {}-chart",
                    stage, chart.id, chart.variables[exceptional_index]
                ),
                variables: chart.variables.clone(),
                polynomial: transformed_polynomial,
                active: true,
                substitutions: display,
                projective_coordinates,
                blowup_coordinates,
            });
        }

        self.charts[chart_position].active = false;
        self.charts.extend(new_charts);
        self.blowups.push(BlowupRecord {
            stage,
            input_chart: chart.id,
            center: CoordinateCenter {
                assignments: change.assignments,
            },
            multiplicity,
            output_charts,
        });
        self.cached_singular_analysis = None;
        Ok(())
    }

    fn blowup_prepared_coordinate_center(
        &mut self,
        chart_position: usize,
        chart: Chart,
        assignments: BTreeMap<String, Rational>,
        indexed_assignments: BTreeMap<usize, Rational>,
        force: bool,
        prefix_substitutions: Vec<String>,
        transform: BlowupTransform,
        stage: usize,
        known_multiplicity: Option<usize>,
    ) -> Result<(), String> {
        if !force {
            validate_center_in_singular_locus(&chart, &indexed_assignments)?;
        }

        let center_indices = indexed_assignments.keys().copied().collect::<Vec<_>>();
        let translated = chart.polynomial.translated_by(&indexed_assignments);
        let multiplicity = known_multiplicity
            .unwrap_or_else(|| translated.center_order(&center_indices).unwrap_or(0));
        if multiplicity == 0 && !force {
            return Err("center is not contained in the hypersurface".to_string());
        }

        let mut output_charts = Vec::new();
        let mut new_charts = Vec::new();
        for exceptional_index in center_indices.iter().copied() {
            let substitutions =
                blowup_substitutions(&chart, &indexed_assignments, exceptional_index);
            let transformed = substitute_zero_centered_coordinate_blowup(
                &translated,
                &center_indices,
                exceptional_index,
            );
            let exceptional_power = transform.exceptional_division_power(multiplicity);
            let transformed_polynomial =
                transformed.divide_by_var_power(exceptional_index, exceptional_power)?;
            let projective_coordinates =
                substitute_polys(&chart.projective_coordinates, &substitutions.polys);
            let mut blowup_coordinates =
                substitute_blowup_coordinates(&chart.blowup_coordinates, &substitutions.polys);
            blowup_coordinates.push(BlowupCoordinateChart {
                stage,
                coordinates: exceptional_projective_coordinates(
                    &chart.variables,
                    &center_indices,
                    exceptional_index,
                ),
            });
            let new_id = self.next_chart_id;
            self.next_chart_id += 1;
            output_charts.push(new_id);
            let mut display = prefix_substitutions.clone();
            display.extend(substitutions.display);
            new_charts.push(Chart {
                id: new_id,
                parent: Some(chart.id),
                label: format!(
                    "blowup {} of chart {}, {}-chart",
                    stage, chart.id, chart.variables[exceptional_index]
                ),
                variables: chart.variables.clone(),
                polynomial: transformed_polynomial,
                active: true,
                substitutions: display,
                projective_coordinates,
                blowup_coordinates,
            });
        }

        self.charts[chart_position].active = false;
        self.charts.extend(new_charts);
        self.blowups.push(BlowupRecord {
            stage,
            input_chart: chart.id,
            center: CoordinateCenter { assignments },
            multiplicity,
            output_charts,
        });
        self.cached_singular_analysis = None;

        Ok(())
    }

    fn next_blowup_stage(&self) -> usize {
        self.blowups
            .iter()
            .map(|record| record.stage)
            .max()
            .unwrap_or(0)
            + 1
    }

    pub fn chart(&self, chart_id: usize) -> Option<&Chart> {
        self.charts.iter().find(|chart| chart.id == chart_id)
    }

    fn component_key(
        &self,
        chart: &Chart,
        ideal_generators: &[Poly],
        dimension: usize,
        coordinate_assignments: &Option<BTreeMap<String, Rational>>,
    ) -> String {
        if dimension == 1 {
            if let Some(key) =
                projective_line_image_key(&chart.projective_coordinates, ideal_generators)
            {
                return format!("projective-line:{key}");
            }
            if let Some(key) = exceptional_line_image_key(chart, ideal_generators) {
                return key;
            }
        }

        if dimension == 0 {
            if let Some(assignments) = coordinate_assignments {
                if let Some(key) = self.global_point_key(chart, assignments) {
                    return key;
                }
            }
        }

        if chart.parent.is_none() && chart.id < self.projective_vars.len() {
            if let Some(projective_key) =
                projective_homogenized_key(ideal_generators, chart, &self.projective_vars)
            {
                return format!("projective:{projective_key}");
            }
        }

        format!("chart:{}:{}", chart.id, groebner_key(ideal_generators))
    }

    fn global_point_key(
        &self,
        chart: &Chart,
        assignments: &BTreeMap<String, Rational>,
    ) -> Option<String> {
        let mut values = Vec::with_capacity(chart.variables.len());
        for variable in &chart.variables {
            values.push(assignments.get(variable)?.clone());
        }

        let mut key = format!(
            "projective-point:{}",
            projective_point_key_from_polys(&chart.projective_coordinates, &values)?
        );
        for blowup_coordinate in &chart.blowup_coordinates {
            if let Some(direction_key) =
                projective_point_key_from_polys(&blowup_coordinate.coordinates, &values)
            {
                key.push_str(&format!(
                    ":exceptional:{}:{direction_key}",
                    blowup_coordinate.stage
                ));
            }
        }
        Some(key)
    }

    fn projective_linear_singular_components(
        &self,
        options: SingularComputationOptions,
    ) -> Option<Vec<SingularComponent>> {
        if !self.blowups.is_empty() || self.projective_vars.len() != 4 {
            return None;
        }
        if self
            .charts
            .iter()
            .any(|chart| chart.parent.is_some() || !chart.active)
        {
            return None;
        }

        let projective_components =
            find_projective_linear_singular_components(&self.initial_polynomial)?;
        if projective_components.is_empty() {
            return None;
        }

        let mut components = Vec::new();
        for projective_component in projective_components {
            let mut charts = Vec::new();
            let mut component_multiplicity = 0usize;
            for chart in self.active_charts() {
                let Some(ideal_generators) =
                    dehomogenize_projective_linear_ideal(&projective_component.generators, chart)
                else {
                    continue;
                };
                if ideal_generators.iter().any(Poly::is_nonzero_constant) {
                    continue;
                }

                let dimension =
                    ideal_dimension(&ideal_generators, chart.variables.len()).unwrap_or_default();
                if dimension != projective_component.dimension {
                    continue;
                }

                let linear_change = linear_center_change(&ideal_generators, &chart.variables);
                let multiplicity = generic_multiplicity_on_component(
                    chart,
                    &ideal_generators,
                    linear_change.as_ref(),
                );
                component_multiplicity = component_multiplicity.max(multiplicity);
                let coordinate_assignments =
                    coordinate_assignments_from_basis(&ideal_generators, &chart.variables);
                let singularity_type = if options.include_singularity_types {
                    coordinate_assignments
                        .as_ref()
                        .filter(|assignments| assignments.len() == chart.variables.len())
                        .map(|assignments| analyze_chart_point_singularity(chart, assignments))
                } else {
                    None
                };

                charts.push(SingularComponentChart {
                    chart_id: chart.id,
                    chart_label: chart.label.clone(),
                    variables: chart.variables.clone(),
                    ideal_generators,
                    dimension,
                    multiplicity,
                    coordinate_assignments,
                    affine_linear_center: linear_change.is_some(),
                    singularity_type,
                });
            }

            if charts.is_empty() {
                continue;
            }

            let mut component = SingularComponent {
                dimension: projective_component.dimension,
                multiplicity: component_multiplicity,
                generic_singularity_type: None,
                generic_singularity_chart: None,
                generic_singularity_point: None,
                charts,
            };
            if options.include_singularity_types && component.dimension == 1 {
                if let Some((chart_id, point, singularity_type)) =
                    generic_line_singularity_type(&component, &self.charts)
                {
                    component.generic_singularity_type = Some(singularity_type);
                    component.generic_singularity_chart = Some(chart_id);
                    component.generic_singularity_point = Some(point);
                }
            }
            components.push(component);
        }

        if components.is_empty() {
            return None;
        }

        components.sort_by(|left, right| {
            right
                .dimension
                .cmp(&left.dimension)
                .then_with(|| left.charts[0].chart_id.cmp(&right.charts[0].chart_id))
                .then_with(|| component_display_key(left).cmp(&component_display_key(right)))
        });
        Some(components)
    }

    fn choose_crepant_center(
        &self,
        analysis: &SingularAnalysis,
    ) -> Option<AutoResolutionCandidate> {
        let mut candidates = Vec::new();

        for (index, component) in analysis.components.iter().enumerate() {
            if let Some(candidate) = self.crepant_component_candidate(index + 1, component) {
                candidates.push(candidate);
            }
        }

        for (index, intersection) in analysis.intersections.iter().enumerate() {
            if let Some(candidate) = self.crepant_intersection_candidate(index + 1, intersection) {
                candidates.push(candidate);
            }
        }

        candidates.sort_by_key(|candidate| candidate.score);
        candidates.into_iter().next()
    }

    fn crepant_component_candidate(
        &self,
        display_index: usize,
        component: &SingularComponent,
    ) -> Option<AutoResolutionCandidate> {
        let mut ambient_codimension = None::<usize>;
        let mut multiplicity = None::<usize>;
        let mut input_chart_ids = Vec::new();

        for chart_component in &component.charts {
            if !component_chart_has_supported_center(chart_component) {
                return None;
            }

            let local_codimension = chart_component
                .variables
                .len()
                .checked_sub(chart_component.dimension)?;
            if local_codimension < 2 {
                return None;
            }
            if !is_crepant_double_cover_branch_blowup(
                local_codimension,
                chart_component.multiplicity,
            ) {
                return None;
            }

            if ambient_codimension.is_some_and(|codimension| codimension != local_codimension) {
                return None;
            }
            if multiplicity.is_some_and(|value| value != chart_component.multiplicity) {
                return None;
            }

            ambient_codimension = Some(local_codimension);
            multiplicity = Some(chart_component.multiplicity);
            input_chart_ids.push(chart_component.chart_id);
        }

        let point_priority = usize::from(component.dimension == 0);
        let kind_priority = if component.dimension == 0 { 0 } else { 2 };
        Some(AutoResolutionCandidate {
            center: AutoResolutionCenter::Component {
                index: display_index,
            },
            center_dimension: component.dimension,
            ambient_codimension: ambient_codimension?,
            multiplicity: multiplicity?,
            input_chart_ids,
            score: (
                point_priority,
                kind_priority,
                component.multiplicity,
                display_index,
            ),
        })
    }

    fn crepant_intersection_candidate(
        &self,
        display_index: usize,
        intersection: &SingularIntersection,
    ) -> Option<AutoResolutionCandidate> {
        let assignments = intersection.coordinate_assignments.as_ref()?;
        let chart = self.chart(intersection.chart_id)?;
        if assignments.len() != chart.variables.len() || assignments.len() < 2 {
            return None;
        }

        let ambient_codimension = chart.variables.len();
        let singularity_type = analyze_chart_point_singularity(chart, assignments);
        if !is_crepant_double_cover_branch_blowup(
            ambient_codimension,
            singularity_type.multiplicity,
        ) {
            return None;
        }

        Some(AutoResolutionCandidate {
            center: AutoResolutionCenter::Intersection {
                index: display_index,
            },
            center_dimension: 0,
            ambient_codimension,
            multiplicity: singularity_type.multiplicity,
            input_chart_ids: vec![intersection.chart_id],
            score: (1, 1, singularity_type.multiplicity, display_index),
        })
    }
}

pub fn analyze_affine_polynomial_at(
    polynomial_text: &str,
    variables: &[String],
    values: &[Rational],
) -> Result<SingularityType, String> {
    if variables.len() != values.len() {
        return Err(format!(
            "expected {} point coordinates, received {}",
            variables.len(),
            values.len()
        ));
    }
    let polynomial = parse_polynomial(polynomial_text, variables)?;
    let assignments = values.iter().cloned().enumerate().collect();
    let maximum_degree = formal_series_degree_bound(&polynomial);
    let mut jet_degree = maximum_degree.min(8);
    loop {
        let translated = polynomial.translated_jet_by(&assignments, jet_degree);
        let singularity = analyze_transverse_singularity_with_bound(&translated, Some(jet_degree));
        let unresolved_corank_one = singularity.multiplicity == 2
            && singularity.quadratic_rank == Some(2)
            && singularity.embedding_dimension == 3
            && !singularity.label.starts_with("cA");
        if !unresolved_corank_one || jet_degree == maximum_degree {
            return Ok(singularity);
        }
        jet_degree = (jet_degree * 2).min(maximum_degree);
    }
}

fn compute_local_singular_components_in_parallel(
    charts: Vec<&Chart>,
    positive_dimensional_only: bool,
) -> Vec<(usize, Vec<Vec<Poly>>)> {
    if charts.len() <= 1 {
        return charts
            .into_iter()
            .map(|chart| {
                let equations = singularity_equations(chart);
                (
                    chart.id,
                    local_singular_components(chart, &equations, positive_dimensional_only),
                )
            })
            .collect();
    }

    let worker_count = 8.min(charts.len());
    let chunk_size = charts.len().div_ceil(worker_count);
    let mut results = std::thread::scope(|scope| {
        let handles = charts
            .chunks(chunk_size)
            .map(|chunk| {
                scope.spawn(move || {
                    chunk
                        .iter()
                        .map(|chart| {
                            let equations = singularity_equations(chart);
                            (
                                chart.id,
                                local_singular_components(
                                    chart,
                                    &equations,
                                    positive_dimensional_only,
                                ),
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .flat_map(|handle| handle.join().expect("singular chart worker panicked"))
            .collect::<Vec<_>>()
    });
    results.sort_by_key(|(chart_id, _)| *chart_id);
    results
}

impl AutoResolutionCandidate {
    fn into_step(self, records: &[BlowupRecord]) -> AutoResolutionStep {
        let blowup_stages = records
            .iter()
            .map(|record| record.stage)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        AutoResolutionStep {
            center: self.center,
            center_dimension: self.center_dimension,
            ambient_codimension: self.ambient_codimension,
            multiplicity: self.multiplicity,
            input_chart_ids: self.input_chart_ids,
            blowup_stages,
            output_chart_ids: records
                .iter()
                .flat_map(|record| record.output_charts.iter().copied())
                .collect(),
        }
    }
}

fn component_chart_has_supported_center(chart_component: &SingularComponentChart) -> bool {
    chart_component
        .coordinate_assignments
        .as_ref()
        .is_some_and(|assignments| assignments.len() >= 2)
        || chart_component.affine_linear_center
}

fn is_crepant_double_cover_branch_blowup(ambient_codimension: usize, multiplicity: usize) -> bool {
    // For a double cover branched along f = 0, blowing up a center of
    // codimension c has discrepancy c - 1 - floor(m / 2).
    ambient_codimension >= 2 && multiplicity / 2 == ambient_codimension - 1
}

fn double_cover_branch_division_power(multiplicity: usize) -> usize {
    multiplicity - (multiplicity % 2)
}

fn substitute_polys(polys: &[Poly], replacements: &[Poly]) -> Vec<Poly> {
    polys
        .iter()
        .map(|poly| poly.substitute(replacements))
        .collect()
}

fn substitute_blowup_coordinates(
    coordinates: &[BlowupCoordinateChart],
    replacements: &[Poly],
) -> Vec<BlowupCoordinateChart> {
    coordinates
        .iter()
        .map(|coordinate_chart| BlowupCoordinateChart {
            stage: coordinate_chart.stage,
            coordinates: substitute_polys(&coordinate_chart.coordinates, replacements),
        })
        .collect()
}

fn exceptional_projective_coordinates(
    variables: &[String],
    center_indices: &[usize],
    exceptional_index: usize,
) -> Vec<Poly> {
    center_indices
        .iter()
        .map(|index| {
            if *index == exceptional_index {
                Poly::constant(variables, Rational::one())
            } else {
                Poly::var(variables, *index)
            }
        })
        .collect()
}

fn projective_line_image_key(coordinate_polys: &[Poly], generators: &[Poly]) -> Option<String> {
    let points = projective_image_points(coordinate_polys, generators)?;
    projective_linear_key_from_points(&points)
}

fn exceptional_line_image_key(chart: &Chart, generators: &[Poly]) -> Option<String> {
    let base_point = projective_point_image_key(&chart.projective_coordinates, generators)?;
    let mut prefix = format!("projective-point:{base_point}");

    for blowup_coordinate in &chart.blowup_coordinates {
        if let Some(line_key) =
            projective_line_image_key(&blowup_coordinate.coordinates, generators)
        {
            return Some(format!(
                "exceptional-line:{prefix}:stage:{}:{line_key}",
                blowup_coordinate.stage
            ));
        }

        if let Some(point_key) =
            projective_point_image_key(&blowup_coordinate.coordinates, generators)
        {
            prefix.push_str(&format!(
                ":exceptional:{}:{point_key}",
                blowup_coordinate.stage
            ));
        }
    }

    None
}

fn projective_point_image_key(coordinate_polys: &[Poly], generators: &[Poly]) -> Option<String> {
    let points = projective_image_points(coordinate_polys, generators)?;
    (points.len() == 1).then(|| projective_point_key_from_normalized(&points[0]))
}

fn projective_image_points(
    coordinate_polys: &[Poly],
    generators: &[Poly],
) -> Option<Vec<Vec<Rational>>> {
    if coordinate_polys.is_empty() {
        return None;
    }

    let samples = affine_component_sample_points(generators)?;
    let mut unique = BTreeSet::new();
    for sample in samples {
        let coordinates = evaluate_polys(coordinate_polys, &sample)?;
        if let Some(normalized) = normalized_projective_coordinates(coordinates) {
            unique.insert(normalized);
        }
    }

    if unique.is_empty() {
        None
    } else {
        Some(unique.into_iter().collect())
    }
}

fn affine_component_sample_points(generators: &[Poly]) -> Option<Vec<Vec<Rational>>> {
    let variable_count = generators.first()?.vars.len();
    let rows = generators
        .iter()
        .map(affine_linear_row)
        .collect::<Option<Vec<_>>>()?;
    let reduced_rows = rref_affine_rows(rows, variable_count)?;
    let pivot_indices = reduced_rows
        .iter()
        .map(|(pivot, _)| *pivot)
        .collect::<Vec<_>>();
    let free_indices = (0..variable_count)
        .filter(|index| !pivot_indices.contains(index))
        .collect::<Vec<_>>();

    let mut samples = Vec::new();
    for parameters in affine_parameter_samples(free_indices.len()) {
        let mut point = vec![Rational::zero(); variable_count];
        for (parameter_index, variable_index) in free_indices.iter().enumerate() {
            point[*variable_index] = parameters[parameter_index].clone();
        }
        for (pivot, row) in &reduced_rows {
            let mut value = -row[variable_count].clone();
            for variable_index in &free_indices {
                if row[*variable_index].is_zero() {
                    continue;
                }
                value = value - row[*variable_index].clone() * point[*variable_index].clone();
            }
            point[*pivot] = value;
        }
        samples.push(point);
    }

    Some(samples)
}

fn projective_point_key_from_polys(polys: &[Poly], values: &[Rational]) -> Option<String> {
    let coordinates = evaluate_polys(polys, values)?;
    let normalized = normalized_projective_coordinates(coordinates)?;
    Some(projective_point_key_from_normalized(&normalized))
}

fn evaluate_polys(polys: &[Poly], values: &[Rational]) -> Option<Vec<Rational>> {
    polys
        .iter()
        .map(|poly| poly.evaluate(values).ok())
        .collect()
}

fn normalized_projective_coordinates(coordinates: Vec<Rational>) -> Option<Vec<Rational>> {
    let first_nonzero = coordinates.iter().find(|value| !value.is_zero())?.clone();
    Some(
        coordinates
            .into_iter()
            .map(|value| value / first_nonzero.clone())
            .collect(),
    )
}

fn projective_point_key_from_normalized(coordinates: &[Rational]) -> String {
    coordinates
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn projective_linear_key_from_points(points: &[Vec<Rational>]) -> Option<String> {
    let coordinate_count = points.first()?.len();
    if points.iter().any(|point| point.len() != coordinate_count) {
        return None;
    }

    let rank = rref_homogeneous_rows(points.to_vec(), coordinate_count).len();
    if rank != 2 {
        return None;
    }
    if coordinate_count == 2 {
        return Some("P1".to_string());
    }

    let annihilator_rows = homogeneous_annihilator_rows(points, coordinate_count)?;
    if annihilator_rows.is_empty() {
        return None;
    }
    Some(homogeneous_rows_key(&annihilator_rows))
}

fn homogeneous_annihilator_rows(
    points: &[Vec<Rational>],
    coordinate_count: usize,
) -> Option<Vec<Vec<Rational>>> {
    if points.iter().any(|point| point.len() != coordinate_count) {
        return None;
    }

    let rref = rref_homogeneous_rows(points.to_vec(), coordinate_count);
    let pivots = rref.iter().map(|(pivot, _)| *pivot).collect::<Vec<_>>();
    let free_indices = (0..coordinate_count)
        .filter(|index| !pivots.contains(index))
        .collect::<Vec<_>>();

    let mut rows = Vec::new();
    for free_index in free_indices {
        let mut row = vec![Rational::zero(); coordinate_count];
        row[free_index] = Rational::one();
        for (pivot, pivot_row) in &rref {
            row[*pivot] = -pivot_row[free_index].clone();
        }
        rows.push(row);
    }

    Some(
        rref_homogeneous_rows(rows, coordinate_count)
            .into_iter()
            .map(|(_, row)| row)
            .collect(),
    )
}

fn rref_homogeneous_rows(
    mut rows: Vec<Vec<Rational>>,
    column_count: usize,
) -> Vec<(usize, Vec<Rational>)> {
    rows.retain(|row| row.len() == column_count && row.iter().any(|value| !value.is_zero()));

    let mut pivot_row = 0usize;
    for column in 0..column_count {
        let Some(source_row) = (pivot_row..rows.len()).find(|row| !rows[*row][column].is_zero())
        else {
            continue;
        };

        rows.swap(pivot_row, source_row);
        let pivot = rows[pivot_row][column].clone();
        for value in &mut rows[pivot_row] {
            *value = value.clone() / pivot.clone();
        }

        let normalized_pivot = rows[pivot_row].clone();
        for row in 0..rows.len() {
            if row == pivot_row || rows[row][column].is_zero() {
                continue;
            }
            let factor = rows[row][column].clone();
            for col in 0..column_count {
                rows[row][col] =
                    rows[row][col].clone() - factor.clone() * normalized_pivot[col].clone();
            }
        }

        pivot_row += 1;
        if pivot_row == rows.len() {
            break;
        }
    }

    let mut reduced_rows = rows
        .into_iter()
        .filter_map(|row| {
            row.iter()
                .position(|value| !value.is_zero())
                .map(|pivot| (pivot, row))
        })
        .collect::<Vec<_>>();
    reduced_rows.sort_by_key(|(pivot, _)| *pivot);
    reduced_rows
}

fn homogeneous_rows_key(rows: &[Vec<Rational>]) -> String {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>()
        .join(";")
}

struct BlowupSubstitutions {
    polys: Vec<Poly>,
    display: Vec<String>,
}

struct LinearCenterChange {
    replacements: Vec<Poly>,
    assignments: BTreeMap<String, Rational>,
    indexed_assignments: BTreeMap<usize, Rational>,
    display: Vec<String>,
}

fn singularity_equations(chart: &Chart) -> Vec<Poly> {
    let mut equations = Vec::with_capacity(chart.variables.len() + 1);
    equations.push(chart.polynomial.clone());
    equations.extend(chart.polynomial.partials());
    equations
}

fn local_singular_components(
    chart: &Chart,
    equations: &[Poly],
    positive_dimensional_only: bool,
) -> Vec<Vec<Poly>> {
    if positive_dimensional_only || should_use_affine_linear_fast_path(chart) {
        if let Some(components) = affine_linear_singular_components(equations, chart) {
            if !components.is_empty() {
                return components;
            }
        }
    }

    if positive_dimensional_only {
        if singular_backend_enabled() {
            if let Some(components) = singular_minimal_primes(equations, &chart.variables, None) {
                return components
                    .into_iter()
                    .filter(|component| {
                        ideal_dimension(component, chart.variables.len()).is_some_and(|dim| dim > 0)
                    })
                    .collect();
            }
            eprintln!(
                "warning: Singular positive-dimensional audit failed on chart {}; treating the affine-linear search as incomplete",
                chart.id
            );
        }
        return Vec::new();
    }

    if singular_backend_enabled() {
        let point_bound = std::env::var("BLOWUP_SINGULAR_POINT_BOUND")
            .ok()
            .and_then(|value| value.parse::<i32>().ok())
            .filter(|bound| (1..=16).contains(bound))
            .unwrap_or(2);
        let (rational_points, sample_limit_hit) =
            enumerate_singular_points(chart, equations, point_bound, 10_000);
        let rational_points = (!sample_limit_hit).then_some(rational_points);
        if let Some(components) =
            singular_minimal_primes(equations, &chart.variables, rational_points.as_deref())
        {
            return components;
        }
        eprintln!(
            "warning: Singular minimal-prime backend failed on chart {}; using the in-process fallback",
            chart.id
        );
    }

    decompose_ideal(equations)
}

fn singular_backend_enabled() -> bool {
    std::env::var("BLOWUP_SINGULAR_BACKEND").is_ok_and(|value| {
        matches!(
            value.trim().to_ascii_lowercase().as_str(),
            "1" | "true" | "yes" | "singular"
        )
    })
}

fn singular_minimal_primes(
    equations: &[Poly],
    variables: &[String],
    rational_points: Option<&[SingularityPoint]>,
) -> Option<Vec<Vec<Poly>>> {
    if variables.is_empty() || equations.is_empty() {
        return None;
    }

    let mut child = Command::new("Singular")
        .arg("-q")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .ok()?;
    let script = singular_minimal_primes_script(equations, variables, rational_points);
    child.stdin.as_mut()?.write_all(script.as_bytes()).ok()?;
    drop(child.stdin.take());
    let output = child.wait_with_output().ok()?;
    if !output.status.success() {
        return None;
    }
    let stdout = String::from_utf8(output.stdout).ok()?;
    parse_singular_minimal_primes(&stdout, variables)
}

fn singular_minimal_primes_script(
    equations: &[Poly],
    variables: &[String],
    rational_points: Option<&[SingularityPoint]>,
) -> String {
    let ideal = equations
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let decomposition = if let Some(points) = rational_points {
        let (point_ideal_declarations, certified_list) = if points.is_empty() {
            (String::new(), "list(ideal(1))".to_string())
        } else {
            let declarations = points
                .iter()
                .enumerate()
                .map(|(point_index, point)| {
                    let generators = variables
                        .iter()
                        .zip(point.values.iter())
                        .map(|(variable, value)| format!("{variable}-({value})"))
                        .collect::<Vec<_>>()
                        .join(",");
                    format!("ideal P{point_index}={generators};")
                })
                .collect::<Vec<_>>()
                .join("\n");
            let names = (0..points.len())
                .map(|index| format!("P{index}"))
                .collect::<Vec<_>>()
                .join(",");
            (declarations, format!("list({names})"))
        };
        format!(
            "{point_ideal_declarations}\n\
             ideal G=std(I);\n\
             int D=vdim(G);\n\
             list L;\n\
             if (D=={}) {{ L={certified_list}; }} else {{ L=minAss(I); }}",
            points.len()
        )
    } else {
        "list L=minAss(I);".to_string()
    };
    format!(
        "LIB \"primdec.lib\";\n\
         ring r=0,({}),dp;\n\
         ideal I={ideal};\n\
         {decomposition}\n\
         print(\"BLOWUP_COUNT=\"+string(size(L)));\n\
         int i,j;\n\
         for (i=1;i<=size(L);i++) {{\n\
           ideal J=L[i];\n\
           print(\"BLOWUP_BEGIN=\"+string(i));\n\
           for (j=1;j<=size(J);j++) {{ print(\"BLOWUP_GEN=\"+string(J[j])); }}\n\
           print(\"BLOWUP_END=\"+string(i));\n\
         }}\n\
         quit;\n",
        variables.join(",")
    )
}

fn parse_singular_minimal_primes(output: &str, variables: &[String]) -> Option<Vec<Vec<Poly>>> {
    let count = output
        .lines()
        .find_map(|line| line.strip_prefix("BLOWUP_COUNT="))?
        .trim()
        .parse::<usize>()
        .ok()?;
    let mut components = Vec::with_capacity(count);
    let mut current = None::<Vec<Poly>>;

    for line in output.lines() {
        if line.starts_with("BLOWUP_BEGIN=") {
            if current.is_some() {
                return None;
            }
            current = Some(Vec::new());
        } else if let Some(generator) = line.strip_prefix("BLOWUP_GEN=") {
            let component = current.as_mut()?;
            component.push(parse_polynomial(generator.trim(), variables).ok()?);
        } else if line.starts_with("BLOWUP_END=") {
            let component = current.take()?;
            if component.is_empty() {
                return None;
            }
            components.push(component);
        }
    }

    (current.is_none() && components.len() == count).then_some(components)
}

fn should_use_affine_linear_fast_path(chart: &Chart) -> bool {
    // Large charts in these examples often have affine-linear singular centers.
    // Try that exact search before falling back to Groebner decomposition.
    chart.polynomial.total_degree() >= 8 || chart.polynomial.terms.len() >= 32
}

fn affine_linear_singular_components(equations: &[Poly], chart: &Chart) -> Option<Vec<Vec<Poly>>> {
    let variables = &chart.variables;
    if !(2..=3).contains(&variables.len()) {
        return Some(Vec::new());
    }

    let coefficient_values = small_projective_search_coefficients();
    let mut components = tangent_reconstructed_affine_lines(chart, equations);
    if !components.is_empty() && sampled_singular_locus_is_covered(chart, equations, &components) {
        return Some(remove_redundant_linear_components(components));
    }
    let mut seen = BTreeSet::new();
    for component in &components {
        seen.insert(groebner_key(component));
    }

    let maximum_codimension = if variables.len() == 3 {
        variables.len() - 1
    } else {
        variables.len()
    };
    for codimension in 2..=maximum_codimension {
        for rows in rref_affine_linear_rows(variables.len(), codimension, &coefficient_values) {
            if !affine_linear_subspace_annihilates(equations, &rows) {
                continue;
            }

            let generators = rows
                .iter()
                .map(|row| affine_row_to_poly(row, variables).monic())
                .collect::<Vec<_>>();
            let Some(basis) = affine_linear_basis(&generators, variables) else {
                continue;
            };
            if basis.is_empty() {
                continue;
            }

            let key = groebner_key(&basis);
            if seen.insert(key) {
                components.push(basis);
                if components.len() > 256 {
                    return None;
                }
            }
        }
    }

    Some(remove_redundant_linear_components(components))
}

fn sampled_singular_locus_is_covered(
    chart: &Chart,
    equations: &[Poly],
    components: &[Vec<Poly>],
) -> bool {
    let (points, _) = enumerate_singular_points(chart, equations, 2, 1_000);
    !points.is_empty()
        && points.iter().all(|point| {
            components.iter().any(|component| {
                component.iter().all(|generator| {
                    generator
                        .evaluate(&point.values)
                        .is_ok_and(|value| value.is_zero())
                })
            })
        })
}

fn tangent_reconstructed_affine_lines(chart: &Chart, equations: &[Poly]) -> Vec<Vec<Poly>> {
    if chart.variables.len() != 3 {
        return Vec::new();
    }

    let (points, _) = enumerate_singular_points(chart, equations, 2, 1_000);
    let mut components = BTreeMap::new();
    for point in points {
        let jacobian_rows = equations
            .iter()
            .map(|equation| {
                equation
                    .partials()
                    .iter()
                    .map(|partial| partial.evaluate(&point.values).ok())
                    .collect::<Option<Vec<_>>>()
            })
            .collect::<Option<Vec<_>>>();
        let Some(jacobian_rows) = jacobian_rows else {
            continue;
        };
        let reduced_jacobian = rref_homogeneous_rows(jacobian_rows, chart.variables.len());
        if reduced_jacobian.len() != 2 {
            continue;
        }
        let jacobian_basis = reduced_jacobian
            .into_iter()
            .map(|(_, row)| row)
            .collect::<Vec<_>>();
        let Some(tangent_directions) =
            homogeneous_annihilator_rows(&jacobian_basis, chart.variables.len())
        else {
            continue;
        };
        let [direction] = tangent_directions.as_slice() else {
            continue;
        };
        let Some(normal_rows) =
            homogeneous_annihilator_rows(std::slice::from_ref(direction), chart.variables.len())
        else {
            continue;
        };

        let affine_rows = normal_rows
            .into_iter()
            .map(|mut row| {
                let constant = -row
                    .iter()
                    .zip(point.values.iter())
                    .fold(Rational::zero(), |sum, (coefficient, value)| {
                        sum + coefficient.clone() * value.clone()
                    });
                row.push(constant);
                row
            })
            .collect::<Vec<_>>();
        if !affine_linear_subspace_annihilates(equations, &affine_rows) {
            continue;
        }
        let generators = affine_rows
            .iter()
            .map(|row| affine_row_to_poly(row, &chart.variables).monic())
            .collect::<Vec<_>>();
        let Some(basis) = affine_linear_basis(&generators, &chart.variables) else {
            continue;
        };
        components.entry(groebner_key(&basis)).or_insert(basis);
    }
    components.into_values().collect()
}

fn rref_affine_linear_rows(
    variable_count: usize,
    codimension: usize,
    coefficient_values: &[Rational],
) -> Vec<Vec<Vec<Rational>>> {
    if codimension == 0 || codimension > variable_count {
        return Vec::new();
    }

    let pivot_sets = combinations(variable_count, codimension);
    let mut row_sets = Vec::new();
    for pivots in pivot_sets {
        let free_indices = (0..variable_count)
            .filter(|index| !pivots.contains(index))
            .collect::<Vec<_>>();
        let slots = pivots
            .iter()
            .enumerate()
            .flat_map(|(row_index, pivot)| {
                free_indices
                    .iter()
                    .copied()
                    .filter(move |free_index| free_index > pivot)
                    .chain(std::iter::once(variable_count))
                    .map(move |column| (row_index, column))
            })
            .collect::<Vec<_>>();
        let mut rows = vec![vec![Rational::zero(); variable_count + 1]; codimension];
        for (row_index, pivot) in pivots.iter().enumerate() {
            rows[row_index][*pivot] = Rational::one();
        }
        fill_affine_row_coefficients(0, &slots, coefficient_values, &mut rows, &mut row_sets);
    }
    row_sets
}

fn fill_affine_row_coefficients(
    slot_index: usize,
    slots: &[(usize, usize)],
    coefficient_values: &[Rational],
    rows: &mut Vec<Vec<Rational>>,
    row_sets: &mut Vec<Vec<Vec<Rational>>>,
) {
    if slot_index == slots.len() {
        row_sets.push(rows.clone());
        return;
    }

    let (row_index, column) = slots[slot_index];
    for value in coefficient_values {
        rows[row_index][column] = value.clone();
        fill_affine_row_coefficients(slot_index + 1, slots, coefficient_values, rows, row_sets);
    }
    rows[row_index][column] = Rational::zero();
}

fn affine_linear_subspace_annihilates(equations: &[Poly], rows: &[Vec<Rational>]) -> bool {
    let variable_count = equations.first().map_or(0, |poly| poly.vars.len());
    let Some((parameter_vars, replacements)) =
        affine_subspace_parameterization(variable_count, rows)
    else {
        return false;
    };

    if !affine_subspace_sample_checks(equations, &replacements, parameter_vars.len()) {
        return false;
    }

    if parameter_vars.len() == 1 {
        let degree_bound = equations.iter().map(Poly::total_degree).max().unwrap_or(0);
        return (0..=degree_bound).all(|value| {
            let parameter = [Rational::from_i128(value as i128)];
            let Some(point) = replacements
                .iter()
                .map(|replacement| replacement.evaluate(&parameter).ok())
                .collect::<Option<Vec<_>>>()
            else {
                return false;
            };
            equations.iter().all(|equation| {
                equation
                    .evaluate(&point)
                    .is_ok_and(|result| result.is_zero())
            })
        });
    }

    equations
        .iter()
        .all(|equation| equation.substitute(&replacements).is_zero())
}

fn affine_subspace_parameterization(
    variable_count: usize,
    rows: &[Vec<Rational>],
) -> Option<(Vec<String>, Vec<Poly>)> {
    let pivots = rows
        .iter()
        .map(|row| {
            row.iter()
                .take(variable_count)
                .position(|value| !value.is_zero())
        })
        .collect::<Option<Vec<_>>>()?;
    let free_indices = (0..variable_count)
        .filter(|index| !pivots.contains(index))
        .collect::<Vec<_>>();

    let parameter_vars = (0..free_indices.len())
        .map(|index| format!("t{index}"))
        .collect::<Vec<_>>();
    let mut replacements = vec![Poly::zero(&parameter_vars); variable_count];
    for (parameter_index, variable_index) in free_indices.iter().enumerate() {
        replacements[*variable_index] = Poly::var(&parameter_vars, parameter_index);
    }
    for (row, pivot) in rows.iter().zip(pivots.iter()) {
        let mut replacement = Poly::constant(&parameter_vars, -row[variable_count].clone());
        for (parameter_index, variable_index) in free_indices.iter().enumerate() {
            if row[*variable_index].is_zero() {
                continue;
            }
            replacement = replacement.sub(
                &Poly::var(&parameter_vars, parameter_index).scale(row[*variable_index].clone()),
            );
        }
        replacements[*pivot] = replacement;
    }

    Some((parameter_vars, replacements))
}

fn affine_subspace_sample_checks(
    equations: &[Poly],
    replacements: &[Poly],
    parameter_count: usize,
) -> bool {
    affine_parameter_samples(parameter_count)
        .into_iter()
        .all(|sample| {
            let Some(point) = replacements
                .iter()
                .map(|replacement| replacement.evaluate(&sample).ok())
                .collect::<Option<Vec<_>>>()
            else {
                return false;
            };
            equations
                .iter()
                .all(|equation| equation.evaluate(&point).is_ok_and(|value| value.is_zero()))
        })
}

fn affine_parameter_samples(parameter_count: usize) -> Vec<Vec<Rational>> {
    match parameter_count {
        0 => vec![Vec::new()],
        1 => [-2, -1, 0, 1, 2]
            .into_iter()
            .map(|value| vec![Rational::from_i128(value)])
            .collect(),
        count => projective_parameter_samples(count)
            .into_iter()
            .chain(std::iter::once(vec![Rational::zero(); count]))
            .collect(),
    }
}

fn remove_redundant_linear_components(components: Vec<Vec<Poly>>) -> Vec<Vec<Poly>> {
    let mut unique = BTreeMap::new();
    for component in components {
        unique.entry(groebner_key(&component)).or_insert(component);
    }
    let components = unique.into_values().collect::<Vec<_>>();

    let mut keep = vec![true; components.len()];
    for left in 0..components.len() {
        for right in 0..components.len() {
            if left == right {
                continue;
            }
            let left_contains_right = ideal_contains_all(&components[left], &components[right]);
            let right_contains_left = ideal_contains_all(&components[right], &components[left]);
            if left_contains_right && !right_contains_left {
                keep[left] = false;
                break;
            }
        }
    }

    let mut minimal = components
        .into_iter()
        .enumerate()
        .filter_map(|(index, component)| keep[index].then_some(component))
        .collect::<Vec<_>>();
    minimal.sort_by(|left, right| {
        let left_dim = ideal_dimension(left, left.first().map_or(0, |poly| poly.vars.len()));
        let right_dim = ideal_dimension(right, right.first().map_or(0, |poly| poly.vars.len()));
        right_dim
            .cmp(&left_dim)
            .then_with(|| groebner_key(left).cmp(&groebner_key(right)))
    });
    minimal
}

fn find_projective_linear_singular_components(
    hypersurface: &Poly,
) -> Option<Vec<ProjectiveLinearComponent>> {
    if hypersurface.vars.len() != 4 {
        return None;
    }

    let equations = hypersurface
        .partials()
        .into_iter()
        .filter(|poly| !poly.is_zero())
        .collect::<Vec<_>>();
    if equations.is_empty() {
        return None;
    }

    let lines = contained_projective_linear_subspaces(&equations, &hypersurface.vars, 2, 64)?;
    let points = contained_projective_linear_subspaces(&equations, &hypersurface.vars, 3, 128)?;
    let mut components = lines
        .into_iter()
        .map(|generators| ProjectiveLinearComponent {
            generators,
            dimension: 1,
        })
        .collect::<Vec<_>>();

    for point_generators in points {
        let Some(point) = projective_point_coordinates(&point_generators) else {
            continue;
        };
        if components
            .iter()
            .any(|component| projective_generators_vanish_at(&component.generators, &point))
        {
            continue;
        }
        components.push(ProjectiveLinearComponent {
            generators: point_generators,
            dimension: 0,
        });
    }

    (!components.is_empty()).then_some(components)
}

fn contained_projective_linear_subspaces(
    equations: &[Poly],
    variables: &[String],
    codimension: usize,
    limit: usize,
) -> Option<Vec<Vec<Poly>>> {
    let coefficient_values = small_projective_search_coefficients();
    let mut components = Vec::new();
    let mut seen = BTreeSet::new();

    for rows in rref_homogeneous_linear_rows(variables.len(), codimension, &coefficient_values) {
        if !linear_subspace_annihilates(equations, &rows) {
            continue;
        }
        let generators = rows
            .iter()
            .map(|row| homogeneous_row_to_poly(row, variables).monic())
            .collect::<Vec<_>>();
        let key = groebner_key(&generators);
        if seen.insert(key) {
            components.push(generators);
            if components.len() > limit {
                return None;
            }
        }
    }

    Some(components)
}

fn small_projective_search_coefficients() -> Vec<Rational> {
    let mut values = BTreeSet::new();
    for denominator in 1..=2 {
        for numerator in -2..=2 {
            values.insert(Rational::new(numerator, denominator));
        }
    }
    values.into_iter().collect()
}

fn rref_homogeneous_linear_rows(
    variable_count: usize,
    codimension: usize,
    coefficient_values: &[Rational],
) -> Vec<Vec<Vec<Rational>>> {
    if codimension == 0 || codimension >= variable_count {
        return Vec::new();
    }

    let pivot_sets = combinations(variable_count, codimension);
    let mut row_sets = Vec::new();
    for pivots in pivot_sets {
        let free_indices = (0..variable_count)
            .filter(|index| !pivots.contains(index))
            .collect::<Vec<_>>();
        let slots = pivots
            .iter()
            .enumerate()
            .flat_map(|(row_index, pivot)| {
                free_indices
                    .iter()
                    .copied()
                    .filter(move |free_index| free_index > pivot)
                    .map(move |free_index| (row_index, free_index))
            })
            .collect::<Vec<_>>();
        let mut rows = vec![vec![Rational::zero(); variable_count]; codimension];
        for (row_index, pivot) in pivots.iter().enumerate() {
            rows[row_index][*pivot] = Rational::one();
        }
        fill_homogeneous_row_coefficients(0, &slots, coefficient_values, &mut rows, &mut row_sets);
    }
    row_sets
}

fn fill_homogeneous_row_coefficients(
    slot_index: usize,
    slots: &[(usize, usize)],
    coefficient_values: &[Rational],
    rows: &mut Vec<Vec<Rational>>,
    row_sets: &mut Vec<Vec<Vec<Rational>>>,
) {
    if slot_index == slots.len() {
        row_sets.push(rows.clone());
        return;
    }

    let (row_index, variable_index) = slots[slot_index];
    for value in coefficient_values {
        rows[row_index][variable_index] = value.clone();
        fill_homogeneous_row_coefficients(
            slot_index + 1,
            slots,
            coefficient_values,
            rows,
            row_sets,
        );
    }
    rows[row_index][variable_index] = Rational::zero();
}

fn combinations(item_count: usize, choose: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = Vec::new();
    combinations_recursive(0, item_count, choose, &mut current, &mut result);
    result
}

fn combinations_recursive(
    start: usize,
    item_count: usize,
    choose: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == choose {
        result.push(current.clone());
        return;
    }
    let remaining = choose - current.len();
    for index in start..=item_count - remaining {
        current.push(index);
        combinations_recursive(index + 1, item_count, choose, current, result);
        current.pop();
    }
}

fn linear_subspace_annihilates(equations: &[Poly], rows: &[Vec<Rational>]) -> bool {
    let variable_count = equations.first().map_or(0, |poly| poly.vars.len());
    let Some((parameter_vars, replacements)) =
        projective_subspace_parameterization(variable_count, rows)
    else {
        return false;
    };

    if parameter_vars.is_empty() {
        return false;
    }
    if !linear_subspace_sample_checks(equations, &replacements, parameter_vars.len()) {
        return false;
    }
    if parameter_vars.len() == 1 {
        return true;
    }

    equations
        .iter()
        .all(|equation| equation.substitute(&replacements).is_zero())
}

fn projective_subspace_parameterization(
    variable_count: usize,
    rows: &[Vec<Rational>],
) -> Option<(Vec<String>, Vec<Poly>)> {
    let pivots = rows
        .iter()
        .map(|row| row.iter().position(|value| !value.is_zero()))
        .collect::<Option<Vec<_>>>()?;
    let free_indices = (0..variable_count)
        .filter(|index| !pivots.contains(index))
        .collect::<Vec<_>>();
    if free_indices.is_empty() {
        return None;
    }

    let parameter_vars = (0..free_indices.len())
        .map(|index| format!("t{index}"))
        .collect::<Vec<_>>();
    let mut replacements = vec![Poly::zero(&parameter_vars); variable_count];
    for (parameter_index, variable_index) in free_indices.iter().enumerate() {
        replacements[*variable_index] = Poly::var(&parameter_vars, parameter_index);
    }
    for (row, pivot) in rows.iter().zip(pivots.iter()) {
        let mut replacement = Poly::zero(&parameter_vars);
        for (parameter_index, variable_index) in free_indices.iter().enumerate() {
            if row[*variable_index].is_zero() {
                continue;
            }
            replacement = replacement.sub(
                &Poly::var(&parameter_vars, parameter_index).scale(row[*variable_index].clone()),
            );
        }
        replacements[*pivot] = replacement;
    }

    Some((parameter_vars, replacements))
}

fn linear_subspace_sample_checks(
    equations: &[Poly],
    replacements: &[Poly],
    parameter_count: usize,
) -> bool {
    projective_parameter_samples(parameter_count)
        .into_iter()
        .all(|sample| {
            let Some(point) = replacements
                .iter()
                .map(|replacement| replacement.evaluate(&sample).ok())
                .collect::<Option<Vec<_>>>()
            else {
                return false;
            };
            equations
                .iter()
                .all(|equation| equation.evaluate(&point).is_ok_and(|value| value.is_zero()))
        })
}

fn projective_parameter_samples(parameter_count: usize) -> Vec<Vec<Rational>> {
    match parameter_count {
        0 => Vec::new(),
        1 => vec![vec![Rational::one()]],
        2 => vec![
            vec![Rational::one(), Rational::zero()],
            vec![Rational::zero(), Rational::one()],
            vec![Rational::one(), Rational::one()],
            vec![Rational::one(), -Rational::one()],
            vec![Rational::from_i128(2), Rational::one()],
        ],
        count => {
            let mut samples = Vec::new();
            for index in 0..count {
                let mut sample = vec![Rational::zero(); count];
                sample[index] = Rational::one();
                samples.push(sample);
            }
            samples.push(vec![Rational::one(); count]);
            samples
        }
    }
}

fn homogeneous_row_to_poly(row: &[Rational], variables: &[String]) -> Poly {
    let terms = row.iter().enumerate().filter_map(|(index, coeff)| {
        if coeff.is_zero() {
            None
        } else {
            let mut exp = vec![0; variables.len()];
            exp[index] = 1;
            Some((exp, coeff.clone()))
        }
    });
    Poly::from_terms(variables, terms)
}

fn dehomogenize_projective_linear_ideal(generators: &[Poly], chart: &Chart) -> Option<Vec<Poly>> {
    let replacements = generators
        .first()?
        .vars
        .iter()
        .map(|variable| {
            if !chart.variables.contains(variable) {
                Some(Poly::constant(&chart.variables, Rational::one()))
            } else {
                let index = chart.variables.iter().position(|name| name == variable)?;
                Some(Poly::var(&chart.variables, index))
            }
        })
        .collect::<Option<Vec<_>>>()?;
    let affine_generators = generators
        .iter()
        .map(|generator| generator.substitute(&replacements))
        .collect::<Vec<_>>();

    affine_linear_basis(&affine_generators, &chart.variables)
}

fn projective_point_coordinates(generators: &[Poly]) -> Option<Vec<Rational>> {
    let variable_count = generators.first()?.vars.len();
    let rows = generators
        .iter()
        .map(affine_linear_row)
        .collect::<Option<Vec<_>>>()?;
    let reduced_rows = rref_affine_rows(rows, variable_count)?;
    if reduced_rows.len() + 1 != variable_count {
        return None;
    }

    let pivots = reduced_rows
        .iter()
        .map(|(pivot, _)| *pivot)
        .collect::<Vec<_>>();
    let free_index = (0..variable_count).find(|index| !pivots.contains(index))?;
    let mut values = vec![Rational::zero(); variable_count];
    values[free_index] = Rational::one();
    for (pivot, row) in reduced_rows {
        let mut value = Rational::zero();
        for (index, coeff) in row.iter().take(variable_count).enumerate() {
            if index == pivot || coeff.is_zero() {
                continue;
            }
            value = value - coeff.clone() * values[index].clone();
        }
        values[pivot] = value;
    }
    Some(values)
}

fn projective_generators_vanish_at(generators: &[Poly], point: &[Rational]) -> bool {
    generators
        .iter()
        .all(|generator| generator.evaluate(point).is_ok_and(|value| value.is_zero()))
}

fn generic_multiplicity_on_component(
    chart: &Chart,
    ideal_generators: &[Poly],
    linear_change: Option<&LinearCenterChange>,
) -> usize {
    if let Some(change) = linear_change {
        let center_indices = change
            .indexed_assignments
            .keys()
            .copied()
            .collect::<Vec<_>>();
        if !center_indices.is_empty() {
            let transformed = chart.polynomial.substitute(&change.replacements);
            return transformed.center_order(&center_indices).unwrap_or(0);
        }
    }

    generic_multiplicity(&chart.polynomial, ideal_generators)
}

fn component_intersection_bases(
    left: &SingularComponentChart,
    right: &SingularComponentChart,
    chart: &Chart,
) -> Vec<Vec<Poly>> {
    let mut generators = left.ideal_generators.clone();
    generators.extend(right.ideal_generators.clone());

    if left.affine_linear_center && right.affine_linear_center {
        if let Some(linear_basis) = affine_linear_basis(&generators, &chart.variables) {
            return vec![linear_basis];
        }
    }

    decompose_ideal(&generators)
}

fn affine_linear_basis(generators: &[Poly], variables: &[String]) -> Option<Vec<Poly>> {
    if generators.is_empty() {
        return None;
    }

    let rows = generators
        .iter()
        .map(affine_linear_row)
        .collect::<Option<Vec<_>>>()?;
    let reduced_rows = rref_affine_rows(rows, variables.len())?;
    if reduced_rows.is_empty() {
        return Some(Vec::new());
    }

    Some(
        reduced_rows
            .into_iter()
            .map(|(_, row)| affine_row_to_poly(&row, variables).monic())
            .collect(),
    )
}

fn affine_row_to_poly(row: &[Rational], variables: &[String]) -> Poly {
    let mut terms = Vec::new();
    for (index, coeff) in row.iter().take(variables.len()).enumerate() {
        if coeff.is_zero() {
            continue;
        }
        let mut exp = vec![0; variables.len()];
        exp[index] = 1;
        terms.push((exp, coeff.clone()));
    }
    if !row[variables.len()].is_zero() {
        terms.push((vec![0; variables.len()], row[variables.len()].clone()));
    }
    Poly::from_terms(variables, terms)
}

fn generic_line_singularity_type(
    component: &SingularComponent,
    charts: &[Chart],
) -> Option<(usize, BTreeMap<String, Rational>, SingularityType)> {
    let mut best = None::<(usize, BTreeMap<String, Rational>, SingularityType)>;

    for component_chart in &component.charts {
        if component_chart.dimension != 1 || !component_chart.affine_linear_center {
            continue;
        }
        let Some(chart) = charts
            .iter()
            .find(|chart| chart.id == component_chart.chart_id)
        else {
            continue;
        };
        let Some(candidate) =
            generic_line_singularity_type_on_chart(chart, &component_chart.ideal_generators)
        else {
            continue;
        };

        let replace = best.as_ref().is_none_or(|(_, _, current_type)| {
            generic_type_score(&candidate.1) < generic_type_score(current_type)
        });
        if replace {
            best = Some((chart.id, candidate.0, candidate.1));
        }
    }

    best
}

fn generic_line_singularity_type_on_chart(
    chart: &Chart,
    generators: &[Poly],
) -> Option<(BTreeMap<String, Rational>, SingularityType)> {
    let rows = generators
        .iter()
        .map(affine_linear_row)
        .collect::<Option<Vec<_>>>()?;
    let reduced_rows = rref_affine_rows(rows, chart.variables.len())?;
    if reduced_rows.len() + 1 != chart.variables.len() {
        return None;
    }

    let pivot_indices = reduced_rows
        .iter()
        .map(|(pivot, _)| *pivot)
        .collect::<Vec<_>>();
    let free_index = (0..chart.variables.len()).find(|index| !pivot_indices.contains(index))?;
    let change = linear_center_change(generators, &chart.variables)?;
    let transformed = chart.polynomial.substitute(&change.replacements);
    let mut best = None::<(BTreeMap<String, Rational>, SingularityType)>;

    for value in generic_parameter_samples() {
        let point = point_on_line(&reduced_rows, &chart.variables, free_index, value.clone());
        let mut specialization = BTreeMap::new();
        specialization.insert(free_index, value);
        let transverse = transformed.specialize(&specialization);
        let transverse = project_to_variables(&transverse, &pivot_indices)?;
        let singularity_type = analyze_transverse_singularity(&transverse);
        if !singularity_type.is_singular {
            continue;
        }

        let replace = best.as_ref().is_none_or(|(_, current_type)| {
            generic_type_score(&singularity_type) < generic_type_score(current_type)
        });
        if replace {
            best = Some((point, singularity_type));
        }
    }

    best
}

fn generic_parameter_samples() -> Vec<Rational> {
    [0, 1, -1, 2, -2, 3, -3, 4, -4]
        .into_iter()
        .map(Rational::from_i128)
        .collect()
}

fn point_on_line(
    reduced_rows: &[(usize, Vec<Rational>)],
    variables: &[String],
    free_index: usize,
    free_value: Rational,
) -> BTreeMap<String, Rational> {
    let mut values = vec![Rational::zero(); variables.len()];
    values[free_index] = free_value;

    for (pivot, row) in reduced_rows {
        let mut value = -row[variables.len()].clone();
        for (index, coeff) in row.iter().take(variables.len()).enumerate() {
            if index == *pivot || coeff.is_zero() {
                continue;
            }
            value = value - coeff.clone() * values[index].clone();
        }
        values[*pivot] = value;
    }

    variables.iter().cloned().zip(values).collect()
}

fn project_to_variables(poly: &Poly, keep_indices: &[usize]) -> Option<Poly> {
    let variables = keep_indices
        .iter()
        .map(|index| poly.vars[*index].clone())
        .collect::<Vec<_>>();
    let mut terms = Vec::new();

    for (exp, coeff) in &poly.terms {
        for index in 0..poly.vars.len() {
            if !keep_indices.contains(&index) && exp[index] != 0 {
                return None;
            }
        }
        let next_exp = keep_indices.iter().map(|index| exp[*index]).collect();
        terms.push((next_exp, coeff.clone()));
    }

    Some(Poly::from_terms(&variables, terms))
}

fn analyze_transverse_singularity(poly: &Poly) -> SingularityType {
    analyze_transverse_singularity_with_bound(poly, None)
}

fn analyze_transverse_singularity_with_bound(
    poly: &Poly,
    formal_degree_bound: Option<usize>,
) -> SingularityType {
    let center_indices = (0..poly.vars.len()).collect::<Vec<_>>();
    let multiplicity = poly.center_order(&center_indices).unwrap_or(0);
    let tangent_cone_degree = (!poly.is_zero()).then_some(multiplicity);
    let tangent_cone = tangent_cone_degree
        .map(|degree| poly.homogeneous_part(degree))
        .unwrap_or_else(|| Poly::zero(&poly.vars));
    let quadratic_rank = if tangent_cone_degree == Some(2) {
        Some(quadratic_form_rank(&tangent_cone))
    } else {
        None
    };
    let mut label =
        classify_singularity(multiplicity, poly.vars.len(), quadratic_rank, &tangent_cone);

    if multiplicity == 2 {
        let local_a_number = match (poly.vars.len(), quadratic_rank) {
            (2, Some(2)) => Some(1),
            (2, Some(1)) => local_bivariate_a_number(poly),
            (3, Some(3)) => Some(1),
            (3, Some(2)) => Some(local_trivariate_a_number_with_bound(
                poly,
                formal_degree_bound.unwrap_or_else(|| formal_series_degree_bound(poly)),
            ))
            .flatten(),
            _ => None,
        };
        if let Some(a_number) = local_a_number.filter(|number| *number > 0) {
            label = format!("cA{a_number}");
        }
    }

    SingularityType {
        label,
        multiplicity,
        tangent_cone,
        tangent_cone_degree,
        quadratic_rank,
        embedding_dimension: poly.vars.len(),
        is_singular: multiplicity >= 2,
    }
}

fn local_bivariate_a_number(poly: &Poly) -> Option<usize> {
    if poly.vars.len() != 2 {
        return None;
    }

    // For a rank-one quadratic germ, choose a variable in the nondegenerate
    // Hessian direction. The formal implicit function theorem gives a unique
    // critical series x = phi(y). The order of f(phi(y), y) is n + 1 for A_n.
    let split_index = (0..2).find(|index| {
        poly.terms
            .get(&quadratic_exponent(*index, 2))
            .is_some_and(|coefficient| !coefficient.is_zero())
    })?;
    let parameter_index = 1 - split_index;
    let split_quadratic_coefficient = poly.terms.get(&quadratic_exponent(split_index, 2))?.clone();
    let implicit_linear_coefficient = Rational::from_i128(2) * split_quadratic_coefficient;

    let degree_bound = poly.total_degree().saturating_sub(1);
    let max_degree = degree_bound
        .saturating_mul(degree_bound)
        .saturating_add(1)
        .clamp(3, 512);
    let gradient = poly.derivative(split_index);
    let mut critical_series = vec![Rational::zero(); max_degree + 1];

    for degree in 1..max_degree {
        let gradient_series = substitute_bivariate_series(
            &gradient,
            split_index,
            parameter_index,
            &critical_series,
            degree,
        );
        critical_series[degree] =
            -gradient_series[degree].clone() / implicit_linear_coefficient.clone();

        let residual_degree = degree + 1;
        let residual = substitute_bivariate_series(
            poly,
            split_index,
            parameter_index,
            &critical_series,
            residual_degree,
        );
        if !residual[residual_degree].is_zero() {
            return Some(residual_degree - 1);
        }
    }
    None
}

fn local_trivariate_a_number(poly: &Poly) -> Option<usize> {
    local_trivariate_a_number_with_bound(poly, formal_series_degree_bound(poly))
}

fn formal_series_degree_bound(poly: &Poly) -> usize {
    let degree_bound = poly.total_degree().saturating_sub(1);
    degree_bound
        .saturating_mul(degree_bound)
        .saturating_add(1)
        .clamp(4, 512)
}

fn local_trivariate_a_number_with_bound(poly: &Poly, max_degree: usize) -> Option<usize> {
    if poly.vars.len() != 3 {
        return None;
    }

    let split_indices = (0..3)
        .flat_map(|left| ((left + 1)..3).map(move |right| [left, right]))
        .find(|[left, right]| {
            let a = quadratic_hessian_entry(poly, *left, *left);
            let b = quadratic_hessian_entry(poly, *left, *right);
            let d = quadratic_hessian_entry(poly, *right, *right);
            !(a * d - b.clone() * b).is_zero()
        })?;
    let parameter_index = (0..3).find(|index| !split_indices.contains(index))?;
    let [left, right] = split_indices;
    let a = quadratic_hessian_entry(poly, left, left);
    let b = quadratic_hessian_entry(poly, left, right);
    let d = quadratic_hessian_entry(poly, right, right);
    let determinant = a.clone() * d.clone() - b.clone() * b.clone();

    let gradients = [poly.derivative(left), poly.derivative(right)];
    let mut left_series = vec![Rational::zero(); max_degree + 1];
    let mut right_series = vec![Rational::zero(); max_degree + 1];

    for degree in 1..max_degree {
        let left_residual = substitute_trivariate_series(
            &gradients[0],
            [left, right],
            parameter_index,
            [&left_series, &right_series],
            degree,
        )[degree]
            .clone();
        let right_residual = substitute_trivariate_series(
            &gradients[1],
            [left, right],
            parameter_index,
            [&left_series, &right_series],
            degree,
        )[degree]
            .clone();
        left_series[degree] = (-d.clone() * left_residual.clone()
            + b.clone() * right_residual.clone())
            / determinant.clone();
        right_series[degree] =
            (b.clone() * left_residual - a.clone() * right_residual) / determinant.clone();

        let residual_degree = degree + 1;
        let residual = substitute_trivariate_series(
            poly,
            [left, right],
            parameter_index,
            [&left_series, &right_series],
            residual_degree,
        );
        if !residual[residual_degree].is_zero() {
            return Some(residual_degree - 1);
        }
    }
    None
}

fn quadratic_hessian_entry(poly: &Poly, left: usize, right: usize) -> Rational {
    let mut exponent = vec![0usize; poly.vars.len()];
    exponent[left] += 1;
    exponent[right] += 1;
    let coefficient = poly
        .terms
        .get(&exponent)
        .cloned()
        .unwrap_or_else(Rational::zero);
    if left == right {
        Rational::from_i128(2) * coefficient
    } else {
        coefficient
    }
}

fn substitute_trivariate_series(
    poly: &Poly,
    series_indices: [usize; 2],
    parameter_index: usize,
    series: [&[Rational]; 2],
    max_degree: usize,
) -> Vec<Rational> {
    let mut result = vec![Rational::zero(); max_degree + 1];
    for (exponent, coefficient) in &poly.terms {
        let left_power = truncated_series_power(series[0], exponent[series_indices[0]], max_degree);
        let right_power =
            truncated_series_power(series[1], exponent[series_indices[1]], max_degree);
        let product = multiply_truncated_series(&left_power, &right_power, max_degree);
        let parameter_power = exponent[parameter_index];
        for (degree, value) in product.into_iter().enumerate() {
            let output_degree = degree + parameter_power;
            if output_degree > max_degree || value.is_zero() {
                continue;
            }
            result[output_degree] = result[output_degree].clone() + coefficient.clone() * value;
        }
    }
    result
}

fn quadratic_exponent(index: usize, power: usize) -> Vec<usize> {
    let mut exponent = vec![0usize; 2];
    exponent[index] = power;
    exponent
}

fn substitute_bivariate_series(
    poly: &Poly,
    series_index: usize,
    parameter_index: usize,
    series: &[Rational],
    max_degree: usize,
) -> Vec<Rational> {
    let mut result = vec![Rational::zero(); max_degree + 1];
    for (exponent, coefficient) in &poly.terms {
        let series_power = truncated_series_power(series, exponent[series_index], max_degree);
        let parameter_power = exponent[parameter_index];
        for (degree, value) in series_power.into_iter().enumerate() {
            let output_degree = degree + parameter_power;
            if output_degree > max_degree || value.is_zero() {
                continue;
            }
            result[output_degree] = result[output_degree].clone() + coefficient.clone() * value;
        }
    }
    result
}

fn truncated_series_power(
    series: &[Rational],
    exponent: usize,
    max_degree: usize,
) -> Vec<Rational> {
    let mut result = vec![Rational::zero(); max_degree + 1];
    result[0] = Rational::one();
    for _ in 0..exponent {
        result = multiply_truncated_series(&result, series, max_degree);
    }
    result
}

fn multiply_truncated_series(
    left: &[Rational],
    right: &[Rational],
    max_degree: usize,
) -> Vec<Rational> {
    let mut result = vec![Rational::zero(); max_degree + 1];
    for (left_degree, left_coefficient) in left.iter().enumerate().take(max_degree + 1) {
        if left_coefficient.is_zero() {
            continue;
        }
        for (right_degree, right_coefficient) in
            right.iter().enumerate().take(max_degree + 1 - left_degree)
        {
            if right_coefficient.is_zero() {
                continue;
            }
            let degree = left_degree + right_degree;
            result[degree] =
                result[degree].clone() + left_coefficient.clone() * right_coefficient.clone();
        }
    }
    result
}

fn generic_type_score(singularity_type: &SingularityType) -> (usize, usize, usize) {
    let ca_number = singularity_type
        .label
        .strip_prefix("cA")
        .and_then(|value| value.parse::<usize>().ok());
    (
        usize::from(ca_number.is_none()),
        singularity_type.multiplicity,
        ca_number.unwrap_or(usize::MAX),
    )
}

fn component_display_key(component: &SingularComponent) -> String {
    component
        .charts
        .iter()
        .map(|chart| {
            format!(
                "{}:{}",
                chart.chart_id,
                chart
                    .ideal_generators
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join("|")
}

fn assignment_key(assignments: &BTreeMap<String, Rational>) -> String {
    assignments
        .iter()
        .map(|(name, value)| format!("{name}={value}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn assignment_key_option(assignments: &Option<BTreeMap<String, Rational>>) -> String {
    assignments.as_ref().map(assignment_key).unwrap_or_default()
}

fn require_full_point_assignment(
    chart: &Chart,
    assignments: &BTreeMap<String, Rational>,
) -> Result<(), String> {
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
    Ok(())
}

fn analyze_chart_point_singularity(
    chart: &Chart,
    assignments: &BTreeMap<String, Rational>,
) -> SingularityType {
    let indexed_assignments = assignments
        .iter()
        .filter_map(|(name, value)| {
            chart
                .polynomial
                .variable_index(name)
                .map(|index| (index, value.clone()))
        })
        .collect::<BTreeMap<_, _>>();
    let translated = chart.polynomial.translated_by(&indexed_assignments);
    let center_indices = (0..chart.variables.len()).collect::<Vec<_>>();
    let multiplicity = translated.center_order(&center_indices).unwrap_or(0);
    let tangent_cone_degree = (!translated.is_zero()).then_some(multiplicity);
    let tangent_cone = tangent_cone_degree
        .map(|degree| translated.homogeneous_part(degree))
        .unwrap_or_else(|| Poly::zero(&chart.variables));
    let quadratic_rank = if tangent_cone_degree == Some(2) {
        Some(quadratic_form_rank(&tangent_cone))
    } else {
        None
    };
    let mut label = classify_singularity(
        multiplicity,
        chart.variables.len(),
        quadratic_rank,
        &tangent_cone,
    );
    if multiplicity == 2 && quadratic_rank == Some(2) && chart.variables.len() == 3 {
        if let Some(a_number) = local_trivariate_a_number(&translated) {
            label = format!("cA{a_number}");
        }
    }

    SingularityType {
        label,
        multiplicity,
        tangent_cone,
        tangent_cone_degree,
        quadratic_rank,
        embedding_dimension: chart.variables.len(),
        is_singular: multiplicity >= 2,
    }
}

fn classify_singularity(
    multiplicity: usize,
    variable_count: usize,
    quadratic_rank: Option<usize>,
    tangent_cone: &Poly,
) -> String {
    if multiplicity == 0 {
        return "point is not on the hypersurface".to_string();
    }
    if multiplicity == 1 {
        return "smooth point".to_string();
    }

    match (multiplicity, variable_count, quadratic_rank) {
        (2, 2, Some(2)) => "ordinary double point (node/A1)".to_string(),
        (2, 2, Some(1)) => "plane curve double point with repeated tangent".to_string(),
        (2, 3, Some(3)) => "ordinary double point (A1)".to_string(),
        (2, 3, Some(2)) => "surface double point with rank-2 quadratic tangent cone".to_string(),
        (2, _, Some(rank)) if rank > 0 => {
            format!("double point with rank-{rank} quadratic tangent cone")
        }
        _ if !tangent_cone.is_zero() => format!("multiplicity-{multiplicity} singularity"),
        _ => "undetermined singularity type".to_string(),
    }
}

fn quadratic_form_rank(poly: &Poly) -> usize {
    let size = poly.vars.len();
    let mut matrix = vec![vec![Rational::zero(); size]; size];
    for (exp, coeff) in &poly.terms {
        if exp.iter().sum::<usize>() != 2 {
            continue;
        }

        let support = exp
            .iter()
            .enumerate()
            .filter_map(|(index, power)| (*power > 0).then_some((index, *power)))
            .collect::<Vec<_>>();
        match support.as_slice() {
            [(index, 2)] => {
                matrix[*index][*index] =
                    matrix[*index][*index].clone() + coeff.clone() * Rational::from_i128(2);
            }
            [(left, 1), (right, 1)] => {
                matrix[*left][*right] = matrix[*left][*right].clone() + coeff.clone();
                matrix[*right][*left] = matrix[*right][*left].clone() + coeff.clone();
            }
            _ => {}
        }
    }
    rational_matrix_rank(matrix)
}

fn rational_matrix_rank(mut matrix: Vec<Vec<Rational>>) -> usize {
    let row_count = matrix.len();
    let column_count = matrix.first().map_or(0, Vec::len);
    let mut rank = 0usize;

    for column in 0..column_count {
        let Some(pivot_row) = (rank..row_count).find(|row| !matrix[*row][column].is_zero()) else {
            continue;
        };
        matrix.swap(rank, pivot_row);
        let pivot = matrix[rank][column].clone();
        for column_index in column..column_count {
            matrix[rank][column_index] = matrix[rank][column_index].clone() / pivot.clone();
        }

        let normalized_pivot = matrix[rank].clone();
        for row in 0..row_count {
            if row == rank || matrix[row][column].is_zero() {
                continue;
            }
            let factor = matrix[row][column].clone();
            for column_index in column..column_count {
                matrix[row][column_index] = matrix[row][column_index].clone()
                    - factor.clone() * normalized_pivot[column_index].clone();
            }
        }
        rank += 1;
        if rank == row_count {
            break;
        }
    }

    rank
}

fn coordinate_assignments_from_basis(
    basis: &[Poly],
    variables: &[String],
) -> Option<BTreeMap<String, Rational>> {
    let mut assignments = BTreeMap::new();
    for generator in basis {
        let (variable_index, value) = coordinate_assignment(generator)?;
        let name = variables[variable_index].clone();
        if let Some(existing) = assignments.insert(name.clone(), value.clone()) {
            if existing != value {
                return None;
            }
        }
    }
    Some(assignments)
}

fn coordinate_assignment(poly: &Poly) -> Option<(usize, Rational)> {
    let mut variable_term = None::<(usize, Rational)>;
    let mut constant = Rational::zero();

    for (exp, coeff) in &poly.terms {
        let degree = exp.iter().sum::<usize>();
        if degree == 0 {
            constant = constant + coeff.clone();
            continue;
        }

        if degree != 1 {
            return None;
        }
        let variable_index = exp.iter().position(|power| *power == 1)?;
        if variable_term
            .replace((variable_index, coeff.clone()))
            .is_some()
        {
            return None;
        }
    }

    let (variable_index, coeff) = variable_term?;
    if coeff.is_zero() {
        return None;
    }
    Some((variable_index, -constant / coeff))
}

fn linear_center_change(generators: &[Poly], variables: &[String]) -> Option<LinearCenterChange> {
    if generators.is_empty() {
        return None;
    }

    let mut rows = Vec::new();
    for generator in generators {
        rows.push(affine_linear_row(generator)?);
    }

    let reduced_rows = rref_affine_rows(rows, variables.len())?;
    if reduced_rows.len() < 2 {
        return None;
    }

    let mut replacements = (0..variables.len())
        .map(|index| Poly::var(variables, index))
        .collect::<Vec<_>>();
    let mut assignments = BTreeMap::new();
    let mut indexed_assignments = BTreeMap::new();
    let mut display = Vec::new();

    for (pivot, row) in reduced_rows {
        let constant = row[variables.len()].clone();
        let mut replacement = Poly::var(variables, pivot).sub(&Poly::constant(variables, constant));
        for (index, coeff) in row.iter().take(variables.len()).enumerate() {
            if index == pivot || coeff.is_zero() {
                continue;
            }
            replacement = replacement.sub(&Poly::var(variables, index).scale(coeff.clone()));
        }

        display.push(format!(
            "affine change: {} -> {}",
            variables[pivot], replacement
        ));
        replacements[pivot] = replacement;
        assignments.insert(variables[pivot].clone(), Rational::zero());
        indexed_assignments.insert(pivot, Rational::zero());
    }

    Some(LinearCenterChange {
        replacements,
        assignments,
        indexed_assignments,
        display,
    })
}

fn affine_linear_row(poly: &Poly) -> Option<Vec<Rational>> {
    let mut row = vec![Rational::zero(); poly.vars.len() + 1];
    for (exp, coeff) in &poly.terms {
        let degree = exp.iter().sum::<usize>();
        if degree == 0 {
            row[poly.vars.len()] = row[poly.vars.len()].clone() + coeff.clone();
        } else if degree == 1 {
            let variable_index = exp.iter().position(|power| *power == 1)?;
            row[variable_index] = row[variable_index].clone() + coeff.clone();
        } else {
            return None;
        }
    }
    Some(row)
}

fn rref_affine_rows(
    mut rows: Vec<Vec<Rational>>,
    variable_count: usize,
) -> Option<Vec<(usize, Vec<Rational>)>> {
    let mut pivot_row = 0usize;
    for column in 0..variable_count {
        let Some(source_row) = (pivot_row..rows.len()).find(|row| !rows[*row][column].is_zero())
        else {
            continue;
        };

        rows.swap(pivot_row, source_row);
        let pivot = rows[pivot_row][column].clone();
        for value in &mut rows[pivot_row] {
            *value = value.clone() / pivot.clone();
        }

        let normalized_pivot = rows[pivot_row].clone();
        for row in 0..rows.len() {
            if row == pivot_row || rows[row][column].is_zero() {
                continue;
            }
            let factor = rows[row][column].clone();
            for col in 0..=variable_count {
                rows[row][col] =
                    rows[row][col].clone() - factor.clone() * normalized_pivot[col].clone();
            }
        }

        pivot_row += 1;
        if pivot_row == rows.len() {
            break;
        }
    }

    let mut reduced_rows = Vec::new();
    for row in rows {
        let pivot = row
            .iter()
            .take(variable_count)
            .position(|value| !value.is_zero());
        match pivot {
            Some(pivot) => reduced_rows.push((pivot, row)),
            None if !row[variable_count].is_zero() => return None,
            None => {}
        }
    }
    reduced_rows.sort_by_key(|(pivot, _)| *pivot);
    Some(reduced_rows)
}

fn projective_homogenized_key(
    ideal_generators: &[Poly],
    chart: &Chart,
    projective_vars: &[String],
) -> Option<String> {
    let unit_index = chart.id;
    let mut affine_to_projective = Vec::new();
    for variable in &chart.variables {
        affine_to_projective.push(projective_vars.iter().position(|name| name == variable)?);
    }

    let homogenized = ideal_generators
        .iter()
        .map(|poly| homogenize_from_chart(poly, unit_index, &affine_to_projective, projective_vars))
        .collect::<Vec<_>>();
    let basis = reduced_groebner_basis(&homogenized);
    Some(groebner_key(&basis))
}

fn homogenize_from_chart(
    poly: &Poly,
    unit_index: usize,
    affine_to_projective: &[usize],
    projective_vars: &[String],
) -> Poly {
    let degree = poly.total_degree();
    let terms = poly.terms.iter().map(|(affine_exp, coeff)| {
        let mut projective_exp = vec![0usize; projective_vars.len()];
        let affine_degree = affine_exp.iter().sum::<usize>();
        projective_exp[unit_index] = degree - affine_degree;
        for (affine_index, power) in affine_exp.iter().enumerate() {
            projective_exp[affine_to_projective[affine_index]] = *power;
        }
        (projective_exp, coeff.clone())
    });
    Poly::from_terms(projective_vars, terms)
}

fn blowup_substitutions(
    chart: &Chart,
    assignments: &BTreeMap<usize, Rational>,
    exceptional_index: usize,
) -> BlowupSubstitutions {
    let exceptional = Poly::var(&chart.variables, exceptional_index);
    let mut polys = Vec::new();
    let mut display = Vec::new();

    for index in 0..chart.variables.len() {
        let variable = Poly::var(&chart.variables, index);
        if let Some(value) = assignments.get(&index) {
            let replacement = if index == exceptional_index {
                Poly::constant(&chart.variables, value.clone()).add(&variable)
            } else {
                Poly::constant(&chart.variables, value.clone()).add(&exceptional.mul(&variable))
            };
            display.push(format!("{} -> {}", chart.variables[index], replacement));
            polys.push(replacement);
        } else {
            display.push(format!(
                "{} -> {}",
                chart.variables[index], chart.variables[index]
            ));
            polys.push(variable);
        }
    }

    BlowupSubstitutions { polys, display }
}

fn substitute_zero_centered_coordinate_blowup(
    poly: &Poly,
    center_indices: &[usize],
    exceptional_index: usize,
) -> Poly {
    let terms = poly.terms.iter().map(|(exponent, coefficient)| {
        let mut transformed_exponent = exponent.clone();
        for center_index in center_indices {
            if *center_index != exceptional_index {
                transformed_exponent[exceptional_index] += exponent[*center_index];
            }
        }
        (transformed_exponent, coefficient.clone())
    });
    Poly::from_terms(&poly.vars, terms)
}

fn validate_center_in_singular_locus(
    chart: &Chart,
    assignments: &BTreeMap<usize, Rational>,
) -> Result<(), String> {
    let hypersurface_restriction = chart.polynomial.specialize(assignments);
    if !hypersurface_restriction.is_zero() {
        return Err(format!(
            "center is not contained in chart {id}'s hypersurface; use --force to override",
            id = chart.id
        ));
    }

    let partials = chart.polynomial.partials();
    for partial in partials {
        if !partial.specialize(assignments).is_zero() {
            return Err(format!(
                "center is not contained in chart {id}'s singular locus; use --force to override",
                id = chart.id
            ));
        }
    }
    Ok(())
}

fn enumerate_singular_points(
    chart: &Chart,
    equations: &[Poly],
    bound: i32,
    limit: usize,
) -> (Vec<SingularityPoint>, bool) {
    let candidates = rational_grid_values(bound);
    let mut current = Vec::new();
    let mut points = Vec::new();
    let mut limit_hit = false;
    enumerate_points_recursive(
        chart,
        equations,
        &candidates,
        limit,
        &mut current,
        &mut points,
        &mut limit_hit,
    );
    (points, limit_hit)
}

fn enumerate_points_recursive(
    chart: &Chart,
    equations: &[Poly],
    candidates: &[Rational],
    limit: usize,
    current: &mut Vec<Rational>,
    points: &mut Vec<SingularityPoint>,
    limit_hit: &mut bool,
) {
    if *limit_hit {
        return;
    }

    if current.len() == chart.variables.len() {
        let is_singular = equations.iter().all(|equation| {
            equation
                .evaluate(current)
                .is_ok_and(|value| value.is_zero())
        });
        if is_singular {
            points.push(SingularityPoint {
                chart_id: chart.id,
                values: current.clone(),
            });
            if points.len() >= limit {
                *limit_hit = true;
            }
        }
        return;
    }

    for candidate in candidates {
        current.push(candidate.clone());
        enumerate_points_recursive(
            chart, equations, candidates, limit, current, points, limit_hit,
        );
        current.pop();
        if *limit_hit {
            break;
        }
    }
}

fn rational_grid_values(bound: i32) -> Vec<Rational> {
    let bound = bound.max(1);
    let mut values = BTreeSet::new();
    for denominator in 1..=bound {
        for numerator in -bound..=bound {
            values.insert(Rational::new(numerator as i128, denominator as i128));
        }
    }
    values.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::{
        AppState, AutoResolutionOptions, AutoResolutionProgress, AutoResolutionStatus,
        parse_singular_minimal_primes, singular_minimal_primes,
    };
    use crate::rational::Rational;
    use std::collections::BTreeMap;

    const STALLING_SURFACE: &str = "x0^4*x1^2*x2^2-2*x0^2*x1^4*x2^2+x1^6*x2^2+4*x0^2*x1^3*x2^3-4*x1^5*x2^3-2*x0^2*x1^2*x2^4+6*x1^4*x2^4-4*x1^3*x2^5+x1^2*x2^6-4*x0^4*x1^2*x2*x3+6*x0^3*x1^3*x2*x3-2*x0^2*x1^4*x2*x3+2*x0*x1^5*x2*x3-2*x1^6*x2*x3-6*x0^3*x1^2*x2^2*x3-2*x0*x1^4*x2^2*x3+8*x1^5*x2^2*x3+2*x0^2*x1^2*x2^3*x3-2*x0*x1^3*x2^3*x3-12*x1^4*x2^3*x3+2*x0*x1^2*x2^4*x3+8*x1^3*x2^4*x3-2*x1^2*x2^5*x3+4*x0^4*x1^2*x3^2-8*x0^3*x1^3*x3^2+5*x0^2*x1^4*x3^2-2*x0*x1^5*x3^2+x1^6*x3^2+2*x0^4*x1*x2*x3^2+4*x0^3*x1^2*x2*x3^2-2*x0*x1^4*x2*x3^2-4*x1^5*x2*x3^2+4*x0^3*x1*x2^2*x3^2-7*x0^2*x1^2*x2^2*x3^2+10*x0*x1^3*x2^2*x3^2+6*x1^4*x2^2*x3^2+2*x0^2*x1*x2^3*x3^2-6*x0*x1^2*x2^3*x3^2-4*x1^3*x2^3*x3^2+x1^2*x2^4*x3^2-4*x0^4*x1*x3^3+6*x0^3*x1^2*x3^3-6*x0^2*x1^3*x3^3+4*x0*x1^4*x3^3-6*x0^3*x1*x2*x3^3+8*x0^2*x1^2*x2*x3^3-8*x0*x1^3*x2*x3^3-2*x0^2*x1*x2^2*x3^3+4*x0*x1^2*x2^2*x3^3+x0^4*x3^4";
    const REPORTED_INTERSECTION_SURFACE: &str = "x0^6*x2^2-4*x0^5*x1*x2^2+6*x0^4*x1^2*x2^2-4*x0^3*x1^3*x2^2+x0^2*x1^4*x2^2-4*x0^5*x2^3+10*x0^4*x1*x2^3-8*x0^3*x1^2*x2^3+2*x0^2*x1^3*x2^3+6*x0^4*x2^4-8*x0^3*x1*x2^4+3*x0^2*x1^2*x2^4-4*x0^3*x2^5+2*x0^2*x1*x2^5+x0^2*x2^6-2*x0^6*x2*x3+6*x0^5*x1*x2*x3-6*x0^4*x1^2*x2*x3+2*x0^3*x1^3*x2*x3+6*x0^5*x2^2*x3-12*x0^4*x1*x2^2*x3+2*x0^3*x1^2*x2^2*x3+6*x0^2*x1^3*x2^2*x3-2*x0*x1^4*x2^2*x3-6*x0^4*x2^3*x3+2*x0^3*x1*x2^3*x3+6*x0^2*x1^2*x2^3*x3-4*x0*x1^3*x2^3*x3+2*x0^3*x2^4*x3+6*x0^2*x1*x2^4*x3-4*x0*x1^2*x2^4*x3-2*x0*x1*x2^5*x3+x0^6*x3^2-2*x0^5*x1*x3^2+x0^4*x1^2*x3^2-2*x0^4*x1*x2*x3^2+6*x0^3*x1^2*x2*x3^2-4*x0^2*x1^3*x2*x3^2-5*x0^4*x2^2*x3^2+16*x0^3*x1*x2^2*x3^2-11*x0^2*x1^2*x2^2*x3^2+x1^4*x2^2*x3^2+6*x0^3*x2^3*x3^2-14*x0^2*x1*x2^3*x3^2+4*x0*x1^2*x2^3*x3^2+2*x1^3*x2^3*x3^2-2*x0^2*x2^4*x3^2+2*x0*x1*x2^4*x3^2+x1^2*x2^4*x3^2-2*x0^5*x3^3+4*x0^4*x1*x3^3-2*x0^3*x1^2*x3^3+4*x0^4*x2*x3^3-8*x0^3*x1*x2*x3^3+2*x0^2*x1^2*x2*x3^3+2*x0*x1^3*x2*x3^3-2*x0^3*x2^2*x3^3+2*x0^2*x1*x2^2*x3^3+2*x0*x1^2*x2^2*x3^3-2*x1^3*x2^2*x3^3+2*x0*x1*x2^3*x3^3-2*x1^2*x2^3*x3^3+x0^4*x3^4-2*x0^3*x1*x3^4+x0^2*x1^2*x3^4-2*x0^3*x2*x3^4+4*x0^2*x1*x2*x3^4-2*x0*x1^2*x2*x3^4+x0^2*x2^2*x3^4-2*x0*x1*x2^2*x3^4+x1^2*x2^2*x3^4";
    const SEVEN_LINE_SURFACE: &str = include_str!("../fixtures/seven_line_surface.txt");

    #[test]
    fn creates_projective_charts() {
        let state = AppState::new(1, "x0*x2^2 - x1^3").unwrap();
        assert_eq!(state.charts.len(), 3);
        assert_eq!(state.degree, 3);
    }

    #[test]
    fn parses_singular_minimal_prime_protocol() {
        let vars = vec!["x".to_string(), "y".to_string()];
        let output = "BLOWUP_COUNT=2\nBLOWUP_BEGIN=1\nBLOWUP_GEN=x\nBLOWUP_GEN=y\nBLOWUP_END=1\nBLOWUP_BEGIN=2\nBLOWUP_GEN=x-1/2\nBLOWUP_END=2\n";
        let components = parse_singular_minimal_primes(output, &vars).unwrap();

        assert_eq!(components.len(), 2);
        assert_eq!(components[0][0].to_string(), "x");
        assert_eq!(components[0][1].to_string(), "y");
        assert_eq!(components[1][0].to_string(), "x - 1/2");
    }

    #[test]
    fn singular_backend_decomposes_a_finite_jacobian_locus_when_available() {
        if std::process::Command::new("Singular")
            .arg("--version")
            .output()
            .is_err()
        {
            return;
        }
        let vars = vec!["x1".to_string(), "x2".to_string()];
        let state = AppState::new(1, "x1^2+x2^2").unwrap();
        let chart = state.chart(0).unwrap();
        let equations = super::singularity_equations(chart);
        let points = vec![super::SingularityPoint {
            chart_id: 0,
            values: vec![Rational::zero(), Rational::zero()],
        }];
        let components = singular_minimal_primes(&equations, &vars, Some(&points)).unwrap();

        assert_eq!(components.len(), 1);
        assert_eq!(components[0].len(), 2);
    }

    #[test]
    fn creates_projective_threefold_charts() {
        let state = AppState::new(3, "x0*x1*x4^7 + x2^9").unwrap();
        assert_eq!(state.charts.len(), 5);
        assert_eq!(state.degree, 9);
    }

    #[test]
    fn blows_up_cusp_origin_on_affine_chart() {
        let mut state = AppState::new(1, "x0*x2^2 - x1^3").unwrap();
        let mut center = BTreeMap::new();
        center.insert("x1".to_string(), Rational::zero());
        center.insert("x2".to_string(), Rational::zero());
        state.blowup_coordinate_center(0, center, false).unwrap();

        assert!(!state.chart(0).unwrap().active);
        assert_eq!(state.blowups.len(), 1);
        assert_eq!(state.active_charts().count(), 4);
    }

    #[test]
    fn blows_up_arbitrary_point_on_later_chart_with_force() {
        let mut state = AppState::new(1, "x0*x2^2 - x1^3").unwrap();
        let mut origin = BTreeMap::new();
        origin.insert("x1".to_string(), Rational::zero());
        origin.insert("x2".to_string(), Rational::zero());
        state
            .blowup_coordinate_center(0, origin.clone(), false)
            .unwrap();

        state.blowup_coordinate_center(3, origin, true).unwrap();

        assert_eq!(state.blowups.len(), 2);
        assert_eq!(state.blowups[1].input_chart, 3);
        assert_eq!(state.blowups[1].multiplicity, 1);
        assert!(!state.chart(3).unwrap().active);
    }

    #[test]
    fn analyzes_plane_curve_double_point_types() {
        let node = AppState::new(1, "x1*x2").unwrap();
        let mut origin = BTreeMap::new();
        origin.insert("x1".to_string(), Rational::zero());
        origin.insert("x2".to_string(), Rational::zero());
        let node_type = node.analyze_point(0, origin.clone()).unwrap();

        assert_eq!(node_type.multiplicity, 2);
        assert_eq!(node_type.quadratic_rank, Some(2));
        assert_eq!(node_type.label, "ordinary double point (node/A1)");

        let cusp = AppState::new(1, "x0*x2^2 - x1^3").unwrap();
        let cusp_type = cusp.analyze_point(0, origin).unwrap();

        assert_eq!(cusp_type.multiplicity, 2);
        assert_eq!(cusp_type.quadratic_rank, Some(1));
        assert_eq!(
            cusp_type.label,
            "plane curve double point with repeated tangent"
        );
    }

    #[test]
    fn classifies_local_trivariate_ca_germs() {
        let variables = ["x", "y", "z"].map(str::to_string);
        let origin = [Rational::zero(), Rational::zero(), Rational::zero()];

        let diagonal =
            super::analyze_affine_polynomial_at("x^2+y^2+z^5", &variables, &origin).unwrap();
        assert_eq!(diagonal.label, "cA4");
        assert_eq!(diagonal.quadratic_rank, Some(2));

        let mixed =
            super::analyze_affine_polynomial_at("(x-z)^2+(y-2*z)^2+z^6", &variables, &origin)
                .unwrap();
        assert_eq!(mixed.label, "cA5");
        assert_eq!(mixed.quadratic_rank, Some(2));
    }

    #[test]
    fn auto_crepant_resolution_resolves_branch_quadruple_point() {
        let mut state = AppState::new(2, "x1^4 + x2^4 + x3^4").unwrap();

        let result = state
            .resolve_crepant_with_progress(AutoResolutionOptions { max_steps: 4 }, |_| {})
            .unwrap();

        assert_eq!(result.status, AutoResolutionStatus::Resolved);
        assert_eq!(result.steps.len(), 1);
        assert_eq!(result.steps[0].center_dimension, 0);
        assert_eq!(result.steps[0].ambient_codimension, 3);
        assert_eq!(result.steps[0].multiplicity, 4);
        assert!(state.singular_components().is_empty());
    }

    #[test]
    fn auto_crepant_resolution_skips_branch_surface_double_point() {
        let mut state = AppState::new(2, "x0*x1 - x2^2").unwrap();

        let result = state
            .resolve_crepant_with_progress(AutoResolutionOptions { max_steps: 4 }, |_| {})
            .unwrap();

        assert_eq!(result.status, AutoResolutionStatus::NoCrepantCenter);
        assert!(result.steps.is_empty());
        assert_eq!(state.blowups.len(), 0);
        assert!(!state.singular_components().is_empty());
    }

    #[test]
    fn double_cover_branch_crepancy_accepts_expected_multiplicities() {
        assert!(!super::is_crepant_double_cover_branch_blowup(3, 2));
        assert!(super::is_crepant_double_cover_branch_blowup(3, 4));
        assert!(super::is_crepant_double_cover_branch_blowup(3, 5));
        assert!(super::is_crepant_double_cover_branch_blowup(2, 2));
        assert!(super::is_crepant_double_cover_branch_blowup(2, 3));
        assert!(!super::is_crepant_double_cover_branch_blowup(2, 4));
    }

    #[test]
    fn auto_crepant_resolution_uses_branch_transform_for_odd_multiplicity() {
        let mut state = AppState::new(2, "x1^5 + x2^5 + x3^5").unwrap();

        let result = state
            .resolve_crepant_with_progress(AutoResolutionOptions { max_steps: 1 }, |_| {})
            .unwrap();

        assert_eq!(result.status, AutoResolutionStatus::NoCrepantCenter);
        assert_eq!(result.steps.len(), 1);
        assert_eq!(result.steps[0].multiplicity, 5);
        let first_output_chart = state.chart(result.steps[0].output_chart_ids[0]).unwrap();
        assert_eq!(first_output_chart.polynomial.center_order(&[0]), Some(1));
    }

    #[test]
    fn auto_crepant_resolution_accepts_branch_double_and_triple_lines() {
        let mut double_line = AppState::new(2, "x0^2*x2 + x1^2*x3").unwrap();
        let double_result = double_line
            .resolve_crepant_with_progress(AutoResolutionOptions { max_steps: 1 }, |_| {})
            .unwrap();

        assert_eq!(double_result.steps.len(), 1);
        assert_eq!(double_result.steps[0].center_dimension, 1);
        assert_eq!(double_result.steps[0].ambient_codimension, 2);
        assert_eq!(double_result.steps[0].multiplicity, 2);

        let mut triple_line = AppState::new(2, "x0^3*x2 + x1^3*x3").unwrap();
        let triple_result = triple_line
            .resolve_crepant_with_progress(AutoResolutionOptions { max_steps: 1 }, |_| {})
            .unwrap();

        assert_eq!(triple_result.steps.len(), 1);
        assert_eq!(triple_result.steps[0].center_dimension, 1);
        assert_eq!(triple_result.steps[0].ambient_codimension, 2);
        assert_eq!(triple_result.steps[0].multiplicity, 3);
    }

    #[test]
    fn auto_crepant_resolution_reports_progress() {
        let mut state = AppState::new(2, "x1^4 + x2^4 + x3^4").unwrap();
        let mut events = Vec::new();

        let result = state
            .resolve_crepant_with_progress(AutoResolutionOptions { max_steps: 4 }, |event| {
                events.push(event);
            })
            .unwrap();

        assert_eq!(result.status, AutoResolutionStatus::Resolved);
        assert!(matches!(
            events.first(),
            Some(AutoResolutionProgress::CheckingSingularLocus {
                completed_steps: 0,
                max_steps: 4,
            })
        ));
        assert!(events.iter().any(|event| matches!(
            event,
            AutoResolutionProgress::BlowingUp {
                step: 1,
                max_steps: 4,
                ..
            }
        )));
        assert!(
            events
                .iter()
                .any(|event| matches!(event, AutoResolutionProgress::StepFinished { step: 1, .. }))
        );
        assert!(matches!(
            events.last(),
            Some(AutoResolutionProgress::SingularLocusFound {
                completed_steps: 1,
                components: 0,
                intersections: 0,
            })
        ));
    }

    #[test]
    fn deduplicates_projective_singular_line_across_charts() {
        let state = AppState::new(2, "x0^2*x2 + x1^2*x3").unwrap();
        let components = state.singular_components();

        assert_eq!(components.len(), 1);
        assert_eq!(components[0].dimension, 1);
        assert_eq!(components[0].multiplicity, 2);
        assert_eq!(components[0].charts.len(), 2);
        assert!(components[0].charts.iter().any(|chart| chart.chart_id == 2));
        assert!(components[0].charts.iter().any(|chart| chart.chart_id == 3));
    }

    #[test]
    fn detects_generic_ca_type_on_singular_line() {
        let state = AppState::new(3, "x0*x1*x4^7 + x2^9").unwrap();
        let components = state.singular_components();
        let line = components
            .iter()
            .find(|component| component.dimension == 1 && component.multiplicity == 2)
            .expect("singular line component");

        let singularity_type = line.generic_singularity_type.as_ref().unwrap();
        assert_eq!(singularity_type.label, "cA8");
        assert_eq!(line.generic_singularity_chart, Some(3));
    }

    #[test]
    fn classifies_supplied_seven_line_surface_locally() {
        let state = AppState::new(2, SEVEN_LINE_SURFACE).unwrap();
        let mut labels = state
            .singular_components()
            .into_iter()
            .map(|component| {
                component
                    .generic_singularity_type
                    .expect("generic line type")
                    .label
            })
            .collect::<Vec<_>>();
        labels.sort();

        assert_eq!(labels, ["cA1", "cA1", "cA3", "cA3", "cA3", "cA4", "cA5"]);
    }

    #[test]
    fn singular_locus_includes_generic_line_analysis_by_default() {
        let state = AppState::new(3, "x0*x1*x4^7 + x2^9").unwrap();
        let analysis = state.singular_locus();
        let line = analysis
            .components
            .iter()
            .find(|component| component.dimension == 1 && component.multiplicity == 2)
            .expect("singular line component");

        assert_eq!(
            line.generic_singularity_type
                .as_ref()
                .map(|singularity_type| singularity_type.multiplicity),
            Some(2)
        );
    }

    #[test]
    fn handles_high_degree_surface_with_linear_singular_locus() {
        let state = AppState::new(2, STALLING_SURFACE).unwrap();
        let analysis = state.singular_analysis();

        assert_eq!(
            analysis
                .components
                .iter()
                .filter(|component| component.dimension == 1)
                .count(),
            8
        );
        assert_eq!(
            analysis
                .components
                .iter()
                .filter(|component| component.dimension == 0)
                .count(),
            1
        );
        assert_eq!(analysis.intersections.len(), 10);
        assert!(
            analysis
                .components
                .iter()
                .filter(|component| component.dimension == 1)
                .all(|component| component.generic_singularity_type.is_some())
        );
    }

    #[test]
    fn handles_reported_intersection_blowup_without_component_duplication() {
        let mut state = AppState::new(2, REPORTED_INTERSECTION_SURFACE).unwrap();
        let before = state.singular_locus();

        assert_eq!(before.components.len(), 11);
        assert_eq!(before.intersections.len(), 9);

        state.blowup_intersection(9, false).unwrap();
        let after = state.singular_locus();

        assert_eq!(after.components.len(), 13);
        assert_eq!(after.intersections.len(), 17);
        assert_eq!(
            after
                .components
                .iter()
                .filter(|component| component.dimension == 1)
                .count(),
            13
        );
        assert!(after.components.iter().any(|component| {
            component.dimension == 1
                && component
                    .charts
                    .iter()
                    .filter(|chart| chart.chart_id >= 4)
                    .count()
                    >= 3
        }));
    }

    #[test]
    fn uses_affine_linear_fast_path_after_high_degree_blowup() {
        let mut state = AppState::new(2, STALLING_SURFACE).unwrap();
        state.blowup_singular_component(1, false).unwrap();
        let chart = state.chart(4).expect("first post-blowup chart");
        let equations = super::singularity_equations(chart);
        let components = super::local_singular_components(chart, &equations, false);

        assert_eq!(components.len(), 7);
        assert!(components.iter().all(|basis| {
            basis.len() >= 2
                && basis
                    .iter()
                    .all(|generator| super::affine_linear_row(generator).is_some())
        }));
    }

    #[test]
    fn finds_rational_intersection_points_of_singular_components() {
        let state = AppState::new(2, "x0*x1*x2").unwrap();
        let intersections = state.singular_component_intersections();

        assert_eq!(intersections.len(), 1);
        assert_eq!(intersections[0].chart_id, 3);
        assert_eq!(intersections[0].component_indices, vec![1, 2, 3]);
        let assignments = intersections[0].coordinate_assignments.as_ref().unwrap();
        assert_eq!(assignments.get("x0"), Some(&Rational::zero()));
        assert_eq!(assignments.get("x1"), Some(&Rational::zero()));
        assert_eq!(assignments.get("x2"), Some(&Rational::zero()));
        assert_eq!(
            intersections[0]
                .singularity_type
                .as_ref()
                .unwrap()
                .multiplicity,
            3
        );
    }

    #[test]
    fn deduplicates_projective_intersection_points_across_charts() {
        let state = AppState::new(2, "x0*x1*(x2-x3)").unwrap();
        let intersections = state.singular_component_intersections();

        assert_eq!(intersections.len(), 1);
        let assignments = intersections[0].coordinate_assignments.as_ref().unwrap();
        assert_eq!(assignments.get("x0"), Some(&Rational::zero()));
        assert_eq!(assignments.get("x1"), Some(&Rational::zero()));
    }

    #[test]
    fn blows_up_singular_component_intersection_point() {
        let mut state = AppState::new(2, "x0*x1*x2").unwrap();
        state.blowup_intersection(1, false).unwrap();

        assert_eq!(state.blowups.len(), 1);
        assert!(!state.chart(3).unwrap().active);
        assert_eq!(state.blowups[0].input_chart, 3);
        assert_eq!(state.blowups[0].output_charts.len(), 3);
    }

    #[test]
    fn blows_up_component_on_each_chart_appearance() {
        let mut state = AppState::new(2, "x0^2*x2 + x1^2*x3").unwrap();
        state.blowup_singular_component(1, false).unwrap();

        assert_eq!(state.blowups.len(), 2);
        assert!(!state.chart(2).unwrap().active);
        assert!(!state.chart(3).unwrap().active);
        assert_eq!(state.active_charts().count(), 6);
    }

    #[test]
    fn deduplicates_line_that_is_not_coordinate_in_every_chart() {
        let state = AppState::new(2, "(x0-x1)^2*x3 + x2^2*x0").unwrap();
        let components = state.singular_components();

        assert_eq!(components.len(), 1);
        assert_eq!(components[0].dimension, 1);
        assert_eq!(components[0].multiplicity, 2);
        assert_eq!(components[0].charts.len(), 3);
        assert!(components[0].charts.iter().any(|chart| {
            chart.chart_id == 3
                && chart.coordinate_assignments.is_none()
                && chart.affine_linear_center
        }));
    }

    #[test]
    fn blows_up_affine_linear_component_chart() {
        let mut state = AppState::new(2, "(x0-x1)^2*x3 + x2^2*x0").unwrap();
        state.blowup_singular_component(1, false).unwrap();

        assert_eq!(state.blowups.len(), 3);
        assert!(!state.chart(0).unwrap().active);
        assert!(!state.chart(1).unwrap().active);
        assert!(!state.chart(3).unwrap().active);
        assert_eq!(state.active_charts().count(), 7);
    }
}
