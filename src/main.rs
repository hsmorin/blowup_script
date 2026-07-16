mod cli;
mod grobner;
mod json;
mod model;
mod parser;
mod poly;
mod rational;

fn main() {
    let result = if std::env::args().nth(1).as_deref() == Some("local-type") {
        run_local_type()
    } else {
        cli::run()
    };
    if let Err(err) = result {
        eprintln!("error: {err}");
        std::process::exit(1);
    }
}

fn run_local_type() -> Result<(), String> {
    let args = std::env::args().skip(2).collect::<Vec<_>>();
    let [variables, coordinates] = args.as_slice() else {
        return Err(
            "usage: blowup-script local-type VAR1,VAR2,... VALUE1,VALUE2,... < polynomial"
                .to_string(),
        );
    };
    let variables = variables
        .split(',')
        .map(str::trim)
        .map(str::to_string)
        .collect::<Vec<_>>();
    let coordinates = coordinates
        .split(',')
        .map(str::trim)
        .map(rational::Rational::parse)
        .collect::<Result<Vec<_>, _>>()?;
    let mut polynomial = String::new();
    std::io::Read::read_to_string(&mut std::io::stdin(), &mut polynomial)
        .map_err(|err| format!("failed to read polynomial: {err}"))?;
    let singularity = model::analyze_affine_polynomial_at(&polynomial, &variables, &coordinates)?;
    println!("label: {}", singularity.label);
    println!("multiplicity: {}", singularity.multiplicity);
    if let Some(rank) = singularity.quadratic_rank {
        println!("quadratic rank: {rank}");
    }
    Ok(())
}
