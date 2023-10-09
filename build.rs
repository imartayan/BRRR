/// This is a hack to support dynamic K values.
/// K values are implemented as a const generic in our code
/// as we expect it to remain constant across executions
/// and benefit from compile-time optimizations.
/// This build script will set the value of K at compile-time
/// from an environment variable, so one can easily build
/// the project "just in time" with the desired K value.
/// This will not re-build if the K value does not change.
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=K");

    let out_dir: std::path::PathBuf = std::env::var("OUT_DIR")
        .expect("Failed to obtain OUT_DIR")
        .into();
    let mut code = Vec::new();

    let k: usize = std::env::var("K")
        .unwrap_or_else(|_| "31".into())
        .parse()
        .expect("Failed to parse K");
    assert!(k % 2 == 1, "K must be odd!");
    code.push(format!("pub const K: usize = {k};"));

    let kmer_bits = 2 * k;
    code.push(format!("pub const KMER_BITS: usize = {kmer_bits};"));

    let kt = select_type(kmer_bits);
    code.push(format!("pub type KT = {kt};"));

    let m: usize = std::env::var("M")
        .unwrap_or_else(|_| "21".into())
        .parse()
        .expect("Failed to parse M");
    assert!(m % 2 == 1, "M must be odd!");
    code.push(format!("pub const M: usize = {m};"));

    let mmer_bits = 2 * m;
    code.push(format!("pub const MMER_BITS: usize = {mmer_bits};"));

    let mt = select_type(mmer_bits);
    code.push(format!("pub type MT = {mt};"));

    std::fs::write(out_dir.join("constants.rs"), code.join("\n"))
        .expect("Failed to write const file");
}

fn select_type(n_bits: usize) -> &'static str {
    match n_bits.next_power_of_two() {
        8 => "u8",
        16 => "u16",
        32 => "u32",
        64 => "u64",
        _ => "u128",
    }
}
