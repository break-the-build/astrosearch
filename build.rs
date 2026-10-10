//! Embed the git revision in the binary so provenance records can name the
//! exact source that produced an output.

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8_lossy(&output.stdout).trim().to_string())
}

fn main() {
    let commit = git(&["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".into());
    let dirty = match git(&["status", "--porcelain", "--untracked-files=no"]) {
        Some(status) if status.is_empty() => "false",
        Some(_) => "true",
        None => "unknown",
    };
    println!("cargo:rustc-env=ASTROSEARCH_GIT_COMMIT={commit}");
    println!("cargo:rustc-env=ASTROSEARCH_GIT_DIRTY={dirty}");
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/index");
    println!("cargo:rerun-if-changed=.git/refs");
}
