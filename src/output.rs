//! CSV rendering of candidates.

use crate::Candidate;

/// CSV includes both display order and rank within each detector class.
pub fn render_csv(candidates: &[Candidate]) -> String {
    render(candidates, None)
}

/// `render_csv` plus `jackknife_score` and `tier` columns: transit candidates
/// whose jackknife score (`Candidate::jackknife_score`) reaches
/// `detection_threshold` are `detection`, other transit candidates
/// `screening`; non-transit classes are left blank.
pub fn render_csv_with_tier(candidates: &[Candidate], detection_threshold: f64) -> String {
    render(candidates, Some(detection_threshold))
}

fn tier(candidate: &Candidate, threshold: f64) -> &'static str {
    if !candidate.kind.is_transit() {
        ""
    } else if candidate.jackknife_score() >= threshold {
        "detection"
    } else {
        "screening"
    }
}

fn render(candidates: &[Candidate], detection_threshold: Option<f64>) -> String {
    use std::collections::BTreeMap;
    use std::fmt::Write;
    let mut output = String::from(
        "rank,class,class_rank,score,period_days,epoch,duration_days,depth_ppm,amplitude_ppm,n_events,odd_even_sigma,white_score,noise_factor,event_evidence,alias_ephemerides,note",
    );
    output.push_str(if detection_threshold.is_some() {
        ",jackknife_score,tier\n"
    } else {
        "\n"
    });
    let mut ranks = BTreeMap::new();
    for (i, c) in candidates.iter().enumerate() {
        let rank = ranks.entry(c.kind).or_insert(0);
        *rank += 1;
        let optional = |v: Option<f64>| v.map(|v| format!("{v:.3}")).unwrap_or_default();
        write!(
            output,
            "{},{},{},{:.5},{:.8},{:.8},{:.8},{:.4},{:.4},{},{},{:.5},{},\"{}\",\"{}\",\"{}\"",
            i + 1,
            c.kind,
            rank,
            c.score,
            c.period,
            c.epoch,
            c.duration,
            c.depth,
            c.amplitude,
            c.n_events,
            optional(c.odd_even_sigma),
            c.white_score,
            optional(c.noise_factor),
            c.event_evidence
                .iter()
                .map(|e| format!(
                    "{:.8}:{:.4}:{:.3}:{}",
                    e.epoch, e.depth, e.significance, e.points
                ))
                .collect::<Vec<_>>()
                .join(";"),
            c.alias_ephemerides
                .iter()
                .map(|a| format!(
                    "{:.8}:{:.8}:{:.8}:{:.3}:{}:{}:{}",
                    a.period,
                    a.epoch,
                    a.duration,
                    a.score,
                    a.n_events,
                    a.event_epochs
                        .iter()
                        .map(|epoch| format!("{epoch:.8}"))
                        .collect::<Vec<_>>()
                        .join("+"),
                    a.alias_type
                ))
                .collect::<Vec<_>>()
                .join(";"),
            c.note.replace('"', "\"\"")
        )
        .unwrap();
        if let Some(threshold) = detection_threshold {
            write!(output, ",{:.5},{}", c.jackknife_score(), tier(c, threshold)).unwrap();
        }
        output.push('\n');
    }
    output
}
