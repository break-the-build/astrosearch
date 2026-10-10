use crate::{math::median, Point};
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

/// Quoted fields, doubled quotes, empty fields, and CRLF are supported.
/// Unterminated or misplaced quotes are reported instead of shifting columns.
fn fields(line: &str, delimiter: char) -> Result<Vec<String>, String> {
    let mut output = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars().peekable();
    let mut quoted = false;
    let mut closed = false;
    while let Some(ch) = chars.next() {
        if quoted {
            if ch == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    field.push('"');
                } else {
                    quoted = false;
                    closed = true;
                }
            } else {
                field.push(ch);
            }
        } else if ch == delimiter {
            output.push(field.trim().to_string());
            field.clear();
            closed = false;
        } else if ch == '"' {
            if !field.trim().is_empty() || closed {
                return Err("misplaced quote in delimited record".into());
            }
            field.clear();
            quoted = true;
        } else if closed && !ch.is_whitespace() {
            return Err("unexpected text after a quoted field".into());
        } else if !closed {
            field.push(ch);
        }
    }
    if quoted {
        return Err("unterminated quoted field; multiline CSV records are unsupported".into());
    }
    output.push(field.trim().to_string());
    Ok(output)
}

pub(crate) fn read_points(path: &str) -> Result<(Vec<Point>, usize), String> {
    let file = File::open(path).map_err(|e| format!("open {path}: {e}"))?;
    parse_reader(BufReader::new(file))
}

fn parse_reader(reader: impl BufRead) -> Result<(Vec<Point>, usize), String> {
    let mut lines = reader.lines().enumerate();
    let (header_line, header) = loop {
        let (index, line) = lines.next().ok_or("empty input")?;
        let line = line.map_err(|e| e.to_string())?;
        let line = line.trim_start_matches('\u{feff}').trim().to_string();
        if !line.is_empty() && !line.starts_with('#') {
            break (index + 1, line);
        }
    };
    let delimiter = if header.contains('\t') { '\t' } else { ',' };
    let names: Vec<String> = fields(&header, delimiter)
        .map_err(|e| format!("line {header_line}: {e}"))?
        .into_iter()
        .map(|v| v.to_ascii_lowercase())
        .collect();
    // Preference follows the supplied key order, rather than CSV column order.
    let column = |keys: &[&str]| {
        keys.iter()
            .find_map(|key| names.iter().position(|v| v == key))
    };
    let time =
        column(&["time", "btjd", "bjd", "t"]).ok_or("header needs a time/BTJD/BJD column")?;
    let flux = column(&["flux", "pdcsap_flux", "sap_flux", "relative_flux"])
        .ok_or("header needs a flux column")?;
    let paired_error = match names[flux].as_str() {
        "pdcsap_flux" => column(&["pdcsap_flux_err"]),
        "sap_flux" => column(&["sap_flux_err"]),
        _ => None,
    };
    let error =
        paired_error.or_else(|| column(&["error", "flux_err", "flux_error", "uncertainty"]));
    let quality = column(&["quality", "quality_flag"]);
    let mut points = Vec::new();
    let mut skipped = 0;
    for (index, line) in lines {
        let line = line.map_err(|e| e.to_string())?;
        if line.trim().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        let record = fields(&line, delimiter).map_err(|e| format!("line {}: {e}", index + 1))?;
        let number = |column: usize| record.get(column).and_then(|v| v.parse::<f64>().ok());
        if quality.is_some_and(|q| number(q) != Some(0.0)) {
            skipped += 1;
            continue;
        }
        match (number(time), number(flux)) {
            (Some(time), Some(flux)) if time.is_finite() && flux.is_finite() && flux > 0.0 => {
                let error = error
                    .and_then(number)
                    .filter(|v| v.is_finite() && *v > 0.0)
                    .unwrap_or(0.0);
                points.push(Point { time, flux, error });
            }
            _ => skipped += 1,
        }
    }
    Ok((points, skipped))
}

pub(crate) fn clean_points(points: Vec<Point>) -> Vec<Point> {
    let mut points: Vec<Point> = points
        .into_iter()
        .filter(|p| p.time.is_finite() && p.flux.is_finite() && p.flux > 0.0)
        .map(|p| Point {
            error: if p.error.is_finite() && p.error > 0.0 {
                p.error
            } else {
                0.0
            },
            ..p
        })
        .collect();
    points.sort_by(|a, b| a.time.total_cmp(&b.time));
    let mut merged = Vec::with_capacity(points.len());
    let mut i = 0;
    while i < points.len() {
        let mut j = i + 1;
        while j < points.len() && points[j].time == points[i].time {
            j += 1;
        }
        merged.push(Point {
            time: points[i].time,
            flux: median(points[i..j].iter().map(|p| p.flux)),
            // Conservatively avoid claiming independent precision for duplicates.
            error: median(points[i..j].iter().map(|p| p.error).filter(|e| *e > 0.0)),
        });
        i = j;
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quoted_csv_uses_preferred_flux_and_matching_uncertainty() {
        let data = "\u{feff}# comment\nSAP_FLUX,SAP_FLUX_ERR,TIME,PDCSAP_FLUX,PDCSAP_FLUX_ERR,QUALITY,label\n999,90,1,10,0.2,\"0\",\"source, A\"\n999,90,2,11,0.3,1,x\n";
        let (points, skipped) = parse_reader(data.as_bytes()).unwrap();
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].flux, 10.0);
        assert_eq!(points[0].error, 0.2);
        assert_eq!(skipped, 1);
        assert_eq!(fields("\"a\"\"b\",,c", ',').unwrap(), ["a\"b", "", "c"]);
    }

    #[test]
    fn invalid_rows_and_uncertainties_do_not_propagate_nan_or_infinity() {
        let data = "time,flux,error\n1,1,inf\n2,-1,0.1\n3,0,0.1\n4,nan,0.1\n5,1,-1\n";
        let (points, skipped) = parse_reader(data.as_bytes()).unwrap();
        assert_eq!(skipped, 3);
        assert_eq!(points.len(), 2);
        assert!(points.iter().all(|p| p.error == 0.0));
        assert!(parse_reader("time,flux\n1,\"1\n".as_bytes()).is_err());
    }

    #[test]
    fn distinct_nearby_times_survive_and_exact_duplicates_merge() {
        let points = clean_points(vec![
            Point {
                time: 1.0,
                flux: 2.0,
                error: 0.1,
            },
            Point {
                time: 1.0,
                flux: 4.0,
                error: 0.1,
            },
            Point {
                time: 1.0 + 5e-10,
                flux: 6.0,
                error: 0.1,
            },
        ]);
        assert_eq!(points.len(), 2);
        assert_eq!(points[0].flux, 3.0);
    }
}
