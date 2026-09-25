use rust_htslib::bam::{self as htsbam, Format, Header, Read, Reader, Writer};
use std::path::{Path, PathBuf};
use tempfile::{NamedTempFile, TempPath};

pub struct PreparedAlignment {
    path: PathBuf,
    _temporary: Option<TempPath>,
}

impl PreparedAlignment {
    pub fn path(&self) -> &Path {
        &self.path
    }
}

fn looks_native(path: &str, format_hint: &str) -> bool {
    if path == "-" {
        return false;
    }
    match format_hint {
        "sam" => true,
        "bam" => true,
        "auto" => {
            let lower = path.to_ascii_lowercase();
            lower.ends_with(".sam") || lower.ends_with(".bam")
        }
        _ => false,
    }
}

pub fn prepare_alignment(
    path: &str,
    threads: u16,
    format_hint: &str,
) -> Result<PreparedAlignment, String> {
    // Keep the fast pure-Rust BAM/SAM path for ordinary files. HTSlib is used
    // only when we need stdin, CRAM, or content-based format autodetection.
    if looks_native(path, format_hint) {
        return Ok(PreparedAlignment {
            path: PathBuf::from(path),
            _temporary: None,
        });
    }

    let mut reader = if path == "-" {
        Reader::from_stdin().map_err(|e| format!("Could not read alignments from stdin: {e}"))?
    } else {
        Reader::from_path(path)
            .map_err(|e| format!("Could not open alignment file '{path}': {e}"))?
    };
    if threads > 1 {
        reader
            .set_threads((threads - 1) as usize)
            .map_err(|e| format!("Could not enable HTSlib reader threads: {e}"))?;
    }

    let header = Header::from_template(reader.header());
    let temp = NamedTempFile::new()
        .map_err(|e| format!("Could not create temporary BAM for input conversion: {e}"))?;
    let temp_path = temp.into_temp_path();

    {
        let mut writer = Writer::from_path(&temp_path, &header, Format::Bam)
            .map_err(|e| format!("Could not create temporary BAM: {e}"))?;
        if threads > 1 {
            writer
                .set_threads((threads - 1) as usize)
                .map_err(|e| format!("Could not enable HTSlib writer threads: {e}"))?;
        }

        for record in reader.records() {
            let record =
                record.map_err(|e| format!("Error while decoding input alignment: {e}"))?;
            writer
                .write(&record)
                .map_err(|e| format!("Error while converting input alignment: {e}"))?;
        }
    }

    Ok(PreparedAlignment {
        path: temp_path.to_path_buf(),
        _temporary: Some(temp_path),
    })
}

pub fn convert_alignment_output(
    sam_path: &Path,
    output_path: &Path,
    format: &str,
    threads: u16,
) -> Result<(), String> {
    let target_format = match format.to_ascii_lowercase().as_str() {
        "sam" => {
            std::fs::copy(sam_path, output_path)
                .map_err(|e| format!("Could not copy SAM output: {e}"))?;
            return Ok(());
        }
        "bam" => Format::Bam,
        other => return Err(format!("Unsupported --samout-format '{other}'")),
    };

    let mut reader =
        Reader::from_path(sam_path).map_err(|e| format!("Could not reopen SAM output: {e}"))?;
    let header = Header::from_template(reader.header());
    let mut writer = Writer::from_path(output_path, &header, target_format)
        .map_err(|e| format!("Could not create BAM output: {e}"))?;

    if threads > 1 {
        reader
            .set_threads((threads - 1) as usize)
            .map_err(|e| format!("Could not enable HTSlib reader threads: {e}"))?;
        writer
            .set_threads((threads - 1) as usize)
            .map_err(|e| format!("Could not enable HTSlib writer threads: {e}"))?;
    }

    for record in reader.records() {
        let record = record.map_err(|e| format!("Could not read SAM output: {e}"))?;
        writer
            .write(&record)
            .map_err(|e| format!("Could not write BAM output: {e}"))?;
    }
    Ok(())
}
