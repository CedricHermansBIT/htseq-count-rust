use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use tempfile::TempPath;

#[cfg(unix)]
use rust_htslib::bam::{Format, Header, Read as HtsRead, Reader, Writer};

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

#[cfg(unix)]
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
    let temp = tempfile::Builder::new()
        .suffix(".bam")
        .tempfile()
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

#[cfg(windows)]
fn write_sniffed_temporary(bytes: &[u8]) -> Result<PreparedAlignment, String> {
    if bytes.starts_with(b"CRAM") {
        return Err(
            "CRAM input is not available in the Windows build because upstream hts-sys              does not currently build natively on Windows. Use SAM/BAM on Windows, or              use the Linux/macOS TallySeq build for CRAM."
                .to_string(),
        );
    }

    // BAM is BGZF and therefore begins with the gzip magic bytes. Ordinary
    // SAM is text. A suffix is enough for the native reader selected later.
    let suffix = if bytes.starts_with(&[0x1f, 0x8b]) {
        ".bam"
    } else {
        ".sam"
    };
    let mut temp = tempfile::Builder::new()
        .suffix(suffix)
        .tempfile()
        .map_err(|e| format!("Could not create temporary alignment file: {e}"))?;
    temp.write_all(bytes)
        .map_err(|e| format!("Could not buffer alignment input: {e}"))?;
    temp.flush()
        .map_err(|e| format!("Could not flush temporary alignment input: {e}"))?;
    let temp_path = temp.into_temp_path();

    Ok(PreparedAlignment {
        path: temp_path.to_path_buf(),
        _temporary: Some(temp_path),
    })
}

#[cfg(windows)]
pub fn prepare_alignment(
    path: &str,
    _threads: u16,
    format_hint: &str,
) -> Result<PreparedAlignment, String> {
    if looks_native(path, format_hint) {
        return Ok(PreparedAlignment {
            path: PathBuf::from(path),
            _temporary: None,
        });
    }

    if path == "-" {
        let mut bytes = Vec::new();
        std::io::stdin()
            .read_to_end(&mut bytes)
            .map_err(|e| format!("Could not read alignments from stdin: {e}"))?;
        return write_sniffed_temporary(&bytes);
    }

    let bytes = std::fs::read(path)
        .map_err(|e| format!("Could not open alignment file '{path}': {e}"))?;
    write_sniffed_temporary(&bytes)
}
