use rust_hdf5::{H5Dataset, H5File, HBool, VarLenUnicode};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct OutputTable {
    /// Number of annotation-derived rows before HTSeq's special __... rows.
    pub real_feature_count: usize,
    pub feature_ids: Vec<String>,
    pub metadata_names: Vec<String>,
    pub metadata_values: Vec<Vec<String>>,
    pub sample_names: Vec<String>,
    /// Sample-major matrix with shape [n_samples, n_features].
    pub values: Vec<f64>,
}

impl OutputTable {
    pub fn n_features(&self) -> usize {
        self.feature_ids.len()
    }

    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }

    pub fn value(&self, sample: usize, feature: usize) -> f64 {
        self.values[sample * self.n_features() + feature]
    }
}

pub fn write_output(
    table: &OutputTable,
    output_path: Option<&str>,
    delimiter: &str,
    append: bool,
    with_header: bool,
    sparse: bool,
) -> Result<(), String> {
    match output_path {
        None => {
            let stdout = std::io::stdout();
            let mut out = stdout.lock();
            write_tabular(&mut out, table, delimiter, with_header)
                .map_err(|e| e.to_string())
        }
        Some(path) => {
            let suffix = Path::new(path)
                .extension()
                .and_then(|ext| ext.to_str())
                .unwrap_or("")
                .to_ascii_lowercase();

            match suffix.as_str() {
                "" | "tsv" | "csv" | "txt" => {
                    let mut options = OpenOptions::new();
                    options.create(true).write(true);
                    if append {
                        options.append(true);
                    } else {
                        options.truncate(true);
                    }
                    let file = options.open(path).map_err(|e| e.to_string())?;
                    let mut out = BufWriter::new(file);
                    write_tabular(&mut out, table, delimiter, with_header)
                        .map_err(|e| e.to_string())
                }
                "mtx" => write_mtx(path, table, sparse),
                "h5ad" => write_h5ad(path, table, sparse),
                "loom" => write_loom(path, table),
                other => Err(format!(
                    "Format not recognized for output count file: {other}"
                )),
            }
        }
    }
}

fn write_tabular<W: Write>(
    out: &mut W,
    table: &OutputTable,
    delimiter: &str,
    with_header: bool,
) -> std::io::Result<()> {
    if with_header {
        write!(out, "{}", delimiter)?;
        for _ in &table.metadata_names {
            write!(out, "{}", delimiter)?;
        }
        for (index, sample) in table.sample_names.iter().enumerate() {
            if index != 0 {
                write!(out, "{}", delimiter)?;
            }
            write!(out, "{}", sample)?;
        }
        writeln!(out)?;
    }

    for feature in 0..table.n_features() {
        write!(out, "{}", table.feature_ids[feature])?;
        for value in &table.metadata_values[feature] {
            write!(out, "{}{}", delimiter, value)?;
        }
        for sample in 0..table.n_samples() {
            write!(out, "{}{}", delimiter, table.value(sample, feature))?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn sidecar_prefix(path: &str) -> String {
    path.strip_suffix(".mtx").unwrap_or(path).to_string()
}

fn write_mtx(path: &str, table: &OutputTable, sparse: bool) -> Result<(), String> {
    let mut out = BufWriter::new(File::create(path).map_err(|e| e.to_string())?);
    if sparse {
        let nnz = table.values.iter().filter(|value| **value != 0.0).count();
        writeln!(out, "%%MatrixMarket matrix coordinate real general")
            .map_err(|e| e.to_string())?;
        writeln!(out, "{} {} {}", table.n_samples(), table.n_features(), nnz)
            .map_err(|e| e.to_string())?;
        for sample in 0..table.n_samples() {
            for feature in 0..table.n_features() {
                let value = table.value(sample, feature) as f32;
                if value != 0.0 {
                    writeln!(out, "{} {} {}", sample + 1, feature + 1, value)
                        .map_err(|e| e.to_string())?;
                }
            }
        }
    } else {
        writeln!(out, "%%MatrixMarket matrix array real general")
            .map_err(|e| e.to_string())?;
        writeln!(out, "{} {}", table.n_samples(), table.n_features())
            .map_err(|e| e.to_string())?;
        // Matrix Market array format is column-major.
        for feature in 0..table.n_features() {
            for sample in 0..table.n_samples() {
                writeln!(out, "{}", table.value(sample, feature) as f32)
                    .map_err(|e| e.to_string())?;
            }
        }
    }

    let prefix = sidecar_prefix(path);
    let mut samples = BufWriter::new(
        File::create(format!("{prefix}_samples.tsv")).map_err(|e| e.to_string())?,
    );
    for sample in &table.sample_names {
        writeln!(samples, "{}", sample).map_err(|e| e.to_string())?;
    }

    let mut features = BufWriter::new(
        File::create(format!("{prefix}_features.tsv")).map_err(|e| e.to_string())?,
    );
    write!(features, "id").map_err(|e| e.to_string())?;
    for name in &table.metadata_names {
        write!(features, "\t{}", name).map_err(|e| e.to_string())?;
    }
    writeln!(features).map_err(|e| e.to_string())?;
    for (index, id) in table.feature_ids.iter().enumerate() {
        write!(features, "{}", id).map_err(|e| e.to_string())?;
        for value in &table.metadata_values[index] {
            write!(features, "\t{}", value).map_err(|e| e.to_string())?;
        }
        writeln!(features).map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn set_dataset_string_attr(
    dataset: &H5Dataset,
    name: &str,
    value: &str,
) -> Result<(), String> {
    let attr = dataset
        .new_attr::<VarLenUnicode>()
        .shape(())
        .create(name)
        .map_err(|e| e.to_string())?;
    attr.write_string(value).map_err(|e| e.to_string())
}

fn write_string_array(
    group: &rust_hdf5::H5Group,
    name: &str,
    values: &[String],
) -> Result<H5Dataset, String> {
    let refs: Vec<&str> = values.iter().map(String::as_str).collect();
    let dataset = group
        .write_vlen_strings(name, &refs)
        .map_err(|e| e.to_string())?;
    set_dataset_string_attr(&dataset, "encoding-type", "string-array")?;
    set_dataset_string_attr(&dataset, "encoding-version", "0.2.0")?;
    Ok(dataset)
}

fn write_categorical_metadata(
    group: &rust_hdf5::H5Group,
    name: &str,
    values: &[String],
    valid_count: usize,
) -> Result<(), String> {
    use std::collections::BTreeSet;

    let categories: Vec<String> = values
        .iter()
        .take(valid_count)
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    let category_index: std::collections::HashMap<&str, i32> = categories
        .iter()
        .enumerate()
        .map(|(index, value)| (value.as_str(), index as i32))
        .collect();

    let mut codes = Vec::with_capacity(values.len());
    for (index, value) in values.iter().enumerate() {
        if index >= valid_count {
            codes.push(-1);
        } else {
            codes.push(*category_index.get(value.as_str()).unwrap());
        }
    }

    let categorical = group.create_group(name).map_err(|e| e.to_string())?;
    categorical
        .set_attr_string("encoding-type", "categorical")
        .map_err(|e| e.to_string())?;
    categorical
        .set_attr_string("encoding-version", "0.2.0")
        .map_err(|e| e.to_string())?;
    categorical
        .set_attr_numeric("ordered", &HBool::from(false))
        .map_err(|e| e.to_string())?;

    write_string_array(&categorical, "categories", &categories)?;

    let codes_ds = categorical
        .new_dataset::<i32>()
        .shape(&[codes.len()])
        .create("codes")
        .map_err(|e| e.to_string())?;
    codes_ds.write_raw(&codes).map_err(|e| e.to_string())?;
    set_dataset_string_attr(&codes_ds, "encoding-type", "array")?;
    set_dataset_string_attr(&codes_ds, "encoding-version", "0.2.0")?;
    Ok(())
}

fn write_h5ad(path: &str, table: &OutputTable, sparse: bool) -> Result<(), String> {
    let file = H5File::create(path).map_err(|e| e.to_string())?;
    file.set_attr_string("encoding-type", "anndata")
        .map_err(|e| e.to_string())?;
    file.set_attr_string("encoding-version", "0.1.0")
        .map_err(|e| e.to_string())?;

    if sparse {
        let mut data = Vec::<f32>::new();
        let mut indices = Vec::<i32>::new();
        let mut indptr = Vec::<i32>::with_capacity(table.n_samples() + 1);
        indptr.push(0);

        for sample in 0..table.n_samples() {
            for feature in 0..table.n_features() {
                let value = table.value(sample, feature);
                if value != 0.0 {
                    data.push(value as f32);
                    indices.push(feature as i32);
                }
            }
            indptr.push(data.len() as i32);
        }

        let x = file.create_group("X").map_err(|e| e.to_string())?;
        x.set_attr_string("encoding-type", "csr_matrix")
            .map_err(|e| e.to_string())?;
        x.set_attr_string("encoding-version", "0.1.0")
            .map_err(|e| e.to_string())?;
        x.set_attr_array_numeric(
            "shape",
            &[table.n_samples() as i64, table.n_features() as i64],
        )
        .map_err(|e| e.to_string())?;

        let data_ds = x
            .new_dataset::<f32>()
            .shape(&[data.len()])
            .create("data")
            .map_err(|e| e.to_string())?;
        data_ds.write_raw(&data).map_err(|e| e.to_string())?;
        set_dataset_string_attr(&data_ds, "encoding-type", "array")?;
        set_dataset_string_attr(&data_ds, "encoding-version", "0.2.0")?;

        let indices_ds = x
            .new_dataset::<i32>()
            .shape(&[indices.len()])
            .create("indices")
            .map_err(|e| e.to_string())?;
        indices_ds
            .write_raw(&indices)
            .map_err(|e| e.to_string())?;
        set_dataset_string_attr(&indices_ds, "encoding-type", "array")?;
        set_dataset_string_attr(&indices_ds, "encoding-version", "0.2.0")?;

        let indptr_ds = x
            .new_dataset::<i32>()
            .shape(&[indptr.len()])
            .create("indptr")
            .map_err(|e| e.to_string())?;
        indptr_ds.write_raw(&indptr).map_err(|e| e.to_string())?;
        set_dataset_string_attr(&indptr_ds, "encoding-type", "array")?;
        set_dataset_string_attr(&indptr_ds, "encoding-version", "0.2.0")?;
    } else {
        let matrix_f32: Vec<f32> =
            table.values.iter().map(|value| *value as f32).collect();
        let x = file
            .new_dataset::<f32>()
            .shape(&[table.n_samples(), table.n_features()])
            .create("X")
            .map_err(|e| e.to_string())?;
        x.write_raw(&matrix_f32).map_err(|e| e.to_string())?;
        set_dataset_string_attr(&x, "encoding-type", "array")?;
        set_dataset_string_attr(&x, "encoding-version", "0.2.0")?;
    }

    let obs = file.create_group("obs").map_err(|e| e.to_string())?;
    obs.set_attr_string("_index", "_index")
        .map_err(|e| e.to_string())?;
    obs.set_attr_string_array("column-order", &[])
        .map_err(|e| e.to_string())?;
    obs.set_attr_string("encoding-type", "dataframe")
        .map_err(|e| e.to_string())?;
    obs.set_attr_string("encoding-version", "0.2.0")
        .map_err(|e| e.to_string())?;
    write_string_array(&obs, "_index", &table.sample_names)?;

    let var = file.create_group("var").map_err(|e| e.to_string())?;
    var.set_attr_string("_index", "_index")
        .map_err(|e| e.to_string())?;
    let metadata_refs: Vec<&str> =
        table.metadata_names.iter().map(String::as_str).collect();
    var.set_attr_string_array("column-order", &metadata_refs)
        .map_err(|e| e.to_string())?;
    var.set_attr_string("encoding-type", "dataframe")
        .map_err(|e| e.to_string())?;
    var.set_attr_string("encoding-version", "0.2.0")
        .map_err(|e| e.to_string())?;
    write_string_array(&var, "_index", &table.feature_ids)?;
    for (column, name) in table.metadata_names.iter().enumerate() {
        let values: Vec<String> = table
            .metadata_values
            .iter()
            .map(|row| row.get(column).cloned().unwrap_or_default())
            .collect();
        // HTSeq builds a pandas DataFrame whose metadata columns are shorter
        // than the ID column by the five special __... rows. AnnData's default
        // writer converts these string columns to categoricals and represents
        // the missing special-row values using categorical code -1.
        write_categorical_metadata(
            &var,
            name,
            &values,
            table.real_feature_count,
        )?;
    }

    for name in ["obsm", "varm", "obsp", "varp", "layers", "uns"] {
        let group = file.create_group(name).map_err(|e| e.to_string())?;
        group
            .set_attr_string("encoding-type", "dict")
            .map_err(|e| e.to_string())?;
        group
            .set_attr_string("encoding-version", "0.1.0")
            .map_err(|e| e.to_string())?;
    }

    file.close().map_err(|e| e.to_string())
}

fn write_loom(path: &str, table: &OutputTable) -> Result<(), String> {
    let file = H5File::create(path).map_err(|e| e.to_string())?;
    file.set_attr_string("LOOM_SPEC_VERSION", "3.0.0")
        .map_err(|e| e.to_string())?;

    // Loom stores features as rows and samples as columns.
    let mut matrix = vec![0.0f32; table.n_features() * table.n_samples()];
    for feature in 0..table.n_features() {
        for sample in 0..table.n_samples() {
            matrix[feature * table.n_samples() + sample] =
                table.value(sample, feature) as f32;
        }
    }
    let dataset = file
        .new_dataset::<f32>()
        .shape(&[table.n_features(), table.n_samples()])
        .create("matrix")
        .map_err(|e| e.to_string())?;
    dataset.write_raw(&matrix).map_err(|e| e.to_string())?;

    let row_attrs = file.create_group("row_attrs").map_err(|e| e.to_string())?;
    let ids: Vec<&str> = table.feature_ids.iter().map(String::as_str).collect();
    row_attrs
        .write_vlen_strings("id", &ids)
        .map_err(|e| e.to_string())?;
    for (column, name) in table.metadata_names.iter().enumerate() {
        let values: Vec<String> = table
            .metadata_values
            .iter()
            .map(|row| row.get(column).cloned().unwrap_or_default())
            .collect();
        let refs: Vec<&str> = values.iter().map(String::as_str).collect();
        row_attrs
            .write_vlen_strings(name, &refs)
            .map_err(|e| e.to_string())?;
    }

    let col_attrs = file.create_group("col_attrs").map_err(|e| e.to_string())?;
    let samples: Vec<&str> = table.sample_names.iter().map(String::as_str).collect();
    col_attrs
        .write_vlen_strings("_index", &samples)
        .map_err(|e| e.to_string())?;

    for group in ["attrs", "layers", "row_graphs", "col_graphs"] {
        file.create_group(group).map_err(|e| e.to_string())?;
    }

    file.close().map_err(|e| e.to_string())
}
