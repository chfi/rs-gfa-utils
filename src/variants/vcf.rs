use bstr::BString;
use std::{
    fmt,
    fmt::{Display, Formatter},
    path::{Path, PathBuf},
};

use chrono::prelude::*;

/// A struct that holds Variants, as defined in the VCF format
#[derive(Debug, PartialEq)]
pub struct VCFRecord {
    pub chromosome: BString,
    pub position: i64,
    pub id: Option<BString>,
    pub reference: BString,
    pub alternate: Option<BString>,
    pub quality: Option<i32>,
    pub filter: Option<BString>,
    pub info: Option<BString>,
    pub format: Option<BString>,
    pub sample_name: Option<BString>,
}

impl VCFRecord {
    pub fn vcf_cmp(&self, other: &VCFRecord) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        let chr_cmp = self.chromosome.cmp(&other.chromosome);
        if let Ordering::Equal = chr_cmp {
            self.position.cmp(&other.position)
        } else {
            chr_cmp
        }
    }
}

impl Display for VCFRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        fn display_field<T: Display>(field: Option<T>) -> String {
            if let Some(x) = field {
                x.to_string()
            } else {
                ".".to_string()
            }
        }

        write!(f, "{}\t", self.chromosome)?;
        write!(f, "{}\t", self.position)?;
        write!(f, "{}\t", display_field(self.id.as_ref()))?;
        write!(f, "{}\t", self.reference)?;
        write!(f, "{}\t", display_field(self.alternate.as_ref()))?;
        write!(f, "{}\t", display_field(self.quality.as_ref()))?;
        write!(f, "{}\t", display_field(self.filter.as_ref()))?;
        write!(f, "{}", display_field(self.info.as_ref()))?;
        if let Some(format) = self.format.as_ref() {
            if let Some(sample) = self.sample_name.as_ref() {
                write!(f, "\t{}", format)?;
                write!(f, "\t{}", sample)?;
            }
        }
        Ok(())
    }
}

pub struct VCFHeader {
    reference: PathBuf,
}

impl VCFHeader {
    pub fn new<T: AsRef<Path>>(path: T) -> Self {
        let reference = path.as_ref().to_owned();
        Self { reference }
    }
}

impl Display for VCFHeader {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let date: DateTime<Utc> = Utc::now();

        writeln!(f, "##fileformat=VCFv4.2")?;
        writeln!(f, "##fileDate={}", date.format("%Y%m%d"))?;
        writeln!(f, "##reference={}", self.reference.display())?;

        writeln!(
            f,
            r#"##INFO=<ID=TYPE,Number=A,Type=String,Description="Type of each allele (snv, ins, del, mnp, clumped)">"#
        )?;

        // writeln!(
        //     f,
        //     r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#
        // )?;

        let header_line: BString = bstr::join(
            "\t",
            [
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                "INFO",
                // "FORMAT",
                // "SampleName",
            ]
            .iter(),
        )
        .into();

        write!(f, "{}", header_line)
    }
}
