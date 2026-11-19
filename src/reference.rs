//! Reference genome builds and the two APOE-defining SNP sites.
//!
//! Coordinates were verified against the UCSC reference sequences (hg19, hg38,
//! hs1): the 340 bp windows around both SNPs are identical in all three builds,
//! and the hs1 positions agree with the UCSC hg38ToHs1 liftOver chain.

use std::fmt;
use std::str::FromStr;

/// Supported reference genome builds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Build {
    /// GRCh37 / hg19.
    Hg19,
    /// GRCh38 / hg38.
    Hg38,
    /// T2T-CHM13v2.0 / hs1.
    Hs1,
}

/// One APOE-defining SNP in a particular build (1-based position, + strand bases).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SnpSite {
    /// dbSNP identifier, e.g. `rs429358`.
    pub rsid: &'static str,
    /// 1-based position on chromosome 19.
    pub pos: u64,
    /// Reference base (+ strand) in this build.
    pub ref_base: u8,
    /// The non-reference base that defines the APOE allele
    /// (`C` at rs429358 marks ε4/ε1, `T` at rs7412 marks ε2/ε1).
    pub alt_base: u8,
}

/// Canonical chromosome label used in report headers, independent of the VCF naming.
pub const CHR19_LABEL: &str = "chr19";

impl Build {
    /// All supported builds, in CLI/GUI display order.
    pub const ALL: [Build; 3] = [Build::Hg19, Build::Hg38, Build::Hs1];

    /// Short identifier used on the command line (`hg19`, `hg38`, `hs1`).
    pub fn id(self) -> &'static str {
        match self {
            Build::Hg19 => "hg19",
            Build::Hg38 => "hg38",
            Build::Hs1 => "hs1",
        }
    }

    /// Human readable description.
    pub fn description(self) -> &'static str {
        match self {
            Build::Hg19 => "GRCh37/hg19",
            Build::Hg38 => "GRCh38/hg38",
            Build::Hs1 => "T2T-CHM13v2.0/hs1",
        }
    }

    /// Length of chromosome 19 in this build (UCSC chromInfo).
    pub fn chr19_length(self) -> u64 {
        match self {
            Build::Hg19 => 59_128_983,
            Build::Hg38 => 58_617_616,
            Build::Hs1 => 61_707_364,
        }
    }

    /// RefSeq and GenBank accessions of chromosome 19 in this build (UCSC chromAlias).
    pub fn chr19_accessions(self) -> &'static [&'static str] {
        match self {
            Build::Hg19 => &["NC_000019.9", "CM000681.1"],
            Build::Hg38 => &["NC_000019.10", "CM000681.2"],
            Build::Hs1 => &["NC_060943.1", "CP068259.2"],
        }
    }

    /// rs429358 (T>C; C marks ε4 and ε1).
    pub fn rs429358(self) -> SnpSite {
        SnpSite {
            rsid: "rs429358",
            pos: match self {
                Build::Hg19 => 45_411_941,
                Build::Hg38 => 44_908_684,
                Build::Hs1 => 47_733_380,
            },
            ref_base: b'T',
            alt_base: b'C',
        }
    }

    /// rs7412 (C>T; T marks ε2 and ε1).
    pub fn rs7412(self) -> SnpSite {
        SnpSite {
            rsid: "rs7412",
            pos: match self {
                Build::Hg19 => 45_412_079,
                Build::Hg38 => 44_908_822,
                Build::Hs1 => 47_733_518,
            },
            ref_base: b'C',
            alt_base: b'T',
        }
    }

    /// Both sites in report order: `[rs429358, rs7412]`.
    pub fn sites(self) -> [SnpSite; 2] {
        [self.rs429358(), self.rs7412()]
    }

    /// Parses a build name, case-insensitively, accepting common aliases
    /// (`hg19`/`GRCh37`/`b37`, `hg38`/`GRCh38`, `hs1`/`CHM13`/`T2T`).
    pub fn from_name(name: &str) -> Option<Build> {
        match name.trim().to_ascii_lowercase().as_str() {
            "hg19" | "grch37" | "b37" | "hs37d5" => Some(Build::Hg19),
            "hg38" | "grch38" | "b38" => Some(Build::Hg38),
            "hs1" | "chm13" | "chm13v2" | "chm13v2.0" | "t2t" | "t2t-chm13" | "t2t-chm13v2"
            | "t2t-chm13v2.0" => Some(Build::Hs1),
            _ => None,
        }
    }

    /// Identifies the build from a declared chromosome 19 length.
    pub fn from_chr19_length(len: u64) -> Option<Build> {
        Build::ALL.into_iter().find(|b| b.chr19_length() == len)
    }

    /// Identifies the build from a versioned chromosome 19 accession
    /// (e.g. `NC_000019.10` -> hg38). Case-insensitive.
    pub fn from_chr19_accession(name: &str) -> Option<Build> {
        Build::ALL.into_iter().find(|b| {
            b.chr19_accessions()
                .iter()
                .any(|acc| acc.eq_ignore_ascii_case(name.trim()))
        })
    }
}

impl fmt::Display for Build {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.id())
    }
}

impl FromStr for Build {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Build::from_name(s)
            .ok_or_else(|| format!("unknown reference build '{s}' (expected hg19, hg38 or hs1)"))
    }
}

/// Returns true if `name` denotes chromosome 19 in any supported naming scheme:
/// `19`, `chr19` (any letter case) or a known RefSeq/GenBank accession of any build.
pub fn is_chr19_name(name: &str) -> bool {
    let name = name.trim();
    name == "19" || name.eq_ignore_ascii_case("chr19") || Build::from_chr19_accession(name).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hg19_sites() {
        assert_eq!(Build::Hg19.rs429358().pos, 45_411_941);
        assert_eq!(Build::Hg19.rs7412().pos, 45_412_079);
    }

    #[test]
    fn hg38_sites() {
        assert_eq!(Build::Hg38.rs429358().pos, 44_908_684);
        assert_eq!(Build::Hg38.rs7412().pos, 44_908_822);
    }

    #[test]
    fn hs1_sites() {
        assert_eq!(Build::Hs1.rs429358().pos, 47_733_380);
        assert_eq!(Build::Hs1.rs7412().pos, 47_733_518);
    }

    #[test]
    fn snp_spacing_is_identical_in_all_builds() {
        for b in Build::ALL {
            assert_eq!(b.rs7412().pos - b.rs429358().pos, 138, "{b}");
            assert_eq!(b.rs429358().ref_base, b'T');
            assert_eq!(b.rs429358().alt_base, b'C');
            assert_eq!(b.rs7412().ref_base, b'C');
            assert_eq!(b.rs7412().alt_base, b'T');
        }
    }

    #[test]
    fn parse_build_names() {
        assert_eq!("hg19".parse::<Build>(), Ok(Build::Hg19));
        assert_eq!("GRCh37".parse::<Build>(), Ok(Build::Hg19));
        assert_eq!("HG38".parse::<Build>(), Ok(Build::Hg38));
        assert_eq!("hs1".parse::<Build>(), Ok(Build::Hs1));
        assert_eq!("T2T-CHM13v2.0".parse::<Build>(), Ok(Build::Hs1));
        assert!("hg12".parse::<Build>().is_err());
    }

    #[test]
    fn build_from_chr19_length_and_accession() {
        assert_eq!(Build::from_chr19_length(59_128_983), Some(Build::Hg19));
        assert_eq!(Build::from_chr19_length(58_617_616), Some(Build::Hg38));
        assert_eq!(Build::from_chr19_length(61_707_364), Some(Build::Hs1));
        assert_eq!(Build::from_chr19_length(1), None);
        assert_eq!(Build::from_chr19_accession("NC_000019.9"), Some(Build::Hg19));
        assert_eq!(Build::from_chr19_accession("NC_000019.10"), Some(Build::Hg38));
        assert_eq!(Build::from_chr19_accession("NC_060943.1"), Some(Build::Hs1));
        assert_eq!(Build::from_chr19_accession("chr19"), None);
    }

    #[test]
    fn chr19_names() {
        for n in ["19", "chr19", "CHR19", "Chr19", "NC_000019.10", "CP068259.2"] {
            assert!(is_chr19_name(n), "{n}");
        }
        for n in ["1", "chr1", "chr19_KI270938v1_alt", "chr190", "NC_000019"] {
            assert!(!is_chr19_name(n), "{n}");
        }
    }
}
