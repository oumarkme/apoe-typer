//! Allele observations at the two SNP sites and the APOE (ε1-ε4) caller.
//!
//! SKELETON: types and signatures are the public contract (see SPEC). Bodies
//! marked `todo!()` are to be implemented.

use std::fmt;

/// What one haplotype carries at a target SNP position.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AlleleObs {
    /// A concrete base at the target position (always stored uppercase, e.g. `b'T'`).
    Base(u8),
    /// The target base is deleted on this haplotype (overlapping deletion,
    /// spanning-deletion `*` allele or symbolic `<DEL>`).
    Deleted,
    /// Called, but not resolvable to a base at the target (symbolic allele,
    /// complex indel, invalid allele index, ...). The text is a short reason.
    Unknown(String),
    /// No call (`.` in GT, GT absent, zero-depth reference block, no data).
    Missing,
}

/// One APOE allele.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ApoeAllele {
    E1,
    E2,
    E3,
    E4,
}

impl ApoeAllele {
    /// 1, 2, 3 or 4.
    pub fn number(self) -> u8 {
        match self {
            ApoeAllele::E1 => 1,
            ApoeAllele::E2 => 2,
            ApoeAllele::E3 => 3,
            ApoeAllele::E4 => 4,
        }
    }

    /// Haplotype (rs429358 base, rs7412 base) -> APOE allele:
    /// (C,T)=ε1, (T,T)=ε2, (T,C)=ε3, (C,C)=ε4; anything else -> None.
    pub fn from_haplotype(rs429358: u8, rs7412: u8) -> Option<ApoeAllele> {
        match (rs429358.to_ascii_uppercase(), rs7412.to_ascii_uppercase()) {
            (b'C', b'T') => Some(ApoeAllele::E1),
            (b'T', b'T') => Some(ApoeAllele::E2),
            (b'T', b'C') => Some(ApoeAllele::E3),
            (b'C', b'C') => Some(ApoeAllele::E4),
            _ => None,
        }
    }
}

impl fmt::Display for ApoeAllele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ε{}", self.number())
    }
}

/// The resolved genotype of one sample at one SNP site.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SiteCall {
    /// Observed alleles, normally exactly 2 (diploid). For "no data" use two `Missing`.
    pub alleles: Vec<AlleleObs>,
    /// True if the GT used `|` separators and the call was not merged from several records.
    pub phased: bool,
    /// FORMAT/PS value when present and not `.`.
    pub phase_set: Option<String>,
    /// Site-level notes WITHOUT the `rsid: ` prefix (the analysis layer adds it).
    pub notes: Vec<String>,
}

impl SiteCall {
    /// Display string, e.g. `T/C`, `C/C`, `./.`, `T/del`, `C/G`, `T/?`.
    /// Alleles are always joined with `/` and ordered T, C, A, G, other letters
    /// (alphabetical), `del`, `?`, `.` - so heterozygotes read `T/C` exactly as in
    /// the lab template, regardless of REF/ALT order or GT order.
    pub fn display(&self) -> String {
        todo!("SiteCall::display")
    }
}

/// Result of combining the rs429358 and rs7412 calls of one sample.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ApoeCall {
    /// Sorted ascending (`[E2, E4]` for ε2/ε4). `None` = not determinable (NA).
    pub genotype: Option<[ApoeAllele; 2]>,
    /// True when both sites are phased (`|`) and in the same phase set.
    pub phased: bool,
    /// Sample-level notes (ambiguity, phase-set mismatch, ...).
    pub notes: Vec<String>,
}

/// Combines the two site calls into an APOE genotype (see SPEC "APOE calling").
pub fn call_apoe(rs429358: &SiteCall, rs7412: &SiteCall) -> ApoeCall {
    let _ = (rs429358, rs7412);
    todo!("call_apoe")
}
