//! Allele observations at the two SNP sites and the APOE (ε1-ε4) caller.
//!
//! Besides the public types this module holds the pure logic of the caller:
//! GT parsing, allele-to-base resolution at a target position, merging of
//! overlapping records and the APOE genotype call.

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
        let mut sorted: Vec<&AlleleObs> = self.alleles.iter().collect();
        sorted.sort_by_key(|a| display_rank(a));
        let mut out = String::new();
        for (i, a) in sorted.iter().enumerate() {
            if i > 0 {
                out.push('/');
            }
            match a {
                AlleleObs::Base(b) => out.push(*b as char),
                AlleleObs::Deleted => out.push_str("del"),
                AlleleObs::Unknown(_) => out.push('?'),
                AlleleObs::Missing => out.push('.'),
            }
        }
        out
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

/// Sort key used by [`SiteCall::display`]: T, C, A, G, other letters
/// (alphabetical), `del`, `?`, `.`.
fn display_rank(a: &AlleleObs) -> (u8, u8) {
    match a {
        AlleleObs::Base(b'T') => (0, 0),
        AlleleObs::Base(b'C') => (1, 0),
        AlleleObs::Base(b'A') => (2, 0),
        AlleleObs::Base(b'G') => (3, 0),
        AlleleObs::Base(b) => (4, *b),
        AlleleObs::Deleted => (5, 0),
        AlleleObs::Unknown(_) => (6, 0),
        AlleleObs::Missing => (7, 0),
    }
}

/// Note added when an unphased double heterozygote cannot be resolved.
pub const NOTE_AMBIGUOUS: &str =
    "unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 (ε1 is very rare)";
/// Note added when both sites are phased but belong to different phase sets.
pub const NOTE_PHASE_SETS_DIFFER: &str = "phase sets differ; treated as unphased";

fn sorted_pair(a: ApoeAllele, b: ApoeAllele) -> [ApoeAllele; 2] {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

/// Returns the two bases of a site when it is exactly two T/C alleles.
fn tc_pair(site: &SiteCall) -> Option<(u8, u8)> {
    match site.alleles.as_slice() {
        [AlleleObs::Base(x), AlleleObs::Base(y)]
            if matches!(x, b'T' | b'C') && matches!(y, b'T' | b'C') =>
        {
            Some((*x, *y))
        }
        _ => None,
    }
}

/// Combines the two site calls into an APOE genotype (see SPEC "APOE calling").
pub fn call_apoe(rs429358: &SiteCall, rs7412: &SiteCall) -> ApoeCall {
    let none = ApoeCall {
        genotype: None,
        phased: false,
        notes: Vec::new(),
    };
    let (Some((a0, a1)), Some((b0, b1))) = (tc_pair(rs429358), tc_pair(rs7412)) else {
        return none;
    };
    let both_phased = rs429358.phased && rs7412.phased;
    let joinable = both_phased && rs429358.phase_set == rs7412.phase_set;
    let hap = ApoeAllele::from_haplotype;

    if joinable {
        return match (hap(a0, b0), hap(a1, b1)) {
            (Some(x), Some(y)) => ApoeCall {
                genotype: Some(sorted_pair(x, y)),
                phased: true,
                notes: Vec::new(),
            },
            _ => none,
        };
    }

    let (Some(p1), Some(p2), Some(q1), Some(q2)) =
        (hap(a0, b0), hap(a1, b1), hap(a0, b1), hap(a1, b0))
    else {
        return none;
    };
    let g1 = sorted_pair(p1, p2);
    let g2 = sorted_pair(q1, q2);
    let mut notes = Vec::new();
    let genotype = if g1 == g2 {
        g1
    } else {
        if both_phased && a0 != a1 && b0 != b1 {
            notes.push(NOTE_PHASE_SETS_DIFFER.to_string());
        }
        notes.push(NOTE_AMBIGUOUS.to_string());
        [ApoeAllele::E2, ApoeAllele::E4]
    };
    ApoeCall {
        genotype: Some(genotype),
        phased: false,
        notes,
    }
}

// ---------------------------------------------------------------------------
// GT parsing
// ---------------------------------------------------------------------------

/// One token of a GT field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GtToken {
    /// `.`
    Missing,
    /// Allele index (0 = REF, i = ALT[i-1]).
    Allele(usize),
    /// Anything that is neither `.` nor an unsigned integer (kept verbatim).
    Invalid(String),
}

/// A parsed GT field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gt {
    /// The GT text exactly as written in the file.
    pub raw: String,
    pub tokens: Vec<GtToken>,
    /// At least one separator and every separator is `|`.
    pub phased: bool,
}

impl Gt {
    /// Number of allele tokens.
    pub fn ploidy(&self) -> usize {
        self.tokens.len()
    }

    /// True when every token is `.`.
    pub fn all_missing(&self) -> bool {
        self.tokens.iter().all(|t| *t == GtToken::Missing)
    }
}

/// Parses a GT string such as `0/1`, `1|0`, `./.`, `1`, `0/0/1`, `10/2`.
pub fn parse_gt(raw: &str) -> Gt {
    let mut tokens = Vec::new();
    let mut separators = 0usize;
    let mut pipes = 0usize;
    let mut start = 0usize;
    for (i, c) in raw.char_indices() {
        if c == '/' || c == '|' {
            separators += 1;
            if c == '|' {
                pipes += 1;
            }
            tokens.push(parse_gt_token(&raw[start..i]));
            start = i + 1;
        }
    }
    tokens.push(parse_gt_token(&raw[start..]));
    Gt {
        raw: raw.to_string(),
        tokens,
        phased: separators > 0 && pipes == separators,
    }
}

fn parse_gt_token(tok: &str) -> GtToken {
    if tok == "." {
        GtToken::Missing
    } else if !tok.is_empty() && tok.bytes().all(|b| b.is_ascii_digit()) {
        match tok.parse::<usize>() {
            Ok(i) => GtToken::Allele(i),
            Err(_) => GtToken::Invalid(tok.to_string()),
        }
    } else {
        GtToken::Invalid(tok.to_string())
    }
}

// ---------------------------------------------------------------------------
// Allele -> base at the target position
// ---------------------------------------------------------------------------

/// Offset of `target` inside a record starting at `pos` (0 when `target < pos`).
fn offset(pos: u64, target: u64) -> usize {
    usize::try_from(target.saturating_sub(pos)).unwrap_or(usize::MAX)
}

/// What the REF allele shows at `target`. Positions past the end of REF
/// (reference blocks, symbolic records) are the build's reference base.
pub fn ref_allele_at(ref_seq: &[u8], pos: u64, target: u64, build_ref: u8) -> AlleleObs {
    match ref_seq.get(offset(pos, target)) {
        Some(b) => AlleleObs::Base(b.to_ascii_uppercase()),
        None => AlleleObs::Base(build_ref),
    }
}

/// What one ALT allele shows at `target` (see SPEC section 4).
pub fn allele_at_target(
    ref_seq: &[u8],
    alt: &[u8],
    pos: u64,
    target: u64,
    build_ref: u8,
) -> AlleleObs {
    let k = offset(pos, target);
    if alt == b"*" {
        return AlleleObs::Deleted;
    }
    if alt.first() == Some(&b'<') {
        if alt.eq_ignore_ascii_case(b"<DEL>") {
            return if target > pos {
                AlleleObs::Deleted
            } else {
                ref_allele_at(ref_seq, pos, pos, build_ref)
            };
        }
        return AlleleObs::Unknown(format!("symbolic allele {}", String::from_utf8_lossy(alt)));
    }
    if alt == b"." {
        return AlleleObs::Unknown("missing ALT allele".to_string());
    }
    if alt.iter().any(|b| *b == b'[' || *b == b']') {
        return AlleleObs::Unknown("breakend allele".to_string());
    }
    if alt.is_empty() || !alt.iter().all(|b| b.is_ascii_alphabetic()) {
        return AlleleObs::Unknown(format!(
            "unsupported allele {}",
            String::from_utf8_lossy(alt)
        ));
    }
    if k >= ref_seq.len() {
        return AlleleObs::Base(build_ref);
    }
    if alt.len() == ref_seq.len() {
        return AlleleObs::Base(alt[k].to_ascii_uppercase());
    }
    // Different length: strip the common prefix and suffix to find the changed core.
    let p = ref_seq
        .iter()
        .zip(alt)
        .take_while(|(r, a)| r.eq_ignore_ascii_case(a))
        .count();
    let s = ref_seq[p..]
        .iter()
        .rev()
        .zip(alt[p..].iter().rev())
        .take_while(|(r, a)| r.eq_ignore_ascii_case(a))
        .count();
    if k < p || k >= ref_seq.len() - s {
        AlleleObs::Base(ref_seq[k].to_ascii_uppercase())
    } else if alt.len() - s == p {
        AlleleObs::Deleted
    } else {
        AlleleObs::Unknown("complex variant overlaps site".to_string())
    }
}

/// Resolves a parsed GT against one record into one [`AlleleObs`] per token.
pub fn resolve_gt(
    gt: &Gt,
    ref_seq: &[u8],
    alts: &[Vec<u8>],
    pos: u64,
    target: u64,
    build_ref: u8,
) -> Vec<AlleleObs> {
    gt.tokens
        .iter()
        .map(|tok| match tok {
            GtToken::Missing => AlleleObs::Missing,
            GtToken::Allele(0) => ref_allele_at(ref_seq, pos, target, build_ref),
            GtToken::Allele(i) => match alts.get(i - 1) {
                Some(alt) => allele_at_target(ref_seq, alt, pos, target, build_ref),
                None => AlleleObs::Unknown(format!("invalid allele index {i}")),
            },
            GtToken::Invalid(t) => AlleleObs::Unknown(format!("invalid allele index {t}")),
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Merging overlapping records
// ---------------------------------------------------------------------------

/// Classification of one record's observation for one sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ObsKind {
    /// No call (all tokens `.`, no GT, zero-depth block).
    Missing,
    /// Ploidy other than 2 (and not all missing).
    NonDiploid,
    /// Diploid with at least one non-missing token.
    Called,
}

/// What one record says about one sample at one target site.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Observation {
    pub kind: ObsKind,
    pub alleles: Vec<AlleleObs>,
    /// GT text as written (empty when there was no GT).
    pub raw_gt: String,
    pub phased: bool,
    pub phase_set: Option<String>,
    /// Explanation for `Missing` observations.
    pub note: Option<String>,
    /// FILTER value of the record when it is neither `PASS` nor `.`.
    pub filter: Option<String>,
}

impl Observation {
    /// A no-call observation with an explanatory note.
    pub fn missing(raw_gt: &str, note: String, filter: Option<String>) -> Observation {
        Observation {
            kind: ObsKind::Missing,
            alleles: vec![AlleleObs::Missing, AlleleObs::Missing],
            raw_gt: raw_gt.to_string(),
            phased: false,
            phase_set: None,
            note: Some(note),
            filter,
        }
    }

    /// Builds an observation from a parsed GT and its resolved alleles.
    pub fn from_gt(
        gt: &Gt,
        alleles: Vec<AlleleObs>,
        phase_set: Option<String>,
        filter: Option<String>,
    ) -> Observation {
        if gt.all_missing() {
            return Observation::missing(&gt.raw, format!("no call (GT {})", gt.raw), filter);
        }
        let kind = if gt.ploidy() == 2 {
            ObsKind::Called
        } else {
            ObsKind::NonDiploid
        };
        Observation {
            kind,
            alleles,
            raw_gt: gt.raw.clone(),
            phased: gt.phased,
            phase_set,
            note: None,
            filter,
        }
    }
}

/// A call with no data at all.
pub fn no_call(note: &str) -> SiteCall {
    SiteCall {
        alleles: vec![AlleleObs::Missing, AlleleObs::Missing],
        phased: false,
        phase_set: None,
        notes: vec![note.to_string()],
    }
}

/// Homozygous reference call used when the site is absent from the VCF.
pub fn assumed_reference(ref_base: u8) -> SiteCall {
    let r = ref_base as char;
    SiteCall {
        alleles: vec![AlleleObs::Base(ref_base), AlleleObs::Base(ref_base)],
        phased: false,
        phase_set: None,
        notes: vec![format!("not in VCF; assumed reference ({r}/{r})")],
    }
}

/// Combines all observations overlapping one site for one sample (SPEC section 5).
pub fn combine_observations(obs: &[Observation], ref_base: u8) -> SiteCall {
    let called: Vec<&Observation> = obs.iter().filter(|o| o.kind == ObsKind::Called).collect();
    if called.is_empty() {
        if let Some(nd) = obs.iter().find(|o| o.kind == ObsKind::NonDiploid) {
            return SiteCall {
                alleles: nd.alleles.clone(),
                phased: false,
                phase_set: None,
                notes: vec![format!("non-diploid genotype ({}); no call", nd.raw_gt)],
            };
        }
        let note = obs
            .iter()
            .find_map(|o| o.note.clone())
            .unwrap_or_else(|| "no call".to_string());
        return no_call(&note);
    }

    let ref_obs = AlleleObs::Base(ref_base);
    let informative: Vec<&Observation> = called
        .iter()
        .copied()
        .filter(|o| o.alleles.iter().any(|a| *a != ref_obs))
        .collect();
    let base = informative.first().copied().unwrap_or(called[0]);
    let mut contributors: Vec<&Observation> = vec![base];
    let mut alleles = base.alleles.clone();
    let mut merged = false;

    for other in informative.iter().skip(1) {
        contributors.push(other);
        let mut seen: Vec<&AlleleObs> = Vec::new();
        for v in other.alleles.iter() {
            if *v == ref_obs || *v == AlleleObs::Missing || seen.contains(&v) {
                continue;
            }
            seen.push(v);
            let in_other = other.alleles.iter().filter(|a| *a == v).count();
            let in_base = alleles.iter().filter(|a| *a == v).count();
            for _ in 0..in_other.saturating_sub(in_base) {
                match alleles
                    .iter_mut()
                    .find(|a| **a == ref_obs || **a == AlleleObs::Missing)
                {
                    Some(slot) => {
                        *slot = v.clone();
                        merged = true;
                    }
                    None => {
                        let unknown = AlleleObs::Unknown("conflicting overlapping records".into());
                        return SiteCall {
                            alleles: vec![unknown.clone(), unknown],
                            phased: false,
                            phase_set: None,
                            notes: vec!["conflicting overlapping records; no call".to_string()],
                        };
                    }
                }
            }
        }
    }

    let mut notes = Vec::new();
    if alleles.contains(&AlleleObs::Missing) {
        notes.push(format!("no call (GT {})", base.raw_gt));
    }
    if alleles.contains(&AlleleObs::Deleted) {
        notes.push("allele deleted by overlapping deletion".to_string());
    }
    let mut seen_bases: Vec<u8> = Vec::new();
    let mut seen_reasons: Vec<&str> = Vec::new();
    for a in &alleles {
        match a {
            AlleleObs::Base(b) if !matches!(b, b'T' | b'C') && !seen_bases.contains(b) => {
                seen_bases.push(*b);
                notes.push(format!("unexpected allele {}", *b as char));
            }
            AlleleObs::Unknown(r) if !seen_reasons.contains(&r.as_str()) => {
                seen_reasons.push(r);
                notes.push(format!("unresolvable allele ({r})"));
            }
            _ => {}
        }
    }
    let mut seen_filters: Vec<&str> = Vec::new();
    for c in &contributors {
        if let Some(f) = &c.filter {
            if !seen_filters.contains(&f.as_str()) {
                seen_filters.push(f);
                notes.push(format!("FILTER={f}"));
            }
        }
    }
    let (phased, phase_set) = if merged {
        (false, None)
    } else {
        (base.phased, base.phase_set.clone())
    };
    SiteCall {
        alleles,
        phased,
        phase_set,
        notes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(x: u8) -> AlleleObs {
        AlleleObs::Base(x)
    }

    fn site(alleles: Vec<AlleleObs>, phased: bool, ps: Option<&str>) -> SiteCall {
        SiteCall {
            alleles,
            phased,
            phase_set: ps.map(str::to_string),
            notes: vec![],
        }
    }

    fn pair(x: u8, y: u8, phased: bool, ps: Option<&str>) -> SiteCall {
        site(vec![b(x), b(y)], phased, ps)
    }

    // ---- GT parsing -------------------------------------------------------

    #[test]
    fn vcf_gt_parsing_unphased_phased_haploid_polyploid() {
        let g = parse_gt("0/1");
        assert_eq!(g.tokens, vec![GtToken::Allele(0), GtToken::Allele(1)]);
        assert!(!g.phased);
        let g = parse_gt("1|0");
        assert!(g.phased);
        assert_eq!(g.ploidy(), 2);
        let g = parse_gt("1");
        assert!(!g.phased);
        assert_eq!(g.ploidy(), 1);
        let g = parse_gt("0/0/1");
        assert_eq!(g.ploidy(), 3);
        // Mixed separators are not phased.
        assert!(!parse_gt("0|1/1").phased);
        assert!(parse_gt("0|1|1").phased);
    }

    #[test]
    fn vcf_gt_parsing_missing_partial_and_multi_digit() {
        let g = parse_gt("./.");
        assert!(g.all_missing());
        let g = parse_gt(".|.");
        assert!(g.all_missing() && g.phased);
        assert!(parse_gt(".").all_missing());
        let g = parse_gt("0/.");
        assert!(!g.all_missing());
        assert_eq!(g.tokens, vec![GtToken::Allele(0), GtToken::Missing]);
        let g = parse_gt("10/12");
        assert_eq!(g.tokens, vec![GtToken::Allele(10), GtToken::Allele(12)]);
        let g = parse_gt("a/-1");
        assert_eq!(
            g.tokens,
            vec![GtToken::Invalid("a".into()), GtToken::Invalid("-1".into())]
        );
        assert_eq!(parse_gt("0/1").raw, "0/1");
    }

    // ---- allele -> base at target ----------------------------------------

    #[test]
    fn snp_allele_at_target_substitution() {
        assert_eq!(allele_at_target(b"T", b"C", 100, 100, b'T'), b(b'C'));
        assert_eq!(allele_at_target(b"T", b"g", 100, 100, b'T'), b(b'G'));
        assert_eq!(ref_allele_at(b"t", 100, 100, b'T'), b(b'T'));
    }

    #[test]
    fn snp_mnp_allele_at_target() {
        // GT>GC at 683: target 684 is the second base.
        assert_eq!(allele_at_target(b"GT", b"GC", 683, 684, b'T'), b(b'C'));
        assert_eq!(allele_at_target(b"GT", b"GC", 683, 683, b'T'), b(b'G'));
        assert_eq!(ref_allele_at(b"GT", 683, 684, b'T'), b(b'T'));
    }

    #[test]
    fn indel_insertion_leaves_target_unchanged() {
        assert_eq!(allele_at_target(b"T", b"TGC", 100, 100, b'T'), b(b'T'));
        // Insertion before the final REF base: the target base is still REF[k].
        assert_eq!(allele_at_target(b"CG", b"CTG", 10, 11, b'G'), b(b'G'));
    }

    #[test]
    fn indel_deletion_upstream_and_after_the_site() {
        // CGTGT>C at 680 deletes 681..684: target 684 is deleted.
        assert_eq!(
            allele_at_target(b"CGTGT", b"C", 680, 684, b'T'),
            AlleleObs::Deleted
        );
        // Target is the padding base itself: unchanged.
        assert_eq!(allele_at_target(b"CGTGT", b"C", 680, 680, b'C'), b(b'C'));
        // Deletion that starts at the site (CG>C at target) keeps the first base.
        assert_eq!(allele_at_target(b"CG", b"C", 822, 822, b'C'), b(b'C'));
        // Deletion after the site: target is the padding base.
        assert_eq!(allele_at_target(b"CG", b"C", 822, 822, b'C'), b(b'C'));
    }

    #[test]
    fn indel_complex_variant_is_unknown() {
        assert_eq!(
            allele_at_target(b"GTG", b"A", 683, 684, b'T'),
            AlleleObs::Unknown("complex variant overlaps site".into())
        );
    }

    #[test]
    fn indel_star_del_and_symbolic_alleles() {
        assert_eq!(
            allele_at_target(b"T", b"*", 100, 100, b'T'),
            AlleleObs::Deleted
        );
        assert_eq!(
            allele_at_target(b"C", b"<DEL>", 100, 150, b'T'),
            AlleleObs::Deleted
        );
        assert_eq!(allele_at_target(b"C", b"<DEL>", 100, 100, b'C'), b(b'C'));
        for alt in [&b"<NON_REF>"[..], b"<*>", b"<INS>"] {
            assert!(matches!(
                allele_at_target(b"C", alt, 100, 150, b'T'),
                AlleleObs::Unknown(_)
            ));
        }
        assert!(matches!(
            allele_at_target(b"C", b"C[chr2:5[", 100, 100, b'C'),
            AlleleObs::Unknown(_)
        ));
        assert!(matches!(
            allele_at_target(b"C", b".", 100, 100, b'C'),
            AlleleObs::Unknown(_)
        ));
    }

    #[test]
    fn gvcf_reference_block_alleles_use_build_base() {
        // The block REF is a single padding base far upstream of the target.
        assert_eq!(ref_allele_at(b"C", 100, 150, b'T'), b(b'T'));
        let gt = parse_gt("0/0");
        let alleles = resolve_gt(&gt, b"C", &[b"<NON_REF>".to_vec()], 100, 150, b'T');
        assert_eq!(alleles, vec![b(b'T'), b(b'T')]);
        let gt = parse_gt("0/1");
        let alleles = resolve_gt(&gt, b"C", &[b"<NON_REF>".to_vec()], 100, 150, b'T');
        assert_eq!(alleles[0], b(b'T'));
        assert!(matches!(alleles[1], AlleleObs::Unknown(_)));
    }

    #[test]
    fn vcf_resolve_gt_invalid_index_and_missing() {
        let gt = parse_gt("./2");
        let alleles = resolve_gt(&gt, b"T", &[b"C".to_vec()], 5, 5, b'T');
        assert_eq!(alleles[0], AlleleObs::Missing);
        assert_eq!(
            alleles[1],
            AlleleObs::Unknown("invalid allele index 2".into())
        );
    }

    // ---- merging ----------------------------------------------------------

    fn obs(gt: &str, alleles: Vec<AlleleObs>) -> Observation {
        Observation::from_gt(&parse_gt(gt), alleles, None, None)
    }

    #[test]
    fn snp_merge_split_multiallelic_records() {
        let a = obs("1/0", vec![b(b'C'), b(b'T')]);
        let c = obs("0/1", vec![b(b'T'), b(b'G')]);
        let call = combine_observations(&[a, c], b'T');
        assert_eq!(call.display(), "C/G");
        assert!(!call.phased);
        assert_eq!(call.notes, vec!["unexpected allele G"]);
    }

    #[test]
    fn snp_merge_duplicate_records_are_idempotent() {
        let a = obs("0|1", vec![b(b'T'), b(b'C')]);
        let call = combine_observations(&[a.clone(), a], b'T');
        assert_eq!(call.display(), "T/C");
        assert!(call.phased, "no replacement happened, phasing is kept");
    }

    #[test]
    fn snp_merge_conflict_gives_unknown_pair() {
        let a = obs("1/1", vec![b(b'C'), b(b'C')]);
        let c = obs("0/1", vec![b(b'T'), b(b'G')]);
        let call = combine_observations(&[a, c], b'T');
        assert_eq!(call.display(), "?/?");
        assert_eq!(call.notes, vec!["conflicting overlapping records; no call"]);
    }

    #[test]
    fn indel_merge_deletion_and_snp() {
        let del = obs("0/1", vec![b(b'T'), AlleleObs::Deleted]);
        let snp = obs("0/1", vec![b(b'T'), b(b'C')]);
        let call = combine_observations(&[del, snp], b'T');
        assert_eq!(call.display(), "C/del");
        assert_eq!(call.notes, vec!["allele deleted by overlapping deletion"]);
    }

    #[test]
    fn indel_merge_ref_only_records_keep_first_called() {
        let a = obs("0|0", vec![b(b'T'), b(b'T')]);
        let c = obs("0/0", vec![b(b'T'), b(b'T')]);
        let call = combine_observations(&[a, c], b'T');
        assert_eq!(call.display(), "T/T");
        assert!(call.phased);
    }

    #[test]
    fn vcf_merge_no_called_observations() {
        let nd = obs("0/0/1", vec![b(b'T'), b(b'T'), b(b'C')]);
        let call = combine_observations(&[nd], b'T');
        assert_eq!(call.display(), "T/T/C");
        assert_eq!(call.notes, vec!["non-diploid genotype (0/0/1); no call"]);

        let m = obs("./.", vec![AlleleObs::Missing, AlleleObs::Missing]);
        let call = combine_observations(&[m], b'T');
        assert_eq!(call.display(), "./.");
        assert_eq!(call.notes, vec!["no call (GT ./.)"]);

        let none = Observation::missing("", "no GT field; no call".into(), None);
        assert_eq!(
            combine_observations(&[none], b'T').notes,
            vec!["no GT field; no call"]
        );
    }

    #[test]
    fn vcf_merge_partial_call_and_unknown_notes() {
        let p = obs("0/.", vec![b(b'T'), AlleleObs::Missing]);
        let call = combine_observations(&[p], b'T');
        assert_eq!(call.display(), "T/.");
        assert_eq!(call.notes, vec!["no call (GT 0/.)"]);

        let u = obs(
            "0/1",
            vec![
                b(b'T'),
                AlleleObs::Unknown("complex variant overlaps site".into()),
            ],
        );
        let call = combine_observations(&[u], b'T');
        assert_eq!(call.display(), "T/?");
        assert_eq!(
            call.notes,
            vec!["unresolvable allele (complex variant overlaps site)"]
        );
    }

    #[test]
    fn vcf_merge_filter_notes_once_per_value() {
        let mut a = obs("0/1", vec![b(b'T'), b(b'C')]);
        a.filter = Some("LowQual".into());
        let mut c = obs("0/1", vec![b(b'T'), b(b'G')]);
        c.filter = Some("LowQual".into());
        let call = combine_observations(&[a, c], b'T');
        assert_eq!(
            call.notes
                .iter()
                .filter(|n| n.starts_with("FILTER="))
                .count(),
            1
        );
    }

    // ---- absent-site rules ------------------------------------------------

    #[test]
    fn wildtype_absent_site_calls() {
        let t = assumed_reference(b'T');
        assert_eq!(t.display(), "T/T");
        assert_eq!(t.notes, vec!["not in VCF; assumed reference (T/T)"]);
        let c = assumed_reference(b'C');
        assert_eq!(c.notes, vec!["not in VCF; assumed reference (C/C)"]);
        let n = no_call("not covered by gVCF; no call");
        assert_eq!(n.display(), "./.");
    }

    // ---- APOE calling -----------------------------------------------------

    fn unphased(a: (u8, u8), b: (u8, u8)) -> ApoeCall {
        call_apoe(&pair(a.0, a.1, false, None), &pair(b.0, b.1, false, None))
    }

    fn g(x: ApoeAllele, y: ApoeAllele) -> Option<[ApoeAllele; 2]> {
        Some([x, y])
    }

    #[test]
    fn unphased_all_nine_combinations() {
        use ApoeAllele::*;
        let (t, c) = (b'T', b'C');
        let cases = [
            ((t, t), (t, t), g(E2, E2)),
            ((t, t), (t, c), g(E2, E3)),
            ((t, c), (t, c), g(E2, E4)),
            ((t, t), (c, c), g(E3, E3)),
            ((t, c), (c, c), g(E3, E4)),
            ((c, c), (c, c), g(E4, E4)),
            ((c, c), (t, t), g(E1, E1)),
            ((t, c), (t, t), g(E1, E2)),
            ((c, c), (t, c), g(E1, E4)),
        ];
        for (a, bb, want) in cases {
            let call = unphased(a, bb);
            assert_eq!(call.genotype, want, "{a:?} {bb:?}");
            assert!(!call.phased);
        }
    }

    #[test]
    fn unphased_double_heterozygote_is_ambiguous() {
        let call = unphased((b'T', b'C'), (b'T', b'C'));
        assert_eq!(call.genotype, g(ApoeAllele::E2, ApoeAllele::E4));
        assert_eq!(call.notes, vec![NOTE_AMBIGUOUS]);
        assert!(NOTE_AMBIGUOUS.contains("ε2/ε4 or ε1/ε3"));
    }

    #[test]
    fn phased_haplotypes_resolve_double_heterozygote() {
        use ApoeAllele::*;
        // 0|1 / 0|1: haplotypes (T,C)=e3 and (C,T)=e1.
        let call = call_apoe(&pair(b'T', b'C', true, None), &pair(b'C', b'T', true, None));
        assert_eq!(call.genotype, g(E1, E3));
        assert!(call.phased && call.notes.is_empty());
        // (T,T)=e2 and (C,C)=e4.
        let call = call_apoe(&pair(b'T', b'C', true, None), &pair(b'T', b'C', true, None));
        assert_eq!(call.genotype, g(E2, E4));
        assert!(call.phased);
    }

    #[test]
    fn phased_same_phase_set_is_joinable() {
        let call = call_apoe(
            &pair(b'C', b'T', true, Some("100")),
            &pair(b'T', b'C', true, Some("100")),
        );
        assert_eq!(call.genotype, g(ApoeAllele::E1, ApoeAllele::E3));
        assert!(call.phased);
    }

    #[test]
    fn phased_different_phase_sets_are_unphased_with_notes() {
        let call = call_apoe(
            &pair(b'C', b'T', true, Some("100")),
            &pair(b'C', b'T', true, Some("200")),
        );
        assert_eq!(call.genotype, g(ApoeAllele::E2, ApoeAllele::E4));
        assert!(!call.phased);
        assert_eq!(call.notes, vec![NOTE_PHASE_SETS_DIFFER, NOTE_AMBIGUOUS]);
        // One side without PS also differs.
        let call = call_apoe(
            &pair(b'C', b'T', true, Some("100")),
            &pair(b'C', b'T', true, None),
        );
        assert_eq!(call.notes.len(), 2);
        // A homozygous site makes the phase set irrelevant: no notes.
        let call = call_apoe(
            &pair(b'C', b'C', true, Some("100")),
            &pair(b'C', b'T', true, Some("200")),
        );
        assert_eq!(call.genotype, g(ApoeAllele::E1, ApoeAllele::E4));
        assert!(!call.phased && call.notes.is_empty());
    }

    #[test]
    fn phased_mixed_with_unphased_is_not_joined() {
        let call = call_apoe(
            &pair(b'T', b'C', true, None),
            &pair(b'T', b'C', false, None),
        );
        assert_eq!(call.genotype, g(ApoeAllele::E2, ApoeAllele::E4));
        assert!(!call.phased);
        assert_eq!(call.notes, vec![NOTE_AMBIGUOUS]);
    }

    #[test]
    fn snp_call_apoe_not_determinable() {
        let good = pair(b'T', b'T', false, None);
        let cases = vec![
            site(vec![AlleleObs::Missing, AlleleObs::Missing], false, None),
            site(vec![b(b'T'), AlleleObs::Missing], false, None),
            site(vec![b(b'T'), AlleleObs::Deleted], false, None),
            site(vec![b(b'T'), AlleleObs::Unknown("x".into())], false, None),
            site(vec![b(b'C'), b(b'G')], false, None),
            site(vec![b(b'C')], false, None),
            site(vec![b(b'T'), b(b'T'), b(b'C')], false, None),
        ];
        for bad in cases {
            for call in [call_apoe(&bad, &good), call_apoe(&good, &bad)] {
                assert_eq!(call.genotype, None);
                assert!(!call.phased);
            }
        }
    }

    #[test]
    fn apoe_allele_helpers() {
        assert_eq!(ApoeAllele::from_haplotype(b'c', b't'), Some(ApoeAllele::E1));
        assert_eq!(ApoeAllele::from_haplotype(b'A', b'T'), None);
        assert_eq!(ApoeAllele::E4.to_string(), "ε4");
        assert_eq!(ApoeAllele::E2.number(), 2);
    }

    // ---- display ordering -------------------------------------------------

    #[test]
    fn snp_display_ordering() {
        let d = |a: Vec<AlleleObs>| site(a, false, None).display();
        assert_eq!(d(vec![b(b'C'), b(b'T')]), "T/C");
        assert_eq!(d(vec![b(b'T'), b(b'C')]), "T/C");
        assert_eq!(d(vec![b(b'G'), b(b'A')]), "A/G");
        assert_eq!(d(vec![b(b'G'), b(b'C')]), "C/G");
        assert_eq!(d(vec![b(b'N'), b(b'G')]), "G/N");
        assert_eq!(d(vec![AlleleObs::Deleted, b(b'C')]), "C/del");
        assert_eq!(
            d(vec![AlleleObs::Unknown("r".into()), AlleleObs::Deleted]),
            "del/?"
        );
        assert_eq!(d(vec![AlleleObs::Missing, b(b'T')]), "T/.");
        assert_eq!(d(vec![AlleleObs::Missing, AlleleObs::Missing]), "./.");
        assert_eq!(d(vec![b(b'C')]), "C");
        assert_eq!(d(vec![b(b'C'), b(b'T'), b(b'T')]), "T/T/C");
    }
}
