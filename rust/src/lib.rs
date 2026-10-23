//! Batched Rust core for PairLoom's per-base hot paths.
//!
//! This crate is a bit-exact port of the Python reference implementations in
//! `pairloom.alignment.build_read_structure`, `pairloom.merge.try_merge`,
//! and `pairloom.overlap.aggregate_overlap`. The Python implementations
//! remain the semantic reference and the unit-test target; the pipeline
//! calls this core for speed. Any divergence is a bug — the benchmark
//! equivalence gate compares outputs byte for byte.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::collections::{BTreeMap, BTreeSet};

const CIGAR_M: u8 = 0;
const CIGAR_I: u8 = 1;
const CIGAR_D: u8 = 2;
const CIGAR_N: u8 = 3;
const CIGAR_S: u8 = 4;
const CIGAR_H: u8 = 5;
const CIGAR_EQUAL: u8 = 7;
const CIGAR_X: u8 = 8;

const UNRESOLVED_QUALITY: u8 = 2;

fn is_query_op(op: u8) -> bool {
    matches!(op, CIGAR_M | CIGAR_I | CIGAR_S | CIGAR_EQUAL | CIGAR_X)
}
fn is_aligned_op(op: u8) -> bool {
    matches!(op, CIGAR_M | CIGAR_EQUAL | CIGAR_X)
}
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Port of `pairloom.cigarops.normalize_cigar`: merge adjacent identical ops.
fn normalize_cigar(cigar: &[(u8, u32)]) -> Vec<(u8, u32)> {
    let mut merged: Vec<(u8, u32)> = Vec::with_capacity(cigar.len());
    for &(op, length) in cigar {
        if let Some(last) = merged.last_mut() {
            if last.0 == op {
                last.1 += length;
                continue;
            }
        }
        merged.push((op, length));
    }
    merged
}

fn validate_cigar(cigar: &[(u8, u32)], seq_len: usize, context: &str) -> PyResult<()> {
    if cigar.is_empty() {
        return Err(data_invariant(format!("empty CIGAR for {context}")));
    }
    let mut total_q: u32 = 0;
    for &(op, length) in cigar {
        if op > 8 || op == 6 {
            return Err(data_invariant(format!(
                "unsupported CIGAR op {op} (len {length}) for {context}"
            )));
        }
        if length == 0 {
            return Err(data_invariant(format!(
                "non-positive CIGAR length {length} (op {op}) for {context}"
            )));
        }
        if is_query_op(op) {
            total_q += length;
        }
    }
    if total_q as usize != seq_len {
        return Err(data_invariant(format!(
            "CIGAR consumes {total_q} query bases but sequence has {seq_len} for {context}"
        )));
    }
    Ok(())
}

fn data_invariant(msg: String) -> PyErr {
    PyErr::new::<pyo3::exceptions::PyException, _>(format!("DataInvariantError: {msg}"))
}

/// A read expressed in reference orientation — port of `ReadStructure`.
#[pyclass]
struct ReadStructureR {
    #[pyo3(get)]
    tid: i32,
    #[pyo3(get)]
    start0: i32,
    #[pyo3(get)]
    end0: i32,
    #[pyo3(get)]
    is_reverse: bool,
    mapq: i32,
    #[pyo3(get)]
    query_name: String,
    /// parallel arrays: ascending reference position, base, phred quality
    positions: Vec<i32>,
    bases: Vec<u8>,
    quals: Vec<u8>,
    deleted: BTreeSet<i32>,
    skipped: BTreeSet<i32>,
    /// boundary -> (inserted bases, phred qualities)
    insertions: BTreeMap<i32, (Vec<u8>, Vec<u8>)>,
    #[pyo3(get)]
    lead_clip_bases: Vec<u8>,
    lead_clip_quals: Vec<u8>,
    #[pyo3(get)]
    trail_clip_bases: Vec<u8>,
    trail_clip_quals: Vec<u8>,
    #[pyo3(get)]
    hard_clipped: bool,
    #[pyo3(get)]
    has_qualities: bool,
    #[pyo3(get)]
    read_group: Option<String>,
    #[pyo3(get)]
    aligned_len: i32,
}

#[pymethods]
impl ReadStructureR {
    #[getter]
    fn mapq(&self) -> i32 {
        self.mapq
    }
    #[setter]
    fn set_mapq(&mut self, value: i32) {
        self.mapq = value;
    }
    #[getter]
    fn n_positions(&self) -> usize {
        self.positions.len()
    }
    fn has_insertions(&self) -> bool {
        !self.insertions.is_empty()
    }
}

/// Port of `pairloom.alignment.build_read_structure`.
#[pyfunction]
#[pyo3(signature = (tid, start0, is_reverse, mapq, query_name, cigar, seq, quals, read_group))]
fn build_structure(
    tid: i32,
    start0: i32,
    is_reverse: bool,
    mapq: i32,
    query_name: String,
    cigar: Vec<(u8, u32)>,
    seq: String,
    quals: Option<Vec<u8>>,
    read_group: Option<String>,
) -> PyResult<ReadStructureR> {
    let seq_up = seq.to_uppercase();
    let seqb = seq_up.as_bytes();
    let has_qualities = quals.is_some();
    let zeros = vec![0u8; seqb.len()];
    let quals_ref: &[u8] = quals.as_deref().unwrap_or(&zeros);
    let context = format!("read {query_name}");
    if quals.is_some() && quals_ref.len() != seqb.len() {
        return Err(data_invariant(format!(
            "quality length {} != sequence length {} for {context}",
            quals_ref.len(),
            seqb.len()
        )));
    }
    validate_cigar(&cigar, seqb.len(), &context)?;
    let cigar = normalize_cigar(&cigar);

    let mut positions: Vec<i32> = Vec::with_capacity(seqb.len());
    let mut bases: Vec<u8> = Vec::with_capacity(seqb.len());
    let mut quals_v: Vec<u8> = Vec::with_capacity(seqb.len());
    let mut deleted: BTreeSet<i32> = BTreeSet::new();
    let mut skipped: BTreeSet<i32> = BTreeSet::new();
    let mut insertions: BTreeMap<i32, (Vec<u8>, Vec<u8>)> = BTreeMap::new();
    let mut lead_bases: Vec<u8> = Vec::new();
    let mut lead_quals: Vec<u8> = Vec::new();
    let mut trail_bases: Vec<u8> = Vec::new();
    let mut trail_quals: Vec<u8> = Vec::new();
    let mut in_trail = false;

    let mut qi = 0usize;
    let mut ri = start0;
    let mut aligned_len: i32 = 0;
    let mut hard_clipped = false;
    for &(op, length) in &cigar {
        let len = length as usize;
        match op {
            op if is_aligned_op(op) => {
                for k in 0..len {
                    positions.push(ri + k as i32);
                    bases.push(seqb[qi + k]);
                    quals_v.push(quals_ref[qi + k]);
                }
                aligned_len += length as i32;
                qi += len;
                ri += length as i32;
                in_trail = true;
            }
            CIGAR_I => {
                insertions.insert(
                    ri,
                    (seqb[qi..qi + len].to_vec(), quals_ref[qi..qi + len].to_vec()),
                );
                qi += len;
            }
            CIGAR_D => {
                for k in 0..length as i32 {
                    deleted.insert(ri + k);
                }
                ri += length as i32;
            }
            CIGAR_N => {
                for k in 0..length as i32 {
                    skipped.insert(ri + k);
                }
                ri += length as i32;
            }
            CIGAR_S => {
                if in_trail {
                    trail_bases.extend_from_slice(&seqb[qi..qi + len]);
                    trail_quals.extend_from_slice(&quals_ref[qi..qi + len]);
                } else {
                    lead_bases.extend_from_slice(&seqb[qi..qi + len]);
                    lead_quals.extend_from_slice(&quals_ref[qi..qi + len]);
                }
                qi += len;
            }
            CIGAR_H => hard_clipped = true,
            _ => unreachable!("validate_cigar rejects other ops"),
        }
    }
    Ok(ReadStructureR {
        tid,
        start0,
        end0: ri,
        is_reverse,
        mapq,
        query_name,
        positions,
        bases,
        quals: quals_v,
        deleted,
        skipped,
        insertions,
        lead_clip_bases: lead_bases,
        lead_clip_quals: lead_quals,
        trail_clip_bases: trail_bases,
        trail_clip_quals: trail_quals,
        hard_clipped,
        has_qualities,
        read_group,
        aligned_len,
    })
}

// ---------------------------------------------------------------------------
// try_merge port
// ---------------------------------------------------------------------------

const R_NON_PRIMARY: &str = "non_primary_records";
const R_HARD_CLIP: &str = "hard_clipped_mate";
const R_NO_QUALITIES: &str = "qualities_missing";
const R_MAPQ: &str = "mapq_below_min";
const R_SKIP: &str = "reference_skip_unsupported";
const R_CONTIG: &str = "different_contig";
const R_ORIENTATION: &str = "orientation_incompatible";
const R_FRAGMENT: &str = "fragment_length_out_of_bounds";
const R_OVERLAP: &str = "insufficient_overlap";
const R_RG: &str = "read_group_conflict";
const R_GAP: &str = "fragment_gap";
const R_INDEL: &str = "indel_conflict";
const R_INSERTION: &str = "insertion_conflict";
const R_STRICT: &str = "strict_conflict";
const R_EMPTY: &str = "empty_consensus";

/// Port of `pairloom.pairing.leftmost_is_read1`.
fn leftmost_is_read1(pos1: i32, pos2: i32, rev1: bool, _rev2: bool, model: &str) -> bool {
    if pos1 != pos2 {
        return pos1 < pos2;
    }
    match model {
        "fr" => !rev1,
        "rf" => rev1,
        _ => true,
    }
}

/// Port of `pairloom.pairing.orientation_ok`.
fn orientation_ok(rev1: bool, rev2: bool, left_is_read1: bool, model: &str) -> PyResult<bool> {
    let (left_rev, right_rev) = if left_is_read1 {
        (rev1, rev2)
    } else {
        (rev2, rev1)
    };
    Ok(match model {
        "fr" => !left_rev && right_rev,
        "rf" => left_rev && !right_rev,
        "ff" => !rev1 && !rev2,
        "rr" => rev1 && rev2,
        "any" => true,
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown library orientation model {other:?}"
            )))
        }
    })
}

/// Port of `pairloom.merge._consensus_base`.
enum Consensus {
    None,
    Base(u8, u8),
    Strict,
}

fn consensus_base(
    in1: Option<(u8, u8)>,
    in2: Option<(u8, u8)>,
    min_quality_difference: u8,
    max_consensus_quality: u8,
    strict: bool,
) -> Consensus {
    match (in1, in2) {
        (None, None) => Consensus::None,
        (Some(x), None) => Consensus::Base(x.0, x.1),
        (None, Some(x)) => Consensus::Base(x.0, x.1),
        (Some((b1, q1)), Some((b2, q2))) => {
            if b1 == b2 {
                if is_acgt(b1) {
                    let q = q1.saturating_add(q2).min(max_consensus_quality).min(93);
                    return Consensus::Base(b1, q);
                }
                return Consensus::Base(b1, q1.min(q2));
            }
            if b1 == b'N' {
                return Consensus::Base(b2, q2);
            }
            if b2 == b'N' {
                return Consensus::Base(b1, q1);
            }
            let (hi, lo) = if q1 >= q2 {
                ((b1, q1), (b2, q2))
            } else {
                ((b2, q2), (b1, q1))
            };
            if hi.1.saturating_sub(lo.1) >= min_quality_difference {
                return Consensus::Base(hi.0, hi.1.min(max_consensus_quality));
            }
            if strict {
                return Consensus::Strict;
            }
            Consensus::Base(b'N', UNRESOLVED_QUALITY)
        }
    }
}

/// Port of `pairloom.merge._lead_owner` / `_trail_owner` selection.
fn lead_owner_is_1(s1: &ReadStructureR, s2: &ReadStructureR) -> bool {
    if s1.start0 != s2.start0 {
        s1.start0 < s2.start0
    } else {
        s1.lead_clip_bases.len() >= s2.lead_clip_bases.len()
    }
}
fn trail_owner_is_1(s1: &ReadStructureR, s2: &ReadStructureR) -> bool {
    if s1.end0 != s2.end0 {
        s1.end0 > s2.end0
    } else {
        s1.trail_clip_bases.len() >= s2.trail_clip_bases.len()
    }
}

fn covers(s: &ReadStructureR, boundary: i32) -> bool {
    s.start0 <= boundary && boundary <= s.end0
}

enum InsertionOutcome {
    None,
    Ins(Vec<u8>, Vec<u8>),
    Failure(&'static str),
}

/// Port of `pairloom.merge._merge_insertion`.
fn merge_insertion(boundary: i32, s1: &ReadStructureR, s2: &ReadStructureR) -> InsertionOutcome {
    let ins1 = s1.insertions.get(&boundary);
    let ins2 = s2.insertions.get(&boundary);
    match (ins1, ins2) {
        (None, None) => InsertionOutcome::None,
        (Some((b1, q1)), Some((b2, q2))) => {
            if b1 == b2 {
                let quals = q1
                    .iter()
                    .zip(q2.iter())
                    .map(|(a, b)| (*a).max(*b))
                    .collect();
                InsertionOutcome::Ins(b1.clone(), quals)
            } else {
                InsertionOutcome::Failure(R_INSERTION)
            }
        }
        (Some(ins), None) => {
            let other = s2;
            if covers(other, boundary) && other.start0 < boundary && boundary < other.end0 {
                return InsertionOutcome::Failure(R_INSERTION);
            }
            InsertionOutcome::Ins(ins.0.clone(), ins.1.clone())
        }
        (None, Some(ins)) => {
            let other = s1;
            if covers(other, boundary) && other.start0 < boundary && boundary < other.end0 {
                return InsertionOutcome::Failure(R_INSERTION);
            }
            InsertionOutcome::Ins(ins.0.clone(), ins.1.clone())
        }
    }
}

/// Successful merge outcome — field-for-field port of `MergeSuccess`.
#[pyclass]
struct MergeOutcomeR {
    #[pyo3(get)]
    ok: bool,
    #[pyo3(get)]
    reason: String,
    #[pyo3(get)]
    start0: i32,
    #[pyo3(get)]
    cigar: Vec<(u8, u32)>,
    #[pyo3(get)]
    seq: String,
    #[pyo3(get)]
    quals: Vec<u8>,
    #[pyo3(get)]
    mapq: i32,
    #[pyo3(get)]
    overlap_len: i32,
    #[pyo3(get)]
    read_group: Option<String>,
}

fn failure(reason: &str) -> MergeOutcomeR {
    MergeOutcomeR {
        ok: false,
        reason: reason.to_string(),
        start0: 0,
        cigar: Vec::new(),
        seq: String::new(),
        quals: Vec::new(),
        mapq: 0,
        overlap_len: 0,
        read_group: None,
    }
}

/// Port of `pairloom.merge.try_merge` (gate order and per-position walk
/// included, exactly).
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(signature = (s1, s2, min_mapq, min_overlap, min_fragment_length, max_fragment_length, min_quality_difference, max_consensus_quality, conflict_policy, library_orientation, merged_mapq, rescue_bypass_mapq, has_non_primary))]
fn merge_pair(
    s1: PyRef<'_, ReadStructureR>,
    s2: PyRef<'_, ReadStructureR>,
    min_mapq: i32,
    min_overlap: i32,
    min_fragment_length: i32,
    max_fragment_length: i32,
    min_quality_difference: u8,
    max_consensus_quality: u8,
    conflict_policy: &str,
    library_orientation: &str,
    merged_mapq: Option<i32>,
    rescue_bypass_mapq: bool,
    has_non_primary: bool,
) -> PyResult<MergeOutcomeR> {
    let (s1, s2) = (&*s1, &*s2);
    let strict = conflict_policy == "strict";
    if has_non_primary {
        return Ok(failure(R_NON_PRIMARY));
    }
    if s1.hard_clipped || s2.hard_clipped {
        return Ok(failure(R_HARD_CLIP));
    }
    if !s1.has_qualities || !s2.has_qualities {
        return Ok(failure(R_NO_QUALITIES));
    }
    if !rescue_bypass_mapq && (s1.mapq < min_mapq || s2.mapq < min_mapq) {
        return Ok(failure(R_MAPQ));
    }
    if !s1.skipped.is_empty() || !s2.skipped.is_empty() {
        return Ok(failure(R_SKIP));
    }
    if s1.tid != s2.tid {
        return Ok(failure(R_CONTIG));
    }
    let left_is_read1 = leftmost_is_read1(
        s1.start0,
        s2.start0,
        s1.is_reverse,
        s2.is_reverse,
        library_orientation,
    );
    if !orientation_ok(s1.is_reverse, s2.is_reverse, left_is_read1, library_orientation)? {
        return Ok(failure(R_ORIENTATION));
    }
    let start = s1.start0.min(s2.start0);
    let end = s1.end0.max(s2.end0);
    let span = end - start;
    if !(min_fragment_length..=max_fragment_length).contains(&span) {
        return Ok(failure(R_FRAGMENT));
    }
    // overlap = |positions1 ∩ positions2| via ascending two-pointer
    let mut overlap: i32 = 0;
    {
        let (mut i, mut j) = (0usize, 0usize);
        while i < s1.positions.len() && j < s2.positions.len() {
            match s1.positions[i].cmp(&s2.positions[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    overlap += 1;
                    i += 1;
                    j += 1;
                }
            }
        }
    }
    if overlap < min_overlap {
        return Ok(failure(R_OVERLAP));
    }
    if s1.read_group.is_some() && s2.read_group.is_some() && s1.read_group != s2.read_group {
        return Ok(failure(R_RG));
    }

    let mut seq_parts: Vec<u8> = Vec::with_capacity((span.max(0) + 8) as usize * 2);
    let mut qual_parts: Vec<u8> = Vec::with_capacity(seq_parts.capacity());
    let mut ops: Vec<(u8, u32)> = Vec::new();

    macro_rules! emit {
        ($op:expr, $length:expr, $bases:expr, $quals:expr) => {
            if $length > 0 {
                if let Some(last) = ops.last_mut() {
                    if last.0 == $op {
                        last.1 += $length;
                    } else {
                        ops.push(($op, $length));
                    }
                } else {
                    ops.push(($op, $length));
                }
                seq_parts.extend_from_slice(&$bases);
                qual_parts.extend_from_slice(&$quals);
            }
        };
    }

    let (lead, trail) = {
        let lead_is_1 = lead_owner_is_1(s1, s2);
        let trail_is_1 = trail_owner_is_1(s1, s2);
        (
            if lead_is_1 { s1 } else { s2 },
            if trail_is_1 { s1 } else { s2 },
        )
    };
    emit!(
        CIGAR_S,
        lead.lead_clip_bases.len() as u32,
        lead.lead_clip_bases,
        lead.lead_clip_quals
    );

    // per-position walk with two-pointer base lookups
    let (mut i1, mut i2) = (0usize, 0usize);
    for rpos in start..end {
        match merge_insertion(rpos, s1, s2) {
            InsertionOutcome::Failure(r) => return Ok(failure(r)),
            InsertionOutcome::Ins(bases, quals) => {
                emit!(CIGAR_I, bases.len() as u32, bases, quals);
            }
            InsertionOutcome::None => {}
        }
        while i1 < s1.positions.len() && s1.positions[i1] < rpos {
            i1 += 1;
        }
        while i2 < s2.positions.len() && s2.positions[i2] < rpos {
            i2 += 1;
        }
        let in1: Option<(u8, u8)> = if i1 < s1.positions.len() && s1.positions[i1] == rpos {
            Some((s1.bases[i1], s1.quals[i1]))
        } else {
            None
        };
        let in2: Option<(u8, u8)> = if i2 < s2.positions.len() && s2.positions[i2] == rpos {
            Some((s2.bases[i2], s2.quals[i2]))
        } else {
            None
        };
        let d1 = s1.deleted.contains(&rpos);
        let d2 = s2.deleted.contains(&rpos);
        if in1.is_none() && in2.is_none() {
            if d1 || d2 {
                emit!(CIGAR_D, 1, [0u8; 0], [0u8; 0]);
                continue;
            }
            return Ok(failure(R_GAP));
        }
        if (in1.is_some() && d2) || (in2.is_some() && d1) {
            return Ok(failure(R_INDEL));
        }
        match consensus_base(in1, in2, min_quality_difference, max_consensus_quality, strict) {
            Consensus::Strict => return Ok(failure(R_STRICT)),
            Consensus::None => return Ok(failure(R_GAP)),
            Consensus::Base(base, qual) => emit!(CIGAR_M, 1, [base], [qual]),
        }
    }
    match merge_insertion(end, s1, s2) {
        InsertionOutcome::Failure(r) => return Ok(failure(r)),
        InsertionOutcome::Ins(bases, quals) => {
            emit!(CIGAR_I, bases.len() as u32, bases, quals);
        }
        InsertionOutcome::None => {}
    }

    emit!(
        CIGAR_S,
        trail.trail_clip_bases.len() as u32,
        trail.trail_clip_bases,
        trail.trail_clip_quals
    );
    if !ops.iter().any(|&(op, _)| op == CIGAR_M || op == CIGAR_I || op == CIGAR_D) {
        return Ok(failure(R_EMPTY));
    }
    let read_group = s1.read_group.clone().or_else(|| s2.read_group.clone());
    Ok(MergeOutcomeR {
        ok: true,
        reason: String::new(),
        start0: start,
        cigar: normalize_cigar(&ops),
        seq: String::from_utf8(seq_parts)
            .map_err(|e| data_invariant(format!("non-UTF8 consensus sequence: {e}")))?,
        quals: qual_parts,
        mapq: merged_mapq.unwrap_or_else(|| s1.mapq.min(s2.mapq)),
        overlap_len: overlap,
        read_group,
    })
}

// ---------------------------------------------------------------------------
// aggregate_overlap port (blob transport)
// ---------------------------------------------------------------------------

/// Aggregate one template's eligible overlap observations.
///
/// Returns raw little-endian blobs (consumed by `merge_aggregates`):
/// (keys u64 = (tid << 32) | rpos, disg u8, refcode u8, altkey u64,
///  altbase u8, altstrand u8). Each key carries dp = 2 (both mates'
/// observations); every alt row is one observation (count 1).
#[pyfunction]
#[pyo3(signature = (s1, s2, tid, ref_seq, ref_offset, min_mapq, min_bq))]
fn aggregate_overlap<'py>(
    py: Python<'py>,
    s1: PyRef<'_, ReadStructureR>,
    s2: PyRef<'_, ReadStructureR>,
    tid: i32,
    ref_seq: &str,
    ref_offset: i32,
    min_mapq: i32,
    min_bq: u8,
) -> PyResult<(
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
)> {
    let (s1, s2) = (&*s1, &*s2);
    let refb = ref_seq.as_bytes();
    let tid64 = (tid as i64 as u64) << 32;
    let mut keys: Vec<u64> = Vec::new();
    let mut disg: Vec<u8> = Vec::new();
    let mut refcodes: Vec<u8> = Vec::new();
    let mut altkeys: Vec<u64> = Vec::new();
    let mut altbases: Vec<u8> = Vec::new();
    let mut altstrands: Vec<u8> = Vec::new();

    if s1.mapq >= min_mapq && s2.mapq >= min_mapq {
        let (mut i, mut j) = (0usize, 0usize);
        while i < s1.positions.len() && j < s2.positions.len() {
            let (p1, p2) = (s1.positions[i], s2.positions[j]);
            match p1.cmp(&p2) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    let (b1, q1) = (s1.bases[i], s1.quals[i]);
                    let (b2, q2) = (s2.bases[j], s2.quals[j]);
                    let ok1 = is_acgt(b1) && q1 >= min_bq;
                    let ok2 = is_acgt(b2) && q2 >= min_bq;
                    if ok1 && ok2 {
                        let rpos = p1;
                        let ref_base = refb
                            .get((rpos - ref_offset) as usize)
                            .copied()
                            .unwrap_or(0);
                        keys.push(tid64 | (rpos as u64 & 0xFFFF_FFFF));
                        let disagrees = b1 != b2;
                        disg.push(u8::from(disagrees));
                        refcodes.push(ref_base);
                        if b1 != ref_base {
                            altkeys.push(tid64 | (rpos as u64 & 0xFFFF_FFFF));
                            altbases.push(b1);
                            altstrands.push(u8::from(s1.is_reverse));
                        }
                        if b2 != ref_base {
                            altkeys.push(tid64 | (rpos as u64 & 0xFFFF_FFFF));
                            altbases.push(b2);
                            altstrands.push(u8::from(s2.is_reverse));
                        }
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
    }
    let keys_b = PyBytes::new(py, &to_le_bytes_u64(&keys));
    let disg_b = PyBytes::new(py, &disg);
    let ref_b = PyBytes::new(py, &refcodes);
    let ak_b = PyBytes::new(py, &to_le_bytes_u64(&altkeys));
    let ab_b = PyBytes::new(py, &altbases);
    let as_b = PyBytes::new(py, &altstrands);
    Ok((keys_b, disg_b, ref_b, ak_b, ab_b, as_b))
}

fn to_le_bytes_u64(v: &[u64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 8);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}
fn from_le_bytes_u64(buf: &[u8]) -> Vec<u64> {
    buf.chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

#[derive(Default)]
struct AggEntry {
    dp: u64,
    disg: bool,
    refcode: u8,
    seen_ref: bool,
    alts: BTreeMap<(u8, u8), u64>, // (base, strand) -> count
}

/// Merge per-template aggregation blobs for one chunk into chunk-level
/// blobs, ready for `FastDisagreementAccumulator.add_chunk`.
#[pyfunction]
#[pyo3(signature = (items,))]
fn merge_aggregates<'py>(
    py: Python<'py>,
    items: Vec<(
        Bound<'py, PyBytes>,
        Bound<'py, PyBytes>,
        Bound<'py, PyBytes>,
        Bound<'py, PyBytes>,
        Bound<'py, PyBytes>,
        Bound<'py, PyBytes>,
    )>,
) -> PyResult<(
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
    Bound<'py, PyBytes>,
)> {
    let mut entries: BTreeMap<u64, AggEntry> = BTreeMap::new();
    for (kb, db, rb, akb, abb, asb) in items {
        let keys = from_le_bytes_u64(kb.as_bytes());
        let disg = db.as_bytes();
        let refs = rb.as_bytes();
        for (idx, &key) in keys.iter().enumerate() {
            let e = entries.entry(key).or_default();
            e.dp += 2;
            e.disg |= disg[idx] != 0;
            if !e.seen_ref {
                e.refcode = refs[idx];
                e.seen_ref = true;
            } else if e.refcode != refs[idx] {
                return Err(data_invariant(format!(
                    "inconsistent reference base at tid {} pos {}",
                    key >> 32,
                    (key & 0xFFFF_FFFF) + 1
                )));
            }
        }
        let akeys = from_le_bytes_u64(akb.as_bytes());
        let abases = abb.as_bytes();
        let astrands = asb.as_bytes();
        for idx in 0..akeys.len() {
            let e = entries.entry(akeys[idx]).or_default();
            *e.alts.entry((abases[idx], astrands[idx])).or_insert(0) += 1;
        }
    }
    // NOTE: alt-only positions (dp = 0) cannot occur: every alt row's
    // position also carries a base observation in the same template.
    let mut keys: Vec<u64> = Vec::with_capacity(entries.len());
    let mut dp: Vec<u64> = Vec::with_capacity(entries.len());
    let mut disg: Vec<u8> = Vec::with_capacity(entries.len());
    let mut refs: Vec<u8> = Vec::with_capacity(entries.len());
    let mut ak: Vec<u64> = Vec::new();
    let mut ab: Vec<u8> = Vec::new();
    let mut astr: Vec<u8> = Vec::new();
    let mut ac: Vec<u64> = Vec::new();
    for (key, e) in entries {
        keys.push(key);
        dp.push(e.dp);
        disg.push(u8::from(e.disg));
        refs.push(e.refcode);
        for ((base, strand), count) in e.alts {
            ak.push(key);
            ab.push(base);
            astr.push(strand);
            ac.push(count);
        }
    }
    Ok((
        PyBytes::new(py, &to_le_bytes_u64(&keys)),
        PyBytes::new(py, &to_le_bytes_u64(&dp)),
        PyBytes::new(py, &disg),
        PyBytes::new(py, &refs),
        PyBytes::new(py, &to_le_bytes_u64(&ak)),
        PyBytes::new(py, &ab),
        PyBytes::new(py, &astr),
        PyBytes::new(py, &to_le_bytes_u64(&ac)),
    ))
}

/// Port of `pairloom.cigarops.compute_md_nm`: MD tag and NM count by
/// comparing bases to the reference (query in BAM storage orientation;
/// `ref_seq` covers exactly [pos0, pos0 + cigar_ref_len)).
#[pyfunction]
#[pyo3(signature = (pos0, cigar, query_seq, ref_seq))]
fn compute_md_nm(
    pos0: i32,
    cigar: Vec<(u8, u32)>,
    query_seq: &str,
    ref_seq: &str,
) -> PyResult<(String, u32)> {
    let need: u32 = cigar
        .iter()
        .filter(|&&(op, _)| is_aligned_op(op) || op == CIGAR_D || op == CIGAR_N)
        .map(|&(_, len)| len)
        .sum();
    let q = query_seq.as_bytes();
    let r = ref_seq.as_bytes();
    if (r.len() as u32) < need {
        return Err(data_invariant(format!(
            "reference slice too short: need {need} bases from pos {pos0}, got {}",
            r.len()
        )));
    }
    let mut md_parts: Vec<u8> = Vec::new();
    let mut run: u32 = 0;
    let mut nm: u32 = 0;
    let mut qi = 0usize;
    let mut ri = 0usize;
    let mut push_run = |md_parts: &mut Vec<u8>, run: &mut u32| {
        let s = run.to_string();
        md_parts.extend_from_slice(s.as_bytes());
        *run = 0;
    };
    for &(op, length) in &cigar {
        match op {
            op if is_aligned_op(op) => {
                for k in 0..length as usize {
                    if q[qi + k] == r[ri + k] {
                        run += 1;
                    } else {
                        push_run(&mut md_parts, &mut run);
                        md_parts.push(r[ri + k]);
                        nm += 1;
                    }
                }
                qi += length as usize;
                ri += length as usize;
            }
            CIGAR_I => {
                nm += length;
                qi += length as usize;
            }
            CIGAR_D => {
                push_run(&mut md_parts, &mut run);
                md_parts.push(b'^');
                md_parts.extend_from_slice(&r[ri..ri + length as usize]);
                nm += length;
                ri += length as usize;
            }
            CIGAR_N => {
                push_run(&mut md_parts, &mut run);
                md_parts.push(b'^');
                md_parts.extend_from_slice(&r[ri..ri + length as usize]);
                ri += length as usize;
            }
            CIGAR_S => qi += length as usize,
            CIGAR_H => {}
            _ => {
                return Err(data_invariant("unsupported CIGAR op in MD recompute".to_string()))
            }
        }
    }
    push_run(&mut md_parts, &mut run);
    Ok((
        String::from_utf8(md_parts)
            .map_err(|e| data_invariant(format!("non-UTF8 MD tag: {e}")))?,
        nm,
    ))
}

/// All ascending start offsets of `probe` in `haystack` (rescue scans).
#[pyfunction]
fn find_all(haystack: &str, probe: &str) -> Vec<i64> {
    if probe.is_empty() {
        return Vec::new();
    }
    haystack
        .match_indices(probe)
        .map(|(i, _)| i as i64)
        .collect()
}

/// Occurrences of `block` inside `window`, mapped by `lo` offset, capped at
/// `limit` + 1 entries so repetitive windows stay bounded (candidate
/// discovery sliding). Returns (positions, hit_limit).
#[pyfunction]
fn slide_window(window: &str, block: &str, lo: i64, limit: usize) -> (Vec<i64>, bool) {
    let hits = find_all(window, block);
    let hit_limit = hits.len() > limit;
    (hits.into_iter().map(|p| lo + p).take(limit).collect(), hit_limit)
}

#[pymodule]
fn pairloom_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(build_structure, m)?)?;
    m.add_function(wrap_pyfunction!(merge_pair, m)?)?;
    m.add_function(wrap_pyfunction!(aggregate_overlap, m)?)?;
    m.add_function(wrap_pyfunction!(merge_aggregates, m)?)?;
    m.add_function(wrap_pyfunction!(compute_md_nm, m)?)?;
    m.add_function(wrap_pyfunction!(find_all, m)?)?;
    m.add_function(wrap_pyfunction!(slide_window, m)?)?;
    m.add("UNRESOLVED_QUALITY", UNRESOLVED_QUALITY)?;
    Ok(())
}
