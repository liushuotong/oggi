//! OGGI anchor DP. The caller retains legacy NumPy sorting and p-value semantics.
use std::slice;

#[repr(C)]
pub struct View {
    scores: *const f64,
    used: *const u64,
    offsets: *const usize,
    paths: *const usize,
    paths_len: usize,
}

pub struct ResultData {
    scores: Vec<f64>,
    used: Vec<u64>,
    offsets: Vec<usize>,
    paths: Vec<usize>,
}

fn calculate(
    x: &[i64],
    y: &[i64],
    grading: &[f64],
    order: &[usize],
    reverse: bool,
    mg1: i64,
    mg2: i64,
    penalty: f64,
) -> ResultData {
    let n = x.len();
    let mut scores = grading.to_vec();
    let mut used = vec![0; n];
    // Immutable links preserve path snapshots even for unsorted forward input.
    let mut links: Vec<(usize, Option<usize>)> = (0..n).map(|i| (i, None)).collect();
    let mut heads: Vec<usize> = (0..n).collect();
    let mut by_x: Vec<usize> = (0..n).collect();
    by_x.sort_by_key(|&i| x[i]);
    let mut rank = vec![0; n];
    for (r, &i) in order.iter().enumerate() {
        rank[i] = r;
    }
    let mut candidates = Vec::new();
    for &i in order {
        let (lo, hi) = if reverse {
            (
                by_x.partition_point(|&j| x[j] <= x[i].saturating_sub(mg1)),
                by_x.partition_point(|&j| x[j] < x[i]),
            )
        } else {
            (
                by_x.partition_point(|&j| x[j] <= x[i]),
                by_x.partition_point(|&j| x[j] < x[i].saturating_add(mg1)),
            )
        };
        candidates.clear();
        if lo < hi {
            candidates.extend(
                by_x[lo..hi]
                    .iter()
                    .copied()
                    .filter(|&j| y[j] > y[i] && y[j] < y[i].saturating_add(mg2)),
            );
        }
        // Range indexing must not alter the historical neighbour visit order.
        candidates.sort_unstable_by_key(|&j| rank[j]);
        let mut old_x = x[i];
        let mut gap = mg2;
        for &j in &candidates {
            let dy = y[j] - y[i];
            if dy > gap && if reverse { x[j] < old_x } else { x[j] > old_x } {
                break;
            }
            let dx = if reverse { x[i] - x[j] } else { x[j] - x[i] };
            let step = grading[j] + (dx as f64 + dy as f64) * penalty;
            let score = step + scores[i];
            if step > 0.0 && scores[j] < score {
                scores[j] = score;
                used[i] += 1;
                used[j] += 1;
                links.push((j, Some(heads[i])));
                heads[j] = links.len() - 1;
                gap = gap.min(dy);
                old_x = x[j];
            }
        }
    }
    let mut paths = Vec::new();
    let mut offsets = vec![0];
    for head in heads {
        let start = paths.len();
        let mut cursor = Some(head);
        while let Some(k) = cursor {
            paths.push(links[k].0);
            cursor = links[k].1;
        }
        paths[start..].reverse();
        offsets.push(paths.len());
    }
    ResultData {
        scores,
        used,
        offsets,
        paths,
    }
}

#[no_mangle]
pub extern "C" fn oggi_subcoli_abi() -> u32 {
    1
}

/// Pointers must reference n elements; order must be a permutation of 0..n.
#[no_mangle]
pub unsafe extern "C" fn oggi_subcoli_score(
    n: usize,
    x: *const i64,
    y: *const i64,
    grading: *const f64,
    order: *const usize,
    reverse: bool,
    mg1: i64,
    mg2: i64,
    penalty: f64,
    view: *mut View,
) -> *mut ResultData {
    if x.is_null() || y.is_null() || grading.is_null() || order.is_null() || view.is_null() {
        return std::ptr::null_mut();
    }
    let result = std::panic::catch_unwind(|| {
        let order = slice::from_raw_parts(order, n);
        let mut seen = vec![false; n];
        for &i in order {
            if i >= n || seen[i] {
                return None;
            }
            seen[i] = true;
        }
        Some(Box::new(calculate(
            slice::from_raw_parts(x, n),
            slice::from_raw_parts(y, n),
            slice::from_raw_parts(grading, n),
            order,
            reverse,
            mg1,
            mg2,
            penalty,
        )))
    });
    match result {
        Ok(Some(data)) => {
            *view = View {
                scores: data.scores.as_ptr(),
                used: data.used.as_ptr(),
                offsets: data.offsets.as_ptr(),
                paths: data.paths.as_ptr(),
                paths_len: data.paths.len(),
            };
            Box::into_raw(data)
        }
        _ => std::ptr::null_mut(),
    }
}

/// Release a handle exactly once, after copying its view.
#[no_mangle]
pub unsafe extern "C" fn oggi_subcoli_free(data: *mut ResultData) {
    if !data.is_null() {
        drop(Box::from_raw(data));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn forward_and_reverse() {
        let forward = calculate(
            &[0, 1, 2],
            &[0, 1, 2],
            &[50.; 3],
            &[0, 1, 2],
            false,
            40,
            40,
            -1.,
        );
        assert_eq!(forward.scores, vec![50., 98., 146.]);
        assert_eq!(
            &forward.paths[forward.offsets[2]..forward.offsets[3]],
            &[0, 1, 2]
        );
        let reverse = calculate(
            &[0, 1, 2],
            &[2, 1, 0],
            &[50.; 3],
            &[2, 1, 0],
            true,
            40,
            40,
            -1.,
        );
        assert_eq!(reverse.scores, vec![146., 98., 50.]);
    }
    #[test]
    fn exclusive_gap_boundary() {
        let data = calculate(&[0, 40], &[0, 1], &[50.; 2], &[0, 1], false, 40, 40, -1.);
        assert_eq!(data.used, vec![0, 0]);
    }
}
