import numpy as np
import pandas as pd


def resolve_backend(backend="auto"):
    """Select a prebuilt native library; never download or compile at runtime."""
    if backend not in ("auto", "python", "rust"):
        raise ValueError("subcoli backend must be auto, python or rust")
    if backend == "python":
        return backend
    library = _native_library()
    if library is not None:
        return "rust"
    if backend == "rust":
        raise RuntimeError("Rust subcoli library not found. Build subcoli_rs with "
                           "cargo build --release, or set OGGI_SUBCOLI_LIB to the library file.")
    return "python"


def _native_library():
    import ctypes as ct
    import os
    import sys
    from pathlib import Path
    explicit = os.environ.get("OGGI_SUBCOLI_LIB")
    name = ("oggi_subcoli.dll" if sys.platform == "win32" else
            "liboggi_subcoli.dylib" if sys.platform == "darwin" else "liboggi_subcoli.so")
    root = Path(__file__).resolve().parent
    candidates = ([Path(explicit)] if explicit else
                  [root / name, root / "subcoli_rs" / "target" / "release" / name])
    for path in candidates:
        if not path.is_file():
            continue
        key = str(path.resolve())
        if key in _NATIVE_CACHE:
            return _NATIVE_CACHE[key]
        lib = ct.CDLL(key)
        lib.oggi_subcoli_abi.restype = ct.c_uint32
        if lib.oggi_subcoli_abi() != 1:
            raise RuntimeError("Unsupported subcoli Rust ABI: " + key)
        lib.oggi_subcoli_score.argtypes = [ct.c_size_t, ct.c_void_p, ct.c_void_p,
            ct.c_void_p, ct.c_void_p, ct.c_bool, ct.c_int64, ct.c_int64,
            ct.c_double, ct.c_void_p]
        lib.oggi_subcoli_score.restype = ct.c_void_p
        lib.oggi_subcoli_free.argtypes = [ct.c_void_p]
        lib.oggi_subcoli_free.restype = None
        _NATIVE_CACHE[key] = lib
        return lib
    if explicit:
        raise FileNotFoundError("OGGI_SUBCOLI_LIB does not exist: " + explicit)
    return None


_NATIVE_CACHE = {}


def run_blocks(options, points, backend="auto"):
    engine = resolve_backend(backend)
    return (RustCollinearity if engine == "rust" else collinearity)(options, points).run()

class collinearity:
    def __init__(self, options, points):
        # Default values
        self.gap_penalty = -1
        self.over_length = 0
        self.mg1 = 40
        self.mg2 = 40
        self.pvalue = 1
        self.over_gap = 3
        self.points = points
        self.p_value = 0
        self.coverage_ratio = 0.8
        
        # Set user-defined options
        for k, v in options:
            setattr(self, str(k), v)

        # Initialize grading and mg values
        self.grading = [50, 40, 25] if not hasattr(self, 'grading') else [int(k) for k in self.grading.split(',')]
        self.mg1, self.mg2 = [40, 40] if not hasattr(self, 'mg') else [int(k) for k in self.mg.split(',')]

        # Convert string values to floats
        self.pvalue = float(self.pvalue)
        self.coverage_ratio = float(self.coverage_ratio)

    def get_matrix(self):
        """Initialize the matrix for the collinearity points."""
        self.points['usedtimes1'] = 0
        self.points['usedtimes2'] = 0
        self.points['times'] = 1
        self.points['score1'] = self.points['grading']
        self.points['score2'] = self.points['grading']
        self.points['path1'] = self.points.index.to_numpy().reshape(len(self.points), 1).tolist()
        self.points['path2'] = self.points['path1']
        self.points_init = self.points.copy()
        self.mat_points = self.points

    def run(self):
        """Run the main collinearity processing."""
        self.get_matrix()
        self.score_matrix()
        data = []

        # Process points for maxPath in the positive direction
        points1 = self.points[['loc1', 'loc2', 'score1', 'path1', 'usedtimes1']].sort_values(by=['score1'], ascending=False)
        points1.drop(index=points1[points1['usedtimes1'] < 1].index, inplace=True)
        points1.columns = ['loc1', 'loc2', 'score', 'path', 'usedtimes']
        
        while (self.over_length >= self.over_gap or len(points1) >= self.over_gap):
            if self.max_path(points1):
                if self.p_value > self.pvalue:
                    continue
                data.append([self.path, self.p_value, self.score])

        # Process points for maxPath in the negative direction
        points2 = self.points[['loc1', 'loc2', 'score2', 'path2', 'usedtimes2']].sort_values(by=['score2'], ascending=False)
        points2.drop(index=points2[points2['usedtimes2'] < 1].index, inplace=True)
        points2.columns = ['loc1', 'loc2', 'score', 'path', 'usedtimes']

        while (self.over_length >= self.over_gap) or (len(points2) >= self.over_gap):
            if self.max_path(points2):
                if self.p_value > self.pvalue:
                    continue
                data.append([self.path, self.p_value, self.score])

        return data

    def score_matrix(self):
        """Calculate the scoring matrix for the points."""
        for index, row, col in self.points[['loc1', 'loc2']].itertuples():
            # Get points within a certain range
            points = self.points[(self.points['loc1'] > row) & 
                                 (self.points['loc2'] > col) & 
                                 (self.points['loc1'] < row + self.mg1) & 
                                 (self.points['loc2'] < col + self.mg2)]
            
            row_i_old, gap = row, self.mg2
            for index_ij, row_i, col_j, grading in points[['loc1', 'loc2', 'grading']].itertuples():
                if col_j - col > gap and row_i > row_i_old:
                    break
                score = grading + (row_i - row + col_j - col) * self.gap_penalty
                score1 = score + self.points.at[index, 'score1']
                if score > 0 and self.points.at[index_ij, 'score1'] < score1:
                    self.points.at[index_ij, 'score1'] = score1
                    self.points.at[index, 'usedtimes1'] += 1
                    self.points.at[index_ij, 'usedtimes1'] += 1
                    self.points.at[index_ij, 'path1'] = self.points.at[index, 'path1'] + [index_ij]
                    gap = min(col_j - col, gap)
                    row_i_old = row_i

        # Reverse processing to handle negative direction
        points_reverse = self.points.sort_values(by=['loc1', 'loc2'], ascending=[False, True])
        for index, row, col in points_reverse[['loc1', 'loc2']].itertuples():
            points = points_reverse[(points_reverse['loc1'] < row) & 
                                    (points_reverse['loc2'] > col) & 
                                    (points_reverse['loc1'] > row - self.mg1) & 
                                    (points_reverse['loc2'] < col + self.mg2)]
            
            row_i_old, gap = row, self.mg2
            for index_ij, row_i, col_j, grading in points[['loc1', 'loc2', 'grading']].itertuples():
                if col_j - col > gap and row_i < row_i_old:
                    break
                score = grading + (row - row_i + col_j - col) * self.gap_penalty
                score2 = score + self.points.at[index, 'score2']
                if score > 0 and self.points.at[index_ij, 'score2'] < score2:
                    self.points.at[index_ij, 'score2'] = score2
                    self.points.at[index, 'usedtimes2'] += 1
                    self.points.at[index_ij, 'usedtimes2'] += 1
                    self.points.at[index_ij, 'path2'] = self.points.at[index, 'path2'] + [index_ij]
                    gap = min(col_j - col, gap)
                    row_i_old = row_i

    def max_path(self, points):
        """Find the maximum path for the given points."""
        if len(points) == 0:
            self.over_length = 0
            return False
        
        # Initialize path score and index
        self.score, self.path_index = points.loc[points.index[0], ['score', 'path']]
        self.path = points[points.index.isin(self.path_index)]
        self.over_length = len(self.path_index)
        
        # Check if the block overlaps with other blocks
        if self.over_length >= self.over_gap and len(self.path) / self.over_length > self.coverage_ratio:
            points.drop(index=self.path.index, inplace=True)
            [loc1_min, loc2_min], [loc1_max, loc2_max] = self.path[['loc1', 'loc2']].agg(['min', 'max']).to_numpy()

            # Calculate p-value
            gap_init = self.points_init[(loc1_min <= self.points_init['loc1']) & 
                                        (self.points_init['loc1'] <= loc1_max) & 
                                        (loc2_min <= self.points_init['loc2']) & 
                                        (self.points_init['loc2'] <= loc2_max)].copy()
            
            self.p_value = self.p_value_estimated(gap_init, loc1_max - loc1_min + 1, loc2_max - loc2_min + 1)
            self.path = self.path.sort_values(by=['loc1'], ascending=[True])[['loc1', 'loc2']]
            return True
        else:
            points.drop(index=points.index[0], inplace=True)
        return False

    def p_value_estimated(self, gap, L1, L2):
        """Estimate p-value based on the given gap and lengths."""
        N1 = gap['times'].sum()
        N = len(gap)
        self.points_init.loc[gap.index, 'times'] += 1
        m = len(self.path)
        a = (1 - self.score / m / self.grading[0]) * (N1 - m + 1) / N * (L1 - m + 1) * (L2 - m + 1) / L1 / L2
        return round(a, 4)


class RustCollinearity(collinearity):
    """Native DP with array-based extraction and legacy score ordering/rounding.

    Unlike the reference class, the input DataFrame is not mutated. Public block
    frames retain its index and the same loc1/loc2 columns.
    """

    def run(self):
        import ctypes as ct
        n = len(self.points)
        if not n:
            return []
        if not self.points.index.is_unique:
            raise ValueError("anchor index must be unique")
        coordinates = self.points[["loc1", "loc2"]].to_numpy()
        if (not np.isfinite(coordinates).all() or
                not np.equal(coordinates, np.floor(coordinates)).all() or
                np.abs(coordinates).max() >= 2**52):
            raise ValueError("anchor coordinates must be finite integers smaller than 2**52")
        if not np.isfinite(self.points["grading"].to_numpy(dtype=float)).all():
            raise ValueError("anchor grading must be finite")
        if self.over_gap < 1 or self.grading[0] <= 0 or not np.isfinite(self.gap_penalty):
            raise ValueError("invalid subcoli scoring parameters")
        lib = _native_library()
        if lib is None:
            raise RuntimeError("Rust subcoli library not found")

        class View(ct.Structure):
            _fields_ = [("scores", ct.POINTER(ct.c_double)), ("used", ct.POINTER(ct.c_uint64)),
                        ("offsets", ct.POINTER(ct.c_size_t)), ("paths", ct.POINTER(ct.c_size_t)),
                        ("paths_len", ct.c_size_t)]

        x = np.ascontiguousarray(coordinates[:, 0], dtype=np.int64)
        y = np.ascontiguousarray(coordinates[:, 1], dtype=np.int64)
        grading = np.ascontiguousarray(self.points["grading"], dtype=np.float64)
        reverse_order = self.points.index.get_indexer(
            self.points.sort_values(["loc1", "loc2"], ascending=[False, True]).index)
        times = np.ones(n, dtype=np.int64)
        data = []
        over_length = self.over_length
        for reverse, order in ((False, np.arange(n)), (True, reverse_order)):
            order = np.ascontiguousarray(order, dtype=np.uintp)
            view = View()
            handle = lib.oggi_subcoli_score(n, x.ctypes.data, y.ctypes.data, grading.ctypes.data,
                order.ctypes.data, reverse, self.mg1, self.mg2, self.gap_penalty, ct.byref(view))
            if not handle:
                raise RuntimeError("Rust subcoli scoring failed")
            try:
                score = np.ctypeslib.as_array(view.scores, shape=(n,)).copy()
                used = np.ctypeslib.as_array(view.used, shape=(n,)).copy()
                offsets = np.ctypeslib.as_array(view.offsets, shape=(n + 1,)).copy()
                paths = np.ctypeslib.as_array(view.paths, shape=(view.paths_len,)).copy()
            finally:
                lib.oggi_subcoli_free(handle)
            if (pd.api.types.is_integer_dtype(self.points["grading"].dtype)
                    and float(self.gap_penalty).is_integer()):
                score = score.astype(np.int64)
            # pandas' default unstable tie order is part of the existing results.
            ranked = pd.Series(score).sort_values(ascending=False).index.to_numpy()
            ranked = ranked[used[ranked] >= 1]
            active = np.zeros(n, dtype=bool)
            active[ranked] = True
            remaining = len(ranked)
            cursor = 0
            while over_length >= self.over_gap or remaining >= self.over_gap:
                while cursor < len(ranked) and not active[ranked[cursor]]:
                    cursor += 1
                if cursor == len(ranked):
                    over_length = 0
                    break
                head = ranked[cursor]
                path = paths[offsets[head]:offsets[head + 1]]
                over_length = len(path)
                membership = np.zeros(n, dtype=bool)
                membership[path] = True
                selected = ranked[active[ranked] & membership[ranked]]
                m = len(selected)
                if over_length >= self.over_gap and m / over_length > self.coverage_ratio:
                    active[selected] = False
                    remaining -= m
                    xmin, xmax = x[selected].min(), x[selected].max()
                    ymin, ymax = y[selected].min(), y[selected].max()
                    inside = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
                    N1, N = times[inside].sum(), int(inside.sum())
                    times[inside] += 1
                    L1, L2 = xmax - xmin + 1, ymax - ymin + 1
                    pv = round((1 - score[head] / m / self.grading[0]) * (N1 - m + 1) / N
                               * (L1 - m + 1) * (L2 - m + 1) / L1 / L2, 4)
                    if pv <= self.pvalue:
                        block = self.points.iloc[selected].sort_values("loc1")[["loc1", "loc2"]].copy()
                        data.append([block, pv, score[head]])
                else:
                    active[head] = False
                    remaining -= 1
        return data
