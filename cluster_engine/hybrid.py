"""Serial HOG refinement with bounded resources, progress and validated checkpoints."""
import hashlib
import json
import time
from pathlib import Path

from .data import partition


def _hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    temporary.replace(path)


def refine_hogs(method, groups, unsupported, params, context, work, sequence_groups):
    args = context['args']
    threads = min(getattr(args, 'hybrid_threads', 1), args.threads)
    timeout = getattr(args, 'hybrid_timeout', 300.0)
    subcontext = dict(context, sequence_threads=threads, command_timeout=timeout)
    candidate_dir = Path(context.get('candidate_directory', work))
    checkpoint_dir = candidate_dir/'hog_checkpoints'
    progress_file = candidate_dir/'hybrid_progress.json'
    # In normal CLI use the fingerprint covers code, tools, sequence file and
    # arguments. Tests/embedded callers without one do not reuse checkpoints.
    fingerprint = context.get('fingerprint')
    source_by_gene = {r['gene_ID']: r['HOG'] for r in
        context.get('manifest', {}).get('orthofinder_source', {}).get('assignments', []) if r['HOG']}
    progress = dict(method=method, total_groups=len(groups), completed_groups=0,
                    cached_groups=0, completed_genes=0, threads_per_hog=threads,
                    seconds_per_hog_limit=timeout, state='running', current=None)
    refined, missing_all = [], list(unsupported)
    started = time.monotonic()
    for i, group in enumerate(groups):
        group = sorted(group)
        hog = source_by_gene.get(group[0], 'unresolved' if group[0] in unsupported else 'HOG')
        label = '%s HOG %d/%d (%s; %d genes; %d thread%s)' % (
            method, i+1, len(groups), hog, len(group), threads, '' if threads == 1 else 's')
        group_work = Path(work)/('hog_%06d' % i)
        progress['current'] = dict(index=i+1, HOG=hog, gene_count=len(group), work_directory=str(group_work))
        progress['elapsed_seconds'] = time.monotonic()-started
        _atomic_json(progress_file, progress)
        print(label + ' starting', flush=True)
        group_started = time.monotonic()
        try:
            key = _hash(dict(fingerprint=fingerprint, method=method, params=params, threads=threads,
                             sequences=[(g, hashlib.sha256(context['seqs'][g].encode()).hexdigest()) for g in group]))
            checkpoint = checkpoint_dir/(key+'.json')
            cached = False
            if fingerprint and checkpoint.exists():
                record = json.loads(checkpoint.read_text(encoding='utf-8'))
                payload = record['result']
                if record['key'] != key or record['result_sha256'] != _hash(payload):
                    raise ValueError('HOG checkpoint checksum mismatch: ' + str(checkpoint))
                sub, missing = payload['groups'], payload['unsupported']
                partition(sub, group)
                if len(set(missing)) != len(missing) or not set(missing) <= set(group):
                    raise ValueError('invalid HOG checkpoint unsupported genes: ' + str(checkpoint))
                cached = True
            else:
                sub, missing = sequence_groups(method.split('-')[1], group, params, subcontext, group_work)
                partition(sub, group)
                if len(set(missing)) != len(missing) or not set(missing) <= set(group):
                    raise ValueError('invalid HOG subdivision unsupported genes')
                if fingerprint:
                    payload = dict(groups=sub, unsupported=missing)
                    _atomic_json(checkpoint, dict(key=key, result=payload, result_sha256=_hash(payload)))
            refined.extend(sub)
            missing_all.extend(missing)
            progress['completed_groups'] = i+1
            progress['completed_genes'] += len(group)
            progress['cached_groups'] += int(cached)
            progress['elapsed_seconds'] = time.monotonic()-started
            progress['current']['status'] = 'cached' if cached else 'complete'
            _atomic_json(progress_file, progress)
            print(label + (' cached' if cached else ' done in %.2fs' % (time.monotonic()-group_started)), flush=True)
        except BaseException as exc:
            progress.update(state='interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed',
                            error=str(exc) or type(exc).__name__, elapsed_seconds=time.monotonic()-started)
            _atomic_json(progress_file, progress)
            raise
    progress.update(state='complete', current=None, elapsed_seconds=time.monotonic()-started)
    _atomic_json(progress_file, progress)
    return refined, sorted(set(missing_all))
