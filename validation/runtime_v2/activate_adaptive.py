"""Activate a verified immutable candidate without rewriting old computations.

The candidate receipt supplies both old and new identities. All three entry
leases and an OS process check are required; capability receipts are reverified.
This updates only the active pointer, shared config and launcher metadata.
"""
from contextlib import ExitStack
from pathlib import Path
import argparse
import os
import time

import resume
import windows_job
from migration import verify_evidence
from state_store import lease, read_json, atomic_json


def activate(home, candidate, receipt_path):
    candidate=read_json(candidate)
    config=candidate['candidate_config']
    prepared=candidate['prepared']
    bundle=Path(prepared['bundle']).resolve()
    if not bundle.is_relative_to(home.resolve()/'runtime/releases'):
        raise ValueError('Candidate bundle escaped this deployment')
    with ExitStack() as stack:
        for path in (Path(os.environ['LOCALAPPDATA'])/'COV/gaussian-runtime-v2.lock',
                     Path(config['work_root'])/'jobs/supervisor.lock', home/'helper-operation.lock'):
            stack.enter_context(lease(path))
        if read_json(home/'runtime/active.json')!=candidate['previous_active'] or read_json(home/'runtime/config.json')!=candidate['previous_config']:
            raise ValueError('Active deployment changed after candidate preparation')
        data=resume.FrozenData(config)
        if windows_job.existing_gaussian(data.gaussian):
            raise RuntimeError('Gaussian still running; refuse activation')
        engine=resume.Engine(config,data,windows_job)
        if engine.identity!=candidate['runtime_identity'] or read_json(engine.root/'binding.json')['runtime_identity']!=engine.identity:
            raise ValueError('Candidate source and imported runtime identity differ')
        for method in ('RPBE1PBE','UPBE1PBE'):engine.require_capability('ram_pause',method)
        imported=read_json(engine.root/'import-summary.json')
        verify_evidence(imported['source_evidence'])
        manifest=read_json(bundle/'manifest.json')
        for name,digest in manifest['files'].items():
            if resume.sha(bundle/'validation'/name)!=digest:
                raise ValueError('Immutable candidate changed: '+name)
        # Each byte that contributes to execution identity must match this code.
        for path in resume.HERE.glob('*.py'):
            if resume.sha(path)!=resume.sha(bundle/'validation/runtime_v2'/path.name):
                raise ValueError('Candidate bundle misses final source: '+path.name)
        if engine.mode()!='pause':raise ValueError('Candidate must be idle and paused')
        receipt=dict(candidate,activation_started_epoch=time.time(),completed=False,
            originals_modified=False,calculation_started=False,source_import_counts=imported['counts'])
        atomic_json(receipt_path,receipt)
        # If interrupted between replacements, the old/new identity mismatch
        # fails closed. The prepared receipt supplies both exact rollback values.
        atomic_json(home/'runtime/config.json',config)
        atomic_json(home/'runtime/active.json',prepared)
        launcher=read_json(home/'runtime/launcher.json')
        launcher.update(helper_version='3.0.0',runtime_identity=engine.identity)
        atomic_json(home/'runtime/launcher.json',launcher)
        receipt.update(completed=True,activated=True,completed_epoch=time.time())
        atomic_json(receipt_path,receipt)
        print('Activated immutable adaptive runtime '+engine.identity,flush=True)
        return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--home',type=Path,required=True)
    parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--receipt',type=Path,required=True)
    args=parser.parse_args();activate(args.home,args.candidate,args.receipt)
