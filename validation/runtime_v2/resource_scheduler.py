"""Atomic host-wide whole-tier reservations, including utilities and RAM holds.

The ledger is durable provenance; OS process and Job identities decide liveness.
No elapsed-time expiry can discard a living writer. A reservation spans a whole
stage and its sequential probes; each native tree also has a durable sublease.
"""
from contextlib import contextmanager
from pathlib import Path
import os
import threading
import time
import uuid

from policy import VERSION, core_budget
from state_store import atomic_json, read_json, lease, LeaseOccupied

_threads = threading.RLock()
_local = threading.local()


class AdmissionDenied(RuntimeError):
    pass


def default_path():
    return Path(os.environ.get('LOCALAPPDATA', str(Path.home()))) / 'COV/adaptive-physical-resources.json'


class Ledger:
    def __init__(self, path, host, owner, alive, job_alive=lambda name: False):
        self.path = Path(path)
        self.host, self.owner = host, owner
        self.alive, self.job_alive = alive, job_alive
        self.budget = core_budget(host['physical_cores'])

    @contextmanager
    def transaction(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with _threads:
            deadline = time.monotonic() + 10
            while True:
                try:
                    guard = lease(self.path.with_suffix('.lock'))
                    guard.__enter__()
                    break
                except LeaseOccupied:
                    if time.monotonic() >= deadline:
                        raise
                    time.sleep(.01)
            try:
                value = read_json(self.path, {'schema': 1, 'active': {}, 'events': []})
                yield value
                atomic_json(self.path, value)
            finally:
                guard.__exit__(None, None, None)

    def event(self, value, kind, **fields):
        value['events'].append(dict(sequence=len(value['events']) + 1,
            epoch=time.time(), event=kind, policy_version=VERSION, **fields))

    def reconcile(self, value):
        for token, row in list(value['active'].items()):
            if self.alive(row['owner']):
                continue
            tree = row.get('tree')
            if tree and (self.job_alive(tree['job_name']) or
                         (tree.get('process') and self.alive(tree['process']))):
                row['state'] = 'orphan_tree_still_alive'
                continue
            self.event(value, 'crash_reconciled_release', token=token, reservation=row)
            del value['active'][token]

    def reserve(self, requested, *, label, budget=None, metadata=None):
        budget = self.budget if budget is None else budget
        if (isinstance(requested, bool) or not isinstance(requested, int) or
                not 1 <= requested <= budget <= self.budget):
            raise AdmissionDenied('Whole request does not fit the physical host budget')
        with self.transaction() as value:
            self.reconcile(value)
            active = value['active']
            # A more conservative existing coordinator's budget remains binding.
            effective = min([budget] + [r['budget'] for r in active.values()])
            if len(active) >= 2 or sum(r['allocated'] for r in active.values()) + requested > effective:
                raise AdmissionDenied('Physical budget or two-job limit occupied')
            token = uuid.uuid4().hex
            row = {'token': token, 'owner': self.owner, 'label': label,
                'requested': requested, 'allocated': requested, 'budget': budget,
                'concurrent_jobs_at_admission': len(active) + 1, 'policy_version': VERSION,
                'host': self.host, 'metadata': metadata or {}, 'state': 'reserved', 'tree': None}
            active[token] = row
            self.event(value, 'allocate', token=token, reservation=row.copy())
        return Reservation(self, token, row)

    def update_tree(self, token, tree):
        with self.transaction() as value:
            row = value['active'][token]
            if tree is not None and row['tree'] is not None and row['tree']['job_name'] != tree['job_name']:
                raise AdmissionDenied('A reservation cannot launch concurrent auxiliary trees')
            row['tree'] = tree
            self.event(value, 'tree_start' if tree else 'tree_exit_confirmed', token=token, tree=tree)

    def hold(self, token, held):
        with self.transaction() as value:
            row = value['active'][token]
            state = 'ram_hold_reservation_retained' if held else 'reserved'
            if row['state'] != state:
                row['state'] = state
                self.event(value, state, token=token,
                    reason='The living tree retains all allocated cores and its job slot throughout RAM hold')

    def release(self, token):
        with self.transaction() as value:
            row = value['active'][token]
            if row['tree'] is not None:
                raise AdmissionDenied('Cannot release before the entire owned tree confirms exit')
            self.event(value, 'release', token=token, reservation=row)
            del value['active'][token]


class Reservation:
    def __init__(self, ledger, token, record):
        self.ledger, self.token, self.record = ledger, token, record

    @contextmanager
    def scope(self):
        if getattr(_local, 'reservation', None) is not None:
            raise AdmissionDenied('Nested resource scope must reuse its existing reservation')
        _local.reservation = self
        try:
            yield self
        finally:
            _local.reservation = None

    def release(self):
        self.ledger.release(self.token)


def current():
    return getattr(_local, 'reservation', None)
