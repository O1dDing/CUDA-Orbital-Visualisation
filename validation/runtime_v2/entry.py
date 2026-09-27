"""Shared launcher, including the original Resume helper command-line interface."""
from pathlib import Path
import argparse
import sys
import resume


def translate(argv):
    if '--legacy' not in argv:
        return argv or ['menu']
    argv = [arg for arg in argv if arg != '--legacy']
    parser = argparse.ArgumentParser(description='COV unified runtime (legacy Resume arguments supported)')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--limit', type=int, default=1)
    parser.add_argument('--workers', type=int, choices=(1, 2), default=1)
    parser.add_argument('--case', action='append', default=[])
    args = parser.parse_args(argv)
    if not args.execute:
        return ['menu']
    result = ['run', '--limit', str(args.limit), '--workers', str(args.workers)]
    for case in args.case:
        result.extend(['--case', case])
    return result


if __name__ == '__main__':
    # The compatibility stub supplies the shared configuration explicitly.
    config, raw = sys.argv[1], sys.argv[2:]
    sys.argv = [str(Path(resume.__file__)), '--config', config] + translate(raw)
    resume.main()
