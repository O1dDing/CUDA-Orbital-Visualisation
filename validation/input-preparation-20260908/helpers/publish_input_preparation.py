"""Publish only reviewed preparation artifacts and verify their exact Git bytes."""
from pathlib import Path
import collections
import ctypes
import hashlib
import json
import shutil
import subprocess
import sys

repo=Path(r'F:\Dev\cov-native-validation-20260905')
workspace=Path(__file__).resolve().parent.parent
base=workspace/'outputs/cov-complete-validation-20260906'
relative='validation/input-preparation-20260908'
target=repo/relative
sys.path.insert(0,str(repo/'tests'))
from validation_process import physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
git=r'F:\Dev\Git\cmd\git.exe'
already_staged=subprocess.check_output([git,'diff','--cached','--name-only'],cwd=repo,text=True).splitlines()
assert all(p.startswith(relative+'/') for p in already_staged)
assert not subprocess.check_output([git,'diff','--name-only'],cwd=repo)
for root in [base,target]:
    path=root/'benchmark-geometry-preparation-v1/README.md'
    text=path.read_text(encoding='utf-8')
    text=text.replace('数据集 CC BY 4.0。','数据集 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode)。')
    path.write_text(text,encoding='utf-8',newline='\n')
for name in ['GMTKN55-v1-LICENSE.txt','GMTKN55-v1-LICENSE.txt.retrieval.json']:
    shutil.copy2(base/'external-primary-sources-v1'/name,target/'external-primary-sources-v1'/name)
for name in ['checkpoint_input_preparation.py','publish_input_preparation.py']:
    shutil.copy2(workspace/'work'/name,target/'helpers'/name)
files=[p for p in sorted(target.rglob('*')) if p.is_file() and p.name!='files-sha256.json']
assert not any(p.suffix.lower() in {'.pdf','.html','.pyc','.exe','.dll'} for p in files)
assert sum(p.stat().st_size for p in files)<16*1024*1024
manifest={str(p.relative_to(target)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(target/'files-sha256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
subprocess.run([git,'add','--force','--',*[relative+'/'+name for name in manifest],relative+'/files-sha256.json'],cwd=repo,check=True)
changed=subprocess.check_output([git,'diff','--cached','--name-only'],cwd=repo,text=True).splitlines()
assert len(changed)==len(manifest)+1 and all(p.startswith(relative+'/') for p in changed)
check=subprocess.run([git,'diff','--cached','--check'],cwd=repo,text=True,capture_output=True)
# Preserve the three reviewed Markdown line endings in the author's exact license.
if check.returncode:
    (base/'input-preparation-whitespace-review.txt').write_text(check.stdout+check.stderr,encoding='utf-8')
    headers=[line for line in check.stdout.splitlines() if not line.startswith('+')]
    expected=[relative+'/external-primary-sources-v1/GMTKN55-v1-LICENSE.txt:'+str(n)+': trailing whitespace.' for n in [12,15,16]]
    assert headers==expected and not check.stderr
    license_file=target/'external-primary-sources-v1/GMTKN55-v1-LICENSE.txt'
    assert hashlib.sha256(license_file.read_bytes()).hexdigest()=='29853614116e554998fb67ca0654c69dadedaf70960344027acdfa8ee42f0a68'
print(json.dumps(dict(files=len(files),bytes=sum(p.stat().st_size for p in files),
                      extensions=dict(collections.Counter(p.suffix for p in files)))),flush=True)
subprocess.run([git,'commit','-m','docs(validation): preserve density controls and external coordinate sources'],cwd=repo,check=True,stdout=subprocess.DEVNULL)
subprocess.run([git,'-c','http.version=HTTP/1.1','push','origin','test/fchk-validation-native'],cwd=repo,check=True)
subprocess.run([sys.executable,str(workspace/'work/verify_publication_git_bytes.py'),
                '--relative',relative,'--receipt',str(base/'github-input-preparation-git-bytes.json')],cwd=workspace,check=True)
