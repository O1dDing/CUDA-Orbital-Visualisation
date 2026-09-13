from pathlib import Path
from collections import Counter
import hashlib
import json
import shutil
import time

workspace=Path(__file__).resolve().parent.parent
source=workspace/'outputs/cov-complete-validation-20260906/external-preparation-v1'
target=Path(r'F:\Dev\cov-native-validation-20260905\validation\external-preparation-20260906')
summary=json.loads((source/'identity-review-v2/summary.json').read_text(encoding='utf-8'))
assert summary['all_metadata_reviews_terminal'] and summary['formal_external_cases']==0
matched=[r for r in summary['records'] if r['lookup_status']=='unique-metadata-match']
composition_new=sum(not r['original_cases_with_same_element_counts'] for r in matched)
families=dict(Counter(r['family'] for r in summary['records']))
table=['| 准备编号 | 候选 | 分子式 / 电荷 | 当前身份状态 | 用途 |',
       '|---|---|---|---|---|']
labels={'unique-metadata-match':'身份字段相符；几何待建','lookup-failed-retained':'名称查询无匹配，保留',
        'ambiguous-or-mismatched-metadata':'数据库字段不符，保留反例','source-coordinate-work-pending':'按原始结构来源继续核对'}
for r in summary['records']:
    table.append('| '+ ' | '.join([r['candidate_id'],r['query'],r['expected_formula']+' / '+str(r['expected_charge']),
                                  labels.get(r['lookup_status'],r['lookup_status']),r['purpose']])+' |')
body=f'''# 外部集合准备：67 个候选，尚未冻结正式集合

54 个候选的公开身份字段经电荷和元素计数核对后相符，其中 {composition_new} 个的元素计数与原库全部 273 案例不同。另有 9 个名称查询无匹配、1 个数据库字段与目标不符、3 个按原始结构论文继续核对。全部原始响应、错误响应和复核过程保留。正式外部分子数仍为 **0**；没有启动 Gaussian，也没有形成可验收的外部 FCHK。

查询按 [PubChem PUG REST 官方接口](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)顺序进行，并限制请求速率。保存 CID、InChIKey、连接 SMILES、分子式、电荷、响应原文和 SHA-256。第一版读取器未支持 PubChem 分子式尾部的电荷表示，第二次复核直接重读原始响应并显式比较电荷，未改写外部数据。

同分子式不等于同一结构。cubane/COT、pentatetraene/spiropentadiene、benzyl/tropylium 分别保留原库碰撞记录，正式入库前还须核对图结构。54 个相符候选的连接 InChIKey 首段之间没有重复；金属配合物的数据库标准化可能拆开配体，此结果不能证明配位连接、几何或电子态。

六甲基钨的名称响应给出了断开的甲基组分及 −6 电荷，与目标中性 W(CH3)6 不符，已作为输入反例保留，不能直接用于 Gaussian。Ni(PMe3)4 名称查询同时返回质子化与中性记录，逐条核对后仅保留中性匹配的身份候选。[Ni(PF3)4 气相结构原始研究](https://pubs.rsc.org/en/content/articlelanding/1970/c2/c29700000595)和 [Rh/PMe3 配合物原始结构研究](https://pubs.rsc.org/en/content/articlelanding/1980/dt/dt9800000511)提供后续坐标核查入口，当前尚未取得并验收这些坐标。

甲基过氧自由基采用 [NIST CCCBDB 的身份标识](https://cccbdb.nist.gov/exp2x.asp?casno=2143580&charge=0)补查；硝普盐阴离子采用 [PubChem 对应阴离子条目](https://pubchem.ncbi.nlm.nih.gov/compound/11963622)补查。保留最初的查无结果或中性错误匹配，不能让名称别名静默改变目标电荷。

正式集合还需补充弱复合物/环境模型，例如氨—氟化氢、HCN—水、甲酸二聚体、甲醇—水及 CO2—水；这五项目前是用途候选，尚未另计入 67 项身份清单。新增分子中的纯/笛卡尔 g 表示作为同一分子的输入变体覆盖，不额外增加分子计数。

后续顺序：核对真实结构或适用初始几何 → 明确电荷、电子态及环境用途 → 对原库和新增集合进行化学图去重 → 明确冻结成员、方法和判据 → 在用户解除 Gaussian 暂停并满足前置条件后计算及正式验证。准备中的失败条目继续保留；不得把身份匹配数当作计算完成或 COV 通过数。

{chr(10).join(table)}

源码数值/布局正在 NUM-FIX-005 固定版本复验，此目录只增加外部准备资料，不改变该轮输入、算法或判据。REF-001 保持用户暂停状态。
'''
(source/'README.md').write_text(body,encoding='utf-8')
record={'snapshot_epoch':time.time(),'candidate_count':67,'identity_counts':summary['counts'],
        'matched_composition_disjoint_from_original':composition_new,'family_counts':families,
        'formal_external_cases':0,'gaussian_started':False,'REF-001':'paused-by-user',
        'preparation_directory':str(source),'review_sha256':hashlib.sha256((source/'identity-review-v2/summary.json').read_bytes()).hexdigest()}
(source/'status.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
target.mkdir(exist_ok=True)
for path in source.rglob('*'):
    if path.is_file():
        out=target/path.relative_to(source);out.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,out)
(target/'.gitattributes').write_text('* -text\n',encoding='utf-8')
for name in ('prepare_external_identities.py','review_external_identities.py','publish_external_preparation.py'):
    shutil.copy2(workspace/'work'/name,target/name)
shutil.copy2(workspace/'outputs/general-fixes-20260906/next-iteration-scope-contract.md',target/'next-iteration-scope-contract.md')
hashes={str(p.relative_to(target)):hashlib.sha256(p.read_bytes()).hexdigest() for p in target.rglob('*')
        if p.is_file() and p.name!='files-sha256.json'}
(target/'files-sha256.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
state_path=workspace/'outputs/cov-complete-validation-20260906/阶段状态.json'
state=json.loads(state_path.read_text(encoding='utf-8'));state['external_preparation']=record;state['snapshot_epoch']=time.time()
temporary=state_path.with_name(state_path.name+'.tmp');temporary.write_text(json.dumps(state,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');temporary.replace(state_path)
checkpoint=workspace/'outputs/cov-complete-validation-20260906/执行检查点.md'
heading='## 外部候选身份准备已完成第一批复核'
previous=checkpoint.read_text(encoding='utf-8')
if heading not in previous:
    checkpoint.write_text(previous+'\n\n'+heading+'\n\n67 个准备候选，54 个身份字段相符，51 个相符候选的元素计数与原库不同。9 项名称无匹配、1 项字段不符、3 项结构源待核对均保留。正式外部数 0、Gaussian 未启动、REF-001 暂停。详见 external-preparation-v1/README.md 和 identity-review-v2/summary.json。\n',encoding='utf-8')
print(json.dumps({'publication_directory':str(target),'files':len(hashes),'matched_composition_disjoint':composition_new}))
