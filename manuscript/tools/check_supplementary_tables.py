"""Checks the supplementary workbook against the manuscript text (usage, from the repository root:
python manuscript/tools/check_supplementary_tables.py <xlsx> <manuscript.md>).

Tables are found by their titles, so the checks do not depend on the numbering; a separate check confirms that each
sheet's title matches the item with the same number in the draft's supplementary list."""
import re
import sys

import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
xlsx, draft = sys.argv[1], sys.argv[2]
sheets = pd.read_excel(xlsx, sheet_name=None, header=None)
text = open(draft, encoding="utf-8").read()
fails = 0
TITLES = {name: str(df.iloc[0, 0]) for name, df in sheets.items() if name != "Contents"}


def check(label, ok):
    global fails
    fails += 0 if ok else 1
    print(("OK   " if ok else "FAIL ") + label)


def by_title(fragment: str) -> tuple[str, pd.DataFrame]:
    hits = [n for n, t in TITLES.items() if fragment in t]
    if len(hits) != 1:
        sys.exit(f"title fragment {fragment!r} matches {hits}")
    return hits[0], pd.read_excel(xlsx, sheet_name=hits[0], header=4)


def words(t: str) -> set:
    return set(re.findall(r"[a-z0-9]+", t.lower())) - {"the", "and", "of", "on", "in", "by", "per", "a", "each", "part"}


# numbering: sheet Sn title vs draft list item n
block = text[text.index("## Supplementary tables\n"):]
block = block[: block.index("\n---")]
items = {int(m.group(1)): m.group(2) for m in re.finditer(r"^(\d+)\. (.*)$", block, flags=re.M)}
for n, item in sorted(items.items()):
    scores = {name: len(words(item) & words(t)) / len(words(item) | words(t)) for name, t in TITLES.items()}
    best = max(scores, key=scores.get)
    check(f"list item {n} ({item[:50]}...) best matches sheet {best}", best == f"S{n}")
check(f"workbook has {len(TITLES)} tables, draft list has {len(items)}", len(TITLES) == len(items))
# Table 2 of the draft against the by-region calibration table
name, reg_t = by_title("Calibration by region")
reg = {"Donor": "donor region", "Acceptor": "acceptor region", "Near intronic": "near intronic",
       "Deep intronic": "deep intronic"}
cols = ["Supporting pathogenic", "Moderate pathogenic", "Strong pathogenic", "Supporting benign", "Moderate benign"]
for line in text.split("\n"):
    m = re.match(r"\| (Donor|Acceptor|Near intronic|Deep intronic) \| (AVI|SpliceAI) \| (.*)\|$", line)
    if not m:
        continue
    cells = [c.strip() for c in m.group(3).split("|")][1:]
    for col, cell in zip(cols, cells):
        sub = reg_t[(reg_t["tool"] == m.group(2)) & (reg_t["region"] == reg[m.group(1)]) & (reg_t["level"] == col.lower())]
        if cell == "not reached":
            check(f"Table 2 {m.group(1)} {m.group(2)} {col}: not reached / absent in {name}", sub.empty)
            continue
        flag = "yes" if cell.endswith(", S") else "no"
        lr = float(cell.split(" (")[0].replace(",", ""))
        ok = (not sub.empty) and sub["supported"].iloc[0] == flag and abs(float(sub["likelihood ratio"].iloc[0]) - lr) <= 0.01 * max(1, lr)
        check(f"Table 2 {m.group(1)} {m.group(2)} {col}: {cell} vs {name} {sub[['likelihood ratio','supported']].values.tolist()}", ok)
# key numbers
name, t = by_title("paired AUROC differences, AVI minus each tool")
x = t[(t["set"] == "main set") & (t["region"].str.startswith("main region")) & (t["comparison"] == "AVI − SpliceAI")]
check(f"{name} main AVI-SpliceAI {x[['AUROC difference','95% low','95% high']].values.tolist()} ~ -0.005 (-0.017, 0.006)",
      round(x["AUROC difference"].iloc[0], 3) == -0.005 and round(x["95% low"].iloc[0], 3) == -0.017)
name, t = by_title("within each evidence level")
x = t[(t["tool"] == "AVI") & (t["level"] == "strong pathogenic")]
check(f"{name} AVI strong {x[['likelihood ratio','95% low','95% high','supported']].values.tolist()} ~ 397 (294-569)",
      round(x["likelihood ratio"].iloc[0]) == 397)
name, t = by_title("thresholds derived on the calibration sample")
x = t[t["tool"] == "CADD"]["pathogenic threshold (score ≥)"].tolist()[:3]
check(f"{name} CADD thresholds {x} = 17.56, 19.18, 23.4", x == [17.56, 19.18, 23.4])
name, t = by_title("pooled estimates and leave-one-dataset-out")
x = t[(t["comparison"] == "AVI − SpliceAI") & t["dataset or pooled estimate"].str.startswith("pooled (")]
check(f"{name} pooled AVI-SpliceAI {x[['AUROC difference','Hartung–Knapp low','Hartung–Knapp high']].values.tolist()}",
      round(x["AUROC difference"].iloc[0], 3) == -0.089)
name, t = by_title("Part 2 datasets")
check(f"{name} analysed counts {t['analysed variants'].tolist()[:6]} sum of main = 334,230",
      int(t[t["role"] == "main"]["analysed variants"].sum()) == 334230)
name, t = by_title("flow from the ClinVar")
check(f"{name} flow contains 87,666", (t == 87666).any().any())
# the draft's references to the numbered tables point at the right content
refs = {"intronic CADD thresholds derived here": "thresholds derived on the calibration sample",
        "held in every supporting analysis": "paired AUROC differences, AVI minus each tool",
        "likelihood ratios of all variants beyond each threshold": "beyond each threshold",
        "definition of intermediate variants are listed in": "Part 2 datasets",
        "Part 2 used six main datasets": "Part 2 datasets",
        "in cumulative form": "beyond each threshold"}
for phrase, title in refs.items():
    i = text.index(phrase)
    num = int(re.search(r"Supplementary Table (\d+)", text[i:i + 400]).group(1))
    name, _ = by_title(title)
    check(f"'{phrase}' cites Supplementary Table {num}, sheet {name}", name == f"S{num}")

# definitions, coverage notes and flags
name, t = by_title("beyond each threshold")
diff = []
name7, t7 = by_title("within each evidence level")
for _, x in t.iterrows():
    y = t7[(t7["tool"] == x["tool"]) & (t7["level"] == x["level"])]
    if not y.empty and y["supported"].iloc[0] != x["supported"]:
        diff.append((x["tool"], x["level"]))
check(f"{name} vs {name7}: verdicts differ only for phyloP 241-way and phastCons 470-way supporting pathogenic {diff}",
      sorted(diff) == [("phastCons 470-way", "supporting pathogenic"), ("phyloP 241-way", "supporting pathogenic")])
i = text.index("likelihood ratios of all variants beyond each threshold")
check("Methods sentence names phyloP 241-way and phastCons 470-way", "phyloP 241-way and phastCons 470-way" in text[i:i + 500])
olds = []
for n, df in sheets.items():
    if n == "Contents":
        continue
    cells = df.iloc[4:].astype(str)
    for bad in ("inf", "unbounded", ">=", "<="):
        if cells.apply(lambda c: c.str.fullmatch(r"-?inf|unbounded") if bad in ("inf", "unbounded") else c.str.contains(bad, regex=False)).any().any():
            olds.append((n, bad))
check(f"no old bound notation in table cells {olds}", not olds)
name, t = by_title("coverage per tool")
check(f"{name} Pangolin rows carry the sample note", t[t["tool"] == "Pangolin"]["note"].str.contains("sample").all())
name, t = by_title("discrimination by set, region and tool")
check(f"{name} has no row without scored variants", not ((t["P/LP"] == 0) & (t["B/LB"] == 0)).any())
check(f"{name} flags {int((t['fewer than 5 P/LP or B/LB'] == 'yes').sum())} small cells",
      (t["fewer than 5 P/LP or B/LB"] == "yes").sum() == ((t["P/LP"] < 5) | (t["B/LB"] < 5)).sum())
contents = sheets["Contents"].astype(str)
terms = set(contents.iloc[:, 0])
needed = set(t["set"]) - {"main set"} | {"main set", "pooled test set", "AVI", "SpliceAI", "Pangolin", "Distance"}
missing = sorted(x for x in needed if not any(x == term or x in term for term in terms))
check(f"Contents defines every set and the key tools {missing}", not missing)

# effect measures, counts not analysed and notation
hy = sum(df.astype(str).apply(lambda col: col.str.contains("AVI - |mutant - WT", regex=True)).sum().sum() for df in sheets.values())
check(f"comparison labels use the minus sign (cells with a hyphen: {hy})", hy == 0)
name, t3 = by_title("Part 2 datasets")
em = t3.loc[t3["dataset"].str.startswith("COMPASS"), "effect measure (used for the Spearman correlation)"].iloc[0]
check(f"{name} COMPASS effect measure is one phrase", ";" not in em and "(cells:" not in em)
sa = t3.set_index("dataset")["set aside before the counts (design or data rules)"].astype(str)
check(f"{name} COMPASS set-aside names the ClinVar and Geuvadis subsets and 1,534 others",
      "ClinVar subset (1,885)" in sa.filter(like="COMPASS").iloc[0] and "1,534" in sa.filter(like="COMPASS").iloc[0])
check(f"{name} OpenSplice set-aside names its 36 duplicate SNVs",
      sa.filter(like="OpenSplice (Quarantani").iloc[0].startswith("36 duplicate SNVs (14 alleles"))
check(f"{name} COMPASS set-aside splits 1,534 into 1,426 reference-PSI and 108 duplicates",
      "1,534 other SNVs: 1,426" in sa.filter(like="COMPASS").iloc[0] and "108 duplicates" in sa.filter(like="COMPASS").iloc[0])
med = t3.set_index("dataset")["removed: median-based estimates without a test"].astype(str)
check(f"{name} MST1R median-based estimates removed = 106", med.filter(like="MST1R").iloc[0] == "106")


def num(v):
    try:
        return int(str(v).split(" ")[0].replace(",", ""))
    except ValueError:
        return 0


bad_rows = []
for _, x in t3.iterrows():
    if str(x["role"]) == "reported separately":
        continue
    expect = num(x["binary-labelled SNVs"]) - num(x["removed: Part 1 new labels"]) - num(x["removed: developers' benchmark data"]) - num(x["removed: median-based estimates without a test"])
    if expect != int(x["analysed variants"]):
        bad_rows.append((x["dataset"][:20], expect, x["analysed variants"]))
check(f"{name} binary - exclusions = analysed in every row {bad_rows}", not bad_rows)
rules = " ".join(t3["label rule (authors' calls or the authors' threshold)"].astype(str)) + " " + " ".join(t3["effect measure (used for the Spearman correlation)"].astype(str))
import re as _re
check(f"{name} rule texts use the minus sign before numbers and in KTS", not _re.search(r"[≤≥<>]\s?-\d|KTS-", rules))
name, t7 = by_title("within each evidence level")
x = t7[(t7["tool"] == "phastCons 470-way") & (t7["level"] == "supporting pathogenic")]["interval method"].iloc[0]
check(f"{name} phastCons supporting pathogenic method = {x!r}", x == "empty score interval")
# no per-variant content
bad = [n for n, df in sheets.items() if df.astype(str).apply(lambda c: c.str.contains(r"\b(?:chr)?\d{1,2}:\d{4,}:[ACGT]", regex=True)).any().any()]
check(f"no allele identifiers in any sheet {bad}", not bad)
print("FAILS:", fails)
