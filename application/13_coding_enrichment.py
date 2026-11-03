"""Coding-consequence enrichment among each method's confident calls.

The seven positive controls in validated_variants.csv are too few to separate
methods. This check needs no truth set and scales to every call: of the
variants a method declares causal at PIP >= 0.95, what share are
protein-altering? The reference is the set of GWAS lead variants, i.e. what
you would have picked with no fine-mapping at all.

Consequences come from the Ensembl VEP REST endpoint, queried BY rsID so no
GRCh37/38 liftover enters, and are cached in vep_consequences.json so the
table reproduces offline.

  python 13_coding_enrichment.py --export ~/Downloads/chol_results_v2/export

Two limits, both of which have to travel with the number:

  1. Protein-altering consequence is a PROXY for causality. It is biased
     against genuinely regulatory causal variants, which at lipid loci are a
     large share of the truth -- SORT1 is the canonical example and it is
     non-coding.
  2. It is partly circular for a functionally informed method. PolyFun's
     `combined` prior is built on baselineLF, which contains explicit coding
     and conserved categories, so that column is partly instructed to upweight
     coding variants. `polyfun_sniff` is not: it sees only the Borzoi PCs,
     the same information fb_xregion's head sees, which makes sniff-vs-
     fb_xregion the comparison worth reading.

The call sets also overlap each other and the lead set heavily (39-43% of each
method's calls ARE its block's lead), so the intervals below are Wilson
binomial intervals for each set on its own and NO test between sets is
offered -- they would not be independent samples.
"""
import argparse, json, math, os, sys, time, urllib.request, urllib.error
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--export", required=True)
ap.add_argument("--cache", default=os.path.join(HERE, "vep_consequences.json"))
ap.add_argument("--refresh", action="store_true", help="re-query Ensembl for ids not cached")
a = ap.parse_args()
E = os.path.expanduser(a.export)

# Sequence Ontology terms that change the protein product.
CODING = {"missense_variant", "stop_gained", "stop_lost", "start_lost",
          "frameshift_variant", "inframe_insertion", "inframe_deletion",
          "splice_acceptor_variant", "splice_donor_variant",
          "protein_altering_variant", "incomplete_terminal_codon_variant"}
ORDER = ["polyfun_combined", "polyfun_sniff", "fb_xregion", "funmap",
         "fb_prefix", "fb_memomean", "beatrice"]

conf = pd.read_csv(os.path.join(E, "confident_calls.csv"))
S    = pd.read_csv(os.path.join(E, "region_method_summary.csv"))
leads = sorted(set(S[S.sig].lead_snp))
need  = sorted(set(conf.SNP) | set(leads))

V = json.load(open(a.cache)) if os.path.exists(a.cache) else {}
missing = [i for i in need if i not in V]
if missing and a.refresh:
    URL, B = "https://rest.ensembl.org/vep/human/id", 60
    for i in range(0, len(missing), B):
        chunk = missing[i:i + B]
        req = urllib.request.Request(
            URL, data=json.dumps({"ids": chunk}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        for _ in range(4):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    for rec in json.load(r):
                        V[rec["input"]] = rec.get("most_severe_consequence", "NA")
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(float(e.headers.get("Retry-After", 2))); continue
                print(f"  HTTP {e.code} on batch {i//B}", file=sys.stderr); break
            except Exception as e:
                print(f"  {type(e).__name__} on batch {i//B}, retrying", file=sys.stderr)
                time.sleep(2)
        print(f"  {len(V)} cached", flush=True)
    json.dump(V, open(a.cache, "w"))
    missing = [i for i in need if i not in V]
if missing:
    print(f"WARNING {len(missing)} of {len(need)} ids have no consequence; "
          f"re-run with --refresh" if not a.refresh else
          f"WARNING {len(missing)} ids still unannotated")

def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p, d = k / n, 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)

conf["csq"] = conf.SNP.map(V)
lead_csq = pd.Series(leads).map(V).dropna()
base = lead_csq.isin(CODING).mean()

print(f"{'set':<20}{'calls':>7}{'prot-alt':>10}{'%':>8}{'95% CI':>18}"
      f"{'x lead rate':>13}{'blocks':>8}")
for m in ORDER:
    g = conf[conf.method == m].dropna(subset=["csq"])
    if g.empty:
        continue
    k, n = int(g.csq.isin(CODING).sum()), len(g)
    lo, hi = wilson(k, n)
    print(f"{m:<20}{n:>7}{k:>10}{100*k/n:>7.1f}%   [{100*lo:>4.1f}, {100*hi:>4.1f}]"
          f"{k/n/base:>12.1f}x{g[g.csq.isin(CODING)].region_id.nunique():>8}")
lo, hi = wilson(int(base * len(lead_csq)), len(lead_csq))
print(f"{'[GWAS leads]':<20}{len(lead_csq):>7}{int(base*len(lead_csq)):>10}"
      f"{100*base:>7.1f}%   [{100*lo:>4.1f}, {100*hi:>4.1f}]{1.0:>12.1f}x")

print("\nshare of each method's calls that ARE its block's lead variant:")
ls = S.set_index(["region_id", "method"]).lead_snp
for m in ORDER:
    g = conf[conf.method == m]
    if g.empty:
        continue
    n = sum(1 for r, s in zip(g.region_id, g.SNP) if ls.get((r, m)) == s)
    print(f"  {m:<20}{n:>4} of {len(g):>4}  ({100*n/len(g):.0f}%)")
