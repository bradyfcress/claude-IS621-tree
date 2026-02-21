#!/usr/bin/env python3
"""Local genome search: download assemblies, predict proteins, BLAST + HMM search.

V3 replaces NCBI protein DB searches (scripts 02/05) with:
  1. Download genome assemblies for selected species from NCBI FTP
  2. Predict proteins with Prodigal (standard mode for complete genomes)
  3. BLAST E. coli references against predicted proteins
  4. hmmsearch with IS110 HMM profiles for divergent IS110 detection
  5. tBLASTn IS621/IS911 against nucleotide for unannotated elements
  6. Validate coverage meets expectations, update metadata.csv

Reads:
    data/taxonomy.csv               (species list from 03v2)
    data/metadata.csv               (existing metadata)
    data/ecoli_proteins.fasta       (5 reference proteins)
    data/hmm/all_is110_pfam.hmm     (IS110 Pfam HMM profiles)
    data/gtdb_cache/*_taxonomy.tsv.gz

Writes:
    data/genomes/{accession}/       (downloaded genomes + predicted proteins)
    data/targeted_cache/genome_search_cache.json
    data/metadata.csv               (overwritten with local search results)
"""

import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlretrieve

from Bio import Entrez, SeqIO

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
TAXONOMY_CSV = Path("data/taxonomy.csv")
METADATA_CSV = Path("data/metadata.csv")
ECOLI_FASTA = Path("data/ecoli_proteins.fasta")
GTDB_DIR = Path("data/gtdb_cache")
HMM_FILE = Path("data/hmm/all_is110_pfam.hmm")
GENOMES_DIR = Path("data/genomes")
CACHE_DIR = Path("data/targeted_cache")
CACHE_FILE = CACHE_DIR / "genome_search_cache.json"
ASSEMBLY_CACHE = CACHE_DIR / "assembly_info.json"

Entrez.email = "is621tree@example.com"
REQUEST_DELAY = 0.4  # seconds between NCBI calls

PROTEINS = ["IS621", "EF_Tu", "SecY", "MuA", "IS911"]
# E-value thresholds: stringent for universal proteins, relaxed for IS elements
EVALUE = {"EF_Tu": "1e-10", "SecY": "0.01", "MuA": "1e-5",
          "IS621": "1e-5", "IS911": "1e-5"}
HMM_EVALUE = 1e-5

# Tools
BLASTP = "blastp"
TBLASTN = "tblastn"
MAKEBLASTDB = "makeblastdb"
HMMSEARCH = os.path.expanduser("~/miniforge3/bin/hmmsearch")
PRODIGAL = "prodigal"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def sanitize_tip_label(name: str) -> str:
    """Make a tree-safe tip label (must match 03v2 logic exactly)."""
    label = re.sub(r"[^A-Za-z0-9_]", "_", name)
    label = re.sub(r"_+", "_", label)
    return label.strip("_")[:60]


def load_cache(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def save_cache(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


# ---------------------------------------------------------------------------
# [1/8] Load data
# ---------------------------------------------------------------------------
def load_data() -> tuple[list[dict], dict[str, str], dict[str, str]]:
    """Load taxonomy, GTDB genome mapping, and reference proteins.

    Returns (taxonomy_rows, tip_to_accession, ref_proteins).
    """
    # Load taxonomy.csv
    with open(TAXONOMY_CSV) as f:
        taxonomy = list(csv.DictReader(f))
    print(f"  Loaded {len(taxonomy)} species from {TAXONOMY_CSV}", flush=True)

    # Parse GTDB to get Tip_ID → assembly accession
    tip_to_accn = {}
    for gz_file in sorted(GTDB_DIR.glob("*_taxonomy.tsv.gz")):
        with gzip.open(gz_file, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) != 2:
                    continue
                genome_id = parts[0]
                ranks = parts[1].split(";")
                if len(ranks) != 7:
                    continue
                species = ranks[6].strip()
                if species.startswith("s__"):
                    species = species[3:]
                if not species:
                    continue
                tip = sanitize_tip_label(species)
                accn = genome_id
                for pfx in ("GB_", "RS_"):
                    if accn.startswith(pfx):
                        accn = accn[len(pfx):]
                tip_to_accn[tip] = accn

    matched = sum(1 for t in taxonomy if t["Tip_ID"] in tip_to_accn)
    print(f"  GTDB accession mapping: {matched}/{len(taxonomy)} species", flush=True)

    # Load reference proteins
    ref_proteins = {}
    for record in SeqIO.parse(str(ECOLI_FASTA), "fasta"):
        for protein in PROTEINS:
            if protein in record.id:
                ref_proteins[protein] = str(record.seq)
    print(f"  Reference proteins: {list(ref_proteins.keys())}", flush=True)

    return taxonomy, tip_to_accn, ref_proteins


# ---------------------------------------------------------------------------
# [2/8] Resolve FTP paths
# ---------------------------------------------------------------------------
def resolve_ftp_paths(
    taxonomy: list[dict],
    tip_to_accn: dict[str, str],
) -> dict[str, dict]:
    """Resolve NCBI FTP paths for each assembly.

    Returns {tip_id: {accession, ftp_path, source}}.
    """
    assembly_cache = load_cache(ASSEMBLY_CACHE)
    results = {}
    to_resolve = []

    for row in taxonomy:
        tip = row["Tip_ID"]
        accn = tip_to_accn.get(tip)
        if not accn:
            continue

        # Check cache
        if tip in assembly_cache and "ftp_path" in assembly_cache[tip]:
            results[tip] = assembly_cache[tip]
            continue
        to_resolve.append((tip, accn))

    print(f"  Cached: {len(results)}, need FTP lookup: {len(to_resolve)}", flush=True)

    for i, (tip, accn) in enumerate(to_resolve):
        try:
            time.sleep(REQUEST_DELAY)
            handle = Entrez.esearch(db="assembly", term=accn, retmax=1)
            search = Entrez.read(handle)
            handle.close()

            if not search["IdList"]:
                print(f"  WARN: No assembly for {accn} ({tip})", flush=True)
                results[tip] = {"accession": accn, "ftp_path": None, "source": "not_found"}
                assembly_cache[tip] = results[tip]
                continue

            time.sleep(REQUEST_DELAY)
            handle = Entrez.esummary(db="assembly", id=search["IdList"][0])
            summary = Entrez.read(handle)
            handle.close()

            doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
            ftp = doc.get("FtpPath_GenBank") or doc.get("FtpPath_RefSeq") or ""

            if ftp:
                results[tip] = {"accession": accn, "ftp_path": ftp, "source": "entrez"}
                assembly_cache[tip] = results[tip]
            else:
                # Try constructing FTP path from accession
                ftp = _construct_ftp_path(accn)
                results[tip] = {"accession": accn, "ftp_path": ftp, "source": "constructed"}
                assembly_cache[tip] = results[tip]

        except Exception as e:
            print(f"  ERROR resolving {accn}: {e}", flush=True)
            ftp = _construct_ftp_path(accn)
            results[tip] = {"accession": accn, "ftp_path": ftp, "source": "fallback"}
            assembly_cache[tip] = results[tip]

        if (i + 1) % 20 == 0:
            save_cache(ASSEMBLY_CACHE, assembly_cache)
            print(f"  Resolved {i + 1}/{len(to_resolve)}...", flush=True)

    save_cache(ASSEMBLY_CACHE, assembly_cache)
    valid = sum(1 for v in results.values() if v.get("ftp_path"))
    print(f"  FTP paths resolved: {valid}/{len(results)}", flush=True)
    return results


def _construct_ftp_path(accession: str) -> str:
    """Construct NCBI FTP path from accession (e.g., GCA_000006765.1)."""
    # NCBI FTP convention: ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/765/
    prefix = accession[:3]  # GCA or GCF
    digits = accession.split("_")[1].split(".")[0]  # 000006765
    # Pad to 9 digits
    digits = digits.zfill(9)
    d1, d2, d3 = digits[:3], digits[3:6], digits[6:9]
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{d1}/{d2}/{d3}/"
    return base


# ---------------------------------------------------------------------------
# [3/8] Download genomes
# ---------------------------------------------------------------------------
def download_genomes(
    ftp_info: dict[str, dict],
) -> dict[str, Path]:
    """Download genome FASTA files.

    Returns {tip_id: genome_dir_path}.
    """
    GENOMES_DIR.mkdir(parents=True, exist_ok=True)
    genome_dirs = {}
    n_downloaded = 0
    n_cached = 0
    n_failed = 0

    items = [(tip, info) for tip, info in ftp_info.items() if info.get("ftp_path")]

    for i, (tip, info) in enumerate(items):
        accn = info["accession"]
        genome_dir = GENOMES_DIR / accn
        genome_fna = genome_dir / "genome.fna.gz"
        protein_faa = genome_dir / "proteins.faa"

        # Check if already downloaded
        if genome_fna.exists() and genome_fna.stat().st_size > 10_000:
            genome_dirs[tip] = genome_dir
            n_cached += 1
            continue
        if protein_faa.exists() and protein_faa.stat().st_size > 100:
            # Already has proteins (from previous run with pre-downloaded proteins)
            genome_dirs[tip] = genome_dir
            n_cached += 1
            continue

        genome_dir.mkdir(parents=True, exist_ok=True)
        ftp_path = info["ftp_path"]

        try:
            # Build download URL from FTP path
            # Entrez returns: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/.../GCA_xxx_AsmName
            # Convert to HTTPS and append _genomic.fna.gz
            https_path = ftp_path.replace("ftp://", "https://")
            asm_name = https_path.rsplit("/", 1)[-1]
            url = f"{https_path}/{asm_name}_genomic.fna.gz"

            time.sleep(REQUEST_DELAY)
            urlretrieve(url, genome_fna)

            if genome_fna.stat().st_size < 1000:
                print(f"  WARN: {accn} genome too small ({genome_fna.stat().st_size} bytes)",
                      flush=True)
                genome_fna.unlink()
                n_failed += 1
                continue

            genome_dirs[tip] = genome_dir
            n_downloaded += 1

            # Also try to get pre-annotated proteins
            protein_url = url.replace("_genomic.fna.gz", "_protein.faa.gz")
            protein_gz = genome_dir / "protein.faa.gz"
            try:
                time.sleep(REQUEST_DELAY)
                urlretrieve(protein_url, protein_gz)
                if protein_gz.stat().st_size > 1000:
                    _decompress_gz(protein_gz, protein_faa)
                protein_gz.unlink(missing_ok=True)
            except Exception:
                protein_gz.unlink(missing_ok=True)

        except Exception as e:
            print(f"  FAIL download {accn}: {e}", flush=True)
            n_failed += 1

        if (i + 1) % 10 == 0:
            print(f"  Progress: {i + 1}/{len(items)} "
                  f"(downloaded={n_downloaded}, cached={n_cached}, failed={n_failed})",
                  flush=True)

    print(f"  Genome download: {n_downloaded} new, {n_cached} cached, {n_failed} failed",
          flush=True)
    print(f"  Total genomes available: {len(genome_dirs)}/{len(items)}", flush=True)

    if n_failed > 10:
        print(f"  WARNING: {n_failed} genomes failed to download!", flush=True)

    return genome_dirs


def _resolve_genome_url(accession: str, base_url: str) -> str:
    """Resolve the full URL for a genome FASTA given a base FTP path.

    For constructed paths, we need to discover the assembly name.
    Uses datasets API as fallback.
    """
    # Try NCBI datasets API to get the direct download info
    # First, try the simple convention: accession + version as directory name
    # E.g., GCA_000006765.2 → GCA_000006765.2_ASM676v2
    # Just use the Entrez-returned path directly since we already have it
    # For constructed paths, try fetching the directory listing

    # Simple approach: use NCBI datasets CLI or just try common patterns
    # Most complete genomes follow: accession_assemblyname
    asm_name = accession  # Start with just accession as guess
    url = f"{base_url}{asm_name}/{asm_name}_genomic.fna.gz"
    return url


def _decompress_gz(gz_path: Path, out_path: Path) -> None:
    """Decompress a .gz file."""
    with gzip.open(gz_path, "rb") as f_in:
        with open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


# ---------------------------------------------------------------------------
# [4/8] Predict proteins with Prodigal
# ---------------------------------------------------------------------------
def predict_proteins(genome_dirs: dict[str, Path]) -> dict[str, Path]:
    """Run Prodigal on genomes that don't have protein files yet.

    Returns {tip_id: protein_faa_path}.
    """
    protein_files = {}
    n_predicted = 0
    n_cached = 0
    n_failed = 0

    for tip, genome_dir in genome_dirs.items():
        protein_faa = genome_dir / "proteins.faa"

        # Check if proteins already exist
        if protein_faa.exists() and protein_faa.stat().st_size > 100:
            protein_files[tip] = protein_faa
            n_cached += 1
            continue

        # Need to decompress genome and run Prodigal
        genome_gz = genome_dir / "genome.fna.gz"
        genome_fna = genome_dir / "genome.fna"

        if not genome_gz.exists() and not genome_fna.exists():
            n_failed += 1
            continue

        try:
            # Decompress if needed
            if not genome_fna.exists():
                _decompress_gz(genome_gz, genome_fna)

            # Run Prodigal (standard mode for complete genomes)
            result = subprocess.run(
                [PRODIGAL, "-i", str(genome_fna), "-a", str(protein_faa), "-p", "single",
                 "-q"],
                capture_output=True, text=True, timeout=300,
            )
            if result.returncode != 0:
                # Try meta mode as fallback
                result = subprocess.run(
                    [PRODIGAL, "-i", str(genome_fna), "-a", str(protein_faa), "-p", "meta",
                     "-q"],
                    capture_output=True, text=True, timeout=300,
                )

            if protein_faa.exists() and protein_faa.stat().st_size > 100:
                protein_files[tip] = protein_faa
                n_predicted += 1
                # Count proteins
                n_proteins = sum(1 for line in open(protein_faa) if line.startswith(">"))
                if n_proteins < 500:
                    print(f"  WARN: {tip} has only {n_proteins} predicted proteins",
                          flush=True)
            else:
                n_failed += 1
                print(f"  FAIL Prodigal: {tip}", flush=True)

            # Clean up decompressed genome to save space
            if genome_gz.exists() and genome_fna.exists():
                genome_fna.unlink()

        except Exception as e:
            print(f"  ERROR Prodigal {tip}: {e}", flush=True)
            n_failed += 1

    print(f"  Prodigal: {n_predicted} new, {n_cached} cached, {n_failed} failed", flush=True)
    print(f"  Total protein files: {len(protein_files)}/{len(genome_dirs)}", flush=True)
    return protein_files


# ---------------------------------------------------------------------------
# [5/8] BLAST all references against predicted proteins
# ---------------------------------------------------------------------------
def blast_references(
    protein_files: dict[str, Path],
    ref_proteins: dict[str, str],
    cache: dict,
) -> dict:
    """BLAST each reference protein against each genome's predicted proteins.

    Updates and returns cache: {tip_id: {protein: identity, ...}}.
    """
    n_searched = 0
    n_cached = 0

    for tip, prot_file in protein_files.items():
        # Check if all proteins already searched
        if tip in cache and all(f"{p}_blastp" in cache[tip] for p in PROTEINS):
            n_cached += 1
            continue

        if tip not in cache:
            cache[tip] = {}

        genome_dir = prot_file.parent
        db_path = genome_dir / "blastdb"

        # Make BLAST database
        try:
            subprocess.run(
                [MAKEBLASTDB, "-in", str(prot_file), "-dbtype", "prot",
                 "-out", str(db_path)],
                capture_output=True, text=True, timeout=60,
            )
        except Exception as e:
            print(f"  ERROR makeblastdb {tip}: {e}", flush=True)
            continue

        # BLAST each reference
        for protein, seq in ref_proteins.items():
            cache_key = f"{protein}_blastp"
            if cache_key in cache[tip]:
                continue

            identity = _run_blastp(seq, db_path, EVALUE.get(protein, "1e-5"))
            cache[tip][cache_key] = identity

        n_searched += 1

        # Clean up BLAST DB files
        for ext in [".pdb", ".phr", ".pin", ".pjs", ".pot", ".psq", ".ptf", ".pto"]:
            (genome_dir / f"blastdb{ext}").unlink(missing_ok=True)

        if n_searched % 5 == 0:
            save_cache(CACHE_FILE, cache)

    save_cache(CACHE_FILE, cache)
    print(f"  BLAST: {n_searched} genomes searched, {n_cached} cached", flush=True)
    return cache


def _run_blastp(query_seq: str, db_path: Path, evalue: str) -> float | None:
    """Run blastp with a single query sequence. Returns best hit identity or None."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as f:
        f.write(f">query\n{query_seq}\n")
        query_file = f.name

    try:
        result = subprocess.run(
            [BLASTP, "-query", query_file, "-db", str(db_path),
             "-evalue", evalue, "-outfmt", "6 sseqid pident evalue bitscore",
             "-max_target_seqs", "1", "-num_threads", "4"],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode == 0 and result.stdout.strip():
            first_line = result.stdout.strip().split("\n")[0]
            fields = first_line.split("\t")
            return round(float(fields[1]), 1)
    except Exception:
        pass
    finally:
        os.unlink(query_file)

    return None


# ---------------------------------------------------------------------------
# [6/8] HMM search + tBLASTn for unannotated elements
# ---------------------------------------------------------------------------
def hmm_and_tblastn_search(
    protein_files: dict[str, Path],
    genome_dirs: dict[str, Path],
    ref_proteins: dict[str, str],
    cache: dict,
) -> dict:
    """Run hmmsearch (IS110 HMM) and tBLASTn (IS621/IS911 vs nucleotide).

    Updates cache with hmm_is110 and tblastn_IS621/tblastn_IS911 keys.
    """
    n_searched = 0
    n_cached = 0

    for tip, prot_file in protein_files.items():
        if tip not in cache:
            cache[tip] = {}

        # Check if HMM + tBLASTn already done
        done = all(k in cache[tip] for k in ["hmm_is110", "tblastn_IS621", "tblastn_IS911"])
        if done:
            n_cached += 1
            continue

        genome_dir = genome_dirs.get(tip)
        if not genome_dir:
            continue

        # --- HMM search for IS110 ---
        if "hmm_is110" not in cache[tip] and HMM_FILE.exists():
            identity = _run_hmmsearch_is110(prot_file, ref_proteins.get("IS621"))
            cache[tip]["hmm_is110"] = identity

        # --- tBLASTn IS621/IS911 against nucleotide ---
        genome_gz = genome_dir / "genome.fna.gz"
        genome_fna = genome_dir / "genome.fna"

        # Decompress if needed
        need_cleanup = False
        if not genome_fna.exists() and genome_gz.exists():
            _decompress_gz(genome_gz, genome_fna)
            need_cleanup = True

        if genome_fna.exists():
            for is_name in ["IS621", "IS911"]:
                tblastn_key = f"tblastn_{is_name}"
                if tblastn_key not in cache[tip] and is_name in ref_proteins:
                    identity = _run_tblastn(ref_proteins[is_name], genome_fna)
                    cache[tip][tblastn_key] = identity

            # Clean up decompressed genome
            if need_cleanup:
                genome_fna.unlink()

        n_searched += 1
        if n_searched % 5 == 0:
            save_cache(CACHE_FILE, cache)

    save_cache(CACHE_FILE, cache)
    print(f"  HMM+tBLASTn: {n_searched} genomes searched, {n_cached} cached", flush=True)
    return cache


def _run_hmmsearch_is110(protein_file: Path, is621_ref_seq: str | None) -> float | None:
    """Run hmmsearch with IS110 profiles against a protein file.

    If a hit is found, BLAST the hit against IS621 reference to get % identity.
    Returns best identity or None.
    """
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tblout", delete=False) as f:
            tblout = f.name

        result = subprocess.run(
            [HMMSEARCH, "--tblout", tblout, "-E", str(HMM_EVALUE),
             "--cpu", "4", str(HMM_FILE), str(protein_file)],
            capture_output=True, text=True, timeout=120,
        )

        if result.returncode != 0:
            return None

        # Parse tblout for hits
        hit_names = []
        with open(tblout) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    hit_names.append(parts[0])

        os.unlink(tblout)

        if not hit_names or not is621_ref_seq:
            return len(hit_names) > 0 or None  # Just report presence/absence

        # Extract hit sequences and BLAST against IS621 reference
        best_identity = None
        # Read protein sequences to find hits
        hit_seqs = {}
        for record in SeqIO.parse(str(protein_file), "fasta"):
            if record.id in hit_names:
                hit_seqs[record.id] = str(record.seq).rstrip("*")

        if not hit_seqs:
            return None

        # BLAST each hit against IS621 reference
        with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as f:
            f.write(f">IS621_ref\n{is621_ref_seq}\n")
            ref_file = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as f:
            for name, seq in hit_seqs.items():
                f.write(f">{name}\n{seq}\n")
            hits_file = f.name

        try:
            # Make DB from IS621 reference
            with tempfile.TemporaryDirectory() as tmpdir:
                db_path = os.path.join(tmpdir, "ref_db")
                subprocess.run(
                    [MAKEBLASTDB, "-in", ref_file, "-dbtype", "prot", "-out", db_path],
                    capture_output=True, timeout=30,
                )
                result = subprocess.run(
                    [BLASTP, "-query", hits_file, "-db", db_path,
                     "-evalue", "1e-3", "-outfmt", "6 qseqid pident evalue",
                     "-max_target_seqs", "1"],
                    capture_output=True, text=True, timeout=60,
                )
                if result.stdout.strip():
                    for line in result.stdout.strip().split("\n"):
                        fields = line.split("\t")
                        identity = float(fields[1])
                        if best_identity is None or identity > best_identity:
                            best_identity = identity
        except Exception:
            pass
        finally:
            os.unlink(ref_file)
            os.unlink(hits_file)

        return round(best_identity, 1) if best_identity else None

    except Exception as e:
        return None


def _run_tblastn(query_seq: str, genome_fna: Path) -> float | None:
    """Run tBLASTn query against a nucleotide genome. Returns best identity or None."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as f:
        f.write(f">query\n{query_seq}\n")
        query_file = f.name

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "nuc_db")
            subprocess.run(
                [MAKEBLASTDB, "-in", str(genome_fna), "-dbtype", "nucl", "-out", db_path],
                capture_output=True, timeout=60,
            )
            result = subprocess.run(
                [TBLASTN, "-query", query_file, "-db", db_path,
                 "-evalue", "1e-5", "-outfmt", "6 sseqid pident evalue bitscore",
                 "-max_target_seqs", "1", "-num_threads", "4"],
                capture_output=True, text=True, timeout=120,
            )
            if result.returncode == 0 and result.stdout.strip():
                first_line = result.stdout.strip().split("\n")[0]
                fields = first_line.split("\t")
                return round(float(fields[1]), 1)
    except Exception:
        pass
    finally:
        os.unlink(query_file)

    return None


# ---------------------------------------------------------------------------
# [7/8] Validation
# ---------------------------------------------------------------------------
def validate_results(
    taxonomy: list[dict],
    cache: dict,
) -> bool:
    """Validate coverage meets biological expectations.

    Returns True if all critical gates pass, False otherwise.
    """
    n_total = len(taxonomy)
    searched = [t for t in taxonomy if t["Tip_ID"] in cache]
    n_searched = len(searched)

    if n_searched == 0:
        print("  FATAL: No genomes were searched!", flush=True)
        return False

    print(f"\n  Validation ({n_searched} genomes searched):", flush=True)

    # For each protein, compute best-of-all-methods identity
    protein_found = defaultdict(int)
    for row in searched:
        tip = row["Tip_ID"]
        entry = cache.get(tip, {})
        for protein in PROTEINS:
            best = _best_identity(entry, protein)
            if best is not None:
                protein_found[protein] += 1

    ok = True
    for protein in PROTEINS:
        n = protein_found[protein]
        pct = 100 * n / n_searched if n_searched else 0
        status = ""
        if protein == "EF_Tu":
            if pct < 95:
                status = " *** FAIL: Expected >=95% ***"
                ok = False
            elif pct < 99:
                status = " (WARN: expected ~100%)"
        elif protein == "SecY":
            if pct < 90:
                status = " *** FAIL: Expected >=90% ***"
                ok = False
            elif pct < 95:
                status = " (WARN: expected >95%, some ultra-reduced genomes may lack)"
        elif protein == "IS621":
            if pct < 30:
                status = " (WARN: expected >=30% with HMM)"
        print(f"    {protein}: {n}/{n_searched} ({pct:.1f}%){status}", flush=True)

    # Domain breakdown
    bacteria = [t for t in searched if t.get("Domain") == "Bacteria"]
    archaea = [t for t in searched if t.get("Domain") == "Archaea"]
    for domain, subset in [("Bacteria", bacteria), ("Archaea", archaea)]:
        if not subset:
            continue
        print(f"\n  {domain} ({len(subset)} species):", flush=True)
        for protein in PROTEINS:
            n = sum(1 for t in subset
                    if _best_identity(cache.get(t["Tip_ID"], {}), protein) is not None)
            pct = 100 * n / len(subset)
            print(f"    {protein}: {n}/{len(subset)} ({pct:.1f}%)", flush=True)

    return ok


def _best_identity(entry: dict, protein: str) -> float | None:
    """Get best identity for a protein across all search methods."""
    candidates = []

    # blastp
    v = entry.get(f"{protein}_blastp")
    if v is not None:
        candidates.append(v)

    # HMM (only for IS621 → hmm_is110)
    if protein == "IS621":
        v = entry.get("hmm_is110")
        if isinstance(v, (int, float)):
            candidates.append(v)

    # tBLASTn
    v = entry.get(f"tblastn_{protein}")
    if isinstance(v, (int, float)):
        candidates.append(v)

    return max(candidates) if candidates else None


# ---------------------------------------------------------------------------
# [8/8] Update metadata.csv
# ---------------------------------------------------------------------------
def update_metadata(
    taxonomy: list[dict],
    cache: dict,
) -> None:
    """Overwrite metadata.csv with local search results."""
    # Load existing metadata for non-identity fields
    existing = {}
    if METADATA_CSV.exists():
        with open(METADATA_CSV) as f:
            for row in csv.DictReader(f):
                existing[row["Tip_ID"]] = row

    rows = []
    for tax_row in taxonomy:
        tip = tax_row["Tip_ID"]
        entry = cache.get(tip, {})
        old = existing.get(tip, {})

        row = {
            "Tip_ID": tip,
            "Phylum": old.get("Phylum", tax_row.get("Phylum", "")),
            "Ecoli_Reference": old.get("Ecoli_Reference", "FALSE"),
        }

        for protein in PROTEINS:
            col = f"{protein}_Identity"
            best = _best_identity(entry, protein)
            if best is not None:
                row[col] = best
            elif row["Ecoli_Reference"] == "TRUE":
                row[col] = 100.0
            else:
                row[col] = "NA"

        rows.append(row)

    # Write
    fields = ["Tip_ID", "Phylum", "Ecoli_Reference",
              "EF_Tu_Identity", "SecY_Identity", "MuA_Identity",
              "IS621_Identity", "IS911_Identity"]
    with open(METADATA_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    # Print before/after comparison
    print(f"\n  Wrote {len(rows)} rows to {METADATA_CSV}", flush=True)
    print(f"\n  Coverage summary:", flush=True)
    for protein in PROTEINS:
        col = f"{protein}_Identity"
        n_hit = sum(1 for r in rows if r[col] != "NA")
        print(f"    {protein}: {n_hit}/{len(rows)}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60, flush=True)
    print("Local genome search (V3)", flush=True)
    print("=" * 60, flush=True)

    # [1/8] Load data
    print("\n[1/8] Loading data...", flush=True)
    taxonomy, tip_to_accn, ref_proteins = load_data()

    # [2/8] Resolve FTP paths
    print("\n[2/8] Resolving FTP paths...", flush=True)
    ftp_info = resolve_ftp_paths(taxonomy, tip_to_accn)

    # [3/8] Download genomes
    print("\n[3/8] Downloading genomes...", flush=True)
    genome_dirs = download_genomes(ftp_info)

    # [4/8] Predict proteins
    print("\n[4/8] Predicting proteins with Prodigal...", flush=True)
    protein_files = predict_proteins(genome_dirs)

    # Load search cache
    cache = load_cache(CACHE_FILE)

    # [5/8] BLAST references
    print("\n[5/8] BLAST reference proteins...", flush=True)
    cache = blast_references(protein_files, ref_proteins, cache)

    # [6/8] HMM + tBLASTn
    print("\n[6/8] HMM search + tBLASTn...", flush=True)
    cache = hmm_and_tblastn_search(protein_files, genome_dirs, ref_proteins, cache)

    # [7/8] Validate
    print("\n[7/8] Validation...", flush=True)
    passed = validate_results(taxonomy, cache)

    # [8/8] Update metadata
    print("\n[8/8] Updating metadata...", flush=True)
    update_metadata(taxonomy, cache)

    if passed:
        print("\n=== ALL VALIDATION GATES PASSED ===", flush=True)
    else:
        print("\n=== WARNING: SOME VALIDATION GATES FAILED ===", flush=True)
        print("Review the output above and investigate.", flush=True)

    print("\nDone! Run the R script to generate the tree:", flush=True)
    print("  Rscript scripts/04_plot_tree.R", flush=True)


if __name__ == "__main__":
    main()
