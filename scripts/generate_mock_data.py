#!/usr/bin/env python3
"""Generate biologically realistic mock data for testing the visualization.

Outputs: data/metadata.csv  and  data/taxonomy.csv
in the same format as the real pipeline (03_build_metadata.py).

Mock species span ~15 major prokaryotic phyla with realistic identity patterns:
  - EF-Tu: universal (45-99%), higher in Pseudomonadota
  - SecY: universal (30-85%), slightly less conserved
  - MuA: patchy (~30% of species), 25-70%
  - IS621: very sparse (~15%), concentrated in Pseudomonadota
"""

import csv
import random
import re
from pathlib import Path

OUTPUT_METADATA = Path("data/metadata.csv")
OUTPUT_TAXONOMY = Path("data/taxonomy.csv")

random.seed(42)

# Phyla with realistic taxonomic structure
# Format: {phylum: [(class, order, family, genus, species_name), ...]}
TAXONOMY = {
    # --- Bacteria ---
    "Pseudomonadota": {  # formerly Proteobacteria
        "Gammaproteobacteria": [
            ("Enterobacterales", "Enterobacteriaceae", "Escherichia", "Escherichia_coli_K12"),
            ("Enterobacterales", "Enterobacteriaceae", "Salmonella", "Salmonella_enterica"),
            ("Enterobacterales", "Enterobacteriaceae", "Klebsiella", "Klebsiella_pneumoniae"),
            ("Enterobacterales", "Yersiniaceae", "Yersinia", "Yersinia_pestis"),
            ("Pseudomonadales", "Pseudomonadaceae", "Pseudomonas", "Pseudomonas_aeruginosa"),
            ("Pseudomonadales", "Moraxellaceae", "Acinetobacter", "Acinetobacter_baumannii"),
            ("Vibrionales", "Vibrionaceae", "Vibrio", "Vibrio_cholerae"),
            ("Xanthomonadales", "Xanthomonadaceae", "Xanthomonas", "Xanthomonas_campestris"),
            ("Legionellales", "Legionellaceae", "Legionella", "Legionella_pneumophila"),
            ("Pasteurellales", "Pasteurellaceae", "Haemophilus", "Haemophilus_influenzae"),
        ],
        "Betaproteobacteria": [
            ("Burkholderiales", "Burkholderiaceae", "Burkholderia", "Burkholderia_cenocepacia"),
            ("Burkholderiales", "Comamonadaceae", "Comamonas", "Comamonas_testosteroni"),
            ("Neisseriales", "Neisseriaceae", "Neisseria", "Neisseria_meningitidis"),
            ("Nitrosomonadales", "Nitrosomonadaceae", "Nitrosomonas", "Nitrosomonas_europaea"),
        ],
        "Alphaproteobacteria": [
            ("Rhizobiales", "Rhizobiaceae", "Rhizobium", "Rhizobium_leguminosarum"),
            ("Rhizobiales", "Bradyrhizobiaceae", "Bradyrhizobium", "Bradyrhizobium_japonicum"),
            ("Caulobacterales", "Caulobacteraceae", "Caulobacter", "Caulobacter_vibrioides"),
            ("Sphingomonadales", "Sphingomonadaceae", "Sphingomonas", "Sphingomonas_paucimobilis"),
            ("Rickettsiales", "Rickettsiaceae", "Rickettsia", "Rickettsia_prowazekii"),
        ],
        "Deltaproteobacteria": [
            ("Myxococcales", "Myxococcaceae", "Myxococcus", "Myxococcus_xanthus"),
            ("Bdellovibrionales", "Bdellovibrionaceae", "Bdellovibrio", "Bdellovibrio_bacteriovorus"),
        ],
        "Epsilonproteobacteria": [
            ("Campylobacterales", "Campylobacteraceae", "Campylobacter", "Campylobacter_jejuni"),
            ("Campylobacterales", "Helicobacteraceae", "Helicobacter", "Helicobacter_pylori"),
        ],
    },
    "Bacillota": {  # formerly Firmicutes
        "Bacilli": [
            ("Bacillales", "Bacillaceae", "Bacillus", "Bacillus_subtilis"),
            ("Bacillales", "Staphylococcaceae", "Staphylococcus", "Staphylococcus_aureus"),
            ("Bacillales", "Listeriaceae", "Listeria", "Listeria_monocytogenes"),
            ("Lactobacillales", "Lactobacillaceae", "Lactobacillus", "Lactobacillus_acidophilus"),
            ("Lactobacillales", "Streptococcaceae", "Streptococcus", "Streptococcus_pneumoniae"),
            ("Lactobacillales", "Enterococcaceae", "Enterococcus", "Enterococcus_faecalis"),
        ],
        "Clostridia": [
            ("Eubacteriales", "Clostridiaceae", "Clostridium", "Clostridium_difficile"),
            ("Eubacteriales", "Lachnospiraceae", "Roseburia", "Roseburia_intestinalis"),
            ("Eubacteriales", "Peptostreptococcaceae", "Peptoclostridium", "Clostridioides_difficile"),
        ],
    },
    "Actinomycetota": {  # formerly Actinobacteria
        "Actinomycetes": [
            ("Mycobacteriales", "Mycobacteriaceae", "Mycobacterium", "Mycobacterium_tuberculosis"),
            ("Mycobacteriales", "Corynebacteriaceae", "Corynebacterium", "Corynebacterium_glutamicum"),
            ("Streptomycetales", "Streptomycetaceae", "Streptomyces", "Streptomyces_coelicolor"),
            ("Propionibacteriales", "Propionibacteriaceae", "Cutibacterium", "Cutibacterium_acnes"),
            ("Micrococcales", "Micrococcaceae", "Micrococcus", "Micrococcus_luteus"),
            ("Bifidobacteriales", "Bifidobacteriaceae", "Bifidobacterium", "Bifidobacterium_longum"),
        ],
    },
    "Bacteroidota": {  # formerly Bacteroidetes
        "Bacteroidia": [
            ("Bacteroidales", "Bacteroidaceae", "Bacteroides", "Bacteroides_fragilis"),
            ("Bacteroidales", "Prevotellaceae", "Prevotella", "Prevotella_copri"),
            ("Bacteroidales", "Porphyromonadaceae", "Porphyromonas", "Porphyromonas_gingivalis"),
            ("Flavobacteriales", "Flavobacteriaceae", "Flavobacterium", "Flavobacterium_johnsoniae"),
            ("Cytophagales", "Cytophagaceae", "Cytophaga", "Cytophaga_hutchinsonii"),
        ],
    },
    "Cyanobacteriota": {
        "Cyanophyceae": [
            ("Synechococcales", "Synechococcaceae", "Synechococcus", "Synechococcus_elongatus"),
            ("Nostocales", "Nostocaceae", "Nostoc", "Nostoc_punctiforme"),
            ("Oscillatoriales", "Oscillatoriaceae", "Oscillatoria", "Oscillatoria_nigroviridis"),
            ("Chroococcales", "Microcystaceae", "Microcystis", "Microcystis_aeruginosa"),
        ],
    },
    "Spirochaetota": {
        "Spirochaetia": [
            ("Spirochaetales", "Spirochaetaceae", "Treponema", "Treponema_pallidum"),
            ("Spirochaetales", "Borreliaceae", "Borrelia", "Borrelia_burgdorferi"),
            ("Leptospirales", "Leptospiraceae", "Leptospira", "Leptospira_interrogans"),
        ],
    },
    "Thermodesulfobacteriota": {
        "Thermodesulfobacteria": [
            ("Thermodesulfobacteriales", "Thermodesulfobacteriaceae", "Thermodesulfobacterium", "Thermodesulfobacterium_commune"),
        ],
    },
    "Deinococcota": {
        "Deinococci": [
            ("Deinococcales", "Deinococcaceae", "Deinococcus", "Deinococcus_radiodurans"),
            ("Thermales", "Thermaceae", "Thermus", "Thermus_thermophilus"),
        ],
    },
    "Chloroflexota": {
        "Chloroflexia": [
            ("Chloroflexales", "Chloroflexaceae", "Chloroflexus", "Chloroflexus_aurantiacus"),
            ("Herpetosiphonales", "Herpetosiphonaceae", "Herpetosiphon", "Herpetosiphon_aurantiacus"),
        ],
    },
    "Planctomycetota": {
        "Planctomycetes": [
            ("Planctomycetales", "Planctomycetaceae", "Planctomyces", "Planctomyces_limnophilus"),
            ("Pirellulales", "Pirellulaceae", "Pirellula", "Pirellula_staleyi"),
        ],
    },
    "Verrucomicrobiota": {
        "Verrucomicrobiae": [
            ("Verrucomicrobiales", "Verrucomicrobiaceae", "Verrucomicrobium", "Verrucomicrobium_spinosum"),
            ("Chthoniobacterales", "Chthoniobacteraceae", "Chthoniobacter", "Chthoniobacter_flavus"),
        ],
    },
    "Fusobacteriota": {
        "Fusobacteriia": [
            ("Fusobacteriales", "Fusobacteriaceae", "Fusobacterium", "Fusobacterium_nucleatum"),
        ],
    },
    "Chlamydiota": {
        "Chlamydiia": [
            ("Chlamydiales", "Chlamydiaceae", "Chlamydia", "Chlamydia_trachomatis"),
        ],
    },
    # --- Archaea ---
    "Euryarchaeota": {
        "Methanobacteria": [
            ("Methanobacteriales", "Methanobacteriaceae", "Methanobacterium", "Methanobacterium_thermoautotrophicum"),
        ],
        "Methanococci": [
            ("Methanococcales", "Methanococcaceae", "Methanococcus", "Methanococcus_jannaschii"),
        ],
        "Halobacteria": [
            ("Halobacteriales", "Halobacteriaceae", "Halobacterium", "Halobacterium_salinarum"),
            ("Haloferacales", "Haloferacaceae", "Haloferax", "Haloferax_volcanii"),
        ],
        "Thermococci": [
            ("Thermococcales", "Thermococcaceae", "Thermococcus", "Thermococcus_kodakarensis"),
            ("Thermococcales", "Thermococcaceae", "Pyrococcus", "Pyrococcus_furiosus"),
        ],
    },
    "Thermoproteota": {  # formerly Crenarchaeota
        "Thermoprotei": [
            ("Sulfolobales", "Sulfolobaceae", "Sulfolobus", "Sulfolobus_acidocaldarius"),
            ("Desulfurococcales", "Desulfurococcaceae", "Desulfurococcus", "Desulfurococcus_mucosus"),
            ("Thermoproteales", "Thermoproteaceae", "Thermoproteus", "Thermoproteus_tenax"),
        ],
    },
    "Asgardarchaeota": {
        "Lokiarchaeia": [
            ("Lokiarchaeales", "Lokiarchaeaceae", "Lokiarchaeum", "Lokiarchaeum_sp_GC14_75"),
        ],
        "Thorarchaeia": [
            ("Thorarchaeales", "Thorarchaeaceae", "Thorarchaeum", "Thorarchaeum_sp_SMTZ_45"),
        ],
    },
}

# Phylogenetic distance from E. coli (rough estimate, 0=same, 1=most distant)
PHYLUM_DISTANCE = {
    "Pseudomonadota": 0.0,
    "Bacillota": 0.4,
    "Actinomycetota": 0.45,
    "Bacteroidota": 0.5,
    "Cyanobacteriota": 0.55,
    "Spirochaetota": 0.5,
    "Thermodesulfobacteriota": 0.6,
    "Deinococcota": 0.55,
    "Chloroflexota": 0.6,
    "Planctomycetota": 0.55,
    "Verrucomicrobiota": 0.5,
    "Fusobacteriota": 0.55,
    "Chlamydiota": 0.5,
    "Euryarchaeota": 0.85,
    "Thermoproteota": 0.9,
    "Asgardarchaeota": 0.8,
}

# Class-level distance adjustments within Pseudomonadota
CLASS_DISTANCE = {
    "Gammaproteobacteria": 0.0,
    "Betaproteobacteria": 0.15,
    "Alphaproteobacteria": 0.2,
    "Deltaproteobacteria": 0.3,
    "Epsilonproteobacteria": 0.25,
}


def generate_identity(base_high: float, base_low: float, distance: float,
                      noise: float = 5.0) -> float:
    """Generate a realistic identity value based on phylogenetic distance."""
    identity = base_high - (base_high - base_low) * distance
    identity += random.gauss(0, noise)
    return max(base_low * 0.8, min(100, round(identity, 1)))


def main():
    taxonomy_rows = []
    metadata_rows = []

    for phylum, classes in TAXONOMY.items():
        domain = "Archaea" if phylum in ("Euryarchaeota", "Thermoproteota", "Asgardarchaeota") else "Bacteria"
        base_dist = PHYLUM_DISTANCE.get(phylum, 0.5)

        for cls, species_list in classes.items():
            for order, family, genus, species in species_list:
                # Calculate phylogenetic distance
                dist = base_dist
                if phylum == "Pseudomonadota":
                    dist += CLASS_DISTANCE.get(cls, 0.2)

                tip_id = species

                # Taxonomy row
                taxonomy_rows.append({
                    "Tip_ID": tip_id,
                    "Domain": domain,
                    "Phylum": phylum,
                    "Class": cls,
                    "Order": order,
                    "Family": family,
                    "Genus": genus,
                    "Species": tip_id,
                })

                # Metadata row
                is_ecoli = species == "Escherichia_coli_K12"

                # EF-Tu: universal, 99% for E.coli, ~50% for distant archaea
                eftu = 100.0 if is_ecoli else generate_identity(99, 50, dist)

                # SecY: universal, slightly less conserved
                secy = 100.0 if is_ecoli else generate_identity(90, 30, dist, noise=6)

                # MuA: patchy, mainly in Proteobacteria and some Firmicutes
                mua = "NA"
                if is_ecoli:
                    mua = 100.0
                elif phylum == "Pseudomonadota" and random.random() < 0.6:
                    mua = generate_identity(70, 30, dist * 0.8)
                elif phylum == "Bacillota" and random.random() < 0.25:
                    mua = generate_identity(50, 25, 0.5)
                elif phylum in ("Actinomycetota", "Bacteroidota") and random.random() < 0.1:
                    mua = generate_identity(40, 25, 0.6)

                # IS621: very sparse, mainly in Enterobacteriaceae and some other Gamma
                is621 = "NA"
                if is_ecoli:
                    is621 = 100.0
                elif phylum == "Pseudomonadota" and cls == "Gammaproteobacteria":
                    if family == "Enterobacteriaceae":
                        is621 = generate_identity(95, 50, dist * 0.5)
                    elif random.random() < 0.3:
                        is621 = generate_identity(60, 25, dist * 1.2)
                elif phylum == "Pseudomonadota" and random.random() < 0.08:
                    is621 = generate_identity(40, 20, 0.6)
                elif random.random() < 0.03:
                    is621 = generate_identity(30, 20, 0.7)

                metadata_rows.append({
                    "Tip_ID": tip_id,
                    "Phylum": phylum,
                    "Ecoli_Reference": "TRUE" if is_ecoli else "FALSE",
                    "EF_Tu_Identity": eftu,
                    "SecY_Identity": secy,
                    "MuA_Identity": mua,
                    "IS621_Identity": is621,
                })

    # Write taxonomy CSV
    OUTPUT_TAXONOMY.parent.mkdir(parents=True, exist_ok=True)
    tax_fields = ["Tip_ID", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    with open(OUTPUT_TAXONOMY, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tax_fields)
        writer.writeheader()
        writer.writerows(taxonomy_rows)

    # Write metadata CSV
    meta_fields = ["Tip_ID", "Phylum", "Ecoli_Reference",
                   "EF_Tu_Identity", "SecY_Identity", "MuA_Identity", "IS621_Identity"]
    with open(OUTPUT_METADATA, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=meta_fields)
        writer.writeheader()
        writer.writerows(metadata_rows)

    print(f"Generated mock data:")
    print(f"  {len(taxonomy_rows)} species across {len(TAXONOMY)} phyla")
    print(f"  Taxonomy: {OUTPUT_TAXONOMY}")
    print(f"  Metadata: {OUTPUT_METADATA}")

    # Summary stats
    n_mua = sum(1 for r in metadata_rows if r["MuA_Identity"] != "NA")
    n_is621 = sum(1 for r in metadata_rows if r["IS621_Identity"] != "NA")
    print(f"  MuA hits: {n_mua}/{len(metadata_rows)} ({100*n_mua/len(metadata_rows):.0f}%)")
    print(f"  IS621 hits: {n_is621}/{len(metadata_rows)} ({100*n_is621/len(metadata_rows):.0f}%)")


if __name__ == "__main__":
    main()
