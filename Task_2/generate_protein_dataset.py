"""
generate_protein_dataset.py

Reads a CSV or Excel file of UniProt accessions from a specified column
(default: "Master Protein Accessions"), deduplicates them, shows a summary,
and then lets the user choose between:
  1) Generating a dataset for ALL unique accessions
  2) Generating a dataset for the FIRST N unique accessions (user-specified N)
  3) Abort without generating a dataset

Uses the EBI Proteins API (https://www.ebi.ac.uk/proteins/api/) to fetch
JSON-based information for each accession, then parses fields based on that
schema to extract:
  - proteinName
  - accessions (primary + secondary)
  - gene
  - status (Reviewed/Unreviewed)
  - organism
  - variantIDs
  - goAnnotations
  - associatedDiseases
  - familiesAndDomains
  - sequence
  - isoforms
  - regions

Usage:
    python generate_protein_dataset.py \
        --input accessions.csv \
        --column "Master Protein Accessions" \
        --output proteins.json \
        [--delay 0.1]
"""

import argparse
import json
import time
import requests
import pandas as pd


REQUEST_URL = "https://www.ebi.ac.uk/proteins/api/proteins/"


def fetch_uniprot_entry(accession: str, sleep_between: float = 0.1) -> dict:
    """
    Fetches a UniProt entry via EBI Proteins API in JSON format for the given
    accession. Returns the parsed JSON as a Python dict; if not found or
    error, returns an empty dictionary ({}).
    """
    url = REQUEST_URL + accession
    try:
        resp = requests.get(url, headers={"Accept": "application/json"},
                            timeout=30)
        if resp.status_code == 200:
            return resp.json()
        else:
            print(f"Warning: EBI API returned status {resp.status_code} for "
                  f" {accession}")
            return {}
    except Exception as e:
        print(f"Error fetching {accession}: {e}")
        return {}
    finally:
        time.sleep(sleep_between)


def _parse_properties(raw) -> dict:
    """
    Given a 'properties' field from EBI JSON, normalizes it into a dictionary
    It might be:
     - a list of dictionaries: [{"key": "...", "value": "..."}, ...]
     - a single dictionary (rare)
     - a semicolon-separated string "key1=value1;key2=value2"
    """
    props = {}
    if isinstance(raw, list):
        for item in raw:
            if isinstance(item, dict):
                k = item.get("key")
                v = item.get("value")
                if k and v is not None:
                    props[k] = v
    elif isinstance(raw, dict):
        for k, v in raw.items():
            props[k] = v
    elif isinstance(raw, str):
        for segment in raw.split(";"):
            if "=" in segment:
                k, v = segment.split("=", 1)
                props[k.strip()] = v.strip()
    return props


def parse_protein_name(entry_json: dict) -> str:
    """
    Extracts the recommended full name from entry_json.ebi, falling back to
    first available. EBI returns names under entry_json["protein"][
    "recommendedName"]["fullName"] or similar.
    """
    prot_block = entry_json.get("protein", {})
    rec = prot_block.get("recommendedName", {})
    fn = rec.get("fullName")
    if fn:
        return fn
    # If no recommendedName, check alternativeNames list
    for alt in prot_block.get("alternativeNames", []):
        fn2 = alt.get("fullName")
        if fn2:
            return fn2
    return ""


def parse_existence_info(entry_json: dict) -> dict:
    """
    Extracts the protein existence level and general entry information from:
      {
        "proteinExistence": <string>,     # e.g. "Evidence at protein level"
        "sourceDatabase": <string>,       # e.g. "Swiss-Prot"
        "created": <YYYY-MM-DD string>,
        "modified": <YYYY-MM-DD string>,
        "version": <int>
      }

    EBI JSON example:
      "proteinExistence": "Evidence at protein level",
      "info": {
        "type": "Swiss-Prot",
        "created": "2007-10-02",
        "modified": "2025-04-09",
        "version": 125
      }
    """
    pe = entry_json.get("proteinExistence", "")
    info = entry_json.get("info", {})

    return {
        "proteinExistence": pe,
        "sourceDatabase":   info.get("type", ""),
        "created":          info.get("created", ""),
        "modified":         info.get("modified", ""),
        "version":          info.get("version")
    }


def parse_accessions(entry_json: dict) -> list:
    """
    Returns a list of dictionaries: { "id": accession, "type":
    "Primary"/"Secondary" }.
    EBI returns .accession (primary) and .secondaryAccessions (list).
    """
    accs = []
    primary = entry_json.get("accession")
    if primary:
        accs.append({"id": primary, "type": "Primary"})
    for sec in entry_json.get("secondaryAccession", []):
        accs.append({"id": sec, "type": "Secondary"})
    return accs


def parse_gene(entry_json: dict) -> list:
    """
    Extracts the first gene name. EBI has entry_json["genes"], each with
    "geneName" object.
    """
    genes = []
    for g in entry_json.get("gene", []):
        name = None
        if isinstance(g.get("name"), dict):
            name = g["name"].get("value")

        syns = []
        for s in g.get("synonyms", []):
            val = s.get("value")
            if val:
                syns.append(val)

        if name:
            genes.append({
                "name": name,
                "synonyms": syns})

    return genes


def parse_status(entry_json: dict) -> list:
    """
    Extracts all top‐level 'comments' from the EBI JSON and normalises them into
    a list of dictionaries with the right sub‐fields for each comment type.

    Supported types:
      - FUNCTION             → extracts text.value and evidences
      - INTERACTION          → passes through the 'interactions' list
      - SUBCELLULAR_LOCATION → extracts locations and topology blocks
      - ALTERNATIVE_PRODUCTS → extracts event and isoforms
      - TISSUE_SPECIFICITY   → extracts text.value and evidences
      - DISEASE              → extracts text.value
      - SIMILARITY           → extracts text.value and evidences
      - SEQUENCE_CAUTION     → extracts conflictType, sequence, text (if
      present), evidences

    Returns a list like:
      [
        { "type": "FUNCTION", "text": [ { "value":…, "evidences":[…] }, … ] },
        { "type": "INTERACTION", "interactions": [ … ] },
        …
      ]
    """
    comments_out = []

    for com in entry_json.get("comments", []):
        ctype = com.get("type")
        out = {"type": ctype}

        if ctype == "FUNCTION":
            out["text"] = com.get("text", [])

        elif ctype == "INTERACTION":
            out["interactions"] = com.get("interactions", [])

        elif ctype == "SUBCELLULAR_LOCATION":
            out["locations"] = com.get("locations", [])

        elif ctype == "ALTERNATIVE_PRODUCTS":
            # event is a list of strings, isoforms a list of dicts
            out["event"]    = com.get("event", [])
            out["isoforms"] = com.get("isoforms", [])

        elif ctype == "TISSUE_SPECIFICITY":
            out["text"] = com.get("text", [])

        elif ctype == "DISEASE":
            out["text"] = com.get("text", [])

        elif ctype == "SIMILARITY":
            out["text"] = com.get("text", [])

        elif ctype == "SEQUENCE_CAUTION":
            # There may be multiple sequence‐caution entries, each with:
            #   conflictType, sequence, text (optional), evidences
            out["conflictType"] = com.get("conflictType")
            out["sequence"]     = com.get("sequence")
            if "text" in com:
                out["text"] = com.get("text")
            out["evidences"] = com.get("evidences", [])

        else:
            # any other comment types you didn’t list—just dump raw
            out.update({k: com[k] for k in com if k not in ("type",)})

        comments_out.append(out)

    return comments_out


def parse_organism(entry_json: dict) -> dict:
    """
    Returns scientific name from entry_json["organism"]["scientificName"].
    """
    org = entry_json.get("organism", {})
    taxonomy = org.get("taxonomy")

    scientific = ""
    common = ""
    for name_obj in org.get("names", []):
        t = name_obj.get("type")
        v = name_obj.get("value")
        if t == "scientific":
            scientific = v
        elif t == "common":
            common = v

    return {
        "taxonomy": taxonomy,
        "scientific": scientific,
        "common": common
    }

def parse_variants(entry_json: dict) -> list:
    """
    Extracts variant IDs from entry_json["features"] where type == "VARIANT".
    """
    variants = []
    for feat in entry_json.get("features", []):
        if feat.get("category") != "VARIANTS":
            continue

        # common fields
        ctype = feat.get("type")
        ftId = feat.get("ftId", "")
        try:
            begin = int(feat.get("begin"))
            end = int(feat.get("end"))
        except (TypeError, ValueError):
            # skip if no numeric range
            continue

        # collect evidences
        evs = []
        for e in feat.get("evidences", []):
            evs.append({
                "code": e.get("code"),
                "source": e.get("source", {}).get("name"),
                "id": e.get("source", {}).get("id")
            })

        base = {
            "type": ctype,
            "ftId": ftId,
            "begin": begin,
            "end": end,
        }
        if evs:
            base["evidences"] = evs

        if ctype == "VAR_SEQ":
            # nothing extra
            pass
        elif ctype == "VARIANT":
            # add the variant‐specific fields
            base["description"] = feat.get("description", "")
            base["alternativeSequence"] = feat.get("alternativeSequence", "")

        variants.append(base)

    return variants


def parse_go_annotations(entry_json: dict) -> list:
    """
    Extracts GO annotations from entry_json["dbReferences"] where type == "GO".
    Each dbReference may have 'properties' describing term/category,
    these are also extracted.
    """
    go_list = []
    for db_ref in entry_json.get("dbReferences", []):
        if db_ref.get("type") == "GO":
            props = _parse_properties(db_ref.get("properties", []))
            go_list.append({
                "id": db_ref.get("id"),
                "term": props.get("term"),
                "source": props.get("source")
            })
    return go_list


def parse_associated_diseases(entry_json: dict) -> list:
    """
    Extracts associated diseases from entry_json["comments"] where category
    == "DISEASE".
    """
    diseases = []
    for com in entry_json.get("comments", []):
        if com.get("type") == "DISEASE":
            # Pull all text blocks
            texts = [txt.get("value") for txt in com.get("text", []) if
                     txt.get("value")]

            # Any evidences attached to each text
            evs = []
            for txt in com.get("text", []):
                for e in txt.get("evidences", []):
                    evs.append({
                        "code": e.get("code"),
                        "source": e.get("source", {}).get("name"),
                        "id": e.get("source", {}).get("id")
                    })

            # A comment‐level "scope" array (if present)
            scope = com.get("scope", [])

            diseases.append({
                "texts": texts,
                "evidences": evs,
                "scope": scope
            })
    return diseases


def parse_families_and_domains(entry_json: dict) -> list:
    """
    Extracts family/domain info from entry_json["dbReferences"] for relevant
    sources.
    """
    fd_list = []
    relevant = {"Pfam", "SMART", "PROSITE", "InterPro", "SUPFAM"}
    for db_ref in entry_json.get("dbReferences", []):
        dbtype = db_ref.get("type")
        if dbtype in relevant:
            props = _parse_properties(db_ref.get("properties", []))
            fd_list.append({
                "id": db_ref.get("id"),
                "type": dbtype,
                "description": props.get("entryName") or props.get("name")
            })
    return fd_list


def parse_sequence(entry_json: dict) -> str:
    """
    Returns canonical sequences from entry_json["sequence"]["sequence"].
    EBI uses nested .sequence object with "sequence" key.
    """
    return entry_json.get("sequence", {}).get("sequence", "")


def parse_isoforms(entry_json: dict) -> list:
    """
    Extracts isoform info from entry_json["comments"] where category ==
    "ALTERNATIVE PRODUCTS".
    """
    iso_list = []
    for com in entry_json.get("comments", []):
        if com.get("type") == "ALTERNATIVE_PRODUCTS":
            for iso in com.get("isoforms", []):
                # ids is already a list
                ids = iso.get("ids", [])

                # name.value
                name = iso.get("name", {}).get("value", "")

                # status
                status = iso.get("sequenceStatus", "")

                # a possible 'sequence' array of ftId strings
                seq = iso.get("sequence", [])

                iso_list.append({
                    "ids": ids,
                    "name": name,
                    "sequenceStatus": status,
                    "sequence": seq
                })
    return iso_list


def parse_regions(entry_json: dict) -> list:
    """
    Collects features with a start/end range from entry_json["features"].
    """
    regions = []
    for feat in entry_json.get("features", []):
        if feat.get("category") == "VARIANTS":
            # skip, parsed by parse_variants()
            continue

        try:
            begin = int(feat.get("begin"))
            end = int(feat.get("end"))
        except (TypeError, ValueError):
            continue

        region = {
            "type": feat.get("type"),
            "category": feat.get("category"),
            "description": feat.get("description", ""),
            "begin": begin,
            "end": end
        }

        evs = []
        for e in feat.get("evidences", []):
            evs.append({
                "code": e.get("code"),
                "source": e.get("source", {}).get("name"),
                "id": e.get("source", {}).get("id")
            })
        if evs:
            region["evidences"] = evs

        regions.append(region)

    return regions


def build_protein_record(accession: str) -> dict:
    """
    Fetches EBI JSON for 'accession' and builds a record with:
      {
        queryAccession,
        proteinName,
        accessions,
        gene,
        status,
        organism,
        variantIDs,
        goAnnotations,
        associatedDiseases,
        familiesAndDomains,
        sequence,
        isoforms,
        regions
      }
    """
    entry = fetch_uniprot_entry(accession)
    if not entry:
        return {"queryAccession": accession, "error": "No data returned or  "
                                                      "accession not found."}

    return {
        "queryAccession": accession,
        "proteinExistence": parse_existence_info(entry),
        "proteinName": parse_protein_name(entry),
        "accessions": parse_accessions(entry),
        "gene": parse_gene(entry),
        "status": parse_status(entry),
        "organism": parse_organism(entry),
        "variantIDs": parse_variants(entry),
        "goAnnotations": parse_go_annotations(entry),
        "associatedDiseases": parse_associated_diseases(entry),
        "familiesAndDomains": parse_families_and_domains(entry),
        "sequence": parse_sequence(entry),
        "isoforms": parse_isoforms(entry),
        "regions": parse_regions(entry),
    }


def read_accessions_any_format(path: str, column_name: str) -> list:
    """
    Attempts to read `path` as CSV (utf-8 then latin1). If header is
    invalid, Excel is tried. Returns a sorted list of unique accession
    strings from column `column_name`, splitting on ";".

    Raises ValueError if column not found.
    """
    encodings = ["utf-8", "ISO-8859-1"]
    df = None
    for enc in encodings:
        try:
            df = pd.read_csv(path, encoding=enc, dtype=str)
            first_col = df.columns[0]
            if isinstance(first_col, str) and first_col.startswith(
                    "PK\x03\x04"):
                raise ValueError("Looks like XLSX, switch to read_excel.")
            break
        except (UnicodeDecodeError, ValueError):
            continue
        except Exception:
            df = None
            break

    if df is None:
        try:
            df = pd.read_excel(path, dtype=str)
        except Exception as e:
            raise ValueError(f"Could not read '{path}' as CSV or Excel: {e}")

    if column_name not in df.columns:
        raise ValueError(
            f'Input must have column "{column_name}". Found: {list(df.columns)}'
        )

    unique_accs = set()
    for cell in df[column_name].dropna().astype(str):
        for part in cell.split(";"):
            acc = part.strip()
            if acc:
                unique_accs.add(acc)

    return sorted(unique_accs)


def format_time(seconds: float) -> str:
    """
    Converts seconds to "Xm Ys" or "Z.Zs".
    """
    if seconds < 60:
        return f"{seconds:.1f} seconds"
    minutes = int(seconds // 60)
    secs = seconds - (minutes * 60)
    return f"{minutes} minute(s) {secs:.1f} seconds"


def main():
    parser = argparse.ArgumentParser(
        description="Generate a JSON dataset of protein annotations "
                    "from UniProt accessions in a specified column."
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to input CSV/Excel file containing a column of UniProt  "
             "accessions."
    )
    parser.add_argument(
        "--column", "-c", default="Master Protein Accessions",
        help='Name of the column that holds UniProt accessions (separate  '
             'multiple IDs by ";"). '
             'Default: "Master Protein Accessions".'
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Path where the JSON output will be written."
    )
    parser.add_argument(
        "--delay", "-d", type=float, default=0.1,
        help="Seconds to wait between API requests (default: 0.1)."
    )
    args = parser.parse_args()

    try:
        accessions = read_accessions_any_format(args.input, args.column)
    except ValueError as ve:
        print(f"Error: {ve}")
        return

    total_acc = len(accessions)
    if total_acc == 0:
        print(f"No accessions found in column '{args.column}'. Exiting.")
        return

    total_seconds_all = total_acc * args.delay
    eta_all = format_time(total_seconds_all)

    print(f"\nFound {total_acc} unique accession(s) in column '{args.column}'.")
    print(f"With a delay of {args.delay:.2f} seconds per accession,")
    print(f"fetching ALL would take approximately {eta_all}.\n")

    print("Choose an option:")
    print(f"  1) Generate dataset for ALL {total_acc} accession(s)")
    print("  2) Generate dataset for the FIRST N accession(s) (you specify N)")
    print("  3) Do NOT generate any dataset (exit)\n")

    choice = input("Enter 1, 2, or 3: ").strip()
    if choice == "3":
        print("Aborted. No data has been fetched or written.")
        return

    if choice == "1":
        to_process = accessions
        print(f"\nProceeding with ALL {total_acc} accession(s).")
    elif choice == "2":
        while True:
            n_str = input(f"Enter N (1–{total_acc}): ").strip()
            if not n_str.isdigit():
                print("→ Invalid: please enter a positive integer.")
                continue
            n = int(n_str)
            if n < 1 or n > total_acc:
                print(f"→ Invalid: N must be between 1 and {total_acc}.")
                continue
            break
        to_process = accessions[:n]
        eta_subset = format_time(n * args.delay)
        print(f"\nProceeding with the first {n} accession(s).")
        print(f"Estimated time for these {n} accession(s): {eta_subset}.\n")
    else:
        print("Invalid choice. Exiting without generating anything.")
        return

    num_to_process = len(to_process)
    all_records = []
    for idx, acc in enumerate(to_process, start=1):
        print(f"[{idx}/{num_to_process}] Processing {acc} ...")
        rec = build_protein_record(acc)
        all_records.append(rec)

    with open(args.output, "w", encoding="utf-8") as outf:
        json.dump({"proteins": all_records}, outf, ensure_ascii=False,
                  indent=4)

    print(f"\nFinished. JSON dataset written to {args.output}")


if __name__ == "__main__":
    main()
