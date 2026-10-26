"""Exact, conservative matching against the bundled Cambridge contaminant snapshot."""
import csv
import hashlib
import io

# <shared-core> Shared exact matching for the app, package and portable exports.
from app.services.copylot_core import contaminant_mask
# </shared-core>

from app.services.reference_db import resolve_reference_db


def contaminant_reference():
    payload = resolve_reference_db("crap_ccp_metadata.csv", "Cambridge contaminant metadata").read_bytes()
    return {
        "database": "Cambridge crap_ccp v1.0",
        "source": "https://zenodo.org/records/15115102",
        "metadataSha256": hashlib.sha256(payload).hexdigest(),
        "records": list(csv.DictReader(io.StringIO(payload.decode("utf-8-sig")))),
        # Portable notebooks redistribute the lookup, so they must carry its
        # attribution and licence notices as well as the scientific records.
        "notices": {name: resolve_reference_db(name, "contaminant licence notice").read_text(encoding="utf-8")
                    for name in ("CAMBRIDGE_NOTICE.md", "CAMBRIDGE_LICENSE.txt", "CC-BY-4.0.txt")},
    }


def filter_peptide_contaminants(data, config, reference):
    return data.loc[~contaminant_mask(data, config, reference)].copy()
