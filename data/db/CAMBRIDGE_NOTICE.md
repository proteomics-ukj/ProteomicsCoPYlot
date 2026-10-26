# Cambridge CCP contaminant database

`crap_ccp.fasta` is preserved byte-for-byte from the following deposit:

Gatto, Laurent (2025). *Proteomics contaminant databases*, version 1.0.
Zenodo, published 31 March 2025. https://doi.org/10.5281/zenodo.15115102
Download: https://zenodo.org/api/records/15115102/files/crap_ccp.fasta/content

The deposit declares CC BY 4.0. See `CC-BY-4.0.txt` and
https://creativecommons.org/licenses/by/4.0/. Redistribution and adaptation,
including commercial use, are permitted subject to those terms.

The Cambridge Centre for Proteomics curated the CCP collection distributed
through camprotR. Credit Charlotte Dawson (maintainer) and Tom Smith (coauthor).
Gatto archived the files extracted from camprotR; he is not credited here as
originator of the Cambridge collection.
https://cambridgecentreforproteomics.github.io/camprotR/articles/crap.html
https://github.com/CambridgeCentreForProteomics/camprotR

The camprotR MIT notice, Copyright (c) 2020 Cambridge Centre for Proteomics,
is retained in `CAMBRIDGE_LICENSE.txt` from commit
197a44685aba5e490cab3877f4ac7ec3ef3e989a.
Underlying UniProt records: credit The UniProt Consortium, *UniProt: the Universal
Protein Knowledgebase*. UniProt applies CC BY 4.0 to copyrightable database parts:
https://www.uniprot.org/help/license
The CCP collection is largely based on GPM cRAP; credit The Global Proteome Machine
Organization for that foundation: https://www.thegpm.org/crap/

The snapshot contains 125 sequences: 123 distinct non-placeholder accessions and
two commercial sequences (Endoproteinase GluC, NEB P8100S; recombinant Lys-C,
Promega V167A) using placeholder accession 000000. Original headers are retained.
The deposit does not specify the UniProt release used to generate this snapshot;
its publication date must not be presented as a UniProt release date.

CoPYlot has made no changes to the FASTA. It is bundled as reference data;
optional contaminant filtering uses its identifiers and the category metadata below.
Filtering is off by default. Listing a protein
here does not establish that it is contamination in a particular experiment.
The data retain their source terms independently of CoPYlot's GPL-3.0 licence.
No upstream endorsement is implied. The data are provided without warranties;
upstream notices and any applicable third-party rights remain in force.

SHA-256: `e9da6a4361fac82331c63fd5d1a9534b48310ebf8136ebc245fd3451a1b230d0`
Upstream MD5: `e56af49e6037c2c1fb84beba4b25a78c`

## Derived category metadata

`crap_ccp_metadata.csv` is a CoPYlot-generated annotation of the unchanged Cambridge
snapshot. The 111 direct category matches come from The Global Proteome Machine Organization,
*cRAP protein sequences*, version 2012.01.01: https://www.thegpm.org/crap/
Copyright notice on source page: Copyright 2004-2011, The Global Proteome Machine
Organization. The source page does not state an explicit redistribution licence.
Only factual identifiers and category assignments are extracted, not its protein
descriptions, webpage or FASTA. No CC BY/MIT permission from GPM is asserted.
Cambridge/UniProt attribution and terms above continue to apply to their fields.

The other 14 assignments are labelled scientific decision. These are CoPYlot
annotations, not GPM or Cambridge assignments. All 125 records have a category;
the two 000000 commercial records remain distinct by Cambridge ID and entry name.
Categories describe potential contamination sources, not instructions to remove
proteins. Filtering requires the user's explicit selection in Filter Settings.
