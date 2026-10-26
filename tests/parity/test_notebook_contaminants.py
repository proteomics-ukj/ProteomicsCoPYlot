import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from _r_runtime import RSCRIPT
from proteomicscopylot.notebook_functions import filter_contaminants, filter_data, normalize_data
from app.services.contaminant_filter import contaminant_mask, contaminant_reference

ROOT = Path(__file__).resolve().parents[2]


def test_package_filter_preserves_indices_and_normalization_basis():
    data = pd.DataFrame({'ProteinNames': ['ALBU_HUMAN', 'ALBU_HUMAN;SAMPLE', 'SAMPLE', 'SAMPLE2'],
                         'A': [100., 10., 20., None], 'B': [200., 12., 22., 23.]}, index=[4, 7, 10, 13])
    meta = pd.DataFrame({'sample': ['A', 'B'], 'condition': ['one', 'one']})
    pd.testing.assert_frame_equal(filter_contaminants(data), data)
    kept = filter_contaminants(filter_data(data, meta, 2), enabled=True)
    assert kept.index.tolist() == [7, 10]
    normalized = normalize_data(data, meta, 'median', fit_rows=kept.index)
    assert normalized.loc[kept.index, 'A'].median() == pytest.approx(normalized.loc[kept.index, 'B'].median())
    expected = data.loc[~contaminant_mask(data, {'enabled': True}, contaminant_reference())]
    pd.testing.assert_frame_equal(filter_contaminants(data, enabled=True), expected)


@pytest.mark.requires_r()
@pytest.mark.parametrize('multi_category', [False, True])
@pytest.mark.parametrize('categories,accessions', [([1, 2, 3, 4, 5], True), ([1, 4], False), ([2], True)])
def test_r_and_python_match_the_app(tmp_path, categories, accessions, multi_category):
    headers = [line for line in (ROOT / 'data/db/crap_ccp.fasta').read_text().splitlines() if line.startswith('>')]
    names = ['ALBU_HUMAN', 'ALBU_HUMAN;SAMPLE', 'ALBU_HUMAN;', 'ALB', 'prefix_ALBU_HUMAN',
             None, 'display name', 'ENDOP_GLUC', 'SAMPLE', 'P02768-2', 'ALBU_HUMAN;SAMPLE'] + headers
    ids = [None, None, None, None, None, 'P02768', 'P02768', '000000', '000000', None, 'P02768'] + [None] * len(headers)
    extra = ['CON__P02768', 'contam_sp|P02768|ALBU_HUMAN', 'CON__P02768;SAMPLE',
             'REV__P02768', 'rev_sp|P02768|ALBU_HUMAN', 'CON__UNKNOWN']
    names += extra
    ids += [None] * len(extra)
    names += ['ENDOP_GLUC', 'SAMPLE']
    ids += ['CON__000000', 'CON__000000']
    data = pd.DataFrame({'ProteinNames': names, 'ProteinAccessions': ids, 'row_id': range(len(names))})
    input_path, output_path = tmp_path / 'input.csv', tmp_path / 'output.csv'
    data.to_csv(input_path, index=False)
    database = pd.read_csv(ROOT / 'data/db/crap_ccp_metadata.csv', dtype=str, keep_default_na=False)
    if multi_category:
        database.loc[database.entry_name == 'ALBU_HUMAN', 'gpm_category_id'] = '1;2'
    db_path = tmp_path / 'metadata.csv'
    database.to_csv(db_path, index=False)
    expected = filter_contaminants(data, enabled=True, categories=categories, match_accessions=accessions, database=database)
    script = f'''
source({json.dumps((ROOT / 'proteomicscopylot/notebook_functions.R').as_posix())})
x <- read.csv({json.dumps(input_path.as_posix())}, colClasses="character", na.strings="")
db <- read.csv({json.dumps(db_path.as_posix())}, colClasses="character", na.strings="")
y <- filter_contaminants(x, enabled=TRUE, categories=c({','.join(map(str, categories))}), match_accessions={'TRUE' if accessions else 'FALSE'}, database=db)
write.csv(y, {json.dumps(output_path.as_posix())}, row.names=FALSE)
'''
    completed = subprocess.run([RSCRIPT, '-e', script], capture_output=True, text=True, timeout=60)
    assert completed.returncode == 0, completed.stderr
    assert pd.read_csv(output_path).row_id.tolist() == expected.row_id.tolist()
