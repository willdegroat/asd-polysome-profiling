from pathlib import Path

import openpyxl
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data' / 'reference' / 'fu2022SupplementaryTables.xlsx'
OUTPUT = ROOT / 'data' / 'geneSets' / 'fuAsd2022.gmt'

workbook = openpyxl.load_workbook(SOURCE, read_only = True, data_only = True)
sheet = workbook['Supplementary Table 11']
rows = sheet.iter_rows(values_only = True)
table = pd.DataFrame(rows, columns = next(rows))
workbook.close()

required = {'gene', 'ASD72', 'ASD185'}
missing = required.difference(table.columns)
if missing:
    raise ValueError(f'Fu table is missing columns: {sorted(missing)}')


def extract_genes(column: str, expected: int) -> list[str]:
    genes = (
        table.loc[table[column].eq(True), 'gene'].astype('string').str.strip()
    )
    if genes.isna().any() or genes.eq('').any():
        raise ValueError(f'{column} contains a missing gene symbol.')
    if len(genes) != expected or genes.nunique() != expected:
        raise ValueError(
            f'{column} must contain {expected} unique gene symbols.'
        )
    return sorted(genes.tolist())


genes72 = extract_genes('ASD72', 72)
genes185 = extract_genes('ASD185', 185)

if not set(genes72).issubset(genes185):
    raise ValueError('The Fu 72-gene tier must be nested in the 185-gene tier.')

source = 'https://doi.org/10.1038/s41588-022-01104-0'
lines = [
    '\t'.join(['Fu ASD FDR <= 0.001', source, *genes72]),
    '\t'.join(['Fu ASD FDR <= 0.05', source, *genes185]),
]
OUTPUT.write_text('\n'.join(lines) + '\n', encoding = 'utf-8')
