import type { AnnotationResultResponse, ContaminantFilterConfig } from "../../lib/types";

export const CONTAMINANT_CATEGORIES = [
  [1, "Laboratory proteins"],
  [2, "Dust/contact proteins"],
  [3, "Molecular weight markers/Standard proteins"],
  [4, "Sigma-Aldrich Universal Protein Standard"],
  [5, "Common viral contaminants"],
] as const;

export const DEFAULT_CONTAMINANTS: ContaminantFilterConfig = {
  enabled: false, categories: [1, 2, 3, 4, 5], matchAccessions: true,
};

export function filterResultSummary(result: AnnotationResultResponse): string {
  const contaminants = result.filter?.contaminants;
  const removal = result.warnings?.map((warning) => /^Contaminant filter: removed ([\d,]+) rows/.exec(warning)).find(Boolean);
  const totalRemoved = Math.max(0, result.log2Rows - result.filteredRows);
  const contaminantRemoved = removal ? Number(removal[1].replaceAll(",", "")) : 0;
  const categories = CONTAMINANT_CATEGORIES.filter(([id]) => contaminants?.categories.includes(id)).map(([, label]) => label);
  if (contaminants?.enabled && !removal) {
    return `${totalRemoved.toLocaleString()} rows filtered; ${result.filteredRows.toLocaleString()} rows remain. Categories: ${categories.join(", ")}.`;
  }
  const featureRemoved = Math.max(0, totalRemoved - contaminantRemoved);
  return `${featureRemoved.toLocaleString()} rows filtered by feature; ${contaminantRemoved.toLocaleString()} rows filtered by contaminant.`
    + (contaminants?.enabled ? ` Categories: ${categories.join(", ")}.` : "");
}
