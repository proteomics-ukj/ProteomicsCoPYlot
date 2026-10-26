import type { ContaminantFilterConfig } from "../../lib/types";

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
