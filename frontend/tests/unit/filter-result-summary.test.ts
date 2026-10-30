import { expect, it } from "vitest";
import { filterResultSummary } from "../../src/features/metadata/contaminantOptions";
import type { AnnotationResultResponse } from "../../src/lib/types";

it("separates feature and contaminant removals without exposing provenance prose", () => {
  const result = { log2Rows: 100, filteredRows: 70,
    warnings: ["Contaminant filter: removed 10 rows after the valid-value filter; metadata SHA-256: secret."],
    filter: { contaminants: { enabled: true, categories: [1, 5] } },
  } as AnnotationResultResponse;
  expect(filterResultSummary(result)).toBe("20 rows filtered by feature; 10 rows filtered by contaminant. Categories: Laboratory proteins, Common viral contaminants.");
});

it("does not invent contaminant counts for legacy results without a count", () => {
  const result = { log2Rows: 100, filteredRows: 70, warnings: [],
    filter: { contaminants: { enabled: true, categories: [1] } },
  } as AnnotationResultResponse;
  expect(filterResultSummary(result)).toBe("30 rows filtered; 70 rows remain. Categories: Laboratory proteins.");
});
