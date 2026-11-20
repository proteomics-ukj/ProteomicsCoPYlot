/**
 * The list of plots added to the report.
 *
 * Add stores the plot as it is; Remove matches on the page selection only and
 * takes every entry that matches; Delete all is limited to one plot type of
 * one workflow.
 */
import { beforeEach, expect, it } from "vitest";
import {
  addReportEntry,
  countReportEntries,
  deleteAllReportEntries,
  getSummaryReportContext,
  hasReportEntryMatching,
  listReportEntries,
  removeMatchingReportEntries,
  syncSummaryReportStateWithDatasets,
} from "../../src/lib/reportState";
import type { CurrentDatasetsResponse, ReportEntry } from "../../src/lib/types";

const volcano = (overrides: Partial<ReportEntry> = {}, request: Record<string, unknown> = {}): ReportEntry => ({
  type: "volcano",
  kind: "protein",
  request: {
    dataSource: "imputed", condition1: "A", condition2: "B", identifier: "workflow", testType: "unpaired",
    pValueThreshold: 0.05, log2fcThreshold: 1, useUncorrected: false, proteinCorrection: false,
    highlightTerms: [], labelTopUp: 0, labelTopDown: 0, ...request,
  },
  render: { type: "Normal", widthCm: 18, heightCm: 14, dpi: 300 },
  advanced: null,
  ...overrides,
});

beforeEach(() => sessionStorage.clear());

it("refuses only a fully identical entry", () => {
  expect(addReportEntry(volcano())).toBe(true);
  expect(addReportEntry(volcano())).toBe(false);
  expect(addReportEntry(volcano({ advanced: { title: "Mine" } }))).toBe(true);
  expect(addReportEntry(volcano({}, { highlightTerms: ["P1"] }))).toBe(true);
  expect(addReportEntry(volcano({ render: { type: "Normal", widthCm: 10, heightCm: 14, dpi: 300 } }))).toBe(true);
  expect(countReportEntries("volcano", "protein")).toBe(4);
});

it("removes every entry with the same selection, whatever its appearance, highlights or size", () => {
  addReportEntry(volcano());
  addReportEntry(volcano({ advanced: { title: "Mine" } }));
  addReportEntry(volcano({}, { highlightTerms: ["P1"], labelTopUp: 5 }));
  addReportEntry(volcano({ render: { type: "Normal", widthCm: 10, heightCm: 14, dpi: 300 } }));
  // A different plot: other label column, other render type, other level.
  addReportEntry(volcano({}, { identifier: "genes" }));
  addReportEntry(volcano({ render: { type: "Interactive" } }));
  addReportEntry(volcano({ kind: "phospho" }));
  expect(removeMatchingReportEntries(volcano({ advanced: { title: "Other" } }, { highlightTerms: ["X"] }))).toBe(4);
  expect(listReportEntries()).toHaveLength(3);
});

it("ignores the alias and the colour of a comparison when removing", () => {
  const pearson = (alias1: string, color: string): ReportEntry => ({
    type: "pearson", kind: "protein", advanced: null, render: { type: "Normal" },
    request: { mode: "Condition", sample1: "", sample2: "", condition1: "A", condition2: "B", alias1, alias2: "", color, header: true },
  });
  addReportEntry(pearson("first", "#111111"));
  addReportEntry(pearson("second", "#222222"));
  expect(removeMatchingReportEntries(pearson("", ""))).toBe(2);
});

it("deletes all entries of one plot type for one workflow only", () => {
  addReportEntry(volcano());
  addReportEntry(volcano({}, { condition2: "C" }));
  addReportEntry(volcano({ kind: "phospho" }));
  addReportEntry({ type: "venn", kind: "protein", request: { mode: "Single", first: "s1", second: "s2" }, render: { type: "Normal" }, advanced: null });
  expect(deleteAllReportEntries("volcano", "protein")).toBe(2);
  expect(countReportEntries("volcano", "phospho")).toBe(1);
  expect(countReportEntries("venn", "protein")).toBe(1);
});

it("carries volcano entries saved before the other plot types existed", () => {
  sessionStorage.setItem("pcopylot.report.volcanoLog.v1", JSON.stringify([
    { kind: "protein", control: false, dataSource: "filtered", condition1: "A", condition2: "B", identifier: "workflow",
      pValueThreshold: 0.05, log2fcThreshold: 1, testType: "unpaired", useUncorrected: false, highlightTerms: [],
      renderType: "Normal", labelTopUp: 2 },
    { kind: "protein", control: true, dataSource: "data", condition1: "A", condition2: "B" },
  ]));
  const entries = listReportEntries();
  expect(entries).toHaveLength(1);
  expect(entries[0]).toMatchObject({ type: "volcano", kind: "protein", render: { type: "Normal" }, advanced: null });
  expect(entries[0].request).toMatchObject({ dataSource: "filtered", labelTopUp: 2 });
  expect(sessionStorage.getItem("pcopylot.report.volcanoLog.v1")).toBeNull();
  expect(getSummaryReportContext().entries).toHaveLength(1);
});

it("finds the plot an ORA is based on by its test settings", () => {
  addReportEntry(volcano({ advanced: { title: "Mine" } }, { highlightTerms: ["P1"] }));
  const keys = ["dataSource", "condition1", "condition2", "testType", "pValueThreshold", "log2fcThreshold", "useUncorrected"];
  const ora = { dataSource: "imputed", condition1: "A", condition2: "B", testType: "unpaired", pValueThreshold: 0.05, log2fcThreshold: 1, useUncorrected: false };
  expect(hasReportEntryMatching("volcano", "protein", ora, keys)).toBe(true);
  expect(hasReportEntryMatching("volcano", "protein", { ...ora, log2fcThreshold: 2 }, keys)).toBe(false);
  expect(hasReportEntryMatching("volcano", "phospho", ora, keys)).toBe(false);
});

const datasets = (peptidePath: string, proteinRows = 10): CurrentDatasetsResponse => ({
  activePackage: "Data", activeMetadataProfile: "Metadata",
  protein: { filename: "protein.tsv", rows: proteinRows, columns: 5 },
  phospho: null, phosprot: null,
  peptide: { filename: "peptide.tsv", path: peptidePath },
} as unknown as CurrentDatasetsResponse);

it("keeps the list when only the peptide file's server path changes, and clears it for another dataset", () => {
  syncSummaryReportStateWithDatasets(datasets("/tmp/a/peptide.tsv"));
  addReportEntry(volcano());
  syncSummaryReportStateWithDatasets(datasets("/tmp/restored/peptide.tsv"));
  expect(listReportEntries()).toHaveLength(1);
  syncSummaryReportStateWithDatasets(datasets("/tmp/restored/peptide.tsv", 11));
  expect(listReportEntries()).toHaveLength(0);
});
