import { act, render, waitFor } from "@testing-library/react";
import { beforeEach, expect, it, vi } from "vitest";
import SingleProteinPage from "../../src/features/single-protein/SingleProteinPage";
import PhosphoSpecificPage from "../../src/features/phospho/PhosphoSpecificPage";
import ComparisonPipelinePage from "../../src/features/comparison/ComparisonPipelinePage";
import * as api from "../../src/lib/api";
import { serveDatasets, silenceNetwork } from "./apiMock";

vi.mock("../../src/lib/api", async (original) => (await import("./apiMock")).mockApiModule(original));

beforeEach(() => { silenceNetwork(); serveDatasets(api); });

it.each(["pearson", "venn"] as const)("does not rewrite valid restored %s selections when options arrive", async (activeTab) => {
  window.sessionStorage.setItem("pcopylot.ui.comparison.pipeline.v1", JSON.stringify({
    kind: "protein",
    pearson: { mode: "Single", sample1: "s1", sample2: "s2", condition1: "A", condition2: "B", alias1: "First", alias2: "Second" },
    venn: { mode: "Single", first: "s1", second: "s2", third: "", alias1: "First", alias2: "Second", alias3: "" },
  }));
  let resolve!: (value: Awaited<ReturnType<typeof api.getComparisonOptions>>) => void;
  vi.mocked(api.getComparisonOptions).mockImplementation(() => new Promise(done => { resolve = done; }));
  vi.mocked(api.getComparisonTable).mockResolvedValue({ rows: [] } as Awaited<ReturnType<typeof api.getComparisonTable>>);
  const writes = vi.spyOn(Storage.prototype, "setItem");
  render(<ComparisonPipelinePage activeTab={activeTab} />);
  await waitFor(() => expect(resolve).toBeTypeOf("function"));
  writes.mockClear();
  await act(async () => resolve({ samples: ["s1", "s2"], conditions: ["A", "B"], sampleLabels: {} } as Awaited<ReturnType<typeof api.getComparisonOptions>>));
  expect(writes.mock.calls.filter(([key]) => key.startsWith("pcopylot.ui.comparison.pipeline.v1"))).toEqual([]);
  expect(api.getComparisonTable).toHaveBeenCalledTimes(1);
  expect(vi.mocked(api.getComparisonTable).mock.calls[0][2]).toMatchObject({ alias1: "First", alias2: "Second" });
});

it.each(["boxplot", "lineplot"] as const)("does not restart the %s table request for already-valid restored selections", async (activeTab) => {
  window.sessionStorage.setItem("pcopylot.ui.single.protein.v1", JSON.stringify({
    kind: "protein", dataSource: "data",
    box: { identifier: "workflow", protein: "P1", conditions: ["B", "A"] },
    line: { identifier: "workflow", proteins: ["P2", "P1"], conditions: ["B", "A"] },
  }));
  let resolve!: (value: Awaited<ReturnType<typeof api.getSingleProteinOptions>>) => void;
  vi.mocked(api.getSingleProteinOptions).mockImplementation(() => new Promise(done => { resolve = done; }));
  vi.mocked(api.getSingleProteinTable).mockResolvedValue({ rows: [] } as Awaited<ReturnType<typeof api.getSingleProteinTable>>);
  render(<SingleProteinPage activeTab={activeTab} />);
  await waitFor(() => expect(resolve).toBeTypeOf("function"));
  await act(async () => resolve({
    proteins: ["P1", "P2"], conditions: ["A", "B"], proteinCount: 2, conditionCount: 2,
    identifier: "workflow", availableIdentifiers: [{ key: "workflow", label: "Protein" }],
    dataSource: "data", availableDataSources: ["data"], defaultDataSource: "data",
  }));
  expect(api.getSingleProteinTable).toHaveBeenCalledTimes(1);
  const call = vi.mocked(api.getSingleProteinTable).mock.calls[0];
  expect(call[2].conditions).toBe('["B","A"]');
  if (activeTab === "lineplot") expect(call[2].proteins).toBe('["P2","P1"]');
  else expect(call[2].protein).toBe("P1");
  expect(call[3]?.signal?.aborted).toBe(false);
});

it("does not schedule another phosphoprotein table request when restored conditions remain valid", async () => {
  window.sessionStorage.setItem("pcopylot.ui.phospho.specific.v1", JSON.stringify({
    dataSource: "data", phosprotRegulation: { condition1: "A", condition2: "B" },
  }));
  let resolve!: (value: Awaited<ReturnType<typeof api.getPhosphoOptions>>) => void;
  vi.mocked(api.getPhosphoOptions).mockImplementation(() => new Promise(done => { resolve = done; }));
  vi.mocked(api.getPhosphoTable).mockResolvedValue({ rows: [] } as Awaited<ReturnType<typeof api.getPhosphoTable>>);
  render(<PhosphoSpecificPage activeTab="phosprotRegulation" />);
  await waitFor(() => expect(resolve).toBeTypeOf("function"));
  await act(async () => resolve({
    conditions: ["A", "B"], samplesByCondition: {}, pairIdsBySample: {}, features: [],
    imputedAvailable: false, availableDataSources: ["data"], defaultDataSource: "data",
  }));
  // Let the existing 650 ms debounce expire so a redundant update cannot hide.
  await act(async () => { await new Promise(done => window.setTimeout(done, 750)); });
  expect(api.getPhosphoTable).toHaveBeenCalledTimes(1);
  expect(vi.mocked(api.getPhosphoTable).mock.calls[0][1]).toMatchObject({ condition1: "A", condition2: "B" });
});
