/**
 * Phosphoprotein volcano plots cannot be added to the report.
 *
 * The report has a protein and a phosphosite level only, so the panel offers
 * no report buttons for the phosphoprotein level, while protein and phospho
 * keep theirs.
 */
import { render, screen, waitFor } from "@testing-library/react";
import { beforeEach, expect, it, vi } from "vitest";
import { VolcanoPanel } from "../../src/features/stats/VolcanoPanel";
import * as api from "../../src/lib/api";
import type { StatisticalOptionsResponse, VolcanoResultResponse } from "../../src/lib/types";
import { serveDatasets, silenceNetwork, statisticalOptions } from "../component/apiMock";

vi.mock("../../src/lib/api", async (original) => (await import("../component/apiMock")).mockApiModule(original));

const options = (kind: string): StatisticalOptionsResponse => statisticalOptions({ kind: kind as StatisticalOptionsResponse["kind"] });

beforeEach(() => {
  silenceNetwork();
  serveDatasets(api);
  vi.mocked(api.getStatisticalOptions).mockImplementation(async (kind) => options(String(kind)));
  vi.mocked(api.runVolcanoAnalysis).mockResolvedValue({
    kind: "protein", sourceUsed: "data", totalRows: 1, upregulatedCount: 0, downregulatedCount: 0,
    notSignificantCount: 1, rows: [{ ProteinNames: "P1", log2FC: 0.1, pval: 0.5 }], warnings: [],
  } as unknown as VolcanoResultResponse);
  sessionStorage.setItem("pcopylot.ui.stats.volcano.v1", JSON.stringify({ condition1: "A", condition2: "B" }));
});

it("offers no report buttons for phosphoprotein volcanoes", async () => {
  render(<VolcanoPanel kindOptions={[{ value: "phosprot", label: "Phosphoprotein" }]} />);
  await screen.findByLabelText("Highlight top N up-regulated");
  expect(screen.queryByRole("button", { name: "Add Volcano To Report" })).toBeNull();
  expect(screen.queryByRole("button", { name: "Remove Volcano From Report" })).toBeNull();
  expect(screen.queryByRole("button", { name: "Delete All From Report" })).toBeNull();
});

it.each(["protein", "phospho"])("keeps the report buttons for %s volcanoes", async (kind) => {
  render(<VolcanoPanel kindOptions={[{ value: kind, label: kind }]} />);
  await waitFor(() => expect(screen.getByRole("button", { name: "Add Volcano To Report" })).toBeTruthy());
  expect(screen.getByRole("button", { name: "Remove Volcano From Report" })).toBeTruthy();
  expect(screen.getByRole("button", { name: "Delete All From Report" })).toBeTruthy();
});
import { fireEvent } from "@testing-library/react";
import { listReportEntries } from "../../src/lib/reportState";

it("uses the same top-N selection for interactive outlines and saved reports", async () => {
  render(<VolcanoPanel kindOptions={[{ value: "protein", label: "Protein" }]} />);
  await screen.findByLabelText("Highlight top N up-regulated");
  fireEvent.change(screen.getByLabelText("Type"), { target: { value: "Interactive" } });
  fireEvent.change(screen.getByLabelText("Highlight top N up-regulated"), { target: { value: "3" } });
  fireEvent.blur(screen.getByLabelText("Highlight top N up-regulated"));
  expect(screen.getByLabelText("Outline highlighted points")).toBeTruthy();
  expect(screen.queryByLabelText("Label highlighted points")).toBeNull();
  await waitFor(() => {
    const urls = vi.mocked(globalThis.fetch).mock.calls.map(([url]) => new URL(String(url), "http://localhost"));
    expect(urls.some(url => url.pathname.endsWith("volcano.html") && url.searchParams.get("labelTopUp") === "3")).toBe(true);
  });
  fireEvent.click(screen.getByRole("button", { name: "Add Volcano To Report" }));
  expect(listReportEntries().some(entry =>
    entry.type === "volcano" && entry.request.labelTopUp === 3 && entry.render.type === "Interactive")).toBe(true);
});
