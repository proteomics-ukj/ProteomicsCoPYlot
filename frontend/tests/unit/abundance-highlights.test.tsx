import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";
import QCPipelinePage from "../../src/features/qc/QCPipelinePage";
import * as api from "../../src/lib/api";
import { serveDatasets, silenceNetwork } from "../component/apiMock";

vi.mock("../../src/lib/api", async (original) => (await import("../component/apiMock")).mockApiModule(original));

describe("abundance highlighting", () => {
  it("restores selection, commits custom input and sends both display switches", async () => {
    silenceNetwork();
    serveDatasets(api);
    vi.mocked(api.getQcPlotOptions).mockResolvedValue({
      conditions: ["A", "B"], imputedAvailable: false, samples: ["a1", "b1"],
      availableDataSources: ["data"], defaultDataSource: "data", centredPhosprotValues: false,
    } as Awaited<ReturnType<typeof api.getQcPlotOptions>>);
    sessionStorage.setItem("pcopylot.ui.qc.pipeline.v1", JSON.stringify({ abundance: {
      topHigh: 2, topLow: 3, search: "AKT1", highlightOutline: false, highlightLabels: true,
    } }));
    render(<QCPipelinePage activeTab="abundance" />);
    const requested = () => vi.mocked(globalThis.fetch).mock.calls.map(([url]) => new URL(String(url), "http://localhost"))
      .filter((url) => /\/abundance\.png$/.test(url.pathname));
    await waitFor(() => expect(requested().length).toBeGreaterThan(0));
    let query = requested().at(-1)!.searchParams;
    expect(query.get("topHigh")).toBe("2");
    expect(query.get("topLow")).toBe("3");
    expect(query.get("highlightOutline")).toBe("false");
    expect(query.get("highlightLabels")).toBe("true");
    const input = screen.getByLabelText("Highlight terms");
    fireEvent.change(input, { target: { value: "MTOR\nAKT1" } });
    fireEvent.blur(input);
    await waitFor(() => expect(requested().at(-1)!.searchParams.get("search")).toBe("MTOR\nAKT1"));
    fireEvent.click(screen.getByLabelText("Outline highlighted points"));
    await waitFor(() => expect(requested().at(-1)!.searchParams.get("highlightOutline")).toBe("true"));
    query = requested().at(-1)!.searchParams;
    expect(query.get("condition")).toBe("All Conditions");
  });
});

it("offers genes, keeps top-N together and puts static size after highlighting", async () => {
  silenceNetwork();
  serveDatasets(api);
  vi.mocked(api.getQcPlotOptions).mockResolvedValue({
    conditions: ["A"], imputedAvailable: false, availableDataSources: ["data"], defaultDataSource: "data",
    centredPhosprotValues: false, abundanceIdentifiers: [{ key: "workflow", label: "Protein Names" }, { key: "genes", label: "Gene Names" }],
  });
  render(<QCPipelinePage activeTab="abundance" />);
  await waitFor(() => expect(api.getQcPlotOptions).toHaveBeenCalled());
  await screen.findByRole("option", { name: "Gene Names" });
  fireEvent.change(screen.getByLabelText("Label column"), { target: { value: "genes" } });
  const high = screen.getByLabelText("Highlight top N highest abundance");
  const low = screen.getByLabelText("Highlight top N lowest abundance");
  expect(high.parentElement?.parentElement).toBe(low.parentElement?.parentElement);
  fireEvent.change(high, { target: { value: "2" } });
  fireEvent.blur(high);
  const terms = screen.getByLabelText("Highlight terms");
  const width = screen.getByLabelText("Width (cm)");
  expect(terms.compareDocumentPosition(width) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
  expect(screen.getByLabelText("Label highlighted points")).toBeTruthy();
  fireEvent.change(screen.getByLabelText("Type"), { target: { value: "Interactive" } });
  expect(screen.queryByLabelText("Label highlighted points")).toBeNull();
  expect(screen.getByLabelText("Outline highlighted points")).toBeTruthy();
  await waitFor(() => {
    const urls = vi.mocked(globalThis.fetch).mock.calls.map(([url]) => new URL(String(url), "http://localhost"));
    expect(urls.some(url => url.pathname.endsWith("abundance-interactive.html") && url.searchParams.get("identifier") === "genes" && url.searchParams.get("topHigh") === "2")).toBe(true);
  });
});
