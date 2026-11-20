/**
 * The report's notebook export belongs to the report it was generated with.
 *
 * A notebook token names one rendered output. The button is enabled only for
 * the output it was captured from, so a notebook can never be exported beside a
 * report it does not reproduce; on the Summary page it sits next to the HTML
 * download and uses the token that came back with that report, and a notebook
 * that could not be built says why.
 *
 * Converted from a node --test file that compiled NotebookDownload.tsx with a
 * fake React, and matched SummaryPage.tsx's JSX with regular expressions; the
 * page is mounted and a report generated here.
 */
import { fireEvent, render, screen, within } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";
import NotebookDownload from "../../src/components/NotebookDownload";
import SummaryPage from "../../src/features/summary/SummaryPage";
import * as api from "../../src/lib/api";
import { serveDatasets, silenceNetwork } from "../component/apiMock";

vi.mock("../../src/lib/api", async (original) => (await import("../component/apiMock")).mockApiModule(original));

const TOKEN = "a".repeat(32);

describe("NotebookDownload", () => {
  it("is enabled only for the captured report", () => {
    const { rerender } = render(<NotebookDownload plotUrl="report-one" outputLabel="report" rendered={{ url: "report-one", token: "snapshot" }} />);
    const button = screen.getByRole("button") as HTMLButtonElement;
    expect(button.disabled).toBe(false);
    expect(button.textContent).toBe("Export notebook");
    expect(button.title).toMatch(/this report/);
    rerender(<NotebookDownload plotUrl="report-one" outputLabel="report" rendered={{ url: "report-two", token: "snapshot" }} />);
    expect((screen.getByRole("button") as HTMLButtonElement).disabled).toBe(true);
    rerender(<NotebookDownload plotUrl="report-one" outputLabel="report" rendered={{ url: "report-one", token: null }} />);
    expect((screen.getByRole("button") as HTMLButtonElement).disabled).toBe(true);
  });
});

describe("the Summary report", () => {
  beforeEach(() => {
    silenceNetwork();
    serveDatasets(api);
  });

  async function generate(response: Partial<Awaited<ReturnType<typeof api.generateSummaryReport>>>) {
    vi.mocked(api.generateSummaryReport).mockResolvedValue({ fileName: "report.html", html: "<html>Report</html>", warnings: [], ...response });
    render(<SummaryPage activeTab="report" />);
    fireEvent.click(screen.getByText("Generate HTML Report"));
    return screen.findByText("Download HTML Report");
  }

  it("puts the notebook export next to the HTML download, armed with the report's own token", async () => {
    const download = await generate({ notebookToken: TOKEN });
    const notebook = download.nextElementSibling as HTMLButtonElement;
    expect(notebook.textContent).toBe("Export notebook");
    expect(notebook.disabled).toBe(false);
    expect(notebook.title).toMatch(/this report/);
    vi.mocked(globalThis.fetch).mockClear();
    fireEvent.click(notebook);
    await vi.waitFor(() => expect(globalThis.fetch).toHaveBeenCalled());
    expect(String(vi.mocked(globalThis.fetch).mock.calls[0][0])).toContain(`/api/plots/notebooks/${TOKEN}.zip`);
  });

  it("disables the export and says why when the report came without a notebook", async () => {
    const download = await generate({ notebookToken: null, notebookError: "The notebook could not be built for this report." });
    expect((download.nextElementSibling as HTMLButtonElement).disabled).toBe(true);
    expect(screen.getByRole("alert").textContent).toBe("The notebook could not be built for this report.");
  });

  it("offers separate selected-profile and all-profile metadata downloads", async () => {
    vi.mocked(api.getSummaryTables).mockResolvedValue({
      protein: { kind: "protein", available: true, rows: 1, columns: 2, columnNames: ["sample", "condition"], table: [{ sample: "run1", condition: "A" }] },
      phospho: { kind: "phospho", available: false, rows: 0, columns: 0, columnNames: [], table: [] }, sampleNames: ["run1"],
    });
    vi.mocked(api.downloadCurrentAnnotationMetadataData).mockResolvedValue({ blob: new Blob(["metadata"]), filename: "metadata.xlsx" });
    render(<SummaryPage activeTab="tables" />);
    const selected = await screen.findByRole("button", { name: "Download selected profile (CSV)" });
    fireEvent.click(selected);
    await vi.waitFor(() => expect(api.downloadCurrentAnnotationMetadataData).toHaveBeenCalledWith("protein", "selected"));
    const all = await screen.findByRole("button", { name: "Download all profiles (Excel)" });
    fireEvent.click(all);
    await vi.waitFor(() => expect(api.downloadCurrentAnnotationMetadataData).toHaveBeenCalledWith("protein", "all"));
  });
});

describe("the report text box", () => {
  beforeEach(() => {
    silenceNetwork();
    serveDatasets(api);
  });

  it("asks the browser for en-US spell checking", async () => {
    // 2026-10-09: the browser's own dictionary underlines typos in the prose.
    render(<SummaryPage activeTab="text" />);
    const box = (await screen.findByLabelText("Introduction Text")) as HTMLTextAreaElement;
    expect(box.getAttribute("spellcheck")).toBe("true");
    expect(box.getAttribute("lang")).toBe("en-US");
    expect(screen.queryByRole("button", { name: "Add / Update" })).toBeNull();
  });

  it("keeps the saved section and draft while dataset options load", async () => {
    sessionStorage.setItem("pcopylot.ui.summary.editor.v1", JSON.stringify({
      selectedPlot: "CoverageProt", position: "Below", editorValue: "Recovered draft",
      textEntries: { Introduction: "Study introduction", CoverageProtBelow: "Older text" },
    }));
    const first = render(<SummaryPage activeTab="text" />);
    await vi.waitFor(() => expect((screen.getByLabelText("Select Plot") as HTMLSelectElement).value).toBe("CoverageProt"));
    expect((screen.getByLabelText("Text", { exact: true }) as HTMLTextAreaElement).value).toBe("Recovered draft");
    first.unmount();
    render(<SummaryPage activeTab="text" />);
    await vi.waitFor(() => expect((screen.getByLabelText("Select Plot") as HTMLSelectElement).value).toBe("CoverageProt"));
    expect((screen.getByLabelText("Text", { exact: true }) as HTMLTextAreaElement).value).toBe("Recovered draft");
    fireEvent.change(screen.getByLabelText("Select Plot"), { target: { value: "Introduction" } });
    expect((screen.getByLabelText("Introduction Text") as HTMLTextAreaElement).value).toBe("Study introduction");
  });

  it("saves text across sections, positions and remounts and includes it in the report", async () => {
    const { unmount } = render(<SummaryPage activeTab="text" />);
    await screen.findByRole("option", { name: "Coverage Plot - Protein Level" });
    fireEvent.change(screen.getByLabelText("Introduction Text"), { target: { value: "Study introduction" } });
    fireEvent.change(screen.getByLabelText("Select Plot"), { target: { value: "CoverageProt" } });
    fireEvent.change(screen.getByLabelText("Text", { exact: true }), { target: { value: "Before coverage\nSecond line" } });
    fireEvent.click(screen.getByRole("button", { name: "Below" }));
    fireEvent.change(screen.getByLabelText("Text", { exact: true }), { target: { value: "After coverage" } });
    fireEvent.click(screen.getByRole("button", { name: "Above" }));
    expect((screen.getByLabelText("Text", { exact: true }) as HTMLTextAreaElement).value).toBe("Before coverage\nSecond line");
    fireEvent.click(screen.getByLabelText("Show stored text entries"));
    const blocks = screen.getAllByRole("article");
    expect(blocks).toHaveLength(3);
    expect(within(blocks[0]).getByRole("heading").textContent).toBe("Introduction");
    expect(within(blocks[1]).getByRole("heading").textContent).toBe("Coverage Plot - Protein Level");
    expect(within(blocks[1]).getByText("Above the plot")).toBeTruthy();
    expect(within(blocks[1]).getByText("Before coverage Second line")).toBeTruthy();
    unmount();

    vi.mocked(api.generateSummaryReport).mockResolvedValue({ fileName: "report.html", html: "<html>Report</html>", warnings: [] });
    render(<SummaryPage activeTab="report" />);
    fireEvent.click(screen.getByText("Generate HTML Report"));
    await screen.findByText("Download HTML Report");
    expect(vi.mocked(api.generateSummaryReport).mock.calls.at(-1)?.[0].textEntries).toEqual({
      Introduction: "Study introduction", CoverageProtAbove: "Before coverage\nSecond line", CoverageProtBelow: "After coverage",
    });
  });

  it("clears saved text without bringing it back on a section switch", async () => {
    render(<SummaryPage activeTab="text" />);
    await screen.findByRole("option", { name: "Coverage Plot - Protein Level" });
    fireEvent.change(screen.getByLabelText("Introduction Text"), { target: { value: "Remove me" } });
    fireEvent.click(screen.getByRole("button", { name: "Delete" }));
    fireEvent.change(screen.getByLabelText("Select Plot"), { target: { value: "CoverageProt" } });
    fireEvent.change(screen.getByLabelText("Select Plot"), { target: { value: "Introduction" } });
    expect((screen.getByLabelText("Introduction Text") as HTMLTextAreaElement).value).toBe("");
    fireEvent.click(screen.getByLabelText("Show stored text entries"));
    expect(screen.getByText("No stored entries yet.")).toBeTruthy();
  });

  it("keeps shared custom sample labels while toggled off and submits them when enabled", async () => {
    vi.mocked(api.getSummaryTables).mockResolvedValue({
      protein: { kind: "protein", available: false, rows: 0, columns: 0, columnNames: [], table: [] },
      phospho: { kind: "phospho", available: false, rows: 0, columns: 0, columnNames: [], table: [] },
      sampleNames: ["run_1", "run_2"],
    });
    const view = render(<SummaryPage activeTab="text" />);
    const options = screen.getAllByRole("option").map((option) => option.textContent);
    expect(options.slice(0, 2)).toEqual(["Introduction", "Sample table"]);
    fireEvent.change(screen.getByLabelText("Select Plot"), { target: { value: "SampleTable" } });
    const label = await screen.findByLabelText("Custom label for run_1") as HTMLInputElement;
    const enabled = screen.getByLabelText("Add custom labels to report sample tables") as HTMLInputElement;
    expect(enabled.checked).toBe(false);
    expect(label.value).toBe("");
    expect(screen.queryByRole("button", { name: "Above" })).toBeNull();
    fireEvent.change(label, { target: { value: "Patient 1" } });
    fireEvent.click(enabled);
    fireEvent.click(enabled);
    view.unmount();
    const restored = render(<SummaryPage activeTab="text" />);
    expect((await screen.findByLabelText("Custom label for run_1") as HTMLInputElement).value).toBe("Patient 1");
    expect((screen.getByLabelText("Custom label for run_2") as HTMLInputElement).value).toBe("");
    expect((screen.getByLabelText("Add custom labels to report sample tables") as HTMLInputElement).checked).toBe(false);
    fireEvent.click(screen.getByLabelText("Add custom labels to report sample tables"));
    vi.mocked(api.generateSummaryReport).mockResolvedValue({ fileName: "report.html", html: "<html>Report</html>", warnings: [] });
    restored.rerender(<SummaryPage activeTab="report" />);
    fireEvent.click(screen.getByText("Generate HTML Report"));
    await screen.findByText("Download HTML Report");
    expect(vi.mocked(api.generateSummaryReport).mock.calls.at(-1)?.[0]).toMatchObject({
      includeSampleLabels: true, sampleLabels: { run_1: "Patient 1" },
    });
  });
});
