/**
 * The Add / Remove / Delete all control every addable plot shares.
 */
import { fireEvent, render, screen } from "@testing-library/react";
import { beforeEach, expect, it } from "vitest";
import ReportEntryControl from "../../src/components/ReportEntryControl";
import { addReportEntry, hasReportEntryMatching, listReportEntries } from "../../src/lib/reportState";
import type { ReportEntry } from "../../src/lib/types";

const venn = (kind: ReportEntry["kind"], second: string, color1 = ""): ReportEntry => ({
  type: "venn", kind, request: { mode: "Condition", first: "A", second, third: "", minPresent: 1, header: true, color1 },
  render: { type: "Normal", widthCm: 15, heightCm: 12, dpi: 300 }, advanced: null,
});

function control(kind: string, entries: ReportEntry[], notice?: () => string | null) {
  return render(
    <ReportEntryControl type="venn" kind={kind} kindLabel="Protein" title="Venn" noun="Venn diagram"
      buildEntries={() => entries} invalidMessage="Please choose a valid selection before adding this plot." notice={notice} />
  );
}

beforeEach(() => sessionStorage.clear());

it("adds, refuses a duplicate, removes by selection and counts the workflow's entries", () => {
  addReportEntry(venn("protein", "B", "#123456"));
  addReportEntry(venn("phospho", "B"));
  control("protein", [venn("protein", "B")]);
  expect(screen.getByText("Saved entries for Protein: 1")).toBeTruthy();
  fireEvent.click(screen.getByRole("button", { name: "Add Venn To Report" }));
  expect(screen.getByText("Venn diagram added to report list.")).toBeTruthy();
  fireEvent.click(screen.getByRole("button", { name: "Add Venn To Report" }));
  expect(screen.getByText("This Venn diagram configuration is already in the report list.")).toBeTruthy();
  expect(screen.getByText("Saved entries for Protein: 2")).toBeTruthy();
  fireEvent.click(screen.getByRole("button", { name: "Remove Venn From Report" }));
  expect(screen.getByText("Removed 2 Venn diagram entries from report list.")).toBeTruthy();
  expect(listReportEntries()).toEqual([venn("phospho", "B")]);
});

it("deletes all entries of this plot for the selected workflow only", () => {
  addReportEntry(venn("protein", "B"));
  addReportEntry(venn("protein", "C"));
  addReportEntry(venn("phospho", "B"));
  control("protein", []);
  fireEvent.click(screen.getByRole("button", { name: "Delete All From Report" }));
  expect(screen.getByText("Saved entries for Protein: 0")).toBeTruthy();
  expect(listReportEntries()).toEqual([venn("phospho", "B")]);
});

it("says so when the plot cannot be added yet", () => {
  control("protein", []);
  fireEvent.click(screen.getByRole("button", { name: "Add Venn To Report" }));
  expect(screen.getByText("Please choose a valid selection before adding this plot.")).toBeTruthy();
});

it("offers nothing for the phosphoprotein level", () => {
  const { container } = control("phosprot", [venn("protein", "B")]);
  expect(container.innerHTML).toBe("");
});

it("shows a notice about the saved entries once it applies", () => {
  const notice = () =>
    hasReportEntryMatching("venn", "protein", { second: "B" }, ["second"]) ? "The volcano this ORA is based on is not in the report." : null;
  control("protein", [venn("protein", "B")], notice);
  expect(screen.queryByRole("status")).toBeNull();
  fireEvent.click(screen.getByRole("button", { name: "Add Venn To Report" }));
  expect(screen.getByRole("status").textContent).toBe("The volcano this ORA is based on is not in the report.");
});
