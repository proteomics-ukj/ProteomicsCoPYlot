import { useState } from "react";
import { fireEvent, render, screen } from "@testing-library/react";
import { expect, it } from "vitest";
import ContaminantFilter from "../../src/features/metadata/ContaminantFilter";
import { DEFAULT_CONTAMINANTS } from "../../src/features/metadata/contaminantOptions";

it("defaults off, selects database categories on enable, and allows combinations but not an empty selection", () => {
  function Harness() {
    const [value, setValue] = useState(DEFAULT_CONTAMINANTS);
    return <ContaminantFilter value={value} onChange={setValue} />;
  }
  render(<Harness />);
  const enabled = screen.getByLabelText("Filter contaminants") as HTMLInputElement;
  expect(enabled.checked).toBe(false);
  expect(screen.queryByText("Contaminant categories")).toBeNull();
  fireEvent.click(enabled);
  const labels = ["Laboratory proteins", "Dust/contact proteins", "Molecular weight markers/Standard proteins",
    "Sigma-Aldrich Universal Protein Standard", "Common viral contaminants", "Upstream-flagged contaminants"];
  for (const label of labels) expect((screen.getByLabelText(label) as HTMLInputElement).checked).toBe(true);
  for (const label of labels.slice(0, 4)) fireEvent.click(screen.getByLabelText(label));
  fireEvent.click(screen.getByLabelText("Upstream-flagged contaminants"));
  expect((screen.getByLabelText("Common viral contaminants") as HTMLInputElement).disabled).toBe(true);
  fireEvent.click(screen.getByLabelText("Laboratory proteins"));
  expect((screen.getByLabelText("Common viral contaminants") as HTMLInputElement).disabled).toBe(false);
  fireEvent.click(enabled);
  fireEvent.click(enabled);
  expect((screen.getByLabelText("Dust/contact proteins") as HTMLInputElement).checked).toBe(false);
});

it("allows upstream flags alone and retains the selection when toggled off and on", () => {
  function Harness() {
    const [value, setValue] = useState(DEFAULT_CONTAMINANTS);
    return <><ContaminantFilter value={value} onChange={setValue} /><output>{JSON.stringify(value.categories)}</output></>;
  }
  render(<Harness />);
  fireEvent.click(screen.getByLabelText("Filter contaminants"));
  const upstream = screen.getByLabelText("Upstream-flagged contaminants") as HTMLInputElement;
  expect(upstream.checked).toBe(true);
  for (const label of ["Laboratory proteins", "Dust/contact proteins", "Molecular weight markers/Standard proteins",
    "Sigma-Aldrich Universal Protein Standard", "Common viral contaminants"]) fireEvent.click(screen.getByLabelText(label));
  expect(screen.getByRole("status").textContent).toBe("[6]");
  expect(upstream.disabled).toBe(true);
  fireEvent.click(screen.getByLabelText("Filter contaminants"));
  fireEvent.click(screen.getByLabelText("Filter contaminants"));
  expect((screen.getByLabelText("Upstream-flagged contaminants") as HTMLInputElement).checked).toBe(true);
});
