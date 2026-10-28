import { act, fireEvent, render, screen, waitFor } from "@testing-library/react";
import { expect, it, vi } from "vitest";
import * as api from "../../src/lib/api";
import AnnotationPage from "../../src/features/metadata/AnnotationPage";
import { serveDatasets, silenceNetwork } from "./apiMock";

vi.mock("../../src/lib/api", async (original) => (await import("./apiMock")).mockApiModule(original));

it("does not carry protein contaminant settings into a new phospho annotation", async () => {
  silenceNetwork();
  serveDatasets(api);
  vi.mocked(api.getCurrentAnnotation).mockRejectedValue(new Error("No annotation"));
  // Flush dataset loading and its filter-reset effect before interacting.
  // The dataset selector can become enabled before that effect has settled.
  await act(async () => { render(<AnnotationPage />); });
  await waitFor(() => expect((screen.getByLabelText("Dataset level") as HTMLSelectElement).disabled).toBe(false));
  fireEvent.click(screen.getByLabelText("Filter contaminants"));
  expect((screen.getByLabelText("Filter contaminants") as HTMLInputElement).checked).toBe(true);
  fireEvent.change(screen.getByLabelText("Dataset level"), { target: { value: "phospho" } });
  await waitFor(() => expect((screen.getByLabelText("Filter contaminants") as HTMLInputElement).checked).toBe(false));
});
