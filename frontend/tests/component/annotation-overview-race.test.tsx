import { act, render, screen, waitFor } from "@testing-library/react";
import { beforeEach, expect, it, vi } from "vitest";
import AnnotationPage from "../../src/features/metadata/AnnotationPage";
import * as api from "../../src/lib/api";
import { emitDatasetStateChanged } from "../../src/lib/datasetEvents";
import { DATASETS, silenceNetwork } from "./apiMock";

vi.mock("../../src/lib/api", async (original) => (await import("./apiMock")).mockApiModule(original));

type Overview = Awaited<ReturnType<typeof api.getPeptideOverview>>;
function deferred() {
  let resolve!: (value: Overview) => void;
  let reject!: (reason: Error) => void;
  const promise = new Promise<Overview>((yes, no) => { resolve = yes; reject = no; });
  return { promise, resolve, reject };
}
function datasets(filename: string | null) {
  return { ...DATASETS, protein: null, phospho: null,
    peptide: filename ? { filename, path: filename, rows: 1, columns: 1 } : null,
  } as Awaited<ReturnType<typeof api.getCurrentDatasets>>;
}
function overview(sample: string) {
  return { availableSamples: [sample], rows: 10, columns: 2 } as Overview;
}
beforeEach(() => {
  silenceNetwork();
  vi.mocked(api.getCurrentDatasets).mockResolvedValue(datasets("old.tsv"));
});

it.each(["success", "failure"])("ignores obsolete %s after a newer overview arrives", async (outcome) => {
  const old = deferred(), current = deferred();
  vi.mocked(api.getPeptideOverview).mockReturnValueOnce(old.promise).mockReturnValueOnce(current.promise);
  render(<AnnotationPage />);
  await waitFor(() => expect(api.getPeptideOverview).toHaveBeenCalledTimes(1));
  const signal = vi.mocked(api.getPeptideOverview).mock.calls[0][0]?.signal;
  act(() => emitDatasetStateChanged(datasets("new.tsv")));
  await waitFor(() => expect(api.getPeptideOverview).toHaveBeenCalledTimes(2));
  await act(async () => current.resolve(overview("new-sample")));
  expect(await screen.findByRole("option", { name: "new-sample" })).toBeTruthy();
  await act(async () => {
    if (outcome === "success") old.resolve(overview("old-sample"));
    else old.reject(new Error("obsolete failure"));
  });
  expect(screen.queryByRole("option", { name: "old-sample" })).toBeNull();
  expect(screen.getByRole("option", { name: "new-sample" })).toBeTruthy();
  expect(screen.queryByText("obsolete failure")).toBeNull();
  expect(signal?.aborted).toBe(true);
});

it.each(["success", "failure"])("cancels an outstanding overview on dataset removal: %s", async (outcome) => {
  const old = deferred();
  vi.mocked(api.getPeptideOverview).mockReturnValueOnce(old.promise);
  render(<AnnotationPage />);
  await waitFor(() => expect(api.getPeptideOverview).toHaveBeenCalledTimes(1));
  const signal = vi.mocked(api.getPeptideOverview).mock.calls[0][0]?.signal;
  act(() => emitDatasetStateChanged(datasets(null)));
  await act(async () => {
    if (outcome === "success") old.resolve(overview("old-sample"));
    else old.reject(new Error("removed dataset failure"));
  });
  expect(screen.queryByRole("option", { name: "old-sample" })).toBeNull();
  expect(screen.queryByText("removed dataset failure")).toBeNull();
  expect(signal?.aborted).toBe(true);
  expect(api.getPeptideOverview).toHaveBeenCalledTimes(1);
});

it("still displays an error from the current overview request", async () => {
  vi.mocked(api.getPeptideOverview).mockRejectedValue(new Error("current request failure"));
  render(<AnnotationPage />);
  expect(await screen.findByText("current request failure")).toBeTruthy();
});

it("cancels the overview when leaving the page", async () => {
  vi.mocked(api.getPeptideOverview).mockReturnValue(deferred().promise);
  const view = render(<AnnotationPage />);
  await waitFor(() => expect(api.getPeptideOverview).toHaveBeenCalledTimes(1));
  const signal = vi.mocked(api.getPeptideOverview).mock.calls[0][0]?.signal;
  view.unmount();
  expect(signal?.aborted).toBe(true);
});
