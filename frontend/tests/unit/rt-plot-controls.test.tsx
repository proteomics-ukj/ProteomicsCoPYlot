import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";
import PeptideLevelPage from "../../src/features/peptide/PeptideLevelPage";
import * as api from "../../src/lib/api";
import { serveDatasets, silenceNetwork } from "../component/apiMock";

vi.mock("../../src/lib/api", async (original) => (await import("../component/apiMock")).mockApiModule(original));

beforeEach(() => {
  silenceNetwork();
  serveDatasets(api);
  vi.mocked(api.getPeptideOverview).mockResolvedValue({
    filename: "peptides.tsv", path: "peptides.tsv", rows: 10, columns: 4, columnNames: ["Sequence", "a1", "b1"],
    availableSamples: ["a1", "b1"], availableProteins: ["P1"], metadataLoaded: true, metadataFilename: "meta.tsv", warnings: [],
  } as unknown as Awaited<ReturnType<typeof api.getPeptideOverview>>);
});

const field = (label: string) => screen.getByLabelText(label) as HTMLInputElement | HTMLSelectElement;

/** The last RT plot the page asked for. */
function rtRequest() {
  const urls = vi.mocked(globalThis.fetch).mock.calls.map(([url]) => String(url))
    .concat(Array.from(document.querySelectorAll("iframe[src]")).map((frame) => frame.getAttribute("src")!))
    .filter((url) => /\/api\/plots\/peptide\/rt/.test(url));
  return urls.length ? new URL(urls.at(-1)!, "http://localhost") : null;
}

async function mount() {
  render(<PeptideLevelPage activeTab="rtPlot" />);
  await screen.findByLabelText("Mode");
  await waitFor(() => expect(rtRequest()).not.toBeNull());
}

describe("the RT plot options", () => {
  it("leave Mode and Bins usable and send them to the static plot", async () => {
    await mount();
    expect(screen.queryByLabelText("Type")).toBeNull();
    expect(field("Random Seed").value).toBe("1");
    expect(field("Mode").value).toBe("Scatter Plot");
    expect(field("Bins").value).toBe("100");
    expect(field("Mode").disabled).toBe(false);
    expect(field("Bins").disabled).toBe(false);
    fireEvent.change(field("Mode"), { target: { value: "Hexbin Plot" } });
    await waitFor(() => expect(rtRequest()!.searchParams.get("method")).toBe("Hexbin Plot"));
    expect(rtRequest()!.pathname).toMatch(/\/api\/plots\/peptide\/rt\.png$/);
    expect(rtRequest()!.searchParams.has("bins")).toBe(true);
  });

  it("sends a custom seed to the static plot", async () => {
    await mount();
    fireEvent.change(field("Random Seed"), { target: { value: "42" } });
    fireEvent.blur(field("Random Seed"));
    await waitFor(() => expect(rtRequest()!.searchParams.get("randomSeed")).toBe("42"));
    expect(screen.queryByText("Open Interactive Plot")).toBeNull();
  });
});
