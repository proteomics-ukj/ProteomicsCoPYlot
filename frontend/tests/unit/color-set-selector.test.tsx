import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";
import ColorsPage from "../../src/features/data/ColorsPage";
import * as api from "../../src/lib/api";
import { COPYLOT_COLORS } from "../../src/lib/conditionColors";
import { MATPLOTLIB_TAB10, UKJ_COLORS } from "../../src/lib/conditionColorSets";

const snapshot = vi.hoisted(() => ({
  availableKinds: ["protein", "phospho"],
  kindOptions: [{ value: "protein", label: "Protein" }, { value: "phospho", label: "Phosphosite" }],
}));
vi.mock("../../src/lib/datasetAvailability", () => ({ useCurrentDatasetsSnapshot: () => snapshot }));
vi.mock("../../src/lib/api", () => ({
  getQcPlotOptions: vi.fn(), getConditionPalette: vi.fn(), setConditionPalette: vi.fn(),
  isAbortError: () => false,
}));

const conditions = Array.from({ length: 41 }, (_, i) => `C${i + 1}`);
let saved: Record<string, Record<string, string>>;

beforeEach(() => {
  vi.clearAllMocks();
  saved = { protein: {}, phospho: {} };
  vi.mocked(api.getQcPlotOptions).mockResolvedValue({ conditions, samples: [] } as unknown as Awaited<ReturnType<typeof api.getQcPlotOptions>>);
  vi.mocked(api.getConditionPalette).mockImplementation(async (kind) => ({ kind, palette: saved[kind] }));
  vi.mocked(api.setConditionPalette).mockImplementation(async (kind, palette) => {
    saved[kind] = { ...palette };
    return { kind, palette };
  });
});

async function mount() {
  const view = render(<ColorsPage />);
  await waitFor(() => expect((screen.getByLabelText("Color set") as HTMLSelectElement).disabled).toBe(false));
  return view;
}

async function choose(value: string) {
  fireEvent.change(screen.getByLabelText("Color set"), { target: { value } });
  await waitFor(() => expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe(value));
}

describe("Colors tab color sets", () => {
  it("saves and reloads UKJ with twenty unique colors and the three logo anchors", async () => {
    expect(UKJ_COLORS).toHaveLength(20);
    expect(new Set(UKJ_COLORS).size).toBe(20);
    expect(UKJ_COLORS.slice(0, 3)).toEqual(["#004899", "#009FE3", "#96C11F"]);
    const view = await mount();
    expect(screen.getByRole("option", { name: "UKJ Colors" })).toBeTruthy();
    await choose("ukj");
    conditions.forEach((condition, i) => expect(saved.protein[condition]).toBe(UKJ_COLORS[i % 20]));
    expect(saved.phospho).toEqual({});
    view.unmount();
    await mount();
    expect((screen.getByLabelText("Color set") as HTMLSelectElement).value).toBe("ukj");
    expect((screen.getByLabelText("C21 color as hex") as HTMLInputElement).value).toBe(UKJ_COLORS[0]);
    fireEvent.click(screen.getByRole("button", { name: "Reset Default Palette" }));
    await waitFor(() => expect((screen.getByLabelText("Color set") as HTMLSelectElement).value).toBe("copylot"));
  });

  it("defaults to CoPYlot and applies, saves and reloads the previous tab10 palette", async () => {
    const view = await mount();
    expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe("copylot");
    expect(screen.getByRole("option", { name: "Proteomics CoPYlot Colors" })).toBeTruthy();
    expect(screen.getByRole("option", { name: "Matplotlib tab10 Colors" })).toBeTruthy();
    expect((screen.getByLabelText("C41 color as hex") as HTMLInputElement).value).toBe(COPYLOT_COLORS[0]);
    await choose("tab10");
    conditions.forEach((condition, i) => expect(saved.protein[condition]).toBe(MATPLOTLIB_TAB10[i % 10]));
    expect(saved.phospho).toEqual({});
    view.unmount();
    await mount();
    expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe("tab10");
    expect((screen.getByLabelText("C11 color as hex") as HTMLInputElement).value).toBe(MATPLOTLIB_TAB10[0]);
  });

  it("preserves custom colors and clears incomplete hex drafts when resetting", async () => {
    saved.protein = { C1: "#123456" };
    await mount();
    expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe("custom");
    expect(api.setConditionPalette).not.toHaveBeenCalled();
    fireEvent.change(screen.getByLabelText("C1 color as hex"), { target: { value: "#12" } });
    fireEvent.click(screen.getByRole("button", { name: "Reset Default Palette" }));
    await waitFor(() => expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe("copylot"));
    expect((screen.getByLabelText("C1 color as hex") as HTMLInputElement).value).toBe(COPYLOT_COLORS[0]);
    conditions.forEach((condition, i) => expect(saved.protein[condition]).toBe(COPYLOT_COLORS[i % 40]));
    fireEvent.change(screen.getByLabelText("C1 color"), { target: { value: "#abcdef" } });
    expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe("custom");
    expect(saved.protein.C1).toBe("#ABCDEF");
  });

  it("keeps the current palette when saving a selected set fails", async () => {
    await mount();
    vi.mocked(api.setConditionPalette).mockRejectedValueOnce(new Error("Save unavailable"));
    fireEvent.change(screen.getByLabelText("Color set"), { target: { value: "tab10" } });
    expect(await screen.findByText("Save unavailable")).toBeTruthy();
    expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe("copylot");
    expect((screen.getByLabelText("C1 color as hex") as HTMLInputElement).value).toBe(COPYLOT_COLORS[0]);
    expect(saved.protein).toEqual({});
  });

  it("copies the selected set to another dataset level", async () => {
    await mount();
    await choose("tab10");
    fireEvent.click(screen.getByRole("button", { name: "Copy Colors" }));
    await waitFor(() => expect(saved.phospho).toEqual(saved.protein));
    fireEvent.change(screen.getByLabelText("Dataset level"), { target: { value: "phospho" } });
    await waitFor(() => expect((screen.getByLabelText("Color set") as HTMLSelectElement).disabled).toBe(false));
    expect((screen.getByLabelText("Color set") as HTMLInputElement).value).toBe("tab10");
  });
});
