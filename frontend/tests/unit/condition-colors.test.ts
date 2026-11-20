import { describe, expect, it } from "vitest";
import { COPYLOT_COLORS, conditionColor } from "../../src/lib/conditionColors";
import { paletteColorOverride } from "../../src/lib/conditionColorSets";

describe("CoPYlot Balanced condition colors", () => {
  it("inherits the palette for legacy defaults but preserves custom plot colors", () => {
    expect(paletteColorOverride(undefined, "#2563eb")).toBe("");
    expect(paletteColorOverride("#2563EB", "#2563eb")).toBe("");
    expect(paletteColorOverride("#123456", "#2563eb")).toBe("#123456");
  });
  it("keeps forty distinct colors and repeats at the boundary", () => {
    expect(COPYLOT_COLORS).toHaveLength(40);
    expect(new Set(COPYLOT_COLORS).size).toBe(40);
    expect(conditionColor(0)).toBe("#2563EB");
    expect(conditionColor(1)).toBe("#FF4B5C");
    for (let i = 0; i < 120; i += 1) {
      expect(conditionColor(i)).toBe(COPYLOT_COLORS[i % 40]);
    }
  });
});
