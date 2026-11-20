import { COPYLOT_COLORS } from "./conditionColors";

// The previous Colors-tab default, in Matplotlib's tab10 order.
export const MATPLOTLIB_TAB10 = [
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
] as const;

// Logo anchors sampled from the official UKJ web logo, followed by derived
// shades with a soft preference for UKJ identity, balanced against separation.
// https://www.uniklinikum-jena.de/skin/medweb2016/images/logo.png
export const UKJ_COLORS = [
  "#004899", "#009FE3", "#96C11F", "#7E391B", "#0C3ADF",
  "#416024", "#244F60", "#D63DAD", "#23BFC7", "#CDA67A",
  "#2A89F4", "#55BE71", "#7C7E07", "#A6C15C", "#418DAA",
  "#407E07", "#4170AA", "#8B68DF", "#6CC3E5", "#2347C7",
] as const;

export const CONDITION_COLOR_SETS = [
  { value: "copylot", label: "Proteomics CoPYlot Colors", colors: COPYLOT_COLORS },
  { value: "ukj", label: "UKJ Colors", colors: UKJ_COLORS },
  { value: "tab10", label: "Matplotlib tab10 Colors", colors: MATPLOTLIB_TAB10 },
] as const;

// Older pages persisted their built-in defaults even without a user override.
export function paletteColorOverride(value: unknown, previousDefault: string): string {
  const color = typeof value === "string" ? value.trim() : "";
  return color.toUpperCase() === previousDefault.toUpperCase() ? "" : color;
}
