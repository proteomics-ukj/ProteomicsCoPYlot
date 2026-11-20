import { IS_VIEWER_MODE } from "./appMode";

/**
 * What a notebook export ships besides its code and outputs.
 *
 * Both are off by default: an export then holds only what cannot be recomputed
 * from the original file. Kept per browser session and saved in a save state,
 * like the other page settings.
 */
export type NotebookExportOptions = { matrices: boolean; imputed: boolean };

const KEY = "pcopylot.ui.notebookExport.v1";

/**
 * Viewers do not get notebook exports (D11). Turning them on later is this
 * constant plus `VIEWER_NOTEBOOK_EXPORTS` in the backend's viewer_policy.py.
 */
export const NOTEBOOK_EXPORTS_ENABLED = !IS_VIEWER_MODE;

export function getNotebookExportOptions(): NotebookExportOptions {
  try {
    const stored = JSON.parse(globalThis.sessionStorage?.getItem(KEY) ?? "null") as Partial<NotebookExportOptions> | null;
    return { matrices: stored?.matrices === true, imputed: stored?.imputed === true };
  } catch {
    return { matrices: false, imputed: false };
  }
}

export function setNotebookExportOptions(options: NotebookExportOptions): void {
  try {
    globalThis.sessionStorage?.setItem(KEY, JSON.stringify({ matrices: options.matrices, imputed: options.imputed }));
  } catch {}
}

/** The query a notebook download carries for the current options ('' when both are off). */
export function notebookExportQuery(): string {
  const { matrices, imputed } = getNotebookExportOptions();
  const parts = [matrices ? "matrices=true" : "", imputed ? "imputed=true" : ""].filter(Boolean);
  return parts.length ? `?${parts.join("&")}` : "";
}
