const UI_STATE_PREFIX = "pcopylot.ui.";

function storageAvailable(): boolean {
  try {
    return typeof globalThis.sessionStorage !== "undefined";
  } catch {
    return false;
  }
}

function fullKey(key: string): string {
  const normalized = String(key || "").trim();
  return normalized.startsWith(UI_STATE_PREFIX) ? normalized : `${UI_STATE_PREFIX}${normalized}`;
}

export function readPersistedUiState<T>(key: string, fallback: T): T {
  if (!storageAvailable()) return fallback;
  const storageKey = fullKey(key);
  try {
    const raw = globalThis.sessionStorage.getItem(storageKey);
    if (!raw) return fallback;
    const parsed = JSON.parse(raw) as T;
    return parsed ?? fallback;
  } catch {
    return fallback;
  }
}

export function writePersistedUiState(key: string, value: unknown): void {
  if (!storageAvailable()) return;
  const storageKey = fullKey(key);
  try {
    globalThis.sessionStorage.setItem(storageKey, JSON.stringify(value));
  } catch {
    // ignore sessionStorage failures
  }
}

// Older sessions kept one settings object with its owning dataset kind.
export function readWorkflowState(key: string, kind: string): Record<string, unknown> {
  const saved = readPersistedUiState<Record<string, unknown> | null>(key + "." + kind, null);
  if (saved) return saved;
  const previous = readPersistedUiState<Record<string, unknown>>(key, {});
  return (previous.kind ?? "protein") === kind ? previous : {};
}

export function writeWorkflowState(key: string, state: Record<string, unknown> & { kind: string }): void {
  const previous = readPersistedUiState<Record<string, unknown> | null>(key, null);
  if (previous && previous.kind !== state.kind) {
    const oldKey = key + "." + (previous.kind ?? "protein");
    if (!readPersistedUiState(oldKey, null)) writePersistedUiState(oldKey, previous);
  }
  writePersistedUiState(key + "." + state.kind, state);
  writePersistedUiState(key, state);
}

export function workflowPageForPlot(plot: string): string | undefined {
  const [section, tab, scope] = plot.split(":");
  if (scope) return undefined;
  if (section === "stats") return "stats." + tab + ".v1";
  if (section === "comparison") return "comparison.pipeline.v1";
  if (section === "singleProtein") return "single.protein.v1";
  if (plot === "data:distribution") return "data.distribution.v1";
}

export function lastPlotWorkflow(plot: string, page: string): string {
  if (plot === "singleProtein:heatmap") return "phospho";
  return readPersistedUiState<{ kind?: string }>(page, {}).kind ?? "protein";
}

export function readQcWorkflowState(kind: string) {
  const state = readWorkflowState("qc.pipeline.v1", kind);
  const owner = readPersistedUiState<{ kind?: string }>("qc.pipeline.v1", {}).kind ?? "protein";
  if (state.advancedV2 || owner !== kind) return state;
  // Early viewer sessions stored appearance separately from the active workflow's controls.
  return { ...state, advancedV2: readPersistedUiState("qc.advanced.v2", {}) } as Record<string, unknown>;
}
