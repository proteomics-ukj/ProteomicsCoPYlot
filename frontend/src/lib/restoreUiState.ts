const ANALYSIS_STORAGE_PREFIXES = ["pcopylot.report.", "pcopylot.ui.", "copylot:advanced"] as const;
// Saved plot styles travel in a save state. Unlike the analysis keys they are
// not cleared when the session is reset.
const STYLE_STORAGE_PREFIXES = ["copylot:savedPlotStyles.", "copylot:plotStyleTemplate.", "copylot:plotStyleSeen."] as const;

function hasPrefix(key: string, prefixes: readonly string[]): boolean {
  return prefixes.some((prefix) => key.startsWith(prefix));
}

function removeKeys(storage: Storage, prefixes: readonly string[]): void {
  const keysToRemove: string[] = [];
  for (let index = 0; index < storage.length; index += 1) {
    const key = storage.key(index);
    if (key && hasPrefix(key, prefixes)) keysToRemove.push(key);
  }
  for (const key of keysToRemove) storage.removeItem(key);
}

function removeCopylotKeys(storage: Storage): void {
  removeKeys(storage, ANALYSIS_STORAGE_PREFIXES);
}

export function clearCopylotClientStorageState(): void {
  try {
    if (typeof globalThis.sessionStorage !== "undefined") {
      removeCopylotKeys(globalThis.sessionStorage);
    }
  } catch {}
}

export function captureCopylotSessionStorageState(): Record<string, string | null> {
  try {
    if (typeof globalThis.sessionStorage === "undefined") return {};
    const payload: Record<string, string | null> = {};
    for (let index = 0; index < globalThis.sessionStorage.length; index += 1) {
      const key = globalThis.sessionStorage.key(index);
      if (!key) continue;
      if (!hasPrefix(key, ANALYSIS_STORAGE_PREFIXES) && !hasPrefix(key, STYLE_STORAGE_PREFIXES)) continue;
      payload[key] = globalThis.sessionStorage.getItem(key);
    }
    return payload;
  } catch {
    return {};
  }
}

export function applyCopylotSessionStorageState(
  payload: Record<string, unknown> | null | undefined,
  replaceStyles = false
): void {
  try {
    if (typeof globalThis.sessionStorage === "undefined") return;
    const entries = payload && typeof payload === "object" ? payload : {};
    removeCopylotKeys(globalThis.sessionStorage);
    // A save state that records styles replaces the current ones; an older one
    // does not know them and leaves them alone.
    if (replaceStyles) removeKeys(globalThis.sessionStorage, STYLE_STORAGE_PREFIXES);
    for (const [key, raw] of Object.entries(entries)) {
      const restorable = hasPrefix(key, ANALYSIS_STORAGE_PREFIXES) || (replaceStyles && hasPrefix(key, STYLE_STORAGE_PREFIXES));
      if (restorable && typeof raw === "string") globalThis.sessionStorage.setItem(key, raw);
    }
  } catch {}
}
