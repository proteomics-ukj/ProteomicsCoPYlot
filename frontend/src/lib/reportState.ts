import type {
  AnnotationKind,
  AdvancedPlotSettings,
  CurrentDatasetsResponse,
  ReportEntry,
  ReportEntryType,
  SummaryReportContext,
} from "./types";
import { restoreVolcanoRenderType } from "./types";
import { getActiveAdvancedPlotSettingsJson, getAdvancedPlotMode, getAdvancedPlotSettingsByKey } from "./appSettings";
import { readPersistedUiState, readWorkflowState, readQcWorkflowState } from "./persistedUiState";
import { IS_VIEWER_MODE } from "./appMode";

const QC_SETTINGS_KEY = "pcopylot.report.qcSettings.v1";
const VOLCANO_LOG_KEY = "pcopylot.report.volcanoLog.v1";
const ENTRIES_KEY = "pcopylot.report.entries.v2";
const REPORT_SCOPE_KEY = "pcopylot.report.datasetScope.v1";

function readJson<T>(key: string, fallback: T): T {
  try {
    const raw = globalThis.sessionStorage?.getItem(key);
    if (!raw) return fallback;
    const parsed = JSON.parse(raw) as T;
    return parsed ?? fallback;
  } catch {
    return fallback;
  }
}

function writeJson(key: string, value: unknown) {
  try {
    globalThis.sessionStorage?.setItem(key, JSON.stringify(value));
  } catch {
    // Ignore storage write failures (private mode, quota, etc.).
  }
}

function readRaw(key: string): string | null {
  try {
    return globalThis.sessionStorage?.getItem(key) ?? null;
  } catch {
    return null;
  }
}

function writeRaw(key: string, value: string) {
  try {
    globalThis.sessionStorage?.setItem(key, value);
  } catch {
    // Ignore storage write failures (private mode, quota, etc.).
  }
}

function removeKey(key: string) {
  try {
    globalThis.sessionStorage?.removeItem(key);
  } catch {
    // Ignore storage write failures (private mode, quota, etc.).
  }
}

function datasetScopeToken(
  datasets: CurrentDatasetsResponse | null
): string {
  if (!datasets) return "none";
  const pkg = `package:${datasets.activePackage || "Data"}`;
  const metadataProfile = `metadata:${datasets.activeMetadataProfile || "Metadata"}`;
  const protein = datasets.protein
    ? `protein:${datasets.protein.filename}:${datasets.protein.rows}:${datasets.protein.columns}`
    : "protein:none";
  const phospho = datasets.phospho
    ? `phospho:${datasets.phospho.filename}:${datasets.phospho.rows}:${datasets.phospho.columns}`
    : "phospho:none";
  const phosprot = datasets.phosprot
    ? `phosprot:${datasets.phosprot.filename}:${datasets.phosprot.rows}:${datasets.phosprot.columns}`
    : "phosprot:none";
  // The file name only: its path on the server changes when a saved session is
  // restored, and a restore must keep the report list it brings back.
  const peptide = datasets.peptide
    ? `peptide:${datasets.peptide.filename}`
    : "peptide:none";
  return [pkg, metadataProfile, protein, phospho, phosprot, peptide].join("|");
}

function clearSummaryReportState() {
  removeKey(VOLCANO_LOG_KEY);
  removeKey(ENTRIES_KEY);
  removeKey(QC_SETTINGS_KEY);
}

export function syncSummaryReportStateWithDatasets(
  datasets: CurrentDatasetsResponse | null
) {
  const token = datasetScopeToken(datasets);
  const previous = readRaw(REPORT_SCOPE_KEY);
  if (previous == null) {
    clearSummaryReportState();
    writeRaw(REPORT_SCOPE_KEY, token);
    return;
  }
  if (previous === token) return;
  clearSummaryReportState();
  writeRaw(REPORT_SCOPE_KEY, token);
}

export function saveQcReportSettings(
  kind: AnnotationKind | "peptide" | "precursor",
  settings: Record<string, unknown>
) {
  const all = readJson<Record<string, Record<string, unknown>>>(
    QC_SETTINGS_KEY,
    {}
  );
  all[kind] = { ...all[kind], ...settings };
  writeJson(QC_SETTINGS_KEY, all);
}

function readAdvancedUiState(key: string): Record<string, AdvancedPlotSettings> {
  const state = readPersistedUiState<Record<string, unknown>>(key, {});
  // `advancedV2` since the panel stopped storing its displayed defaults; the
  // QC and Completeness pages write nothing under the old name any more.
  const advanced = state.advancedV2;
  return advanced && typeof advanced === "object" && !Array.isArray(advanced)
    ? (advanced as Record<string, AdvancedPlotSettings>)
    : {};
}

// What decides which plot an entry is, per type. Remove matches on these and
// on Normal/Interactive; highlight options, aliases, colours, the figure size
// and the advanced appearance do not make a different plot.
const SELECTION_KEYS: Record<ReportEntryType, readonly string[]> = {
  volcano: ["dataSource", "condition1", "condition2", "identifier", "testType", "pValueThreshold",
    "log2fcThreshold", "useUncorrected", "proteinCorrection"],
  pearson: ["mode", "sample1", "sample2", "condition1", "condition2", "header"],
  venn: ["mode", "first", "second", "third", "minPresent", "header"],
  anova: ["source", "dataSource", "identifier", "pValueThreshold", "useUncorrected", "showDendrograms",
    "showLegend", "showRowLabels", "nClusters", "clusterMethod", "selectedCluster", "postHoc",
    "postHocCondition1", "postHocCondition2"],
  ora: ["direction", "source", "dataSource", "condition1", "condition2", "condition1Control", "condition2Control",
    "pValueThreshold", "log2fcThreshold", "testType", "useUncorrected", "selectedCluster", "nClusters",
    "clusterMethod", "adjPValueCutoff", "topN", "minTermSize", "maxTermSize", "backgroundMode", "organism"],
  single_box: ["protein", "conditions", "identifier", "dataSource", "outliers", "dots", "sampleN", "randomSeed",
    "header", "legend"],
  single_line: ["proteins", "conditions", "includeId", "showNaLow", "identifier", "dataSource", "header", "legend"],
};

function selectionOf(entry: ReportEntry): string {
  return JSON.stringify([
    entry.type,
    entry.kind,
    entry.render.type,
    SELECTION_KEYS[entry.type].map((key) => entry.request[key] ?? null),
  ]);
}

function sortedJson(value: unknown): string {
  return JSON.stringify(value, (_key, item) =>
    item && typeof item === "object" && !Array.isArray(item)
      ? Object.fromEntries(Object.entries(item as Record<string, unknown>).sort(([a], [b]) => a.localeCompare(b)))
      : item
  );
}

// Volcano entries saved before the other plot types could be added.
function convertVolcanoLog(): ReportEntry[] {
  const saved = readJson<Record<string, unknown>[]>(VOLCANO_LOG_KEY, []);
  return saved
    .filter((entry) => !entry.control && (entry.kind === "protein" || entry.kind === "phospho"))
    .map((entry) => {
      const request = { ...entry };
      for (const key of ["kind", "control", "renderType", "condition1Control", "condition2Control"]) delete request[key];
      return {
        type: "volcano",
        kind: entry.kind as ReportEntry["kind"],
        request,
        render: { type: restoreVolcanoRenderType(entry.renderType) },
        advanced: null,
      };
    });
}

export function listReportEntries(): ReportEntry[] {
  const entries = readJson<ReportEntry[]>(ENTRIES_KEY, []);
  if (readRaw(VOLCANO_LOG_KEY) == null) return entries;
  const merged = [...convertVolcanoLog(), ...entries];
  writeJson(ENTRIES_KEY, merged);
  removeKey(VOLCANO_LOG_KEY);
  return merged;
}

/** `only` narrows a plot type to one of its views (the single clusters of an ANOVA, say). */
export function countReportEntries(type: ReportEntryType, kind: string, only?: (entry: ReportEntry) => boolean): number {
  return listReportEntries().filter((entry) => entry.type === type && entry.kind === kind && (!only || only(entry))).length;
}

/** The advanced appearance of the plot on screen, as an entry stores it. */
export function currentAdvancedAppearance(): AdvancedPlotSettings | null {
  try {
    const parsed = JSON.parse(getActiveAdvancedPlotSettingsJson() || "null");
    return parsed && typeof parsed === "object" && !Array.isArray(parsed) ? (parsed as AdvancedPlotSettings) : null;
  } catch {
    return null;
  }
}

export function addReportEntry(entry: ReportEntry): boolean {
  const current = listReportEntries();
  const identity = sortedJson(entry);
  if (current.some((existing) => sortedJson(existing) === identity)) return false;
  writeJson(ENTRIES_KEY, [...current, entry]);
  return true;
}

export function removeMatchingReportEntries(entry: ReportEntry): number {
  const current = listReportEntries();
  const selection = selectionOf(entry);
  const next = current.filter((existing) => selectionOf(existing) !== selection);
  writeJson(ENTRIES_KEY, next);
  return current.length - next.length;
}

export function deleteAllReportEntries(type: ReportEntryType, kind: string, only?: (entry: ReportEntry) => boolean): number {
  const current = listReportEntries();
  const next = current.filter((entry) => !(entry.type === type && entry.kind === kind && (!only || only(entry))));
  writeJson(ENTRIES_KEY, next);
  return current.length - next.length;
}

/** Whether the report holds a plot of this type with the same values for `keys`. */
export function hasReportEntryMatching(
  type: ReportEntryType,
  kind: string,
  request: Record<string, unknown>,
  keys: readonly string[]
): boolean {
  return listReportEntries().some(
    (entry) =>
      entry.type === type && entry.kind === kind && keys.every((key) => (entry.request[key] ?? null) === (request[key] ?? null))
  );
}

export function getSummaryReportContext(): SummaryReportContext {
  const qcSettings = readJson<Record<string, Record<string, unknown>>>(
    QC_SETTINGS_KEY,
    {}
  );
  const workflowPlotSettings: NonNullable<SummaryReportContext["workflowPlotSettings"]> = {};
  for (const kind of ["protein", "phospho", "phosprot", "peptide", "precursor"]) {
    const qc = readQcWorkflowState(kind);
    const completeness = readWorkflowState("completeness.v1", kind);
    workflowPlotSettings[kind] = {
      qcAppearance: (qc.advancedV2 ?? {}) as Record<string, AdvancedPlotSettings>,
      completenessAppearance: (completeness.advancedV2 ?? {}) as Record<string, AdvancedPlotSettings>,
      missingPlot: (completeness.missingPlot ?? {}) as Record<string, unknown>,
    };
  }
  return {
    qcSettings,
    workflowPlotSettings,
    entries: listReportEntries(),
    completenessSettings: readPersistedUiState<{ missingPlot?: Record<string, unknown> }>("completeness.v1", {}).missingPlot ?? {},
    advancedPlotMode: getAdvancedPlotMode(),
    advancedPlotSettingsByKey: getAdvancedPlotSettingsByKey(),
    // In viewer mode the QC page writes its advanced settings to their own
    // session key instead of into the persisted page state, so reading only
    // "qc.pipeline.v1" produced a report whose QC plots carried none of the
    // styling -- renamed conditions included -- that the viewer was looking at.
    qcAdvancedSettingsByTab: IS_VIEWER_MODE
      ? readPersistedUiState<Record<string, AdvancedPlotSettings>>("qc.advanced.v2", {})
      : readAdvancedUiState("qc.pipeline.v1"),
    completenessAdvancedSettingsByTab: readAdvancedUiState("completeness.v1"),
  };
}
