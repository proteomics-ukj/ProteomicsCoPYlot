import { useWorkflow, withWorkflowSettings } from "../../lib/workflowSettings";
// Volcano panel: ordinary, paired and control-adjusted comparisons.

import { SizeRow } from "../../components/PlotFields";
import { restorePlotDimensionCm, restorePlotDpi } from "../../lib/plotLimits";
import { coercePlotDownloadFormat, type PlotDownloadFormat } from "../../lib/plotDownload";
import { HighlightControls } from "../../components/HighlightTermsField";
import { useEffect, useMemo, useState, type ReactNode } from "react";
import PairingEditor from "../../components/PairingEditor";
import { buildPlotUrl, getStatisticalOptions, getVolcanoNotebookTable, isAbortError, pairingIsReady, runVolcanoAnalysis, runVolcanoControlAnalysis } from "../../lib/api";
import type { RenderedNotebook } from "../../components/NotebookDownload";
import { useDatasetRevision } from "../../lib/datasetAvailability";
import ReportEntryControl from "../../components/ReportEntryControl";
import { currentAdvancedAppearance } from "../../lib/reportState";
import { useAdvancedPlotSignal } from "../../lib/useAdvancedPlotSignal";
import { useDebouncedValue } from "../../lib/useDebouncedValue";
import { readWorkflowState, writeWorkflowState } from "../../lib/persistedUiState";
import type { ContextMode, AnnotationKind, StatsIdentifier, StatsTestType, VolcanoResultResponse, ReportEntry, DataSource, VolcanoRenderType, StatisticalOptionsResponse } from "../../lib/types";
import { VOLCANO_RENDER_TYPE_OPTIONS } from "../../lib/types";
import { PROTEIN_CORRECTION_LABEL, useProteinCorrection } from "../../lib/proteinCorrection";
import { useDataSourceChoice } from "../../lib/useDataSourceChoice";
import { restoreTestType, CONTROL_TEST_TYPE_OPTIONS, controlTestType, useTestTypeOptions, availableTestType } from "../../lib/statisticalTests";
import SelectField from "../../components/SelectField";
import SectionCard from "../../components/SectionCard";
import CheckboxField from "../../components/CheckboxField";
import NumberField from "../../components/NumberField";
import { keepCondition, keepIdentifier, useStatisticalOptionsLoader } from "./useStatisticalOptions";
import { ImageSection, InfoSection, Notice, PlotFrame, SummarySection, TableSection, ConditionPicker } from "./statsShared";
import { conditionOptions, isPairedTestType, removeSourceUsageWarnings, restorePValueThreshold, sanitizeFilenamePart } from "./statsHelpers";
import LocalizationCutoffField from "../../components/LocalizationCutoffField";


function volcanoExportFilename(
  sourceMode: "volcano" | "volcano_control",
  kind: string,
  condition1: string,
  condition2: string
): string {
  const prefix = sourceMode === "volcano_control" ? "volcano_control" : "volcano";
  const first = sanitizeFilenamePart(condition1);
  const second = sanitizeFilenamePart(condition2);
  if (!first || !second) return `${prefix}_${kind}.csv`;
  return `${prefix}_${kind}_${first}_vs_${second}.csv`;
}

function VolcanoPanelView({
  kindOptions,
  onContextModeChange,
  advancedPlotOptions,
}: {
  kindOptions: { value: string; label: string }[];
  onContextModeChange?: (mode: ContextMode | null) => void;
  advancedPlotOptions?: ReactNode;
}) {
  const advancedPlotSignal = useAdvancedPlotSignal();
  const [kind, setKind] = useWorkflow<AnnotationKind>();
  const persisted = useMemo(() => readWorkflowState("stats.volcano.v1", kind), [kind]);
  const revision = useDatasetRevision();
  const [sourceMode, setSourceMode] = useState<"volcano" | "volcano_control">(() => {
    const value = String(persisted.sourceMode ?? "volcano").trim();
    return value === "volcano_control" ? "volcano_control" : "volcano";
  });
  const {
    dataSource,
    applyChoices: applyDataSourceChoices,
    chooseDataSource,
  } = useDataSourceChoice(persisted.dataSource);
  const [options, setOptions] = useState<StatisticalOptionsResponse | null>(null);
  const [condition1, setCondition1] = useState(() => String(persisted.condition1 ?? ""));
  const [condition2, setCondition2] = useState(() => String(persisted.condition2 ?? ""));
  const [condition1Control, setCondition1Control] = useState(() => String(persisted.condition1Control ?? ""));
  const [condition2Control, setCondition2Control] = useState(() => String(persisted.condition2Control ?? ""));
  const [identifier, setIdentifier] = useState<StatsIdentifier>(() => {
    const value = String(persisted.identifier ?? "workflow").trim();
    return value === "genes" ? "genes" : "workflow";
  });
  const [selectedTestType, setTestType] = useState<StatsTestType>(() => {
    return restoreTestType(persisted.testType);
  });
  const { options: testTypeOptions, moderatedAvailable } = useTestTypeOptions();
  const testType = sourceMode === "volcano_control"
    ? controlTestType(selectedTestType)
    : availableTestType(selectedTestType, moderatedAvailable);
  const [pValueThreshold, setPValueThreshold] = useState(() =>
    restorePValueThreshold(persisted.pValueThreshold)
  );
  const [log2fcThreshold, setLog2fcThreshold] = useState(() => Number(persisted.log2fcThreshold ?? 1));
  const [useUncorrected, setUseUncorrected] = useState(() => Boolean(persisted.useUncorrected ?? false));
  const [highlightText, setHighlightText] = useState(() => String(persisted.highlightText ?? ""));
  // What the highlight terms draw: a ring (default) and/or, on the static
  // figure, a text label placed like the top hits.
  const [highlightOutline, setHighlightOutline] = useState(() => persisted.highlightOutline !== false);
  const [highlightLabels, setHighlightLabels] = useState(() => Boolean(persisted.highlightLabels ?? false));
  // Interactive by default, because this panel is for exploring and the hover
  // text is how you find out which point you are looking at. Normal is the
  // figure you take away.
  const [renderType, setRenderType] = useState<VolcanoRenderType>(() =>
    String(persisted.renderType ?? "Interactive").trim() === "Normal" ? "Normal" : "Interactive"
  );
  const [labelTopUp, setLabelTopUp] = useState(() => Number(persisted.labelTopUp ?? 0));
  const [labelTopDown, setLabelTopDown] = useState(() => Number(persisted.labelTopDown ?? 0));
  // Size and export belong to the static figure only. The interactive plot
  // fills its frame and exports itself through Plotly's camera button, so
  // neither is offered -- nor sent -- while it is selected.
  const [widthCm, setWidthCm] = useState(() => restorePlotDimensionCm(persisted.widthCm, 18));
  const [heightCm, setHeightCm] = useState(() => restorePlotDimensionCm(persisted.heightCm, 14));
  const [dpi, setDpi] = useState(() => restorePlotDpi(persisted.dpi, 300));
  const [downloadFormat, setDownloadFormat] = useState<PlotDownloadFormat>(() =>
    coercePlotDownloadFormat(persisted.downloadFormat)
  );
  const [error, setError] = useState<string | null>(null);
  const [renderedNotebook, setRenderedNotebook] = useState<RenderedNotebook | null>(null);
  const [notebookTable, setNotebookTable] = useState<{ url: string; token: string; data: VolcanoResultResponse } | null>(null);
  const [resultState, setResultState] = useState<{
    key: string;
    data: VolcanoResultResponse;
  } | null>(null);

  const kindAvailable = kindOptions.some((option) => option.value === kind);
  const proteinCorrection = useProteinCorrection(
    kind, options, persisted.proteinCorrection, dataSource, chooseDataSource);
  const volcanoDataSourceOptions = proteinCorrection.dataSourceOptions;

  useStatisticalOptionsLoader(
    kind,
    kindAvailable,
    (data) => {
      setOptions(data);
      applyDataSourceChoices(data);
      setCondition1(keepCondition(data.availableConditions));
      setCondition2(keepCondition(data.availableConditions));
      setCondition1Control(keepCondition(data.availableConditions));
      setCondition2Control(keepCondition(data.availableConditions));
      setIdentifier(keepIdentifier(data));
    },
    setError,
  );


  const controlSelectionValid =
    sourceMode !== "volcano_control" ||
    new Set([condition1, condition1Control, condition2, condition2Control].filter(Boolean)).size === 4;
  const canRun = Boolean(
    kindAvailable &&
      condition1 &&
      condition2 &&
      condition1 !== condition2 &&
      controlSelectionValid &&
      (!isPairedTestType(testType) ||
        pairingIsReady(
          condition1,
          condition2,
          options?.samplesByCondition ?? {},
          options?.pairIdsBySample ?? {}
        ))
  );
  const highlightTerms = useMemo(() => highlightText.split(/\s+/).map((term) => term.trim()).filter(Boolean), [highlightText]);
  const requestPayload = useMemo(() => {
    if (!canRun) return null;
    return {
      source: sourceMode,
      kind,
      dataSource,
      condition1,
      condition2,
      condition1Control,
      condition2Control,
      identifier,
      pValueThreshold,
      log2fcThreshold,
      testType,
      useUncorrected,
      proteinCorrection: proteinCorrection.active,
      highlightTerms,
    };
  }, [
    canRun,
    sourceMode,
    kind,
    dataSource,
    condition1,
    condition2,
    condition1Control,
    condition2Control,
    identifier,
    pValueThreshold,
    log2fcThreshold,
    testType,
    useUncorrected,
    proteinCorrection.active,
    highlightTerms,
  ]);

  useEffect(() => {
    if (!onContextModeChange) return;
    const mode: ContextMode = {
      source: sourceMode,
      kind: kind,
      dataSource: dataSource,
      condition1: condition1,
      condition2: condition2,
      condition1Control: condition1Control,
      condition2Control: condition2Control,
      identifier: identifier,
      testType: testType,
      pValueThreshold: pValueThreshold,
      log2fcThreshold: log2fcThreshold,
      useUncorrected: useUncorrected,
      proteinCorrection: proteinCorrection.active,
      type: renderType,
    };
    onContextModeChange(mode);
  }, [
    onContextModeChange,
    sourceMode,
    kind,
    dataSource,
    condition1,
    condition2,
    condition1Control,
    condition2Control,
    identifier,
    testType,
    pValueThreshold,
    log2fcThreshold,
    useUncorrected,
    proteinCorrection.active,
    renderType,
  ]);
  const debouncedRequestPayload = useDebouncedValue(requestPayload, 650);
  const currentRequestKey = requestPayload ? JSON.stringify({ revision, ...requestPayload }) : "";
  const result = currentRequestKey && resultState?.key === currentRequestKey ? resultState.data : null;

  // The plot as it is now, captured at the moment it is added, so two
  // comparisons in one report can render differently.
  function buildReportEntries(): ReportEntry[] {
    if (!canRun || (kind !== "protein" && kind !== "phospho")) return [];
    return [{
      type: "volcano",
      kind,
      request: {
        dataSource,
        condition1,
        condition2,
        identifier,
        pValueThreshold,
        log2fcThreshold,
        testType,
        useUncorrected,
        proteinCorrection: proteinCorrection.active,
        highlightTerms,
        labelTopUp,
        labelTopDown,
        highlightOutline,
        highlightLabels,
      },
      render: { type: renderType, widthCm, heightCm, dpi },
      advanced: currentAdvancedAppearance(),
    }];
  }

  useEffect(() => {
    writeWorkflowState("stats.volcano.v1", {
      kind,
      sourceMode,
      dataSource,
      condition1,
      condition2,
      condition1Control,
      condition2Control,
      identifier,
      testType,
      pValueThreshold,
      log2fcThreshold,
      useUncorrected,
      proteinCorrection: proteinCorrection.enabled,
      highlightText,
      highlightOutline,
      highlightLabels,
      renderType,
      labelTopUp,
      labelTopDown,
      widthCm,
      heightCm,
      dpi,
      downloadFormat,
    });
  }, [
    kind,
    sourceMode,
    dataSource,
    condition1,
    condition2,
    condition1Control,
    condition2Control,
    identifier,
    testType,
    pValueThreshold,
    log2fcThreshold,
    useUncorrected,
    proteinCorrection.enabled,
    highlightText,
    highlightOutline,
    highlightLabels,
    renderType,
    labelTopUp,
    labelTopDown,
    widthCm,
    heightCm,
    dpi,
    downloadFormat,
  ]);

  useEffect(() => {
    if (!debouncedRequestPayload) {
      return;
    }

    const controller = new AbortController();

    const run = async () => {
      try {
        const { source, condition1Control, condition2Control, ...request } = debouncedRequestPayload;
        const response = source === "volcano_control"
          ? await runVolcanoControlAnalysis({ ...request, condition1Control, condition2Control }, { signal: controller.signal })
          : await runVolcanoAnalysis(request, { signal: controller.signal });
        if (controller.signal.aborted) return;
        setResultState({ key: JSON.stringify({ revision, ...debouncedRequestPayload }), data: response });
        setError(null);
      } catch (err) {
        if (isAbortError(err)) return;
        setResultState(null);
        setError(err instanceof Error ? err.message : "Failed to run volcano analysis");
      }
    };

    run();

    return () => {
      controller.abort();
    };
  }, [
    debouncedRequestPayload,
    revision,
  ]);

  const plotUrl = useMemo(() => {
    void advancedPlotSignal;
    // Also a recompute trigger rather than a value: buildPlotUrl reads the
    // dataset-scope nonce, so a new revision has to rebuild this URL.
    void revision;
    if (!debouncedRequestPayload) return "";
    const isStatic = renderType === "Normal";
    const extension = isStatic ? "png" : "html";
    // The label counts only mean something to the static renderer; the
    // interactive one names a point when you hover it.
    const shared = {
      condition1: debouncedRequestPayload.condition1,
      condition2: debouncedRequestPayload.condition2,
      dataSource: debouncedRequestPayload.dataSource,
      identifier: debouncedRequestPayload.identifier,
      pValueThreshold: debouncedRequestPayload.pValueThreshold,
      log2fcThreshold: debouncedRequestPayload.log2fcThreshold,
      testType: debouncedRequestPayload.testType,
      useUncorrected: debouncedRequestPayload.useUncorrected,
      proteinCorrection: debouncedRequestPayload.proteinCorrection,
      highlightTerms: debouncedRequestPayload.highlightTerms.join(","),
      highlightOutline,
      labelTopUp, labelTopDown,
      ...(isStatic ? { highlightLabels, widthCm, heightCm, dpi } : {}),
    };
    if (debouncedRequestPayload.source === "volcano_control") {
      return buildPlotUrl(
        `/api/plots/stats/${debouncedRequestPayload.kind}/volcano-control.${extension}`,
        {
          ...shared,
          condition1Control: debouncedRequestPayload.condition1Control,
          condition2Control: debouncedRequestPayload.condition2Control,
        }
      );
    }
    return buildPlotUrl(
      `/api/plots/stats/${debouncedRequestPayload.kind}/volcano.${extension}`,
      shared
    );
  }, [
    debouncedRequestPayload,
    advancedPlotSignal,
    revision,
    renderType,
    labelTopUp,
    labelTopDown,
    highlightOutline,
    highlightLabels,
    widthCm,
    heightCm,
    dpi,
  ]);

  useEffect(() => {
    if (renderedNotebook?.url !== plotUrl || !renderedNotebook.token) return;
    const controller = new AbortController();
    const token = renderedNotebook.token;
    getVolcanoNotebookTable(kind, token, controller.signal).then(data => {
      if (!controller.signal.aborted) setNotebookTable({ url: plotUrl, token, data });
    }).catch(cause => {
      if (!controller.signal.aborted && !isAbortError(cause)) {
        setError(cause instanceof Error ? cause.message : "Failed to load the captured volcano table");
      }
    });
    return () => controller.abort();
  }, [kind, plotUrl, renderedNotebook]);
  const capturedTable = notebookTable?.url === plotUrl
    && notebookTable.token === renderedNotebook?.token ? notebookTable : null;
  const displayedResult = capturedTable?.data ?? result;

  return (
    <div className="space-y-6">
      <SectionCard title="Volcano Plot">
        <div className="grid max-w-2xl gap-4 sm:grid-cols-2">
          <SelectField label="Dataset level" value={kind} onChange={(value) => setKind(value as AnnotationKind)} options={kindOptions} />
          <LocalizationCutoffField kind={kind} />
        </div>
      </SectionCard>

      <SectionCard title="Options">
        <div className="grid gap-4 lg:grid-cols-2">
          <SelectField
            label="Source"
            value={sourceMode}
            onChange={(value) => setSourceMode(value as "volcano" | "volcano_control")}
            options={[
              { value: "volcano", label: "Volcano" },
              { value: "volcano_control", label: "Volcano Control" },
            ]}
          />
          <SelectField
            label="Data Source"
            value={dataSource}
            onChange={(value) => chooseDataSource(value as DataSource)}
            options={volcanoDataSourceOptions}
          />
        </div>
        <ConditionPicker
          conditions={conditionOptions(options)}
          condition1={{ value: condition1, onChange: setCondition1 }}
          condition2={{ value: condition2, onChange: setCondition2 }}
          controls={sourceMode === "volcano_control" ? {
            condition1: { value: condition1Control, onChange: setCondition1Control },
            condition2: { value: condition2Control, onChange: setCondition2Control },
          } : undefined}
        />

        <div className="mt-4 grid gap-4 sm:grid-cols-2">
          <SelectField
            label="Type"
            value={renderType}
            onChange={(value) => setRenderType(value as VolcanoRenderType)}
            options={VOLCANO_RENDER_TYPE_OPTIONS}
          />
          <SelectField label="Label column" value={identifier} onChange={(value) => setIdentifier(value as StatsIdentifier)} options={(options?.availableIdentifiers ?? []).map((item) => ({ value: item.key, label: item.label }))} />
        </div>
        <div className="mt-4 grid gap-4 sm:grid-cols-3">
          <SelectField label="Test type" value={testType} onChange={(value) => setTestType(value as StatsTestType)} options={sourceMode === "volcano_control" ? CONTROL_TEST_TYPE_OPTIONS : testTypeOptions} />
          <NumberField label="P-value threshold" value={pValueThreshold} onChange={setPValueThreshold} step={0.001} />
          <NumberField label="log2 FC threshold" value={log2fcThreshold} onChange={setLog2fcThreshold} step={0.1} />
        </div>
        {isPairedTestType(testType) ? (
          <PairingEditor
            kind={kind}
            condition1={condition1}
            condition2={condition2}
            samplesByCondition={options?.samplesByCondition ?? {}}
            pairIdsBySample={options?.pairIdsBySample ?? {}}
            onSaved={async () => setOptions(await getStatisticalOptions(kind))}
          />
        ) : null}

        <div className="mt-4">
          <div className="grid gap-4 sm:grid-cols-2 xl:grid-cols-3">
            <CheckboxField label="Use uncorrected p-values" checked={useUncorrected} onChange={setUseUncorrected} />
            {proteinCorrection.offered ? (
              <CheckboxField
                label={PROTEIN_CORRECTION_LABEL}
                checked={proteinCorrection.enabled}
                onChange={proteinCorrection.setEnabled}
              />
            ) : null}
          </div>
        </div>

        <div className="mt-4">
          <HighlightControls terms={highlightText} onTerms={setHighlightText}
            topFirst={labelTopUp} topSecond={labelTopDown}
            firstLabel="Highlight top N up-regulated" secondLabel="Highlight top N down-regulated"
            onFirst={setLabelTopUp} onSecond={setLabelTopDown}
            outline={highlightOutline} onOutline={setHighlightOutline}
            labels={highlightLabels} onLabels={setHighlightLabels} interactive={renderType === "Interactive"} />
        </div>
        {renderType === "Normal" ? (
          <div className="mt-4 max-w-3xl">
            <SizeRow
              widthCm={widthCm}
              heightCm={heightCm}
              dpi={dpi}
              onWidthChange={setWidthCm}
              onHeightChange={setHeightCm}
              onDpiChange={setDpi}
            />
          </div>
        ) : null}

        {sourceMode === "volcano" ? (
          <ReportEntryControl
            type="volcano"
            kind={kind}
            kindLabel={kindOptions.find((option) => option.value === kind)?.label ?? kind}
            title="Volcano"
            noun="volcano plot"
            buildEntries={buildReportEntries}
            invalidMessage="Please choose valid conditions before adding this volcano plot."
          />
        ) : null}

      <Notice error={error} warnings={removeSourceUsageWarnings(options?.warnings)} />
      </SectionCard>
      {advancedPlotOptions}
      {!canRun && options ? (
        <InfoSection message="Please select two different conditions to generate the plot." />
      ) : null}
      {canRun && result && displayedResult ? (
        <>
          <SummarySection items={[{ label: "Rows", value: String(displayedResult.totalRows) }, { label: "Upregulated", value: String(displayedResult.upregulatedCount) }, { label: "Downregulated", value: String(displayedResult.downregulatedCount) }, { label: "Not Significant", value: String(displayedResult.notSignificantCount) }]} />
          {renderType === "Normal" ? (
            <ImageSection
              title="Volcano Plot"
              url={plotUrl}
              notebook={renderedNotebook}
              onRendered={(url, token) => setRenderedNotebook({ url, token })}
              download={{
                filename: volcanoExportFilename(sourceMode, kind, condition1, condition2).replace(/\.csv$/, ".png"),
                format: downloadFormat,
                onFormatChange: setDownloadFormat,
              }}
            />
          ) : (
            <PlotFrame title="Volcano Plot" url={plotUrl} height={760} notebook={renderedNotebook}
              onRendered={(url, token) => setRenderedNotebook({ url, token })} />
          )}
          <Notice error={null} warnings={removeSourceUsageWarnings(Array.from(new Set([...(displayedResult.warnings ?? []), ...(result.warnings ?? [])])))} />
          <TableSection
            title="Result Table"
            rows={displayedResult.rows}
            plotUrl={plotUrl}
            notebook={capturedTable ? { url: capturedTable.url, token: capturedTable.token } : null}
            filename={volcanoExportFilename(sourceMode, kind, condition1, condition2)}
          />
          {displayedResult.untestedRows && displayedResult.untestedRows.length > 0 ? (
            <TableSection
              title="Rows not tested"
              rows={displayedResult.untestedRows}
              plotUrl={plotUrl}
              notebook={capturedTable ? { url: capturedTable.url, token: capturedTable.token } : null}
              filename={volcanoExportFilename(sourceMode, kind, condition1, condition2).replace(/\.csv$/, "_not_tested.csv")}
            />
          ) : null}
        </>
      ) : null}
    </div>
  );
}

export const VolcanoPanel = withWorkflowSettings("stats.volcano.v1", VolcanoPanelView, props => props.kindOptions);
