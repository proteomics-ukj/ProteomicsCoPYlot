import { useEffect, useId, useState } from "react";
import NumberField from "./NumberField";
import CheckboxField from "./CheckboxField";
import { AssistantTextarea } from "./AssistantFields";

export default function HighlightTermsField({
  label,
  value,
  onChange,
  rows = 3,
  placeholder,
}: {
  label: string;
  value: string;
  onChange: (value: string) => void;
  rows?: number;
  placeholder?: string;
}) {
  const fieldId = useId();
  const [draft, setDraft] = useState(value);

  useEffect(() => {
    setDraft(value);
  }, [value]);

  const commit = () => {
    onChange(draft);
  };

  return (
    <div>
      <label htmlFor={fieldId} className="mb-2 block text-sm font-medium text-slate-700">{label}</label>
      <AssistantTextarea id={fieldId} data-assistant-permission="view"
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={(e) => {
          if (e.key === "Enter" && !e.shiftKey) {
            e.preventDefault();
            commit();
          }
        }}
        rows={rows}
        className="w-full rounded-xl border border-slate-300 bg-white px-3 py-2 text-sm outline-none focus:border-slate-900"
        placeholder={placeholder}
      />
    </div>
  );
}


export function HighlightControls({ terms, onTerms, topFirst, topSecond, firstLabel, secondLabel,
  onFirst, onSecond, outline, onOutline, labels, onLabels, interactive }: {
  terms: string; onTerms: (value: string) => void;
  topFirst: number; topSecond: number; firstLabel: string; secondLabel: string;
  onFirst: (value: number) => void; onSecond: (value: number) => void;
  outline: boolean; onOutline: (value: boolean) => void;
  labels: boolean; onLabels: (value: boolean) => void; interactive: boolean;
}) {
  return <div className="max-w-3xl space-y-4" data-assistant-appearance="true">
    <div className="grid grid-cols-2 gap-4">
      <NumberField label={firstLabel} value={topFirst} step={1}
        onChange={(value) => onFirst(Math.max(0, Math.round(value)))} />
      <NumberField label={secondLabel} value={topSecond} step={1}
        onChange={(value) => onSecond(Math.max(0, Math.round(value)))} />
    </div>
    <HighlightTermsField label="Highlight terms" value={terms} onChange={onTerms}
      placeholder="Space-separated labels to highlight" />
    {(terms.trim() || topFirst > 0 || topSecond > 0) ?
      <div className="grid gap-4 sm:grid-cols-2">
        <CheckboxField label="Outline highlighted points" checked={outline} onChange={onOutline} />
        {!interactive ? <CheckboxField label="Label highlighted points" checked={labels} onChange={onLabels} /> : null}
      </div> : null}
  </div>;
}
