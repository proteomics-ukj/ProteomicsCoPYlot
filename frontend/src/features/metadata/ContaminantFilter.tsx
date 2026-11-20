import type { ContaminantFilterConfig } from "../../lib/types";
import { CONTAMINANT_CATEGORIES, DEFAULT_CONTAMINANTS } from "./contaminantOptions";

export default function ContaminantFilter({ value = DEFAULT_CONTAMINANTS, onChange }: {
  value?: ContaminantFilterConfig;
  onChange: (value: ContaminantFilterConfig) => void;
}) {
  return <div className="mt-5 rounded-xl border border-slate-200 p-4">
    <label className="flex items-center gap-2 font-medium text-slate-900">
      <input type="checkbox" checked={value.enabled}
        onChange={(e) => onChange({ ...value, enabled: e.target.checked })} />
      Filter contaminants
    </label>
    {value.enabled && <div className="mt-3 space-y-3 text-sm text-slate-700">
      <fieldset className="space-y-2">
        <legend className="mb-2 font-medium">Contaminant categories</legend>
        {CONTAMINANT_CATEGORIES.map(([id, label]) => <label key={id} className="flex items-center gap-2">
          <input type="checkbox" checked={value.categories.includes(id)}
            disabled={value.categories.length === 1 && value.categories.includes(id)}
            onChange={(e) => onChange({ ...value, categories: e.target.checked
              ? [...value.categories, id].sort() : value.categories.filter((category) => category !== id) })} />
          {label}
        </label>)}
      </fieldset>
      <details>
        <summary className="cursor-pointer font-medium">Advanced matching</summary>
        <label className="mt-2 flex items-center gap-2">
          <input type="checkbox" checked={value.matchAccessions}
            onChange={(e) => onChange({ ...value, matchAccessions: e.target.checked })} />
          Also match exact protein accessions when available
        </label>
      </details>
    </div>}
  </div>;
}
