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
      <p>Cambridge cRAP: all categories are selected by default. Choose any combination below.
        Rows are removed after the valid-value filter only when every protein in the group matches.
        Mixed groups and missing identifiers stay. Gene names are never matched.</p>
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
        <p className="mt-2">Uses protein entry names or full protein names, with accessions as an additional
          match or fallback. MaxQuant CON__ and FragPipe contam_ prefixes are ignored for matching;
          search-engine flags alone do not remove rows. No partial matches or automatic isoform merging.</p>
      </details>
      <p>Protein settings also apply to peptide parent proteins. Phosphosite settings are independent.
        Review the selected categories: some listed proteins may be genuine components of your sample.</p>
    </div>}
  </div>;
}
