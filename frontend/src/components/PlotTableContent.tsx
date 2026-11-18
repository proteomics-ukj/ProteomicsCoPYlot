import type { ReactNode } from "react";

export default function PlotTableContent({ loading, error, renderTable }: {
  loading: boolean;
  error: string | null;
  renderTable: () => ReactNode;
}) {
  return loading ? (
    <div className="rounded-xl border border-slate-200 bg-slate-50 px-4 py-6 text-sm text-slate-500">
      Loading table...
    </div>
  ) : error ? (
    <div className="rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-sm text-red-700">
      {error}
    </div>
  ) : renderTable();
}
