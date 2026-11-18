import type { ReactNode } from "react";

export default function PlotOptionsContent({ loading, error, renderOptions }: {
  loading: boolean;
  error: string | null;
  renderOptions: () => ReactNode;
}) {
  return loading ? (
    <div className="mt-4 rounded-xl border border-slate-200 bg-slate-50 px-4 py-6 text-sm text-slate-500">
      Loading options...
    </div>
  ) : error ? (
    <div className="mt-4 rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-sm text-red-700">
      {error}
    </div>
  ) : (
    <div className="mt-4">{renderOptions()}</div>
  );
}
