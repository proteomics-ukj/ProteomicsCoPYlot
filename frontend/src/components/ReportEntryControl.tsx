import { useState } from "react";
import { IS_VIEWER_MODE } from "../lib/appMode";
import {
  addReportEntry,
  countReportEntries,
  deleteAllReportEntries,
  removeMatchingReportEntries,
} from "../lib/reportState";
import type { ReportEntry, ReportEntryType } from "../lib/types";

type Props = {
  type: ReportEntryType;
  kind: string;
  kindLabel: string;
  /** Button wording, e.g. "Volcano". */
  title: string;
  /** Message wording, e.g. "volcano plot". */
  noun: string;
  /** The plots on screen as entries; none while they cannot be added. */
  buildEntries: () => ReportEntry[];
  invalidMessage: string;
  /** A notice about the saved entries, read whenever the list changes. */
  notice?: () => string | null;
  /** The view of this plot type the control is for; its count and Delete All cover only that view. */
  only?: (entry: ReportEntry) => boolean;
};

export default function ReportEntryControl({
  type,
  kind,
  kindLabel,
  title,
  noun,
  buildEntries,
  invalidMessage,
  notice,
  only,
}: Props) {
  const [version, setVersion] = useState(0);
  const [message, setMessage] = useState<string | null>(null);

  // The report has a protein and a phosphosite level only.
  if (IS_VIEWER_MODE || (kind !== "protein" && kind !== "phospho")) return null;

  void version;
  const count = countReportEntries(type, kind, only);
  const noticeText = notice ? notice() : null;
  const removedEntries = (removed: number) => `${noun} entr${removed === 1 ? "y" : "ies"}`;
  const sentence = noun.charAt(0).toUpperCase() + noun.slice(1);

  function handleAdd() {
    const entries = buildEntries();
    if (!entries.length) {
      setMessage(invalidMessage);
      return;
    }
    const added = entries.map((entry) => addReportEntry(entry)).some(Boolean);
    setVersion((value) => value + 1);
    setMessage(added ? `${sentence} added to report list.` : `This ${noun} configuration is already in the report list.`);
  }

  function handleRemove() {
    const removed = buildEntries().reduce((total, entry) => total + removeMatchingReportEntries(entry), 0);
    setVersion((value) => value + 1);
    setMessage(
      removed > 0
        ? `Removed ${removed} ${removedEntries(removed)} from report list.`
        : `No matching ${noun} entry found in report list.`
    );
  }

  function handleDeleteAll() {
    const removed = deleteAllReportEntries(type, kind, only);
    setVersion((value) => value + 1);
    setMessage(
      removed > 0
        ? `Removed ${removed} ${removedEntries(removed)} from report list.`
        : `No ${noun} entry found in report list.`
    );
  }

  return (
    <div className="mt-4 rounded-xl border border-slate-200 bg-slate-50 p-4">
      <div className="flex flex-wrap items-center gap-2">
        <button
          type="button"
          onClick={handleAdd}
          className="rounded-xl bg-slate-900 px-3 py-2 text-sm font-medium text-white transition hover:bg-slate-800"
        >
          Add {title} To Report
        </button>
        <button
          type="button"
          onClick={handleRemove}
          className="rounded-xl border border-slate-300 bg-white px-3 py-2 text-sm text-slate-700 transition hover:bg-slate-50"
        >
          Remove {title} From Report
        </button>
        <button
          type="button"
          onClick={handleDeleteAll}
          className="rounded-xl border border-slate-300 bg-white px-3 py-2 text-sm text-slate-700 transition hover:bg-slate-50"
        >
          Delete All From Report
        </button>
        <span className="text-sm text-slate-600">
          Saved entries for {kindLabel}: {count}
        </span>
      </div>
      {noticeText ? (
        <div role="status" className="mt-3 rounded-xl border border-amber-200 bg-amber-50 px-3 py-2 text-sm text-amber-900">
          {noticeText}
        </div>
      ) : null}
      {message ? <div className="mt-3 text-sm text-slate-700">{message}</div> : null}
    </div>
  );
}
