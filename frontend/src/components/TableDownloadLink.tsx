import { collectColumns } from "../lib/tableData";
import { rowsToXlsx, XLSX_MIME } from "../lib/xlsx";

export type DownloadableTable = { rows: Record<string, unknown>[]; columns?: string[] };

/** A table as an Excel file. The workbook is written when the link is clicked, not while the page renders. */
export default function TableDownloadLink({ table, filename, label = "Download table (.xlsx)" }: {
  table: DownloadableTable;
  filename: string;
  label?: string;
}) {
  function prepare(link: HTMLAnchorElement) {
    const columns = table.columns ?? collectColumns(table.rows);
    const url = URL.createObjectURL(new Blob([rowsToXlsx(table.rows, columns) as BlobPart], { type: XLSX_MIME }));
    link.href = url;
    setTimeout(() => URL.revokeObjectURL(url), 0);
  }

  return (
    <a href="#" download={filename} onClick={(event) => prepare(event.currentTarget)}
      className="rounded-xl border border-slate-300 bg-white px-3 py-2 text-sm text-slate-700 transition hover:bg-slate-50">
      {label}
    </a>
  );
}
