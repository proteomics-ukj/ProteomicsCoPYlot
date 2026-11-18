export default function CsvDownloadLink({ href, filename }: { href: string; filename: string }) {
  return (
    <a href={href} download={filename}
      className="rounded-xl border border-slate-300 bg-white px-3 py-2 text-sm text-slate-700 transition hover:bg-slate-50">
      Download CSV
    </a>
  );
}
