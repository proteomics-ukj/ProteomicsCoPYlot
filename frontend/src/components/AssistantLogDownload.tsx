import { useState, useSyncExternalStore } from "react";
import type { AssistantLog } from "../lib/assistantLog";

export default function AssistantLogDownload({ log }: { log: AssistantLog }) {
  const count = useSyncExternalStore(log.subscribe, log.count);
  const [open, setOpen] = useState(false);
  const download = () => {
    const url = URL.createObjectURL(new Blob([log.text()], { type: "text/plain;charset=utf-8" }));
    const link = document.createElement("a");
    link.href = url;
    link.download = `copylot-ai-conversation-${log.started.replace(/[:.]/g, "-")}.txt`;
    link.click();
    window.setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
  return <details onToggle={event => setOpen(event.currentTarget.open)} className="border-b border-slate-200 px-4 py-2 text-xs">
    <summary className="cursor-pointer">Conversation log ({count} events)</summary>
    <p className="my-2">Review your messages, AI responses, function calls and changes before sharing. The log may contain research data. It stays in this browser tab until you download it; closing or reloading the page clears it.</p>
    <button type="button" onClick={download} disabled={!count} className="rounded border px-3 py-1 disabled:opacity-50">Download conversation log</button>
    {open && <pre aria-label="Conversation log preview" className="mt-2 max-h-64 overflow-auto whitespace-pre-wrap break-words">{log.text()}</pre>}
  </details>;
}
