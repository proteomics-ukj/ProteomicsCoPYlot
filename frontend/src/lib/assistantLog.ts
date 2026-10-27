export type LogEvent = { time: string; kind: string; data: unknown };

const excluded = /api.?key|authorization|password|secret|credential|signature/i;

/** Session-only record. No browser storage or automatic network upload. */
export class AssistantLog {
  private events: LogEvent[] = [];
  private secrets = new Set<string>();
  private listeners = new Set<() => void>();
  readonly started = new Date().toISOString();
  subscribe = (listener: () => void) => { this.listeners.add(listener); return () => { this.listeners.delete(listener); }; };
  count = () => this.events.length;

  exclude(...values: (string | undefined)[]) {
    values.filter(Boolean).forEach(value => this.secrets.add(value!));
    // Also remove a credential if it appeared in an earlier user message.
    this.events = this.events.map(event => ({ ...event, data: this.clean(event.data) }));
  }

  private clean(value: unknown): unknown {
    if (typeof value === "string") {
      let text = value;
      for (const secret of this.secrets) text = text.split(secret).join("[redacted]");
      return text.replace(/data:image\/[^\s"']+/g, "[image excluded]");
    }
    if (Array.isArray(value)) return value.map(item => this.clean(item));
    if (value && typeof value === "object") return Object.fromEntries(Object.entries(value).map(([key, item]) =>
      [key, excluded.test(key) || ["image", "plotImage", "ollamaUrl", "localServerUrl"].includes(key) ? "[excluded]" : this.clean(item)]));
    return value;
  }

  record = (kind: string, data: unknown) => {
    const cleaned = this.clean(data);
    // The progress line is re-rendered every 1.5 s while a model thinks; a
    // repeat of the previous event is a tick, not a change, and is not logged.
    const last = this.events[this.events.length - 1];
    if (kind === "Displayed assistant response / progress" && last && last.kind === kind
      && JSON.stringify(last.data) === JSON.stringify(cleaned)) return;
    this.events.push({ time: new Date().toISOString(), kind, data: cleaned });
    this.listeners.forEach(listener => listener());
  };

  text() {
    return ["CoPYlot AI conversation log — format 1", `CoPYlot frontend version: ${version}`, `Session started: ${this.started}`,
      "Times are UTC. This log contains messages, analysis context, function calls and results, proposals, approvals, applied changes and errors.",
      "API keys, configured endpoint addresses and image payloads are excluded. Provider-supplied reasoning is included when available; private reasoning that the provider does not return is unavailable. Messages/results may contain sensitive research information: inspect before sharing.",
      "A proposed action is not an applied action. Applied changes are recorded separately. Stopped or disconnected requests attempt bounded recovery in the background; later or unreachable server activity may be unavailable.",
      "Kept in browser memory until this page closes/reloads; no automatic upload. Model changes retain the log and add a reset event.",
      ...this.events.map((event, index) => `\n${index + 1}. ${event.time} — ${event.kind}\n${typeof event.data === "string" ? event.data : JSON.stringify(event.data, null, 2)}`),
      "\nEnd of current log snapshot.", ""].join("\n");
  }
}
import { version } from "../../package.json";
