import { afterEach, expect, it, vi } from "vitest";
import { getAssistantPartialTrace } from "../../src/lib/api";

afterEach(() => { vi.useRealTimers(); vi.unstubAllGlobals(); });

it.each(["connection", "body"])("bounds recovery when the %s stalls", async phase => {
  vi.useFakeTimers();
  const fetcher = vi.fn<(url: unknown, options: RequestInit) => Promise<Response>>(() => phase === "connection"
    ? new Promise<Response>(() => {})
    : Promise.resolve({ ok: true, json: () => new Promise(() => {}) } as unknown as Response));
  vi.stubGlobal("fetch", fetcher);
  const result = getAssistantPartialTrace("trace-12345");
  await vi.advanceTimersByTimeAsync(2000);
  expect(await result).toEqual([]);
  expect(fetcher.mock.calls[0][1].signal?.aborted).toBe(true);
});

it("returns recovered events without waiting for the deadline", async () => {
  const events = [{ time: "now", kind: "function returned", data: "done" }];
  vi.stubGlobal("fetch", vi.fn().mockResolvedValue({ ok: true, json: async () => ({ trace: events }) }));
  expect(await getAssistantPartialTrace("trace-12345")).toEqual(events);
});

it("keeps a recovery failure from replacing the original request outcome", async () => {
  vi.stubGlobal("fetch", vi.fn().mockResolvedValue(new Response("{}", { status: 422 })));
  expect(await getAssistantPartialTrace("trace-12345")).toEqual([]);
});
