import { act, fireEvent, render, screen } from "@testing-library/react";
import { expect, it, vi } from "vitest";
import AssistantLogDownload from "../../src/components/AssistantLogDownload";
import { AssistantLog } from "../../src/lib/assistantLog";

it("previews text safely and downloads the complete current log", async () => {
  const log = new AssistantLog();
  render(<AssistantLogDownload log={log} />);
  expect((screen.getByRole("button", { name: "Download conversation log" }) as HTMLButtonElement).disabled).toBe(true);
  act(() => log.record("User message", "<img src=x onerror=alert(1)>"));
  const details = screen.getByText(/Conversation log \(1 events\)/).closest("details")!;
  details.open = true;
  fireEvent(details, new Event("toggle"));
  expect(screen.getByLabelText("Conversation log preview").textContent).toContain("<img src=x");
  expect(details.querySelector("img")).toBeNull();
  const create = vi.fn(() => "blob:log");
  vi.stubGlobal("URL", { createObjectURL: create, revokeObjectURL: vi.fn() });
  const click = vi.spyOn(HTMLAnchorElement.prototype, "click").mockImplementation(() => {});
  fireEvent.click(screen.getByRole("button", { name: "Download conversation log" }));
  expect(create.mock.calls).toHaveLength(1);
  const blob = (create.mock.calls as unknown as [Blob][])[0][0];
  const content = await new Promise(resolve => {
    const reader = new FileReader();
    reader.onload = () => resolve(reader.result);
    reader.readAsText(blob);
  });
  expect(content).toBe(log.text());
  expect(click.mock.instances[0].download).toMatch(/^copylot-ai-conversation-.*\.txt$/);
  vi.unstubAllGlobals();
});
