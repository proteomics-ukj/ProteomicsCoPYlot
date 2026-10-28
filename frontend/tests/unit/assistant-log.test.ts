// @vitest-environment node
import { expect, it } from "vitest";
import { AssistantLog } from "../../src/lib/assistantLog";

it("exports immutable, full-length events with credentials and images excluded", () => {
  const log = new AssistantLog();
  log.exclude("my-secret");
  const input = { message: "x".repeat(40000), apiKey: "another-secret", nested: { text: "my-secret", plotImage: "pixels" } };
  log.record("User message", input);
  input.message = "mutated";
  log.record("Applied action", "Changed plot title.");
  const text = log.text();
  expect(text).toContain("x".repeat(40000));
  for (const secret of ["my-secret", "another-secret", "pixels", "mutated"]) expect(text).not.toContain(secret);
  expect(text.indexOf("User message")).toBeLessThan(text.indexOf("Applied action"));
});

it("redacts a newly entered key from earlier messages and retains multiple conversations", () => {
  const log = new AssistantLog();
  log.record("User message", "accidentally pasted credential-value");
  log.exclude("credential-value");
  log.record("Conversation reset", "new model");
  log.record("User message", "next request");
  expect(log.count()).toBe(3);
  expect(log.text()).not.toContain("credential-value");
  expect(log.text()).toContain("next request");
});

it("logs a changed progress line once, not every poll tick", () => {
  const log = new AssistantLog();
  log.record("Displayed assistant response / progress", "Thinking… (12 s)");
  log.record("Displayed assistant response / progress", "Thinking… (12 s)");
  log.record("Displayed assistant response / progress", "Thinking… (14 s)");
  log.record("Displayed assistant response / progress", "Thinking… (14 s)");
  log.record("Applied action", "Changed plot title.");
  expect(log.count()).toBe(3);
  expect(log.text()).toContain("(14 s)");
});


it("retains repeated actual actions and function events", () => {
  const log = new AssistantLog();
  for (let i = 0; i < 2; i++) log.record("Applied action / search result", "Added tile.");
  for (let i = 0; i < 2; i++) log.record("function started", { name: "read_results" });
  expect(log.count()).toBe(4);
});
