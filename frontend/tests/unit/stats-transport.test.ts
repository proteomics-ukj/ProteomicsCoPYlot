import { afterEach, expect, it, vi } from "vitest";
import { runVolcanoAnalysis } from "../../src/lib/api";
import { collectColumns, rowsToQuotedCsv } from "../../src/lib/tableData";

afterEach(() => vi.unstubAllGlobals());

it("restores every stats cell and produces the same downloadable CSV", async () => {
  const columns = ["label", "sample,1", "__proto__", 'sample"2'];
  const values = [["protein α", null, "A", 1.0000000000000002], ["β", 0, "B", -3]];
  const rows = values.map(row => Object.fromEntries(columns.map((key, i) => [key, row[i]])));
  const fetch = vi.fn().mockResolvedValueOnce(new Response(JSON.stringify({
    rows: { columns, values }, untestedRows: [], warnings: ["notice"],
  }), { headers: { "Content-Type": "application/vnd.copylot.stats+json" } }));
  vi.stubGlobal("fetch", fetch);
  const result = await runVolcanoAnalysis({ kind: "protein", condition1: "A", condition2: "B" });
  expect(result.rows).toEqual(rows);
  expect(Object.keys(result.rows[0])).toEqual(columns);
  expect(result.warnings).toEqual(["notice"]);
  expect(rowsToQuotedCsv(result.rows, collectColumns(result.rows))).toEqual(rowsToQuotedCsv(rows, columns));
});
