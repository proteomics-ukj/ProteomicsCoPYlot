import { describe, expect, it } from "vitest";
import { rowsToXlsx } from "../../src/lib/xlsx";

// The workbook is a ZIP whose parts are stored uncompressed, so they can be read back as text.
function parts(bytes: Uint8Array): Record<string, string> {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const decoder = new TextDecoder();
  const found: Record<string, string> = {};
  let offset = 0;
  while (view.getUint32(offset, true) === 0x04034b50) {
    const size = view.getUint32(offset + 18, true);
    const nameLength = view.getUint16(offset + 26, true);
    const name = decoder.decode(bytes.subarray(offset + 30, offset + 30 + nameLength));
    found[name] = decoder.decode(bytes.subarray(offset + 30 + nameLength, offset + 30 + nameLength + size));
    offset += 30 + nameLength + size;
  }
  return found;
}

describe("table downloads as Excel workbooks", () => {
  const rows = [
    { Sample: "WT_1", Number: 7410, Share: 0.25, Note: "a < b & c", Flag: true },
    { Sample: "GRK_1", Number: null, Share: Number.NaN, Note: "", Flag: false },
  ];
  const bytes = rowsToXlsx(rows, ["Sample", "Number", "Share", "Note", "Flag"], "Coverage");

  it("is a workbook with one sheet", () => {
    const files = parts(bytes);
    expect(Object.keys(files)).toEqual([
      "[Content_Types].xml", "_rels/.rels", "xl/workbook.xml", "xl/_rels/workbook.xml.rels", "xl/worksheets/sheet1.xml",
    ]);
    expect(files["xl/workbook.xml"]).toContain('<sheet name="Coverage" sheetId="1" r:id="rId1"/>');
    // The directory at the end names every part.
    expect(new DataView(bytes.buffer).getUint16(bytes.length - 22 + 10, true)).toBe(5);
  });

  it("keeps numbers as numbers and text as text", () => {
    const sheet = parts(bytes)["xl/worksheets/sheet1.xml"];
    expect(sheet).toContain('<c r="A1" t="inlineStr"><is><t xml:space="preserve">Sample</t></is></c>');
    expect(sheet).toContain('<c r="B2"><v>7410</v></c>');
    expect(sheet).toContain('<c r="C2"><v>0.25</v></c>');
    expect(sheet).toContain("a &lt; b &amp; c");
    expect(sheet).toContain('<c r="E2" t="b"><v>1</v></c>');
  });

  it("leaves missing values empty and writes non-finite numbers as text", () => {
    const sheet = parts(bytes)["xl/worksheets/sheet1.xml"];
    expect(sheet).not.toContain('r="B3"');
    expect(sheet).not.toContain('r="D3"');
    expect(sheet).toContain('<c r="C3" t="inlineStr"><is><t xml:space="preserve">NaN</t></is></c>');
  });

  it("names columns past Z", () => {
    const wide = Object.fromEntries(Array.from({ length: 28 }, (_, index) => [`c${index}`, index]));
    expect(parts(rowsToXlsx([wide], Object.keys(wide)))["xl/worksheets/sheet1.xml"]).toContain('<c r="AB2"><v>27</v></c>');
  });
});
