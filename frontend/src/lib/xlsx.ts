// A table as an Excel workbook (.xlsx), written here so that a table download
// opens in a spreadsheet as it is: numbers are numbers and text is text, with
// no separator or decimal-mark guessing as in CSV. One sheet, no styles; the
// file is a ZIP of five small XML parts, stored without compression.

export const XLSX_MIME = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet";

const encoder = new TextEncoder();

const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n += 1) {
    let c = n;
    for (let k = 0; k < 8; k += 1) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

function crc32(bytes: Uint8Array): number {
  let crc = 0xffffffff;
  for (let i = 0; i < bytes.length; i += 1) crc = CRC_TABLE[(crc ^ bytes[i]) & 0xff] ^ (crc >>> 8);
  return (crc ^ 0xffffffff) >>> 0;
}

function zip(files: [string, string][]): Uint8Array {
  const parts: Uint8Array[] = [];
  const central: Uint8Array[] = [];
  let offset = 0;
  for (const [name, text] of files) {
    const nameBytes = encoder.encode(name);
    const data = encoder.encode(text);
    const crc = crc32(data);
    const header = new DataView(new ArrayBuffer(30));
    header.setUint32(0, 0x04034b50, true);
    header.setUint16(4, 20, true); // version needed
    header.setUint16(6, 0x0800, true); // UTF-8 names
    header.setUint16(8, 0, true); // stored
    header.setUint16(10, 0, true); // time
    header.setUint16(12, 0x21, true); // date: 1980-01-01
    header.setUint32(14, crc, true);
    header.setUint32(18, data.length, true);
    header.setUint32(22, data.length, true);
    header.setUint16(26, nameBytes.length, true);
    header.setUint16(28, 0, true);
    const entry = new DataView(new ArrayBuffer(46));
    entry.setUint32(0, 0x02014b50, true);
    entry.setUint16(4, 20, true);
    entry.setUint16(6, 20, true);
    entry.setUint16(8, 0x0800, true);
    entry.setUint16(10, 0, true);
    entry.setUint16(12, 0, true);
    entry.setUint16(14, 0x21, true);
    entry.setUint32(16, crc, true);
    entry.setUint32(20, data.length, true);
    entry.setUint32(24, data.length, true);
    entry.setUint16(28, nameBytes.length, true);
    entry.setUint32(42, offset, true);
    parts.push(new Uint8Array(header.buffer), nameBytes, data);
    central.push(new Uint8Array(entry.buffer), nameBytes);
    offset += 30 + nameBytes.length + data.length;
  }
  const centralSize = central.reduce((total, part) => total + part.length, 0);
  const end = new DataView(new ArrayBuffer(22));
  end.setUint32(0, 0x06054b50, true);
  end.setUint16(8, files.length, true);
  end.setUint16(10, files.length, true);
  end.setUint32(12, centralSize, true);
  end.setUint32(16, offset, true);
  const all = [...parts, ...central, new Uint8Array(end.buffer)];
  const out = new Uint8Array(all.reduce((total, part) => total + part.length, 0));
  let position = 0;
  for (const part of all) {
    out.set(part, position);
    position += part.length;
  }
  return out;
}

function escapeXml(text: string): string {
  let out = "";
  for (const character of text) {
    const code = character.codePointAt(0) ?? 0;
    // Characters XML 1.0 cannot carry.
    if ((code < 0x20 && code !== 9 && code !== 10 && code !== 13) || code === 0xfffe || code === 0xffff) continue;
    out += character === "&" ? "&amp;" : character === "<" ? "&lt;" : character === ">" ? "&gt;" : character;
  }
  return out;
}

function columnLetters(index: number): string {
  let letters = "";
  for (let n = index + 1; n > 0; n = Math.floor((n - 1) / 26)) letters = String.fromCharCode(65 + ((n - 1) % 26)) + letters;
  return letters;
}

function cell(reference: string, value: unknown): string {
  if (value == null || value === "") return "";
  if (typeof value === "number") {
    // Excel has no NaN or infinity; they stay readable as text.
    if (Number.isFinite(value)) return `<c r="${reference}"><v>${value}</v></c>`;
    value = String(value);
  }
  if (typeof value === "boolean") return `<c r="${reference}" t="b"><v>${value ? 1 : 0}</v></c>`;
  const text = typeof value === "object" ? JSON.stringify(value) : String(value);
  // A cell holds at most 32,767 characters.
  return `<c r="${reference}" t="inlineStr"><is><t xml:space="preserve">${escapeXml(text.slice(0, 32767))}</t></is></c>`;
}

/** The rows as one worksheet, with `columns` as the header row. */
export function rowsToXlsx(rows: Record<string, unknown>[], columns: string[], sheetName = "Table"): Uint8Array {
  const lines = [columns.map((column, index) => cell(`${columnLetters(index)}1`, column)).join("")];
  rows.forEach((row, rowIndex) => {
    lines.push(columns.map((column, index) => cell(`${columnLetters(index)}${rowIndex + 2}`, row[column])).join(""));
  });
  const sheet = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    + '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'
    + lines.map((line, index) => `<row r="${index + 1}">${line}</row>`).join("")
    + "</sheetData></worksheet>";
  const name = escapeXml(sheetName.replace(/[\\/?*[\]:]/g, " ").slice(0, 31) || "Table").replace(/"/g, "&quot;");
  return zip([
    ["[Content_Types].xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
      + '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
      + '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
      + '<Default Extension="xml" ContentType="application/xml"/>'
      + '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
      + '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
      + "</Types>"],
    ["_rels/.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
      + '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
      + '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
      + "</Relationships>"],
    ["xl/workbook.xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
      + '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
      + 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
      + `<sheets><sheet name="${name}" sheetId="1" r:id="rId1"/></sheets></workbook>`],
    ["xl/_rels/workbook.xml.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
      + '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
      + '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
      + "</Relationships>"],
    ["xl/worksheets/sheet1.xml", sheet],
  ]);
}
