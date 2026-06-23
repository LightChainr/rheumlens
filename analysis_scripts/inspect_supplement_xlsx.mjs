import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const input = process.argv[2];
if (!input) throw new Error("usage: inspect_supplement_xlsx.mjs workbook.xlsx");
const blob = await FileBlob.load(input);
const workbook = await SpreadsheetFile.importXlsx(blob);
const overview = await workbook.inspect({
  kind: "sheet,table",
  include: "id,name,values",
  maxChars: 30000,
  tableMaxRows: 12,
  tableMaxCols: 16,
  tableMaxCellChars: 120,
});
console.log(overview.ndjson);
