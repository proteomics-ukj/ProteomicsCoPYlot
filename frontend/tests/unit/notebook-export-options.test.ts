import { beforeEach, describe, expect, it } from "vitest";
import {
  getNotebookExportOptions,
  notebookExportQuery,
  setNotebookExportOptions,
} from "../../src/lib/notebookExportOptions";
import { captureCopylotSessionStorageState } from "../../src/lib/restoreUiState";

describe("the notebook export options", () => {
  beforeEach(() => sessionStorage.clear());

  it("are off by default, so an export ships no data matrix", () => {
    expect(getNotebookExportOptions()).toEqual({ matrices: false, imputed: false });
    expect(notebookExportQuery()).toBe("");
  });

  it("are sent with the download and travel in a save state", () => {
    setNotebookExportOptions({ matrices: true, imputed: false });
    expect(notebookExportQuery()).toBe("?matrices=true");
    setNotebookExportOptions({ matrices: true, imputed: true });
    expect(notebookExportQuery()).toBe("?matrices=true&imputed=true");
    expect(Object.keys(captureCopylotSessionStorageState())).toContain("pcopylot.ui.notebookExport.v1");
  });

  it("read an unreadable stored value as off", () => {
    sessionStorage.setItem("pcopylot.ui.notebookExport.v1", "{broken");
    expect(getNotebookExportOptions()).toEqual({ matrices: false, imputed: false });
  });
});
