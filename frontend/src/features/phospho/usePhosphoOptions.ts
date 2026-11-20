import { useEffect } from "react";
import { getPhosphoOptions, isAbortError } from "../../lib/api";
import type { PhosphoOptionsResponse } from "../../lib/types";

export function usePhosphoOptionsLoader(
  hasPhosphoDataset: boolean,
  applyDataSourceChoices: (options: PhosphoOptionsResponse) => void,
  setOptions: (options: PhosphoOptionsResponse | null) => void,
  setOptionsLoading: (loading: boolean) => void,
  setOptionsError: (error: string | null) => void,
) {
  useEffect(() => {
    if (!hasPhosphoDataset) {
      setOptions(null);
      setOptionsLoading(false);
      setOptionsError("No phospho dataset loaded.");
      return;
    }
    const controller = new AbortController();
    setOptionsLoading(true);
    setOptionsError(null);
    getPhosphoOptions({ signal: controller.signal })
      .then((payload) => {
        setOptions(payload);
        applyDataSourceChoices(payload);
      })
      .catch((err) => {
        if (isAbortError(err)) return;
        setOptions(null);
        setOptionsError(err instanceof Error ? err.message : "Failed to load phospho options");
      })
      .finally(() => {
        if (controller.signal.aborted) return;
        setOptionsLoading(false);
      });
    return () => {
      controller.abort();
    };
  }, [hasPhosphoDataset, applyDataSourceChoices, setOptions, setOptionsLoading, setOptionsError]);
}
