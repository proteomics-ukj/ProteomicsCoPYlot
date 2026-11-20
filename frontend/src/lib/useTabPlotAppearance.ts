import { useEffect, useRef, type Dispatch, type SetStateAction } from "react";
import { getAdvancedPlotMode } from "./appSettings";
import { registerStyleStore, styleTemplate, withStyle } from "./plotStyles";
import type { AdvancedPlotSettings } from "./types";

export function useAppearanceModeSync(setEnabled: (enabled: boolean) => void) {
  useEffect(() => {
    const refresh = () => setEnabled(getAdvancedPlotMode());
    window.addEventListener("focus", refresh);
    window.addEventListener("storage", refresh);
    window.addEventListener("copylot:advancedPlotModeChanged", refresh);
    return () => {
      window.removeEventListener("focus", refresh);
      window.removeEventListener("storage", refresh);
      window.removeEventListener("copylot:advancedPlotModeChanged", refresh);
    };
  }, [setEnabled]);
}

export function useAppearanceRenderDelay<T>(settings: T, setRendered: (settings: T) => void) {
  useEffect(() => {
    const timeout = window.setTimeout(() => setRendered(settings), 350);
    return () => window.clearTimeout(timeout);
  }, [settings, setRendered]);
}

type TabSettings<T extends string> = Partial<Record<T, AdvancedPlotSettings>>;

export function useTabStyleStore<T extends string>(
  storeId: string, activeTab: T, advancedByTab: TabSettings<T>,
  setAdvancedByTab: Dispatch<SetStateAction<TabSettings<T>>>,
) {
  const styleStoreState = useRef({ tab: activeTab, byTab: advancedByTab });
  useEffect(() => { styleStoreState.current = { tab: activeTab, byTab: advancedByTab }; }, [activeTab, advancedByTab]);
  useEffect(() => registerStyleStore(storeId, {
    active: () => styleStoreState.current.byTab[styleStoreState.current.tab] ?? styleTemplate() ?? {},
    setActive: (settings) => setAdvancedByTab((prev) => ({ ...prev, [styleStoreState.current.tab]: settings })),
    setAll: (style) => {
      setAdvancedByTab((prev) => Object.fromEntries(Object.entries(prev).map(
        ([tab, value]) => [tab, withStyle(value as AdvancedPlotSettings | undefined, style)],
      )) as TabSettings<T>);
      return Object.keys(styleStoreState.current.byTab).length;
    },
  }), [storeId, setAdvancedByTab]);
}
